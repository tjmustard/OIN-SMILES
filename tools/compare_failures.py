"""Build input-vs-generated overlay XYZ files for round-trip sweep failures.

For each failed molecule in a sweep results directory that has a *saved* generated
structure, this writes a two-frame XYZ (frame 1 = original input, frame 2 = generated
"failed" structure). Both frames are translated so the metal sits at the origin and
rotated into the generator's canonical coordination frame for that molecule's geometry
(SPL / OCT / TET / TBP / TPL / TPY / SPY / LIN / SQA / PBP / TCT), so every file opens
in a consistent, comparable orientation. The generated frame is additionally snapped
onto the input so a correct pair overlays and a wrong pair visibly diverges -- which is
exactly what makes the accuracy easy to eyeball in a GUI (VMD/ASE/Avogadro/Jmol read a
multi-frame XYZ as a short trajectory).

It reuses the sweep's on-disk generated structures (it does NOT re-run the generator),
so the geometry shown is exactly the one the sweep produced, and the result is
deterministic. Failures that crashed before any 3D existed (encode/generation stage) are
skipped and recorded in ``skipped_manifest.json``.

Reuses the encoder's own geometry machinery: ``TEMPLATES`` (ideal axis positions),
``classify_coordination_geometry``, and ``OINDiscreteAligner._map_to_template`` from
``oinsmiles.utils.oin_aligner``.
"""

import argparse
import glob
import itertools
import json
import os
import re
import sys

import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

# --- Reuse the encoder's geometry templates + matcher (see src/oinsmiles/utils/oin_aligner.py)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

from rdkit import Chem  # noqa: E402

from oinsmiles.core.constants import TRANSITION_METALS  # noqa: E402
from oinsmiles.utils.oin_aligner import (  # noqa: E402
    TEMPLATES,
    OINDiscreteAligner,
    classify_coordination_geometry,
)

_PT = Chem.GetPeriodicTable()
_ALIGNER = OINDiscreteAligner(0, [])  # pure w.r.t. instance state; used for _map_to_template
_TM_SYMBOLS = set(TRANSITION_METALS)
_DONOR_CUTOFF_FACTOR = 1.3  # bonded if dist < factor * (Rcov_metal + Rcov_atom)


# --------------------------------------------------------------------------- I/O
def parse_xyz(path):
    """Parse a standard XYZ file -> (symbols: list[str], coords: (N,3) float array)."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    symbols, coords = [], []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.asarray(coords, dtype=float)


def write_multiframe_xyz(path, frames):
    """Write frames = [(comment, symbols, coords), ...] as one multi-frame XYZ."""
    chunks = []
    for comment, symbols, coords in frames:
        chunks.append(str(len(symbols)))
        chunks.append(comment)
        for sym, (x, y, z) in zip(symbols, coords):
            chunks.append(f"{sym:<2s} {x:14.8f} {y:14.8f} {z:14.8f}")
    with open(path, "w") as fh:
        fh.write("\n".join(chunks) + "\n")


# ------------------------------------------------------------------- geometry glue
def _rcov(symbol):
    return _PT.GetRcovalent(_PT.GetAtomicNumber(symbol))


def find_metal_index(symbols, coords, metal_symbol=None):
    """Return the metal atom index.

    Prefers ``metal_symbol`` (the element named in the OIN tag, authoritative). Falls
    back to any transition-metal element. With several candidates, picks the one with
    the most close contacts (the coordination center). Returns None if none found.
    """
    if metal_symbol:
        cand = [i for i, s in enumerate(symbols) if s == metal_symbol]
    else:
        cand = [i for i, s in enumerate(symbols) if s in _TM_SYMBOLS]
    if not cand:
        return None
    if len(cand) == 1:
        return cand[0]
    # Multiple metals: choose the one with the most neighbours within 2.8 A.
    dmat = cdist(coords[cand], coords)
    return cand[int(np.argmax((dmat < 2.8).sum(axis=1)))]


def detect_donors(symbols, coords, metal_idx):
    """Indices of atoms covalently bonded to the metal (covalent-radius cutoff).

    Used only as a diagnostic and as the no-geometry-tag fallback; the primary donor
    set for alignment is ``closest_heavy_indices`` (see build_overlay).
    """
    rm = _rcov(symbols[metal_idx])
    mpos = coords[metal_idx]
    donors = []
    for i, (sym, pos) in enumerate(zip(symbols, coords)):
        if i == metal_idx:
            continue
        if np.linalg.norm(pos - mpos) < _DONOR_CUTOFF_FACTOR * (rm + _rcov(sym)):
            donors.append(i)
    return donors


def closest_heavy_indices(symbols, coords, metal_idx, n):
    """Indices of the n non-hydrogen atoms closest to the metal (the primary sphere).

    The coordination number is authoritative from the OIN geometry (slot count), so we
    take exactly that many nearest heavy atoms as donors. This keeps the input and
    generated structures aligned to the *same* template consistently, and handles
    haptic ligands sensibly (one nearest ring atom per coordinated face). Metal-centred
    coords assumed only for distance ordering, which is translation-invariant anyway.
    """
    heavy = [i for i, s in enumerate(symbols) if i != metal_idx and s != "H"]
    if not heavy:
        return []
    dist = np.linalg.norm(coords[heavy] - coords[metal_idx], axis=1)
    return [heavy[k] for k in np.argsort(dist)[:n]]


def _pai_rotation(heavy_coords):
    """Deterministic principal-axis rotation (fallback when no template fits).

    heavy_coords are metal-centred. Returns a Rotation whose frame is the inertia
    principal axes, ordered by eigenvalue and sign-fixed for stability. Residual
    axis/sign ambiguity is resolved later by the discrete refine step.
    """
    if len(heavy_coords) < 3:
        return Rotation.from_matrix(np.eye(3))
    tensor = heavy_coords.T @ heavy_coords
    _, vecs = np.linalg.eigh(tensor)  # ascending eigenvalues
    axes = vecs[:, ::-1]  # largest spread first
    # Deterministic sign: make the largest-magnitude projection positive per axis.
    for k in range(3):
        proj = heavy_coords @ axes[:, k]
        if proj[int(np.argmax(np.abs(proj)))] < 0:
            axes[:, k] *= -1
    if np.linalg.det(axes) < 0:
        axes[:, 2] *= -1
    return Rotation.from_matrix(axes.T)  # rows are principal axes -> maps world into frame


def canonical_rotation(coords, symbols, donor_indices, metal_idx, geo_code):
    """Rotation putting a metal-centred structure into its canonical geometry frame.

    Returns (Rotation, method) where method is "template" (donors snapped onto the
    geometry's ideal axes) or "pai" (principal-axis fallback for haptic / CN-mismatch
    cases where the detected donor count doesn't fit a template).
    """
    donor_vecs = coords[donor_indices] - coords[metal_idx]
    for code in (geo_code, classify_coordination_geometry(donor_vecs) if donor_indices else None):
        template = TEMPLATES.get(code) if code else None
        if template is not None and len(donor_indices) == len(template):
            virtual = [{"coords": v} for v in donor_vecs]
            _, _, rot = _ALIGNER._map_to_template(virtual, np.asarray(template, dtype=float))
            return rot, "template"
    heavy = coords[[i for i, s in enumerate(symbols) if s != "H"]] - coords[metal_idx]
    return _pai_rotation(heavy), "pai"


def _octahedral_rotations():
    """The 24 proper rotations of the octahedron (signed permutation matrices, det +1)."""
    rots = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            m = np.zeros((3, 3))
            for col, (row, sign) in enumerate(zip(perm, signs)):
                m[row, col] = sign
            if abs(np.linalg.det(m) - 1.0) < 1e-9:
                rots.append(Rotation.from_matrix(m))
    return rots


def _build_candidate_rotations():
    """Discrete rotation set for snapping the generated frame onto the input frame.

    The 24 octahedral rotations cover axis relabelling / sign ambiguity (OCT, SPL, TET,
    SQA, PBP); a 15-degree z-axis sweep covers the principal-axis (C3/C4/C5) ambiguity of
    trigonal / axial geometries (TBP, TPL, TPY, SPY, LIN). Together they resolve the
    template's residual symmetry without per-geometry point-group tables.
    """
    cands = _octahedral_rotations()
    for k in range(1, 24):
        cands.append(Rotation.from_euler("z", k * 15.0, degrees=True))
    return cands


_CANDIDATE_ROTATIONS = _build_candidate_rotations()


def _element_aware_chamfer(a_coords, a_syms, b_coords, b_syms):
    """Mean nearest-neighbour distance from A to B, matching only like elements."""
    total, count = 0.0, 0
    a_syms = np.asarray(a_syms)
    b_syms = np.asarray(b_syms)
    for el in set(a_syms):
        a_sel = a_coords[a_syms == el]
        b_sel = b_coords[b_syms == el]
        if len(b_sel) == 0:
            total += 1e3 * len(a_sel)  # element present in gen but not input -> large penalty
        else:
            total += cdist(a_sel, b_sel).min(axis=1).sum()
        count += len(a_sel)
    return total / max(count, 1)


def refine_onto_input(gen_coords, gen_syms, input_coords, input_syms):
    """Pick the candidate rotation of the generated frame that best overlays the input.

    Both are already canonically oriented and metal-centred; this removes the residual
    template-symmetry ambiguity by minimising an element-aware Chamfer over heavy atoms.
    """
    gh = gen_syms != "H" if isinstance(gen_syms, np.ndarray) else np.array(gen_syms) != "H"
    ih = np.array(input_syms) != "H"
    g_heavy, g_hsym = gen_coords[gh], np.asarray(gen_syms)[gh]
    i_heavy, i_hsym = input_coords[ih], np.asarray(input_syms)[ih]
    if len(g_heavy) == 0 or len(i_heavy) == 0:
        return Rotation.from_matrix(np.eye(3))
    best_rot, best_score = _CANDIDATE_ROTATIONS[0], float("inf")
    for rot in _CANDIDATE_ROTATIONS:
        score = _element_aware_chamfer(rot.apply(g_heavy), g_hsym, i_heavy, i_hsym)
        if score < best_score:
            best_score, best_rot = score, rot
    return best_rot


# ------------------------------------------------------------------- report parsing
_GEO_RE = re.compile(r"^\[([A-Za-z]+)_([A-Za-z0-9]+)\]")


def parse_metal_and_geo(oin_string):
    """(metal_symbol, geo_code) from an OIN string's leading metal tag, or (None, None)."""
    if not oin_string:
        return None, None
    m = _GEO_RE.match(oin_string.strip())
    if not m:
        return None, None
    return m.group(1), m.group(2)


def classify_failure(error):
    """Coarse failure class from the report's error string."""
    e = error or ""
    for needle, label in (
        ("String mismatch", "string_mismatch"),
        ("Atom count mismatch", "atom_count_mismatch"),
        ("High RMSD", "high_rmsd"),
        ("RMSD mapping failed", "rmsd_mapping_failed"),
        ("TimeoutException", "timeout"),
        ("XYZToSMILES failed", "encode_failed"),
        ("Generation/Verification failed", "generation_failed"),
        ("UncoordinatedFragment", "uncoordinated_fragment"),
    ):
        if needle in e:
            return label
    return "other" if e else "unknown"


# ------------------------------------------------------------------------ pipeline
def build_overlay(report, results_dir):
    """Build the two aligned frames for one failed molecule.

    Returns (frames, row) on success or (None, row) on failure/skip; ``row`` always
    carries diagnostic fields for the manifest.
    """
    mol = report["molecule"]
    row = {
        "molecule": mol,
        "failure_class": classify_failure(report.get("error")),
        "rmsd": (report.get("metrics") or {}).get("rmsd"),
        "smiles_1": report.get("smiles_1"),
        "smiles_2": report.get("smiles_2"),
    }

    input_path = report.get("input_xyz")
    gen_path = os.path.join(results_dir, "structures", f"{mol}_generated.xyz")
    if not input_path or not os.path.exists(input_path) or not os.path.exists(gen_path):
        row["error"] = "missing input or generated xyz"
        return None, row

    in_syms, in_coords = parse_xyz(input_path)
    gen_syms, gen_coords = parse_xyz(gen_path)

    metal_sym, geo_in = parse_metal_and_geo(report.get("smiles_1"))
    _, geo_gen = parse_metal_and_geo(report.get("smiles_2"))

    in_metal = find_metal_index(in_syms, in_coords, metal_sym)
    gen_metal = find_metal_index(gen_syms, gen_coords, metal_sym)
    if in_metal is None or gen_metal is None:
        row["error"] = "metal not found"
        return None, row

    # Recentre both so the metal is at the origin.
    in_coords = in_coords - in_coords[in_metal]
    gen_coords = gen_coords - gen_coords[gen_metal]

    # Covalent detection is a diagnostic and the fallback for missing geometry tags.
    cov_in = detect_donors(in_syms, in_coords, in_metal)
    cov_gen = detect_donors(gen_syms, gen_coords, gen_metal)
    geo_code = geo_in or classify_coordination_geometry(in_coords[cov_in] - in_coords[in_metal])

    # Primary donor set: the geometry's slot-count nearest heavy atoms, so input and
    # generated snap to the SAME template consistently.
    slots = len(TEMPLATES[geo_code]) if geo_code in TEMPLATES else None
    if slots:
        in_donors = closest_heavy_indices(in_syms, in_coords, in_metal, slots)
        gen_donors = closest_heavy_indices(gen_syms, gen_coords, gen_metal, slots)
    else:
        in_donors, gen_donors = cov_in, cov_gen

    r_in, m_in = canonical_rotation(in_coords, in_syms, in_donors, in_metal, geo_code)
    in_coords = r_in.apply(in_coords)

    r_gen, m_gen = canonical_rotation(gen_coords, gen_syms, gen_donors, gen_metal, geo_code)
    gen_coords = r_gen.apply(gen_coords)

    # Snap the generated frame onto the (already canonical) input frame.
    gen_coords = refine_onto_input(gen_coords, gen_syms, in_coords, in_syms).apply(gen_coords)

    method = "template" if (m_in == "template" and m_gen == "template") else "fallback"
    row.update(
        {
            "geom_in": geo_in,
            "geom_gen": geo_gen,
            "geom_used": geo_code,
            "donors_used": len(in_donors),
            "slots_expected": slots,
            "cov_cn_in": len(cov_in),
            "cov_cn_gen": len(cov_gen),
            "alignment": method,
            "n_atoms_in": len(in_syms),
            "n_atoms_gen": len(gen_syms),
        }
    )
    rmsd_str = f"{row['rmsd']:.3f}" if isinstance(row["rmsd"], (int, float)) else "n/a"
    frames = [
        (
            f"INPUT {mol} | geom={geo_code} | {os.path.basename(input_path)}",
            in_syms,
            in_coords,
        ),
        (
            f"GENERATED {mol} | geom_in={geo_in} geom_gen={geo_gen} | "
            f"fail={row['failure_class']} | rmsd={rmsd_str} | align={method}",
            gen_syms,
            gen_coords,
        ),
    ]
    return frames, row


# ------------------------------------------------------------------- selection glue
def load_failed_reports(results_dir):
    """All individual reports with status == 'failed', sorted by molecule name."""
    reports = []
    for path in sorted(glob.glob(os.path.join(results_dir, "individual_reports", "*.json"))):
        with open(path) as fh:
            rep = json.load(fh)
        if rep.get("status") == "failed":
            reports.append(rep)
    return reports


def has_structure(results_dir, mol):
    return os.path.exists(os.path.join(results_dir, "structures", f"{mol}_generated.xyz"))


def diverse_sample(reports, n):
    """Pick n reports spread across (geometry, failure-class) groups, deterministically."""
    groups = {}
    for rep in reports:
        _, geo = parse_metal_and_geo(rep.get("smiles_1"))
        key = (geo, classify_failure(rep.get("error")))
        groups.setdefault(key, []).append(rep)
    ordered_keys = sorted(groups, key=lambda k: (str(k[0]), str(k[1])))
    picked, idx = [], 0
    while len(picked) < n and any(groups[k] for k in ordered_keys):
        key = ordered_keys[idx % len(ordered_keys)]
        if groups[key]:
            picked.append(groups[key].pop(0))
        idx += 1
    return picked[:n]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--results-dir",
        type=str,
        default="../tmCAT-tmPHOTO_xyz_dataset/results-v042-quick-rerun",
        help="Sweep results dir (must contain individual_reports/ and structures/).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="../results-v042-failure-overlays",
        help="Directory to write <mol>_compare.xyz + manifest.json + skipped_manifest.json.",
    )
    parser.add_argument(
        "--only",
        type=str,
        default=None,
        help="Comma-separated molecule names (e.g. ATUROX_comp_0) to process, bypassing "
        "sampling/limit filters. For targeted repro of specific cases.",
    )
    parser.add_argument(
        "--sample",
        type=int,
        default=None,
        help="Process N failures spread across geometry types and failure classes "
        "(deterministic). Use for the initial format-QA batch, e.g. --sample 10.",
    )
    parser.add_argument("--limit", type=int, default=None, help="Cap number of molecules.")
    parser.add_argument(
        "--shard",
        type=str,
        default=None,
        help="I:N — process only the I-th of N deterministic slices of the molecule list "
        "(e.g. --shard 2:5). For parallel workers.",
    )
    args = parser.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    failed = load_failed_reports(results_dir)
    eligible = [r for r in failed if has_structure(results_dir, r["molecule"])]
    skipped = [
        {"molecule": r["molecule"], "failure_class": classify_failure(r.get("error"))}
        for r in failed
        if not has_structure(results_dir, r["molecule"])
    ]
    print(
        f"{len(failed)} failed; {len(eligible)} with a saved structure (eligible), "
        f"{len(skipped)} without (skipped)."
    )

    if args.only:
        wanted = {n.strip().removesuffix(".xyz") for n in args.only.split(",") if n.strip()}
        eligible = [r for r in eligible if r["molecule"] in wanted]
        print(f"Only mode: matched {len(eligible)} of {len(wanted)} requested molecules.")
    elif args.sample:
        eligible = diverse_sample(eligible, args.sample)
        print(f"Sample mode: selected {len(eligible)} diverse molecules.")

    if args.shard:
        try:
            shard_i, shard_n = (int(x) for x in args.shard.split(":"))
        except ValueError:
            parser.error(f"--shard must be I:N (got {args.shard!r})")
        if not (1 <= shard_i <= shard_n):
            parser.error(f"--shard index out of range: {args.shard!r}")
        eligible = eligible[shard_i - 1 :: shard_n]
        print(f"Shard {shard_i}/{shard_n}: {len(eligible)} molecules in this slice.")

    if args.limit:
        eligible = eligible[: args.limit]

    manifest, errors = [], []
    for rep in eligible:
        try:
            frames, row = build_overlay(rep, results_dir)
        except Exception as exc:  # keep the batch going; record the failure
            errors.append({"molecule": rep["molecule"], "error": f"{type(exc).__name__}: {exc}"})
            continue
        if frames is None:
            errors.append(row)
            continue
        write_multiframe_xyz(os.path.join(output_dir, f"{rep['molecule']}_compare.xyz"), frames)
        manifest.append(row)

    with open(os.path.join(output_dir, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    with open(os.path.join(output_dir, "skipped_manifest.json"), "w") as fh:
        json.dump({"no_structure": skipped, "build_errors": errors}, fh, indent=2)

    n_fallback = sum(1 for r in manifest if r.get("alignment") == "fallback")
    print(
        f"Wrote {len(manifest)} overlay files to {output_dir} "
        f"({n_fallback} used the PAI fallback, {len(errors)} build errors)."
    )


if __name__ == "__main__":
    main()
