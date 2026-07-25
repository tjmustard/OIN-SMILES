"""Lane 3 diagnostic: attribute every eta winding marker to the tier that decided it.

Answers, per molecule and per haptic slot:

* which of the three heading-atom tiers fired (content-canonical / geometric fallback /
  symmetric-ligand override);
* what ``_orientation_free_from_mol`` said (True / False / **None** = undecidable);
* whether ``_canonical_heading_atom`` / ``_topological_heading_atom`` could resolve;
* the emitted winding character and the metal->centroid axis it was measured against.

It also runs the two experiments that decide whether a winding difference is *encoder
drift* or *a real structural difference*:

``--rotate N``   re-encode the input under N random **proper** rotations. The winding sign
                 math (``cross(v_star, v_next) . axis``) is rotation-equivariant, so a
                 change here is a genuine embedding dependence -- the thing Lane 3 exists
                 to remove.
``--mirror``     re-encode the **reflected** input. Winding is *designed* to flip for a
                 mirror image (that is the stereochemistry it carries), so if a
                 generated structure's string equals the mirrored input's string, the
                 generator produced the enantiomer and the encoder reported it correctly.

Usage:
    PYTHONPATH=src python tools/winding_diag.py --xyz A.xyz [--xyz B.xyz] --rotate 4 --mirror
    PYTHONPATH=src python tools/winding_diag.py --results-dir DIR --only MOL1,MOL2
"""

import argparse
import contextlib
import glob
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.utils import oin_aligner as AL  # noqa: E402

SLOT_RE = re.compile(r"\{(\d+)([<>^]?)\}")


@contextlib.contextmanager
def _silence_fds():
    with open(os.devnull, "w") as devnull:
        old_out, old_err = os.dup(1), os.dup(2)
        try:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            os.dup2(old_out, 1)
            os.dup2(old_err, 2)
            os.close(old_out)
            os.close(old_err)


# --------------------------------------------------------------------------- tracing


TRACE: list = []


def install_trace():
    """Wrap the aligner's heading/winding helpers so each decision records its tier."""
    cls = AL.OINDiscreteAligner

    orig_permute = cls._permute_and_serialize
    orig_topological = cls._topological_heading_atom.__func__
    orig_canonical = cls._canonical_heading_atom.__func__
    orig_winding = cls._determine_winding

    def permute(self, slot_assignment, tmpl_vectors, geometry_name=None, alignment_rotation=None):
        rec = {
            "geometry": geometry_name,
            "gated_out": not (
                geometry_name
                and alignment_rotation is not None
                and geometry_name in AL.TEMPLATE_SPECS
            ),
            "items": [],
            "winding_calls": [],
        }
        TRACE.append(rec)
        if slot_assignment:
            # Identical fragments share their local atom indices, so an item is only
            # identified by (slot, constituent set) -- never by the index set alone.
            for slot_idx, atom in enumerate(slot_assignment):
                if atom is None:
                    continue
                cons = atom.get("constituent_indices") or [atom.get("local_idx")]
                if len(cons) < 2:
                    continue
                smiles = atom["chem_id"][1]
                free_raw = AL._winding_is_orientation_free(smiles, tuple(sorted(cons)))
                rec["items"].append(
                    {
                        "pre_slot": slot_idx,
                        "rank": atom.get("rank"),
                        "constituents": sorted(cons),
                        "smiles": smiles,
                        "n_constituents": len(cons),
                        "in_SYMMETRIC_LIGANDS": smiles in AL.SYMMETRIC_LIGANDS,
                        "orientation_free_raw": free_raw,
                        "item_orientation_free": bool(free_raw),
                        "strict_canonical_heading": orig_canonical(cls, smiles, cons),
                        "topological_heading": orig_topological(cls, smiles, cons),
                        "ring_signature": cls._canonical_ring_signature(smiles),
                    }
                )
        return orig_permute(self, slot_assignment, tmpl_vectors, geometry_name, alignment_rotation)

    def determine_winding(
        self, grp_coords, star_idx, constituent_indices, slot_z, slot_x_ref, alignment_rotation=None
    ):
        out = orig_winding(
            self, grp_coords, star_idx, constituent_indices, slot_z, slot_x_ref, alignment_rotation
        )
        if TRACE:
            axis = np.asarray(grp_coords, dtype=float).mean(axis=0)
            cons = sorted(constituent_indices)
            match = next(
                (it for it in TRACE[-1]["items"] if it["constituents"] == cons),
                None,
            )
            # Mirrors the tier routing in `_permute_and_serialize` step 4: every eta
            # group now resolves topologically, and `_topological_heading_atom` reports
            # WHICH of its own three sub-tiers answered.
            if match is None:
                tier = "?"
            elif match["in_SYMMETRIC_LIGANDS"]:
                tier = "3-symmetric-override"
            elif match["strict_canonical_heading"] is not None:
                tier = "1a-strict-canonical-rank"
            elif match["topological_heading"] is not None:
                tier = "1b-topological-fallback"
            else:
                tier = "2-GEOMETRIC-FALLBACK"
            TRACE[-1]["winding_calls"].append(
                {
                    "constituents": cons,
                    "heading": star_idx,
                    "winding": out,
                    "tier": tier,
                    "axis": [round(float(v), 4) for v in axis],
                    "ring_signature": (match or {}).get("ring_signature"),
                    "orientation_free_raw": (match or {}).get("orientation_free_raw"),
                    "strict_canonical_heading": (match or {}).get("strict_canonical_heading"),
                    "topological_heading": (match or {}).get("topological_heading"),
                    "in_SYMMETRIC_LIGANDS": (match or {}).get("in_SYMMETRIC_LIGANDS"),
                }
            )
        return out

    cls._permute_and_serialize = permute
    cls._determine_winding = determine_winding


# --------------------------------------------------------------------------- geometry


def read_xyz(path):
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    comment = lines[1] if len(lines) > 1 else ""
    syms, xyz = [], []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        syms.append(parts[0])
        xyz.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return syms, np.array(xyz), comment


def write_xyz(path, syms, xyz, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n{comment}\n")
        for s, (x, y, z) in zip(syms, xyz):
            f.write(f"{s} {x:.10f} {y:.10f} {z:.10f}\n")


def random_rotation(rng):
    """Uniform random proper rotation (QR of a Gaussian matrix, det forced +1)."""
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


# --------------------------------------------------------------------------- reporting


def winding_profile(oin):
    """Ordered [(slot, char)] for every slot marker carrying a winding character."""
    return [(int(m.group(1)), m.group(2)) for m in SLOT_RE.finditer(oin or "") if m.group(2)]


def encode(conv, path):
    TRACE.clear()
    try:
        with _silence_fds():
            oin = conv.convert(path)
    except Exception as e:  # noqa: BLE001
        return None, f"{type(e).__name__}: {e}", []
    return oin, None, [dict(r) for r in TRACE]


def summarize_trace(trace):
    """The winding decisions of the LAST _permute_and_serialize call (the emitted one)."""
    if not trace:
        return []
    return trace[-1]["winding_calls"]


def run_one(conv, name, paths, rotate, mirror, scratch, rng):
    print("=" * 100)
    print(name)
    results = {}
    for label, path in paths.items():
        oin, err, trace = encode(conv, path)
        results[label] = oin
        print(f"  [{label}] {'ERROR ' + err if err else oin}")
        print(f"      winding: {winding_profile(oin)}")
        for s in summarize_trace(trace):
            print(
                f"      eta n={len(s['constituents']):2d} tier={s['tier']} "
                f"free_raw={s['orientation_free_raw']} sym={s['in_SYMMETRIC_LIGANDS']} "
                f"strict={s['strict_canonical_heading']} topo={s['topological_heading']} "
                f"head={s['heading']} wind={s['winding']}  ring={s['ring_signature']}"
            )

    base_path = paths.get("input") or next(iter(paths.values()))
    syms, xyz, comment = read_xyz(base_path)

    if rotate:
        rot_strings = set()
        rot_windings = set()
        for i in range(rotate):
            R = random_rotation(rng)
            p = os.path.join(scratch, f"{name}_rot{i}.xyz")
            write_xyz(p, syms, xyz @ R.T, comment)
            oin, err, _ = encode(conv, p)
            rot_strings.add(oin if not err else f"ERR:{err}")
            rot_windings.add(tuple(winding_profile(oin)))
        print(
            f"  ROTATION x{rotate}: distinct strings={len(rot_strings)} "
            f"distinct winding profiles={len(rot_windings)}"
        )
        for w in sorted(rot_windings, key=str):
            print(f"      {w}")
        results["_rotation_invariant"] = len(rot_strings) == 1
        results["_rotation_winding_invariant"] = len(rot_windings) == 1
        results["_rotation_matches_input"] = rot_strings == {results.get("input")}

    if mirror:
        p = os.path.join(scratch, f"{name}_mirror.xyz")
        mx = xyz.copy()
        mx[:, 0] *= -1.0
        write_xyz(p, syms, mx, comment)
        oin, err, _ = encode(conv, p)
        print(f"  MIRRORED input: {'ERROR ' + err if err else oin}")
        print(f"      winding: {winding_profile(oin)}")
        results["mirror"] = oin
        if "generated" in results:
            print(f"      mirror == generated ? {oin == results['generated']}")
            print(
                "      mirror winding == generated winding ? "
                f"{winding_profile(oin) == winding_profile(results.get('generated'))}"
            )
        if "input" in results:
            print(f"      mirror == input ? {oin == results['input']}")
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--xyz", action="append", default=[], help="XYZ file (repeatable)")
    ap.add_argument("--results-dir", help="Sweep dir: use each report's input_xyz + generated")
    ap.add_argument("--only", help="Comma-separated molecule ids (with --results-dir)")
    ap.add_argument("--rotate", type=int, default=0, help="N random proper rotations of the input")
    ap.add_argument("--mirror", action="store_true", help="Also encode the reflected input")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--scratch", default=None, help="Where to write rotated/mirrored XYZ")
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args()

    scratch = args.scratch or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", ".winding_diag"
    )
    os.makedirs(scratch, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    install_trace()
    conv = XYZToSMILES()

    jobs = []
    for p in args.xyz:
        jobs.append((os.path.basename(p), {"input": p}))
    if args.results_dir:
        src = os.path.abspath(args.results_dir)
        wanted = {m.strip() for m in (args.only or "").split(",") if m.strip()}
        for rp in sorted(glob.glob(os.path.join(src, "individual_reports", "*.json"))):
            with open(rp) as f:
                rep = json.load(f)
            mol = rep.get("molecule")
            if wanted and mol not in wanted:
                continue
            paths = {}
            if rep.get("input_xyz") and os.path.exists(rep["input_xyz"]):
                paths["input"] = rep["input_xyz"]
            for cand in (
                os.path.join(src, "structures", f"{mol}_generated.xyz"),
                os.path.join(src, "test_failures", mol, "last_generated.xyz"),
                rep.get("generated_xyz") or "",
            ):
                if cand and os.path.exists(cand):
                    paths["generated"] = cand
                    break
            if paths:
                jobs.append((mol, paths))

    out = {}
    for name, paths in jobs:
        out[name] = run_one(conv, name, paths, args.rotate, args.mirror, scratch, rng)

    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"\nWrote {args.json_out}")


if __name__ == "__main__":
    main()
