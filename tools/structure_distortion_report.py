"""Quantify how physically distorted MetalloGen's generated 3D structures are.

Research tool seeding v0.4.3 (structure-quality wave). For every generated structure in a
sweep results dir it scores physical distortion from four independent angles, then rolls
them into a single MPO (multi-parameter optimization) quality score in [0, 1] so structures
are rankable and the worst cases surface for prioritization. Read-only: it reuses the
sweep's on-disk structures and reports, never re-runs the generator, and never edits
MetalloGen.

The four angles (each metric mapped to a Derringer-Suich desirability in [0,1] declared as
minimize / maximize / target-range with bounds; see METRIC_SPECS):
  * 3D geometry  - steric clashes (vdW-overlap + MetalloGen's own covalent gate), bond-length
                   deviation (metal-ligand vs bond_lengths table, organic vs covalent sum),
                   bond-angle strain.
  * graph        - perceived coordination number vs the input's, perceived bond-count vs the
                   input's, perception health (did the sweep re-encode succeed).
  * reference    - coordination-sphere RMSD to input (stored) + full heavy-atom overlay
                   divergence to input (exposes ligand-body distortion the sphere RMSD misses).
  * relax proxy  - ligands-only UFF relaxation displacement + energy drop ("how far would it
                   move under optimization"); rdkit-only, metal-free, best-effort.

Cohorts: the generated structures (headline) and the real input crystal structures (control,
for calibration/contrast). Bounds and weights live in DEFAULT_CONFIG, overridable via
--config <json>. Emits distortion_metrics.json (per molecule) + distortion_report.md.
"""

import argparse
import copy
import glob
import json
import math
import os
import sys

import numpy as np
from scipy.spatial.distance import pdist, squareform

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../tests/integration")))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))  # sibling tools (compare_failures)

import compare_failures as cf  # noqa: E402  (reuse its geometry primitives)
from rdkit import (
    Chem,  # noqa: E402
    RDLogger,  # noqa: E402
)
from rdkit.Chem import AllChem  # noqa: E402

from oinsmiles.core.constants import TRANSITION_METALS  # noqa: E402
from oinsmiles.generation.oin_parser import OINParser  # noqa: E402
from oinsmiles.generator3d.bond_lengths import bond_length  # noqa: E402
from oinsmiles.utils.oin_aligner import TEMPLATES  # noqa: E402

RDLogger.DisableLog("rdApp.*")
_PT = Chem.GetPeriodicTable()
_TM = set(TRANSITION_METALS)


# --------------------------------------------------------------------- configuration
# Each metric: aggregation angle, desirability direction, and bounds L (=ideal, d=1) -> U
# (=unacceptable, d=0). Bounds are calibrated to the real-input population (see the docs
# report) plus physical caps; override via --config <json>. "min"/"max" are one-sided ramps.
EPS = 0.01  # desirability floor at MPO aggregation so no single angle zeroes the score
# `informational` metrics are computed and reported but EXCLUDED from the MPO. Two were
# demoted after sample QA proved they do not discriminate MetalloGen distortion (see the docs
# report): `bondlen_dev` is inverted (the generator places near-ideal bond lengths, so it
# scores *better* than real crystal inputs — bonds are not the distortion), and the ligands-only
# UFF `relax` proxy is confounded (an all-single-bond, metal-free FF relaxes real inputs just as
# far, so it measures bond-order-model mismatch, not distortion; a valid version needs xtb/MACE).
DEFAULT_CONFIG = {
    "angle_weights": {"geometry": 0.45, "graph": 0.30, "reference": 0.25},
    "metrics": {
        # 3D geometry (scored)
        "clash_vdw": {"angle": "geometry", "dir": "min", "L": 0.0, "U": 6.0, "w": 1.0},
        "clash_severe": {"angle": "geometry", "dir": "min", "L": 0.0, "U": 2.0, "w": 1.0},
        "worst_overlap": {"angle": "geometry", "dir": "max", "L": 0.60, "U": 0.90, "w": 1.0},
        "angle_strain": {"angle": "geometry", "dir": "min", "L": 4.0, "U": 20.0, "w": 1.0},
        # graph / connectivity (scored)
        "cn_divergence": {"angle": "graph", "dir": "min", "L": 0.0, "U": 2.0, "w": 1.0},
        "bondcount_divergence": {"angle": "graph", "dir": "min", "L": 0.0, "U": 6.0, "w": 1.0},
        "perception_ok": {"angle": "graph", "dir": "max", "L": 0.0, "U": 1.0, "w": 1.0},
        # reference vs input (scored)
        "coord_rmsd": {"angle": "reference", "dir": "min", "L": 0.30, "U": 1.50, "w": 1.0},
        "full_divergence": {"angle": "reference", "dir": "min", "L": 0.50, "U": 3.00, "w": 1.0},
        # informational only (not in MPO) — see note above
        "bondlen_dev": {
            "angle": "geometry",
            "dir": "min",
            "L": 0.05,
            "U": 0.25,
            "w": 1.0,
            "informational": True,
        },
        "uff_displacement": {
            "angle": "relax",
            "dir": "min",
            "L": 0.20,
            "U": 1.50,
            "w": 1.0,
            "informational": True,
        },
        "uff_energy_per_atom": {
            "angle": "relax",
            "dir": "min",
            "L": 0.0,
            "U": 5.0,
            "w": 1.0,
            "informational": True,
        },
    },
}
# Which metrics are intrinsic (computable from one structure alone) vs comparative (need the
# matched input) -- drives which cohorts get which metrics.
INTRINSIC = {"clash_vdw", "clash_severe", "worst_overlap", "bondlen_dev", "angle_strain"}


def load_config(path):
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    if path:
        with open(path) as fh:
            override = json.load(fh)
        for k, v in override.get("angle_weights", {}).items():
            cfg["angle_weights"][k] = v
        for m, spec in override.get("metrics", {}).items():
            cfg["metrics"].setdefault(m, {}).update(spec)
    return cfg


def desirability(x, spec):
    """Map a raw metric value to a desirability in [0,1] per its direction and bounds."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    L, U, d = spec["L"], spec["U"], spec["dir"]
    if d == "min":
        if x <= L:
            return 1.0
        if x >= U:
            return 0.0
        return (U - x) / (U - L)
    if d == "max":
        if x >= U:
            return 1.0
        if x <= L:
            return 0.0
        return (x - L) / (U - L)
    raise ValueError(f"unknown direction {d!r}")


def aggregate(desir, cfg):
    """(mpo, angle_scores) from per-metric desirabilities.

    Within an angle: weighted arithmetic mean of available components. Across angles:
    weighted geometric mean (a fully-failed angle strongly penalizes), desirabilities
    floored at EPS, weights renormalized over the angles that have any data.
    """
    specs = cfg["metrics"]
    angle_scores = {}
    for angle in cfg["angle_weights"]:
        comps = [
            (specs[m]["w"], desir[m])
            for m in specs
            if specs[m]["angle"] == angle
            and not specs[m].get("informational")
            and desir.get(m) is not None
        ]
        if comps:
            wsum = sum(w for w, _ in comps)
            angle_scores[angle] = sum(w * v for w, v in comps) / wsum
    if not angle_scores:
        return None, {}
    num = den = 0.0
    for angle, score in angle_scores.items():
        w = cfg["angle_weights"][angle]
        num += w * math.log(max(score, EPS))
        den += w
    return math.exp(num / den), angle_scores


# ------------------------------------------------------------------ geometry primitives
def _rcov(s):
    try:
        return _PT.GetRcovalent(_PT.GetAtomicNumber(s))
    except Exception:
        return 0.7


def _rvdw(s):
    try:
        return _PT.GetRvdw(_PT.GetAtomicNumber(s))
    except Exception:
        return 1.7


def _adjacency(sym, xs, factor=1.3):
    """Distance/covalent-radius bond perception. Returns (bool adj, distance matrix)."""
    D = squareform(pdist(xs))
    np.fill_diagonal(D, 9e9)
    rc = np.array([_rcov(s) for s in sym])
    return D < factor * (rc[:, None] + rc[None, :]), D


def geometry_metrics(sym, xs, metal_symbol=None):
    """Intrinsic single-structure geometry metrics (+ perceived CN and bond count)."""
    n = len(sym)
    metal_idx = cf.find_metal_index(sym, xs, metal_symbol)
    adj, D = _adjacency(sym, xs)
    rv = np.array([_rvdw(s) for s in sym])
    overlap = D / (rv[:, None] + rv[None, :])
    # steric clashes: non-bonded, non-geminal (1-3) pairs inside vdW contact
    A = adj.astype(int)
    geminal = (A @ A) > 0
    nonbond = ~(adj | geminal | np.eye(n, dtype=bool))
    iu = np.triu_indices(n, 1)
    nb = nonbond[iu]
    ov = overlap[iu]
    clash_vdw = int(((ov < 0.75) & nb).sum())
    clash_severe = int(((ov < 0.60) & nb).sum())
    worst_overlap = float(ov[nb].min()) if nb.any() else 1.0
    # bond-length deviation (fraction): metal-ligand vs curated table, organic vs covalent sum
    devs = []
    bi, bj = np.where(np.triu(adj, 1))
    for i, j in zip(bi, bj):
        d = D[i, j]
        if i == metal_idx or j == metal_idx:
            lig = sym[j] if i == metal_idx else sym[i]
            ms = sym[metal_idx]
            ideal = bond_length(ms, lig) or (_rcov(ms) + _rcov(lig))
        else:
            ideal = _rcov(sym[i]) + _rcov(sym[j])
        if ideal > 0:
            devs.append(abs(d - ideal) / ideal)
    bondlen_dev = float(np.median(devs)) if devs else 0.0
    # bond-angle strain at light (non-metal) centers, ideal from #bonded neighbours
    strain = []
    for c in range(n):
        if c == metal_idx or sym[c] in _TM:
            continue
        nbrs = [k for k in range(n) if adj[c, k]]
        if len(nbrs) < 2:
            continue
        ideal = {2: 180.0, 3: 120.0, 4: 109.47, 5: 90.0, 6: 90.0}.get(len(nbrs), 109.47)
        for a in range(len(nbrs)):
            for b in range(a + 1, len(nbrs)):
                v1 = xs[nbrs[a]] - xs[c]
                v2 = xs[nbrs[b]] - xs[c]
                cs = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9)
                strain.append(abs(math.degrees(math.acos(max(-1.0, min(1.0, cs)))) - ideal))
    angle_strain = float(np.mean(strain)) if strain else 0.0
    perceived_cn = int(adj[metal_idx].sum()) if metal_idx is not None else 0
    n_bonds = int(adj[iu].sum())
    return {
        "clash_vdw": clash_vdw,
        "clash_severe": clash_severe,
        "worst_overlap": worst_overlap,
        "bondlen_dev": bondlen_dev,
        "angle_strain": angle_strain,
        "perceived_cn": perceived_cn,
        "n_bonds": n_bonds,
        "n_atoms": n,
    }


def full_divergence(in_sym, in_xs, gen_sym, gen_xs, metal_symbol, geo_code):
    """Symmetric element-aware heavy-atom overlay divergence (A) between input and generated.

    Aligns both to the geometry's canonical frame and snaps generated onto input (reusing
    compare_failures' alignment), then averages the two directional Chamfer distances.
    """
    in_m = cf.find_metal_index(in_sym, in_xs, metal_symbol)
    gen_m = cf.find_metal_index(gen_sym, gen_xs, metal_symbol)
    if in_m is None or gen_m is None:
        return None
    in_c = in_xs - in_xs[in_m]
    gen_c = gen_xs - gen_xs[gen_m]
    slots = len(TEMPLATES[geo_code]) if geo_code in TEMPLATES else None
    if slots:
        in_d = cf.closest_heavy_indices(in_sym, in_c, in_m, slots)
        gen_d = cf.closest_heavy_indices(gen_sym, gen_c, gen_m, slots)
    else:
        in_d = cf.detect_donors(in_sym, in_c, in_m)
        gen_d = cf.detect_donors(gen_sym, gen_c, gen_m)
    in_c = cf.canonical_rotation(in_c, in_sym, in_d, in_m, geo_code)[0].apply(in_c)
    gen_c = cf.canonical_rotation(gen_c, gen_sym, gen_d, gen_m, geo_code)[0].apply(gen_c)
    gen_c = cf.refine_onto_input(gen_c, gen_sym, in_c, in_sym).apply(gen_c)
    gh = np.array(gen_sym) != "H"
    ih = np.array(in_sym) != "H"
    if not gh.any() or not ih.any():
        return None
    g_xyz, g_s = gen_c[gh], np.asarray(gen_sym)[gh]
    i_xyz, i_s = in_c[ih], np.asarray(in_sym)[ih]
    return 0.5 * (
        cf._element_aware_chamfer(g_xyz, g_s, i_xyz, i_s)
        + cf._element_aware_chamfer(i_xyz, i_s, g_xyz, g_s)
    )


def uff_relax_proxy(sym, xs):
    """Ligands-only UFF relaxation: (displacement_rmsd, energy_drop_per_atom) or (None, None).

    Strips the metal (rdkit UFF has no TM params), perceives organic single bonds from
    coordinates, minimizes, and reports how far atoms moved and the per-atom energy drop.
    Best-effort: perception/sanitize/typing failures on distorted ligands return (None, None).
    """
    try:
        keep = [i for i, s in enumerate(sym) if s not in _TM]
        if len(keep) < 4:
            return None, None
        sub = xs[keep]
        subsym = [sym[i] for i in keep]
        D = squareform(pdist(sub))
        np.fill_diagonal(D, 9e9)
        rc = np.array([_rcov(s) for s in subsym])
        bond = D < 1.3 * (rc[:, None] + rc[None, :])
        rw = Chem.RWMol()
        for s in subsym:
            rw.AddAtom(Chem.Atom(s))
        bi, bj = np.where(np.triu(bond, 1))
        for a, b in zip(bi, bj):
            rw.AddBond(int(a), int(b), Chem.BondType.SINGLE)
        m = rw.GetMol()
        conf = Chem.Conformer(m.GetNumAtoms())
        for a in range(len(keep)):
            conf.SetAtomPosition(a, tuple(float(v) for v in sub[a]))
        m.AddConformer(conf)
        flags = (
            Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        )
        Chem.SanitizeMol(m, sanitizeOps=flags)
        ff = AllChem.UFFGetMoleculeForceField(m)
        if ff is None:
            return None, None
        before = m.GetConformer().GetPositions().copy()
        e0 = ff.CalcEnergy()
        ff.Minimize(maxIts=200)
        after = m.GetConformer().GetPositions()
        e1 = ff.CalcEnergy()
        disp = float(np.sqrt((((after - before) ** 2).sum(axis=1)).mean()))
        return disp, float((e0 - e1) / m.GetNumAtoms())
    except Exception:
        return None, None


# ----------------------------------------------------------------- per-molecule analysis
def analyze(report, results_dir, cfg, do_uff=True):
    """Full per-molecule row: raw metrics (gen + input), desirabilities, angle scores, MPO."""
    mol = report["molecule"]
    row = {"molecule": mol, "status": report.get("status")}
    input_path = report.get("input_xyz")
    gen_path = os.path.join(results_dir, "structures", f"{mol}_generated.xyz")
    if not input_path or not os.path.exists(input_path) or not os.path.exists(gen_path):
        row["error"] = "missing input or generated xyz"
        return row

    metal_sym, geo_code = cf.parse_metal_and_geo(report.get("smiles_1"))
    row["metal"] = metal_sym
    row["geom"] = geo_code
    row["failure_class"] = cf.classify_failure(report.get("error")) if report.get("error") else None

    in_sym, in_xs = cf.parse_xyz(input_path)
    gen_sym, gen_xs = cf.parse_xyz(gen_path)

    try:
        gm = geometry_metrics(gen_sym, gen_xs, metal_sym)
        im = geometry_metrics(in_sym, in_xs, metal_sym)
    except Exception as exc:
        row["error"] = f"geometry_metrics: {type(exc).__name__}: {exc}"
        return row
    row["gen"] = gm
    row["input"] = im

    # formal coordination number from the intended OIN (slot count)
    formal_cn = None
    try:
        parsed = OINParser().parse(report["smiles_1"])
        formal_cn = len({v.slot for v in parsed.vectors})
    except Exception:
        pass
    row["formal_cn"] = formal_cn

    # comparative (generated-vs-input) raw metrics
    raw = dict(gm)  # start from generated intrinsic values
    raw["cn_divergence"] = abs(gm["perceived_cn"] - im["perceived_cn"])
    raw["bondcount_divergence"] = abs(gm["n_bonds"] - im["n_bonds"])
    raw["perception_ok"] = 1.0 if report.get("smiles_2") else 0.0
    r = (report.get("metrics") or {}).get("rmsd")
    raw["coord_rmsd"] = r if (isinstance(r, (int, float)) and r < 900) else None
    try:
        raw["full_divergence"] = full_divergence(
            in_sym, in_xs, gen_sym, gen_xs, metal_sym, geo_code
        )
    except Exception:
        raw["full_divergence"] = None
    if do_uff:
        raw["uff_displacement"], raw["uff_energy_per_atom"] = uff_relax_proxy(gen_sym, gen_xs)
    else:
        raw["uff_displacement"] = raw["uff_energy_per_atom"] = None

    desir = {m: desirability(raw.get(m), cfg["metrics"][m]) for m in cfg["metrics"]}
    mpo, angle_scores = aggregate(desir, cfg)
    row["raw"] = {m: raw.get(m) for m in cfg["metrics"]}
    row["desirability"] = desir
    row["angle_scores"] = angle_scores
    row["mpo"] = mpo
    # input-cohort MPO uses intrinsic desirabilities only (calibration contrast)
    in_desir = {m: desirability(im.get(m), cfg["metrics"][m]) for m in INTRINSIC}
    row["input_geometry_score"] = aggregate({**{m: None for m in cfg["metrics"]}, **in_desir}, cfg)[
        1
    ].get("geometry")
    return row


# --------------------------------------------------------------------------- reporting
def _pct(vals, p):
    v = [x for x in vals if isinstance(x, (int, float))]
    return float(np.percentile(v, p)) if v else float("nan")


def _median(vals):
    return _pct(vals, 50)


def build_report_md(rows, cfg, results_dir):
    """Compose the human-readable markdown distortion report."""
    scored = [r for r in rows if r.get("mpo") is not None]
    n = len(scored)
    L = []
    L.append("# MetalloGen structure-distortion report\n")
    L.append(
        f"Population: `{os.path.basename(results_dir)}` — {n} generated structures scored "
        f"(of {len(rows)} attempted).\n"
    )

    # headline metric table: generated vs input
    L.append("## Headline metrics (generated vs. real input)\n")
    L.append("| Metric | Generated median | Generated p95 | Input median | Input p95 |")
    L.append("|---|---|---|---|---|")
    label = {
        "clash_vdw": "vdW clashes (<0.75)",
        "clash_severe": "severe clashes (<0.60)",
        "worst_overlap": "worst overlap ratio",
        "bondlen_dev": "bond-length dev (frac)",
        "angle_strain": "angle strain (deg)",
    }
    for m in ["clash_vdw", "clash_severe", "worst_overlap", "bondlen_dev", "angle_strain"]:
        g = [r["gen"][m] for r in scored if "gen" in r]
        i = [r["input"][m] for r in scored if "input" in r]
        L.append(
            f"| {label[m]} | {_median(g):.2f} | {_pct(g, 95):.2f} | "
            f"{_median(i):.2f} | {_pct(i, 95):.2f} |"
        )
    clashing = 100 * np.mean([r["gen"]["clash_vdw"] > 0 for r in scored])
    in_clashing = 100 * np.mean([r["input"]["clash_vdw"] > 0 for r in scored])
    L.append(
        f"\n**{clashing:.0f}%** of generated structures have ≥1 vdW steric clash "
        f"(vs **{in_clashing:.0f}%** of inputs).\n"
    )

    # MPO + angle-score distribution
    L.append("## MPO quality score\n")
    L.append("| Score | median | p25 | p10 |")
    L.append("|---|---|---|---|")
    mpos = [r["mpo"] for r in scored]
    L.append(
        f"| MPO (overall) | {_median(mpos):.3f} | {_pct(mpos, 25):.3f} | {_pct(mpos, 10):.3f} |"
    )
    for angle in cfg["angle_weights"]:
        a = [
            r["angle_scores"].get(angle) for r in scored if r["angle_scores"].get(angle) is not None
        ]
        L.append(
            f"| {angle} sub-score | {_median(a):.3f} | {_pct(a, 25):.3f} | {_pct(a, 10):.3f} |"
        )
    L.append(
        "\nAngle weights (MPO): "
        + ", ".join(f"{a} {w:g}" for a, w in cfg["angle_weights"].items())
        + ".\n"
    )

    # informational metrics — computed but excluded from the MPO (do not discriminate)
    L.append("## Informational metrics (excluded from MPO)\n")
    L.append("These were computed but proved **not** to discriminate MetalloGen distortion.\n")
    gbd = [r["gen"]["bondlen_dev"] for r in scored]
    ibd = [r["input"]["bondlen_dev"] for r in scored]
    ud = [
        r["raw"].get("uff_displacement")
        for r in scored
        if r["raw"].get("uff_displacement") is not None
    ]
    na = 100 * np.mean([r["raw"].get("uff_displacement") is None for r in scored])
    L.append("| Metric | Generated median | Input median | Verdict |")
    L.append("|---|---|---|---|")
    L.append(
        f"| bond-length dev (frac) | {_median(gbd):.3f} | {_median(ibd):.3f} | "
        "inverted — generator places near-ideal bonds; bonds are not the distortion |"
    )
    uff_cell = f"{_median(ud):.2f}" if ud else "(not run)"
    uff_note = (
        f"confounded (all-single-bond FF); N/A rate {na:.0f}%; needs xtb/MACE to be valid"
        if ud
        else "skipped (--no-uff)"
    )
    L.append(f"| ligands-only UFF displacement (Å) | {uff_cell} | — | {uff_note} |")
    L.append("")

    # stratified MPO
    for key, title in (("geom", "geometry"), ("metal", "metal")):
        L.append(f"## MPO by {title}\n")
        L.append(f"| {title} | n | median MPO | % clashing | mean clashes |")
        L.append("|---|---|---|---|---|")
        groups = {}
        for r in scored:
            groups.setdefault(r.get(key), []).append(r)
        stats = []
        for g, rs in groups.items():
            if len(rs) < 8:
                continue
            stats.append(
                (
                    g,
                    len(rs),
                    _median([x["mpo"] for x in rs]),
                    100 * np.mean([x["gen"]["clash_vdw"] > 0 for x in rs]),
                    np.mean([x["gen"]["clash_vdw"] for x in rs]),
                )
            )
        for g, ng, mm, pc, mc in sorted(stats, key=lambda x: x[2]):
            L.append(f"| {g} | {ng} | {mm:.3f} | {pc:.0f}% | {mc:.2f} |")
        L.append("")

    # worst offenders (v0.4.3 priority queue)
    L.append("## Worst 25 by MPO (v0.4.3 priority queue)\n")
    L.append("| molecule | metal | geom | MPO | clashes | angle° | full-div Å | coord-RMSD |")
    L.append("|---|---|---|---|---|---|---|---|")
    for r in sorted(scored, key=lambda x: x["mpo"])[:25]:
        fd = r["raw"].get("full_divergence")
        cr = r["raw"].get("coord_rmsd")
        L.append(
            f"| {r['molecule']} | {r.get('metal')} | {r.get('geom')} | {r['mpo']:.3f} | "
            f"{r['gen']['clash_vdw']} | {r['gen']['angle_strain']:.1f} | "
            f"{fd:.2f} | {cr if cr is None else round(cr, 2)} |"
        )

    L.append("\n## Reproducing\n")
    L.append("```")
    L.append(
        f"uv run python tools/structure_distortion_report.py --results-dir {results_dir} "
        "--output-dir <out>"
    )
    L.append("```")
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------------- CLI
def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--results-dir",
        default="../tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042",
        help="Sweep results dir (individual_reports/ + structures/). Default: capstone non-quick.",
    )
    parser.add_argument(
        "--output-dir",
        default="../results-v043-distortion",
        help="Where to write distortion_metrics.json + distortion_report.md.",
    )
    parser.add_argument("--config", default=None, help="JSON overriding metric bounds / weights.")
    parser.add_argument("--only", default=None, help="Comma-separated molecule names.")
    parser.add_argument(
        "--sample",
        type=int,
        default=None,
        help="Score N structures spread across (geometry, failure-class).",
    )
    parser.add_argument("--limit", type=int, default=None, help="Cap number of structures.")
    parser.add_argument("--shard", default=None, help="I:N deterministic slice for parallel runs.")
    parser.add_argument(
        "--no-uff", action="store_true", help="Skip the ligands-only UFF relaxation proxy (faster)."
    )
    args = parser.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    cfg = load_config(args.config)

    reports = []
    for path in sorted(glob.glob(os.path.join(results_dir, "individual_reports", "*.json"))):
        with open(path) as fh:
            rep = json.load(fh)
        if rep.get("status") in ("success", "failed") and os.path.exists(
            os.path.join(results_dir, "structures", f"{rep['molecule']}_generated.xyz")
        ):
            reports.append(rep)
    print(f"{len(reports)} structures with saved geometry to score.")

    if args.only:
        wanted = {n.strip().removesuffix(".xyz") for n in args.only.split(",") if n.strip()}
        reports = [r for r in reports if r["molecule"] in wanted]
        print(f"Only mode: matched {len(reports)} of {len(wanted)}.")
    elif args.sample:
        reports = cf.diverse_sample(reports, args.sample)
        print(f"Sample mode: selected {len(reports)}.")
    if args.shard:
        try:
            i, nsh = (int(x) for x in args.shard.split(":"))
        except ValueError:
            parser.error(f"--shard must be I:N (got {args.shard!r})")
        if not (1 <= i <= nsh):
            parser.error(f"--shard index out of range: {args.shard!r}")
        reports = reports[i - 1 :: nsh]
        print(f"Shard {i}/{nsh}: {len(reports)} structures.")
    if args.limit:
        reports = reports[: args.limit]

    rows, errors = [], []
    for k, rep in enumerate(reports):
        row = analyze(rep, results_dir, cfg, do_uff=not args.no_uff)
        if row.get("mpo") is None:
            errors.append({"molecule": row["molecule"], "error": row.get("error", "unscored")})
        else:
            rows.append(row)
        if (k + 1) % 500 == 0:
            print(f"  scored {k + 1}/{len(reports)} ...")

    with open(os.path.join(output_dir, "distortion_metrics.json"), "w") as fh:
        json.dump({"config": cfg, "rows": rows, "errors": errors}, fh, indent=1)
    with open(os.path.join(output_dir, "distortion_report.md"), "w") as fh:
        fh.write(build_report_md(rows, cfg, results_dir))

    if rows:
        mpos = [r["mpo"] for r in rows]
        print(
            f"Scored {len(rows)} ({len(errors)} unscored). Median MPO {np.median(mpos):.3f}. "
            f"Report + JSON in {output_dir}"
        )


if __name__ == "__main__":
    main()
