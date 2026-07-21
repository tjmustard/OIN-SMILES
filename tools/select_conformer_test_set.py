#!/usr/bin/env python3
"""Curate a size-stratified conformer-invariance test set from the passing sweep.

Reads the round-trip sweep summaries (``results-v0.4.0`` quick pass, optionally
intersected with the ``results-capstone-v042`` accuracy-clean pass), keeps only
molecules that round-trip cleanly, and samples ~N of them so the set mirrors the
heavy-atom-count distribution of the full passing pool. The two Pt fixtures
(CisPlatin / TransPlatin) are always included as small, fast anchors. Selected
structures are copied into ``tests/fixtures/conformer_set/`` (the tmCAT-tmPHOTO
dataset itself is gitignored, so the subset must be committed) alongside a
``manifest.json`` and a human-readable ``README.md``.

Selection is deterministic given ``--seed``: no wall-clock or unseeded randomness.

Example (run from anywhere; point --dataset-dir at the gitignored dataset in the
main checkout, since a fresh worktree will not contain it)::

    python tools/select_conformer_test_set.py --n 30 --seed 42 \
        --dataset-dir /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FIXTURES = REPO / "tests" / "fixtures"
DEFAULT_OUT = FIXTURES / "conformer_set"

# The two mandatory anchors (already tracked under tests/fixtures/).
FIXTURE_ANCHORS = [
    ("CisPlatin", FIXTURES / "CisPlatin.xyz"),
    ("TransPlatin", FIXTURES / "TransPlatin.xyz"),
]

# Heavy-atom strata (upper edge inclusive). The final bin is capped by --max-heavy.
BINS = [
    ("<=20", 0, 20),
    ("21-30", 21, 30),
    ("31-40", 31, 40),
    ("41-50", 41, 50),
    ("51-75", 51, 75),
    ("76+", 76, 10_000),
]

# Standard atomic weights for the elements that occur in tmCAT-tmPHOTO (organics
# + the 29 transition/main-group metals). Kept as a local table so the selector
# stays pure-stdlib (no rdkit import needed just to sum masses).
ATOMIC_WEIGHT = {
    "H": 1.008,
    "He": 4.003,
    "Li": 6.941,
    "Be": 9.012,
    "B": 10.811,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Ne": 20.180,
    "Na": 22.990,
    "Mg": 24.305,
    "Al": 26.982,
    "Si": 28.086,
    "P": 30.974,
    "S": 32.065,
    "Cl": 35.453,
    "Ar": 39.948,
    "K": 39.098,
    "Ca": 40.078,
    "Sc": 44.956,
    "Ti": 47.867,
    "V": 50.942,
    "Cr": 51.996,
    "Mn": 54.938,
    "Fe": 55.845,
    "Co": 58.933,
    "Ni": 58.693,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.640,
    "As": 74.922,
    "Se": 78.960,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.468,
    "Sr": 87.620,
    "Y": 88.906,
    "Zr": 91.224,
    "Nb": 92.906,
    "Mo": 95.960,
    "Tc": 98.000,
    "Ru": 101.07,
    "Rh": 102.906,
    "Pd": 106.42,
    "Ag": 107.868,
    "Cd": 112.411,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.600,
    "I": 126.904,
    "Xe": 131.293,
    "Cs": 132.905,
    "Ba": 137.327,
    "La": 138.905,
    "Ce": 140.116,
    "Hf": 178.490,
    "Ta": 180.948,
    "W": 183.840,
    "Re": 186.207,
    "Os": 190.230,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.967,
    "Hg": 200.590,
    "Tl": 204.383,
    "Pb": 207.200,
    "Bi": 208.980,
}

_METAL_RE = re.compile(r"^\[([A-Z][a-z]?)")
_CHARGE_RE = re.compile(r"Charge:\s*(-?\d+)")
_MULT_RE = re.compile(r"Multiplicity:\s*(\d+)")
_SPINS_RE = re.compile(r"Potential Spins:\s*\[([^\]]*)\]")


def parse_metal(smiles: str | None) -> str | None:
    if not smiles:
        return None
    m = _METAL_RE.match(smiles)
    return m.group(1) if m else None


def parse_xyz(path: Path):
    """Return (n_total, n_heavy, mass, charge, mult) from an .xyz file.

    Handles both dataset headers (``... | Charge: 0 | Potential Spins: [1, 3, 5]``)
    and fixture headers (``... Charge: 0, Multiplicity: 1``).
    """
    lines = path.read_text().splitlines()
    n_total = int(lines[0].strip())
    header = lines[1] if len(lines) > 1 else ""
    n_heavy = 0
    mass = 0.0
    unknown = set()
    for ln in lines[2 : 2 + n_total]:
        parts = ln.split()
        if len(parts) < 4:
            continue
        el = parts[0]
        if el != "H":
            n_heavy += 1
        if el in ATOMIC_WEIGHT:
            mass += ATOMIC_WEIGHT[el]
        else:
            unknown.add(el)
    if unknown:
        print(f"  WARN: {path.name}: no mass for {sorted(unknown)}", file=sys.stderr)
    charge = 0
    cm = _CHARGE_RE.search(header)
    if cm:
        charge = int(cm.group(1))
    mult = 1
    mm = _MULT_RE.search(header)
    if mm:
        mult = int(mm.group(1))
    else:
        sm = _SPINS_RE.search(header)
        if sm:
            spins = [int(x) for x in re.findall(r"-?\d+", sm.group(1))]
            if spins:
                mult = min(spins)  # lowest multiplicity = default ground state guess
    return n_total, n_heavy, round(mass, 3), charge, mult


def bin_for(heavy: int) -> str:
    for name, lo, hi in BINS:
        if lo <= heavy <= hi:
            return name
    return BINS[-1][0]


def load_success(summary_path: Path) -> dict:
    """Return {molecule: record} for status==success rows of a sweep summary."""
    with summary_path.open() as f:
        recs = json.load(f)
    return {r["molecule"]: r for r in recs if r.get("status") == "success"}


def largest_remainder(fractions: dict[str, float], total: int) -> dict[str, int]:
    raw = {k: v * total for k, v in fractions.items()}
    base = {k: int(v) for k, v in raw.items()}
    used = sum(base.values())
    remainder = total - used
    order = sorted(fractions, key=lambda k: (-(raw[k] - base[k]), k))
    for k in order[:remainder]:
        base[k] += 1
    return base


def allocate(bin_counts: dict[str, int], n_dataset: int) -> dict[str, int]:
    """Proportional allocation of n_dataset picks across non-empty bins, floor 1."""
    present = {b: c for b, c in bin_counts.items() if c > 0}
    if not present:
        return {}
    total = sum(present.values())
    fractions = {b: c / total for b, c in present.items()}
    alloc = largest_remainder(fractions, n_dataset)
    # Enforce a floor of 1 per present bin, then trim from the largest bins to
    # keep the sum == n_dataset. Never exceed a bin's available candidate count.
    for b in present:
        alloc[b] = max(1, alloc.get(b, 0))
        alloc[b] = min(alloc[b], present[b])

    # Rebalance to hit the exact target.
    def _resize(target):
        while sum(alloc.values()) > target:
            b = max(alloc, key=lambda k: (alloc[k], present[k]))
            if alloc[b] <= 1:
                break
            alloc[b] -= 1
        while sum(alloc.values()) < target:
            candidates = [b for b in present if alloc[b] < present[b]]
            if not candidates:
                break
            b = max(candidates, key=lambda k: (present[k] - alloc[k], present[k]))
            alloc[b] += 1

    _resize(n_dataset)
    return alloc


def quality_key(rec: dict, clean: set, seed: int):
    """Deterministic ordering: accuracy-clean first, then fast, then small."""
    mol = rec["molecule"]
    elapsed = (rec.get("metrics") or {}).get("elapsed_s", 1e9)
    tie = hashlib.sha1(f"{seed}:{mol}".encode()).hexdigest()
    return (0 if mol in clean else 1, elapsed, rec.get("_total", 1e9), tie)


def pick_bin(cands: list, quota: int, metal_count: Counter, cap: int):
    """Greedily pick `quota` candidates preferring rarer-so-far metals (<= cap)."""
    remaining = list(cands)
    chosen = []
    while len(chosen) < quota and remaining:
        avail = [c for c in remaining if metal_count[c["metal"]] < cap]
        pool = avail if avail else remaining  # relax the cap if nothing fits
        best = min(pool, key=lambda c: (metal_count[c["metal"]], c["_q"]))
        chosen.append(best)
        metal_count[best["metal"]] += 1
        remaining.remove(best)
    return chosen


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--dataset-dir",
        type=Path,
        default=REPO / "tmCAT-tmPHOTO_xyz_dataset",
        help="Path to the (gitignored) tmCAT-tmPHOTO dataset in the main checkout.",
    )
    ap.add_argument("--v040-summary", type=Path, default=None)
    ap.add_argument("--capstone-summary", type=Path, default=None)
    ap.add_argument(
        "--pool",
        choices=["both", "v040"],
        default="both",
        help="'both' intersects the quick pass with the accuracy-clean pass.",
    )
    ap.add_argument("--n", type=int, default=30, help="Target total set size (incl. 2 fixtures).")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument(
        "--max-heavy", type=int, default=90, help="Exclude candidates above this many heavy atoms."
    )
    ap.add_argument("--max-per-metal", type=int, default=3)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    ds = args.dataset_dir
    v040 = args.v040_summary or ds / "results-v0.4.0" / "summary_roundtrip.json"
    capstone = args.capstone_summary or ds / "results-capstone-v042" / "summary_roundtrip.json"
    if not v040.exists():
        print(
            f"ERROR: sweep summary not found: {v040}\n"
            f"Point --dataset-dir at the gitignored dataset in the main checkout.",
            file=sys.stderr,
        )
        return 2

    success = load_success(v040)
    clean = set()
    if args.pool == "both" and capstone.exists():
        clean = set(load_success(capstone))
        pool = {m: r for m, r in success.items() if m in clean}
    else:
        if args.pool == "both":
            print(f"WARN: {capstone} missing; falling back to v040-only pool.", file=sys.stderr)
        pool = success
    print(f"pool: {len(pool)} molecules ({args.pool}); accuracy-clean set: {len(clean)}")

    # Build candidates: one component per refcode (best quality), heavy <= max_heavy.
    by_refcode: dict[str, dict] = {}
    fixture_names = {n for n, _ in FIXTURE_ANCHORS}
    for mol, rec in pool.items():
        if mol in fixture_names:
            continue
        xyz = Path(rec["input_xyz"])
        if not xyz.exists():
            # Remap onto --dataset-dir by subset/basename if the recorded path moved.
            subset = xyz.parent.name
            alt = ds / subset / xyz.name
            if alt.exists():
                xyz = alt
            else:
                continue
        try:
            n_total, n_heavy, mass, charge, mult = parse_xyz(xyz)
        except Exception as exc:  # noqa: BLE001 - skip unparseable, keep going
            print(f"  skip {mol}: {exc}", file=sys.stderr)
            continue
        if n_heavy > args.max_heavy:
            continue
        rec = dict(rec)
        rec["_total"] = n_total
        cand = {
            "molecule": mol,
            "src": str(xyz),
            "subset": xyz.parent.name,
            "metal": parse_metal(rec.get("smiles_1")) or "?",
            "n_atoms_total": n_total,
            "n_heavy": n_heavy,
            "mass": mass,
            "charge": charge,
            "mult": mult,
            "bin": bin_for(n_heavy),
            "elapsed_s": (rec.get("metrics") or {}).get("elapsed_s"),
            "source": "both" if mol in clean else "v040",
            "_q": quality_key(rec, clean, args.seed),
        }
        refcode = mol.split("_comp_")[0]
        prev = by_refcode.get(refcode)
        if prev is None or cand["_q"] < prev["_q"]:
            by_refcode[refcode] = cand
    candidates = list(by_refcode.values())
    print(f"candidates after dedupe/cap: {len(candidates)}")

    # Distribution of the eligible pool (for proportional allocation + the report).
    bin_counts = Counter(c["bin"] for c in candidates)
    n_dataset = max(0, args.n - len(FIXTURE_ANCHORS))
    alloc = allocate(dict(bin_counts), n_dataset)
    print("allocation (dataset picks per bin):", dict(alloc))

    # Metal budget starts with the fixtures (both Pt).
    metal_count: Counter = Counter()
    selected = []
    for name, path in FIXTURE_ANCHORS:
        n_total, n_heavy, mass, charge, mult = parse_xyz(path)
        metal = parse_metal(None) or "Pt"  # anchors are cisplatin/transplatin
        metal_count[metal] += 1
        selected.append(
            {
                "molecule": name,
                "path": str(path.relative_to(REPO)),
                "subset": "fixture",
                "metal": metal,
                "n_atoms_total": n_total,
                "n_heavy": n_heavy,
                "mass": mass,
                "charge": charge,
                "mult": mult,
                "bin": bin_for(n_heavy),
                "elapsed_s": None,
                "source": "fixture",
            }
        )

    by_bin: dict[str, list] = defaultdict(list)
    for c in candidates:
        by_bin[c["bin"]].append(c)

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, _, _ in BINS:
        quota = alloc.get(name, 0)
        if quota <= 0:
            continue
        cands = sorted(by_bin.get(name, []), key=lambda c: c["_q"])
        chosen = pick_bin(cands, quota, metal_count, args.max_per_metal)
        for c in chosen:
            dst = out_dir / f"{c['molecule']}.xyz"
            shutil.copy(c["src"], dst)
            selected.append(
                {
                    "molecule": c["molecule"],
                    "path": str(dst.relative_to(REPO)),
                    "subset": c["subset"],
                    "metal": c["metal"],
                    "n_atoms_total": c["n_atoms_total"],
                    "n_heavy": c["n_heavy"],
                    "mass": c["mass"],
                    "charge": c["charge"],
                    "mult": c["mult"],
                    "bin": c["bin"],
                    "elapsed_s": c["elapsed_s"],
                    "source": c["source"],
                }
            )

    selected.sort(key=lambda s: (s["n_heavy"], s["molecule"]))
    cmd = (
        f"python tools/select_conformer_test_set.py --n {args.n} --seed {args.seed} "
        f"--pool {args.pool} --max-heavy {args.max_heavy} --max-per-metal {args.max_per_metal}"
    )
    manifest = {
        "generated_with": cmd,
        "seed": args.seed,
        "pool": args.pool,
        "n_requested": args.n,
        "n_selected": len(selected),
        "structures": selected,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    write_readme(out_dir, selected, bin_counts, alloc, cmd, args)
    print(f"\nSelected {len(selected)} structures -> {out_dir}")
    achieved = Counter(s["bin"] for s in selected)
    print("achieved bins:", dict(achieved))
    print("metals:", dict(Counter(s["metal"] for s in selected)))
    return 0


def write_readme(out_dir, selected, bin_counts, alloc, cmd, args):
    total_pool = sum(bin_counts.values()) or 1
    achieved = Counter(s["bin"] for s in selected)
    metals = Counter(s["metal"] for s in selected)
    lines = [
        "# Conformer-invariance test set",
        "",
        "Size-stratified sample of transition-metal complexes that round-trip cleanly",
        "(XYZ -> OIN -> XYZ) in the sweep, used to verify that multiple conformers of a",
        "structure collapse to the same canonical OIN-SMILES. The two Pt fixtures",
        "(CisPlatin/TransPlatin) are always included as small, fast anchors and live in",
        "`tests/fixtures/` (referenced in place, not copied here).",
        "",
        "Regenerate deterministically with (point `--dataset-dir` at the gitignored",
        "dataset in the main checkout):",
        "",
        "```",
        cmd,
        "```",
        "",
        f"- pool: `{args.pool}` (v0.4.0 quick success"
        + (" ∩ capstone-v042 accuracy-clean)" if args.pool == "both" else ")"),
        f"- seed: `{args.seed}`  ·  max heavy atoms: `{args.max_heavy}`  ·  max per metal: `{args.max_per_metal}`",
        f"- selected: **{len(selected)}** structures",
        "",
        "## Heavy-atom strata (pool share vs. selected)",
        "",
        "| bin | pool share | selected |",
        "|---|---|---|",
    ]
    for name, _, _ in BINS:
        share = 100 * bin_counts.get(name, 0) / total_pool
        lines.append(f"| {name} | {share:.1f}% | {achieved.get(name, 0)} |")
    lines += [
        "",
        f"## Metal coverage ({len(metals)} distinct)",
        "",
        ", ".join(f"{m} × {c}" for m, c in metals.most_common()),
        "",
        "## Structures",
        "",
        "| molecule | metal | heavy | total | mass | bin | charge | mult | source |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for s in selected:
        lines.append(
            f"| {s['molecule']} | {s['metal']} | {s['n_heavy']} | {s['n_atoms_total']} | "
            f"{s['mass']:.1f} | {s['bin']} | {s['charge']} | {s['mult']} | {s['source']} |"
        )
    (out_dir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
