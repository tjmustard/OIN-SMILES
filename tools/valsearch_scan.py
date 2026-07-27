#!/usr/bin/env python
"""Pre-scan: how many ligand perceptions exceed ``_VALENCE_COMBO_CAP``, for free.

Q1's first half. Deliberately does **not** call ``AC2BO``: the over-cap decision is a
function of the adjacency matrix alone, so ``get_basic_mol`` + ``MetalDisconnector`` +
``GetAdjacencyMatrix`` + ``possible_valences`` answers it in milliseconds per molecule.
That is what makes it affordable to ask the rate question at all while the release sweep
is saturating the host.

Reports per molecule: number of ligand fragments, how many are over cap, and the
combo size of the worst one. Molecules with an over-cap fragment are the candidate
population for Q2 (``tools/valsearch_budget_ab.py``).

    $V tools/valsearch_scan.py --dataset <dir> --n 100 --seed 1 --out scan.json
"""

import argparse
import json
import random
import re
import sys
import time
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdmolops
from rdkit.Chem.MolStandardize import rdMolStandardize

from oinsmiles.utils.perception_core import (
    _VALENCE_COMBO_CAP,
    possible_valences,
    valence_combo_size,
)
from oinsmiles.utils.perception_tmc import (
    TRANSITION_METALS_NUM,
    MetalNon_Hg,
    get_basic_mol,
    params,
)

RDLogger.DisableLog("rdApp.*")

_CHARGE_RE = re.compile(r"Charge:\s*(-?\d+)")


def read_charge(path):
    with open(path) as fh:
        fh.readline()
        comment = fh.readline()
    match = _CHARGE_RE.search(comment)
    return int(match.group(1)) if match else 0


def scan_one(path):
    """Per-ligand-fragment over-cap report for one xyz, without running the search."""
    charge = read_charge(path)
    mol, _ = get_basic_mol(str(path), charge)
    for a in mol.GetAtoms():
        a.SetIntProp("__origIdx", a.GetIdx())

    mdis = rdMolStandardize.MetalDisconnector(params)
    mdis.SetMetalNon(Chem.MolFromSmarts(MetalNon_Hg))
    frags = mdis.Disconnect(mol)
    frag_mols = rdmolops.GetMolFrags(frags, asMols=True)

    fragments = []
    for i, frag in enumerate(frag_mols):
        atoms = [a.GetAtomicNum() for a in frag.GetAtoms()]
        if any(n in TRANSITION_METALS_NUM for n in atoms):
            continue
        AC = Chem.rdmolops.GetAdjacencyMatrix(frag)
        AC_valence = list(AC.sum(axis=1))
        # The ladder's first arm uses allow_carbenes=True; the no-carbene arm can only
        # shrink a valence list, so True is the arm that decides "over cap".
        vll = possible_valences(AC_valence, atoms, allow_carbenes=True)
        size = valence_combo_size(vll)
        fragments.append(
            {
                "frag": i,
                "n_atoms": frag.GetNumAtoms(),
                "combo_size": size,
                "over_cap": size > _VALENCE_COMBO_CAP,
            }
        )
    return {
        "n_atoms": mol.GetNumAtoms(),
        "n_lig_frags": len(fragments),
        "n_over_cap": sum(1 for f in fragments if f["over_cap"]),
        "max_combo_size": max((f["combo_size"] for f in fragments), default=0),
        "fragments": fragments,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, help="dir containing cat/ and photo/")
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None)
    ap.add_argument(
        "--mols",
        default=None,
        help="comma-separated refcodes or @file, instead of a random sample. Use this to "
        "scan a cohort stratified by encode time rather than the corpus at large.",
    )
    args = ap.parse_args()

    root = Path(args.dataset)
    files = sorted(root.glob("cat/*.xyz")) + sorted(root.glob("photo/*.xyz"))
    if not files:
        # A worktree does not have the gitignored dataset; fail loudly rather than
        # reporting a serene 0/0 (which cost a full A/B run once already).
        sys.exit(f"FATAL: no .xyz under {root}/cat or {root}/photo -- pass --dataset")
    print(f"corpus: {len(files)} xyz files under {root}", flush=True)

    if args.mols:
        raw = (
            Path(args.mols[1:]).read_text().split()
            if args.mols.startswith("@")
            else args.mols.split(",")
        )
        wanted = [x.strip() for x in raw if x.strip()]
        by_stem = {p.stem: p for p in files}
        sample, missing = [], []
        for name in wanted:
            hit = by_stem.get(name) or by_stem.get(f"{name}_comp_0")
            (sample.append(hit) if hit else missing.append(name))
        if missing:
            print(f"NOT FOUND ({len(missing)}): {', '.join(sorted(missing)[:12])}", flush=True)
        if not sample:
            sys.exit("FATAL: none of the requested refcodes exist in the corpus")
    else:
        rng = random.Random(args.seed)
        sample = rng.sample(files, min(args.n, len(files)))

    results = {}
    errors = {}
    t0 = time.time()
    for idx, path in enumerate(sample, 1):
        # Key by "<subdir>/<stem>", NOT by stem: cat/ and photo/ share refcodes, and a
        # stem-keyed dict silently collapsed 8 of a 2000-molecule sample into 1992,
        # quietly shrinking the denominator of the rate this tool exists to report.
        name = f"{path.parent.name}/{path.stem}"
        try:
            results[name] = scan_one(path)
        except Exception as exc:  # noqa: BLE001 - a scan must not stop on one bad file
            errors[name] = f"{type(exc).__name__}: {exc}"
        if idx % 25 == 0:
            print(f"  {idx}/{len(sample)} ({time.time() - t0:.1f}s)", flush=True)
    collapsed = len(sample) - len(results) - len(errors)
    if collapsed:
        print(f"WARNING: {collapsed} sampled paths collapsed on a duplicate key", flush=True)

    over = {k: v for k, v in results.items() if v["n_over_cap"]}
    total_frags = sum(v["n_lig_frags"] for v in results.values())
    over_frags = sum(v["n_over_cap"] for v in results.values())

    print()
    print(f"sample                       : {len(sample)} molecules (seed {args.seed})")
    print(f"scanned OK                   : {len(results)}   errors: {len(errors)}")
    print(f"ligand fragments             : {total_frags}")
    print(
        f"  over _VALENCE_COMBO_CAP    : {over_frags}"
        f"  ({100.0 * over_frags / total_frags if total_frags else 0:.1f}%)"
    )
    print(
        f"molecules with >=1 over cap  : {len(over)}"
        f"  ({100.0 * len(over) / len(results) if results else 0:.1f}%)"
    )
    print(f"scan wall                    : {time.time() - t0:.1f}s")
    if over:
        print("\nover-cap molecules (ligand frag atoms):")
        for name in sorted(over):
            frs = [f for f in over[name]["fragments"] if f["over_cap"]]
            print(f"  {name:26s} {[f['n_atoms'] for f in frs]}")
    if errors:
        print("\nerrors:")
        for name, msg in sorted(errors.items())[:10]:
            print(f"  {name:26s} {msg}")

    if args.out:
        Path(args.out).write_text(
            json.dumps(
                {
                    "seed": args.seed,
                    "n_requested": args.n,
                    "sample": [p.stem for p in sample],
                    "results": results,
                    "errors": errors,
                },
                indent=1,
            )
        )
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
