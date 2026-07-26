#!/usr/bin/env python
"""How many corpus molecules does the AC pruning loop silently corrupt?

The 34 `boron_cluster` molecules are the ones where amputating the cage happens to
produce something `get_lig_mol` cannot perceive at all. `VEJXOZ_comp_0` shows the
other outcome: the pruning loop deleted 6 of its 12 B-B cage bonds (50%) and the
encoder then invented a **C=B double bond** to balance the valences, producing an
OIN that is *self-consistent and round-trips* while describing the wrong molecule.
That failure mode is invisible to `canonical_roundtrip_key`, because the corrupted
encode is compared against its own corrupted mol.

So the blast radius is not 34. This scans the corpus for the deltahedral motif and
reports how many molecules have cage bonds deleted, split by whether they are in
the known encode_fail 34 or are currently counted as *passing*.

Cheap by construction: a text scan filters to files with >=3 boron first, and the
adjacency matrix is only computed for those. No encoding, no generation.

Usage:
    PYTHONPATH=src python tools/boron_blast_radius.py --dataset-dir <abs> \
        [--reports-dir <abs>/results-capstone-v042/individual_reports]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from collections import Counter

import numpy as np


def boron_count_fast(path: str) -> int:
    """Count 'B ' element tokens without building a mol."""
    n = 0
    try:
        with open(path) as fh:
            for i, line in enumerate(fh):
                if i < 2:
                    continue
                tok = line.split()
                if tok and tok[0] == "B":
                    n += 1
    except Exception:
        return 0
    return n


def analyse(path: str, tolerance: float = 0.5) -> dict | None:
    from rdkit import Chem

    from oinsmiles.utils.xyz2mol_local import boron_cage_vertices, read_xyz_file, xyz2AC_obabel

    atoms, _charge, xyz = read_xyz_file(path)
    pt = Chem.GetPeriodicTable()
    coords = np.asarray(xyz, dtype=float)
    n = len(atoms)
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    rcov = np.array([pt.GetRcovalent(int(z)) for z in atoms])
    AC_raw = ((d <= rcov[:, None] + rcov[None, :] + tolerance) & ~np.eye(n, dtype=bool)).astype(int)

    cage = boron_cage_vertices(list(atoms), AC_raw)
    if not cage:
        return None

    AC_pruned = np.asarray(xyz2AC_obabel(atoms, xyz, tolerance=tolerance)[0])
    b = [i for i, z in enumerate(atoms) if z == 5]
    bb_raw = sum(int(AC_raw[i][j]) for k, i in enumerate(b) for j in b[k + 1 :])
    bb_pruned = sum(int(AC_pruned[i][j]) for k, i in enumerate(b) for j in b[k + 1 :])
    return {
        "nB": len(b),
        "n_cage_vertices": len(cage),
        "BB_raw": bb_raw,
        "BB_pruned": bb_pruned,
        "BB_deleted": bb_raw - bb_pruned,
    }


def main() -> None:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from boron_characterize import BORON_COHORT

    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset-dir", required=True)
    ap.add_argument("--reports-dir", default=None)
    ap.add_argument("--out", default="tools/boron_blast_radius.json")
    args = ap.parse_args()

    passing = set()
    if args.reports_dir:
        for rp in glob.glob(os.path.join(args.reports_dir, "*.json")):
            try:
                with open(rp) as fh:
                    dd = json.load(fh)
            except Exception:
                continue
            if dd.get("smiles_1"):
                passing.add(
                    (dd.get("molecule") or os.path.basename(rp).replace(".json", "")).replace(
                        "_comp_0", ""
                    )
                )
        print(f"{len(passing)} molecules with a stored (non-null) OIN")

    files = []
    for sub in ("cat", "photo"):
        files += sorted(glob.glob(os.path.join(args.dataset_dir, sub, "*.xyz")))
    print(f"{len(files)} xyz files in cat/ + photo/")

    candidates = [p for p in files if boron_count_fast(p) >= 3]
    print(f"{len(candidates)} have >=3 boron (the only ones that can carry a cage)")

    known = set(BORON_COHORT)
    rows, tally = [], Counter()
    for p in candidates:
        name = os.path.basename(p).replace(".xyz", "").replace("_comp_0", "")
        try:
            r = analyse(p)
        except Exception as e:  # noqa: BLE001
            tally["analyse_error"] += 1
            rows.append({"mol": name, "error": f"{type(e).__name__}: {str(e)[:80]}"})
            continue
        if r is None:
            tally["boron_but_no_cage_motif"] += 1
            continue
        r["mol"] = name
        r["in_encode_fail_34"] = name in known
        r["counted_as_passing"] = name in passing
        rows.append(r)
        if r["BB_deleted"] == 0:
            tally["cage_intact_already"] += 1
        elif r["in_encode_fail_34"]:
            tally["cage_amputated__known_encode_fail"] += 1
        elif r["counted_as_passing"]:
            tally["cage_amputated__COUNTED_AS_PASSING"] += 1
        else:
            tally["cage_amputated__other_status"] += 1

    print("\n=== TALLY ===")
    for k, v in tally.most_common():
        print(f"{v:5d}  {k}")
    silent = [
        r
        for r in rows
        if r.get("BB_deleted") and not r.get("in_encode_fail_34") and r.get("counted_as_passing")
    ]
    if silent:
        tot_raw = sum(r["BB_raw"] for r in silent)
        tot_del = sum(r["BB_deleted"] for r in silent)
        print(
            f"\nSILENTLY CORRUPTED, counted as passing: {len(silent)} molecules, "
            f"{tot_del}/{tot_raw} B-B cage bonds deleted "
            f"({100.0 * tot_del / tot_raw:.1f}%)"
        )
        for r in sorted(silent, key=lambda r: -r["BB_deleted"])[:25]:
            print(
                f"  {r['mol']:12s} nB={r['nB']:3} B-B {r['BB_pruned']:3}/{r['BB_raw']:3} "
                f"deleted={r['BB_deleted']:3}"
            )
    with open(args.out, "w") as fh:
        json.dump({"tally": dict(tally), "rows": rows}, fh, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
