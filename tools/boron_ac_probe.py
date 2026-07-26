#!/usr/bin/env python
"""Is the cage lost to geometry, or to a valence RULE?

``xyz2AC_obabel`` builds the adjacency matrix from a covalent-radius distance
criterion, then runs a pruning loop that deletes an atom's longest bonds while
its connectivity exceeds ``max(atomic_valence[Z])``.  For boron that cap is 4,
while a closo/nido cage vertex has 5-6 neighbours.  So the pruning loop shatters
every cage it perceives.

This probe compares, per molecule:
  * B degrees from the *raw distance criterion* (no pruning),
  * B degrees from the encoder's *actual* AC (pruned),
  * how many B-B edges the pruning loop deletes,
so the reframe is a measurement, not an inference.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np


def find_xyz(mol: str, dataset_dir: str) -> str | None:
    for sub in ("cat", "photo"):
        p = os.path.join(dataset_dir, sub, f"{mol}_comp_0.xyz")
        if os.path.exists(p):
            return p
    return None


def raw_ac(atoms, xyz, tolerance=0.5):
    """The distance criterion ONLY -- no valence pruning."""
    from rdkit import Chem

    pt = Chem.GetPeriodicTable()
    n = len(atoms)
    coords = np.asarray(xyz, dtype=float)
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    rcov = np.array([pt.GetRcovalent(int(z)) for z in atoms])
    thr = rcov[:, None] + rcov[None, :] + tolerance
    AC = ((d <= thr) & ~np.eye(n, dtype=bool)).astype(int)
    return AC, d


def probe(mol_name: str, dataset_dir: str, tolerance: float = 0.5) -> dict:
    import warnings

    warnings.filterwarnings("ignore")
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    from oinsmiles.utils.xyz2mol_local import read_xyz_file, xyz2AC_obabel

    path = find_xyz(mol_name, dataset_dir)
    if path is None:
        return {"mol": mol_name, "err": "NO_FILE"}
    atoms, _charge, xyz = read_xyz_file(path)

    AC_raw, dmat = raw_ac(atoms, xyz, tolerance)
    AC_pruned, _mol = xyz2AC_obabel(atoms, xyz, tolerance=tolerance)
    AC_pruned = np.asarray(AC_pruned)

    b = [i for i, z in enumerate(atoms) if z == 5]
    deg_raw = sorted(int(AC_raw[i].sum()) for i in b)
    deg_pruned = sorted(int(AC_pruned[i].sum()) for i in b)

    def count_bb(AC):
        return int(sum(AC[i, j] for ii, i in enumerate(b) for j in b[ii + 1 :]))

    bb_raw, bb_pruned = count_bb(AC_raw), count_bb(AC_pruned)

    # every B-B edge the pruning loop deleted, with its length
    deleted = []
    for ii, i in enumerate(b):
        for j in b[ii + 1 :]:
            if AC_raw[i, j] and not AC_pruned[i, j]:
                deleted.append(round(float(dmat[i, j]), 3))
    # every deleted edge, any element pair
    deleted_any = 0
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if AC_raw[i, j] and not AC_pruned[i, j]:
                deleted_any += 1

    return {
        "mol": mol_name,
        "nB": len(b),
        "B_deg_raw": deg_raw,
        "B_deg_pruned": deg_pruned,
        "max_B_deg_raw": max(deg_raw) if b else 0,
        "max_B_deg_pruned": max(deg_pruned) if b else 0,
        "BB_raw": bb_raw,
        "BB_pruned": bb_pruned,
        "BB_deleted": bb_raw - bb_pruned,
        "deleted_BB_lengths": sorted(deleted),
        "deleted_edges_total": deleted_any,
    }


def main() -> None:
    from boron_characterize import BORON_COHORT

    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset-dir", required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--tolerance", type=float, default=0.5)
    ap.add_argument("--out", default="tools/boron_ac_probe.json")
    args = ap.parse_args()
    todo = args.only.split(",") if args.only else BORON_COHORT
    res = []
    for m in todo:
        r = probe(m, args.dataset_dir, args.tolerance)
        res.append(r)
        print(
            f"{m:10s} nB={r.get('nB'):3} maxdeg raw={r.get('max_B_deg_raw'):2} "
            f"pruned={r.get('max_B_deg_pruned'):2}  B-B raw={r.get('BB_raw'):3} "
            f"pruned={r.get('BB_pruned'):3} DELETED={r.get('BB_deleted'):3}",
            flush=True,
        )
    with open(args.out, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main()
