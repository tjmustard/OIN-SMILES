"""Locate where an OIN string's hydrogen count diverges from the input XYZ.

Runs the real encoder and reports, per heavy atom, the H count the input
geometry supports versus the H count the encoder ends up serializing.  The
`atom_count` hard-fail class is 100% hydrogen (see docs/agentic-notes/v0.4.5/ATOM_COUNT_v0.4.5.md),
so this is the whole question.
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

import numpy as np
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

from oinsmiles.utils.xyz2mol import get_oin_string, get_tmc_mol  # noqa: E402

COV = {
    "H": 0.31,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
}


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[0].split()[0])
    els, xyz = [], []
    for ln in lines[2 : 2 + n]:
        p = ln.split()
        els.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return els, np.array(xyz)


def input_h_on_heavy(els, xyz):
    """For each heavy atom index, how many H are bonded to it in the raw XYZ."""
    n = len(els)
    counts = collections.Counter()
    for i in range(n):
        if els[i] != "H":
            continue
        best, bd = None, 1e9
        for j in range(n):
            if j == i or els[j] == "H":
                continue
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if d < bd:
                best, bd = j, d
        # 1.35 A covers C-H (1.09), N-H (1.01), O-H (0.96), S-H (1.34), B-H (1.19)
        if best is not None and bd < 1.45:
            counts[best] += 1
    return counts


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("molecules", nargs="+")
    ap.add_argument(
        "--dataset",
        default="/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset",
    )
    args = ap.parse_args()

    root = Path(args.dataset)
    for name in args.molecules:
        hits = list(root.glob(f"*/{name}.xyz"))
        if not hits:
            print(f"{name}: NO INPUT XYZ FOUND")
            continue
        p = hits[0]
        els, xyz = read_xyz(p)
        n_in = len(els)
        h_on = input_h_on_heavy(els, xyz)

        try:
            tmc, coords = get_tmc_mol(str(p), 0, with_stereo=True)
        except Exception as e:
            print(f"{name}: get_tmc_mol failed: {type(e).__name__}: {e}")
            continue
        try:
            oin = get_oin_string(tmc, coords)
        except Exception as e:
            print(f"{name}: get_oin_string failed: {type(e).__name__}: {e}")
            continue

        # what the encoder's own mol thinks, per atom, keyed by original XYZ index
        print(f"\n=== {name}  input={n_in} atoms ({sum(1 for e in els if e == 'H')} H) ===")
        print(f"    oin: {oin}")
        rows = []
        for a in tmc.GetAtoms():
            if a.GetAtomicNum() == 1:
                continue
            oi = a.GetIntProp("__origIdx") if a.HasProp("__origIdx") else None
            hs_in = h_on.get(oi, 0) if oi is not None else None
            nbrH = sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() == 1)
            tot = a.GetTotalNumHs()
            exp = a.GetNumExplicitHs()
            imp = a.GetNumImplicitHs()
            enc = nbrH + tot
            if hs_in is not None and enc != hs_in:
                rows.append(
                    (
                        a.GetIdx(),
                        oi,
                        a.GetSymbol(),
                        a.GetIsAromatic(),
                        hs_in,
                        nbrH,
                        exp,
                        imp,
                        tot,
                        a.GetTotalValence(),
                        a.GetFormalCharge(),
                        a.GetNumRadicalElectrons(),
                    )
                )
        if not rows:
            print("    tmc_mol H attribution: matches input on every heavy atom")
        else:
            print("    idx orig sym arom  inH nbrH exp imp tot val chg rad   <-- DIVERGENT")
            for r in rows:
                print(
                    f"    {r[0]:3d} {str(r[1]):>4s} {r[2]:>3s} {str(r[3]):>5s} "
                    f"{r[4]:4d} {r[5]:4d} {r[6]:3d} {r[7]:3d} {r[8]:3d} {r[9]:3d} "
                    f"{r[10]:3d} {r[11]:3d}"
                )


if __name__ == "__main__":
    sys.exit(main())
