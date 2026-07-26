"""Is a 0-H carbon a missing hydrogen, or an under-assigned double bond?

The GAIN half of the `atom_count` class rests on the encoder recording carbons with
0 H. Whether that is *correct* decides whether the class is a defect at all, and the
input geometry answers it without any chemistry model:

* a carbon with 3 heavy neighbours and no H that is **planar** (small out-of-plane
  displacement, angles summing to ~360 deg) is **sp2** -- so a double bond was
  under-assigned by perception and no hydrogen is missing;
* the same carbon **pyramidal** (out-of-plane ~0.3-0.4 A, angles ~109 deg) is **sp3**,
  so the crystal structure simply never located its hydrogen and the input XYZ is
  incomplete.

The two verdicts point at opposite fixes, so they must not be conflated.
"""

from __future__ import annotations

import argparse
import collections
import glob
from pathlib import Path

import numpy as np

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
    "Se": 1.20,
    "Sn": 1.39,
}
DEFAULT_METAL_R = 1.45


def radius(sym: str) -> float:
    return COV.get(sym, DEFAULT_METAL_R)


def read(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[0].split()[0])
    els, xyz = [], []
    for ln in lines[2 : 2 + n]:
        p = ln.split()
        els.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return els, np.array(xyz)


def neighbours(els, xyz, tol=0.45):
    n = len(els)
    nb = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(xyz[i] - xyz[j]) < radius(els[i]) + radius(els[j]) + tol:
                nb[i].append(j)
                nb[j].append(i)
    return nb


def out_of_plane(centre, others) -> float:
    """Distance from `centre` to the best-fit plane of exactly three `others`."""
    a, b, c = others
    normal = np.cross(b - a, c - a)
    norm = np.linalg.norm(normal)
    if norm < 1e-9:
        return 0.0
    return abs(float(np.dot(centre - a, normal / norm)))


def angle_sum(centre, others) -> float:
    tot = 0.0
    for i in range(len(others)):
        for j in range(i + 1, len(others)):
            u = others[i] - centre
            v = others[j] - centre
            cos = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
            tot += np.degrees(np.arccos(max(-1.0, min(1.0, cos))))
    return tot


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("molecules", nargs="+")
    ap.add_argument(
        "--dataset",
        default="/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset",
    )
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    grand: collections.Counter = collections.Counter()
    for name in args.molecules:
        hits = glob.glob(f"{args.dataset}/*/{name}.xyz")
        if not hits:
            print(f"{name}: NOT FOUND")
            continue
        els, xyz = read(Path(hits[0]))
        nb = neighbours(els, xyz)
        rows = []
        for i, sym in enumerate(els):
            if sym != "C":
                continue
            h = [j for j in nb[i] if els[j] == "H"]
            heavy = [j for j in nb[i] if els[j] != "H"]
            if h or len(heavy) not in (2, 3):
                continue  # only 2- and 3-coordinate 0-H carbons are ambiguous
            if len(heavy) == 3:
                oop = out_of_plane(xyz[i], [xyz[j] for j in heavy])
                asum = angle_sum(xyz[i], [xyz[j] for j in heavy])
                verdict = (
                    "3-coord sp2 (double bond under-assigned; 0 H is CORRECT)"
                    if oop < 0.15
                    else "3-coord sp3 (1 H MISSING from input)"
                )
            else:
                oop = float("nan")
                asum = angle_sum(xyz[i], [xyz[j] for j in heavy])
                # One angle only. ~180 = sp (alkyne/nitrile, 0 H correct);
                # ~120 = sp2, so one H is missing; ~109 = sp3, two are missing.
                if asum > 155.0:
                    verdict = "2-coord sp (alkyne/nitrile; 0 H is CORRECT)"
                elif asum > 114.0:
                    verdict = "2-coord sp2 (1 H MISSING from input)"
                else:
                    verdict = "2-coord sp3 (2 H MISSING from input)"
            grand[verdict] += 1
            rows.append((i, oop, asum, verdict, [els[j] for j in heavy]))
        n_h = sum(1 for e in els if e == "H")
        print(f"\n{name}: {len(els)} atoms, {n_h} H, {len(rows)} ambiguous 0-H carbons")
        if not args.quiet:
            for i, oop, asum, verdict, nbs in rows[:20]:
                print(
                    f"   C{i:<4d} out-of-plane={oop:6.3f} A  angle-sum={asum:6.1f} deg  "
                    f"nbrs={','.join(nbs)}  -> {verdict}"
                )
            if len(rows) > 20:
                print(f"   ... {len(rows) - 20} more")

    print("\n=== verdict totals over all molecules ===")
    for k, v in grand.most_common():
        print(f"  {v:5d}  {k}")


if __name__ == "__main__":
    main()
