"""Does ``OIN_CANONICAL_DONOR_FOLD`` ever make a token reflection-invariant? (v0.4.11 Lane 2)

THE RISK THIS EXISTS TO CATCH
=============================
``OIN_CANONICAL_SLOTS`` folds only the geometry's **proper-rotation** group, because an
improper operation maps a structure to its mirror image. v0.4.11's donor fold deliberately
widens past that boundary -- within one fragment, one symmetry class, one colour -- so the
question "did the widening collapse an enantiomer?" is not a formality.

A fixture alone has never been sufficient here. The Y2 axial wave sorted a tie-break on a
stereochemical **sign**, which made the emitted token reflection-invariant; every local test
passed and only a corpus-wide mirror audit caught it. So this runs the whole comparison at
corpus scale.

WHAT IT MEASURES
================
For each molecule: encode the structure and its **mirror** (a z-axis reflection), once with
the fold OFF and once with it ON.

    OFF_distinct   the shipped encoder already separates the enantiomers
    ON_distinct    the encoder with the fold separates them

    REGRESSION  <=>  OFF_distinct and not ON_distinct

That implication is the whole verdict, and it is the right one: a pair the shipped encoder
*already* folds (a metal Delta/Lambda whose descriptor is not emitted -- ``OIN_EMIT_METAL_CONFIG``
is held off, so ZUMNEC's enantiomers share a comparison key today) is a pre-existing gap and
belongs to v0.4.16, not to this lever. This tool must not report it as a regression, and it
must not let the lever hide behind it either -- hence the four-way table, not a pass count.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/mirror_audit_donor_fold.py \\
        --dataset <dir> --n 300 [--seed 42] [--out-json audit.json]
"""

import argparse
import json
import os
import random
import sys
import tempfile
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.oin.compare import canonical_roundtrip_key  # noqa: E402

LEVER = "OIN_CANONICAL_DONOR_FOLD"


def read_xyz(path):
    with open(path) as fh:
        lines = fh.readlines()
    n = int(lines[0].split()[0])
    syms, coords = [], []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        syms.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return syms, np.asarray(coords, dtype=float)


def write_xyz(path, syms, coords):
    with open(path, "w") as fh:
        fh.write(f"{len(syms)}\n\n")
        for sym, xyz in zip(syms, coords):
            fh.write(f"{sym:<3} {xyz[0]:>16.10f} {xyz[1]:>16.10f} {xyz[2]:>16.10f}\n")


def _encode(path, fold):
    """Encode ``path`` with the fold forced on or off. ``None`` if the encoder declines."""
    os.environ[LEVER] = "1" if fold else "0"
    try:
        return XYZToSMILES().convert(path)
    except Exception:  # noqa: BLE001  -- an unencodable molecule is not a mirror verdict
        return None


def _key(s):
    if s is None:
        return None
    try:
        return canonical_roundtrip_key(s)
    except Exception:  # noqa: BLE001
        return None


def audit_one(path, tmpdir):
    """``(verdict, detail)`` for one molecule."""
    syms, coords = read_xyz(path)
    mirror = coords * np.array([1.0, 1.0, -1.0])
    mpath = os.path.join(tmpdir, "mirror_" + os.path.basename(path))
    write_xyz(mpath, syms, mirror)

    row = {}
    for label, fold in (("off", False), ("on", True)):
        a, b = _encode(path, fold), _encode(mpath, fold)
        if a is None or b is None:
            return "encode_failed", {}
        row[label] = {
            "raw_distinct": a != b,
            "key_distinct": _key(a) != _key(b),
            "self": a,
            "mirror": b,
        }

    off, on = row["off"], row["on"]
    if off["raw_distinct"] and not on["raw_distinct"]:
        return "REGRESSION_raw_collapsed", row
    if off["key_distinct"] and not on["key_distinct"]:
        return "REGRESSION_key_collapsed", row
    if not off["raw_distinct"]:
        return "achiral_or_preexisting_fold", row
    if on["raw_distinct"] and not off["raw_distinct"]:
        return "improved", row
    return "distinct_both_arms", row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dataset", required=True, help="dir searched recursively for *.xyz")
    ap.add_argument("--n", type=int, default=300, help="molecules to sample (0 = all)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-json", help="write the per-molecule verdicts here")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    paths = []
    for root, _dirs, files in os.walk(args.dataset):
        paths.extend(os.path.join(root, f) for f in files if f.endswith(".xyz"))
    paths.sort()
    if not paths:
        sys.exit(f"no .xyz found under {args.dataset}")
    if args.n and args.n < len(paths):
        paths = random.Random(args.seed).sample(paths, args.n)
        paths.sort()

    tmpdir = tempfile.mkdtemp(prefix="oin-mirror-audit-")
    tally, rows, regressions = Counter(), [], []
    for i, p in enumerate(paths, 1):
        verdict, detail = audit_one(p, tmpdir)
        tally[verdict] += 1
        name = os.path.basename(p)[:-4]
        rows.append({"molecule": name, "verdict": verdict})
        if verdict.startswith("REGRESSION"):
            regressions.append((name, detail))
        if args.verbose or i % 50 == 0:
            print(f"  [{i}/{len(paths)}] {name:28} {verdict}", flush=True)

    print(f"\n=== mirror audit, {LEVER}: {len(paths)} molecules ===")
    for k, n in tally.most_common():
        print(f"  {k:32} {n:5}  {100 * n / len(paths):5.1f}%")

    n_reg = sum(v for k, v in tally.items() if k.startswith("REGRESSION"))
    if n_reg:
        print(f"\n!! {n_reg} REGRESSION(S) -- the fold made a token reflection-invariant:")
        for name, detail in regressions[:10]:
            print(f"  {name}")
            print(f"    off self  : {detail['off']['self']}")
            print(f"    off mirror: {detail['off']['mirror']}")
            print(f"    on  self  : {detail['on']['self']}")
            print(f"    on  mirror: {detail['on']['mirror']}")
    else:
        print("\nCLEAN: no molecule lost a mirror distinction it had with the lever off.")

    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump({"tally": dict(tally), "records": rows}, fh, indent=1)
        print(f"wrote {args.out_json}")
    sys.exit(1 if n_reg else 0)


if __name__ == "__main__":
    main()
