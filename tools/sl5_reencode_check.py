#!/usr/bin/env python
"""SL5 byte-identity gate: re-encode currently-passing molecules, diff vs stored OIN.

The SL5 encoder changes (AC2BO valence-combo cap, lig_checks resonance cap, boron typed
error) must keep every currently-encodable structure byte-identical. This re-encodes a
sample of the v0.4.2 capstone molecules whose stored ``smiles_1`` is non-null and compares
the fresh OIN against it. Any mismatch on a previously-passing molecule is a regression.

Because the caps only change behaviour for *large* ligands, ``--min-atoms`` targets the
risk set directly; run once over the large set and once as a broad random control.

Usage:
    PYTHONPATH=src python tools/sl5_reencode_check.py \
        --dataset-dir <abs> --reports-dir <abs>/results-capstone-v042/individual_reports \
        [--min-atoms 70] [--sample 800] [--seed 0]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import random
import signal
import warnings

warnings.filterwarnings("ignore")


def find_xyz(mol: str, dataset_dir: str) -> str | None:
    for sub in ("cat", "photo"):
        p = os.path.join(dataset_dir, sub, f"{mol}.xyz")
        if os.path.exists(p):
            return p
    return None


def atom_count(xyz_path: str) -> int:
    with open(xyz_path) as fh:
        try:
            return int(fh.readline().split()[0])
        except Exception:
            return 0


class _Timeout(Exception):
    pass


def _handler(signum, frame):
    raise _Timeout()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset-dir", required=True)
    ap.add_argument("--reports-dir", required=True)
    ap.add_argument("--min-atoms", type=int, default=0)
    ap.add_argument("--sample", type=int, default=0, help="0 = all matching")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--per-mol-timeout", type=int, default=120)
    args = ap.parse_args()

    from oinsmiles import XYZToSMILES

    signal.signal(signal.SIGALRM, _handler)

    # Collect (mol, stored_oin) for molecules that encoded (smiles_1 non-null).
    targets = []
    for rp in glob.glob(os.path.join(args.reports_dir, "*.json")):
        try:
            rec = json.load(open(rp))
        except Exception:
            continue
        stored = rec.get("smiles_1")
        mol = rec.get("molecule") or os.path.splitext(os.path.basename(rp))[0]
        if not stored:
            continue
        xyz = rec.get("input_xyz")
        if not xyz or not os.path.exists(xyz):
            xyz = find_xyz(mol, args.dataset_dir)
        if xyz is None or not os.path.exists(xyz):
            continue
        if args.min_atoms and atom_count(xyz) < args.min_atoms:
            continue
        targets.append((mol, stored, xyz))

    rng = random.Random(args.seed)
    rng.shuffle(targets)
    if args.sample:
        targets = targets[: args.sample]

    identical = diffs = now_fail = errors = 0
    diff_list, fail_list = [], []
    for i, (mol, stored, xyz) in enumerate(targets):
        signal.alarm(args.per_mol_timeout)
        try:
            fresh = XYZToSMILES().convert(xyz)
            signal.alarm(0)
            if fresh == stored:
                identical += 1
            else:
                diffs += 1
                diff_list.append(mol)
        except _Timeout:
            signal.alarm(0)
            now_fail += 1
            fail_list.append(f"{mol}:TIMEOUT")
        except Exception as e:  # noqa: BLE001
            signal.alarm(0)
            now_fail += 1
            fail_list.append(f"{mol}:{type(e).__name__}")
        if (i + 1) % 200 == 0:
            print(f"  ...{i + 1}/{len(targets)} done", flush=True)

    print(
        f"\nchecked={len(targets)}  identical={identical}  "
        f"DIFFERENT={diffs}  now_fail={now_fail}  errors={errors}",
        flush=True,
    )
    if diff_list:
        print("REGRESSION (OIN changed):", ", ".join(diff_list[:40]), flush=True)
    if fail_list:
        print("REGRESSION (now fails):", ", ".join(fail_list[:40]), flush=True)
    if not diff_list and not fail_list:
        print("PASS: every re-encoded molecule is byte-identical.", flush=True)


if __name__ == "__main__":
    main()
