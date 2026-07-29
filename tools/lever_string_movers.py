#!/usr/bin/env python
"""The EXACT population a lever can affect, derived from coordinates. (v0.4.14)

WHY THIS IS THE POPULATION, AND WHY THAT MAKES A SWEEP UNNECESSARY
=================================================================
For an encoder-side lever, a molecule's round-trip verdict can only change if the lever changes a
string. Concretely:

* if ``encode(input.xyz)`` is byte-identical in both arms, the generator receives the **same input
  string**, and generation is seeded (``MetalloGenAdapter(seed=42)``), so it returns the same
  structure; and
* if ``encode(that structure)`` is also byte-identical, the verdict is identical.

So the molecules where **neither** encode moves are unchanged **by construction**, and an A/B
restricted to the movers is not a sample — it is **exact**. That is what lets a release re-measure a
canonicalization lever for the cost of the movers instead of a ~55 CPU-h corpus sweep.

⚠ **THIS IS NOT THE SAME AS `fold_key_invariance.py`'s LICENCE, AND THE DIFFERENCE IS THE POINT.**
That tool asks whether the lever moves a comparison KEY, and a `GENERATOR_NEUTRAL` verdict was read
as licensing an offline re-score over FROZEN structures. It does not: `accept_fn` decides by key, so
key-invariance bounds *acceptance*, but the generator's input is the OIN **string**, so a slot
relabeling changes `ParsedOIN`, the CoordMap and the pool itself. This tool takes the string as the
unit precisely because the string is what the generator consumes.

⚠ **DERIVE THE POPULATION FROM COORDINATES, NOT FROM A FROZEN SWEEP'S STORED STRINGS.** A stored
``smiles_1`` was emitted by whatever encoder ran at freeze time; several levers have been promoted
since. Re-encoding the input XYZ under today's encoder is the only way the mover set describes the
code that actually ships. v0.4.14's first pass used stored strings and could not have detected a
molecule whose *current* emitted string differs.

WHAT IT EMITS
=============
Two lists, because they are at risk in opposite directions and want different sample sizes:

``input_moved``   the lever changes ``encode(input.xyz)`` -> generation changes -> anything can
                  happen. Contains both the potential gains and the potential losses.
``candidates``    ``input_moved`` plus molecules whose input is unchanged (so the generated
                  structure is unchanged) but which are worth carrying anyway when the caller has a
                  stored structure to re-encode.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/lever_string_movers.py \\
        --cohort-dir <MAIN>/tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k \\
        --lever OIN_RESONANCE_DONOR_FOLD --out-json movers.json --out-list movers.txt
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def _encode(path):
    from oinsmiles import XYZToSMILES

    try:
        return XYZToSMILES().convert(path)
    except Exception:  # noqa: BLE001 -- an unencodable input cannot move, and is counted
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort-dir", required=True)
    ap.add_argument("--lever", required=True)
    ap.add_argument(
        "--holding",
        action="append",
        default=[],
        help="lever forced ON in BOTH arms, repeatable (e.g. the fold a widening depends on)",
    )
    ap.add_argument("--limit", type=int)
    ap.add_argument("--out-json")
    ap.add_argument("--out-list", help="newline-separated mover names, for generator_ab_honest.py")
    args = ap.parse_args()

    if not os.path.isdir(args.cohort_dir):
        sys.exit(f"🔴 --cohort-dir {args.cohort_dir} is not a directory")
    paths = sorted(glob.glob(os.path.join(args.cohort_dir, "*.xyz")))
    if args.limit:
        paths = paths[: args.limit]
    if not paths:
        sys.exit(f"🔴 REFUSING: no *.xyz under {args.cohort_dir}")

    for h in args.holding:
        os.environ[h] = "1"

    moved, unmoved, failed = [], [], []
    for i, p in enumerate(paths, 1):
        mol = os.path.splitext(os.path.basename(p))[0]
        os.environ[args.lever] = "0"
        off = _encode(p)
        os.environ[args.lever] = "1"
        on = _encode(p)
        if off is None or on is None:
            failed.append(mol)
        elif off != on:
            moved.append(mol)
        else:
            unmoved.append(mol)
        if i % 250 == 0:
            print(f"  [{i}/{len(paths)}] moved={len(moved)} failed={len(failed)}", flush=True)

    n = len(paths)
    print(f"\n=== {args.lever}: exact mover population over {n} molecules ===")
    print(f"  holding ON in both arms : {' + '.join(args.holding) or '(none)'}")
    print(f"  input string MOVED      : {len(moved)}   <-- the exact at-risk population")
    print(f"  unchanged               : {len(unmoved)}")
    print(f"  encode failed either arm: {len(failed)}")

    # 🔴 A zero here is not a small effect -- it means the lever never fired, which is also what a
    # mis-spelled lever name prints. Refuse rather than hand the caller an empty cohort.
    if not moved:
        sys.exit(
            f"\n🔴 REFUSING: {args.lever} moved 0 of {n} emitted strings. Either it is not wired, "
            "the name is misspelled, or a --holding dependency is missing."
        )

    if args.out_list:
        with open(args.out_list, "w") as fh:
            fh.write("\n".join(moved) + "\n")  # trailing newline: a missing one merged two names
        print(f"\nwrote {args.out_list} ({len(moved)} names)")
    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(
                {
                    "lever": args.lever,
                    "holding": args.holding,
                    "cohort_dir": args.cohort_dir,
                    "n_total": n,
                    "moved": moved,
                    "n_moved": len(moved),
                    "n_unmoved": len(unmoved),
                    "encode_failed": failed,
                },
                fh,
                indent=1,
            )
        print(f"wrote {args.out_json}")


if __name__ == "__main__":
    main()
