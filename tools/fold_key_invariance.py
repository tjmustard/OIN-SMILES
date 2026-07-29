"""Is the donor fold GENERATOR-NEUTRAL? (v0.4.13)

THE QUESTION THIS DECIDES
=========================
v0.4.13 promotes ``OIN_CANONICAL_DONOR_FOLD`` + ``OIN_FOLD_PARITY_VETO`` to default-ON, and the
promotion's headline can be produced one of two ways:

* **offline re-score** of the frozen corpus's stored strings -- minutes, and the A/B confound
  dissolves because the SAME stored geometry goes through both arms; or
* **a fresh ~55 CPU-h generator sweep** -- which re-runs a STOCHASTIC generator, so the
  difference between two independent runs contaminates the measurement. This project's standing
  trap: *never A/B by re-running a stochastic harness.*

The offline route is exact **if and only if the promotion cannot change what the generator
returns.** The fold is encoder-side, but that is not sufficient on its own: the generator's
``accept_fn`` (``_reencode_key_matches``) accepts a conformer by comparing
``canonical_roundtrip_key(re-encoded OIN)`` against the target key. If the fold changed a key,
acceptance would change, a different conformer would be returned, and the offline re-score would
be describing structures the shipped code would never produce.

So the condition is precisely: **does folding ever change the round-trip comparison key?**

v0.4.11 measured this on 992 strings and read 0. This re-measures it on every string in the
frozen corpus, because "0 of 992" is a sample and the promotion is permanent.

The veto does not need its own arm here. It only ever returns the folded string or the
rotation-only one; if those two share a key, then every path through ``resolve`` shares it too.

READ THE DENOMINATOR
====================
A broken version of this -- one where the lever never actually toggled, or where every string
failed to parse and was skipped -- would print ``0 differing`` and look like a clean pass. So the
output states how many strings were compared, how many the fold actually MOVED (which must be
non-zero, or the lever is not wired), and how many were skipped.

ANY ENCODER-SIDE CANONICALIZATION LEVER (v0.4.14)
=================================================
The question generalizes without change: a lever that never moves a comparison key cannot move
``accept_fn``'s verdict. ``--lever`` names the lever under test and ``--holding`` names levers
forced ON in **both** arms.

``--holding`` is not a convenience. ``OIN_RESONANCE_DONOR_FOLD`` only widens a candidate set
that ``OIN_CANONICAL_DONOR_FOLD`` creates, so measuring it against a fold-OFF baseline would
report the *fold's* movement as the widening's -- a larger, more attractive, and wrong number.
Hold the fold ON in both arms and the delta is the widening's alone.

Usage
-----
    PYTHONPATH=src .venv/bin/python tools/fold_key_invariance.py \
        --sweep tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest

    PYTHONPATH=src .venv/bin/python tools/fold_key_invariance.py \
        --sweep ... --lever OIN_RESONANCE_DONOR_FOLD --holding OIN_CANONICAL_DONOR_FOLD
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from contextlib import ExitStack, contextmanager

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

from oinsmiles.oin.canonical_slots import canonicalize_oin_slots  # noqa: E402
from oinsmiles.oin.compare import canonical_roundtrip_key  # noqa: E402

FOLD = "OIN_CANONICAL_DONOR_FOLD"


@contextmanager
def forced(name, value):
    """Force a lever for the duration, restoring whatever was there before."""
    prev = os.environ.get(name)
    os.environ[name] = "1" if value else "0"
    try:
        yield
    finally:
        if prev is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = prev


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--out-json")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--lever", default=FOLD, help=f"lever under test (default {FOLD})")
    ap.add_argument(
        "--holding",
        action="append",
        default=[],
        help="lever forced ON in BOTH arms, repeatable (see module docstring)",
    )
    args = ap.parse_args()

    reports = sorted(glob.glob(os.path.join(args.sweep, "individual_reports", "*.json")))
    if args.limit:
        reports = reports[: args.limit]
    if not reports:
        sys.exit(f"FATAL: no individual_reports/*.json under {args.sweep}")

    n_cmp = n_moved = n_key_diff = n_skipped = 0
    key_diffs, movers = [], []

    for path in reports:
        try:
            rep = json.load(open(path))
        except (OSError, json.JSONDecodeError):
            n_skipped += 1
            continue
        mol = rep.get("molecule") or os.path.basename(path)[:-5]
        for field in ("smiles_1", "smiles_2_indep"):
            s = rep.get(field)
            if not s:
                continue
            try:
                with ExitStack() as held:
                    for name in args.holding:
                        held.enter_context(forced(name, True))
                    with forced(args.lever, False):
                        off = canonicalize_oin_slots(s)
                        k_off = canonical_roundtrip_key(off)
                    with forced(args.lever, True):
                        on = canonicalize_oin_slots(s)
                        k_on = canonical_roundtrip_key(on)
            except Exception:  # noqa: BLE001 -- an unparseable string is not evidence either way
                n_skipped += 1
                continue
            n_cmp += 1
            if off != on:
                n_moved += 1
                movers.append(f"{mol}:{field}")
            if k_off != k_on:
                n_key_diff += 1
                key_diffs.append(f"{mol}:{field}")

    held = " + ".join(args.holding) or "(none)"
    print(f"# key-invariance of {args.lever} over {args.sweep}")
    print(f"# levers held ON in BOTH arms: {held}\n")
    print(f"  strings compared          : {n_cmp}")
    print(f"  strings the LEVER MOVED   : {n_moved}   <-- must be > 0 or the lever is not wired")
    print(f"  strings whose KEY CHANGED : {n_key_diff}")
    print(f"  skipped (unparseable)     : {n_skipped}")

    print()
    if n_moved == 0:
        print(f"  🔴 {args.lever} MOVED NOTHING. This measured nothing -- it is not active.")
        verdict = "INCONCLUSIVE"
    elif n_key_diff == 0:
        print("  ✅ GENERATOR-NEUTRAL. The fold changes strings but never the comparison key, so")
        print("     accept_fn's verdict is identical with the lever on or off, the generator")
        print("     returns the same conformers, and an OFFLINE RE-SCORE of the frozen corpus is")
        print("     EXACT. A fresh generator sweep would add stochastic noise and nothing else.")
        verdict = "GENERATOR_NEUTRAL"
    else:
        print(f"  🔴 NOT generator-neutral: {n_key_diff} keys change. accept_fn would decide")
        print("     differently, so the offline re-score would describe conformers the shipped")
        print("     code never returns. A FULL SWEEP IS REQUIRED.")
        print(f"     first: {key_diffs[:8]}")
        verdict = "SWEEP_REQUIRED"

    if args.out_json:
        json.dump(
            {
                "sweep": args.sweep,
                "lever": args.lever,
                "holding": args.holding,
                "verdict": verdict,
                "n_compared": n_cmp,
                "n_moved": n_moved,
                "n_key_changed": n_key_diff,
                "n_skipped": n_skipped,
                "key_diffs": key_diffs[:200],
                "movers_sample": movers[:200],
            },
            open(args.out_json, "w"),
            indent=2,
        )
        print(f"\nwrote {args.out_json}")
    print(f"\n#DONE {n_cmp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
