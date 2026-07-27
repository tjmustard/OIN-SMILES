#!/usr/bin/env python
"""Measure the gap between `accept_fn`'s predicate and the harness's success predicate.

Why this exists
---------------
The eta runtime tail was twice attributed by reasoning about *which stage* `accept_fn` sees,
and both attributions were imprecise. A cProfile run on HIDCIH_comp_1 (96 s generate) showed
`accept_fn` evaluated 26 times and returning False every time, on a molecule the sweep scores
as a success. So the pool fills to completion even though a scoring conformer exists inside it.

The two predicates are NOT the same test:

  harness success  :  canonical_roundtrip_key(oin_in)
                        == canonical_roundtrip_key(get_oin_string(gen.mol, coords))
  accept_fn        :  fast contract-mol prefilter (may only REJECT)
                      AND  canonical_roundtrip_key(XYZToSMILES().convert(xyz)) == target

`accept_fn` is therefore a SUBSET of harness-success: it adds an independent re-perception
that the score never asks for. Any conformer where the cheap test passes and the strict one
fails is a conformer the pool fill could have stopped on and did not.

What it records
---------------
Per conformer reaching the predicate, in fill order: the cheap verdict, the strict verdict,
and the wall-clock at which it was seen. The reportable number is the first index at which
each verdict turns True -- the difference between them, times the per-attempt cost, is the
recoverable time.

This measures; it changes nothing. `_reencode_key_matches` is monkeypatched to record and then
return False unconditionally, so the pool always fills and every conformer is observed. That
means the run is NOT a normal generate() and its output structure is meaningless -- only the
verdict trace is.

Usage:
    python tools/probe_accept_gap.py <mol.xyz> [--timeout 300] [--json out.json]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.generation import metallogen_adapter as MA  # noqa: E402
from oinsmiles.oin.compare import canonical_roundtrip_key  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz")
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--json")
    ap.add_argument(
        "--stop-after-both",
        action="store_true",
        help="stop the fill once BOTH verdicts have fired once (cohort-affordable)",
    )
    args = ap.parse_args()

    name = os.path.basename(args.xyz)
    oin_in = XYZToSMILES().convert(args.xyz)
    target = canonical_roundtrip_key(oin_in)

    trace: list[dict] = []
    t_start = time.monotonic()

    real_fast = MA._reencode_oin_fast
    real_full = MA._reencode_oin
    real_build = MA.build_contract_mol
    real_predicate = MA._reencode_key_matches

    def patched(parsed, m, target_key, cmol=None, require_no_stretch=False, cache=None, **kw):
        """Record both verdicts on this conformer, then always reject so the pool fills.

        ``**kw`` absorbs predicate kwargs this probe does not model (``independent_confirm``).
        It must not raise on an unexpected one: ``_file_and_maybe_stop`` swallows any exception
        from ``accept_fn`` and logs it at debug, so a signature mismatch here would not fail
        loudly -- it would silently report "no conformer ever matched" for every molecule.
        """
        rec: dict = {"i": len(trace), "t": round(time.monotonic() - t_start, 2)}
        # The key the adapter is testing against must be the key we computed from the input,
        # or the trace measures a different question than the one asked.
        if target_key != target:
            rec["TARGET_KEY_MISMATCH"] = True
        try:
            c = cmol if cmol is not None else real_build(parsed, m)
            fast = real_fast(c) if c is not None else None
            rec["fast_encoded"] = fast is not None
            rec["cheap_match"] = bool(
                fast is not None and canonical_roundtrip_key(fast) == target_key
            )
        except Exception as e:
            rec["cheap_match"] = False
            rec["cheap_error"] = f"{type(e).__name__}: {e}"
        try:
            full = real_full(m)
            rec["strict_encoded"] = full is not None
            rec["strict_match"] = bool(
                full is not None and canonical_roundtrip_key(full) == target_key
            )
        except Exception as e:
            rec["strict_match"] = False
            rec["strict_error"] = f"{type(e).__name__}: {e}"
        rec["t_after"] = round(time.monotonic() - t_start, 2)
        trace.append(rec)
        if args.stop_after_both:
            # Both first-hit indices are known once each verdict has fired once; stopping
            # there costs about what an unpatched run costs (it also stops at the first
            # strict hit), which keeps a cohort run affordable. Without this the pool fills
            # to target_pool and every conformer is observed -- more complete, ~2x slower.
            if any(r.get("cheap_match") for r in trace) and any(
                r.get("strict_match") for r in trace
            ):
                return True
        return False  # force a full pool fill so every conformer is observed

    MA._reencode_key_matches = patched
    gen = MA.OIN3DGeneratorMetallogen(
        optimizer=None, ensemble_size=1, timeout=args.timeout, ff_params=None
    )
    err = None
    try:
        gen.generate(oin_in)
    except Exception as e:
        err = f"{type(e).__name__}: {e}"
    finally:
        MA._reencode_key_matches = real_predicate
    total = time.monotonic() - t_start

    def first(key):
        for r in trace:
            if r.get(key):
                return r
        return None

    fc, fs = first("cheap_match"), first("strict_match")
    out = {
        "molecule": name,
        "oin_in": oin_in,
        "generate_error": err,
        "total_s": round(total, 2),
        "n_conformers_seen": len(trace),
        "first_cheap_match": fc,
        "first_strict_match": fs,
        "n_cheap_match": sum(1 for r in trace if r.get("cheap_match")),
        "n_strict_match": sum(1 for r in trace if r.get("strict_match")),
        "n_cheap_only": sum(1 for r in trace if r.get("cheap_match") and not r.get("strict_match")),
        "trace": trace,
    }

    print(f"=== {name} ===")
    print(f"  conformers seen : {len(trace)}   total {total:.1f}s")
    print(
        f"  cheap  matches  : {out['n_cheap_match']}   first at {fc['i'] if fc else '-'}"
        f" (t={fc['t'] if fc else '-'}s)"
    )
    print(
        f"  strict matches  : {out['n_strict_match']}   first at {fs['i'] if fs else '-'}"
        f" (t={fs['t'] if fs else '-'}s)"
    )
    print(f"  cheap-only (recoverable early exits): {out['n_cheap_only']}")
    if fc and not fs:
        print("  => VERDICT: a scoring conformer exists that accept_fn REJECTS.")
    elif fc and fs and fc["i"] < fs["i"]:
        print(f"  => VERDICT: cheap fires {fs['i'] - fc['i']} conformers earlier than strict.")
    elif not fc:
        print("  => VERDICT: no conformer matched even the cheap test (not an acceptance gap).")
    else:
        print("  => VERDICT: predicates agree on the first hit; no gap on this molecule.")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)
        print(f"  wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
