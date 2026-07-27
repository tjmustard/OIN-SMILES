#!/usr/bin/env python
"""Summarize a cohort of `probe_accept_gap.py` runs into the four cases that matter.

The cohort exists because one molecule (HIDCIH_comp_1) produced the hypothesis and a second
(YIYGAP_comp_0) immediately contradicted it. The useful output is not an average -- it is the
population split, because each case implies a different fix:

  GAP        cheap fires strictly earlier than strict
             -> OIN_ACCEPT_SCORED recovers (first_strict_i - first_cheap_i) conformers
  NO_GAP     both fire at the same index -> the lever changes nothing here
  CHEAP_ONLY cheap fires, strict NEVER does
             -> the molecule only ever passes because the SCORE is the cheap test; the lever
                is what makes it exit early, and there is nothing independently confirmable
  PREFILTER_VETO  strict fires but cheap does not
             -> production accepts NEITHER: `_reencode_key_matches` returns False on a cheap
                mismatch before the strict test ever runs, so these conformers are unreachable
                today. This is the one case that is a potential ACCURACY defect rather than a
                latency one, and it runs opposite to the lever's hypothesis.
  DEAD       neither fires -> not an acceptance-gap molecule at all

Usage:  python tools/summarize_accept_gap.py <dir-of-json> [--md out.md]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import statistics


def classify(r: dict) -> str:
    fc, fs = r.get("first_cheap_match"), r.get("first_strict_match")
    if fc and fs:
        return "GAP" if fc["i"] < fs["i"] else ("NO_GAP" if fc["i"] == fs["i"] else "STRICT_FIRST")
    if fc and not fs:
        return "CHEAP_ONLY"
    if fs and not fc:
        return "PREFILTER_VETO"
    return "DEAD"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--md")
    args = ap.parse_args()

    rows = []
    for f in sorted(glob.glob(os.path.join(args.dir, "*.json"))):
        try:
            r = json.load(open(f))
        except Exception:
            continue
        if "n_conformers_seen" not in r:
            continue
        fc, fs = r.get("first_cheap_match"), r.get("first_strict_match")
        rows.append(
            {
                "mol": r["molecule"].replace(".xyz", ""),
                "case": classify(r),
                "seen": r["n_conformers_seen"],
                "n_cheap": r["n_cheap_match"],
                "n_strict": r["n_strict_match"],
                "cheap_i": fc["i"] if fc else None,
                "cheap_t": fc["t"] if fc else None,
                "strict_i": fs["i"] if fs else None,
                "strict_t": fs["t"] if fs else None,
                "total_s": r["total_s"],
            }
        )

    order = ["CHEAP_ONLY", "GAP", "PREFILTER_VETO", "NO_GAP", "STRICT_FIRST", "DEAD"]
    rows.sort(key=lambda r: (order.index(r["case"]), -(r["total_s"] or 0)))

    hdr = (
        f"{'molecule':22s} {'case':15s} {'seen':>4} {'cheap':>6} {'strict':>7} "
        f"{'c@i':>5} {'s@i':>5} {'c_t':>7} {'s_t':>7} {'total':>8}"
    )
    lines = [hdr, "-" * len(hdr)]
    for r in rows:
        lines.append(
            f"{r['mol']:22s} {r['case']:15s} {r['seen']:4d} {r['n_cheap']:6d} {r['n_strict']:7d} "
            f"{str(r['cheap_i']):>5} {str(r['strict_i']):>5} "
            f"{str(r['cheap_t']):>7} {str(r['strict_t']):>7} {r['total_s']:8.1f}"
        )

    counts: dict[str, int] = {}
    for r in rows:
        counts[r["case"]] = counts.get(r["case"], 0) + 1
    lines += ["", f"n = {len(rows)}"]
    for c in order:
        if counts.get(c):
            lines.append(f"  {c:15s} {counts[c]:3d}  ({100 * counts[c] / len(rows):.0f}%)")

    # Recoverable time: for molecules where a scored-success exists, how much of the observed
    # wall-clock came AFTER it? That is what the lever can return, and no more.
    rec = [
        (r["mol"], r["total_s"] - r["cheap_t"], r["total_s"])
        for r in rows
        if r["cheap_t"] is not None and r["total_s"] is not None
    ]
    if rec:
        saved = [s for _, s, _ in rec]
        lines += [
            "",
            "RECOVERABLE WALL-CLOCK (time spent after the first scored-success appeared):",
            f"  molecules with a scored-success : {len(rec)}",
            f"  median recoverable             : {statistics.median(saved):.1f}s",
            f"  total recoverable              : {sum(saved):.1f}s of "
            f"{sum(t for _, _, t in rec):.1f}s observed",
        ]
        worst = sorted(rec, key=lambda x: -x[1])[:6]
        for m, s, t in worst:
            lines.append(f"    {m:22s} {s:7.1f}s of {t:7.1f}s")

    text = "\n".join(lines)
    print(text)
    if args.md:
        with open(args.md, "w") as fh:
            fh.write("```\n" + text + "\n```\n")
        print(f"\nwrote {args.md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
