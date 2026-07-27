"""Summarise a ``tools/budget_bound_ab.py`` run: overrun, epsilon, and what the bound cost.

Reports four things, and refuses to conflate them:

1. **The overrun ratio, per arm.** ``spent_s / budget_s``. This is the number the release
   exists to move, and it is only meaningful per *generation call* -- never from
   ``metrics.elapsed_s``, which sums up to three separately SIGKILLed harness attempts.
2. **eps = max(spent_s) - budget** on the bounded arm, over rows that were **not** killed
   by the outer watchdog. A killed row's ``spent_s`` is the *kill time*, i.e. a lower bound
   on its true spend; averaging it in would understate exactly the defect being measured.
   Killed rows are counted and reported separately, never folded into eps.
3. **What the bound cost.** Molecules that produced a structure with the lever OFF and none
   with it ON are *converted* late successes, not regressions. They are counted apart from
   molecules that fail in both arms, because a tighter bound turning a late success into an
   honest failure is the intended behaviour and will otherwise read as damage.
4. **Byte-identity.** For molecules that completed in BOTH arms, the generated XYZ length is
   compared as a cheap change-detector. A bound may change *which* molecules finish; it must
   not change *what* a finishing molecule produces.

Usage:
    PYTHONPATH=src .venv/bin/python tools/budget_bound_report.py --ab <ab.jsonl>
"""

from __future__ import annotations

import argparse
import json
import statistics as st
import sys


def load(path):
    rows = {"OFF": {}, "ON": {}}
    n = 0
    done = None
    for line in open(path):
        line = line.strip()
        if line.startswith("#DONE"):
            done = int(line.split()[1])
            continue
        if not line.startswith("{"):
            continue
        r = json.loads(line)
        rows[r["arm"]][r["molecule"]] = r
        n += 1
    # The sentinel's denominator is checked BEFORE anything is trusted -- a truncated
    # file must never be able to look like a complete measurement.
    if done is None:
        sys.exit(f"error: no #DONE sentinel in {path} -- the run was killed or is still going")
    if done != n:
        sys.exit(f"error: #DONE {done} but {n} rows parsed -- refusing to trust a truncated run")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ab", required=True)
    args = ap.parse_args()
    rows = load(args.ab)

    common = sorted(set(rows["OFF"]) & set(rows["ON"]))
    print(f"molecules measured in both arms: {len(common)}\n")

    budget = rows["ON"][common[0]]["budget_s"]
    print(f"requested budget: {budget}s\n")

    print(
        f"{'arm':4s} {'n':>4s} {'killed':>7s} {'median':>8s} {'max':>9s} {'max ratio':>10s} {'>budget':>8s}"
    )
    for arm in ("OFF", "ON"):
        rs = [rows[arm][m] for m in common]
        live = [r for r in rs if r.get("outcome") != "killed"]
        spent = [r["spent_s"] for r in live]
        killed = len(rs) - len(live)
        over = sum(1 for s in spent if s > budget)
        print(
            f"{arm:4s} {len(rs):4d} {killed:7d} {st.median(spent):7.1f}s {max(spent):8.1f}s "
            f"{max(spent) / budget:9.2f}x {over:8d}"
        )

    on_live = [rows["ON"][m] for m in common if rows["ON"][m].get("outcome") != "killed"]
    eps = max(r["spent_s"] for r in on_live) - budget
    worst = max(on_live, key=lambda r: r["spent_s"])
    print(
        f"\neps = max(spent) - budget = {eps:+.1f}s   (worst: {worst['molecule']} "
        f"at {worst['spent_s']:.1f}s, encode {worst.get('encode_s', 0):.1f}s)"
    )
    n_killed_on = sum(1 for m in common if rows["ON"][m].get("outcome") == "killed")
    if n_killed_on:
        print(
            f"  ⚠ {n_killed_on} bounded row(s) hit the OUTER kill -- eps excludes them, and "
            "their true spend is unknown (>= the kill time). The bound did not hold there."
        )

    def produced(r):
        return r.get("outcome") == "ok"

    both = [m for m in common if produced(rows["OFF"][m]) and produced(rows["ON"][m])]
    lost = [m for m in common if produced(rows["OFF"][m]) and not produced(rows["ON"][m])]
    gained = [m for m in common if not produced(rows["OFF"][m]) and produced(rows["ON"][m])]
    neither = [m for m in common if not produced(rows["OFF"][m]) and not produced(rows["ON"][m])]

    print("\nwhat the bound cost -- CONVERTED late successes, counted apart from real failures:")
    print(f"  produced in both arms          : {len(both)}")
    print(f"  OFF produced, ON did not       : {len(lost)}   <- converted, INTENDED")
    print(f"  ON produced, OFF did not       : {len(gained)}")
    print(f"  neither arm produced           : {len(neither)}  <- real failures, unchanged")
    if lost:
        print(f"    converted: {lost}")
    if gained:
        print(f"    gained (OFF was killed by the outer watchdog before finishing): {gained}")

    diff = [m for m in both if rows["OFF"][m].get("xyz_len") != rows["ON"][m].get("xyz_len")]
    print(f"\nbyte-identity over the {len(both)} that finished in both arms:")
    if diff:
        print(f"  ⚠ {len(diff)} differ in generated XYZ length: {diff}")
    else:
        print(
            "  all identical in generated XYZ length -- the bound changed WHICH molecules "
            "finish, not WHAT a finishing molecule produces"
        )

    off_cpu = sum(rows["OFF"][m]["spent_s"] for m in common)
    on_cpu = sum(rows["ON"][m]["spent_s"] for m in common)
    print(
        f"\nCPU over this cohort: OFF {off_cpu / 60:.1f} min -> ON {on_cpu / 60:.1f} min "
        f"({100 * (off_cpu - on_cpu) / off_cpu:+.1f}%)"
    )
    print(
        "  ⚠ the OFF arm ran under an outer kill, so its true cost is a LOWER BOUND and "
        "the real saving is larger than this."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
