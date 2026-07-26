"""Per-molecule transition report for the v0.4.5 re-baseline sweep.

Answers the question the frozen `bucket_report.json` can no longer answer honestly: **what does
today's code actually do**, molecule by molecule, against what the 2026-07-15 capstone snapshot
recorded.

That snapshot is unreliable four independent ways — stale (29 `rmsd_gate` molecules were fixed a
week *after* it was taken), misattributed (`HOCVAY`/`WEFZAL` are generation-side deaths bucketed as
encoder failures), hiding a regression (`XOSTUW_comp_0` passed then and fails now), and
understating `atom_count`. So a single aggregate percentage from it is not trustworthy; only a
**per-molecule transition** is.

The cohort is deliberately two populations, and they answer different questions:

* the **436 gap molecules** — every non-passing row in the capstone. How many now pass?
* a **500-molecule passing guard**, seed 42. Did anything REGRESS? This is the half that a
  gap-only re-run cannot see, and it is where `XOSTUW`-class regressions hide.

Usage:
    PYTHONPATH=src .venv/bin/python tools/rebaseline_report.py \\
        --sweep tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-rebaseline
"""

import argparse
import glob
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

PASSING = {"byte_exact", "key_equal"}
CAPSTONE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..",
    "tmCAT-tmPHOTO_xyz_dataset",
    "results-capstone-v042",
    "bucket_report.json",
)


def load_new(sweep):
    """{molecule: (status, error)} from the sweep's individual reports."""
    out = {}
    for f in glob.glob(os.path.join(sweep, "individual_reports", "*.json")):
        try:
            r = json.load(open(f))
        except Exception:  # noqa: BLE001
            continue
        out[r.get("molecule") or os.path.basename(f)[:-5]] = (r.get("status"), r.get("error") or "")
    return out


def main():
    ap = argparse.ArgumentParser(description="Per-molecule re-baseline transitions.")
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--cohort", default="spec/handoffs/v0.4.5/rebaseline_cohort.json")
    ap.add_argument("--out")
    args = ap.parse_args()

    old_rows = json.load(open(os.path.abspath(CAPSTONE)))
    old = {r["molecule"]: r["bucket"] for r in old_rows}
    new = load_new(args.sweep)
    if not new:
        sys.exit(f"error: no individual_reports under {args.sweep}")

    coh = json.load(open(args.cohort)) if os.path.exists(args.cohort) else {}
    gap_set = set(coh.get("gap", []))
    guard_set = set(coh.get("guard", []))

    # `pending_g-xtb` means the tier-2 g-xTB refinement is QUEUED OR STILL RUNNING -- not that
    # xtb is absent (xtb 6.7.1 is installed; an earlier version of this comment said otherwise
    # and that error was published as a premature pass rate once already). It is NOT a failure
    # and NOT a pass, so it gets its own column. While any row is `pending`, the sweep is still
    # in flight: check `systemctl --user is-active v045-rebaseline`, never the report count.
    def verdict(status):
        if status == "success":
            return "pass"
        if status == "pending_g-xtb":
            return "pending"
        return "fail"

    rows = []
    tally = Counter()
    for m, (status, err) in sorted(new.items()):
        was_pass = old.get(m) in PASSING
        v = verdict(status)
        pop = "gap" if m in gap_set else ("guard" if m in guard_set else "unknown")
        if v == "pending":
            key = f"{pop}:pending"
        elif was_pass and v == "pass":
            key = f"{pop}:still_pass"
        elif was_pass and v == "fail":
            key = f"{pop}:REGRESSED"
        elif not was_pass and v == "pass":
            key = f"{pop}:FIXED"
        else:
            key = f"{pop}:still_fail"
        tally[key] += 1
        rows.append(
            {
                "molecule": m,
                "old_bucket": old.get(m),
                "new_status": status,
                "population": pop,
                "transition": key.split(":")[1],
                "error": err[:160],
            }
        )

    n = len(new)
    print(f"{'=' * 72}\nRE-BASELINE — {n} molecules re-run against today's code\n{'=' * 72}")
    for pop, label in (
        ("gap", "GAP (was failing — how many now pass?)"),
        ("guard", "GUARD (was passing — did anything regress?)"),
        ("unknown", "unclassified"),
    ):
        keys = [k for k in tally if k.startswith(pop + ":")]
        if not keys:
            continue
        tot = sum(tally[k] for k in keys)
        print(f"\n{label}   n={tot}")
        for t in ("FIXED", "REGRESSED", "still_pass", "still_fail", "pending"):
            k = f"{pop}:{t}"
            if tally.get(k):
                mark = "  <<<" if t in ("FIXED", "REGRESSED") else ""
                print(f"    {t:11} {tally[k]:5}  ({100 * tally[k] / tot:5.1f}%){mark}")

    reg = [r for r in rows if r["transition"] == "REGRESSED"]
    fixed = [r for r in rows if r["transition"] == "FIXED"]
    print(f"\n{'-' * 72}")
    print(f"FIXED     : {len(fixed)}")
    print(f"REGRESSED : {len(reg)}   <-- these are the ones that matter most")
    for r in reg[:20]:
        print(f"    {r['molecule']:22} was={r['old_bucket']:16} now={r['new_status']}")
        if r["error"]:
            print(f"        {r['error'][:120]}")
    if len(reg) > 20:
        print(f"    ... and {len(reg) - 20} more")

    pend = sum(v for k, v in tally.items() if k.endswith(":pending"))
    if pend:
        print(f"\nNOTE: {pend} molecules are `pending_g-xtb` — tier-2 g-xTB queued or IN FLIGHT.")
        print("      Counted separately: calling them passes OR failures would be dishonest.")
        print("      The sweep is NOT finished while this is non-zero — poll the systemd unit,")
        print("      not the report count, before quoting any number from this run.")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        json.dump({"n": n, "tally": dict(tally), "rows": rows}, open(args.out, "w"), indent=2)
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

# ---------------------------------------------------------------------------------------
# ⚠ TIMEOUT-BUDGET ARTIFACT — read before calling anything a regression.
#
# This sweep ran at --mol-timeout 300. The capstone baseline it is compared against ran at
# 1800. So a molecule that legitimately needs 300-1800s PASSED then and CANNOT pass now,
# regardless of code quality, and it shows up here as a REGRESSION.
#
# Measured on the first five this report flagged: all 5 had capstone elapsed_s of 303.6,
# 331.5, 344.6, 303.3 and 1117.2 seconds. 5/5 already exceeded 300s while passing. Zero
# were correctness regressions.
#
# This is the same artifact that produced v0.4.4's 11 phantom "regressions" -- all 300s
# timeouts, zero wrong answers -- and the warning is in run_sweep.sh's own header. I wrote
# that warning and then walked into it. So: for every REGRESSED row, check the capstone
# elapsed_s before believing it. A regression whose old elapsed_s exceeds the new budget is
# an artifact, not a defect.
# ---------------------------------------------------------------------------------------
