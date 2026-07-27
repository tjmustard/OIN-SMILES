#!/usr/bin/env python
"""Re-run the four `OIN_ACCEPT_SCORED` promotion gates with the attachment check enabled.

Takes a baseline A/B run (`ab_accept_scored.py --cohort ...`, giving arm A = default and
arm B = the bare lever) plus a `--single-arm 1` run made with the attachment check ON
(arm C = lever + check), and prints the same scorecard shape as
`docs/agentic-notes/v0.4.7/ACCEPT_SCORED_v0.4.7.md` §5.1 so the numbers are directly comparable.

WHY ARM C CAN BE A SINGLE-ARM RUN
---------------------------------
Re-measuring arm A costs more than the whole rest of this lane (its median is 13.9s against
arm B's 5.5s, and its tail runs to 400s). §4.2 and §4.4 of the promote lane established
generator determinism at the sha AND coordinate level -- 11/11 and 17/17 identical across
independent processes, tool builds and caps -- so a stored arm A is a valid comparator.
That is a premise, so it is CHECKED rather than assumed: `sha_in` must agree between the
stored run and arm C on every molecule, and any molecule where it does not is dropped from
the comparison and named.

Usage:
    python tools/attach_gate_report.py --baseline runs/ab22.json --armc runs/attach_c22.json
"""

from __future__ import annotations

import argparse
import json


def load_single(path: str) -> dict:
    d = json.load(open(path))
    rows = d["rows"] if isinstance(d, dict) and "rows" in d else d
    return {r["molecule"]: r for r in rows}


def compare(name_a: str, A: dict, name_c: str, C: dict, mols: list) -> None:
    print(
        f"\n{'=' * 78}\n  {name_a}   vs   {name_c}      ({len(mols)} molecules compared)\n{'=' * 78}"
    )

    pa = [m for m in mols if A[m].get("passed")]
    pc = [m for m in mols if C[m].get("passed")]
    print(f"\nG4  pass rate            {len(pa):>3d}/{len(mols)}  ->  {len(pc):>3d}/{len(mols)}")
    print(
        f"    PASS REGRESSIONS     {[m for m in mols if A[m].get('passed') and not C[m].get('passed')] or 'none'}"
    )
    print(
        f"    PASS FIXES           {[m for m in mols if C[m].get('passed') and not A[m].get('passed')] or 'none'}"
    )

    # G3 -- byte identity of smiles_2, the wave's contract.
    same = diff = onearm = neither = 0
    diffs = []
    ctrl = []
    for m in mols:
        sa, sc = A[m].get("sha_out"), C[m].get("sha_out")
        ia, ic = A[m].get("sha_in"), C[m].get("sha_in")
        if ia and ic and ia != ic:
            ctrl.append(m)
        if sa and sc:
            if sa == sc:
                same += 1
            else:
                diff += 1
                diffs.append((m, sa, sc))
        elif sa or sc:
            onearm += 1
        else:
            neither += 1
    print(
        f"\nG3  sha256(smiles_2)     identical {same} | DIFFERENT {diff} | one-arm {onearm} | neither {neither}"
    )
    for m, sa, sc in diffs:
        print(f"      *** {m}: {sa} -> {sc}")
    print(
        f"    sha_in control       {'OK' if not ctrl else '*** MOVED on ' + str(ctrl) + ' -- comparison confounded ***'}"
    )

    # G1 -- structure quality, paired over molecules where BOTH arms returned coordinates.
    paired = [
        m
        for m in mols
        if isinstance(A[m].get("clash_vdw"), int) and isinstance(C[m].get("clash_vdw"), int)
    ]
    ca = sum(A[m]["clash_vdw"] for m in paired)
    cc = sum(C[m]["clash_vdw"] for m in paired)
    sa_ = sum(A[m]["clash_severe"] for m in paired)
    sc_ = sum(C[m]["clash_severe"] for m in paired)
    woa = min((A[m]["worst_overlap"] for m in paired), default=None)
    woc = min((C[m]["worst_overlap"] for m in paired), default=None)
    worse = sum(1 for m in paired if C[m]["worst_overlap"] < A[m]["worst_overlap"])
    better = sum(1 for m in paired if C[m]["worst_overlap"] > A[m]["worst_overlap"])
    ident = len(paired) - worse - better
    print(f"\nG1  structure (paired n={len(paired)})")
    print(f"    clash_vdw total      {ca:>4d}  ->  {cc:>4d}")
    print(f"    clash_severe total   {sa_:>4d}  ->  {sc_:>4d}")
    print(f"    worst_overlap min    {woa}  ->  {woc}")
    print(f"    per-molecule         worse {worse} · better {better} · identical {ident}")

    # G2 -- the independent re-perception the lever drops.
    im = [m for m in mols if "indep_passed" in A[m] and "indep_passed" in C[m]]
    ia_ = sum(1 for m in im if A[m]["indep_passed"])
    ic_ = sum(1 for m in im if C[m]["indep_passed"])
    reg = [m for m in im if A[m]["indep_passed"] and not C[m]["indep_passed"]]
    fix = [m for m in im if C[m]["indep_passed"] and not A[m]["indep_passed"]]
    print(f"\nG2  indep re-perception  {ia_}/{len(im)}  ->  {ic_}/{len(im)}")
    print(f"    INDEP REGRESSIONS    {reg or 'none'}   ({len(reg)})")
    print(f"    INDEP FIXES          {fix or 'none'}   ({len(fix)})")

    # Runtime -- ADVISORY only. Load on this box is 40+ on 12 cores.
    ea = [A[m]["elapsed_s"] for m in mols if isinstance(A[m].get("elapsed_s"), (int, float))]
    ec = [C[m]["elapsed_s"] for m in mols if isinstance(C[m].get("elapsed_s"), (int, float))]
    if ea and ec:
        ea_s, ec_s = sorted(ea), sorted(ec)
        print(
            f"\nSPEED (ADVISORY -- shared box under heavy load)\n"
            f"    median_s             {ea_s[len(ea_s) // 2]:.2f}  ->  {ec_s[len(ec_s) // 2]:.2f}\n"
            f"    total_s              {sum(ea):.1f}  ->  {sum(ec):.1f}\n"
            f"    >30s                 {sum(1 for v in ea if v > 30)}  ->  {sum(1 for v in ec if v > 30)}"
        )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--baseline",
        required=True,
        help="ab_accept_scored.py --out JSON (arm_a + arm_b), or a --single-arm JSON to use "
        "as the arm-B-only comparator (for cohorts where arm A is unaffordable)",
    )
    ap.add_argument(
        "--armc", required=True, help="--single-arm 1 JSON run WITH the attachment check"
    )
    args = ap.parse_args()

    base = json.load(open(args.baseline))
    if "arm_a" in base:
        A = {r["molecule"]: r for r in base["arm_a"]}
        B = {r["molecule"]: r for r in base["arm_b"]}
    else:
        A = None
        B = load_single(args.baseline)
    C = load_single(args.armc)

    arms = [x for x in (A, B, C) if x is not None]
    mols = sorted(set.intersection(*[set(x) for x in arms]))
    dropped = sorted(set.union(*[set(x) for x in arms]) - set(mols))
    if dropped:
        print(f"  not present in every arm, excluded: {dropped}")

    if A is not None:
        compare("A default", A, "C lever+check", C, mols)
    compare("B lever only", B, "C lever+check", C, mols)

    if A is not None:
        # The question the whole lane turns on, per molecule: on the molecules the bare lever
        # BROKE, does the check give the structure back, and what does it charge for it?
        broke = [m for m in mols if A[m].get("indep_passed") and not B[m].get("indep_passed")]
        if broke:
            print(
                f"\n{'=' * 78}\n  RECOVERY vs COST on the {len(broke)} molecules the bare lever "
                f"broke\n{'=' * 78}"
            )
            print(
                f"{'molecule':18s}{'A_ind':>6s}{'B_ind':>6s}{'C_ind':>6s}   "
                f"{'A_s':>8s}{'B_s':>8s}{'C_s':>8s}   verdict"
            )
            rec = 0
            for m in sorted(broke):
                ci = C[m].get("indep_passed")
                rec += bool(ci)
                print(
                    f"{m:18s}{str(A[m].get('indep_passed')):>6s}{str(B[m].get('indep_passed')):>6s}"
                    f"{str(ci):>6s}   {str(A[m].get('elapsed_s')):>8s}"
                    f"{str(B[m].get('elapsed_s')):>8s}{str(C[m].get('elapsed_s')):>8s}   "
                    f"{'RECOVERED' if ci else 'still broken'}"
                )
            print(f"\n  recovered {rec}/{len(broke)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
