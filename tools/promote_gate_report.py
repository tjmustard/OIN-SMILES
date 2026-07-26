#!/usr/bin/env python
"""Render the four `OIN_ACCEPT_SCORED` promotion gates from `ab_accept_scored.py` output.

Reads one or more `--out` JSONs and prints the tables that go in
`docs/ACCEPT_SCORED_v0.4.7.md`. Kept separate from the measuring script so the tables can be
regenerated without re-running anything -- at load 39 a re-run is not free and is not
reproducible in wall-clock terms anyway.

Every number here is load-independent by construction: sha equality, class counts, clash
counts, overlap ratios. Seconds are printed only where labelled ADVISORY.

Usage:
    python tools/promote_gate_report.py --ab spec/handoffs/v0.4.7/runs/ab22.json \
        [--control spec/handoffs/v0.4.7/runs/aprime.json]
"""

from __future__ import annotations

import argparse
import json
import statistics


def _idx(rows: list[dict]) -> dict[str, dict]:
    return {r["molecule"]: r for r in rows}


def _cell(v) -> str:
    return "-" if v is None else str(v)


def gate3(a: list[dict], b: list[dict], control: list[dict] | None) -> None:
    """G3 -- byte identity of `smiles_2`, the wave's actual contract."""
    A, B = _idx(a), _idx(b)
    C = _idx(control) if control else {}
    print("## G3 - sha256(smiles_2) byte identity (load-independent)\n")
    hdr = "| molecule | sha_A(smiles_2) | sha_B(smiles_2) | equal? | sha_in control |"
    if C:
        hdr += " sha_A'(smiles_2) | A==A'? |"
    print(hdr)
    print("|---" * (hdr.count("|") - 1) + "|")
    same = diff = one = neither = 0
    moved = []
    for m in sorted(set(A) | set(B)):
        ra, rb = A.get(m, {}), B.get(m, {})
        sa, sb = ra.get("sha_out"), rb.get("sha_out")
        if sa and sb:
            eq = "YES" if sa == sb else "**NO**"
            same += sa == sb
            diff += sa != sb
        elif sa or sb:
            eq = "one-arm"
            one += 1
        else:
            eq = "neither"
            neither += 1
        ctrl = "ok" if ra.get("sha_in") == rb.get("sha_in") else "**MOVED**"
        if ra.get("sha_in") != rb.get("sha_in"):
            moved.append(m)
        row = f"| {m} | {_cell(sa)} | {_cell(sb)} | {eq} | {ctrl} |"
        if C:
            sc = C.get(m, {}).get("sha_out")
            aeq = "-" if not (sa and sc) else ("YES" if sa == sc else "**NO**")
            row += f" {_cell(sc)} | {aeq} |"
        print(row)
    print(
        f"\n**byte-identical {same} | DIFFERENT {diff} | only one arm produced a string {one} "
        f"| neither {neither}**"
    )
    if moved:
        print(f"\n*** CONTROL VIOLATION: sha_in moved on {moved} -- run is confounded ***")


def gate12(a: list[dict], b: list[dict]) -> None:
    """G1 (quality) and G2 (independent re-perception) side by side."""
    print("\n## G1 - structure quality, and G2 - the dropped independent test\n")
    print(
        "| arm | pass (circular) | indep pass / measured | credited by circular ONLY | "
        "clash_vdw total | mols with any clash | severe | worst_overlap min | "
        "worst_overlap median | median s (ADVISORY) |"
    )
    print("|---" * 10 + "|")
    for rows, label in ((a, "A default"), (b, "B scored")):
        cl = [r["clash_vdw"] for r in rows if isinstance(r.get("clash_vdw"), int)]
        sev = [r["clash_severe"] for r in rows if isinstance(r.get("clash_severe"), int)]
        wo = [r["worst_overlap"] for r in rows if isinstance(r.get("worst_overlap"), float)]
        el = [r["elapsed_s"] for r in rows if isinstance(r.get("elapsed_s"), (int, float))]
        im = sum(1 for r in rows if "indep_passed" in r)
        ip = sum(1 for r in rows if r.get("indep_passed"))
        co = sum(1 for r in rows if r.get("passed") and not r.get("indep_passed"))
        print(
            f"| {label} | {sum(1 for r in rows if r.get('passed'))}/{len(rows)} | {ip}/{im} | "
            f"{co} | {sum(cl)} | {sum(1 for v in cl if v > 0)}/{len(cl)} | {sum(sev)} | "
            f"{round(min(wo), 4) if wo else '-'} | "
            f"{round(statistics.median(wo), 4) if wo else '-'} | "
            f"{round(statistics.median(el), 2) if el else '-'} |"
        )

    A, B = _idx(a), _idx(b)
    print("\n### Per-molecule quality delta (only where both arms produced coordinates)\n")
    print("| molecule | clash_vdw A | clash_vdw B | worst_overlap A | worst_overlap B | verdict |")
    print("|---" * 6 + "|")
    worse = better = tie = 0
    for m in sorted(set(A) & set(B)):
        ra, rb = A[m], B[m]
        if not isinstance(ra.get("clash_vdw"), int) or not isinstance(rb.get("clash_vdw"), int):
            continue
        wa, wb = ra["worst_overlap"], rb["worst_overlap"]
        if rb["clash_vdw"] > ra["clash_vdw"] or wb < wa - 1e-9:
            v = "B worse"
            worse += 1
        elif rb["clash_vdw"] < ra["clash_vdw"] or wb > wa + 1e-9:
            v = "B better"
            better += 1
        else:
            v = "identical"
            tie += 1
        print(f"| {m} | {ra['clash_vdw']} | {rb['clash_vdw']} | {wa} | {wb} | {v} |")
    print(f"\n**B worse {worse} | B better {better} | identical {tie}**")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ab", required=True)
    ap.add_argument("--control", help="--single-arm JSON re-running arm A (determinism check)")
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    d = json.load(open(args.ab))
    a, b = d["arm_a"], d["arm_b"]
    c = json.load(open(args.control))["rows"] if args.control else None
    if args.title:
        print(f"# {args.title}\n")
    gate3(a, b, c)
    gate12(a, b)
    print(f"\nPASS REGRESSIONS (A pass -> B fail): {d.get('regressions') or 'none'}")
    print(f"PASS FIXES       (A fail -> B pass): {d.get('fixes') or 'none'}")
    print(f"INDEP REGRESSIONS (A -> B): {d.get('indep_regressions') or 'none'}")
    print(f"INDEP FIXES       (A -> B): {d.get('indep_fixes') or 'none'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
