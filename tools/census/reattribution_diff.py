"""Two attribution tables, one question: where did each molecule's fault GO? Reads only.

    python tools/census/reattribution_diff.py --old <census table.tsv> --new <new table.tsv> \
        --movers <union movers.txt> [--out report.json]

The census (``attribution_table.py``) gives every molecule exactly one fault. After a shipped
default changes, the table is stale for every molecule the change touched, and the two headline
gaps -- self-consistent and VERIFIED -- still lean on it. This prints:

  1. the new FAIL partition beside the old, in molecules and in points;
  2. the VERIFIED gap: FAIL + the passes that carry a fault (metric false passes);
  3. the transition matrix old fault -> new fault, SPLIT by mover / non-mover.

WHY THE SPLIT IS THE CONTROL. A non-mover gets the same two strings from the encoder in both
sweeps, so any fault it changes is NOT the lever: it is the generator's run-to-run variation, or an
instrument reading a boundary case differently. That number is the noise floor of attribution
itself, and a lane smaller than it cannot be chosen from this table. A broken join (wrong key,
wrong table) shows up here first -- as thousands of non-movers changing fault.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

PT = 100.0 / 5000


def _table(path: Path):
    with open(path, newline="") as fh:
        rows = {r["molecule"]: r for r in csv.DictReader(fh, delimiter="\t")}
    if len(rows) != 5000:
        sys.exit(f"ABORT: {path} has {len(rows)} rows, not 5000")
    return rows


def _owner(fault: str) -> str:
    return "NONE" if fault == "NONE" else fault.split("_")[0]


def _partition(rows, outcome):
    return Counter(r["fault"] for r in rows.values() if r["outcome"] == outcome)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--old", type=Path, required=True)
    ap.add_argument("--new", type=Path, required=True)
    ap.add_argument("--movers", type=Path, required=True, help="one molecule per line")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    old, new = _table(args.old), _table(args.new)
    if set(old) != set(new):
        sys.exit("ABORT: the two tables do not hold the same 5,000 molecules")
    movers = {ln.strip().removesuffix(".xyz") for ln in args.movers.read_text().split()}
    if not movers or not movers <= set(old):
        sys.exit(f"ABORT: mover list has {len(movers)} names, {len(movers - set(old))} unknown")

    report = {"n": len(old), "n_movers": len(movers)}
    print(
        f"DENOMINATOR {len(old)} molecules   movers {len(movers)}   non-movers {len(old) - len(movers)}"
    )

    # ---- 1. FAIL partition ----------------------------------------------------------------
    fo, fn = _partition(old, "FAIL"), _partition(new, "FAIL")
    print(
        f"\n1. FAIL partition    old {sum(fo.values())} = {sum(fo.values()) * PT:.2f} pts"
        f"    new {sum(fn.values())} = {sum(fn.values()) * PT:.2f} pts"
    )
    print(f"   {'fault':18s} {'old':>5s} {'pts':>6s}   {'new':>5s} {'pts':>6s}   {'delta':>6s}")
    for f in sorted(set(fo) | set(fn), key=lambda k: -fn.get(k, 0)):
        print(
            f"   {f:18s} {fo.get(f, 0):5d} {fo.get(f, 0) * PT:6.2f}   {fn.get(f, 0):5d} "
            f"{fn.get(f, 0) * PT:6.2f}   {fn.get(f, 0) - fo.get(f, 0):+6d}"
        )
    oo, on = Counter(), Counter()
    for f, c in fo.items():
        oo[_owner(f)] += c
    for f, c in fn.items():
        on[_owner(f)] += c
    print(
        "   by owner:  "
        + "   ".join(
            f"{k} {oo.get(k, 0)}->{on.get(k, 0)} ({on.get(k, 0) * PT:.2f} pts)"
            for k in sorted(set(oo) | set(on), key=lambda k: -on.get(k, 0))
        )
    )
    report["fail_old"], report["fail_new"] = dict(fo), dict(fn)
    report["fail_by_owner_old"], report["fail_by_owner_new"] = dict(oo), dict(on)

    # ---- 2. the VERIFIED gap ---------------------------------------------------------------
    po, pn = _partition(old, "PASS"), _partition(new, "PASS")
    fpo = {f: c for f, c in po.items() if f != "NONE"}
    fpn = {f: c for f, c in pn.items() if f != "NONE"}
    print(
        f"\n2. passes that carry a fault (metric false passes)   old {sum(fpo.values())}   "
        f"new {sum(fpn.values())}"
    )
    for f in sorted(set(fpo) | set(fpn), key=lambda k: -fpn.get(k, 0)):
        print(
            f"   {f:18s} {fpo.get(f, 0):5d}   {fpn.get(f, 0):5d}   {fpn.get(f, 0) - fpo.get(f, 0):+6d}"
        )
    for tag, p, f in (("old", po, fo), ("new", pn, fn)):
        n_pass, n_ver = sum(p.values()), p.get("NONE", 0)
        print(
            f"   {tag}: self-consistent {n_pass} = {n_pass * PT:.2f}%   VERIFIED {n_ver} = "
            f"{n_ver * PT:.2f}%   gap to 100: {sum(f.values()) * PT:.2f} / {(5000 - n_ver) * PT:.2f}"
        )
    vg = Counter()
    for r in new.values():
        if r["fault"] != "NONE":
            vg[_owner(r["fault"])] += 1
    print(
        "   NEW verified gap by owner (FAIL + false pass): "
        + "   ".join(f"{k} {c} ({c * PT:.2f})" for k, c in vg.most_common())
    )
    report["false_pass_old"], report["false_pass_new"] = fpo, fpn
    report["verified_gap_by_owner_new"] = dict(vg)

    # ---- 3. transitions, split by mover ----------------------------------------------------
    print("\n3. transitions old fault -> new fault")
    for name, keep in (
        ("MOVERS", lambda m: m in movers),
        ("NON-MOVERS", lambda m: m not in movers),
    ):
        t = Counter((old[m]["fault"], new[m]["fault"]) for m in old if keep(m))
        n = sum(t.values())
        same = sum(c for (a, b), c in t.items() if a == b)
        print(
            f"   {name}: {n} molecules, {same} keep their fault, {n - same} change "
            f"({100 * (n - same) / n:.1f}%)"
        )
        for (a, b), c in sorted(
            ((k, c) for k, c in t.items() if k[0] != k[1]), key=lambda x: -x[1]
        )[:14]:
            print(f"      {c:4d}  {a:18s} -> {b}")
        oc = Counter((old[m]["outcome"], new[m]["outcome"]) for m in old if keep(m))
        print(
            f"      outcome: PASS->FAIL {oc.get(('PASS', 'FAIL'), 0)}   FAIL->PASS "
            f"{oc.get(('FAIL', 'PASS'), 0)}"
        )
        report[f"transitions_{name.lower().replace('-', '_')}"] = {
            f"{a} -> {b}": c for (a, b), c in t.items() if a != b
        }
        report[f"changed_{name.lower().replace('-', '_')}"] = n - same

    if args.out:
        args.out.write_text(json.dumps(report, indent=1, sort_keys=True))
        print(f"\nwrote {args.out.name}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
