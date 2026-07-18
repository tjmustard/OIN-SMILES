#!/usr/bin/env python3
"""Compare two round-trip result sets and report per-molecule status transitions.

Joins a *base* (control) and *candidate* result set on the ``molecule`` key and
classifies each shared molecule's transition:

    fail -> pass   (candidate fixed it; a code-attributable fix when base is a matched control)
    pass -> fail   (candidate regressed it)
    fail -> fail   (still failing)
    pass -> pass   (still passing)

Transitions are grouped by each molecule's *original* failure class (optional
``--classes`` map). "pass" means ``status == "success"``; anything else is "fail".

Each of ``--base`` / ``--cand`` may be either:
  * a results directory (uses ``individual_reports/*.json``, else ``summary_roundtrip.json``), or
  * a ``summary_roundtrip.json`` file (a JSON list of report dicts).

Stdlib only — no chemistry stack required.
"""

import argparse
import glob
import json
import os
import sys
from collections import defaultdict


def load_statuses(path):
    """Return {molecule: status} from a results dir or a summary json file."""
    reports = []
    if os.path.isdir(path):
        indiv = os.path.join(path, "individual_reports")
        if os.path.isdir(indiv):
            for fn in glob.glob(os.path.join(indiv, "*.json")):
                try:
                    reports.append(json.load(open(fn)))
                except Exception as e:  # noqa: BLE001
                    print(f"warning: skipping unreadable {fn}: {e}", file=sys.stderr)
        else:
            summ = os.path.join(path, "summary_roundtrip.json")
            if not os.path.exists(summ):
                sys.exit(
                    f"error: {path} has neither individual_reports/ nor summary_roundtrip.json"
                )
            reports = json.load(open(summ))
    elif os.path.isfile(path):
        reports = json.load(open(path))
    else:
        sys.exit(f"error: {path} is not a file or directory")

    out = {}
    for r in reports:
        mol = r.get("molecule")
        if mol is None:
            continue
        out[mol] = "pass" if r.get("status") == "success" else "fail"
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--base", required=True, help="Control result set (dir or summary json)")
    ap.add_argument("--cand", required=True, help="Candidate result set (dir or summary json)")
    ap.add_argument("--classes", help="Optional molecule->class JSON map for grouping")
    ap.add_argument("--out", help="Optional path to write the full result as JSON")
    args = ap.parse_args()

    base = load_statuses(args.base)
    cand = load_statuses(args.cand)
    classes = {}
    if args.classes:
        classes = json.load(open(args.classes))

    shared = sorted(set(base) & set(cand))
    only_base = set(base) - set(cand)
    only_cand = set(cand) - set(base)
    if only_base or only_cand:
        print(
            f"warning: base has {len(only_base)} molecules absent from cand, "
            f"cand has {len(only_cand)} absent from base; comparing the "
            f"{len(shared)}-molecule intersection.",
            file=sys.stderr,
        )

    # per-class transition tallies
    TRANSITIONS = ("fail->pass", "pass->fail", "fail->fail", "pass->pass")
    per_class = defaultdict(lambda: dict.fromkeys(TRANSITIONS, 0))
    fixes, regressions = [], []
    for mol in shared:
        t = f"{base[mol]}->{cand[mol]}"
        cls = classes.get(mol, "unknown")
        per_class[cls][t] += 1
        if t == "fail->pass":
            fixes.append((mol, cls))
        elif t == "pass->fail":
            regressions.append((mol, cls))

    totals = dict.fromkeys(TRANSITIONS, 0)
    for row in per_class.values():
        for t in TRANSITIONS:
            totals[t] += row[t]

    # ---- report ----
    width = max((len(c) for c in per_class), default=5)
    hdr = f"{'class':<{width}}  {'fix':>5} {'regr':>5} {'stillF':>7} {'stillP':>7}"
    print(hdr)
    print("-" * len(hdr))
    for cls, row in sorted(per_class.items(), key=lambda kv: (-kv[1]["fail->pass"], kv[0])):
        print(
            f"{cls:<{width}}  {row['fail->pass']:>5} {row['pass->fail']:>5} "
            f"{row['fail->fail']:>7} {row['pass->pass']:>7}"
        )
    print("-" * len(hdr))
    print(
        f"{'TOTAL':<{width}}  {totals['fail->pass']:>5} {totals['pass->fail']:>5} "
        f"{totals['fail->fail']:>7} {totals['pass->pass']:>7}"
    )

    print(f"\nproven fixes (fail->pass): {len(fixes)}")
    for mol, cls in sorted(fixes):
        print(f"  + {mol}  [{cls}]")
    print(f"\nregressions (pass->fail): {len(regressions)}")
    for mol, cls in sorted(regressions):
        print(f"  - {mol}  [{cls}]")

    print(
        f"\nHEADLINE: fixes = {totals['fail->pass']} "
        f"(of {len(shared)} compared), regressions = {totals['pass->fail']}"
    )

    if args.out:
        json.dump(
            {
                "n_compared": len(shared),
                "only_base": sorted(only_base),
                "only_cand": sorted(only_cand),
                "totals": totals,
                "per_class": per_class,
                "fixes": fixes,
                "regressions": regressions,
            },
            open(args.out, "w"),
            indent=2,
        )
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
