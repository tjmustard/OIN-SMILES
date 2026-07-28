"""Freeze v0.4.13's authoritative table from the offline veto arm (v0.4.13).

WHAT THIS PRODUCES, AND WHAT IT IS NOT
======================================
`results-v0.4.13-honest/` is a **derived table**, not a sweep. It is derived exactly the way
`results-v0.4.8-honest` was derived from `results-v0.4.6-sweep`: by re-scoring stored artifacts
under a new configuration, with the generator NOT re-run.

That is legitimate here for a reason that was measured rather than assumed
(`tools/fold_key_invariance.py`): the promotion is **generator-neutral**. Over all 9669 strings
in the frozen corpus the fold moves 1019 of them and changes **0** comparison keys, and
`accept_fn` decides by key -- so the generator returns bit-identical conformers with the levers
on or off. A fresh 55 CPU-h sweep would re-run a stochastic generator and contaminate the
measurement with run-to-run variation; it would be worse evidence, not better.

The bucket counts come from `fold_transition_sim.py --arm veto`, which classified every mover
with the same `bucket_of` logic the authoritative report uses, over a complete denominator:
393 movers, **0 excluded for drift, 0 excluded as unavailable**.

Usage
-----
    PYTHONPATH=src .venv/bin/python tools/freeze_v0413_table.py \
        --base <main>/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest \
        --transition docs/agentic-notes/v0.4.13/fold_transition_veto.json \
        --out <main>/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.13-honest
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ORDER = [
    "byte_exact",
    "key_equal",
    "facmer_divergent",
    "structural",
    "hard_fail",
    "encode_fail",
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--transition", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--extra", type=Path, action="append", default=[])
    args = ap.parse_args()

    tr = json.loads(args.transition.read_text())
    if tr.get("arm") != "veto":
        sys.exit(f"FATAL: transition file is arm={tr.get('arm')!r}, expected 'veto'")
    if tr.get("excluded_unavailable") or tr.get("excluded_drift"):
        sys.exit(
            "FATAL: the veto arm excluded molecules "
            f"({len(tr.get('excluded_unavailable', []))} unavailable, "
            f"{len(tr.get('excluded_drift', []))} drift). A table frozen from a partial arm "
            "understates the veto's cost and must not be published as authoritative."
        )

    after, n = tr["after"], tr["n"]
    # Collapse the key_equal sub-split back into a top-level bucket for the headline table.
    buckets, sub = {}, {}
    for k, v in after.items():
        if "/" in k:
            top, s = k.split("/", 1)
            buckets[top] = buckets.get(top, 0) + v
            sub[s] = v
        else:
            buckets[k] = buckets.get(k, 0) + v
    if sum(buckets.values()) != n:
        sys.exit(f"FATAL: buckets sum to {sum(buckets.values())}, expected {n}")

    args.out.mkdir(parents=True, exist_ok=True)
    be = buckets.get("byte_exact", 0)

    lines = [
        "# v0.4.13 authoritative table — DERIVED, not swept",
        "",
        f"Derived from `{args.base.name}` by `tools/fold_transition_sim.py --arm veto`.",
        "**The generator was NOT re-run** — see `FROZEN.md` for why that is exact here.",
        "",
        "## Buckets",
        "",
        "| bucket | count | % |",
        "|---|---:|---:|",
    ]
    for b in ORDER:
        if b in buckets:
            lines.append(f"| {b} | {buckets[b]} | {100 * buckets[b] / n:.2f}% |")
    lines += [
        f"| **total** | **{n}** | **100.00%** |",
        "",
        "### key_equal sub-split",
        "",
        "| subclass | count |",
        "|---|---:|",
    ]
    for s, v in sorted(sub.items()):
        lines.append(f"| {s} | {v} |")
    lines += [
        "",
        "## Movement vs the v0.4.8 baseline",
        "",
        "| | |",
        "|---|---|",
        f"| `byte_exact` before | {tr['byte_exact_before']} ({100 * tr['byte_exact_before'] / n:.2f}%) |",
        f"| `byte_exact` after | **{be} ({100 * be / n:.2f}%)** |",
        f"| points | **{tr['points']:+.2f}** |",
        f"| moved in a BAD direction | **{tr['bad_direction']}** |",
        f"| excluded (drift / unavailable) | "
        f"**{len(tr.get('excluded_drift', []))} / {len(tr.get('excluded_unavailable', []))}** |",
        "",
    ]
    for k, v in tr["transitions"].items():
        lines.append(f"- `{k}` — **{v}**")
    (args.out / "bucket_report_PASS1_authoritative.md").write_text("\n".join(lines) + "\n")

    json.dump(
        {
            "derived_from": str(args.base),
            "method": "fold_transition_sim.py --arm veto (generator NOT re-run)",
            "n": n,
            "buckets": buckets,
            "key_equal_subsplit": sub,
            "byte_exact_before": tr["byte_exact_before"],
            "byte_exact_after": be,
            "points": tr["points"],
            "bad_direction": tr["bad_direction"],
            "excluded_drift": len(tr.get("excluded_drift", [])),
            "excluded_unavailable": len(tr.get("excluded_unavailable", [])),
        },
        open(args.out / "bucket_report_PASS1_authoritative.json", "w"),
        indent=2,
    )

    (args.out / "SOURCE").write_text(str(args.base) + "\n")
    for extra in args.extra:
        if extra.exists():
            shutil.copy2(extra, args.out / extra.name)

    print(f"froze {args.out}")
    print(f"  byte_exact {tr['byte_exact_before']} -> {be}  ({tr['points']:+.2f} pts)")
    print(f"  bad_direction {tr['bad_direction']}, excluded 0/0")
    print(f"\n#DONE {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
