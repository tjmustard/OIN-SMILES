"""Six-bucket round-trip classifier for the v0.4.4 fac/mer-aware canonical key (SL0).

Reads ``individual_reports/*.json`` from a results directory and classifies every
molecule with the NEW ``oin.compare`` key, so it can measure the accuracy the fac/mer
key reclaims (benign canonicalization the old key already collapsed but that reads as a
byte mismatch) AND the accuracy it newly catches (fac/mer superpositions the old key
waved through). This is the measurement instrument every v0.4.4 swimlane gates against.

The six mutually-exclusive buckets (first match wins):

1. ``encode_fail``    -- ``smiles_1`` is null (encoding the INPUT XYZ failed; SL5 worklist).
2. ``hard_fail``      -- ``status != success`` and NOT a string mismatch (no-conformer /
                         timeout / RMSD / crash; SL4 worklist).
3. ``byte_exact``     -- ``smiles_1 == smiles_2``.
4. ``key_equal``      -- not byte-exact but the new key agrees (benign canonicalization --
                         the win reclaimed). Sub-split ``slot_renumber`` / ``fragment_reorder``
                         / ``winding_star_drift`` / ``rdkit_canonical`` for diagnostics.
5. ``facmer_divergent`` -- key differs, but same metal and same ligand-body multiset, so only
                         the donor ARRANGEMENT differs (fac/mer, cis/trans, winding face).
                         The honest cost of the fac/mer key: reported explicitly.
6. ``structural``     -- key differs and the metal or ligand-body multiset differs
                         (bond-order / aromaticity / connectivity error).

Buckets 3-6 cover every row with both strings present (all successes plus the stored
``String mismatch`` failures); a string-mismatch that the new key now collapses lands in
``key_equal`` -- a reclaimed false failure.

Usage:
    PYTHONPATH=src python tools/roundtrip_bucket_report.py \
        --results-dir tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042
"""

import argparse
import glob
import json
import os
import re
import sys
from collections import Counter
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from oinsmiles.oin.compare import (  # noqa: E402
    canonical_roundtrip_key,
    normalize_oin_for_comparison,
)

_ETA_RE = re.compile(r"\{\d+[<>]\}")
_SLOT_RE = re.compile(r"\{(\d+)([<>^]?)\}")


def _is_eta(smiles):
    return bool(_ETA_RE.search(smiles or ""))


def _key(smiles):
    """New fac/mer-aware key, or None if it cannot be computed."""
    try:
        return canonical_roundtrip_key(smiles)
    except Exception:
        return None


def _key_equal_subclass(s1, s2):
    """Diagnostic sub-label for a byte-different but key-equal (benign) pair."""
    n1 = normalize_oin_for_comparison((s1 or "").strip())
    n2 = normalize_oin_for_comparison((s2 or "").strip())
    frags1 = [f for f in n1.split(".") if f]
    frags2 = [f for f in n2.split(".") if f]
    if sorted(frags1) == sorted(frags2) and frags1 != frags2:
        return "fragment_reorder"
    blank_num = lambda s: _SLOT_RE.sub(r"{\2}", s)  # blank slot number, keep winding  # noqa: E731
    if blank_num(n1) == blank_num(n2):
        return "slot_renumber"
    strip_wind = lambda s: re.sub(r"\{(\d+)[<>^]\}", r"{\1}", s)  # noqa: E731
    if strip_wind(n1) == strip_wind(n2):
        return "winding_star_drift"
    return "rdkit_canonical"


def classify(rep):
    """Return (bucket, subclass_or_None) for one report dict."""
    s1 = rep.get("smiles_1")
    s2 = rep.get("smiles_2")
    status = rep.get("status")
    error = rep.get("error") or ""

    if s1 is None:
        return "encode_fail", None
    if status != "success" and not error.startswith("String mismatch"):
        # no-conformer / timeout / high-RMSD / crash -- string round-trip not the issue here
        return "hard_fail", None
    if s2 is None:
        # defensive: a string-mismatch row should carry both, but guard anyway
        return "hard_fail", None

    if s1 == s2:
        return "byte_exact", None

    k1, k2 = _key(s1), _key(s2)
    if k1 is None or k2 is None:
        return "structural", None
    if k1 == k2:
        return "key_equal", _key_equal_subclass(s1, s2)

    # key differs: arrangement-only (same metal + same ligand bodies) vs structural
    metal1, _sig1, bodies1 = k1
    metal2, _sig2, bodies2 = k2
    if metal1 == metal2 and bodies1 == bodies2:
        return "facmer_divergent", None
    return "structural", None


def _percentile(values, q):
    if not values:
        return None
    xs = sorted(values)
    if len(xs) == 1:
        return xs[0]
    pos = (len(xs) - 1) * q / 100.0
    lo = int(pos)
    frac = pos - lo
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * frac


def _pct_row(label, values):
    if not values:
        return f"| {label} | 0 | - | - | - | - | - |"
    return (
        f"| {label} | {len(values)} | "
        f"{_percentile(values, 50):.1f} | {_percentile(values, 90):.1f} | "
        f"{_percentile(values, 95):.1f} | {_percentile(values, 99):.1f} | "
        f"{max(values):.1f} |"
    )


def _load_reports(indiv_dir, only, sample, shard, limit):
    paths = sorted(glob.glob(os.path.join(indiv_dir, "*.json")))
    if shard:
        i, n = shard
        paths = [p for idx, p in enumerate(paths) if idx % n == i]
    rows = []
    for fp in paths:
        try:
            with open(fp) as f:
                rep = json.load(f)
        except Exception:
            continue
        if only and rep.get("molecule") not in only:
            continue
        rows.append(rep)
    if sample:
        rows = rows[:: max(1, len(rows) // sample)][:sample]
    if limit:
        rows = rows[:limit]
    return rows


def build_report(rows, indiv_dir):
    results = []
    for rep in rows:
        bucket, subclass = classify(rep)
        results.append(
            {
                "molecule": rep.get("molecule"),
                "bucket": bucket,
                "subclass": subclass,
                "status": rep.get("status"),
                "eta": _is_eta(rep.get("smiles_1")),
                "elapsed_s": (rep.get("metrics") or {}).get("elapsed_s"),
                "smiles_1": rep.get("smiles_1"),
                "smiles_2": rep.get("smiles_2"),
                "error": rep.get("error"),
            }
        )

    total = len(results)
    counts = Counter(r["bucket"] for r in results)
    subcounts = Counter(r["subclass"] for r in results if r["bucket"] == "key_equal")
    order = [
        "byte_exact",
        "key_equal",
        "facmer_divergent",
        "structural",
        "hard_fail",
        "encode_fail",
    ]

    elapsed_all = [r["elapsed_s"] for r in results if isinstance(r["elapsed_s"], (int, float))]
    elapsed_eta = [
        r["elapsed_s"] for r in results if r["eta"] and isinstance(r["elapsed_s"], (int, float))
    ]

    lines = [
        "# v0.4.4 Round-Trip Bucket Report",
        "",
        f"Generated {datetime.now().isoformat(timespec='seconds')} by "
        "tools/roundtrip_bucket_report.py",
        f"from {total} individual reports in `{indiv_dir}`,",
        "classified with the v0.4.4 fac/mer-aware `oin.compare` key.",
        "",
        "## Buckets",
        "",
        "| bucket | count | % |",
        "|---|---:|---:|",
    ]
    for b in order:
        n = counts.get(b, 0)
        lines.append(f"| {b} | {n} | {100 * n / max(total, 1):.2f}% |")
    lines.append(f"| **total** | **{total}** | **100.00%** |")

    if subcounts:
        lines += ["", "### key_equal sub-split (benign canonicalization reclaimed)", ""]
        lines += ["| subclass | count |", "|---|---:|"]
        for sub, n in subcounts.most_common():
            lines.append(f"| {sub} | {n} |")

    lines += [
        "",
        "## elapsed_s percentiles",
        "",
        "| subset | n | p50 | p90 | p95 | p99 | max |",
        "|---|---:|---:|---:|---:|---:|---:|",
        _pct_row("overall", elapsed_all),
        _pct_row("eta subset", elapsed_eta),
    ]

    # Worklists every downstream swimlane consumes.
    for bucket, title in [
        ("facmer_divergent", "fac/mer-divergent (newly-caught isomer errors -- SL2/SL3 target)"),
        ("hard_fail", "hard-fail (SL4 worklist)"),
        ("encode_fail", "encode-fail (SL5 worklist)"),
        ("structural", "structural mismatch"),
    ]:
        members = [r for r in results if r["bucket"] == bucket]
        lines += ["", f"## {title} ({len(members)})", ""]
        for r in sorted(members, key=lambda x: x["molecule"] or ""):
            lines.append(f"- `{r['molecule']}`")

    return results, lines, {"total": total, "counts": dict(counts), "subcounts": dict(subcounts)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--results-dir",
        required=True,
        help="Results dir containing individual_reports/",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Where to write the report (default: the results dir).",
    )
    parser.add_argument("--only", default=None, help="Comma-separated molecule ids to include.")
    parser.add_argument("--sample", type=int, default=None, help="Evenly sample N reports.")
    parser.add_argument("--limit", type=int, default=None, help="Cap to the first N reports.")
    parser.add_argument("--shard", default=None, help="I:N -- process shard I of N.")
    args = parser.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    indiv_dir = os.path.join(results_dir, "individual_reports")
    output_dir = os.path.abspath(args.output_dir) if args.output_dir else results_dir
    only = set(args.only.split(",")) if args.only else None
    shard = tuple(int(x) for x in args.shard.split(":")) if args.shard else None

    rows = _load_reports(indiv_dir, only, args.sample, shard, args.limit)
    results, lines, summary = build_report(rows, indiv_dir)

    os.makedirs(output_dir, exist_ok=True)
    md_path = os.path.join(output_dir, "bucket_report.md")
    with open(md_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    json_path = os.path.join(output_dir, "bucket_report.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=1)

    print(f"{summary['total']} reports classified.")
    for b in [
        "byte_exact",
        "key_equal",
        "facmer_divergent",
        "structural",
        "hard_fail",
        "encode_fail",
    ]:
        n = summary["counts"].get(b, 0)
        print(f"  {n:6d}  {b:18s}  {100 * n / max(summary['total'], 1):6.2f}%")
    print(f"Wrote {md_path}")
    print(f"Wrote {json_path}")


if __name__ == "__main__":
    main()
