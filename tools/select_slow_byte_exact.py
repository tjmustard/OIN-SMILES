"""Select the slowest byte-exact round-trips from a sweep results dir.

Selection predicate over ``<results-dir>/individual_reports/*.json``:

    report["status"] == "success"
    report["tier_passed"] == "UFF_1"
    report["smiles_1"] and report["smiles_2"]
    report["smiles_1"].strip() == report["smiles_2"].strip()   # byte-exact, NOT key-equal
    isinstance(report["metrics"]["elapsed_s"], (int, float))

Sorted by ``metrics.elapsed_s`` DESC, top N (default 100).

WHY ``metrics.elapsed_s`` AND NOT A TOP-LEVEL ``elapsed_s``
============================================================
``elapsed_s`` lives **inside** each report's ``metrics`` dict. Reading
``report["elapsed_s"]`` from the top level does not raise -- it silently returns
whatever ``.get`` default you supplied (typically ``0``) for *every single row*,
which then sorts as a tie and produces a cohort that is an artifact of file-listing
order, not slowness. That exact trap has cost this project a full analysis before,
so this tool indexes ``metrics`` explicitly and treats a missing/malformed
``metrics.elapsed_s`` as EXCLUDED (not defaulted to 0).

Bands (fixed on the *selection* elapsed_s, i.e. the value this tool sorts on):
    A  > 200s
    B  100-200s
    C  < 100s

eta flag: ``re.search(r"\\{\\d+[<>]\\}", smiles_1)`` -- a haptic OIN, so the cohort
is visibly not just simple octahedral/square-planar complexes.

Emits JSON ``{names, elapsed_s, eta, band}`` (one object per selected molecule,
newline-delimited is NOT used here -- see ``--out`` for the single JSON document
written) to stdout, followed by a ``#DONE <n>`` sentinel line.

Exit codes
==========
Exits non-zero if the input directory yields 0 reports -- an empty corpus that
prints a serene "0/0" and exits 0 has already destroyed one full A/B run in this
project (see MEMORY.md). Also exits non-zero if 0 reports pass the selection
predicate.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/select_slow_byte_exact.py \\
        --results-dir /path/to/results-v0.4.5-rebaseline \\
        --n 100 --out /tmp/cohort_candidates.json
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys

HAPTIC = re.compile(r"\{\d+[<>]\}")

BAND_A_MIN = 200.0  # >200s
BAND_B_MIN = 100.0  # 100-200s


def band_for(elapsed_s: float) -> str:
    if elapsed_s > BAND_A_MIN:
        return "A"
    if elapsed_s >= BAND_B_MIN:
        return "B"
    return "C"


def load_reports(results_dir: str) -> list[dict]:
    indiv = os.path.join(results_dir, "individual_reports")
    if not os.path.isdir(indiv):
        sys.exit(f"error: {indiv} not found")
    reports = []
    for fn in sorted(os.listdir(indiv)):
        if not fn.endswith(".json"):
            continue
        path = os.path.join(indiv, fn)
        try:
            with open(path) as f:
                r = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            print(f"warning: could not parse {path}: {e}", file=sys.stderr)
            continue
        r["_basename"] = fn[: -len(".json")] + ".xyz"
        reports.append(r)
    return reports


def select(reports: list[dict]) -> list[dict]:
    """Apply the byte-exact selection predicate; return rows with derived fields."""
    selected = []
    for r in reports:
        if r.get("status") != "success":
            continue
        if r.get("tier_passed") != "UFF_1":
            continue
        s1 = r.get("smiles_1")
        s2 = r.get("smiles_2")
        if not s1 or not s2:
            continue
        if s1.strip() != s2.strip():
            continue
        metrics = r.get("metrics")
        if not isinstance(metrics, dict):
            continue
        elapsed = metrics.get("elapsed_s")
        if not isinstance(elapsed, (int, float)) or isinstance(elapsed, bool):
            continue
        selected.append(
            {
                "name": r["_basename"],
                "elapsed_s": float(elapsed),
                "eta": bool(HAPTIC.search(s1)),
                "band": band_for(float(elapsed)),
            }
        )
    selected.sort(key=lambda row: row["elapsed_s"], reverse=True)
    return selected


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--results-dir", required=True, help="Sweep results dir (has individual_reports/)"
    )
    ap.add_argument("--n", type=int, default=100, help="Top N by elapsed_s DESC (default 100)")
    ap.add_argument("--out", default=None, help="Optional path to also write the JSON payload")
    args = ap.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    reports = load_reports(results_dir)
    if len(reports) == 0:
        sys.exit(
            f"error: 0 reports found under {results_dir}/individual_reports -- "
            "refusing to proceed with an empty corpus"
        )

    selected = select(reports)
    if len(selected) == 0:
        sys.exit(f"error: 0 / {len(reports)} reports passed the byte-exact selection predicate")

    top = selected[: args.n]

    payload = {
        "results_dir": results_dir,
        "total_reports": len(reports),
        "passed_predicate": len(selected),
        "n_requested": args.n,
        "n_selected": len(top),
        "names": [row["name"] for row in top],
        "elapsed_s": {row["name"]: row["elapsed_s"] for row in top},
        "eta": {row["name"]: row["eta"] for row in top},
        "band": {row["name"]: row["band"] for row in top},
    }

    out_text = json.dumps(payload, indent=2)
    print(out_text)
    if args.out:
        with open(os.path.abspath(args.out), "w") as f:
            f.write(out_text)
        print(f"# wrote {os.path.abspath(args.out)}", file=sys.stderr)

    band_counts = {"A": 0, "B": 0, "C": 0}
    eta_count = 0
    for row in top:
        band_counts[row["band"]] += 1
        if row["eta"]:
            eta_count += 1
    print(
        f"# total_reports={len(reports)} passed_predicate={len(selected)} "
        f"selected={len(top)} bands={band_counts} eta={eta_count}",
        file=sys.stderr,
    )
    print(f"#DONE {len(top)}")


if __name__ == "__main__":
    main()
