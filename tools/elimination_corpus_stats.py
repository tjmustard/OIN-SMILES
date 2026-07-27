#!/usr/bin/env python3
"""Corpus-level statistics for the v0.4.3 round-trip failure elimination study.

Read-only. Mines existing sweep output; never runs the generator.

Produces the three framing results the elimination study depends on:

  E0  Failure partition by pipeline stage reached (buckets A-E), plus the
      quick vs non-quick adjudication that sizes the ``--quick`` budget artifact.
  E2b Flake floor: how often a repeated run of the same molecule flips its gate
      outcome. Every effect size in the study is compared against this.
  P0  Provenance: schema drift between result sets and the disposition of rows
      that are neither pass nor fail (``pending_g-xtb``).

Buckets are assigned by how far down the pipeline a failure got, which is
recoverable from the stored strings alone:

  A_encoder        no ``smiles_1``      -> XYZ->OIN never emitted a string, so
                                          neither MetalloGen nor the OIN->m-SMILES
                                          handoff was ever reached.
  B_no_conformers  ``smiles_1`` only, "any conformers" in the error
  C_timeout        ``smiles_1`` only, timeout in the error
  D_gen_other      ``smiles_1`` only, anything else
  E_gate           both strings present -> a real string-identity or RMSD
                                          disagreement to adjudicate.

Usage:
    uv run python tools/elimination_corpus_stats.py
    uv run python tools/elimination_corpus_stats.py --out results/e0_stats.json
"""

from __future__ import annotations

import argparse
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

REPO = Path(__file__).resolve().parent.parent
DATASET = REPO / "tmCAT-tmPHOTO_xyz_dataset"

QUICK_SET = DATASET / "results-v0.4.0"
NONQUICK_SET = DATASET / "results-capstone-v042"
DETERMINISM_DIR = Path.home() / "capstone-v042" / "determinism"

BUCKETS = ("A_encoder", "B_no_conformers", "C_timeout", "D_gen_other", "E_gate")


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def load_summary(results_dir: Path) -> list[dict[str, Any]]:
    """Load a sweep's summary rows, tolerating a partially-written JSON tail."""
    path = results_dir / "summary_roundtrip.json"
    if not path.exists():
        raise SystemExit(f"no summary at {path}")
    text = path.read_text()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        # A sweep killed mid-write leaves a truncated array; salvage what parsed.
        cut = text.rfind("}")
        return json.loads(text[: cut + 1] + "]")


def classify_bucket(row: dict[str, Any]) -> str:
    """Assign a failed row to the pipeline stage it reached."""
    if not row.get("smiles_1"):
        return "A_encoder"
    if row.get("smiles_2"):
        return "E_gate"
    err = row.get("error") or ""
    if "any conformers" in err:
        return "B_no_conformers"
    if "imeout" in err or "exceeded" in err:
        return "C_timeout"
    return "D_gen_other"


def is_pass(row: dict[str, Any]) -> bool:
    return row.get("status") == "success"


def is_pending(row: dict[str, Any]) -> bool:
    return str(row.get("status", "")).startswith("pending")


# --------------------------------------------------------------------------
# E0 - partition + adjudication
# --------------------------------------------------------------------------
def partition(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    rows = list(rows)
    failed = [r for r in rows if not is_pass(r) and not is_pending(r)]
    buckets = Counter(classify_bucket(r) for r in failed)

    # Within bucket A, separate a timeout from a genuine perception failure:
    # they have different owners (budget vs perception) and must not be merged.
    a_rows = [r for r in failed if classify_bucket(r) == "A_encoder"]
    a_split = Counter(
        "encode_timeout" if "imeout" in (r.get("error") or "") else "encode_perception"
        for r in a_rows
    )

    return {
        "n_rows": len(rows),
        "n_pass": sum(1 for r in rows if is_pass(r)),
        "n_pending": sum(1 for r in rows if is_pending(r)),
        "n_fail": len(failed),
        "buckets": {b: buckets.get(b, 0) for b in BUCKETS},
        "bucket_pct_of_failures": {
            b: round(100.0 * buckets.get(b, 0) / len(failed), 2) for b in BUCKETS
        }
        if failed
        else {},
        "A_encoder_split": dict(a_split),
    }


def adjudicate(quick: list[dict], nonquick: list[dict]) -> dict[str, Any]:
    """Re-judge quick-mode outcomes against the non-quick sweep on the overlap.

    The two sweeps differ in generator budget, so a quick failure that passes
    non-quick was never a correctness defect -- it was the 30s wall.
    """
    nq = {r["molecule"]: r for r in nonquick}
    table: dict[str, Counter] = defaultdict(Counter)
    flips_pass_to_fail: list[str] = []

    for r in quick:
        other = nq.get(r["molecule"])
        if other is None:
            continue
        if is_pending(other):
            outcome = "pending"
        else:
            outcome = "pass" if is_pass(other) else "fail"
        key = "PASS_quick" if is_pass(r) else classify_bucket(r)
        table[key][outcome] += 1
        if is_pass(r) and outcome == "fail":
            flips_pass_to_fail.append(r["molecule"])

    fail_keys = [k for k in table if k != "PASS_quick"]
    adjudicated = sum(table[k]["pass"] + table[k]["fail"] for k in fail_keys)
    rescued = sum(table[k]["pass"] for k in fail_keys)

    return {
        "overlap": sum(sum(c.values()) for c in table.values()),
        "by_quick_outcome": {k: dict(v) for k, v in sorted(table.items())},
        "quick_failures_adjudicated": adjudicated,
        "quick_failures_rescued_by_budget": rescued,
        "rescue_rate_pct": round(100.0 * rescued / adjudicated, 2) if adjudicated else None,
        "quick_pass_to_nonquick_fail": len(flips_pass_to_fail),
        "quick_pass_to_nonquick_fail_examples": sorted(flips_pass_to_fail)[:20],
    }


# --------------------------------------------------------------------------
# E2b - flake floor
# --------------------------------------------------------------------------
def flake_floor(det_dir: Path) -> dict[str, Any]:
    """Measure run-to-run outcome instability from repeated single-molecule runs.

    Directories are named ``<molecule>_r<N>``; replicates of one molecule differ
    only by run index, so any disagreement in outcome is generator stochasticity
    rather than a code or input difference.
    """
    if not det_dir.exists():
        return {"available": False, "reason": f"{det_dir} not found"}

    reps: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for sub in sorted(det_dir.iterdir()):
        if not sub.is_dir() or "_r" not in sub.name:
            continue
        mol = sub.name.rsplit("_r", 1)[0]
        for report in (sub / "individual_reports").glob("*.json"):
            try:
                reps[mol].append(json.loads(report.read_text()))
            except (json.JSONDecodeError, OSError):
                continue

    multi = {m: rs for m, rs in reps.items() if len(rs) > 1}
    disagreeing, rmsd_spreads = [], []
    for mol, rs in multi.items():
        outcomes = {is_pass(r) for r in rs}
        if len(outcomes) > 1:
            disagreeing.append(mol)
        vals = [
            r["metrics"]["rmsd"]
            for r in rs
            if isinstance(r.get("metrics"), dict)
            and isinstance(r["metrics"].get("rmsd"), (int, float))
            and r["metrics"]["rmsd"] < 900
        ]
        if len(vals) > 1:
            rmsd_spreads.append(max(vals) - min(vals))

    return {
        "available": True,
        "molecules_with_replicates": len(multi),
        "replicate_counts": dict(Counter(len(rs) for rs in multi.values())),
        "molecules_flipping_outcome": len(disagreeing),
        "flake_floor_pct": round(100.0 * len(disagreeing) / len(multi), 2) if multi else None,
        "flipping_examples": sorted(disagreeing)[:20],
        "rmsd_spread_median": round(statistics.median(rmsd_spreads), 4) if rmsd_spreads else None,
        "rmsd_spread_max": round(max(rmsd_spreads), 4) if rmsd_spreads else None,
        "n_rmsd_spreads": len(rmsd_spreads),
    }


# --------------------------------------------------------------------------
# P0 - provenance
# --------------------------------------------------------------------------
def provenance(name: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    fields: Counter = Counter()
    for r in rows:
        fields.update(r.keys())
    n = len(rows)
    return {
        "set": name,
        "n_rows": n,
        "commit_ids": dict(Counter(r.get("commit_id") for r in rows).most_common(5)),
        "rdkit_versions": dict(Counter(r.get("rdkit_version") for r in rows).most_common(5)),
        "quick": dict(Counter(r.get("quick") for r in rows)),
        "xtb_available": dict(Counter(r.get("xtb_available") for r in rows)),
        "optimizer_effective": dict(Counter(r.get("optimizer_effective") for r in rows)),
        "statuses": dict(Counter(r.get("status") for r in rows)),
        "fields_present_on_all_rows": sorted(k for k, c in fields.items() if c == n),
        "fields_missing_on_some_rows": {k: n - c for k, c in sorted(fields.items()) if c < n},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick-dir", type=Path, default=QUICK_SET)
    ap.add_argument("--nonquick-dir", type=Path, default=NONQUICK_SET)
    ap.add_argument("--determinism-dir", type=Path, default=DETERMINISM_DIR)
    ap.add_argument("--out", type=Path, help="write the full result as JSON")
    args = ap.parse_args()

    quick = load_summary(args.quick_dir)
    nonquick = load_summary(args.nonquick_dir)

    result = {
        "E0_partition_quick": partition(quick),
        "E0_partition_nonquick": partition(nonquick),
        "E0_adjudication": adjudicate(quick, nonquick),
        "E2b_flake_floor": flake_floor(args.determinism_dir),
        "P0_provenance": [
            provenance(args.quick_dir.name, quick),
            provenance(args.nonquick_dir.name, nonquick),
        ],
    }

    print(json.dumps(result, indent=2, default=str))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2, default=str))
        print(f"\nwrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
