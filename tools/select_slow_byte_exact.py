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

SELECTION IS FROM **ONE** DIR. TWO DIRS ARE NOT A ROBUSTNESS FILTER.
=====================================================================
An earlier version of this cohort intersected two "independent" results dirs
under a "slow in both runs screens out contention flukes" rule. That rule
presumes both dirs are repeat measurements of the *same population* -- they were
not (one was a partial of the frozen seed-42 5k cohort, the other a differently
-drawn "gap ∪ guard" re-baseline cohort), so intersecting them was a near-total
sample destruction (~25,000-molecule universe -> ~100 shared names -> 62 passing
both), not a filter. **Select top-N from a single results dir.** A second dir may
be used purely for *corroboration* (see ``--corroborate-with``) -- reporting a
second elapsed_s observation where one happens to exist -- never as a filter that
molecules must additionally survive.

``--corroborate-with`` (OPTIONAL, INFORMATIONAL ONLY)
=======================================================
When given a second results dir, each selected molecule's ``metrics.elapsed_s``
is looked up there too *if it passes the same full predicate* in that dir, and
recorded under ``corroboration_elapsed_s`` (``null`` if absent or it didn't pass
there -- absence is fine and expected, not an error). Independently, if the
corroboration dir has ANY report for that name with a non-null ``smiles_1``
(regardless of tier/status), it is compared byte-for-byte against the primary
source's ``smiles_1`` as an encoder-determinism sanity check and recorded under
``corroboration_smiles1_match`` (``true``/``false``/``null``) -- a ``false`` here
would be a genuine anomaly worth investigating (same input molecule, different
encoded string) and is warned about on stderr, but does NOT remove the molecule
from the selection: corroboration reports, it does not gate.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/select_slow_byte_exact.py \\
        --results-dir /path/to/results-v0.4.5-sweep-partial-2697mols \\
        --n 100 --corroborate-with /path/to/results-v0.4.5-rebaseline \\
        --out /tmp/cohort_candidates.json
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


def _passes_full_predicate(r: dict) -> float | None:
    """Return elapsed_s if ``r`` passes the full selection predicate, else None."""
    if r.get("status") != "success":
        return None
    if r.get("tier_passed") != "UFF_1":
        return None
    s1 = r.get("smiles_1")
    s2 = r.get("smiles_2")
    if not s1 or not s2:
        return None
    if s1.strip() != s2.strip():
        return None
    metrics = r.get("metrics")
    if not isinstance(metrics, dict):
        return None
    elapsed = metrics.get("elapsed_s")
    if not isinstance(elapsed, (int, float)) or isinstance(elapsed, bool):
        return None
    return float(elapsed)


def corroborate(
    names: list[str], primary_smiles1: dict[str, str], corroboration_dir: str
) -> tuple[dict[str, float | None], dict[str, bool | None]]:
    """Look up ``names`` in a second results dir. Informational only -- never a filter.

    Returns ``(elapsed_by_name, smiles1_match_by_name)``:
      * ``elapsed_by_name[name]`` is the corroboration dir's ``metrics.elapsed_s`` if that
        report ALSO passes the full predicate there, else ``None`` (absence/non-pass is
        fine and expected, not an error).
      * ``smiles1_match_by_name[name]`` is a byte-for-byte comparison of ``smiles_1``
        against the primary source, IF the corroboration dir has any report for that name
        with a non-null ``smiles_1`` (regardless of tier/status) -- an encoder-determinism
        sanity check, independent of whether the round trip passed there. ``None`` if the
        corroboration dir has no report, or that report has no ``smiles_1``, for this name.
    """
    reports_by_name = {r["_basename"]: r for r in load_reports(corroboration_dir)}
    elapsed_by_name: dict[str, float | None] = {}
    match_by_name: dict[str, bool | None] = {}
    for name in names:
        r = reports_by_name.get(name)
        if r is None:
            elapsed_by_name[name] = None
            match_by_name[name] = None
            continue
        elapsed_by_name[name] = _passes_full_predicate(r)
        s1 = r.get("smiles_1")
        if not s1:
            match_by_name[name] = None
            continue
        match = s1.strip() == primary_smiles1[name].strip()
        match_by_name[name] = match
        if not match:
            print(
                f"warning: {name} smiles_1 DIFFERS between primary and corroboration "
                f"source -- possible encoder non-determinism, investigate",
                file=sys.stderr,
            )
    return elapsed_by_name, match_by_name


def select(reports: list[dict]) -> list[dict]:
    """Apply the byte-exact selection predicate; return rows with derived fields."""
    selected = []
    for r in reports:
        elapsed = _passes_full_predicate(r)
        if elapsed is None:
            continue
        s1 = r["smiles_1"]
        selected.append(
            {
                "name": r["_basename"],
                "elapsed_s": elapsed,
                "eta": bool(HAPTIC.search(s1)),
                "band": band_for(elapsed),
                "smiles_1": s1.strip(),
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
    ap.add_argument(
        "--corroborate-with",
        default=None,
        help=(
            "Optional second results dir. INFORMATIONAL ONLY -- records a second "
            "elapsed_s observation where one happens to exist, and cross-checks "
            "smiles_1 for encoder determinism. Never filters the selection."
        ),
    )
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

    corroboration_count = None
    if args.corroborate_with:
        corroboration_dir = os.path.abspath(args.corroborate_with)
        primary_smiles1 = {row["name"]: row["smiles_1"] for row in top}
        names = [row["name"] for row in top]
        corr_elapsed, corr_match = corroborate(names, primary_smiles1, corroboration_dir)
        payload["corroboration_source"] = corroboration_dir
        payload["corroboration_elapsed_s"] = corr_elapsed
        payload["corroboration_smiles1_match"] = corr_match
        corroboration_count = sum(1 for v in corr_elapsed.values() if v is not None)

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
    corr_msg = (
        f" corroborated={corroboration_count}/{len(top)}" if corroboration_count is not None else ""
    )
    print(
        f"# total_reports={len(reports)} passed_predicate={len(selected)} "
        f"selected={len(top)} bands={band_counts} eta={eta_count}{corr_msg}",
        file=sys.stderr,
    )
    print(f"#DONE {len(top)}")


if __name__ == "__main__":
    main()
