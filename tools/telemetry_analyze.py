#!/usr/bin/env python3
"""Reproduce the Phase-2 firing-rate tables in docs/FALSIFICATION_v0.4.3_ELIMINATION.md §8b.

Reads the committed telemetry run (``docs/data/v0.4.3/telemetry_events.json``, produced by
``tools/telemetry_run.py``) and its stratified sample definition, and recomputes:

  * per-site firing rate by stratum (fraction of molecules where the site fired >=1x),
  * odds ratio vs the S1 clean-pass baseline with a Haldane-Anscombe (+0.5) correction and
    a Wald 95% CI, decided as ENRICHED / depleted / no difference,
  * mean firings per molecule (intensity), and
  * generation outcomes by stratum.

The Haldane-Anscombe correction is what lets an odds ratio be defined even when a cell is
zero (several sites never fire in a given stratum), at the cost of shrinking the estimate
toward 1 -- which is the conservative direction for an elimination study.

Usage:
    uv run python tools/telemetry_analyze.py
    uv run python tools/telemetry_analyze.py --events <path> --strata <path>
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "docs" / "data" / "v0.4.3"

# Sample strata in reporting order; the first is the baseline every OR is taken against.
STRATA_ORDER = ["S1_pass_clean", "S2_pass_distorted", "S3_bucketE", "S4_bucketBD"]
BASELINE = "S1_pass_clean"
ENRICH_LO = 1.5  # OR CI lower bound above this => ENRICHED (pre-registered threshold)


def odds_ratio(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    """OR of exposure (c/d vs a/b) with Haldane-Anscombe correction and Wald 95% CI."""
    A, B, C, D = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_ = (C / D) / (A / B)
    se = math.sqrt(1 / A + 1 / B + 1 / C + 1 / D)
    return or_, math.exp(math.log(or_) - 1.96 * se), math.exp(math.log(or_) + 1.96 * se)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--events", type=Path, default=DATA / "telemetry_events.json")
    ap.add_argument("--strata", type=Path, default=DATA / "telemetry_strata.json")
    args = ap.parse_args()

    rows = json.loads(args.events.read_text())["rows"]
    strata = json.loads(args.strata.read_text())
    mols = {k: [m for m in v if m in rows] for k, v in strata.items()}
    order = [k for k in STRATA_ORDER if k in mols] + [k for k in mols if k not in STRATA_ORDER]
    sites = sorted({s for r in rows.values() for s in r.get("sites", {})})

    def fired(stratum: str, site: str) -> int:
        return sum(1 for m in mols[stratum] if site in rows[m].get("sites", {}))

    print(f"n = {len(rows)} molecules over {len(order)} strata\n")
    print(f"{'site':<40}" + "".join(f"{k.split('_')[0]:>9}" for k in order))
    print("-" * (40 + 9 * len(order)))
    for site in sites:
        line = f"{site:<40}"
        for k in order:
            n = len(mols[k])
            line += f"{100 * fired(k, site) / n if n else 0:8.1f}%"
        print(line)

    print(f"\nODDS RATIO vs {BASELINE} (Haldane-Anscombe, Wald 95% CI)")
    a = None
    for site in sites:
        hits = fired(BASELINE, site)
        a, b = hits, len(mols[BASELINE]) - hits
        print(f"\n  {site}")
        for k in order:
            if k == BASELINE:
                continue
            c, d = fired(k, site), len(mols[k]) - fired(k, site)
            or_, lo, hi = odds_ratio(a, b, c, d)
            verdict = "ENRICHED" if lo > ENRICH_LO else ("depleted" if hi < 1 else "no difference")
            print(f"    vs {k:<20} OR={or_:6.2f}  CI[{lo:5.2f},{hi:7.2f}]  {verdict}")

    print("\nMEAN FIRINGS PER MOLECULE (intensity)")
    for site in sites:
        vals = "".join(
            f"{sum(rows[m].get('sites', {}).get(site, 0) for m in mols[k]) / len(mols[k]):8.1f}"
            for k in order
        )
        print(f"  {site:<40}{vals}")

    print("\nOUTCOMES BY STRATUM")
    for k in order:
        print(f"  {k:<20} {dict(Counter(rows[m].get('outcome') for m in mols[k]))}")


if __name__ == "__main__":
    main()
