#!/usr/bin/env python3
"""Phase 2: run the stratified sample with silent-degradation telemetry on.

Drives ``OIN3DGeneratorMetallogen`` directly rather than going through
``tools/test_dataset_roundtrip.py``. Two reasons: the harness has no telemetry
hook and editing it would change the shared measurement instrument mid-study,
and the round-trip outcome for every molecule in the sample is already recorded
in the sweep, so re-deriving it here would add cost without adding information.

Generator settings are copied verbatim from the harness's PASS-1 configuration
(``optimizer=None, ensemble_size=1, timeout=300, ff_params=None`` for non-quick)
so the firing rates correspond to the outcomes already on disk.

The point of the run is *not* to count how often a fallback fires. It is to
compare firing rates between strata: a site that fires just as often on clean
passing structures as on failures is benign, however alarming it looks in
isolation. S1 (cleanest passes) is the baseline that makes the rest readable.

Usage:
    uv run python tools/telemetry_run.py --strata ~/elimination-v043/telemetry_strata.json \
        --out ~/elimination-v043/telemetry_events.json
    uv run python tools/telemetry_run.py --limit 5          # smoke test
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

# Must be set before the generator imports, so every worker sees it.
os.environ["OIN_TELEMETRY"] = "1"

from oinsmiles.generation import _telemetry  # noqa: E402
from oinsmiles.generation.metallogen_adapter import (  # noqa: E402
    OIN3DGeneratorMetallogen as OIN3DGenerator,
)

REPO = Path(__file__).resolve().parent.parent
CAPSTONE = REPO / "tmCAT-tmPHOTO_xyz_dataset" / "results-capstone-v042"

# Verbatim from tools/test_dataset_roundtrip.py PASS 1, non-quick.
GEN_KWARGS = {"optimizer": None, "ensemble_size": 1, "timeout": 300, "ff_params": None}


def run_molecule(generator: OIN3DGenerator, oin: str) -> dict[str, Any]:
    """Generate once with telemetry collecting; report which sites fired."""
    started = time.monotonic()
    with _telemetry.collecting():
        try:
            generator.generate(oin)
            outcome, error = "generated", None
        except Exception as exc:
            outcome, error = "raised", f"{type(exc).__name__}: {exc}"
        counts = _telemetry.counts()
    return {
        "outcome": outcome,
        "error": error,
        "sites": counts,
        "elapsed_s": round(time.monotonic() - started, 2),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--strata", type=Path, default=Path.home() / "elimination-v043" / "telemetry_strata.json"
    )
    ap.add_argument(
        "--out", type=Path, default=Path.home() / "elimination-v043" / "telemetry_events.json"
    )
    ap.add_argument("--limit", type=int, help="cap molecules per stratum (smoke test)")
    ap.add_argument("--resume", action="store_true", help="skip molecules already in --out")
    args = ap.parse_args()

    if not _telemetry.enabled():
        raise SystemExit("OIN_TELEMETRY did not take effect; refusing to run a blind sweep")

    strata: dict[str, list[str]] = json.loads(args.strata.read_text())
    done: dict[str, Any] = {}
    if args.resume and args.out.exists():
        done = json.loads(args.out.read_text()).get("rows", {})
        print(f"resuming: {len(done)} molecules already recorded")

    generator = OIN3DGenerator(**GEN_KWARGS)
    rows: dict[str, Any] = dict(done)
    structures = CAPSTONE / "structures"

    for stratum, molecules in strata.items():
        todo = [m for m in (molecules[: args.limit] if args.limit else molecules) if m not in rows]
        print(f"\n=== {stratum}: {len(todo)} to run ({len(molecules)} in stratum) ===", flush=True)
        for i, mol in enumerate(todo, 1):
            oin_path = structures / f"{mol}.oin"
            if not oin_path.exists():
                rows[mol] = {"stratum": stratum, "outcome": "no_oin"}
                continue
            result = run_molecule(generator, oin_path.read_text().strip())
            result["stratum"] = stratum
            rows[mol] = result
            fired = ",".join(result["sites"]) or "-"
            print(
                f"  [{i}/{len(todo)}] {mol:<24} {result['outcome']:<10} {result['elapsed_s']:>7}s  {fired}",
                flush=True,
            )
            # Persist as we go: a long run must survive an interruption.
            args.out.write_text(json.dumps({"rows": rows, "gen_kwargs": GEN_KWARGS}, indent=1))

    print("\n" + "=" * 74)
    print("PER-SITE FIRING RATE BY STRATUM")
    print("=" * 74)
    by_stratum: dict[str, list] = {}
    for mol, r in rows.items():
        by_stratum.setdefault(r.get("stratum", "?"), []).append(r)

    all_sites = sorted({s for r in rows.values() for s in r.get("sites", {})})
    if not all_sites:
        print("  no site fired anywhere -- every probed fallback is dormant on this sample")
    for site in all_sites:
        print(f"\n  {site}")
        for stratum in sorted(by_stratum):
            rs = by_stratum[stratum]
            n = len(rs)
            hit = sum(1 for r in rs if site in r.get("sites", {}))
            print(f"    {stratum:<20} {hit:4d}/{n:<4d}  {100.0 * hit / n if n else 0:5.1f}%")

    print("\n  outcomes by stratum:")
    for stratum in sorted(by_stratum):
        print(f"    {stratum:<20} {dict(Counter(r.get('outcome') for r in by_stratum[stratum]))}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
