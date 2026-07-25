"""Aggregate injectivity probes + a dataset population-at-risk scan into a report.

Two jobs:

1. **Curated probes** -- run the named Y1 twin probes (the confirmed blind-spot fixtures)
   and roll their verdicts into a confusion-matrix summary + a falsification verdict.
2. **Population-at-risk** -- sample the dataset, run the full twin probe on each, and
   estimate what fraction of *chiral* structures the encoder collides with its mirror.
   This is the blast radius: how much of the corpus a chirality blind spot can silently
   mis-pass through the round trip.

Deterministic: no timestamps, seeded sampling -- re-running regenerates byte-identical
``injectivity_metrics.json`` and ``report.md``.

Run:
  PYTHONPATH=$PWD/src python -m tools.injectivity.report --probes
  PYTHONPATH=$PWD/src python -m tools.injectivity.report --population 300 --dataset <dir>
"""

from __future__ import annotations

import argparse
import json
import random
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from .twin_collision import (
    VERDICT_DISTINGUISHED,
    VERDICT_ENCODER_BLIND,
    VERDICT_INVARIANT_OK,
    VERDICT_KEY_BLIND,
    VERDICT_OVER_SENSITIVE,
    ProbeOutcome,
    TwinProbe,
    probe_mirror,
    run_probes,
)

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "results-injectivity-y1"
DEFAULT_DATASET = Path("/home/tjmustard/Documents/GitHub/tmCat-tmPhoto/tmCAT-tmPHOTO_xyz_dataset")
SAMPLE_SEED = 42


def default_probe_set() -> list[TwinProbe]:
    """The curated Y1 mirror-twin probes -- fixtures whose chirality is known by construction.

    Extended by each probe swimlane as it lands its fixtures. Only fixtures that exist in
    this checkout are returned, so a partial worktree still runs.
    """
    fx = REPO / "tests" / "fixtures"
    candidates = [
        TwinProbe("CisPlatin (achiral control)", str(fx / "CisPlatin.xyz"), "none (achiral)"),
        TwinProbe("fac-Ir(ppy)3", str(fx / "fac-Ir(ppy)3.xyz"), "metal Δ/Λ"),
        TwinProbe("PdCl2-R-BINAP", str(fx / "PdCl2-R-BINAP.xyz"), "axial (atropisomer)"),
        TwinProbe("POJJOP (metal-bound 2° amine)", str(fx / "POJJOP.xyz"), "amine flip R→S (P3)"),
        # probe swimlanes drop additional fixtures here as they land.
    ]
    return [p for p in candidates if Path(p.xyz_path).exists()]


@dataclass
class Summary:
    outcomes: list[ProbeOutcome]

    def counts(self) -> dict:
        return dict(Counter(o.verdict for o in self.outcomes))

    def collisions(self) -> list[ProbeOutcome]:
        return [o for o in self.outcomes if o.is_collision]

    def verdict_h0(self) -> str:
        """H0 = 'a passing round-trip proves the OIN is lossless'. Any collision refutes it."""
        return "REFUTED" if self.collisions() else "SURVIVES"


def _iter_dataset(dataset: Path) -> list[Path]:
    return sorted(p for sub in ("cat", "photo") for p in (dataset / sub).glob("*.xyz"))


#: Loud caveat carried in the metrics + report so the number is never cited as a rate.
_POP_CAVEAT = (
    "UPPER BOUND ONLY. The oracle here is a RIGID mirror-superposition test, which is valid "
    "only for rigid molecules. On flexible crystal-conformer structures the mirror is a "
    "non-superimposable CONFORMER (not a configurational enantiomer), so 'distinct' is "
    "massively over-reported and 'collision' with it -- because a correctly "
    "conformer-invariant encoder gives the same string -- is inflated accordingly. This is "
    "NOT a losslessness-failure rate. A configurational oracle (Wave 3) is required for a "
    "trustworthy per-axis population. The curated rigid fixtures in the probe table above are "
    "the sound Y1 result; this scan only bounds the blast radius from above."
)


def population_scan(dataset: Path, n: int) -> dict:
    """Sample n structures deterministically; run the full mirror probe on each.

    See ``_POP_CAVEAT``: the rigid oracle conflates conformation with configuration at
    dataset scale, so every fraction below is an inflated UPPER BOUND, not a rate.
    """
    files = _iter_dataset(dataset)
    if not files:
        return {"error": f"no xyz files under {dataset}", "sampled": 0}
    n = min(n, len(files))
    sample = random.Random(SAMPLE_SEED).sample(files, n)
    verdicts: Counter = Counter()
    errors = 0
    rigid_distinct = 0
    for path in sample:
        try:
            o = probe_mirror(path)
        except Exception:
            errors += 1
            continue
        verdicts[o.verdict] += 1
        if o.oracle_distinct:
            rigid_distinct += 1
    collisions = verdicts[VERDICT_ENCODER_BLIND] + verdicts[VERDICT_KEY_BLIND]
    return {
        "caveat": _POP_CAVEAT,
        "dataset": str(dataset),
        "corpus_size": len(files),
        "sampled": n,
        "errors": errors,
        "rigid_mirror_distinct": rigid_distinct,
        "rigid_mirror_distinct_fraction_UPPER_BOUND": round(rigid_distinct / n, 4) if n else 0.0,
        "collisions_UPPER_BOUND": collisions,
        "collision_fraction_of_sample_UPPER_BOUND": round(collisions / n, 4) if n else 0.0,
        "verdict_breakdown": dict(verdicts),
    }


def _matrix_line(o: ProbeOutcome) -> str:
    flag = {
        VERDICT_ENCODER_BLIND: "🔴 encoder-blind (total)",
        VERDICT_KEY_BLIND: "🟠 key-blind (batch FP)",
        VERDICT_DISTINGUISHED: "🟢 distinguished",
        VERDICT_INVARIANT_OK: "⚪ invariant ok",
        VERDICT_OVER_SENSITIVE: "🟡 over-sensitive (FN)",
    }.get(o.verdict, o.verdict)
    return (
        f"| {o.name} | {o.operator} | {o.oracle_distinct} ({o.oracle_rmsd} Å) "
        f"| {o.raw_equal} | {o.key_equal} | {flag} |"
    )


def render_markdown(summary: Summary, population: dict | None) -> str:
    lines = [
        "# Injectivity Probe Report (Y1)",
        "",
        "**Status:** measurement only -- no encoder code changed.",
        "",
        "Auto-generated by `tools/injectivity/report.py`. Every probe encodes a base",
        "structure and its z-mirror twin and asks whether the encoder / round-trip key can",
        "still tell two genuinely different isomers apart. No 3D generator is invoked.",
        "",
        f'## H0 ("a passing round-trip proves the OIN is lossless"): **{summary.verdict_h0()}**',
        "",
        f"Collisions found: **{len(summary.collisions())}** "
        f"(a collision = oracle-distinct isomers the round-trip key calls equal).",
        "",
        "| fixture | operator | oracle-distinct | raw_equal | key_equal | verdict |",
        "|---|---|---|---|---|---|",
    ]
    lines += [_matrix_line(o) for o in summary.outcomes]
    lines += [
        "",
        "Legend: 🔴 total encoder blindness (byte-identical strings for distinct isomers) ·",
        "🟠 key-blind (raw differs only by non-reproducible slot artifact; the batch harness",
        "gates on the key, so it passes) · 🟢 encoder injective on this axis · ⚪ mirror is the",
        "same isomer (correct reflection-invariance) · 🟡 over-sensitive (a false negative).",
    ]
    if population:
        lines += [
            "",
            "## Population-at-risk (dataset sample)",
            "",
            f"> ⚠️ **{population.get('caveat', '')}**",
            "",
            "```json",
            json.dumps(population, indent=2),
            "```",
            "",
            "The fractions above are conformation-inflated UPPER BOUNDS, not losslessness-failure",
            "rates -- see the caveat. A trustworthy per-axis population needs the configurational",
            "oracle deferred to Wave 3 (the UU hunt).",
        ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--probes", action="store_true", help="run the curated Y1 twin probes")
    ap.add_argument("--population", type=int, default=0, help="sample N dataset structures")
    ap.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)

    if not args.probes and not args.population:
        args.probes = True  # default action

    outcomes = run_probes(default_probe_set()) if args.probes else []
    summary = Summary(outcomes)
    population = population_scan(args.dataset, args.population) if args.population else None

    args.out.mkdir(parents=True, exist_ok=True)
    metrics = {
        "h0_verdict": summary.verdict_h0(),
        "verdict_counts": summary.counts(),
        "n_collisions": len(summary.collisions()),
        "probes": [o.to_dict() for o in outcomes],
        "population_at_risk": population,
    }
    (args.out / "injectivity_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    (args.out / "report.md").write_text(render_markdown(summary, population))

    print(f"H0: {summary.verdict_h0()}  |  collisions: {len(summary.collisions())}")
    for o in outcomes:
        print(f"  [{o.verdict:14}] {o.name}")
    if population:
        print(f"  population: {json.dumps(population)}")
    print(f"wrote {args.out}/injectivity_metrics.json and report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
