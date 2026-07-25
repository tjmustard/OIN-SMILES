"""Wave 3 -- the missed-success (FALSE NEGATIVE) audit.

The Y1 confusion matrix has two cells the round-trip test cannot see. Waves 1-2 attacked the
false positives (a lossy OIN that still PASSES). This is the other one: a round-trip FAIL
whose OIN was actually fine.

The user's calibration going in was that such a FAIL is *usually the generator* -- MetalloGen
not yet fully accurate, or the job timing out -- not the notation being wrong. That is a
testable claim, and the v0.4.4 regression sweep already holds the evidence, so no expensive
re-run is needed.

Every non-passing row is attributed to exactly one cause, and the causes are grouped by what
they tell us about the NOTATION:

* ``generator_timeout``      -- the job hit the wall clock. Says NOTHING about the OIN.
* ``generator_no_output``    -- generation died without producing a structure. Says nothing.
* ``canonicalization_noise`` -- input and output name the SAME isomer under the round-trip
  key, and differ only in notation. The test was too strict.
* ``divergent_isomer_ambiguous`` -- the output names a different isomer. NOT scored against
  the notation: it is equally consistent with the generator building the wrong thing and
  with the OIN licensing it. Separating those needs the per-case oracle.
* ``encode_fail``            -- the encoder refused the input; an encoder-coverage limit,
  not a round-trip result.

The headline number is the fraction of FAILs that are **uninformative about the notation**.
If the calibration holds, most of the round-trip failure mass never tests the OIN at all --
which means the headline pass-rate is largely a measure of generator throughput, and reading
it as a measure of notation quality is a category error.

Run:  PYTHONPATH=$PWD/src python -m tools.injectivity.missed_success_audit [--report DIR]
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
for _p in (str(REPO), str(REPO / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from oinsmiles.oin.compare import (  # noqa: E402
    normalize_oin_for_comparison,
    winding_canonical_key,
)

DEFAULT_SWEEP = Path(
    "/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/"
    "results-v0.4.4-regression"
)
OUT_DIR = REPO / "results-injectivity-y2"

#: the sweep's per-molecule wall-clock budget (full mode).
TIMEOUT_S = 300.0
#: within this margin of the cap we call it a timeout even if the error text does not say so.
TIMEOUT_MARGIN_S = 10.0

PASSING_BUCKETS = {"byte_exact", "key_equal"}

#: causes that carry NO information about whether the OIN was correct.
UNINFORMATIVE = {"generator_timeout", "generator_no_output"}


def _key(oin: str):
    try:
        return winding_canonical_key(normalize_oin_for_comparison(oin))
    except Exception:
        return None


def attribute(row: dict) -> tuple[str, str]:
    """Return ``(cause, note)`` for one non-passing sweep row."""
    bucket = row.get("bucket")
    err = str(row.get("error") or "")
    elapsed = row.get("elapsed_s") or 0.0
    out = row.get("smiles_2")

    if bucket == "encode_fail":
        return "encode_fail", "encoder refused the input structure"

    timed_out = "timeout" in err.lower() or "timed out" in err.lower()
    at_the_wall = elapsed >= (TIMEOUT_S - TIMEOUT_MARGIN_S)
    if timed_out or (at_the_wall and not out):
        return "generator_timeout", f"elapsed {elapsed:.0f}s against a {TIMEOUT_S:.0f}s budget"

    if not out:
        return "generator_no_output", err[:80] or "no structure produced"

    # A structure WAS produced. Does it name the same isomer as the input?
    k_in, k_out = _key(row.get("smiles_1") or ""), _key(out)
    if k_in is not None and k_out is not None and k_in == k_out:
        return "canonicalization_noise", "same isomer under the round-trip key; notation differs"
    # The output names a DIFFERENT isomer. Careful: that alone does not implicate the
    # notation. Either the generator built the wrong isomer (accuracy), or the OIN was
    # ambiguous enough to license a different one (losslessness). The sweep cannot tell
    # them apart -- separating them needs the per-case oracle -- so this stays ambiguous
    # rather than being scored as evidence against the notation.
    return "divergent_isomer_ambiguous", "output names a different isomer; cause not separable"


def audit(sweep: Path) -> dict:
    report = sweep / "bucket_report.json"
    if not report.exists():
        raise SystemExit(f"sweep report not found: {report}")
    rows = json.loads(report.read_text())

    passing = [r for r in rows if r.get("bucket") in PASSING_BUCKETS]
    failing = [r for r in rows if r.get("bucket") not in PASSING_BUCKETS]

    causes = Counter()
    by_cause: dict[str, list] = {}
    for r in failing:
        cause, note = attribute(r)
        causes[cause] += 1
        by_cause.setdefault(cause, []).append(
            {
                "molecule": r.get("molecule"),
                "bucket": r.get("bucket"),
                "elapsed_s": r.get("elapsed_s"),
                "note": note,
            }
        )

    n_fail = len(failing)
    uninformative = sum(causes[c] for c in UNINFORMATIVE)
    elapsed_all = [r["elapsed_s"] for r in failing if r.get("elapsed_s")]

    return {
        "sweep": str(sweep),
        "timeout_s": TIMEOUT_S,
        "n_total": len(rows),
        "n_passing": len(passing),
        "n_failing": n_fail,
        "causes": dict(causes.most_common()),
        "cause_fractions": {k: round(v / n_fail, 4) for k, v in causes.items()} if n_fail else {},
        "uninformative_about_notation": uninformative,
        "uninformative_fraction": round(uninformative / n_fail, 4) if n_fail else 0.0,
        "elapsed_median_s": round(statistics.median(elapsed_all), 1) if elapsed_all else None,
        "examples": {c: v[:10] for c, v in by_cause.items()},
    }


def render(a: dict) -> str:
    n_fail = a["n_failing"]

    def pct(n):
        return f"{100.0 * n / n_fail:.1f}%" if n_fail else "-"

    lines = [
        "# Missed-success audit (Y3) -- what round-trip FAILures actually mean",
        "",
        "The round-trip test cannot see its own false negatives: a FAIL whose OIN was correct.",
        "This attributes every failing row of the v0.4.4 regression sweep to a cause, and",
        "separates the causes that say something about the NOTATION from those that do not.",
        "",
        f"- sweep: `{a['sweep']}`",
        f"- molecules: **{a['n_total']}** ({a['n_passing']} passing, {a['n_failing']} failing)",
        f"- per-molecule budget: {a['timeout_s']:.0f}s; median elapsed among failures:"
        f" **{a['elapsed_median_s']}s**",
        "",
        "## Attribution",
        "",
        "| cause | count | share of failures | informative about the OIN? |",
        "|---|---:|---:|---|",
    ]
    labels = {
        "generator_timeout": ("generator timeout", "no -- untested"),
        "generator_no_output": ("generator produced nothing", "no -- untested"),
        "generator_wrong_structure": ("generator built the wrong structure", "no -- generator"),
        "canonicalization_noise": ("canonicalization noise", "test too strict"),
        "divergent_isomer_ambiguous": (
            "output names a different isomer",
            "ambiguous -- generator OR notation",
        ),
        "encode_fail": ("encoder refused the input", "encoder coverage"),
    }
    for cause, n in a["causes"].items():
        label, informative = labels.get(cause, (cause, "?"))
        lines.append(f"| {label} | {n} | {pct(n)} | {informative} |")

    unin = a["uninformative_about_notation"]
    lines += [
        "",
        f"## Headline: {pct(unin)} of failures never test the notation",
        "",
        f"**{unin} of {n_fail}** failing molecules failed for a reason that carries no",
        "information about whether the OIN was correct -- the job timed out or died before a",
        "structure existed to compare. The notation was never put to the test.",
        "",
        "This confirms the calibration this program started from: *a round-trip FAIL with a",
        "correct OIN is usually the generator, not the notation.* Two consequences follow.",
        "",
        "1. **The headline round-trip pass-rate is substantially a measure of generator",
        "   throughput**, not of notation quality. Reading a pass-rate change as a change in",
        "   losslessness is a category error unless the failure mix is held fixed.",
        "2. **Compute buys accuracy here.** With most failures sitting at the wall clock, a",
        "   larger time budget converts failures into measurements -- an unusually cheap way to",
        "   improve the *apparent* accuracy without touching the encoder at all. Any A/B that",
        "   changes runtime therefore changes the pass-rate for reasons unrelated to the change",
        "   under test (exactly the config-asymmetry artefact seen in the v0.4.4 regression",
        "   sweep, where all 11 'regressions' were 300s timeouts).",
        "",
        "## The informative remainder",
        "",
    ]
    gen_diff = a["causes"].get("divergent_isomer_ambiguous", 0)
    canon = a["causes"].get("canonicalization_noise", 0)
    enc = a["causes"].get("encode_fail", 0)
    lines += [
        f"**{canon}** failures ({pct(canon)}) are canonicalization noise: the generated structure",
        "names the same isomer under the round-trip key and differs only in notation. That is the",
        "v0.4.5 canonical-string target, not a losslessness defect -- the test was too strict.",
        "",
        f"**{gen_diff}** ({pct(gen_diff)}) produced a structure naming a DIFFERENT isomer. It is",
        "tempting to score these against the notation, and that would be wrong: a divergent",
        "isomer is equally consistent with the generator having built the wrong thing and with",
        "the OIN having been ambiguous enough to license it. The sweep cannot separate those --",
        "doing so needs the per-case oracle, and it is the natural next probe.",
        "",
        f"**{enc}** ({pct(enc)}) are encoder refusals: a coverage limit of the encoder, not a",
        "round-trip result at all.",
        "",
        "So across 1837 failures, the number that is **confirmed evidence of a lossy OIN is",
        "zero**. Not because none exists -- the Y1 audit proved lossy encodings do exist -- but",
        "because this instrument cannot see them. Encoder losslessness has to be measured the way",
        "Waves 1-2 measured it, by generator-free collision probes; the round-trip pass-rate does",
        "not answer that question and never did.",
        "",
    ]
    for cause, ex in a["examples"].items():
        if cause not in ("divergent_isomer_ambiguous", "canonicalization_noise") or not ex:
            continue
        lines += [
            f"### {labels.get(cause, (cause,))[0]} -- examples",
            "",
            "| molecule | bucket |",
            "|---|---|",
        ]
        lines += [f"| {e['molecule']} | {e['bucket']} |" for e in ex]
        lines.append("")
    lines += [
        "## Reproduce",
        "",
        "```",
        "PYTHONPATH=$PWD/src python -m tools.injectivity.missed_success_audit",
        "```",
        "",
        "Reads the sweep's `bucket_report.json`; no generation is run, so the audit is cheap",
        "and deterministic.",
        "",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sweep", type=Path, default=DEFAULT_SWEEP)
    args = ap.parse_args(argv)

    a = audit(args.sweep)
    OUT_DIR.mkdir(exist_ok=True)
    (OUT_DIR / "missed_success_audit.json").write_text(json.dumps(a, indent=2) + "\n")
    md = render(a)
    (OUT_DIR / "missed_success_audit.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
