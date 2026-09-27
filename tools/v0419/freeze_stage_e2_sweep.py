"""Stage the E2 lane's CANDIDATE sweep (v0.4.19 defaults + OIN_N_VALENCE_2 + OIN_CANONICAL_RESONANCE)
for ``/freeze-measurements``: the sweep's post-processing under the levers and the row-level
comparison with the e2f projection. Reads only.

Reuses ``freeze_stage_release.py``'s plan runner (scrub, THEN gzip with mtime=0, check what was
WRITTEN, #DONE trailers on every .jsonl/.tsv, per-file cap) with this sweep's plan, prefix
``v0419_e2sweep_`` and release ``v0.4.19-e2-candidate-sweep``.

    <main>/.venv/bin/python tools/v0419/freeze_stage_e2_sweep.py      # writes <sweep>/freeze/
    <main>/.venv/bin/python tools/harvest_measurements.py --release v0.4.19-e2-candidate-sweep \\
        --from <sweep>/freeze --dry-run

NOT staged: ``individual_reports/`` and ``structures/`` (bulk), logs, ``.shippedreader/``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import freeze_stage_release as fsr  # noqa: E402

SWEEP = fsr.DATA / "results-v0.4.19-e2-candidate-sweep"
S = "v0419_e2sweep_"

fsr.STAGE = SWEEP / "freeze"
fsr.PLAN = [
    (SWEEP / "run_config.json", S + "run_config.json", False),
    (SWEEP / "two_numbers.txt", S + "two_numbers.txt", False),
    (SWEEP / "bucket_report_honest.md", S + "bucket_report_honest.md", False),
    (SWEEP / "candidate_vs_projection.txt", S + "vs_projection.txt", False),
    (SWEEP / "candidate_vs_projection.json", S + "vs_projection.json", False),
    # per-molecule honest buckets + the ruler + parse-back UNDER THE LEVERS: gunzip the three and
    # tools/v0417/sweep_two_numbers.two_numbers() re-derives both headline numbers
    (SWEEP / "bucket_report_honest.json", S + "bucket_report_honest.json.gz", True),
    (SWEEP / "g_verdict.jsonl", S + "g_verdict.jsonl.gz", True),
    (SWEEP / "parseback.jsonl", S + "parseback.jsonl.gz", True),
    (SWEEP / "parseback_shippedreader.jsonl", S + "parseback_shippedreader.jsonl.gz", True),
    (SWEEP / "parseback_gen.jsonl", S + "parseback_gen.jsonl.gz", True),
]
# optional: the generator-reach table, written by hand after the sweep if it was computed
_REACH = SWEEP / "unchanged_string_structure_movers.json"
if _REACH.exists():
    fsr.PLAN.append((_REACH, S + "unchanged_string_movers.json", False))

if __name__ == "__main__":
    fsr.main()
