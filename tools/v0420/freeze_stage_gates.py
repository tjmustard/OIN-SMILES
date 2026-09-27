"""Stage the v0.4.20 RELEASE gate evidence for ``/freeze-measurements``: ARM 1 at 66 (shipped, the
control, and the two one-lever-at-0 attribution runs) and the ``/refreeze-goldens`` run (ARM 2 from
a FULL run). Reads only.

The release sweep is not staged again: the E2 candidate sweep IS the v0.4.20 release sweep and is
already frozen (``measurements/v0.4.19-e2-candidate-sweep/``). The census on it is staged by
``tools/census/stage_reattribution.py --prefix v0420_rcensus_`` into its own release.

Reuses ``tools/v0419/freeze_stage_release.py``'s plan runner (scrub, THEN gzip with mtime=0, check
what was WRITTEN, #DONE trailers on every .jsonl/.tsv, per-file cap) with this plan, prefix
``v0420_rarm2_`` and release ``v0.4.20-release-gates``.

    <main>/.venv/bin/python tools/v0420/freeze_stage_gates.py      # writes <arm2 dir>/freeze/
    <main>/.venv/bin/python tools/harvest_measurements.py --release v0.4.20-release-gates \\
        --from <arm2 dir>/freeze --dry-run
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v0419"))

import freeze_stage_release as fsr  # noqa: E402

ARM2 = fsr.DATA / "results-v0.4.20-arm2-refreeze"
A = "v0420_rarm2_"

fsr.STAGE = ARM2 / "freeze"
fsr.PLAN = [
    # --- ARM 1 at 66: shipped, the control (both at 0), and one lever at 0 at a time ---------------
    (ARM2 / "arm1_on.tsv", A + "arm1_on.tsv", False),
    (ARM2 / "arm1_off.tsv", A + "arm1_off_control.tsv", False),
    (ARM2 / "arm1_n2off.tsv", A + "arm1_n_valence_2_off.tsv", False),
    (ARM2 / "arm1_resoff.tsv", A + "arm1_canonical_resonance_off.tsv", False),
    # --- /refreeze-goldens: a FULL gate run, and why each re-frozen row moved -----------------------
    (ARM2 / "field2_audit_v047.jsonl", A + "field2_audit_v047.jsonl.gz", True),
    (ARM2 / "field2_audit_v049.jsonl", A + "field2_audit_v049.jsonl.gz", True),
    (ARM2 / "field2_audit_v047.txt", A + "field2_audit_v047.txt", False),
    (ARM2 / "field2_audit_v049.txt", A + "field2_audit_v049.txt", False),
    (ARM2 / "diff.txt", A + "diff_full_gate_vs_goldens.txt", False),
    (ARM2 / "splice.txt", A + "splice_dry_run.txt", False),
    (ARM2 / "on_rows.tsv", A + "rows_on_full_gate.tsv.gz", True),
    (ARM2 / "off_rows.tsv", A + "rows_off_control.tsv", False),
    (ARM2 / "verify_rows.tsv", A + "rows_verify_real_gate.tsv", False),
    (ARM2 / "comment_block.txt", A + "golden_comment_block.txt", False),
]
# DOKROM, the one row killed with the levers on only, timed alone (budget or hang)
for arm in ("on", "off", "n2off", "resoff"):
    fsr.PLAN.append((ARM2 / "dokrom_alone" / f"{arm}.tsv", A + f"dokrom_solo_{arm}.tsv", False))
# optional: the stale-cause run exists only if step 1 found stale rows
_STALE = ARM2 / "stale_cause.jsonl"
if _STALE.exists():
    fsr.PLAN.append((_STALE, A + "stale_cause.jsonl", False))

if __name__ == "__main__":
    fsr.main()
