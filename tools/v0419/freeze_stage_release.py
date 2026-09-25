"""Stage the v0.4.19 RELEASE evidence for ``/freeze-measurements``: the candidate sweep (= the
release sweep), its post-processing under the levers, the row-level comparison with the projection,
and the ``/refreeze-goldens`` run (ARM 1 at 65 and ARM 2 from a FULL run). Scrub, THEN gzip
(mtime=0); two prefixes, ONE release (``v0.4.19-release-sweep``), because the harvester rebuilds a
release's README index from the current picks only. Reads only.

    <main>/.venv/bin/python tools/v0419/freeze_stage_release.py        # writes <sweep>/freeze/
    <main>/.venv/bin/python tools/harvest_measurements.py --release v0.4.19-release-sweep \\
        --from <sweep>/freeze --dry-run

NOT staged: ``individual_reports/`` and ``structures/`` (bulk), logs, ``.shippedreader/``.
"""

from __future__ import annotations

import gzip
import re
import shutil
import sys
from pathlib import Path

DATA = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
SWEEP = DATA / "results-v0.4.19-candidate-sweep"
ARM2 = DATA / "results-v0.4.19-arm2-refreeze"
MOVERS = DATA / "results-v0.4.19-movers14"
STAGE = SWEEP / "freeze"
PER_FILE_CAP = 512 * 1024
S, A = "v0419_rsweep_", "v0419_rarm2_"
_CHECKOUT = re.compile(r"/home/[^/\s\"']+/Documents/GitHub/([A-Za-z0-9._-]+)")
_HOME = re.compile(r"/home/[^/\s\"']+/")
_SCRATCH = re.compile(r"/tmp/claude-\d+/[^\s\"')]*")
_TMPFILE = re.compile(r"/tmp/tmp[A-Za-z0-9_]+")

PLAN = [
    # --- the sweep: the new baseline of record ---------------------------------------------------
    (SWEEP / "run_config.json", S + "run_config.json", False),
    (SWEEP / "two_numbers.txt", S + "two_numbers.txt", False),
    (SWEEP / "bucket_report_honest.md", S + "bucket_report_honest.md", False),
    (SWEEP / "candidate_vs_projection.txt", S + "vs_projection.txt", False),
    (SWEEP / "candidate_vs_projection.json", S + "vs_projection.json", False),
    (SWEEP / "unchanged_string_structure_movers.json", S + "unchanged_string_movers.json", False),
    # per-molecule honest buckets + the ruler + parse-back UNDER THE LEVERS: gunzip the three and
    # tools/v0417/sweep_two_numbers.two_numbers() re-derives both headline numbers
    (SWEEP / "bucket_report_honest.json", S + "bucket_report_honest.json.gz", True),
    (SWEEP / "g_verdict.jsonl", S + "g_verdict.jsonl.gz", True),
    (SWEEP / "parseback.jsonl", S + "parseback.jsonl.gz", True),
    # the same strings read by the SHIPPED (pre-promotion) reader: the OIN_H_FAITHFUL coupling
    (SWEEP / "parseback_shippedreader.jsonl", S + "parseback_shippedreader.jsonl.gz", True),
    (SWEEP / "parseback_gen.jsonl", S + "parseback_gen.jsonl.gz", True),
    # --- ARM 1 at 65: shipped and the control ------------------------------------------------------
    (ARM2 / "arm1_on.tsv", A + "arm1_on.tsv", False),
    (ARM2 / "arm1_off.tsv", A + "arm1_off_control.tsv", False),
    # --- /refreeze-goldens: a FULL gate run, 26 of 425 rows, and why --------------------------------
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


def scrub(text: str) -> str:
    text = _CHECKOUT.sub(r"<checkout:\1>", text)
    text = _SCRATCH.sub("<SCRATCH>", text)
    text = _TMPFILE.sub("<TMPFILE>", text)
    return _HOME.sub("<HOME>/", text)


def main():
    if STAGE.exists():
        shutil.rmtree(STAGE)  # a stale staged file would be harvested as if it were current
    STAGE.mkdir(parents=True)
    total = 0
    for src, name, gz in PLAN:
        if not src.exists():
            sys.exit(f"ABORT: {src} is missing -- refusing to stage a partial release")
        raw = src.read_text()
        if src.suffix in (".jsonl", ".tsv") and not (
            raw.rstrip().splitlines()[-1].startswith("#DONE")
        ):
            sys.exit(f"ABORT: {src} has no #DONE trailer -- an unfinished run must not be frozen")
        clean = scrub(raw)
        if "/home/" in clean or "/tmp/claude" in clean:
            sys.exit(f"ABORT: a local path survives scrubbing in {src}")
        dest = STAGE / name
        if gz:
            with open(dest, "wb") as raw_fh:
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, mtime=0) as fh:
                    fh.write(clean.encode())
            back = gzip.open(dest, "rt").read()  # check what was WRITTEN, not what was intended
            if back != clean or "/home/" in back or "/tmp/claude" in back:
                sys.exit(f"ABORT: {dest} does not decompress to the scrubbed text")
        else:
            dest.write_text(clean)
        size = dest.stat().st_size
        total += size
        print(f"  {size:8d}  {name}{'   (scrubbed)' if clean != raw else ''}")
        if size > PER_FILE_CAP:
            sys.exit(
                f"ABORT: {name} is {size} bytes -- over the harvester's cap, it would be SKIPPED"
            )
    print(f"staged {len(PLAN)} files, {total / 1024:.0f} KB -> {STAGE}")


if __name__ == "__main__":
    main()
