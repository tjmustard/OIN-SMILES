"""Stage the v0.4.18 PROMOTION evidence for ``/freeze-measurements``: the sweep of record and the
``/refreeze-goldens`` run. Scrub, THEN gzip (mtime=0), two prefixes, one release. Reads only.

Same discipline as ``tools/v0418/freeze_stage.py`` (the lane's A/B, release ``v0.4.18-l2``). This is
a SEPARATE release, ``v0.4.18-sweep``, on purpose: the harvester rebuilds a release's README index
from the CURRENT PICKS ONLY, so adding files to an existing release erases the earlier ones from it.

NOT staged: ``individual_reports/`` and ``structures/`` (bulk -- the sweep's per-molecule rows go
through ``tools/freeze_sweep_extract.py --release v0.4.18-sweep``), logs, and the launcher's default
SCORED bucket report (the number is the HONEST one).

    <main>/.venv/bin/python tools/v0418/freeze_stage_sweep.py        # writes <sweep>/freeze/
"""

from __future__ import annotations

import gzip
import re
import shutil
import sys
from pathlib import Path

DATA = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
SWEEP = DATA / "results-v0.4.18-sweep"
ARM2 = DATA / "results-v0.4.18-arm2-refreeze"
STAGE = SWEEP / "freeze"
PER_FILE_CAP = 512 * 1024
S, A = "v0418_sweep_", "v0418_arm2_"

_CHECKOUT = re.compile(r"/home/[^/\s\"']+/Documents/GitHub/([A-Za-z0-9._-]+)")
_HOME = re.compile(r"/home/[^/\s\"']+/")
_SCRATCH = re.compile(r"/tmp/claude-\d+/[^\s\"')]*")
_TMPFILE = re.compile(r"/tmp/tmp[A-Za-z0-9_]+")

PLAN = [
    # --- the sweep: the new baseline of record ---------------------------------------------------
    (SWEEP / "RUN.md", S + "RUN.md", False),
    (SWEEP / "run_config.json", S + "run_config.json", False),
    (SWEEP / "two_numbers.txt", S + "two_numbers.txt", False),
    (SWEEP / "sweep_vs_ab.txt", S + "vs_ab.txt", False),
    (SWEEP / "sweep_vs_ab.json", S + "vs_ab.json", False),
    (SWEEP / "bucket_report_honest.md", S + "bucket_report_honest.md", False),
    # per-molecule honest buckets: with the two files below this re-derives BOTH headline numbers
    # from the frozen tree alone (tools/v0417/sweep_two_numbers.two_numbers on the gunzipped three)
    (SWEEP / "bucket_report_honest.json", S + "bucket_report_honest.json.gz", True),
    (SWEEP / "g_verdict.jsonl", S + "g_verdict.jsonl.gz", True),
    (SWEEP / "parseback.jsonl", S + "parseback.jsonl.gz", True),
    # --- /refreeze-goldens: a FULL gate run, 0 of 425 rows, and why ---------------------------------
    (ARM2 / "field2_audit_v047.jsonl", A + "field2_audit_v047.jsonl.gz", True),
    (ARM2 / "field2_audit_v049.jsonl", A + "field2_audit_v049.jsonl.gz", True),
    (ARM2 / "field2_audit_v047.txt", A + "field2_audit_v047.txt", False),
    (ARM2 / "field2_audit_v049.txt", A + "field2_audit_v049.txt", False),
    (ARM2 / "diff.txt", A + "diff_full_gate_vs_goldens.txt", False),
    (ARM2 / "on_rows.tsv", A + "rows_on_full_gate.tsv.gz", True),
    (ARM2 / "off_rows.tsv", A + "rows_off_control.tsv", False),
    (ARM2 / "lever_fired.txt", A + "lever_fired_but_unseen.txt", False),
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
