"""Stage v0.4.17's measurements for ``/freeze-measurements``: scrub, THEN gzip, under one prefix.

WHY A STAGING STEP AT ALL
-------------------------
``tools/harvest_measurements.py`` scrubs local paths out of TEXT files on the way into the public
``measurements/`` tree. A ``.gz`` is binary to it and is copied **verbatim**. So anything large
enough to need compressing must be scrubbed BEFORE it is compressed, and its decompressed contents
checked afterwards -- the census did this inside ``attribution_table.py --gz``; this is that step
for this lane.

Every staged file but one is named ``v0417_*`` (``bucket_report_honest.md`` keeps the standard
name every baseline is diffed by) so ONE allowlist family covers the set (the harvester's
allowlist has silently dropped a release's load-bearing file four times for being one pattern too
narrow). Same-named sources from different runs -- the sweep's ``g_verdict.jsonl`` and each A/B
arm's -- get distinct names here, so nothing can overwrite anything in the release directory.

NOT staged, on purpose: ``individual_reports/`` and ``structures/`` (bulk; the sweep's rows go
through ``tools/freeze_sweep_extract.py`` instead), ``*.log``, the launcher's SCORED
``bucket_report.*`` (the honest one is the number), and ``g_verdict_control-mirror.jsonl``, which
is a copy of a file the census already froze.

    <main>/.venv/bin/python tools/v0417/freeze_stage.py          # writes <exactfold>/freeze/
"""

from __future__ import annotations

import gzip
import re
import shutil
import sys
from pathlib import Path

DATA = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
SWEEP, XF = DATA / "results-v0.4.17-sweep", DATA / "results-v0.4.17-exactfold"
STAGE = XF / "freeze"
PER_FILE_CAP = 512 * 1024  # the harvester's; a staged file over it would be silently skipped

_CHECKOUT = re.compile(r"/home/[^/\s\"']+/Documents/GitHub/([A-Za-z0-9._-]+)")
_HOME = re.compile(r"/home/[^/\s\"']+/")
_SCRATCH = re.compile(r"/tmp/claude-\d+/[^\s\"')]*")
_TMPFILE = re.compile(r"/tmp/tmp[A-Za-z0-9_]+")

#: (source, staged path under STAGE, gzip?). Two subdirectories because the precedent is two
#: releases -- ``v0.4.14`` (instruments) and ``v0.4.14-sweep`` (the baseline) -- and the harvester
#: takes one ``--from`` per release.
PLAN = [
    # --- the sweep: the new baseline of record -------------------------------------------------
    (SWEEP / "RUN.md", "sweep/v0417_sweep_RUN.md", False),
    (SWEEP / "two_numbers.txt", "sweep/v0417_sweep_two_numbers.txt", False),
    (SWEEP / "run_config.json", "sweep/v0417_sweep_run_config.json", False),
    (SWEEP / "bucket_report_honest.md", "sweep/bucket_report_honest.md", False),
    (SWEEP / "g_verdict.jsonl", "sweep/v0417_sweep_g_verdict.jsonl.gz", True),
    (SWEEP / "parseback.jsonl", "sweep/v0417_sweep_parseback.jsonl.gz", True),
    # --- L1a: the offline gate, the premise check, the live full-cohort run ----------------------
    (XF / "autofold_audit.json", "l1/v0417_autofold_audit.json.gz", True),
    (XF / "automorphism_extension_check.json", "l1/v0417_automorphism_extension_check.json", False),
    (XF / "e_selfconsistency_exact.jsonl", "l1/v0417_e_selfconsistency_exact.jsonl.gz", True),
    # --- L1b: the population (sample membership) and the generator A/B ---------------------------
    (XF / "exact_fold_movers.txt", "l1/v0417_pop_input_side_movers.txt", False),
    (XF / "generated_side_movers.txt", "l1/v0417_pop_generated_side_movers.txt", False),
    (XF / "ab_union_movers.txt", "l1/v0417_pop_ab_union_movers.txt", False),
    (XF / "generated_side_movers.jsonl", "l1/v0417_generated_side_movers.jsonl.gz", True),
    (XF / "generator_ab_report.txt", "l1/v0417_generator_ab_report.txt", False),
    (XF / "ab_gains.txt", "l1/v0417_ab_gains.txt", False),
    (XF / "ab_losses.txt", "l1/v0417_ab_losses.txt", False),
    (XF / "ab_off/bucket_report_honest.json", "l1/v0417_ab_off_bucket_report_honest.json.gz", True),
    (XF / "ab_on/bucket_report_honest.json", "l1/v0417_ab_on_bucket_report_honest.json.gz", True),
    (XF / "ab_off/g_verdict.jsonl", "l1/v0417_ab_off_g_verdict.jsonl.gz", True),
    (XF / "ab_on/g_verdict.jsonl", "l1/v0417_ab_on_g_verdict.jsonl.gz", True),
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
    total, n_scrubbed = 0, 0
    for src, name, gz in PLAN:
        if not src.exists():
            sys.exit(f"ABORT: {src} is missing -- refusing to stage a partial release")
        raw = src.read_text()
        if src.suffix == ".jsonl" and not raw.rstrip().splitlines()[-1].startswith("#DONE"):
            sys.exit(f"ABORT: {src} has no #DONE trailer -- an unfinished run must not be frozen")
        clean = scrub(raw)
        n_scrubbed += clean != raw
        if "/home/" in clean or "/tmp/claude" in clean:
            sys.exit(f"ABORT: a local path survives scrubbing in {src}")
        dest = STAGE / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if gz:
            # mtime=0 and no embedded filename: the index records a sha256 per file, and "rerun the
            # tool" must reproduce it. A default gzip header carries the wall clock.
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
        flag = (
            "  <-- OVER the harvester's per-file cap: it would be SKIPPED"
            if size > PER_FILE_CAP
            else ""
        )
        print(f"  {size:8d}  {name}{'   (scrubbed)' if clean != raw else ''}{flag}")
        if size > PER_FILE_CAP:
            sys.exit(f"ABORT: {name} is {size} bytes")
    print(
        f"staged {len(PLAN)} files, {total / 1024:.0f} KB, {n_scrubbed} needed scrubbing -> {STAGE}"
    )


if __name__ == "__main__":
    main()
