"""Stage the re-run census for ``/freeze-measurements``: scrub, THEN gzip, one prefix. Reads only.

Same discipline as ``tools/v0417/freeze_stage.py`` and for the same reasons: the harvester copies a
``.gz`` VERBATIM into a PUBLIC tree, so text is scrubbed before it is compressed and the
decompressed bytes are checked afterwards; gzip is written with ``mtime=0`` so a re-stage is
byte-identical; and every file carries one prefix (``v0418_census_``) so one allowlist family and
one PROVENANCE block cover the set. The harvester keys provenance on the FILENAME alone, and the
census's own ``attribution_table.tsv.gz`` is described there as the v0.4.14 sweep's -- reusing that
name for a different sweep's table would publish a false provenance line.

NOT staged: the v0.4.17 sweep's ``g_verdict`` / ``parseback`` (frozen in ``measurements/
v0.4.17-sweep/``), the exact-fold ``e_selfconsistency`` (``measurements/v0.4.17/``), and the two
input-only census files (``measurements/v0.4.17-census/``). The README says where each one is.

    <main>/.venv/bin/python tools/census/stage_reattribution.py       # writes <reattr>/freeze/
"""

from __future__ import annotations

import gzip
import re
import shutil
import sys
from pathlib import Path

DATA = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
# v0.4.18 release: the same layout under a different directory and prefix
#   stage_reattribution.py --results results-v0.4.18-release-census --prefix v0418_rcensus_
_a = sys.argv[1:]
R = DATA / (_a[_a.index("--results") + 1] if "--results" in _a else "results-v0.4.17-reattribution")
STAGE = R / "freeze"
PER_FILE_CAP = 512 * 1024
P = _a[_a.index("--prefix") + 1] if "--prefix" in _a else "v0418_census_"

_CHECKOUT = re.compile(r"/home/[^/\s\"']+/Documents/GitHub/([A-Za-z0-9._-]+)")
_HOME = re.compile(r"/home/[^/\s\"']+/")
_SCRATCH = re.compile(r"/tmp/claude-\d+/[^\s\"')]*")
_TMPFILE = re.compile(r"/tmp/tmp[A-Za-z0-9_]+")

PLAN = [
    (R / "attribution_table.tsv", P + "attribution_table.tsv.gz", True),
    (R / "attribution_summary.json", P + "attribution_summary.json", False),
    (R / "attribution.txt", P + "attribution.txt", False),
    (R / "control" / "attribution_control.txt", P + "control_on_the_record_sweep.txt", False),
    (R / "reattribution_diff.txt", P + "reattribution_diff.txt", False),
    (R / "reattribution_diff.json", P + "reattribution_diff.json", False),
    (R / "parseback_gen.jsonl", P + "parseback_gen.jsonl.gz", True),
    (R / "pflags.jsonl", P + "pflags.jsonl.gz", True),
    (R / "attach_class_audit.json", P + "attach_class_audit.json.gz", True),
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
        if src.suffix == ".jsonl" and not raw.rstrip().splitlines()[-1].startswith("#DONE"):
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
