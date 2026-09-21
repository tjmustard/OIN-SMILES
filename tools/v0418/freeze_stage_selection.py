"""Stage the SELECTION lane (OIN_ETA_RETARGET) for ``/freeze-measurements``. Reads only.

Reuses ``freeze_stage.py``'s scrub -> gzip(mtime=0) -> verify-what-was-written machinery with this
lane's directory, prefix and arms. Its OWN release, ``v0.4.18-selection``: the harvester rebuilds a
release's index from the current picks only.

⚠ ``ab_off`` here was BUILT from the v0.4.18 sweep of record (symlinks), not run -- so the
"noise floor" the verifier prints (OFF == sweep on every comparable structure) is true by
construction and measures nothing. The determinism evidence is in v0.4.18-l2 and v0.4.18-sweep.

    <main>/.venv/bin/python tools/v0418/freeze_stage_selection.py
    <main>/.venv/bin/python tools/v0418/freeze_stage_selection.py --verify <main>/measurements/v0.4.18-selection
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))

import freeze_stage as fs  # noqa: E402

fs.R = fs.DATA / "results-v0.4.18-eta-selection"
fs.RECORD = fs.DATA / "results-v0.4.18-sweep"
fs.STAGE = fs.R / "freeze"
fs.P = "v0418_sel_"
fs.ARMS = ("off", "retarget")
EXPECTED = {"retarget": (49, 4, 43, 5)}  # self gains, losses; verified gains, losses
TEXT = [
    "selection_sim.txt",
    "selection_sim.json",
    "retarget_probe.txt",
    "probe5_classes.json",
    "eta_ab_report_retarget.txt",
    "eta_ab_report_retarget.json",
]


def stage():
    if fs.STAGE.exists():
        shutil.rmtree(fs.STAGE)  # a stale staged file would be harvested as if it were current
    fs.STAGE.mkdir(parents=True)
    total = n = 0
    for f in TEXT:
        if not (fs.R / f).exists():
            sys.exit(f"ABORT: {fs.R / f} is missing -- refusing to stage a partial release")
        total += fs._write(fs.P + f, (fs.R / f).read_text(), False)
        n += 1
    for a in fs.ARMS:
        d = fs.R / f"ab_{a}"
        gv = (d / "g_verdict.jsonl").read_text()
        if not gv.rstrip().splitlines()[-1].startswith("#DONE"):
            sys.exit(f"ABORT: {d}/g_verdict.jsonl has no #DONE trailer")
        total += fs._write(f"{fs.P}ab_{a}_g_verdict.jsonl.gz", gv, True)
        total += fs._write(
            f"{fs.P}ab_{a}_bucket_report_honest.json.gz",
            (d / "bucket_report_honest.json").read_text(),
            True,
        )
        n += 2
    total += fs._write(fs.P + "ab_rows.tsv.gz", fs._ab_rows(), True)
    commits = "".join(
        f"{a}\t{(fs.R / f'ab_{a}' / 'AB_COMMIT').read_text().strip()}\n" for a in fs.ARMS
    )
    total += fs._write(fs.P + "ab_commits.tsv", "arm\tlaunched_from\n" + commits, False)
    print(f"staged {n + 2} files, {total / 1024:.0f} KB -> {fs.STAGE}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--verify", type=Path, help="a frozen release dir; re-derive the headline")
    args = ap.parse_args()
    if args.verify:
        fs.verify(args.verify, expected=EXPECTED, noise=(1072, 1072))
    else:
        stage()


if __name__ == "__main__":
    main()
