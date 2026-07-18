#!/usr/bin/env python3
"""Emit a frozen full backlog report at each 1000-molecule milestone.

Runs alongside the v0.4.0 round-trip accumulator, after
``tools/group_v041_backlog.py`` has refreshed ``V0.4.1_ACCURACY_BACKLOG.md`` and
``case_registry.json``. Each time coverage crosses a clean multiple of the interval
(1000, 2000, 3000, ... -- never 1003, 2003), it freezes a snapshot into
``<output-dir>/reports/``:

* ``REPORT_<M>.md``            -- the tiered backlog as of that milestone
* ``case_registry_<M>.json``  -- the full per-molecule registry as of that milestone

Milestones already emitted are tracked in ``.v041_milestones.json`` so reruns and
harness restarts never double-emit. The report label is always the round milestone
number; the actual coverage at emission (a few molecules past M, since the accumulator
runs continuously and we poll periodically) is recorded inside the header.

First activation note: there are no historical per-molecule snapshots to backfill, so
identical copies for every past thousand would be misleading. On first run we seed all
past thousands as "already accounted for" and emit only the most recent milestone; every
subsequent thousand is then captured live as it is crossed. Historical *counts* remain
available per tick in ``V0.4.1_TREND.tsv``.

Strictly read-only on the round-trip harness's outputs: reads ``case_registry.json`` /
``V0.4.1_ACCURACY_BACKLOG.md`` and writes only under ``reports/`` plus its own state
file. Safe to run repeatedly alongside a live ``--continue`` accumulator.

Usage:
    python tools/milestone_report.py --output-dir <results-dir> [--interval 1000]
"""

import argparse
import json
import os
import shutil
from datetime import datetime


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--output-dir", required=True, help="Results dir with case_registry.json")
    ap.add_argument("--interval", type=int, default=1000, help="Molecules per milestone")
    args = ap.parse_args()
    output_dir = os.path.abspath(args.output_dir)
    interval = args.interval

    with open(os.path.join(output_dir, "case_registry.json")) as f:
        rows = json.load(f)
    coverage = len(rows)

    state_path = os.path.join(output_dir, ".v041_milestones.json")
    first_activation = not os.path.exists(state_path)
    emitted = set()
    if not first_activation:
        try:
            with open(state_path) as f:
                emitted = set(json.load(f).get("emitted", []))
        except Exception:
            emitted = set()

    milestones = list(range(interval, coverage + 1, interval))  # 1000, 2000, ...
    if first_activation and milestones:
        # No historical per-molecule snapshots exist to backfill past thousands with;
        # seed everything but the most recent as accounted-for, emit the latest now.
        emitted = set(milestones[:-1])

    to_emit = [m for m in milestones if m not in emitted]

    if not to_emit:
        with open(state_path, "w") as f:
            json.dump({"emitted": sorted(emitted)}, f)
        print(f"milestone: coverage {coverage}, none new")
        return

    reports_dir = os.path.join(output_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)

    backlog_path = os.path.join(output_dir, "V0.4.1_ACCURACY_BACKLOG.md")
    backlog = ""
    if os.path.exists(backlog_path):
        with open(backlog_path) as f:
            backlog = f.read()
    reg_path = os.path.join(output_dir, "case_registry.json")

    ts = datetime.now().isoformat(timespec="seconds")
    for m in to_emit:
        header = (
            f"# v0.4.1 Milestone Report — {m} molecules\n\n"
            f"Milestone: **{m}** (interval {interval}). "
            f"Coverage at emission: **{coverage}**. Generated {ts}.\n\n"
            "Frozen snapshot of the tiered accuracy backlog at this milestone. "
            "See `../V0.4.1_TREND.tsv` for the continuous per-tick history and "
            f"`case_registry_{m}.json` for the full per-molecule registry.\n\n"
            "---\n\n"
        )
        with open(os.path.join(reports_dir, f"REPORT_{m}.md"), "w") as f:
            f.write(header + backlog)
        try:
            shutil.copyfile(reg_path, os.path.join(reports_dir, f"case_registry_{m}.json"))
        except Exception:
            pass
        emitted.add(m)

    with open(state_path, "w") as f:
        json.dump({"emitted": sorted(emitted)}, f)

    print(f"milestone: coverage {coverage}, emitted " + ", ".join(f"REPORT_{m}" for m in to_emit))


if __name__ == "__main__":
    main()
