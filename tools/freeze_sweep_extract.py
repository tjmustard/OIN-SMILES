#!/usr/bin/env python
"""Freeze a whole 5000-molecule sweep into ONE tracked, git-durable file. (v0.4.16)

WHY THIS EXISTS
===============
``tmCAT-tmPHOTO_xyz_dataset/`` is gitignored **in its entirety**, so every sweep this project has
ever run lives exactly one ``rm -rf`` -- or one ``Kulik_TMC_Dataset`` branch switch, which silently
deletes 26,232 files -- from oblivion. ``measurements/`` was built to stop that, and it works, but
it is an ALLOWLIST of *summary* artifacts: bucket reports, A/B tallies, population lists. It
deliberately does not carry the per-molecule rows.

That gap is real. Every re-analysis this project performs starts from
``individual_reports/*.json``:

* ``roundtrip_bucket_report.py``  -> the six-bucket table and the honest/scored delta
* ``attach_class_audit.py``       -> the attachment split AND v0.4.16's key-component taxonomy
* every "re-derive the gap decomposition on the fresh sweep" step in every CLOSEOUT

If the reports are gone, none of that can be re-run and a later release cannot check whether a
block moved -- it can only re-quote prose. v0.4.11's mirror-audit JSON was lost exactly this way.

WHY NOT JUST COMMIT THE DIRECTORY
=================================
Measured on ``results-v0.4.14-sweep``:

    whole directory            268 MB
      individual_reports/       20 MB   (5000 files)
      structures/               42 MB   (9736 XYZ files)
    this extract, gzipped     0.26 MB

**1000x smaller than the directory and 77x smaller than the reports alone**, because the bulk is
per-file JSON framing and the geometries. Committing 268 MB per release into a repo whose ``.git``
is already 787 MB is not a durability strategy, it is a different failure.

WHAT IT CARRIES, AND WHAT IT DOES NOT
=====================================
Carried: exactly the fields the analyses above read -- ``molecule``, ``status``, ``error``,
``smiles_1``, ``smiles_2`` (scored), ``smiles_2_indep`` (honest), ``indep_key_match``,
``metrics.elapsed_s``, and the ``coordination`` verdict (``intact`` / ``boundary_only``).

🔴 **NOT carried: the geometries.** ``structures/`` is 42 MB of XYZ and this extract cannot
reproduce any analysis that re-reads coordinates -- a mirror audit, a clash re-score, an
attachment re-computation from atoms. Those need the original directory. **The extract is a
re-analysis substrate, not an archive**, and a release that quotes it must not claim otherwise.

DETERMINISM
===========
Rows sorted by molecule, compact separators, gzip ``mtime=0``. Re-running on the same sweep
produces a BYTE-IDENTICAL file, so a diff means the sweep changed rather than the tool.

Usage
=====
    V=$PWD/.venv/bin/python
    $V tools/freeze_sweep_extract.py --results-dir <sweep> --release v0.4.16 [--check]
"""

from __future__ import annotations

import argparse
import glob
import gzip
import hashlib
import json
import os
import sys
from pathlib import Path

#: Read from ``metrics``, which is NESTED. A top-level read of ``elapsed_s`` silently yields 0 --
#: verified 5000/5000 -- and it is a SUM over up to three SIGKILLed harness attempts.
_FIELDS = (
    "molecule",
    "status",
    "error",
    "smiles_1",
    "smiles_2",
    "smiles_2_indep",
    "indep_key_match",
)


def extract_row(report: dict) -> dict:
    coord = report.get("coordination")
    coord = coord if isinstance(coord, dict) else {}
    row = {k: report.get(k) for k in _FIELDS}
    row["elapsed_s"] = (report.get("metrics") or {}).get("elapsed_s")
    row["intact"] = coord.get("intact")
    row["boundary_only"] = coord.get("boundary_only")
    return row


def build(results_dir: Path) -> tuple[bytes, int]:
    paths = sorted(glob.glob(str(results_dir / "individual_reports" / "*.json")))
    if not paths:
        sys.exit(f"🔴 REFUSING: no individual_reports/*.json under {results_dir}")
    rows = []
    for path in paths:
        try:
            rows.append(extract_row(json.loads(Path(path).read_text())))
        except (json.JSONDecodeError, OSError) as exc:
            # An unreadable report is DATA -- recorded as a row with its error, never dropped.
            # Dropping it would silently shrink the denominator, which is how a rate stops being
            # reproducible.
            rows.append({"molecule": Path(path).stem, "status": None, "error": f"unreadable:{exc}"})
    rows.sort(key=lambda r: r.get("molecule") or "")
    raw = "\n".join(json.dumps(r, separators=(",", ":"), sort_keys=True) for r in rows)
    buf = gzip.compress(raw.encode(), mtime=0)
    return buf, len(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--release", required=True, help="e.g. v0.4.16")
    ap.add_argument("--expect", type=int, default=None, help="refuse unless the row count matches")
    ap.add_argument(
        "--check", action="store_true", help="verify an existing extract instead of writing"
    )
    args = ap.parse_args()

    # Always the MAIN checkout: a worktree's files die on `git worktree remove`, which is how a
    # previous release lost its measurements.
    main_checkout = (
        Path(os.popen("git rev-parse --git-common-dir").read().strip() or ".git").resolve().parent
    )
    out_dir = main_checkout / "measurements" / args.release
    out = out_dir / f"sweep_extract_{args.release}.jsonl.gz"

    payload, n = build(args.results_dir.resolve())
    digest = hashlib.sha256(payload).hexdigest()

    # 🔴 COMPLETENESS IS A ROW COUNT, NOT AN EXIT STATUS. A sweep that died mid-run still leaves a
    # plausible, short directory and this tool would happily freeze it.
    if args.expect is not None and n != args.expect:
        sys.exit(f"🔴 REFUSING: {n} rows, expected {args.expect}. A short sweep is not a sweep.")

    if args.check:
        if not out.exists():
            sys.exit(f"🔴 no extract at {out}")
        same = hashlib.sha256(out.read_bytes()).hexdigest() == digest
        print(f"{'✅ byte-identical' if same else '🔴 DIFFERS'}: {out}")
        return 0 if same else 1

    out_dir.mkdir(parents=True, exist_ok=True)
    out.write_bytes(payload)
    print(f"wrote {out}")
    print(f"  rows        {n}")
    print(f"  size        {len(payload) / 1048576:.2f} MB gzipped")
    print(f"  sha256      {digest}")
    print(f"  source      {args.results_dir.name}")
    print(
        "\n⚠ This is a RE-ANALYSIS SUBSTRATE, not an archive. It carries no geometries, so a\n"
        "  mirror audit / clash re-score / attachment recomputation still needs the original dir."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
