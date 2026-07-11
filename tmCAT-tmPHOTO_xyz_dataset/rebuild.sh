#!/bin/bash
# Rebuild the summary_roundtrip.json file from the individual reports.

cd "$(dirname "$0")"/.. || exit 1
uv run python tools/rebuild_summary.py --output-dir tmCAT-tmPHOTO_xyz_dataset/20260707-results
