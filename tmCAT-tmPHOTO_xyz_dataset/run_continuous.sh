#!/bin/bash
# Continuously run roundtrip testing on single random un-processed molecules

cd "$(dirname "$0")"/.. || exit 1

while true; do 
    uv run python tools/test_dataset_roundtrip.py \
        --dataset-dir tmCAT-tmPHOTO_xyz_dataset \
        --quick \
        --output-dir tmCAT-tmPHOTO_xyz_dataset/20260707-results \
        --limit 1 \
        --random \
        --continue
done
