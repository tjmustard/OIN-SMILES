#!/bin/bash
# Ensure you have activated the conda environment before running this script:
# uv sync

# If running from source without installing, uncomment the following line:
# export PYTHONPATH=$PYTHONPATH:$(pwd)/src

# Create artifacts directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
ARTIFACTS_DIR="verification_artifacts_${TIMESTAMP}"
mkdir -p "$ARTIFACTS_DIR"

echo "Artifacts will be saved to: $ARTIFACTS_DIR"

echo "Running Integration Tests (Real Life Examples)..."
uv run python tests/integration/verify_xyz_to_oin.py --output-dir "$ARTIFACTS_DIR" | tee "${ARTIFACTS_DIR}/integration_log.txt"

echo "Running Phase 1 Verification..."
uv run python tests/integration/verify_phase1.py | tee "${ARTIFACTS_DIR}/phase1_log.txt"

echo "Running Round-Trip Verification..."
uv run python tests/integration/verify_roundtrip.py --output-dir "$ARTIFACTS_DIR" | tee "${ARTIFACTS_DIR}/roundtrip_log.txt"

#echo "Running DG Strategy Comparison..."
#uv run python tests/integration/compare_dg_strategies.py --output-dir "$ARTIFACTS_DIR" | tee "${ARTIFACTS_DIR}/comparison_log.txt"

echo "Verification complete. Artifacts in $ARTIFACTS_DIR"
