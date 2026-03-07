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

echo "Running Phase 1 Verification..."
python tests/integration/verify_phase1.py | tee "${ARTIFACTS_DIR}/phase1_log.txt"

echo "Verification complete. Artifacts in $ARTIFACTS_DIR"
