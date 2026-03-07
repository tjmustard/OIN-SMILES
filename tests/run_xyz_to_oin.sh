#!/bin/bash
# Run fast integration tests (limit 22 examples)

# Create artifacts directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
ARTIFACTS_DIR="verification_artifacts_OIN_${TIMESTAMP}"
mkdir -p "$ARTIFACTS_DIR"

echo "Artifacts will be saved to: $ARTIFACTS_DIR"

echo "XYZ to OIN Verification (Limit 22)..."
uv run python tests/integration/verify_xyz_to_oin.py --limit 22 | tee "${ARTIFACTS_DIR}/integration_log.txt"

echo "Verification complete. Artifacts in $ARTIFACTS_DIR"
