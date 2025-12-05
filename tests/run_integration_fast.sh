#!/bin/bash
# Run fast integration tests (limit 5 examples)

# Create artifacts directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
ARTIFACTS_DIR="verification_artifacts_FAST_${TIMESTAMP}"
mkdir -p "$ARTIFACTS_DIR"

echo "Artifacts will be saved to: $ARTIFACTS_DIR"

echo "Phase 1 (Single Example)..."
uv run python tests/integration/verify_phase1.py | tee "${ARTIFACTS_DIR}/phase1_log.txt"

echo "Round-Trip Verification (Limit 5)..."
uv run python tests/integration/verify_roundtrip.py --limit 5 --output-dir "$ARTIFACTS_DIR" | tee "${ARTIFACTS_DIR}/roundtrip_log.txt"

echo "XYZ to OIN Verification (Limit 5)..."
uv run python tests/integration/verify_xyz_to_oin.py --limit 5 | tee "${ARTIFACTS_DIR}/integration_log.txt"

echo "Verification complete. Artifacts in $ARTIFACTS_DIR"
