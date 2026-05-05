#!/bin/bash
# Run fast integration tests (limit 5 examples)

# Create artifacts directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
ARTIFACTS_DIR="verification_artifacts_FAST_${TIMESTAMP}"
mkdir -p "$ARTIFACTS_DIR"

echo "Artifacts will be saved to: $ARTIFACTS_DIR"

echo "XYZ to OIN Verification (Limit 4)..."
uv run python tests/integration/verify_xyz_to_oin.py --limit 4 --output-dir "$ARTIFACTS_DIR" | tee "${ARTIFACTS_DIR}/integration_log.txt"

echo "Phase 1 (Single Example)..."
uv run python tests/integration/verify_phase1.py | tee "${ARTIFACTS_DIR}/phase1_log.txt"

echo "Round-Trip Verification (Limit 4)..."
uv run python tests/integration/verify_roundtrip.py --limit 4 --output-dir "$ARTIFACTS_DIR" | tee "${ARTIFACTS_DIR}/roundtrip_log.txt"

#echo "DG Strategy Comparison (Limit 4)..."
#uv run python tests/integration/compare_dg_strategies.py --limit 4 --output-dir "$ARTIFACTS_DIR" | tee "${ARTIFACTS_DIR}/comparison_log.txt"

echo "Verification complete. Artifacts in $ARTIFACTS_DIR"
