#!/bin/bash
# Run all verification tests including tmQM examples

# Create artifacts directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
ARTIFACTS_DIR="verification_artifacts_ALL_${TIMESTAMP}"
mkdir -p "$ARTIFACTS_DIR"

EXTRA_ARGS=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --optimizer) OPTIMIZER="$2"; EXTRA_ARGS="$EXTRA_ARGS --optimizer $OPTIMIZER"; shift ;;
        --ff-preset) FF_PRESET="$2"; EXTRA_ARGS="$EXTRA_ARGS --ff-preset $FF_PRESET"; shift ;;
        --ensemble-size) ENSEMBLE_SIZE="$2"; EXTRA_ARGS="$EXTRA_ARGS --ensemble-size $ENSEMBLE_SIZE"; shift ;;
        --cpu) EXTRA_ARGS="$EXTRA_ARGS --cpu" ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

echo "Artifacts will be saved to: $ARTIFACTS_DIR"

echo "Running Integration Tests (ALL, including tmQM)..."
uv run python tests/integration/verify_xyz_to_oin.py --include-tmqm --output-dir "$ARTIFACTS_DIR" | tee "${ARTIFACTS_DIR}/integration_log.txt"

echo "Running Phase 1 Verification..."
uv run python tests/integration/verify_phase1.py $EXTRA_ARGS | tee "${ARTIFACTS_DIR}/phase1_log.txt"

echo "Running Round-Trip Verification (ALL)..."
uv run python tests/integration/verify_roundtrip.py --output-dir "$ARTIFACTS_DIR" $EXTRA_ARGS | tee "${ARTIFACTS_DIR}/roundtrip_log.txt"

#echo "Running DG Strategy Comparison (ALL)..."
#uv run python tests/integration/compare_dg_strategies.py --include-tmqm --output-dir "$ARTIFACTS_DIR" | tee "${ARTIFACTS_DIR}/comparison_log.txt"

echo "Verification complete. Artifacts in $ARTIFACTS_DIR"
