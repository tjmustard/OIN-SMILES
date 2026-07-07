#!/bin/bash
# Run all verification tests including tmQM examples

# Create artifacts directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
ARTIFACTS_DIR="verification_artifacts_ALL_${TIMESTAMP}"
mkdir -p "$ARTIFACTS_DIR"

EXTRA_ARGS=""

show_help() {
    echo "Usage: bash tests/run_verification_ALL.sh [OPTIONS]"
    echo ""
    echo "Runs all verification tests including tmQM examples."
    echo ""
    echo "Options:"
    echo "  -h, --help            Show this help message and exit"
    echo "  --optimizer NAME      Post-FF optimizer (e.g. mace-omol-0-extra-large-1024)"
    echo "  --ff-preset PRESET    Force field preset (e.g. uff)"
    echo "  --ensemble-size N     Number of conformers to generate and optimize (default: 1)"
    echo "  --only NAME           Only run tests whose name matches this string"
    echo "  --cpu                 Force CPU execution"
    echo ""
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) show_help; exit 0 ;;
        --optimizer) OPTIMIZER="$2"; EXTRA_ARGS="$EXTRA_ARGS --optimizer $OPTIMIZER"; shift ;;
        --ff-preset) FF_PRESET="$2"; EXTRA_ARGS="$EXTRA_ARGS --ff-preset $FF_PRESET"; shift ;;
        --ensemble-size) ENSEMBLE_SIZE="$2"; EXTRA_ARGS="$EXTRA_ARGS --ensemble-size $ENSEMBLE_SIZE"; shift ;;
        --only) ONLY="$2"; shift ;;
        --cpu) EXTRA_ARGS="$EXTRA_ARGS --cpu" ;;
        *) echo "Unknown parameter passed: $1"; show_help; exit 1 ;;
    esac
    shift
done

echo "Artifacts will be saved to: $ARTIFACTS_DIR"

if [ -n "$ONLY" ]; then
    ONLY_ARG="--only \"$ONLY\""
else
    ONLY_ARG=""
fi

echo "Running Integration Tests (ALL, including tmQM)..."
eval uv run python tests/integration/verify_xyz_to_oin.py --include-tmqm --output-dir \"$ARTIFACTS_DIR\" $ONLY_ARG \| tee \"${ARTIFACTS_DIR}/integration_log.txt\"

echo "Running Phase 1 Verification..."
eval uv run python tests/integration/verify_phase1.py $EXTRA_ARGS $ONLY_ARG \| tee \"${ARTIFACTS_DIR}/phase1_log.txt\"

echo "Running Round-Trip Verification (ALL)..."
eval uv run python tests/integration/verify_roundtrip.py --output-dir \"$ARTIFACTS_DIR\" $ONLY_ARG $EXTRA_ARGS \| tee \"${ARTIFACTS_DIR}/roundtrip_log.txt\"

#echo "Running DG Strategy Comparison (ALL)..."
#uv run python tests/integration/compare_dg_strategies.py --include-tmqm --output-dir "$ARTIFACTS_DIR" | tee "${ARTIFACTS_DIR}/comparison_log.txt"

echo "Verification complete. Artifacts in $ARTIFACTS_DIR"
