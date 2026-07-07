#!/bin/bash
# Ensure you have activated the conda environment before running this script:
# uv sync

# If running from source without installing, uncomment the following line:
# export PYTHONPATH=$PYTHONPATH:$(pwd)/src

# Create artifacts directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
ARTIFACTS_DIR="verification_artifacts_${TIMESTAMP}"
mkdir -p "$ARTIFACTS_DIR"

EXTRA_ARGS=""

show_help() {
    echo "Usage: bash tests/run_verification.sh [OPTIONS]"
    echo ""
    echo "Runs the standard integration tests and verification suite."
    echo ""
    echo "Options:"
    echo "  -h, --help            Show this help message and exit"
    echo "  --optimizer NAME      Post-FF optimizer (e.g. xtb)"
    echo "  --ff-preset PRESET    Force field preset (e.g. uff)"
    echo "  --ensemble-size N     Number of conformers to generate and optimize (default: 1)"
    echo "  --limit N             Limit the number of integration tests to run"
    echo "  --only NAME           Only run tests whose name matches this string (e.g. 'TiCat3')"
    echo "  --cpu                 Force CPU execution"
    echo ""
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) show_help; exit 0 ;;
        --optimizer) OPTIMIZER="$2"; EXTRA_ARGS="$EXTRA_ARGS --optimizer $OPTIMIZER"; shift ;;
        --ff-preset) FF_PRESET="$2"; EXTRA_ARGS="$EXTRA_ARGS --ff-preset $FF_PRESET"; shift ;;
        --ensemble-size) ENSEMBLE_SIZE="$2"; EXTRA_ARGS="$EXTRA_ARGS --ensemble-size $ENSEMBLE_SIZE"; shift ;;
        --limit) LIMIT="$2"; shift ;;
        --only) ONLY="$2"; shift ;;
        --cpu) EXTRA_ARGS="$EXTRA_ARGS --cpu" ;;
        *) echo "Unknown parameter passed: $1"; show_help; exit 1 ;;
    esac
    shift
done

echo "Artifacts will be saved to: $ARTIFACTS_DIR"

if [ -n "$LIMIT" ]; then
    LIMIT_ARG="--limit $LIMIT"
else
    LIMIT_ARG=""
fi

if [ -n "$ONLY" ]; then
    ONLY_ARG="--only \"$ONLY\""
else
    ONLY_ARG=""
fi

echo "Running Integration Tests (Real Life Examples)..."
eval uv run python tests/integration/verify_xyz_to_oin.py --output-dir \"$ARTIFACTS_DIR\" $LIMIT_ARG $ONLY_ARG \| tee \"${ARTIFACTS_DIR}/integration_log.txt\"

echo "Running Phase 1 Verification..."
eval uv run python tests/integration/verify_phase1.py $EXTRA_ARGS $ONLY_ARG \| tee \"${ARTIFACTS_DIR}/phase1_log.txt\"

echo "Running Round-Trip Verification..."
eval uv run python tests/integration/verify_roundtrip.py --output-dir \"$ARTIFACTS_DIR\" $LIMIT_ARG $ONLY_ARG $EXTRA_ARGS \| tee \"${ARTIFACTS_DIR}/roundtrip_log.txt\"

#echo "Running DG Strategy Comparison..."
#uv run python tests/integration/compare_dg_strategies.py --output-dir "$ARTIFACTS_DIR" | tee "${ARTIFACTS_DIR}/comparison_log.txt"

echo "Verification complete. Artifacts in $ARTIFACTS_DIR"
