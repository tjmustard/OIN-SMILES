#!/bin/bash
#
# Download the MACE foundation-model weights used by the opt-in MACE optimizer
# (`OIN3DGenerator(optimizer="mace-omol-0-extra-large-1024")`) and register the
# path in a repo-root .env file.
#
# The MACE-omol-0 "extra large" checkpoint is published as a public asset on the
# ACEsuit/mace-foundations GitHub release, so it can be fetched without any login
# (unlike the Hugging-Face-gated facebook/OMol25 checkpoints).
#
# Usage:
#   bash tools/install_mace_weights.sh
#
# Idempotent: re-running skips the download when the file is already present.

set -euo pipefail

MODEL_NAME="MACE-omol-0-extra-large-1024.model"
URL="https://github.com/ACEsuit/mace-foundations/releases/download/mace_omol_0/${MODEL_NAME}"
ENV_VAR="MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH"
MIN_BYTES=$((100 * 1024 * 1024)) # sanity floor: real checkpoint is ~400 MB

# Resolve the repo root (this script lives in <repo>/tools/).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MODEL_DIR="$REPO_ROOT/models/mace"
MODEL_PATH="$MODEL_DIR/$MODEL_NAME"
REL_PATH="models/mace/$MODEL_NAME"
ENV_FILE="$REPO_ROOT/.env"

file_size() { # portable stat
    wc -c <"$1" | tr -d ' '
}

mkdir -p "$MODEL_DIR"

if [ -f "$MODEL_PATH" ] && [ "$(file_size "$MODEL_PATH")" -ge "$MIN_BYTES" ]; then
    echo "✓ $MODEL_NAME already present ($(($(file_size "$MODEL_PATH") / 1024 / 1024)) MB) — skipping download."
else
    echo "Downloading $MODEL_NAME (~400 MB) from ACEsuit/mace-foundations ..."
    TMP="$(mktemp)"
    if command -v curl >/dev/null 2>&1; then
        curl -fL --progress-bar "$URL" -o "$TMP"
    elif command -v wget >/dev/null 2>&1; then
        wget -O "$TMP" "$URL"
    else
        echo "Error: neither curl nor wget is available. Install one, or download manually:" >&2
        echo "  $URL" >&2
        echo "  -> $MODEL_PATH" >&2
        rm -f "$TMP"
        exit 1
    fi
    if [ "$(file_size "$TMP")" -lt "$MIN_BYTES" ]; then
        echo "Error: downloaded file is smaller than expected — the release asset may have moved." >&2
        echo "Check: $URL" >&2
        rm -f "$TMP"
        exit 1
    fi
    mv "$TMP" "$MODEL_PATH"
    echo "✓ Saved to $MODEL_PATH"
fi

# Register the path in .env (read automatically via python-dotenv at runtime).
touch "$ENV_FILE"
if grep -q "^${ENV_VAR}=" "$ENV_FILE" 2>/dev/null; then
    # Replace the existing line in place (portable sed).
    tmp_env="$(mktemp)"
    grep -v "^${ENV_VAR}=" "$ENV_FILE" >"$tmp_env" || true
    mv "$tmp_env" "$ENV_FILE"
fi
echo "${ENV_VAR}=\"${REL_PATH}\"" >>"$ENV_FILE"
echo "✓ Set ${ENV_VAR} in $ENV_FILE"

echo
echo "Done. Install the ML stack (if you haven't) and run with the MACE optimizer:"
echo "  uv sync --extra mace"
echo '  uv run --extra mace oin-smiles oin2xyz "<OIN string>" --optimizer mace-omol-0-extra-large-1024'
