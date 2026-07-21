#!/usr/bin/env bash
# Install CREST (https://github.com/crest-lab/crest) — the OPTIONAL external tool used
# by the conformer-invariance cross-check (tools/conformer_invariance_crest.py). CREST is
# NOT a dependency of OIN-SMILES; nothing else needs it.
#
# Strategy: prefer a conda-forge environment (micromamba/mamba/conda) because that also
# pulls CREST's tblite/dftd4 backends; fall back to the prebuilt static Linux binary from
# the GitHub releases. On success it prints the directory to add to PATH.
#
# Usage:  bash tools/install_crest.sh
set -euo pipefail

ENVNAME=crest
REPO="$(cd "$(dirname "$0")/.." && pwd)"

if command -v crest >/dev/null 2>&1; then
    echo "crest already on PATH: $(command -v crest)"
    crest --version 2>/dev/null | grep -i version || true
    exit 0
fi

# --- Path 1: a conda-forge environment ---------------------------------------------
for MGR in micromamba mamba conda; do
    if command -v "$MGR" >/dev/null 2>&1; then
        echo "Using $MGR to install crest from conda-forge into env '$ENVNAME'..."
        if "$MGR" env list 2>/dev/null | grep -qw "$ENVNAME"; then
            echo "  env '$ENVNAME' already exists; reusing it."
        else
            "$MGR" create -y -n "$ENVNAME" -c conda-forge crest
        fi
        BIN="$("$MGR" run -n "$ENVNAME" bash -c 'command -v crest')"
        echo "Installed: $BIN"
        "$MGR" run -n "$ENVNAME" crest --version 2>/dev/null | grep -i version || true
        echo
        echo "Add it to PATH for this shell:"
        echo "    export PATH=\"$(dirname "$BIN"):\$PATH\""
        echo "(tools/run_conformer_crest_sweep.sh discovers this env automatically.)"
        exit 0
    fi
done

# --- Path 2: prebuilt static binary from GitHub releases ---------------------------
echo "No conda package manager found; downloading a prebuilt binary..."
DEST="$REPO/.crest/bin"
mkdir -p "$DEST"
API="https://api.github.com/repos/crest-lab/crest/releases/latest"
URL="$(curl -fsSL "$API" \
    | grep browser_download_url \
    | grep -Ei 'ubuntu|linux' \
    | grep -Ei '\.tar\.(xz|gz)' \
    | head -1 | cut -d'"' -f4)"
if [ -z "${URL:-}" ]; then
    echo "ERROR: could not locate a Linux CREST release asset. Install manually from"
    echo "       https://github.com/crest-lab/crest/releases and put 'crest' on PATH." >&2
    exit 1
fi
echo "Downloading $URL"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
curl -fsSL "$URL" -o "$TMP/crest.tar"
tar -xf "$TMP/crest.tar" -C "$TMP"
CBIN="$(find "$TMP" -type f -name crest | head -1)"
[ -n "$CBIN" ] || { echo "ERROR: no 'crest' binary in the archive." >&2; exit 1; }
install -m 0755 "$CBIN" "$DEST/crest"
echo "Installed: $DEST/crest"
"$DEST/crest" --version 2>/dev/null | grep -i version || true
echo
echo "Add it to PATH for this shell:"
echo "    export PATH=\"$DEST:\$PATH\""
