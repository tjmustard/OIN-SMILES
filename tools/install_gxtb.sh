#!/bin/bash

set -e

echo "Detecting OS and architecture..."
OS=$(uname -s)
ARCH=$(uname -m)

if [ "$OS" = "Linux" ] && [ "$ARCH" = "x86_64" ]; then
    URL="https://raw.githubusercontent.com/grimme-lab/g-xtb/main/binaries/xtb-6.7.1-gxtb-140526-linux-x86_64.tar.xz"
    FILE="xtb.tar.xz"
    EXTRACT_CMD="tar -xf"
elif [ "$OS" = "Darwin" ] && [ "$ARCH" = "x86_64" ]; then
    URL="https://raw.githubusercontent.com/grimme-lab/g-xtb/main/binaries/xtb-6.7.1-gxtb-140526-macos-x86_64.tar.gz"
    FILE="xtb.tar.gz"
    EXTRACT_CMD="tar -xzf"
elif [ "$OS" = "Darwin" ] && [ "$ARCH" = "arm64" ]; then
    URL="https://raw.githubusercontent.com/grimme-lab/g-xtb/main/binaries/xtb-6.7.1-gxtb-140526-macos-arm64.tar.gz"
    FILE="xtb.tar.gz"
    EXTRACT_CMD="tar -xzf"
else
    echo "Error: Unsupported OS or architecture ($OS $ARCH). Please download manually from https://github.com/grimme-lab/g-xtb/releases"
    exit 1
fi

echo "Downloading g-xTB binary for $OS ($ARCH)..."

# Create a temporary directory
TMP_DIR=$(mktemp -d)
ORIG_DIR=$PWD

cd "$TMP_DIR"

# Download the archive
curl -sL "$URL" -o "$FILE"

# Extract the archive
echo "Extracting..."
$EXTRACT_CMD "$FILE"

# Determine the target installation path
# If we are in a uv virtualenv, put it in .venv/bin
if [ -d "$ORIG_DIR/.venv/bin" ]; then
    TARGET_DIR="$ORIG_DIR/.venv/bin"
elif [ -n "$VIRTUAL_ENV" ]; then
    TARGET_DIR="$VIRTUAL_ENV/bin"
else
    TARGET_DIR="$HOME/.local/bin"
    mkdir -p "$TARGET_DIR"
fi

echo "Installing to $TARGET_DIR..."
# We will copy all binaries from bin/ to the target directory
cp -a xtb-*/bin/* "$TARGET_DIR/"
# We should also copy share/ if it exists because xtb sometimes needs parameter files (though gxtb says no external param file is needed).
if [ -d xtb-*/share ]; then
    SHARE_DIR="$(dirname "$TARGET_DIR")/share"
    mkdir -p "$SHARE_DIR"
    cp -a xtb-*/share/* "$SHARE_DIR/" 2>/dev/null || true
fi

# Clean up
cd "$ORIG_DIR"
rm -rf "$TMP_DIR"

echo "g-xTB successfully installed to $TARGET_DIR."
echo "You can verify it by running: xtb --version"
