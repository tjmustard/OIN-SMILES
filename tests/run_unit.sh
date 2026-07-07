#!/bin/bash
# Run all unit tests
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: bash tests/run_unit.sh"
    echo ""
    echo "Runs all unit tests in the tests/unit directory."
    echo ""
    echo "Options:"
    echo "  -h, --help            Show this help message and exit"
    echo ""
    exit 0
fi

uv run python -m unittest discover tests/unit
