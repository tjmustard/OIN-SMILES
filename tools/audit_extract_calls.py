#!/usr/bin/env python3
"""
Audit tool: Find all call sites of extract_oin_constraints.

This tool enumerates every call site of extract_oin_constraints() in src/
and tests/ to verify that the rename from extract_oin_constraints to
_extract_oin_constraints is complete.

Usage:
    python tools/audit_extract_calls.py

Output:
    - Per-line references to both prefixed (_extract_oin_constraints) and
      unprefixed (extract_oin_constraints) versions
    - Summary count of each
    - Warning if any unprefixed references remain after the rename (Task 3+)
"""

import re
import sys
from pathlib import Path
from typing import List, Tuple

def find_references(pattern: str) -> List[Tuple[Path, int, str]]:
    """Find all references to a function pattern in src/ and tests/.

    Returns list of (file_path, line_number, line_content) tuples.
    """
    results = []
    search_paths = [Path("src"), Path("tests")]

    for search_path in search_paths:
        if not search_path.exists():
            continue

        for py_file in search_path.rglob("*.py"):
            try:
                with open(py_file) as f:
                    lines = f.readlines()

                for i, line in enumerate(lines, 1):
                    if re.search(pattern, line):
                        results.append((py_file, i, line.rstrip()))

            except Exception as e:
                print(f"WARNING: Could not read {py_file}: {e}", file=sys.stderr)

    return results


def main():
    print("\n" + "=" * 100)
    print("AUDIT: Call Sites of extract_oin_constraints()")
    print("=" * 100 + "\n")

    # Find both prefixed (_extract_oin_constraints) and unprefixed versions
    unprefixed = find_references(r'\bextract_oin_constraints\b')
    prefixed = find_references(r'\b_extract_oin_constraints\b')

    # Print unprefixed (should be zero after rename)
    if unprefixed:
        print("UNPREFIXED (extract_oin_constraints) — SHOULD BE ZERO AFTER TASK 3:")
        for file_path, line_num, line_content in unprefixed:
            snippet = line_content[:80].replace('\n', '')
            print(f"  {file_path}:{line_num}: {snippet}")
        print()

    # Print prefixed (_extract_oin_constraints) — expected to exist
    if prefixed:
        print("PREFIXED (_extract_oin_constraints) — Expected call sites:")
        for file_path, line_num, line_content in prefixed:
            snippet = line_content[:80].replace('\n', '')
            print(f"  {file_path}:{line_num}: {snippet}")
        print()

    # Summary
    print("=" * 100)
    print(f"SUMMARY: {len(prefixed)} prefixed, {len(unprefixed)} unprefixed")
    if len(unprefixed) > 0:
        print("⚠️  WARNING: Unprefixed references remain; rename may be incomplete.")
        return 1
    else:
        print("✓ All references are properly prefixed (or zero total).")
        return 0
    print("=" * 100 + "\n")


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
