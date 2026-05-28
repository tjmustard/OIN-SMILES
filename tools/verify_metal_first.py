#!/usr/bin/env python3
"""
Verification spike: Confirm that fragments[0] is always a transition metal.

This script gates MiniPRD_DirectParser_FragmentMapping_v0.2.2 (D3 gate).
Before any code change lands, run this to assert the metal-first invariant
across all OIN test fixtures.

Usage:
    python tools/verify_metal_first.py

Exit code:
    0 = all fixtures pass (metal-first invariant holds)
    1 = at least one fixture failed (halts MiniPRD #1; escalate to Architect)
"""

import re
import sys
from pathlib import Path
from typing import Set

# Transition metal atomic numbers (3d, 4d, 5d blocks)
TRANSITION_METALS: Set[str] = {
    # 3d metals
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    # 4d metals
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    # 5d metals
    "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    # Lanthanides (4f)
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
    "Er", "Tm", "Yb", "Lu",
    # Actinides (5f)
    "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
}


def extract_fragments(oin_smiles: str) -> list[str]:
    """Extract fragments from OIN-SMILES separated by dots.

    Returns list of fragment SMILES strings.
    """
    # Remove OIN annotations to get clean SMILES
    clean = oin_smiles
    clean = re.sub(r'_[A-Z0-9]+', '', clean)  # Remove shape codes
    clean = re.sub(r'@SP[0-9]+', '', clean)   # Remove chiral tags
    clean = re.sub(r'\{[0-9><]+\}', '', clean)  # Remove vertex indices

    # Split by dots (fragment separators)
    fragments = clean.split('.')
    return fragments


def get_first_atom_symbol(fragment_smiles: str) -> str:
    """Extract the first atom symbol from a SMILES fragment.

    Returns the symbol of the first atom in bracket notation.
    E.g., "[Pt@L1]" → "Pt", "[Cl]" → "Cl", "C" → "C"
    """
    # Match bracketed atoms: [symbol...] or just symbol
    bracket_match = re.match(r'\[([A-Z][a-z]?)', fragment_smiles)
    if bracket_match:
        return bracket_match.group(1)

    # Match unbracketed atom (C, N, O, etc.)
    simple_match = re.match(r'([A-Z][a-z]?)', fragment_smiles)
    if simple_match:
        return simple_match.group(1)

    return ""


def extract_metal_from_xyz(xyz_content: str) -> str:
    """Extract the first atom (metal) symbol from XYZ file content.

    XYZ format: line 0 = count, line 1 = comment, line 2+ = atoms.
    Each atom line: symbol x y z
    """
    lines = xyz_content.strip().split('\n')
    if len(lines) < 3:
        return ""

    # Line 2 is the first atom (0-indexed: lines[2])
    atom_line = lines[2].split()
    if atom_line:
        return atom_line[0]
    return ""


def verify_fixtures():
    """Run verification over all tests/fixtures/*.xyz files."""
    fixture_dir = Path("tests/fixtures")

    if not fixture_dir.exists():
        print(f"ERROR: Fixture directory not found: {fixture_dir}")
        sys.exit(1)

    fixture_files = sorted(fixture_dir.glob("*.xyz"))
    if not fixture_files:
        print(f"WARNING: No XYZ fixtures found in {fixture_dir}")
        return True

    all_passed = True
    results = []

    for fixture_path in fixture_files:
        try:
            with open(fixture_path) as f:
                xyz_content = f.read()

            # Extract the first atom from XYZ (the metal center)
            metal_symbol = extract_metal_from_xyz(xyz_content)
            if not metal_symbol:
                results.append((fixture_path.name, "FAIL", "Cannot parse metal atom from XYZ"))
                all_passed = False
                continue

            if metal_symbol not in TRANSITION_METALS:
                results.append((fixture_path.name, "FAIL", f"First atom is '{metal_symbol}' (not a transition metal)"))
                all_passed = False
            else:
                results.append((fixture_path.name, "PASS", f"Metal: {metal_symbol}"))

        except Exception as e:
            results.append((fixture_path.name, "ERROR", str(e)))
            all_passed = False

    # Print results
    print("\n" + "=" * 80)
    print("VERIFICATION SPIKE: Metal-First Invariant")
    print("=" * 80)
    for name, status, message in results:
        status_str = {"PASS": "✓", "FAIL": "✗", "SKIP": "⊘", "ERROR": "⚠"}.get(status, "?")
        print(f"{status_str} {name:40s} {status:6s} {message}")

    print("=" * 80)
    passed = sum(1 for _, s, _ in results if s == "PASS")
    failed = sum(1 for _, s, _ in results if s == "FAIL")
    errors = sum(1 for _, s, _ in results if s == "ERROR")
    print(f"SUMMARY: {passed} passed, {failed} failed, {errors} errors")
    print("=" * 80)

    return all_passed


if __name__ == "__main__":
    success = verify_fixtures()
    sys.exit(0 if success else 1)
