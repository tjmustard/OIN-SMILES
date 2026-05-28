"""Test helpers for chirality validation and OIN parsing utilities."""

import re
from pathlib import Path
from rdkit import Chem


def extract_ligand_smiles(oin_string: str) -> str:
    """Extract ligand SMILES from OIN string by removing metal and slot markers.

    Test-only utility: strips the metal fragment and slot markers ({0}, {1}, etc.)
    from an OIN string to isolate the ligand SMILES for independent RDKit CIP
    oracle validation.

    Parameters
    ----------
    oin_string : str
        OIN string in v3.6 inline format, e.g.:
        "[Pt@SP1_SPL].[N@@H3]{0}.[Cl]{1}.[Cl]{2}.[N@@H3]{3}"

    Returns
    -------
    str
        Non-canonical SMILES containing only ligand fragments with @/@@ markers
        preserved, e.g.: "N@@H3.Cl.Cl.N@@H3"
    """
    # Remove metal fragment: [Pt@SP1_SPL] → ""
    metal_regex = re.compile(r"\[([A-Z][a-z]?)(?:@[A-Z0-9_]+)?\]")
    no_metal = metal_regex.sub("", oin_string)

    # Remove slot markers: {0}, {1>, etc. → ""
    slot_regex = re.compile(r"\{(\d+)[><]?\}")
    no_slots = slot_regex.sub("", no_metal)

    # Strip leading/trailing dots
    result = no_slots.strip(".")

    return result


def get_chiral_atom(mol: Chem.Mol, atomic_num: int) -> Chem.Atom | None:
    """Find the first P or N atom in a molecule (assumes one chiral center per ligand).

    Test helper: returns the first atom matching the given atomic number.
    For multi-chiral ligands, this will return only the first match.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit Mol object
    atomic_num : int
        Atomic number (7 for N, 15 for P)

    Returns
    -------
    Chem.Atom | None
        The first atom with the given atomic number, or None if not found.
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == atomic_num:
            return atom
    return None


def get_fixture_path(filename: str) -> Path:
    """Return full path to a fixture file in tests/fixtures/."""
    return Path(__file__).parent / "fixtures" / filename
