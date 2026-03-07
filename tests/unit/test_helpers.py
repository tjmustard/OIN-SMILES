"""Test-only helpers for OIN-SMILES chiral test suite.

DO NOT import these from src/oinsmiles — they are test utilities only.
"""
import re


def extract_ligand_smiles(oin_string: str) -> str:
    """Strip the metal fragment and slot markers from an OIN-SMILES string.

    Returns the dot-separated ligand SMILES suitable for RDKit CIP oracle
    testing.  The metal fragment is identified by its ``_GEO`` suffix
    (e.g. ``[Pt@SP1_SPL]``, ``[Fe_LIN]``).  Slot markers ``{N}``,
    ``{N>}``, ``{N<}`` are removed.

    Parameters
    ----------
    oin_string:
        Full OIN-SMILES string (V3.0 inline format).

    Returns
    -------
    str
        Ligand-only SMILES with slot markers stripped.

    Examples
    --------
    >>> extract_ligand_smiles("[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}")
    '[Cl].[Cl].N.N'
    """
    # Remove slot markers: {0}, {0>}, {1<}, etc.
    clean = re.sub(r'\{\d+[><]?\}', '', oin_string)
    # Split on fragment separator and drop the metal fragment
    frags = clean.split('.')
    return '.'.join(f for f in frags if not re.search(r'_[A-Z]{3}', f))
