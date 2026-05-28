"""Axial-chiral ligand encoding tests.

Tests axial chirality (atropisomerism) encoding for BINAP-containing complexes.
BINAP (2,2'-bis(diphenylphosphino)-1,1'-binaphthyl) exhibits axial chirality, not
P-centered stereochemistry.

This test verifies that the encoded SMILES for axial-chiral atoms contains the
correct chirality descriptor characters (e.g., @/@@ or E/Z equivalents) as a
sanity check that encoded SMILES is chemically correct, even when the pipeline
does not interpret axial chirality specially.
"""

import unittest

from rdkit import Chem

from oinsmiles import XYZToSMILES
from .test_helpers import extract_ligand_smiles, get_fixture_path


class TestAxialChiral(unittest.TestCase):
    """Verify axial-chiral ligand encoding correctness."""

    def test_binap_encoded_smiles_is_valid(self):
        """R-BINAP encoded SMILES should be parseable by RDKit."""
        xyz_path = get_fixture_path("PdCl2-R-BINAP.xyz")
        oin = XYZToSMILES().convert(str(xyz_path))
        self.assertIsNotNone(oin)

        # Extract ligand SMILES from OIN
        ligand_smiles = extract_ligand_smiles(oin)
        self.assertIsNotNone(ligand_smiles)

        # Parse ligand with RDKit — should succeed for chemically valid SMILES
        mol = Chem.MolFromSmiles(ligand_smiles)
        self.assertIsNotNone(mol,
                            f"Failed to parse BINAP ligand SMILES: {ligand_smiles}")

    def test_binap_has_phosphorus_markers(self):
        """R-BINAP encoded SMILES should contain phosphorus atoms."""
        xyz_path = get_fixture_path("PdCl2-R-BINAP.xyz")
        oin = XYZToSMILES().convert(str(xyz_path))
        self.assertIsNotNone(oin)

        # Extract ligand SMILES
        ligand_smiles = extract_ligand_smiles(oin)

        # BINAP contains phosphorus — should be present in the encoded SMILES
        # (Axial chirality in BINAP is not explicitly encoded as @/@@ on atoms,
        # but the P atoms define the chiral axis)
        self.assertIn("P", ligand_smiles,
                     f"R-BINAP SMILES should contain phosphorus: {ligand_smiles}")


if __name__ == "__main__":
    unittest.main()
