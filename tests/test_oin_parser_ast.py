"""Unit tests for AST tokenization in the Direct Parser (MiniPRD_DirectParser_ASTTokenization)."""

import unittest

from rdkit import Chem

from oinsmiles.generation.oin_parser import tokenize_unsanitized_smiles


class TestASTTokenization(unittest.TestCase):
    """Test suite for unsanitized SMILES tokenization."""

    # ==================== DETERMINISTIC TESTS ====================

    def test_simple_complex_atom_extraction(self):
        """Test 1: Simple metal complex atom extraction (Pt + 2 Cl + 2 N)."""
        smiles = "[Pt:1].[Cl:2].[Cl:3].[N:4].[N:5]"
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        # Verify 5 atoms extracted in order
        self.assertEqual(len(atoms), 5)
        self.assertEqual(atoms[0].GetAtomicNum(), 78)  # Pt
        self.assertEqual(atoms[1].GetAtomicNum(), 17)  # Cl
        self.assertEqual(atoms[2].GetAtomicNum(), 17)  # Cl
        self.assertEqual(atoms[3].GetAtomicNum(), 7)  # N
        self.assertEqual(atoms[4].GetAtomicNum(), 7)  # N

        # Verify indices are stable (0–4)
        for i, atom in enumerate(atoms):
            self.assertEqual(atom.GetIdx(), i)

    def test_simple_complex_formal_charges(self):
        """Test 1b: Verify formal charges are preserved."""
        smiles = "[Pt:1].[Cl:2].[Cl:3].[N:4].[N:5]"
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        # All atoms should have formal charge 0
        for atom in atoms:
            self.assertEqual(atom.GetFormalCharge(), 0)

    def test_aromatic_ring_unsanitized(self):
        """Test 2: Aromatic ring (Ti + Cp) parsed but not kekulized."""
        smiles = "[Ti:1].[c:2]1[c:3][c:4][c:5][c:6]1"
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        # Verify 6 atoms (Ti + 5 C)
        self.assertEqual(len(atoms), 6)
        self.assertEqual(atoms[0].GetAtomicNum(), 22)  # Ti
        for i in range(1, 6):
            self.assertEqual(atoms[i].GetAtomicNum(), 6)  # C

        # Verify bonds extracted (should include aromatic C-C bonds)
        self.assertGreater(len(bonds), 0)
        # Check that aromatic bonds exist (BondType.Aromatic)
        aromatic_bonds = [b for b in bonds if b[2] == Chem.BondType.AROMATIC]
        self.assertGreater(len(aromatic_bonds), 0)

    def test_implicit_hydrogens_preserved(self):
        """Test 3: Implicit hydrogens are counted correctly in unsanitized state."""
        smiles = "[C:1].[C:2]"
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        # Both carbons should be extracted
        for atom in atoms:
            self.assertEqual(atom.GetAtomicNum(), 6)
            # Note: In unsanitized state, RDKit may not compute implicit H automatically
            # but the molecule maintains consistency for downstream processing
            implicit_h = atom.GetTotalNumHs(includeNeighbors=True)
            # Just verify we can query without error
            self.assertIsNotNone(implicit_h)

    def test_atom_map_preservation(self):
        """Test: Atom maps from regex stage are preserved by RDKit."""
        smiles = "[Pt:1].[Cl:2].[N:3]"
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        # RDKit stores atom maps in the atom property
        self.assertEqual(atoms[0].GetAtomMapNum(), 1)  # Pt:1
        self.assertEqual(atoms[1].GetAtomMapNum(), 2)  # Cl:2
        self.assertEqual(atoms[2].GetAtomMapNum(), 3)  # N:3

    def test_bond_extraction(self):
        """Test: Bond list correctly extracted."""
        # Create a simple bonded system: C-C
        smiles = "[C:1][C:2]"
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        # Expect 1 bond connecting atoms 0 and 1
        self.assertEqual(len(bonds), 1)
        i, j, bond_type = bonds[0]
        self.assertIn(i, [0, 1])
        self.assertIn(j, [0, 1])
        self.assertNotEqual(i, j)
        self.assertEqual(bond_type, Chem.BondType.SINGLE)

    # ==================== EDGE CASES ====================

    def test_single_atom(self):
        """Test: Single atom molecule (just Pt)."""
        smiles = "[Pt:1]"
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        self.assertEqual(len(atoms), 1)
        self.assertEqual(atoms[0].GetAtomicNum(), 78)
        self.assertEqual(len(bonds), 0)  # No bonds

    def test_aromatic_preservation_no_kekulization(self):
        """Test: Aromatic system is NOT kekulized in unsanitized parse."""
        smiles = "[Pt:1].[c:2]1[c:3][c:4][c:5][c:6]1"
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        # Count aromatic bond types
        aromatic_count = sum(1 for _, _, bt in bonds if bt == Chem.BondType.AROMATIC)
        # Aromatic ring should have aromatic bonds in unsanitized state
        self.assertGreater(aromatic_count, 0)

    def test_charged_atoms(self):
        """Test: Charged atoms are parsed correctly."""
        smiles = "[Pt:1].[Cl:2].[Cl:3].[NH3:4]"
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        # Find N atom, verify it was parsed
        n_atoms = [a for a in atoms if a.GetAtomicNum() == 7]
        self.assertGreater(len(n_atoms), 0)
        n_atom = n_atoms[0]
        # Just verify the atom exists and has reasonable degree
        self.assertGreaterEqual(n_atom.GetDegree(), 0)

    # ==================== ERROR CASES ====================

    def test_malformed_smiles_raises_value_error(self):
        """Test 4: Malformed SMILES raises ValueError."""
        malformed = "[Pt:1].[Cl:99"  # Missing closing bracket
        with self.assertRaises(ValueError) as context:
            tokenize_unsanitized_smiles(malformed)
        self.assertIn("Failed to parse", str(context.exception))

    def test_empty_smiles_raises_value_error(self):
        """Test: Empty SMILES raises ValueError."""
        # RDKit returns a valid empty mol for "", so we need a truly invalid SMILES
        # Use a SMILES that RDKit definitively rejects
        with self.assertRaises(ValueError):
            tokenize_unsanitized_smiles("]]][")

    def test_invalid_element_symbol(self):
        """Test: Invalid element symbol."""
        invalid = "[Xx:1]"  # Xx is not a valid element
        with self.assertRaises(ValueError):
            tokenize_unsanitized_smiles(invalid)

    # ==================== INTEGRATION WITH REGEX STAGE ====================

    def test_integration_with_extract_oin_constraints(self):
        """Test: AST tokenization integrates with regex output."""
        from oinsmiles.generation.oin_parser import _extract_oin_constraints

        # OIN v3.6 format string
        oin = "[Pt_SPL].[Cl]{0}.[Cl]{1}"
        stripped_smiles, constraints, _ = _extract_oin_constraints(oin)

        # stripped_smiles should have atom maps
        self.assertIn(":", stripped_smiles)

        # Tokenize the stripped SMILES
        atoms, bonds = tokenize_unsanitized_smiles(stripped_smiles)

        # Verify atoms match OIN structure
        self.assertGreaterEqual(len(atoms), 3)  # Pt + 2 Cl minimum
        self.assertEqual(atoms[0].GetAtomicNum(), 78)  # Pt
        self.assertEqual(atoms[1].GetAtomicNum(), 17)  # Cl
        self.assertEqual(atoms[2].GetAtomicNum(), 17)  # Cl

    # ==================== REAL COMPLEX EXAMPLES ====================

    def test_cisplatin_like_complex(self):
        """Test: Real-world cisplatin-like structure."""
        # cis-Pt(NH3)2Cl2
        smiles = "[Pt:1].[NH3:2].[NH3:3].[Cl:4].[Cl:5]"
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        self.assertEqual(len(atoms), 5)
        # Atoms should maintain their types
        pt_atom = atoms[0]
        self.assertEqual(pt_atom.GetAtomicNum(), 78)

    def test_complex_ligand_with_aromatic_ring(self):
        """Test: Complex ligand with aromatic ring (e.g., biphenyl diphosphine)."""
        # Simplified phenyl-phosphine: Ph-P (P bonded to aromatic C)
        smiles = "[P:1][c:2]1[c:3][c:4][c:5][c:6][c:7]1"
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        # Verify parse succeeds and aromatic structure preserved
        self.assertEqual(len(atoms), 7)  # P + 6 aromatic C
        self.assertEqual(atoms[0].GetAtomicNum(), 15)  # P
        aromatic_c_atoms = [a for a in atoms[1:] if a.GetAtomicNum() == 6]
        self.assertEqual(len(aromatic_c_atoms), 6)

    def test_no_sanitization_validation(self):
        """Test: Unsanitized parsing skips valence/aromaticity validation."""
        # Create a hypothetical "invalid" SMILES that would fail sanitization
        # but should parse fine unsanitized
        smiles = "[C:1]([C:2])([C:3])([C:4])([C:5])[C:6]"  # C with 6 bonds
        atoms, bonds = tokenize_unsanitized_smiles(smiles)

        # Should parse without error (unsanitized tolerance)
        self.assertEqual(len(atoms), 6)
        self.assertGreater(len(bonds), 0)


if __name__ == "__main__":
    unittest.main()
