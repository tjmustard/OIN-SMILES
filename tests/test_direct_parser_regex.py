"""Unit tests for OIN regex preprocessor (_extract_oin_constraints) — OIN v3.6 format."""

import unittest
from src.oinsmiles.generation.oin_parser import _extract_oin_constraints


class TestExtractOINConstraintsV36(unittest.TestCase):
    """Test suite for _extract_oin_constraints function (OIN v3.6 format)."""

    def test_platinum_square_planar_basic(self):
        """Test: Platinum square planar with shape and vertex indices."""
        # Input: [Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}
        oin_smiles = "[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
        stripped, constraints, _ = _extract_oin_constraints(oin_smiles)

        # Verify constraints
        self.assertIn(0, constraints)
        self.assertEqual(constraints[0]['shape'], 'SPL')
        self.assertEqual(constraints[0]['chiral_tag'], '@SP1')
        # vertex_indices values are fragment ranks (metal = rank 0, ligands start at 1); list position = slot
        self.assertEqual(sorted(constraints[0]['vertex_indices']), [1, 2, 3, 4])

        # Verify stripped SMILES is clean
        self.assertNotIn('_', stripped)
        self.assertNotIn('{', stripped)
        self.assertNotIn('}', stripped)
        self.assertNotIn('@SP', stripped)

        # Verify atom maps are present
        self.assertIn(':', stripped)

    def test_iron_linear_with_heading_markers(self):
        """Test: Iron linear with Cp ligands and heading markers (winding direction)."""
        # Input: [Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1
        oin_smiles = "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
        stripped, constraints, _ = _extract_oin_constraints(oin_smiles)

        # Verify shape
        self.assertIn(0, constraints)
        self.assertEqual(constraints[0]['shape'], 'LIN')

        # Verify vertex indices (1 and 2, extracted from heading atoms {0>} and {1>})
        self.assertEqual(sorted(constraints[0]['vertex_indices']), [1, 2])

        # Verify stripped SMILES preserves aromatic notation (with atom maps added)
        self.assertIn('cH', stripped)  # aromatic CH present (may have atom maps)
        self.assertNotIn('_', stripped)
        self.assertNotIn('>', stripped)
        self.assertNotIn('<', stripped)

    def test_no_shape_code(self):
        """Test: OIN without shape code (only vertex indices)."""
        # Input: [Pd].[Cl]{0}.[Cl]{1}
        oin_smiles = "[Pd].[Cl]{0}.[Cl]{1}"
        stripped, constraints, _ = _extract_oin_constraints(oin_smiles)

        # Verify vertex indices are extracted even without shape
        self.assertIn(0, constraints)
        self.assertEqual(constraints[0]['vertex_indices'], [1, 2])

    def test_no_chiral_tag(self):
        """Test: OIN without chiral tag (shape only)."""
        # Input: [Pt_SPL].[Cl]{0}.[Cl]{1}
        oin_smiles = "[Pt_SPL].[Cl]{0}.[Cl]{1}"
        stripped, constraints, _ = _extract_oin_constraints(oin_smiles)

        # Verify shape is extracted
        self.assertEqual(constraints[0]['shape'], 'SPL')
        # Chiral tag should not be present
        self.assertNotIn('chiral_tag', constraints[0])
        # Vertex indices should be present
        self.assertEqual(constraints[0]['vertex_indices'], [1, 2])

    def test_various_shapes(self):
        """Test: Various polyhedral shape codes."""
        test_cases = [
            ("[Pt_SPL].[Cl]{0}", 'SPL'),
            ("[Pt_SQP].[Cl]{0}", 'SQP'),
            ("[Fe_OC].[Cl]{0}", 'OC'),
            ("[Ir_TPY].[Cl]{0}", 'TPY'),
            ("[Ti_TET].[Cl]{0}", 'TET'),
            ("[Mo_PBP].[Cl]{0}", 'PBP'),
        ]

        for oin_smiles, expected_shape in test_cases:
            with self.subTest(shape=expected_shape):
                stripped, constraints, _ = _extract_oin_constraints(oin_smiles)
                self.assertEqual(constraints[0]['shape'], expected_shape)

    def test_various_chiral_tags(self):
        """Test: Various stereochemistry codes."""
        test_cases = [
            "[Pt@SP1_SPL].[Cl]{0}",
            "[Pt@SP2_SPL].[Cl]{0}",
            "[Fe@SP1_OC].[Cl]{0}",
        ]

        for oin_smiles in test_cases:
            with self.subTest(oin=oin_smiles):
                stripped, constraints, _ = _extract_oin_constraints(oin_smiles)
                self.assertIn('chiral_tag', constraints[0])
                self.assertTrue(constraints[0]['chiral_tag'].startswith('@SP'))

    def test_heading_marker_clockwise(self):
        """Test: Heading marker for clockwise winding (>)."""
        # Input: [Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}1
        oin_smiles = "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}1"
        stripped, constraints, _ = _extract_oin_constraints(oin_smiles)

        # Heading markers should be stripped
        self.assertNotIn('>', stripped)
        # Vertex indices should still be extracted (1 and 2)
        self.assertEqual(sorted(constraints[0]['vertex_indices']), [1, 2])

    def test_heading_marker_counterclockwise(self):
        """Test: Heading marker for counterclockwise winding (<)."""
        # Input: [Fe_LIN].[cH]{0<}1[cH]{0}[cH]{0}1.[cH]{1<}1[cH]{1}[cH]{1}1
        oin_smiles = "[Fe_LIN].[cH]{0<}1[cH]{0}[cH]{0}1.[cH]{1<}1[cH]{1}[cH]{1}1"
        stripped, constraints, _ = _extract_oin_constraints(oin_smiles)

        # Heading markers should be stripped
        self.assertNotIn('<', stripped)
        # Vertex indices should still be extracted
        self.assertEqual(sorted(constraints[0]['vertex_indices']), [1, 2])

    def test_multiple_vertex_indices(self):
        """Test: Fragment with many vertex indices."""
        # Input: [Pt@SP1_SQP].[Cl]{0}.[Cl]{1}.[N]{2}.[N]{3}.[C]{4}.[C]{5}
        oin_smiles = "[Pt@SP1_SQP].[Cl]{0}.[Cl]{1}.[N]{2}.[N]{3}.[C]{4}.[C]{5}"
        stripped, constraints, _ = _extract_oin_constraints(oin_smiles)

        # Verify all vertex indices are captured
        expected_indices = [1, 2, 3, 4, 5, 6]
        self.assertEqual(sorted(constraints[0]['vertex_indices']), expected_indices)

    def test_atom_map_insertion(self):
        """Test: RDKit atom maps are inserted for tracking."""
        # Input: [Pt_SPL].[Cl]
        oin_smiles = "[Pt_SPL].[Cl]"
        stripped, constraints, _ = _extract_oin_constraints(oin_smiles)

        # Verify atom maps are inserted
        self.assertIn('[Pt:1]', stripped)
        self.assertIn('[Cl:2]', stripped)

    def test_complex_cyclic_smiles_preservation(self):
        """Test: Complex cyclic SMILES features are preserved."""
        # Input with aromatic rings in ligands
        oin_smiles = "[Fe_OC].[c-]{0}1[cH]{0}[cH]{0}[cH]{0}[c-]{0}1.[Cl]{1}"
        stripped, constraints, _ = _extract_oin_constraints(oin_smiles)

        # Verify aromatic SMILES notation is preserved (with atom maps added)
        self.assertIn('c-', stripped)  # aromatic carbon with negative charge
        self.assertIn('cH', stripped)  # aromatic CH
        # But vertex indices should be removed
        self.assertNotIn('{', stripped)
        self.assertNotIn('}', stripped)


class TestRoundTripCompatibilityV36(unittest.TestCase):
    """Test integration scenarios where preprocessor output feeds to next stages."""

    def test_output_is_valid_rdkit_input(self):
        """Test: Stripped SMILES with atom maps can be parsed by RDKit."""
        from rdkit import Chem

        test_cases = [
            "[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
            "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1",
            "[Pt_SQP].[Cl]{0}.[Cl]{1}",
        ]

        for oin_smiles in test_cases:
            with self.subTest(oin=oin_smiles):
                stripped, constraints, _ = _extract_oin_constraints(oin_smiles)
                # Should be parseable by RDKit (unsanitized)
                mol = Chem.MolFromSmiles(stripped, sanitize=False)
                self.assertIsNotNone(mol, f"Failed to parse: {stripped}")

    def test_constraints_dict_structure(self):
        """Test: Constraints dict has expected structure for downstream use."""
        oin_smiles = "[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
        stripped, constraints, _ = _extract_oin_constraints(oin_smiles)

        # Verify structure
        self.assertIsInstance(constraints, dict)
        for atom_idx, constraint_dict in constraints.items():
            self.assertIsInstance(atom_idx, int)
            self.assertIsInstance(constraint_dict, dict)
            # Expected keys (not all required)
            for key in constraint_dict:
                self.assertIn(key, ['shape', 'chiral_tag', 'vertex_indices'])

    def test_vertex_indices_are_integers(self):
        """Test: Extracted vertex indices are valid integers."""
        oin_smiles = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
        stripped, constraints, _ = _extract_oin_constraints(oin_smiles)

        if 'vertex_indices' in constraints[0]:
            for idx in constraints[0]['vertex_indices']:
                self.assertIsInstance(idx, int)
                self.assertGreaterEqual(idx, 0)


if __name__ == '__main__':
    unittest.main()
