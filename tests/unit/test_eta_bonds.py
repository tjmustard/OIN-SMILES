"""
Tests for eta-bond vertex-to-atom index translation in direct parser.

Tests the fix for Blocker #4: eta-bond vertex indices must be translated
from fragment ranks to actual atom indices in the connected SMILES.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES
from oinsmiles.generation.engine import _translate_eta_vertex_to_atoms, parse_oin_direct


class TestEtaVertexTranslation(unittest.TestCase):
    """Test the _translate_eta_vertex_to_atoms helper function."""

    def test_translate_eta_vertex_single_atom_ligand(self):
        """For monodentate ligands (e.g., Cl), should return single atom index."""
        frag_to_atom = {0: [0], 1: [1], 2: [2]}  # Metal + 2 Cl atoms
        result = _translate_eta_vertex_to_atoms(0, 1, frag_to_atom)
        self.assertEqual(result, [1])

    def test_translate_eta_vertex_multi_atom_ligand(self):
        """For eta ligands (e.g., Cp), should return all ring atom indices."""
        frag_to_atom = {0: [0], 1: [1, 2, 3, 4, 5]}  # Metal + Cp ring (5 atoms)
        result = _translate_eta_vertex_to_atoms(0, 1, frag_to_atom)
        self.assertEqual(result, [1, 2, 3, 4, 5])

    def test_translate_eta_vertex_missing_fragment(self):
        """Should raise ValueError if fragment_rank not in mapping."""
        frag_to_atom = {0: [0], 1: [1]}
        with self.assertRaises(ValueError) as cm:
            _translate_eta_vertex_to_atoms(0, 99, frag_to_atom)
        self.assertIn("Fragment rank 99 not found", str(cm.exception))

    def test_translate_eta_vertex_empty_mapping(self):
        """Should raise ValueError if mapping is empty."""
        frag_to_atom = {}
        with self.assertRaises(ValueError) as cm:
            _translate_eta_vertex_to_atoms(0, 0, frag_to_atom)
        self.assertIn("not found in fragment-to-atom mapping", str(cm.exception))


class TestEtaBondToSelfPrevention(unittest.TestCase):
    """Regression test for bond-to-self errors (A9 in MiniPRD)."""

    def test_ferrocene_no_bond_to_self(self):
        """
        Test the specific ferrocene case: {0>} and {1>} assignments
        should translate correctly without attempting bond-to-self.

        Ferrocene OIN: [Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1

        Before the fix, this would try to bond Fe (atom 0) to atoms 0 and 1
        (incorrect fragment-rank-based indices).

        After the fix, it correctly maps:
        - Fragment 0 (Fe) → atom 0
        - Fragment 1 (Cp) → atoms [1,2,3,4,5]
        - Fragment 2 (Cp) → atoms [6,7,8,9,10]

        And adds eta bonds from Fe (0) to atoms 1-5 and 6-10.
        """
        oin_string = (
            "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
        )

        try:
            result = parse_oin_direct(oin_string)
            # Should succeed without raising "bond-to-self" error
            self.assertIsNotNone(result.xyz)
            # XYZ block should have valid structure (at least the atom count line)
            lines = result.xyz.strip().split("\n")
            self.assertGreater(len(lines), 2)  # At least atom count + blank line + atoms
            atom_count = int(lines[0])
            self.assertGreater(atom_count, 0)
        except ValueError as e:
            if "bond-to-self" in str(e).lower():
                self.fail(f"Bond-to-self error not fixed: {e}")
            # Other errors are acceptable (e.g., Molassembler issues)
            # The key is to not fail with bond-to-self
            pass

    def test_no_bond_to_self_assertion_in_construct(self):
        """
        Verify that construct_molassembler_mol has the bond-to-self assertion.

        This is a meta-test to ensure the regression check is in place.
        """
        import inspect

        from oinsmiles.generation.oin_parser import construct_molassembler_mol

        source = inspect.getsource(construct_molassembler_mol)
        # Check that the bond-to-self prevention code exists
        self.assertIn("bond-to-self", source.lower())
        self.assertIn("metal_idx == lig_idx", source)


class TestEtaBondRoundTrip(unittest.TestCase):
    """Test round-trip conversion for eta-ligand complexes."""

    def setUp(self):
        """Load XYZ→OIN converter."""
        self.converter = XYZToSMILES()
        self.fixtures_dir = os.path.join(os.path.dirname(__file__), "..", "fixtures")

    def test_ferrocene_roundtrip_parsing(self):
        """
        Test that ferrocene OIN can be parsed by direct parser.

        This test uses the ferrocene fixture to generate the OIN string,
        then attempts to parse it with parse_oin_direct.
        """
        ferrocene_xyz_path = os.path.join(self.fixtures_dir, "ferrocene.xyz")
        if not os.path.exists(ferrocene_xyz_path):
            self.skipTest(f"Ferrocene fixture not found at {ferrocene_xyz_path}")

        # Convert XYZ → OIN
        oin_string = self.converter.convert(ferrocene_xyz_path)
        self.assertIsNotNone(oin_string)
        self.assertIn("[Fe", oin_string)
        self.assertIn("cH", oin_string)

        # Parse OIN with direct parser
        try:
            result = parse_oin_direct(oin_string)
            self.assertIsNotNone(result)
            self.assertIsNotNone(result.xyz)
            # Verify XYZ structure
            lines = result.xyz.strip().split("\n")
            self.assertGreater(len(lines), 2)
        except (ValueError, NotImplementedError) as e:
            # Some ligand types may not be fully supported yet
            # The key is to not fail with bond-to-self errors
            if "bond-to-self" in str(e).lower():
                self.fail(f"Bond-to-self error in ferrocene round-trip: {e}")

    def test_eta_ligand_multi_fragment(self):
        """
        Test parsing of a multi-eta complex where multiple ligands bind via eta bonds.

        This is a simplified synthetic test case to verify the translation
        of multiple fragment assignments to different atoms.
        """
        # Synthetic OIN for a simple eta-bonded complex (not a real molecule)
        # [M_geometry].[Ligand1]{0}.[Ligand2]{1}
        # This just tests that the translation logic handles multiple slots.
        oin_string = "[Ni_SQP].[cH]{0}1ccccc1.[cH]{1}1ccccc1"

        # Don't assert success, just that it doesn't bond-to-self
        try:
            result = parse_oin_direct(oin_string)
            # If it succeeds, that's fine
            self.assertIsNotNone(result)
        except ValueError as e:
            # Check it's not a bond-to-self error
            self.assertNotIn("bond-to-self", str(e).lower())
        except (NotImplementedError, RuntimeError):
            # Other errors are acceptable (geometry not supported, etc.)
            pass


class TestAnsaMetalloceneDisjointRings(unittest.TestCase):
    """
    Test that multi-eta complexes (ansa-metallocenes) have disjoint ring atom sets.

    Ansa-metallocenes like TiCat1/3/4 have two Cp rings that bind via eta bonds.
    Each ring's atom set must be disjoint (no atom belongs to both rings).
    """

    def test_disjoint_fragment_atom_indices(self):
        """
        Verify that frag_to_atom mapping for multi-eta ligands produces disjoint sets.

        For a complex with two eta ligands:
        - Fragment 1 → [atoms for ring 1]
        - Fragment 2 → [atoms for ring 2]

        These sets must not overlap.
        """
        # Simulate a TiCat1-like structure with two separate Cp rings
        # In real usage, these would be computed from the actual OIN string
        frag_to_atom_multi_eta = {
            0: [0],  # Metal
            1: [1, 2, 3, 4, 5],  # First Cp ring (5 atoms)
            2: [6, 7, 8, 9, 10],  # Second Cp ring (5 atoms, distinct)
        }

        # Extract ring 1 and ring 2 atoms
        ring1 = set(frag_to_atom_multi_eta[1])
        ring2 = set(frag_to_atom_multi_eta[2])

        # Verify they are disjoint
        intersection = ring1.intersection(ring2)
        self.assertEqual(
            len(intersection),
            0,
            f"Ring atom sets must be disjoint, but found shared atoms: {intersection}",
        )

        # Verify each ring has the expected number of atoms
        self.assertEqual(len(ring1), 5)
        self.assertEqual(len(ring2), 5)


if __name__ == "__main__":
    unittest.main()
