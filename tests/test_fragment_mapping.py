"""Unit tests for fragment mapping in Direct Parser (MiniPRD_DirectParser_FragmentMapping_v0.2.2)."""

import unittest
from src.oinsmiles.generation.oin_parser import _extract_oin_constraints


class TestFragmentMapping(unittest.TestCase):
    """Test suite for fragment-to-atom mapping returned by _extract_oin_constraints."""

    def test_atom_ordering_deterministic(self):
        """Test: Fragment mapping is byte-identical across multiple calls (determinism)."""
        cisplatin_oin = "[Pt_SPL].[Cl]{0}.[Cl]{1}"

        # Call 10 times and collect fragment mappings
        mappings = []
        for _ in range(10):
            _, _, frag_to_atom = _extract_oin_constraints(cisplatin_oin)
            mappings.append(frag_to_atom)

        # All mappings should be identical
        for mapping in mappings[1:]:
            self.assertEqual(mapping, mappings[0], "Fragment mapping changed across calls (non-deterministic)")

    def test_cisplatin_mapping_correct(self):
        """Test: Cisplatin OIN has correct fragment-to-atom mapping."""
        cisplatin_oin = "[Pt_SPL].[Cl]{0}.[Cl]{1}"
        _, constraints, frag_to_atom = _extract_oin_constraints(cisplatin_oin)

        # Should have 3 fragments: Pt, Cl, Cl
        self.assertEqual(len(frag_to_atom), 3)

        # Fragment 0 (metal) should map to atom 0
        self.assertIn(0, frag_to_atom)
        self.assertEqual(frag_to_atom[0], [0])

        # Fragment 1 (first Cl) should map to exactly one atom
        self.assertIn(1, frag_to_atom)
        self.assertEqual(len(frag_to_atom[1]), 1)

        # Fragment 2 (second Cl) should map to exactly one atom
        self.assertIn(2, frag_to_atom)
        self.assertEqual(len(frag_to_atom[2]), 1)

        # All atom indices should be unique and cover 0-2
        all_atom_indices = []
        for atoms in frag_to_atom.values():
            all_atom_indices.extend(atoms)
        self.assertEqual(sorted(all_atom_indices), [0, 1, 2])

    def test_polydentate_mapping_ethanediamine(self):
        """Test: Bidentate ethylenediamine ligand maps to multiple atoms."""
        # cis-PtCl2(en): Pt, 2x Cl, and ethylenediamine (C2N2H8)
        cis_pt_en_oin = "[Pt_SPL].[Cl]{0}.[Cl]{1}.NCC{2}N{3}"

        _, _, frag_to_atom = _extract_oin_constraints(cis_pt_en_oin)

        # Should have 4 fragments
        self.assertEqual(len(frag_to_atom), 4)

        # Fragment 0 (Pt) has 1 atom
        self.assertEqual(len(frag_to_atom[0]), 1)

        # Fragment 1 (Cl) has 1 atom
        self.assertEqual(len(frag_to_atom[1]), 1)

        # Fragment 2 (Cl) has 1 atom
        self.assertEqual(len(frag_to_atom[2]), 1)

        # Fragment 3 (en: NCC{2}N{3}) should have ≥4 atoms (N, C, C, N minimum from connectivity)
        # Actually, looking at "NCC{2}N{3}", this is N-C-C-N unbracketed, so 4 atoms
        self.assertGreaterEqual(len(frag_to_atom[3]), 4)

    def test_no_gaps_in_atom_indices(self):
        """Test: Fragment mapping covers contiguous atom indices (0, 1, 2, ...)."""
        # Generic multi-fragment OIN
        oin_smiles = "[Pt_OCT].[NH3]{0}.[NH3]{1}.[Cl]{2}.[Cl]{3}.[Cl]{4}.[Cl]{5}"
        _, _, frag_to_atom = _extract_oin_constraints(oin_smiles)

        # Collect all atom indices
        all_indices = []
        for atoms in frag_to_atom.values():
            all_indices.extend(atoms)

        # Should be contiguous from 0 to N-1
        all_indices_sorted = sorted(all_indices)
        expected = list(range(len(all_indices)))
        self.assertEqual(all_indices_sorted, expected, "Atom indices have gaps or duplicates")

    def test_metal_always_at_fragment_zero(self):
        """Test: Metal is always in fragment 0 per project invariant."""
        test_cases = [
            "[Pt_SPL].[Cl]{0}.[Cl]{1}",
            "[Fe_TBP].[CO]{0}.[CO]{1}.[CO]{2}",
            "[Ir_OCT].[Cl]{0}.[NH3]{1}.[NH3]{2}.[NH3]{3}.[NH3]{4}.[NH3]{5}",
        ]

        for oin_smiles in test_cases:
            _, _, frag_to_atom = _extract_oin_constraints(oin_smiles)
            # Fragment 0 should exist
            self.assertIn(0, frag_to_atom, f"Fragment 0 missing in {oin_smiles}")
            # Fragment 0 should have at least one atom (the metal)
            self.assertGreater(len(frag_to_atom[0]), 0, f"Fragment 0 has no atoms in {oin_smiles}")


if __name__ == "__main__":
    unittest.main()
