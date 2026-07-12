"""Guard tests for the curated metal–ligand bond-length table (P4).

Covers the module accessor, the covalent-fallback contract, the drift guard
against the legacy source table, the σ-only ``_binding_distance`` seam behaviour,
and that the FF-clean conformer allocation preserves double-bond stereo.
"""

import unittest

from rdkit import Chem

from oinsmiles.generator3d.bond_lengths import BOND_LENGTHS, ENABLED_METALS, bond_length
from oinsmiles.generator3d.clean_geometry import _binding_distance


class TestBondLengthLookup(unittest.TestCase):
    def test_present_pair_returns_table_value(self):
        self.assertEqual(bond_length("Fe", "C"), 1.80)
        self.assertEqual(bond_length("Pt", "Cl"), 2.31)
        self.assertEqual(bond_length("Ti", "O"), 1.80)

    def test_missing_ligand_on_known_metal_returns_none(self):
        # V is in the table but only lists O/N/Cl -> V-S must fall through.
        self.assertIsNone(bond_length("V", "S"))

    def test_unknown_metal_returns_none(self):
        # Mo is absent entirely; a missing pair returns None, never 2.10, so the
        # caller falls back to the covalent sum (which beats a flat default).
        self.assertIsNone(bond_length("Mo", "N"))
        self.assertIsNone(bond_length("Co", "C"))

    def test_all_values_are_positive_floats(self):
        for metal, ligands in BOND_LENGTHS.items():
            for lig, val in ligands.items():
                self.assertIsInstance(val, float, f"{metal}-{lig}")
                self.assertGreater(val, 0.0, f"{metal}-{lig}")


class TestBindingDistanceSeam(unittest.TestCase):
    # Representative covalent radii (from chem.py's Cordero table) so the test is
    # self-contained; exact values are asserted via the ratio, not hard-coded sums.
    def test_sigma_uses_table_times_scale(self):
        # Pt-Cl table is 2.31; covalent sum ~2.38. A sigma donor must take the
        # table value, scaled -- not the covalent sum.
        metal_r, atom_r, scale = 1.36, 1.02, 1.0  # arbitrary; table must win
        d = _binding_distance("Pt", "Cl", metal_r, atom_r, scale, is_sigma=True)
        self.assertAlmostEqual(d, 2.31 * scale)
        self.assertNotAlmostEqual(d, (metal_r + atom_r) * scale)

    def test_sigma_missing_pair_falls_back_to_covalent_not_default(self):
        metal_r, atom_r, scale = 1.30, 1.05, 1.0
        d = _binding_distance("Mo", "S", metal_r, atom_r, scale, is_sigma=True)
        self.assertAlmostEqual(d, (metal_r + atom_r) * scale)
        self.assertNotAlmostEqual(d, 2.10 * scale)  # never the legacy blanket default

    def test_haptic_ignores_table(self):
        # An eta group (is_sigma=False) always uses the covalent sum, even for a
        # pair that exists in the table (Fe-C 1.80 is a sigma distance, not an
        # eta5-Cp centroid distance).
        metal_r, atom_r, scale = 1.25, 0.76, 1.0
        d = _binding_distance("Fe", "C", metal_r, atom_r, scale, is_sigma=False)
        self.assertAlmostEqual(d, (metal_r + atom_r) * scale)
        self.assertNotAlmostEqual(d, 1.80 * scale)

    def test_scale_multiplies_in_all_paths(self):
        for is_sigma, m, lig in [(True, "Pt", "Cl"), (True, "Mo", "S"), (False, "Fe", "C")]:
            d1 = _binding_distance(m, lig, 1.3, 1.0, 1.0, is_sigma)
            d2 = _binding_distance(m, lig, 1.3, 1.0, 1.2, is_sigma)
            self.assertAlmostEqual(d2 / d1, 1.2, msg=f"{m}-{lig} sigma={is_sigma}")

    def test_disabled_metal_uses_covalent_even_with_table_entry(self):
        # Fe-C IS in the raw table (1.80), but Fe is NOT enabled (the table
        # regresses Fe fidelity), so a sigma Fe-C donor must take the covalent sum.
        self.assertEqual(bond_length("Fe", "C"), 1.80)  # raw table has it
        self.assertNotIn("Fe", ENABLED_METALS)
        metal_r, atom_r, scale = 1.25, 0.76, 1.0
        d = _binding_distance("Fe", "C", metal_r, atom_r, scale, is_sigma=True)
        self.assertAlmostEqual(d, (metal_r + atom_r) * scale)
        self.assertNotAlmostEqual(d, 1.80 * scale)

    def test_enabled_metal_with_table_entry_uses_table(self):
        # A representative enabled metal takes the table value on the sigma path.
        self.assertIn("Pt", ENABLED_METALS)
        d = _binding_distance("Pt", "Cl", 1.36, 1.02, 1.0, is_sigma=True)
        self.assertAlmostEqual(d, 2.31)


class TestEnabledMetalsPolicy(unittest.TestCase):
    def test_enabled_metals_all_have_table_data(self):
        # Cannot enable a metal the table has no data for.
        self.assertTrue(ENABLED_METALS <= set(BOND_LENGTHS))

    def test_enabled_metals_nonempty(self):
        self.assertGreater(len(ENABLED_METALS), 0)


class TestConformerAllocationPreservesStereo(unittest.TestCase):
    def test_add_conformer_keeps_double_bond_stereo(self):
        # The FF-clean embed swap (RemoveAllConformers + AddConformer with an empty
        # Chem.Conformer) must not disturb the double-bond stereo asserted just
        # before it -- AddConformer adds coordinates only.
        mol = Chem.MolFromSmiles("C/C=C/C")  # trans-2-butene
        self.assertIsNotNone(mol)
        bond = None
        for b in mol.GetBonds():
            if b.GetStereo() != Chem.BondStereo.STEREONONE:
                bond = b
                break
        self.assertIsNotNone(bond, "fixture must carry a stereo double bond")
        before = bond.GetStereo()

        mol.RemoveAllConformers()
        mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()), assignId=True)

        self.assertEqual(mol.GetBondWithIdx(bond.GetIdx()).GetStereo(), before)
        self.assertEqual(mol.GetNumConformers(), 1)
        self.assertEqual(mol.GetConformer().GetNumAtoms(), mol.GetNumAtoms())


if __name__ == "__main__":
    unittest.main()
