"""Guard tests for the P8 BINAP / ``ff.clean_geometry`` speed regression fix.

P4 added ``Pd`` to ``ENABLED_METALS``, which pinned PdCl₂-BINAP's σ-P donors at the
short curated target (Pd–P 2.30 Å vs the 2.46 Å covalent sum). The short, hard
``ff_clean`` position constraint over-strained the UFF clean, quadrupling its failure
count (3 → 12) and pushing generation 34 s → 57 s (P8, seed 42, deterministic).

The fix routes **only** the Pd–P σ pair back to the covalent sum via
``SHORT_PIN_EXEMPT_PAIRS`` / ``sigma_table_applies``, leaving every other enabled
pair — including Pd–Cl, Pd–N, Pd–O — on the curated table so their P4 fidelity gain
is preserved. These tests lock the seam (``_binding_distance``) so a future edit to
the exempt set is deliberate, not accidental.
"""

import unittest

from oinsmiles.generator3d.bond_lengths import (
    ENABLED_METALS,
    SHORT_PIN_EXEMPT_PAIRS,
    bond_length,
    sigma_table_applies,
)
from oinsmiles.generator3d.clean_geometry import _binding_distance

# Representative covalent radii (chem.py Cordero table). Exact sums are asserted via
# the covalent formula, never hard-coded, so a radius-table change cannot silently
# break these.
_PD_R, _P_R, _CL_R, _N_R = 1.39, 1.07, 1.02, 0.71


class TestShortPinExemptPolicy(unittest.TestCase):
    def test_pd_p_is_the_exempt_pair(self):
        self.assertIn(("Pd", "P"), SHORT_PIN_EXEMPT_PAIRS)

    def test_pd_cl_is_not_exempt(self):
        # Pd–Cl is not the speed culprit (with Pd–P exempt, ff_clean fail stays 3),
        # so it keeps its P4 fidelity gain.
        self.assertNotIn(("Pd", "Cl"), SHORT_PIN_EXEMPT_PAIRS)

    def test_pd_stays_enabled(self):
        # The fix narrows a single pair, not the whole metal: Pd–N/O/Cl still help.
        self.assertIn("Pd", ENABLED_METALS)

    def test_sigma_table_applies_truth_table(self):
        self.assertFalse(sigma_table_applies("Pd", "P"))  # exempt
        self.assertTrue(sigma_table_applies("Pd", "Cl"))  # kept
        self.assertTrue(sigma_table_applies("Pd", "N"))  # kept
        self.assertTrue(sigma_table_applies("Pt", "P"))  # only Pd–P is exempt
        self.assertFalse(sigma_table_applies("Fe", "C"))  # Fe not enabled


class TestBindingDistanceExemptPair(unittest.TestCase):
    def test_pd_p_sigma_falls_back_to_covalent(self):
        # The regression fix: a σ Pd–P donor must take the covalent sum, NOT the
        # 2.30 Å table value that over-strained ff_clean.
        d = _binding_distance("Pd", "P", _PD_R, _P_R, 1.0, is_sigma=True)
        self.assertAlmostEqual(d, (_PD_R + _P_R) * 1.0)
        self.assertNotAlmostEqual(d, 2.30)  # never the table value
        self.assertEqual(bond_length("Pd", "P"), 2.30)  # raw table still has it

    def test_pd_cl_sigma_still_uses_table(self):
        # Kept on the table -> its P4 fidelity gain survives.
        d = _binding_distance("Pd", "Cl", _PD_R, _CL_R, 1.0, is_sigma=True)
        self.assertAlmostEqual(d, 2.30)
        self.assertNotAlmostEqual(d, (_PD_R + _CL_R) * 1.0)

    def test_pd_n_sigma_still_uses_table(self):
        d = _binding_distance("Pd", "N", _PD_R, _N_R, 1.0, is_sigma=True)
        self.assertAlmostEqual(d, 2.10)

    def test_other_enabled_metal_unaffected(self):
        # Only Pd–P is exempt; Pt–Cl (an enabled pair) still takes the table.
        d = _binding_distance("Pt", "Cl", 1.36, _CL_R, 1.0, is_sigma=True)
        self.assertAlmostEqual(d, 2.31)

    def test_pd_p_haptic_always_covalent(self):
        # Haptic (is_sigma=False) was already covalent; exempting Pd–P must not
        # change that (it stays covalent, not the table value).
        d = _binding_distance("Pd", "P", _PD_R, _P_R, 1.0, is_sigma=False)
        self.assertAlmostEqual(d, (_PD_R + _P_R) * 1.0)

    def test_scale_multiplies_on_exempt_path(self):
        # The swept `scale` (pool-diversity source) must keep multiplying even on the
        # exempt covalent path, or the conformer pool collapses.
        d1 = _binding_distance("Pd", "P", _PD_R, _P_R, 1.0, is_sigma=True)
        d2 = _binding_distance("Pd", "P", _PD_R, _P_R, 1.2, is_sigma=True)
        self.assertAlmostEqual(d2 / d1, 1.2)

    def test_scale_multiplies_on_kept_table_path(self):
        d1 = _binding_distance("Pd", "Cl", _PD_R, _CL_R, 1.0, is_sigma=True)
        d2 = _binding_distance("Pd", "Cl", _PD_R, _CL_R, 1.2, is_sigma=True)
        self.assertAlmostEqual(d2 / d1, 1.2)


if __name__ == "__main__":
    unittest.main()
