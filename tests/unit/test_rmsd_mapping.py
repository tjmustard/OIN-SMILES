"""Guard tests for the coordination-sphere RMSD metric (tests/integration/rmsd_utils.py).

Every case here is a regression floor for a defect that made a chemically-correct
round-trip report as catastrophically bad geometry in the tmCAT/tmPHOTO baseline:

- 32 rows returned sentinel 999 ("no metal") because Y and Sc were missing from a
  hand-copied metal list. See test_metal_detection.
- 30 rows returned 996/997 because the input-side coordination sphere was chosen with a
  distance cutoff that both misses real long bonds and admits non-donors. See
  test_composition_matching.
"""

import os
import sys
import unittest

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D
from scipy.spatial.transform import Rotation

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../integration")))

import rmsd_utils
from rmsd_utils import (
    _extract_coordination_sphere,
    _find_metal,
    _select_matched_sphere,
    calculate_tmc_rmsd,
    calculate_tmc_rmsd_detailed,
)

SENTINEL_FLOOR = 900.0


def _complex(metal_sym, donors, bonded=True):
    """Build a metal at the origin with `donors` = [(symbol, (x, y, z)), ...].

    bonded=False mimics Chem.MolFromXYZFile: coordinates but no connectivity, which is
    exactly what the harness hands in as the *input* structure.
    """
    rw = Chem.RWMol()
    metal_idx = rw.AddAtom(Chem.Atom(metal_sym))
    conf = Chem.Conformer(1 + len(donors))
    conf.SetAtomPosition(metal_idx, Point3D(0.0, 0.0, 0.0))
    for sym, (x, y, z) in donors:
        idx = rw.AddAtom(Chem.Atom(sym))
        if bonded:
            rw.AddBond(metal_idx, idx, Chem.BondType.DATIVE)
        conf.SetAtomPosition(idx, Point3D(float(x), float(y), float(z)))
    mol = rw.GetMol()
    mol.AddConformer(conf, assignId=True)
    return mol


class TestMetalDetection(unittest.TestCase):
    """Y (Z=39) and Sc (Z=21) were absent from a duplicated metal list (TD-005)."""

    def test_yttrium_and_scandium_are_metals(self):
        # All 32 sentinel-999 rows in the baseline were Y or Sc complexes.
        for symbol in ("Y", "Sc"):
            mol = _complex(symbol, [("Cl", (2.4, 0, 0))])
            self.assertEqual(_find_metal(mol), 0, f"{symbol} must be recognised as a metal")

    def test_metal_list_is_imported_not_copied(self):
        from oinsmiles.core.constants import TRANSITION_METALS_NUM

        self.assertIs(rmsd_utils.TRANSITION_METALS_NUM, TRANSITION_METALS_NUM)

    def test_yttrium_complex_yields_a_real_rmsd(self):
        donors = [("Cl", (2.6, 0, 0)), ("Cl", (0, 2.6, 0)), ("O", (0, 0, 2.3))]
        rmsd, reason = calculate_tmc_rmsd_detailed(
            _complex("Y", donors, bonded=False),
            _complex("Y", donors),
            mol2_bonded=_complex("Y", donors),
        )
        self.assertIsNone(reason)
        self.assertAlmostEqual(rmsd, 0.0, places=6)

    def test_no_metal_reports_a_reason_not_a_number(self):
        organic = _complex("C", [("O", (1.2, 0, 0))])
        rmsd, reason = calculate_tmc_rmsd_detailed(organic, organic, mol2_bonded=organic)
        self.assertIsNone(rmsd)
        self.assertIn("no transition metal", reason)

    def test_metal_element_mismatch_is_caught(self):
        # Aligning a Pd sphere onto a Ni sphere is a mapping failure, not a bad geometry.
        pd = _complex("Pd", [("N", (2.1, 0, 0))], bonded=False)
        ni = _complex("Ni", [("N", (2.1, 0, 0))])
        rmsd, reason = calculate_tmc_rmsd_detailed(pd, ni, mol2_bonded=ni)
        self.assertIsNone(rmsd)
        self.assertIn("metal element differs", reason)


class TestCompositionMatching(unittest.TestCase):
    """The input sphere follows the generated sphere's bonded composition."""

    def test_long_apical_bond_is_not_dropped(self):
        # DAPZIF: real Pd-N bonds at 2.11 and 2.57 A; the old cutoff was 2.54 A, so the
        # apical donor was lost by 0.03 A and the metric returned 996.
        donors = [("N", (2.11, 0, 0)), ("N", (0, 0, 2.57))]
        gen = _complex("Pd", donors)
        rmsd, reason = calculate_tmc_rmsd_detailed(
            _complex("Pd", donors, bonded=False), gen, mol2_bonded=gen
        )
        self.assertIsNone(reason)
        self.assertAlmostEqual(rmsd, 0.0, places=6)

    def test_spurious_near_atom_is_not_admitted(self):
        # ROJXIY: 8 carbons sit inside the cutoff but only 7 are bonded donors. The
        # input must contribute exactly the 7 nearest, not all 8.
        gen_donors = [("C", (1.76, 0, 0)), ("C", (0, 1.91, 0))]
        input_donors = gen_donors + [("C", (0, 0, 2.19))]  # non-donor, inside any cutoff
        gen = _complex("Mn", gen_donors)
        rmsd, reason = calculate_tmc_rmsd_detailed(
            _complex("Mn", input_donors, bonded=False), gen, mol2_bonded=gen
        )
        self.assertIsNone(reason)
        self.assertAlmostEqual(rmsd, 0.0, places=6)

    def test_selection_picks_the_nearest_atoms_of_each_element(self):
        mol = _complex(
            "Mn", [("C", (1.8, 0, 0)), ("C", (0, 2.2, 0)), ("C", (0, 0, 2.0))], bonded=False
        )
        coords = np.array(
            [mol.GetConformer().GetAtomPosition(i) for i in range(mol.GetNumAtoms())], dtype=float
        )
        sphere, reason, _code = _select_matched_sphere(mol, coords, 0, {"Mn": 1, "C": 2})
        self.assertIsNone(reason)
        # The 2.2 A carbon is the odd one out; the 1.8 and 2.0 carbons are selected.
        picked = sorted(float(np.linalg.norm(p)) for p in sphere["C"])
        self.assertEqual([round(d, 2) for d in picked], [1.80, 2.00])

    def test_zero_regression_when_counts_already_matched(self):
        """A sphere the old cutoff already got right must be selected identically.

        This is why the fix cannot regress a currently-passing complex: when the cutoff
        admitted exactly k atoms of element el, the k *nearest* atoms of el are that
        same set.
        """
        donors = [("N", (2.05, 0, 0)), ("N", (0, 2.05, 0)), ("Cl", (0, 0, 2.30))]
        # A far-away carbon that no cutoff would ever have admitted.
        mol_in = _complex("Pd", donors + [("C", (0, 0, -7.0))], bonded=False)
        coords = np.array(
            [mol_in.GetConformer().GetAtomPosition(i) for i in range(mol_in.GetNumAtoms())],
            dtype=float,
        )
        legacy = _extract_coordination_sphere(mol_in, coords, 0, use_bonds=False)
        matched, reason, _code = _select_matched_sphere(
            mol_in, coords, 0, {el: len(v) for el, v in legacy.items()}
        )
        self.assertIsNone(reason)
        for element, positions in legacy.items():
            np.testing.assert_allclose(
                np.sort(np.asarray(positions), axis=0),
                np.sort(np.asarray(matched[element]), axis=0),
                err_msg=f"selection diverged from the legacy cutoff for {element}",
            )


class TestCeiling(unittest.TestCase):
    """The ceiling is what stops composition-matching from papering over a real defect."""

    def test_absent_donor_is_reported_not_substituted(self):
        # The input genuinely lacks a second N donor; the nearest is 4.6 A away. Matching
        # must refuse rather than reach out and grab it.
        gen = _complex("Pd", [("N", (2.11, 0, 0)), ("N", (0, 0, 2.05))])
        mol_in = _complex("Pd", [("N", (2.11, 0, 0)), ("N", (0, 0, 4.6))], bonded=False)
        rmsd, reason = calculate_tmc_rmsd_detailed(mol_in, gen, mol2_bonded=gen)
        self.assertIsNone(rmsd)
        self.assertIn("beyond the", reason)
        self.assertIn("ceiling", reason)

    def test_element_entirely_missing_from_input(self):
        gen = _complex("Pd", [("N", (2.05, 0, 0)), ("Cl", (0, 0, 2.3))])
        mol_in = _complex("Pd", [("N", (2.05, 0, 0))], bonded=False)
        rmsd, reason = calculate_tmc_rmsd_detailed(mol_in, gen, mol2_bonded=gen)
        self.assertIsNone(rmsd)
        self.assertIn("only 0 Cl", reason)

    def test_ceiling_admits_a_stretched_but_real_bond(self):
        # ABETIK's Zr-C(allyl) overshoots r_cov(Zr) + r_cov(C) by ~0.29 A.
        reach = rmsd_utils._rcov(40) + rmsd_utils._rcov(6) + 0.29
        gen = _complex("Zr", [("C", (reach, 0, 0))])
        mol_in = _complex("Zr", [("C", (reach, 0, 0))], bonded=False)
        rmsd, reason = calculate_tmc_rmsd_detailed(mol_in, gen, mol2_bonded=gen)
        self.assertIsNone(reason)
        self.assertAlmostEqual(rmsd, 0.0, places=6)


class TestRmsdKernel(unittest.TestCase):
    """Mean RMSD, never max-per-atom (project convention)."""

    def test_identical_structures_score_zero(self):
        donors = [("N", (2.0, 0, 0)), ("N", (0, 2.0, 0)), ("Cl", (0, 0, 2.3))]
        gen = _complex("Fe", donors)
        rmsd, reason = calculate_tmc_rmsd_detailed(
            _complex("Fe", donors, bonded=False), gen, mol2_bonded=gen
        )
        self.assertIsNone(reason)
        self.assertAlmostEqual(rmsd, 0.0, places=9)

    def test_score_is_the_mean_not_the_worst_atom(self):
        # One donor displaced by 0.3 A, over a 4-atom sphere (Fe + 3 Cl), and no rotation
        # can improve on the identity => mean RMSD = 0.3 / sqrt(4) = 0.15. A max-per-atom
        # metric would report 0.30 here.
        gen_donors = [("Cl", (2.4, 0, 0)), ("Cl", (0, 2.4, 0)), ("Cl", (0, 0, 2.4))]
        displaced = [("Cl", (2.4, 0, 0)), ("Cl", (0, 2.4, 0)), ("Cl", (0, 0, 2.7))]
        gen = _complex("Fe", gen_donors)
        rmsd, reason = calculate_tmc_rmsd_detailed(
            _complex("Fe", displaced, bonded=False), gen, mol2_bonded=gen
        )
        self.assertIsNone(reason)
        self.assertAlmostEqual(rmsd, 0.15, places=6)

    def test_robust_path_beats_greedy_on_a_rotated_ring(self):
        """>5 atoms in an element group route to the anchor/ICP search.

        Greedy assigns atoms in the un-rotated frame, so it mis-pairs a symmetric ring
        whenever the generated complex is rotated relative to the input. The robust
        search is floored by the greedy estimate and so can only ever lower the score.
        """
        ring = np.array(
            [
                [2.3 * np.cos(t), 2.3 * np.sin(t), 1.5]
                for t in np.linspace(0, 2 * np.pi, 6, endpoint=False)
            ]
        )
        rotated = Rotation.from_euler("y", 110, degrees=True).apply(ring)
        metal = np.array([[0.0, 0.0, 0.0]])

        greedy = rmsd_utils._compute_greedy_rmsd(
            {"Fe": metal, "C": ring}, {"Fe": metal, "C": rotated}
        )
        robust = rmsd_utils._compute_robust_rmsd(
            {"Fe": metal, "C": ring}, {"Fe": metal, "C": rotated}
        )

        self.assertGreater(greedy, 1.0, "this rotation should defeat greedy assignment")
        self.assertLessEqual(robust, greedy)
        self.assertLess(robust, 1e-6, "the ring is rigidly rotated, so the true RMSD is 0")

    def test_large_element_group_routes_through_the_metric(self):
        ring = [
            ("C", (2.3 * np.cos(t), 2.3 * np.sin(t), 1.5))
            for t in np.linspace(0, 2 * np.pi, 6, endpoint=False)
        ]
        gen = _complex("Fe", ring)
        rmsd, reason = calculate_tmc_rmsd_detailed(
            _complex("Fe", ring, bonded=False), gen, mol2_bonded=gen
        )
        self.assertIsNone(reason)
        self.assertLess(rmsd, 1e-6)


class TestLegacyWrapper(unittest.TestCase):
    """Three integration scripts consume a bare float; they must keep working."""

    def test_success_returns_a_float(self):
        donors = [("N", (2.0, 0, 0)), ("Cl", (0, 0, 2.3))]
        gen = _complex("Fe", donors)
        value = calculate_tmc_rmsd(_complex("Fe", donors, bonded=False), gen, mol2_bonded=gen)
        self.assertIsInstance(float(value), float)
        self.assertLess(value, 1.0)

    def test_mapping_failure_returns_a_sentinel_above_the_floor(self):
        gen = _complex("Pd", [("N", (2.11, 0, 0)), ("N", (0, 0, 2.05))])
        mol_in = _complex("Pd", [("N", (2.11, 0, 0)), ("N", (0, 0, 4.6))], bonded=False)
        value = calculate_tmc_rmsd(mol_in, gen, mol2_bonded=gen)
        self.assertGreaterEqual(value, SENTINEL_FLOOR)

    def test_missing_metal_returns_the_no_metal_sentinel(self):
        organic = _complex("C", [("O", (1.2, 0, 0))])
        self.assertEqual(calculate_tmc_rmsd(organic, organic, mol2_bonded=organic), 999.0)


if __name__ == "__main__":
    unittest.main()
