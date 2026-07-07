"""Unit tests for geometry-code-aware conformer selection.

Covers the three layers of the BDNN square-plane stability fix:
  * ``classify_coordination_geometry`` -- pure donor-vector -> OIN code matcher.
  * ``_perceive_geo_code`` -- metal/donor discovery + haptic gate on a contract mol.
  * ``_select_by_geometry`` -- prefer the target-geometry conformer over the
    lowest-energy one, with a non-regressive lowest-energy fallback.
"""

import os
import sys
import types
import unittest
from unittest import mock

from rdkit import Chem
from rdkit.Geometry import Point3D

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation import metallogen_adapter as MA
from oinsmiles.utils.oin_aligner import (
    classify_coordination_geometry,
    coordination_geometry_fit,
)

# Ideal template donor directions (metal at origin), scaled off unit length to
# confirm the matcher is scale-invariant.
_SPL = [[2, 0, 0], [0, 2, 0], [-2, 0, 0], [0, -2, 0]]
_TPY = [[0, 0, 2], [0, 2, 0], [1.7320508, -1, 0], [-1.7320508, -1, 0]]
_TET = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
_OCT = [[2, 0, 0], [-2, 0, 0], [0, 2, 0], [0, -2, 0], [0, 0, 2], [0, 0, -2]]
# A mildly puckered square-plane: still classifies SPL, but a looser fit than _SPL.
_SPL_PUCKERED = [[2, 0, 0.4], [0, 2, -0.4], [-2, 0, 0.4], [0, -2, -0.4]]


def _make_complex(metal_sym, donor_positions, donor_sym="Cl"):
    """A minimal RDKit mol: one metal at origin with DATIVE bonds to donors."""
    rw = Chem.RWMol()
    mi = rw.AddAtom(Chem.Atom(metal_sym))
    conf = Chem.Conformer(1 + len(donor_positions))
    conf.SetAtomPosition(mi, Point3D(0.0, 0.0, 0.0))
    for x, y, z in donor_positions:
        di = rw.AddAtom(Chem.Atom(donor_sym))
        rw.AddBond(mi, di, Chem.BondType.DATIVE)
        conf.SetAtomPosition(di, Point3D(float(x), float(y), float(z)))
    m = rw.GetMol()
    m.AddConformer(conf, assignId=True)
    return m


class TestClassifyCoordinationGeometry(unittest.TestCase):
    def test_square_planar(self):
        self.assertEqual(classify_coordination_geometry(_SPL), "SPL")

    def test_trigonal_pyramidal(self):
        self.assertEqual(classify_coordination_geometry(_TPY), "TPY")

    def test_tetrahedral(self):
        self.assertEqual(classify_coordination_geometry(_TET), "TET")

    def test_octahedral(self):
        self.assertEqual(classify_coordination_geometry(_OCT), "OCT")

    def test_haptic_returns_none(self):
        # 10 coplanar donors: no discrete template of that coordination number.
        import numpy as np

        many = [
            [np.cos(t), np.sin(t), 0.0]
            for t in np.linspace(0, 2 * np.pi, 10, endpoint=False)
        ]
        self.assertIsNone(classify_coordination_geometry(many))


class TestCoordinationGeometryFit(unittest.TestCase):
    def test_ideal_fit_is_near_zero(self):
        self.assertLess(coordination_geometry_fit(_SPL, "SPL"), 1e-6)

    def test_puckered_fits_worse_than_ideal(self):
        # Both are SPL, but the puckered one is a looser fit -- the signal the
        # selector ranks on.
        self.assertEqual(classify_coordination_geometry(_SPL_PUCKERED), "SPL")
        self.assertGreater(
            coordination_geometry_fit(_SPL_PUCKERED, "SPL"),
            coordination_geometry_fit(_SPL, "SPL"),
        )

    def test_unknown_code_is_infinite(self):
        self.assertEqual(coordination_geometry_fit(_SPL, "ZZZ"), float("inf"))


class TestPerceiveGeoCode(unittest.TestCase):
    def test_square_planar(self):
        self.assertEqual(MA._perceive_geo_code(_make_complex("Pd", _SPL), 4), "SPL")

    def test_trigonal_pyramidal(self):
        self.assertEqual(MA._perceive_geo_code(_make_complex("Pd", _TPY), 4), "TPY")

    def test_octahedral(self):
        self.assertEqual(MA._perceive_geo_code(_make_complex("Fe", _OCT), 6), "OCT")

    def test_gate_rejects_mismatched_coordination_number(self):
        # 4 donors but we asked for 6 -> eta/haptic guard -> None (no selection).
        self.assertIsNone(MA._perceive_geo_code(_make_complex("Pd", _SPL), 6))

    def test_no_metal_returns_none(self):
        self.assertIsNone(MA._perceive_geo_code(_make_complex("C", _SPL), 4))


class TestSelectByGeometry(unittest.TestCase):
    """``build_contract_mol`` is patched to identity so the real ``_perceive_geo_code``
    classifies the synthetic RDKit mols we pass in as the conformer pool."""

    def _parsed(self, geo_code):
        return types.SimpleNamespace(geo_code=geo_code)

    def test_prefers_target_over_lower_energy(self):
        tpy = _make_complex("Pd", _TPY)  # index 0 == "lowest energy"
        spl = _make_complex("Pd", _SPL)
        with mock.patch.object(MA, "build_contract_mol", side_effect=lambda p, m: m):
            chosen, cmol = MA._select_by_geometry(self._parsed("SPL"), [tpy, spl])
        self.assertIs(chosen, spl)
        self.assertIs(cmol, spl)

    def test_prefers_cleaner_geometry_over_lower_energy_same_code(self):
        # Both classify as SPL; the lower-energy one (index 0) is puckered, the
        # higher-energy one is a clean square-plane. Rank by template fit -> clean.
        puckered = _make_complex("Pd", _SPL_PUCKERED)  # "lowest energy"
        clean = _make_complex("Pd", _SPL)
        with mock.patch.object(MA, "build_contract_mol", side_effect=lambda p, m: m):
            chosen, _ = MA._select_by_geometry(self._parsed("SPL"), [puckered, clean])
        self.assertIs(chosen, clean)

    def test_sqp_normalizes_to_spl(self):
        tpy = _make_complex("Pd", _TPY)
        spl = _make_complex("Pd", _SPL)
        with mock.patch.object(MA, "build_contract_mol", side_effect=lambda p, m: m):
            chosen, _ = MA._select_by_geometry(self._parsed("SQP"), [tpy, spl])
        self.assertIs(chosen, spl)

    def test_falls_back_to_lowest_energy_when_no_match(self):
        tpy_a = _make_complex("Pd", _TPY)
        tpy_b = _make_complex("Pd", _TPY)
        with mock.patch.object(MA, "build_contract_mol", side_effect=lambda p, m: m):
            chosen, _ = MA._select_by_geometry(self._parsed("SPL"), [tpy_a, tpy_b])
        self.assertIs(chosen, tpy_a)

    def test_unknown_geometry_code_falls_back(self):
        first = _make_complex("Pd", _SPL)
        second = _make_complex("Pd", _TPY)
        with mock.patch.object(MA, "build_contract_mol", side_effect=lambda p, m: m):
            chosen, _ = MA._select_by_geometry(self._parsed("ZZZ"), [first, second])
        self.assertIs(chosen, first)


if __name__ == "__main__":
    unittest.main()
