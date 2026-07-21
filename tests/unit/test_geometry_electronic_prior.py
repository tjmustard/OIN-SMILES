"""Electronic (ligand-field) geometry prior for the OIN metal-geometry label.

When a coordination sphere is geometrically ambiguous (a distorted low-CN sphere that
fits two templates about equally poorly), the discrete OIN geometry label (e.g. Ni_SPL
vs Ni_TET) is broken by the metal's d-electron count instead of by knife-edge geometry.
This keeps the label stable across force-field-distorted conformers of the same molecule
(e.g. a square-planar d8 Ni(II) that GFN-FF puckers toward tau4~0.5). See
oin_aligner._electronic_geometry_override.
"""

import unittest
from pathlib import Path

import numpy as np
from rdkit import Chem

from oinsmiles import XYZToSMILES
from oinsmiles.utils.oin_aligner import (
    _electronic_geometry_override,
    classify_coordination_geometry,
    coordination_geometry_fit,
    metal_d_electron_count,
)
from oinsmiles.utils.xyz2mol import get_tmc_mol

FIX = Path(__file__).resolve().parents[1] / "fixtures"

# Real donor vectors from a GFN-FF-distorted BEPCAC (Ni) conformer: both templates fit
# poorly and nearly-equally (SPL 0.614 / TET 0.621) -- the ambiguity band.
DISTORTED_D8_DONORS = [
    [1.7110, -0.3269, 0.6943],
    [-1.6337, 0.9068, 0.1931],
    [-0.8193, -1.7090, 0.3421],
    [0.7184, 1.1034, -1.4081],
]


class MetalDElectronCountTest(unittest.TestCase):
    def _dn(self, fname):
        mol = get_tmc_mol(str(FIX / fname), 0)[0]
        return metal_d_electron_count(mol)

    def test_d8(self):
        self.assertEqual(self._dn("CisPlatin.xyz"), 8)  # Pt(II)

    def test_d0(self):
        self.assertEqual(self._dn("TiCl4.xyz"), 0)  # Ti(IV)

    def test_d6(self):
        self.assertEqual(self._dn("Ferrocene.xyz"), 6)  # Fe(II)

    def test_non_dblock_returns_none(self):
        self.assertIsNone(metal_d_electron_count(Chem.MolFromSmiles("CCO")))


class OverrideDecisionTest(unittest.TestCase):
    AMBIGUOUS = {"SPL": 0.70, "TET": 0.58}  # both poor, close -> genuinely between
    DECISIVE_SPL = {"SPL": 0.14, "TET": 1.22}  # clean square planar

    def test_d8_ambiguous_overrides_to_spl(self):
        self.assertEqual(_electronic_geometry_override(4, "TET", self.AMBIGUOUS, 8), "SPL")

    def test_d0_ambiguous_overrides_to_tet(self):
        self.assertEqual(
            _electronic_geometry_override(4, "SPL", {"SPL": 0.60, "TET": 0.66}, 0), "TET"
        )

    def test_decisive_geometry_not_overridden(self):
        # A clean square plane is never flipped, even for d8 (pref == current).
        self.assertIsNone(_electronic_geometry_override(4, "SPL", self.DECISIVE_SPL, 8))

    def test_decisive_tet_not_overridden(self):
        # Winner fits well (min rmsd below the distortion floor) -> trust geometry.
        self.assertIsNone(_electronic_geometry_override(4, "TET", {"SPL": 0.20, "TET": 0.10}, 8))

    def test_unknown_dcount_no_override(self):
        self.assertIsNone(_electronic_geometry_override(4, "TET", self.AMBIGUOUS, None))

    def test_dcount_without_prior_cell(self):
        # d6 CN4 has no strong spin-independent preference in the table.
        self.assertIsNone(_electronic_geometry_override(4, "TET", self.AMBIGUOUS, 6))

    def test_pref_equals_current_is_noop(self):
        self.assertIsNone(_electronic_geometry_override(4, "SPL", self.AMBIGUOUS, 8))


class DistortedGeometryTest(unittest.TestCase):
    def test_real_distorted_donors_are_ambiguous_and_resolve_to_spl(self):
        vecs = [np.array(v) for v in DISTORTED_D8_DONORS]
        rmsd_spl = coordination_geometry_fit(vecs, "SPL")
        rmsd_tet = coordination_geometry_fit(vecs, "TET")
        # Geometry alone is ambiguous: both poor, nearly tied.
        self.assertGreater(min(rmsd_spl, rmsd_tet), 0.35)
        self.assertLess(abs(rmsd_spl - rmsd_tet), 0.25)
        self.assertIn(classify_coordination_geometry(vecs), ("SPL", "TET"))
        # d8 resolves to square planar regardless of the geometric coin-flip.
        rmsds = {"SPL": rmsd_spl, "TET": rmsd_tet}
        self.assertEqual(_electronic_geometry_override(4, "TET", rmsds, 8), "SPL")
        # ...and without the d-count, the label stays at the (unstable) geometric winner.
        self.assertIsNone(_electronic_geometry_override(4, "TET", rmsds, None))


class CleanCaseRegressionTest(unittest.TestCase):
    """The prior must not perturb decisive geometries."""

    def test_cisplatin_still_spl(self):
        self.assertIn("Pt_SPL", XYZToSMILES().convert(str(FIX / "CisPlatin.xyz")))

    def test_bepcac_real_still_spl(self):
        oin = XYZToSMILES().convert(str(FIX / "conformer_set" / "BEPCAC_comp_0.xyz"))
        self.assertIn("Ni_SPL", oin)


if __name__ == "__main__":
    unittest.main()
