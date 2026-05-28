"""Regression stability tests: verify OIN encoding stability post-refactor.

Tests all 6 baseline complexes (cisplatin, transplatin, cis-PtCl₂(en), ferrocene,
fac/mer-Ir(ppy)₃) still encode to their known OIN baselines after the full
ChiralEncoding + MolassemblerAdapter refactor.
"""

import unittest
from pathlib import Path

from oinsmiles import XYZToSMILES
from .test_helpers import get_fixture_path


# Baseline OIN strings for each regression test case.
# These are the known correct outputs from the current v0.2.0+ implementation.
BASELINE_OINS = {
    "cisplatin": "[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
    "transplatin": "[Pt@SP2_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}",
    "cis_ptcl2en": "[Pt@SP1_SPL].[NH2]{0}CC[NH2]{1}.[Cl]{2}.[Cl]{3}",
    "ferrocene": "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1",
    "fac_irppy3": "[Ir@OH10_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1",
    "mer_irppy3": "[Ir@OH2_OCT].c{0}1ccccc1-c1ccccn{3}1.c{1}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{4}1",
}


class TestRegressionStability(unittest.TestCase):
    """Verify OIN encoding stability for baseline complexes post-refactor."""

    def test_cisplatin_stability(self):
        """Cisplatin should encode to baseline OIN."""
        xyz_path = get_fixture_path("cisplatin.xyz")
        result = XYZToSMILES().convert(str(xyz_path))
        self.assertEqual(result, BASELINE_OINS["cisplatin"])

    def test_transplatin_stability(self):
        """Transplatin should encode to baseline OIN."""
        xyz_path = get_fixture_path("transplatin.xyz")
        result = XYZToSMILES().convert(str(xyz_path))
        self.assertEqual(result, BASELINE_OINS["transplatin"])

    def test_cis_ptcl2en_stability(self):
        """cis-PtCl₂(en) should encode to baseline OIN."""
        xyz_path = get_fixture_path("cis_ptcl2en.xyz")
        result = XYZToSMILES().convert(str(xyz_path))
        self.assertEqual(result, BASELINE_OINS["cis_ptcl2en"])

    def test_ferrocene_stability(self):
        """Ferrocene should encode to baseline OIN."""
        xyz_path = get_fixture_path("ferrocene.xyz")
        result = XYZToSMILES().convert(str(xyz_path))
        self.assertEqual(result, BASELINE_OINS["ferrocene"])

    def test_fac_irppy3_stability(self):
        """fac-Ir(ppy)₃ should encode to baseline OIN."""
        xyz_path = get_fixture_path("fac_irppy3.xyz")
        result = XYZToSMILES().convert(str(xyz_path))
        self.assertEqual(result, BASELINE_OINS["fac_irppy3"])

    def test_mer_irppy3_stability(self):
        """mer-Ir(ppy)₃ should encode to baseline OIN."""
        xyz_path = get_fixture_path("mer_irppy3.xyz")
        result = XYZToSMILES().convert(str(xyz_path))
        self.assertEqual(result, BASELINE_OINS["mer_irppy3"])


if __name__ == "__main__":
    unittest.main()
