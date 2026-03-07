import unittest
import os
from oinsmiles.core.translator import XYZToSMILES

class TestRegressionStability(unittest.TestCase):
    def setUp(self):
        self.translator = XYZToSMILES()
        self.integration_dir = "tests/integration"
        
    def check_stability(self, xyz_filename, expected_oin):
        path = os.path.join(self.integration_dir, xyz_filename)
        if not os.path.exists(path):
            self.skipTest(f"Fixture {xyz_filename} not found")
        
        oin = self.translator.convert(path)
        self.assertEqual(oin, expected_oin, f"OIN Mismatch for {xyz_filename}")

    def test_cisplatin(self):
        self.check_stability("CisPlatin.xyz", "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}")

    def test_transplatin(self):
        self.check_stability("TransPlatin.xyz", "[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}")

    def test_fac_irppy3(self):
        # Expected OIN from current stable perception (including @b axial tags)
        expected = "[Ir_OCT].[CH]1[CH][CH]N{0}C(C2CCC[CH][C]{3}2)[CH]1.[CH]1[CH][CH]N{5}C(C2CCC[CH][C]{1}2)[CH]1.[CH]1[CH][CH]N{2}C(C2CCC[CH][C]{4}2)[CH]1|b:5-6:STEREOATROP_CW;25-26:STEREOATROP_CW;45-46:STEREOATROP_CW"
        self.check_stability("fac-Ir(ppy)3.xyz", expected)

    def test_mer_irppy3(self):
        # mer-Ir(ppy)3 with improved axial chirality (@b tags)
        expected = "[Ir_OCT].[CH]1[CH]CN{0}C(C2[CH]CCC[C]{3}2)[CH]1.[CH]1[CH]CN{1}C(C2CCCC[C]{5}2)[CH]1.[CH]1[CH]CC(C2[CH]CC[CH][C]{2}2)N{4}[CH]1|b:5-6:STEREOATROP_CCW;25-26:STEREOATROP_CCW;45-46:STEREOATROP_CW"
        self.check_stability("mer-Ir(ppy)3.xyz", expected)

    def test_cis_ptcl2_en(self):
        self.check_stability("Cis-PtCl2(en).xyz", "[Pt_SPL].C(C[NH2]{0})[NH2]{1}.[Cl]{2}.[Cl]{3}")

    def test_ferrocene(self):
        # Ferrocene captures full ring vectors
        expected = "[Fe_LIN].[CH]{0}1[CH]{0}[CH]{0>}[CH]{0}[CH]{0}1.[CH]{1>}1[CH]{1}[CH]{1}[CH]{1}[CH]{1}1"
        self.check_stability("Ferrocene.xyz", expected)

if __name__ == "__main__":
    unittest.main()
