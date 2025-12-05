import unittest
from oinsmiles.oin.writer import OINWriter

class TestOINWriter(unittest.TestCase):
    def setUp(self):
        self.writer = OINWriter()

    def test_write_simple(self):
        smiles = "[Cl].[Pt]"
        coords = [(0, 1.0, 0.0, 0.0), (1, 0.0, 0.0, 0.0)]
        output = self.writer.write(smiles, coords)
        self.assertEqual(output, "[Cl].[Pt] |w:0:1,0,0;1:0,0,0|")

    def test_write_all_tags(self):
        smiles = "[NH3].[Pt]"
        coords = [(0, 2.0, 0.0, 0.0), (1, 0.0, 0.0, 0.0)]
        dative = [(0, 1)]
        haptic = ["0:1.2.3"]
        geometry = "SP_cis"
        
        output = self.writer.write(smiles, coords, dative_bonds=dative, haptic_groups=haptic, geometry_label=geometry)
        
        self.assertIn("w:0:2,0,0;1:0,0,0", output)
        self.assertIn("d:0.1", output)
        self.assertIn("m:0:1.2.3", output)
        self.assertIn("g:SP_cis", output)
        self.assertIn("|", output)

if __name__ == '__main__':
    unittest.main()
