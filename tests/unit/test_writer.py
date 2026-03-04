import unittest
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.oin.writer import OINWriter

class TestOINWriter(unittest.TestCase):
    def setUp(self):
        self.writer = OINWriter()

    def test_write_simple(self):
        smiles = "[Pt].[Cl]"
        coords = [(1, 2.0, 0.0, 0.0)] # Ligand at x=2
        oin = self.writer.write(smiles, coords)
        
        expected = "[Pt].[Cl] |w:1:2,0,0|"
        self.assertIn("w:1:2,0,0", oin)
        self.assertTrue(oin.startswith(smiles))

    def test_write_complex_tags(self):
        smiles = "[Pt].[NH3]"
        coords = [(1, 0.0, 2.0, 0.0)]
        dative = [(1, 0)] # Ligand 1 -> Metal 0
        geo = "SP_cis"
        
        oin = self.writer.write(smiles, coords, dative_bonds=dative, geometry_label=geo)
        
        self.assertIn("g:SP_cis", oin)
        self.assertIn("d:1.0", oin)
        self.assertIn("w:1:0,2,0", oin)

    def test_write_no_tags(self):
        smiles = "C"
        oin = self.writer.write(smiles, [])
        self.assertEqual(oin, "C")

if __name__ == '__main__':
    unittest.main()
