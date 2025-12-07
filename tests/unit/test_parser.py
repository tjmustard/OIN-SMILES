import unittest
from oinsmiles.oin.parser import OINParser

class TestOINParser(unittest.TestCase):
    def setUp(self):
        self.parser = OINParser()

    def test_parse_simple(self):
        oin_string = "[Cl].[Pt]([Cl])([NH3])[NH3] |w:0:1.5,0,0;1:-1.5,0,0;2:0,1.5,0;3:0,-1.5,0;4:0,0,0|"
        smiles, tags = self.parser.parse(oin_string)
        self.assertEqual(smiles, "[Cl].[Pt]([Cl])([NH3])[NH3]")
        self.assertIn('w', tags)
        self.assertEqual(tags['w'], "0:1.5,0,0;1:-1.5,0,0;2:0,1.5,0;3:0,-1.5,0;4:0,0,0")

    def test_parse_multiple_tags(self):
        oin_string = "[NH3].[Pt] |w:0:2.0,2.0,0;1:0,0,0|d:0.1|g:SP_cis|"
        smiles, tags = self.parser.parse(oin_string)
        self.assertEqual(smiles, "[NH3].[Pt]")
        self.assertIn('w', tags)
        self.assertIn('d', tags)
        self.assertIn('g', tags)
        self.assertEqual(tags['d'], "0.1")
        self.assertEqual(tags['g'], "SP_cis")

    def test_parse_coordinates(self):
        w_tag = "0:1.0,2.0,3.0;1:4.0,5.0,6.0"
        coords = self.parser.parse_coordinates(w_tag)
        self.assertEqual(len(coords), 2)
        self.assertEqual(coords[0], (0, 1.0, 2.0, 3.0))
        self.assertEqual(coords[1], (1, 4.0, 5.0, 6.0))



if __name__ == '__main__':
    unittest.main()
