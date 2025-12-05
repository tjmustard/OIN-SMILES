import unittest
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.generation.oin_parser import OINParser

class TestOINParser(unittest.TestCase):
    def test_parse_simple(self):
        oin = "[Pt].[Cl] |w:1:1,0,0|"
        parser = OINParser()
        parsed = parser.parse(oin)
        self.assertEqual(parsed.smiles, "[Pt].[Cl]")
        self.assertEqual(len(parsed.vectors), 1)
        self.assertEqual(parsed.vectors[0].atom_idx, 1)
        self.assertEqual(parsed.vectors[0].fragment_idx, 1) # [Pt] is 0, [Cl] is 1
        self.assertEqual(parsed.vectors[0].atom_in_fragment_idx, 0)

    def test_parse_complex(self):
        oin = "[Pt].[NH2]CC[NH2] |w:1:1,0,0;4:0,1,0|"
        # [Pt] -> 0
        # [NH2]CC[NH2] -> 1. Atoms: N(0), C(1), C(2), N(3). Global: 1, 2, 3, 4.
        parser = OINParser()
        parsed = parser.parse(oin)
        
        self.assertEqual(len(parsed.vectors), 2)
        
        v1 = parsed.vectors[0]
        self.assertEqual(v1.atom_idx, 1)
        self.assertEqual(v1.fragment_idx, 1)
        self.assertEqual(v1.atom_in_fragment_idx, 0)
        
        v2 = parsed.vectors[1]
        self.assertEqual(v2.atom_idx, 4)
        self.assertEqual(v2.fragment_idx, 1)
        self.assertEqual(v2.atom_in_fragment_idx, 3)

if __name__ == '__main__':
    unittest.main()
