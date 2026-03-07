import unittest
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.oin.inline import OINInlineHandler

class TestOINInlineHandler(unittest.TestCase):
    def test_parse_inline_simple(self):
        inline = "[Pt_SPL].[Cl]{0}"
        smiles, geo, vectors = OINInlineHandler.parse_inline_string(inline)
        
        self.assertEqual(geo, "SPL")
        # Need to handle standardizing SMILES. [Pt].[Cl]
        self.assertIn("[Pt]", smiles) 
        self.assertIn("Cl", smiles)
        
        # Vectors: LigandRank (1), AtomIdx (0), Slot (0)
        self.assertEqual(len(vectors), 1)
        self.assertEqual(vectors[0], (1, 0, 0))

    def test_parse_inline_complex(self):
        # [Pt_SPL].N{0}
        inline = "[Pt_SPL].N{0}"
        smiles, geo, vectors = OINInlineHandler.parse_inline_string(inline)
        
        self.assertEqual(geo, "SPL")
        self.assertEqual(smiles, "[Pt].N")
        self.assertEqual(vectors[0], (1, 0, 0, False, 1))

    def test_generate_inline_simple(self):
        # Mocking input V2.4 OIN string
        # [Pt].[Cl] |g:SPL|w:1.0:0|
        
        oin_v2 = "[Pt].[Cl] |g:SPL|w:1.0:0|"
        inline = OINInlineHandler.generate_inline_string(oin_v2)
        
        # Expected: [Pt_SPL].[Cl]{0}
        self.assertIn("[Pt_SPL]", inline)
        self.assertIn("[Cl]{0}", inline)

    def test_generate_inline_unbracketed(self):
        # [Pt].N |g:SPL|w:1.0:0|
        oin_v2 = "[Pt].N |g:SPL|w:1.0:0|"
        inline = OINInlineHandler.generate_inline_string(oin_v2)
        
        # Expected: [Pt_SPL].N{0} (Unbracketed)
        self.assertIn("N{0}", inline)

if __name__ == '__main__':
    unittest.main()
