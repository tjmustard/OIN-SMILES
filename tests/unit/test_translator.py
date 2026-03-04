import unittest
from oinsmiles.core.translator import SMILESToXYZ
from oinsmiles.core.graph import BondType

class TestTranslator(unittest.TestCase):
    def setUp(self):
        self.converter = SMILESToXYZ()

    def test_convert_to_graph(self):
        # [Cl] (0), [Pt] (1). Metal is 1. Ligand is 0.
        # v tag: 1.0 -> Metal 1, Ligand 0.
        oin_string = "[Cl].[Pt] |v:1.0:2.0,0.0,0.0;1:0.0,0.0,0.0|"
        graph = self.converter.convert(oin_string)
        
        self.assertEqual(len(graph.atoms), 2)
        # Check coordinates (approximate float comparison)
        self.assertAlmostEqual(graph.atoms[0].coords[0], 2.0)
        
        # Check dative bond
        # Atom 0 should have a dative bond to Atom 1
        self.assertIn(1, graph.atoms[0].bonds)
        self.assertEqual(graph.atoms[0].bonds[1], BondType.DATIVE)

if __name__ == '__main__':
    unittest.main()
