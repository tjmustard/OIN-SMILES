import unittest
from oinsmiles.core.translator import SMILESToXYZ
from oinsmiles.core.graph import BondType

class TestTranslator(unittest.TestCase):
    def setUp(self):
        self.converter = SMILESToXYZ()

    def test_convert_to_graph(self):
        oin_string = "[Cl].[Pt] |w:0:2.0,0.0,0.0;1:0.0,0.0,0.0|d:0.1|"
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
