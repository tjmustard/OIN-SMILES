import unittest
import sys
import os
import numpy as np

# Add src to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.core.graph import TMCGraph, AtomNode, BondType

class TestTMCGraph(unittest.TestCase):
    def setUp(self):
        self.graph = TMCGraph()

    def test_add_atom(self):
        idx = self.graph.add_atom("Pt", (0.0, 0.0, 0.0))
        self.assertEqual(idx, 0)
        self.assertEqual(len(self.graph.atoms), 1)
        self.assertEqual(self.graph.atoms[0].element, "Pt")
        self.assertEqual(self.graph.atoms[0].coords, (0.0, 0.0, 0.0))

        idx2 = self.graph.add_atom("Cl", (2.0, 0.0, 0.0))
        self.assertEqual(idx2, 1)
        self.assertEqual(len(self.graph.atoms), 2)

    def test_add_bond(self):
        idx1 = self.graph.add_atom("Pt", (0.0, 0.0, 0.0))
        idx2 = self.graph.add_atom("Cl", (2.0, 0.0, 0.0))
        
        self.graph.add_bond(idx1, idx2, BondType.SINGLE)
        
        self.assertEqual(len(self.graph.bonds), 1)
        self.assertEqual(self.graph.bonds[0], (idx1, idx2, BondType.SINGLE))

if __name__ == '__main__':
    unittest.main()
