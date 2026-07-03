import os
import sys
import unittest

# Add src to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.core.graph import BondType, TMCGraph


class TestTMCGraph(unittest.TestCase):
    def setUp(self):
        self.graph = TMCGraph()

    def test_add_atom(self):
        idx = self.graph.add_atom("Pt", (0.0, 0.0, 0.0))
        self.assertEqual(idx, 0)
        self.assertEqual(len(self.graph.atoms), 1)
        self.assertEqual(self.graph.atoms[0].symbol, "Pt")
        self.assertEqual(self.graph.atoms[0].coords, (0.0, 0.0, 0.0))

        idx2 = self.graph.add_atom("Cl", (2.0, 0.0, 0.0))
        self.assertEqual(idx2, 1)
        self.assertEqual(len(self.graph.atoms), 2)

    def test_add_bond(self):
        idx1 = self.graph.add_atom("Pt", (0.0, 0.0, 0.0))
        idx2 = self.graph.add_atom("Cl", (2.0, 0.0, 0.0))

        self.graph.add_bond(idx1, idx2, BondType.COVALENT)

        atom1 = self.graph.atoms[idx1]
        atom2 = self.graph.atoms[idx2]

        self.assertIn(atom2, atom1.neighbors)
        self.assertIn(atom1, atom2.neighbors)
        self.assertEqual(atom1.bonds[idx2], BondType.COVALENT)
        self.assertEqual(atom2.bonds[idx1], BondType.COVALENT)

    def test_set_metal_center(self):
        idx = self.graph.add_atom("Pt", (0.0, 0.0, 0.0))
        self.graph.set_metal_center(idx)
        self.assertEqual(self.graph.metal_index, idx)

    def test_calculate_relative_coords(self):
        self.graph.add_atom("A", (0.0, 0.0, 0.0))
        self.graph.add_atom("B", (2.0, 0.0, 0.0))

        # Centroid should be (1.0, 0.0, 0.0)
        self.graph.calculate_relative_coords()

        rel_A = self.graph.relative_coords[0]
        rel_B = self.graph.relative_coords[1]

        self.assertAlmostEqual(rel_A[0], -1.0)
        self.assertAlmostEqual(rel_B[0], 1.0)


if __name__ == "__main__":
    unittest.main()
