"""Guards for SL1 (v0.4.4 acceptance-gate): the stretched-bond metric.

``clash.stretched_bond_count`` is the complement of ``clash.vdw_clash_count``: the vdW gate
flags *non-bonded* pairs pulled too CLOSE, this flags *bonded* pairs pulled too FAR. It exists
because the generate-until-key-exact early-exit needs to reject the "methyl threaded through a
phenyl / interlocked rings tangled a bond out to ~4 Angstrom" pathology, which no clash or
fusion check sees.

Each test fails against pre-SL1 code: ``stretched_bond_count`` / ``mol_stretched_bond_count``
did not exist. The cross-check asserts the two metrics are genuinely complementary -- a
stretched bond is invisible to the clash count and a vdW clash is invisible to the stretch
count -- so neither can be substituted for the other.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

import numpy as np

from oinsmiles.generator3d import clash

C = 6  # carbon Z


class TestStretchedBondCount(unittest.TestCase):
    """``stretched_bond_count`` on hand-built fixtures (no embedding needed)."""

    def test_normal_bond_not_stretched(self):
        # Two carbons at a normal C-C distance (~1.5 A) that ARE bonded: ratio well under 1.5.
        positions = [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)]
        adjacency = [[0, 1], [1, 0]]
        n, worst = clash.stretched_bond_count(positions, [C, C], adjacency)
        self.assertEqual(n, 0)
        # sanity: worst_ratio is the (single) bonded ratio, ~1.5/(2*rcov(C)) < 1.0
        self.assertLess(worst, 1.5)

    def test_stretched_bond_flagged(self):
        # Same bonded pair pulled out to 4.0 A -- ratio ~2.6, well over the 1.5 bar.
        positions = [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0)]
        adjacency = [[0, 1], [1, 0]]
        n, worst = clash.stretched_bond_count(positions, [C, C], adjacency)
        self.assertEqual(n, 1)
        expected = 4.0 / (2 * clash._rcov(C))
        self.assertAlmostEqual(worst, expected, places=6)
        self.assertGreater(worst, clash.STRETCH_RATIO)

    def test_far_pair_but_not_bonded_is_ignored(self):
        # Two atoms 4 A apart but NOT in the adjacency: a stretched bond needs a real bond.
        positions = [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0)]
        adjacency = [[0, 0], [0, 0]]
        n, worst = clash.stretched_bond_count(positions, [C, C], adjacency)
        self.assertEqual(n, 0)
        self.assertEqual(worst, 0.0)  # no bonds -> 0.0 sentinel

    def test_directed_adjacency_is_symmetrized(self):
        # Only the upper entry set: the metric ORs adj with adj.T, so it still sees the bond.
        positions = [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0)]
        adjacency = [[0, 1], [0, 0]]
        n, _worst = clash.stretched_bond_count(positions, [C, C], adjacency)
        self.assertEqual(n, 1)

    def test_shape_mismatch_degrades_to_zero(self):
        positions = [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0)]
        n, worst = clash.stretched_bond_count(positions, [C, C], adjacency=[[0, 1, 0]])
        self.assertEqual((n, worst), (0, 0.0))

    def test_none_adjacency_degrades_to_zero(self):
        positions = [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0)]
        self.assertEqual(clash.stretched_bond_count(positions, [C, C], None), (0, 0.0))

    def test_complementary_to_vdw_clash(self):
        """A stretched bond and a vdW clash are disjoint by construction.

        Fixture: atoms 0-1 are a bonded C-C stretched to 4 A (a stretched bond, but too far
        apart to be a vdW clash); atoms 2-3 are a non-bonded C-C at 2.2 A (a vdW clash --
        outside the 1.3*sum(rcov)~1.98 A bond-perception window but inside the 0.75 overlap
        cutoff -- and not a bond, so not a stretch). The two clusters are 10 A apart so no
        cross-pair interferes.
        """
        positions = [
            (0.0, 0.0, 0.0),  # 0
            (4.0, 0.0, 0.0),  # 1  (bonded to 0, stretched)
            (10.0, 0.0, 0.0),  # 2
            (12.2, 0.0, 0.0),  # 3  (non-bonded to 2, vdW clash)
        ]
        atomic_numbers = [C, C, C, C]
        adjacency = [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
        ]

        n_stretched, _worst = clash.stretched_bond_count(positions, atomic_numbers, adjacency)
        clash_vdw, _severe, _wo = clash.vdw_clash_count(positions, atomic_numbers)

        # Each metric fires exactly once, on its own pair.
        self.assertEqual(n_stretched, 1, "the bonded 0-1 pair is the only stretched bond")
        self.assertEqual(clash_vdw, 1, "the non-bonded 2-3 pair is the only vdW clash")

        # Complementarity: the stretched bond alone has no clash; the clash alone has no stretch.
        stretch_only = clash.stretched_bond_count(
            positions[:2], atomic_numbers[:2], [[0, 1], [1, 0]]
        )
        clash_only, _s, _w = clash.vdw_clash_count(positions[:2], atomic_numbers[:2])
        self.assertEqual((stretch_only[0], clash_only), (1, 0))

        clash_pair, _s2, _w2 = clash.vdw_clash_count(positions[2:], atomic_numbers[2:])
        stretch_pair = clash.stretched_bond_count(
            positions[2:], atomic_numbers[2:], [[0, 0], [0, 0]]
        )
        self.assertEqual((clash_pair, stretch_pair[0]), (1, 0))


class TestMolStretchedBondCount(unittest.TestCase):
    """``mol_stretched_bond_count`` duck-typing + graceful degradation."""

    def test_reads_atom_list_and_adjacency(self):
        class _Atom:
            def __init__(self, xyz, z):
                self._xyz = xyz
                self._z = z

            def get_coordinate(self):
                return self._xyz

            def get_atomic_number(self):
                return self._z

        class _Mol:
            atom_list = [_Atom((0.0, 0.0, 0.0), C), _Atom((4.0, 0.0, 0.0), C)]

            def get_adj_matrix(self):
                return np.array([[0, 1], [1, 0]])

        self.assertEqual(clash.mol_stretched_bond_count(_Mol()), 1)

    def test_bare_object_degrades_to_zero(self):
        # No atom_list -> AttributeError -> 0 (accept rather than drop the candidate).
        self.assertEqual(clash.mol_stretched_bond_count(object()), 0)

    def test_missing_adjacency_degrades_to_zero(self):
        class _Atom:
            def get_coordinate(self):
                return (0.0, 0.0, 0.0)

            def get_atomic_number(self):
                return C

        class _Mol:
            atom_list = [_Atom(), _Atom()]

            def get_adj_matrix(self):
                return None

        self.assertEqual(clash.mol_stretched_bond_count(_Mol()), 0)


if __name__ == "__main__":
    unittest.main()
