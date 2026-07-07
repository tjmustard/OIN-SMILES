"""Unit tests for eta-ring winding generalization.

The XYZ<->OIN winding machinery must work for *any* number of eta ligands
(2, 3, 4, ...) and *any* hapticity (eta-2/eta-3/eta-5/eta-6), not just the
bis-indenyl ansa-metallocenes in the curated fixtures. These tests exercise the
two load-bearing, fixture-free pieces:

* ``oin.winding.signed_circulation`` -- the per-ring sign, defined for >= 3
  ring atoms (so eta-3 allyl works) and degenerate (default forward) for < 3.
* ``generation.metallogen_adapter._eta_winding_multiset`` -- the order-
  independent per-ring winding multiset used by winding-aware conformer
  selection and (mirrored) by the round-trip harness, which must scale to N
  rings and cleanly separate rac / meso / enantiomer.
"""

import unittest

import numpy as np

from oinsmiles.generation.metallogen_adapter import _eta_winding_multiset
from oinsmiles.oin.winding import signed_circulation


class TestSignedCirculationHapticity(unittest.TestCase):
    """signed_circulation across hapticities."""

    AXIS = np.array([0.0, 0.0, 1.0])  # metal -> centroid, outward (+z)

    def test_eta3_allyl_has_a_definite_winding(self):
        """An eta-3 allyl (3 donor atoms) has a well-defined winding sign --
        eta-3 ligands are a first-class case, not a degenerate one."""
        allyl = np.array([[-1.0, 0.0, 1.0], [0.0, 0.6, 1.0], [1.0, 0.0, 1.0]])
        sign = signed_circulation(allyl, 0, self.AXIS)
        self.assertIn(sign, (">", "<"))
        # Reversing the atom traversal flips the sign (it encodes a real
        # rotational sense, not a constant).
        self.assertNotEqual(sign, signed_circulation(allyl[::-1].copy(), 0, self.AXIS))

    def test_eta2_is_degenerate_forward(self):
        """An eta-2 group (2 atoms) has no rotational sense -> forward default."""
        edge = np.array([[-1.0, 0.0, 1.0], [1.0, 0.0, 1.0]])
        self.assertEqual(signed_circulation(edge, 0, self.AXIS), ">")

    def test_eta5_ring_winding_is_order_sensitive(self):
        """An eta-5 ring has a definite, traversal-direction-dependent sign."""
        ring = np.array(
            [
                [np.cos(t), np.sin(t), 1.0]
                for t in np.linspace(0.0, 2.0 * np.pi, 5, endpoint=False)
            ]
        )
        forward = signed_circulation(ring, 0, self.AXIS)
        reverse = signed_circulation(ring[::-1].copy(), 0, self.AXIS)
        self.assertIn(forward, (">", "<"))
        self.assertNotEqual(forward, reverse)

    def test_axis_direction_flips_sign(self):
        """Flipping the outward axis (which face points at the metal) flips the
        winding -- the property winding exists to encode."""
        ring = np.array(
            [
                [np.cos(t), np.sin(t), 1.0]
                for t in np.linspace(0.0, 2.0 * np.pi, 5, endpoint=False)
            ]
        )
        up = signed_circulation(ring, 0, self.AXIS)
        down = signed_circulation(ring, 0, -self.AXIS)
        self.assertNotEqual(up, down)


class TestEtaWindingMultiset(unittest.TestCase):
    """_eta_winding_multiset scales to N rings and separates stereoisomers."""

    def test_empty_or_none(self):
        self.assertEqual(_eta_winding_multiset(None), [])
        self.assertEqual(_eta_winding_multiset(""), [])
        self.assertEqual(_eta_winding_multiset("[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"), [])

    def test_scales_to_three_and_four_eta_rings(self):
        self.assertEqual(
            _eta_winding_multiset("[M_OCT].a{0>}.b{1<}.c{2>}.X{3}"), ["<", ">", ">"]
        )
        self.assertEqual(
            _eta_winding_multiset("[M].a{0<}.b{1<}.c{2>}.d{3>}"),
            ["<", "<", ">", ">"],
        )

    def test_rac_meso_enantiomer_are_distinct_multisets(self):
        rac = _eta_winding_multiset("a{0>}b{1>}")
        meso = _eta_winding_multiset("a{0<}b{1>}")
        rac_ent = _eta_winding_multiset("a{0<}b{1<}")
        self.assertEqual(rac, [">", ">"])
        self.assertEqual(meso, ["<", ">"])
        self.assertEqual(rac_ent, ["<", "<"])
        # rac vs meso is a diastereomer difference; rac vs its enantiomer is an
        # enantiomer difference. Both must be visible (multisets differ).
        self.assertNotEqual(rac, meso)
        self.assertNotEqual(rac, rac_ent)

    def test_multiset_is_order_independent(self):
        """The two equivalent labelings of a meso structure ({0<}{1>} and
        {0>}{1<}) collapse to the same multiset -- they are one molecule."""
        self.assertEqual(
            _eta_winding_multiset("a{0<}b{1>}"),
            _eta_winding_multiset("a{0>}b{1<}"),
        )


if __name__ == "__main__":
    unittest.main()
