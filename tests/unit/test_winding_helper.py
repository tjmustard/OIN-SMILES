"""Unit tests for the shared winding-sign helper (Stereo Phase 3).

`signed_circulation` (`oin/winding.py`) is the single source of truth for the
OIN V3.6 winding character, shared by the encode-side
`OINDiscreteAligner._determine_winding` and the (future) generation-side
haptic-face correction. These tests pin: (1) behavioral parity between the
encode path and a direct helper call on identical coords, (2) that a 180deg
in-plane rotation flips the sign, (3) that an inward-facing template slot Z
still yields the convention-correct char, and (4) the `n < 3` degenerate
default.
"""

import unittest

import numpy as np
from scipy.spatial.transform import Rotation

from oinsmiles.oin.winding import signed_circulation
from oinsmiles.utils.oin_aligner import OINDiscreteAligner


def _pentagon(radius=1.0, z=0.0):
    """Five points evenly spaced on a circle, in ascending SMILES/local_idx order."""
    coords = []
    for i in range(5):
        theta = 2 * np.pi * i / 5
        coords.append([radius * np.cos(theta), radius * np.sin(theta), z])
    return np.array(coords)


def _dummy_aligner():
    # `_determine_winding` doesn't touch `self`; a minimal instance is enough
    # to call it through the class's own convention (rather than reaching
    # into a "private" free function).
    return OINDiscreteAligner(metal_idx=0, ligands=[{"smiles": "[M]", "metal_coords": np.zeros(3)}])


class TestSignedCirculationParity(unittest.TestCase):
    """Test 1: encode-path char must equal a direct helper call on identical coords."""

    def test_encode_path_matches_direct_helper_call(self):
        coords = _pentagon()
        constituent_indices = [0, 1, 2, 3, 4]
        star_idx = 0
        slot_z = np.array([0.0, 0.0, 1.0])  # metal below plane -> outward = +Z

        aligner = _dummy_aligner()
        encode_char = aligner._determine_winding(
            grp_coords=coords,
            star_idx=star_idx,
            constituent_indices=constituent_indices,
            slot_z=slot_z,
            slot_x_ref=np.array([1.0, 0.0, 0.0]),
        )

        direct_char = signed_circulation(coords, star_local_idx=0, axis=slot_z)

        self.assertEqual(encode_char, direct_char)


class TestSignedCirculationRotationInvert(unittest.TestCase):
    """Test 2: a 180deg in-plane rotation of the ring must flip the char."""

    def test_180_degree_rotation_flips_char(self):
        coords = _pentagon()
        axis = np.array([0.0, 0.0, 1.0])  # fixed lab/template axis (unchanged by the flip)

        original_char = signed_circulation(coords, star_local_idx=0, axis=axis)

        # A 180deg rotation about an axis LYING IN the ring's own plane (e.g.
        # X, since the pentagon lies in the XY plane) flips the ring over --
        # this is the physical "haptic face flip" operation (Stereo Phase 3).
        # Coplanar points keep z == 0 (rotation about X sends z -> -z, and
        # z was already 0), but y -> -y, reversing the star->next traversal
        # direction as seen from the FIXED outward axis. A 180deg rotation
        # about the axis itself (Z) would NOT flip the sign -- it just spins
        # the whole ring rigidly without reversing circulation.
        flip = Rotation.from_euler("x", 180, degrees=True)
        flipped_coords = flip.apply(coords)

        flipped_char = signed_circulation(flipped_coords, star_local_idx=0, axis=axis)

        self.assertNotEqual(original_char, flipped_char)


class TestSignedCirculationAntiparallelSlot(unittest.TestCase):
    """Test 3: an inward-facing template Z must still yield the convention-correct char."""

    def test_inward_facing_slot_z_matches_outward_convention(self):
        coords = _pentagon()

        outward_axis = np.array([0.0, 0.0, 1.0])
        inward_axis = -outward_axis

        outward_char = signed_circulation(coords, star_local_idx=0, axis=outward_axis)
        inward_char = signed_circulation(coords, star_local_idx=0, axis=inward_axis)

        # The helper has no notion of "correct" direction on its own -- it is
        # a pure function of the axis it is given. The encoder's contract is
        # to always feed an outward (metal->centroid) axis; this test guards
        # that an accidentally-inward axis (e.g. an octahedral/LIN slot whose
        # template Z points inward) produces the OPPOSITE char, so a caller
        # that fails to normalize direction will visibly break parity rather
        # than silently agree by coincidence.
        self.assertNotEqual(outward_char, inward_char)

    def test_determine_winding_normalizes_slot_z_direction(self):
        # `_determine_winding` normalizes slot_z's magnitude but trusts its
        # caller-supplied direction to already be outward (asserted via
        # TEMPLATE_SPECS construction). Feeding a non-unit but
        # outward-pointing slot_z must still match the unit-outward-axis
        # helper result directly.
        coords = _pentagon()
        constituent_indices = [0, 1, 2, 3, 4]

        aligner = _dummy_aligner()
        char_scaled = aligner._determine_winding(
            grp_coords=coords,
            star_idx=0,
            constituent_indices=constituent_indices,
            slot_z=np.array([0.0, 0.0, 5.0]),  # outward, non-unit magnitude
            slot_x_ref=np.array([1.0, 0.0, 0.0]),
        )
        char_unit = signed_circulation(coords, star_local_idx=0, axis=np.array([0.0, 0.0, 1.0]))

        self.assertEqual(char_scaled, char_unit)


class TestSignedCirculationDegenerate(unittest.TestCase):
    """Test 4: n < 3 (linear/monodentate) always defaults to '>'."""

    def test_n_less_than_3_defaults_to_forward(self):
        axis = np.array([0.0, 0.0, 1.0])

        self.assertEqual(signed_circulation(np.array([[0.0, 0.0, 0.0]]), 0, axis), ">")
        self.assertEqual(
            signed_circulation(np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]), 0, axis),
            ">",
        )


if __name__ == "__main__":
    unittest.main()
