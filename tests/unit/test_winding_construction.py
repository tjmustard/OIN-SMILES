"""SL2 oin-direct-winding: deterministic eta-ring winding construction (Part B).

The rigid haptic placer (`embed._place_haptic`) turns the placed ring over with one
180 deg proper rotation when its measured winding disagrees with the OIN's -- so an
eta complex generated under `oin_direct` reproduces the requested winding by
construction (no 16-wide search pool). Covered here:
 - unit level: `_place_haptic` constructs BOTH requested winding chars,
 - end to end: a symmetric ring (ferrocene) and load-bearing rings (halide-face Cp
   diastereomers, a chirality-witness ring) round-trip key-exact with the winding
   multiset preserved -- the load-bearing cases would fail if the flip picked the
   wrong face or used a reflection (which would corrupt the sp3 witness).
"""

import math
import unittest

import numpy as np
from rdkit import Chem

from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generation.metallogen_adapter import _eta_winding_multiset
from oinsmiles.generator3d.embed import _place_haptic
from oinsmiles.oin.compare import canonical_roundtrip_key
from oinsmiles.oin.winding import signed_circulation
from oinsmiles.utils.perception_tmc import get_oin_string

FERROCENE = "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
# Two load-bearing halide-face Cp rings (each ring's face is distinguishable, so its
# '<'/'>' is geometric, not notation): a real winding-construction test.
FE_HALIDE_FACE = (
    "[Fe_LIN].Oc{0}1[cH]{0<}c{0}(Cl)c{0}(Br)c{0}1I.Oc{1}1[cH]{1<}c{1}(I)c{1}(Br)c{1}1Cl"
)
# One load-bearing ring ('<') carrying an sp3 chirality witness + one symmetric ring.
FE_CHIRAL_WITNESS = (
    "[Fe_LIN].[C@H](F)(Cl)c{0<}1[cH]{0}c{0}(Cl)c{0}(Br)c{0}1I."
    "[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
)


def _reencode(result):
    lines = result.xyz.splitlines()
    n = int(lines[0].strip())
    coords = np.array([[float(x) for x in lines[i].split()[1:4]] for i in range(2, 2 + n)])
    assert result.mol is not None, "expected a non-None contract mol"
    return get_oin_string(Chem.Mol(result.mol), coords)


class TestPlaceHapticConstructsWinding(unittest.TestCase):
    """_place_haptic returns a ring whose winding equals the requested char."""

    def _pentagon(self):
        n = 5
        angs = [2 * math.pi * k / n for k in range(n)]
        return np.array([[math.cos(a), math.sin(a), 0.0] for a in angs])

    def test_both_chars_constructed(self):
        P = self._pentagon()
        face = [0, 1, 2, 3, 4]
        v_slot = np.array([0.0, 0.0, 1.0])
        for target in (">", "<"):
            placed = _place_haptic(P.copy(), face, v_slot, 1.5, 0.7, 1.0, winding=(face, 0, target))
            ring = placed[face]
            got = signed_circulation(ring, 0, ring.mean(axis=0))
            self.assertEqual(got, target, f"failed to construct winding {target}")

    def test_no_winding_arg_is_noop_path(self):
        # Omitting winding leaves the arbitrary-face placement (byte-identical path).
        P = self._pentagon()
        face = [0, 1, 2, 3, 4]
        a = _place_haptic(P.copy(), face, np.array([0.0, 0.0, 1.0]), 1.5, 0.7, 1.0)
        b = _place_haptic(P.copy(), face, np.array([0.0, 0.0, 1.0]), 1.5, 0.7, 1.0, winding=None)
        np.testing.assert_allclose(a, b)


class TestEtaRoundTripWinding(unittest.TestCase):
    """End-to-end: oin_direct reproduces the requested winding, key-exact."""

    @classmethod
    def setUpClass(cls):
        cls.gen = OIN3DGenerator(
            engine="metallogen", optimizer="ff", ff_params={"oin_direct": True}
        )

    def _assert_roundtrip(self, oin):
        oin2 = _reencode(self.gen.generate(oin))
        self.assertEqual(
            _eta_winding_multiset(oin),
            _eta_winding_multiset(oin2),
            f"winding multiset drifted: {oin2}",
        )
        self.assertEqual(
            canonical_roundtrip_key(oin),
            canonical_roundtrip_key(oin2),
            f"round-trip key mismatch: {oin2}",
        )

    def test_ferrocene_symmetric(self):
        self._assert_roundtrip(FERROCENE)

    def test_halide_face_load_bearing(self):
        self._assert_roundtrip(FE_HALIDE_FACE)

    def test_chirality_witness_ring(self):
        # A reflection instead of a proper rotation would flip the sp3 [C@H] witness
        # and change the key -- so a green key here also proves chirality is preserved.
        self._assert_roundtrip(FE_CHIRAL_WITNESS)


if __name__ == "__main__":
    unittest.main()
