"""OIN_RC1_PROPAGATE (v0.4.19 serializer lane): the aligner's eta rank swap must reach the writer.

`OINDiscreteAligner` step 3b (RC1) re-ranks "same-mass" eta fragments by content -- and same-mass
means the same FIRST BINDING ATOM, so a Cp and an allyl both bound through carbon qualify. Each
w-tag entry keeps the local indices of the fragment it was computed on, but `get_oin_string` reads
the entry's rank against the un-permuted fragment list, so the Cp receives the allyl's three
markers and the allyl three of the Cp's five (the surplus is dropped silently at inline.py).
The string then carries a different coordination sphere from the input (census E1_GRAPH /
SPHERE_DIFF), and the generator builds it faithfully -- a byte-exact FALSE pass when it does.

TULTAX_comp_0 (Ru, Cp + allyl + picolinate, 32 atoms) is the smallest of the 14 release-sweep
strings the lever repairs. The lever is measured and held OFF; both arms are pinned here so a
default flip has to change this file.
"""

from __future__ import annotations

import os
import re
import sys
import unittest
from unittest import mock

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.oin.levers import _DEFAULT_ON, _HELD_OFF  # noqa: E402

TULTAX = os.path.join(os.path.dirname(__file__), "..", "fixtures", "TULTAX_comp_0.xyz")
_MARK = re.compile(r"\{(\d+)[<>^]?\}")


def _markers_per_fragment(oin: str) -> dict[str, int]:
    """{fragment SMILES: number of slot markers}, for every non-metal fragment."""
    out = {}
    for frag in oin.split(".")[1:]:
        out[frag] = len(_MARK.findall(frag))
    return out


def _eta_fragments(oin: str) -> list[tuple[str, int]]:
    """(fragment, marker count) for the fragments whose markers share ONE slot (haptic groups)."""
    eta = []
    for frag, n in _markers_per_fragment(oin).items():
        slots = set(_MARK.findall(frag))
        if n > 1 and len(slots) == 1:
            eta.append((frag, n))
    return eta


class TestRegistry(unittest.TestCase):
    def test_held_off_with_a_reason(self):
        self.assertIn("OIN_RC1_PROPAGATE", _HELD_OFF)
        self.assertNotIn("OIN_RC1_PROPAGATE", _DEFAULT_ON)
        self.assertIn("E1_GRAPH/SPHERE_DIFF", _HELD_OFF["OIN_RC1_PROPAGATE"])


class TestCpPlusAllyl(unittest.TestCase):
    def test_shipped_default_swaps_the_markers(self):
        # Written as the defect, so a default flip is visible here and not only in a sweep.
        with mock.patch.dict(os.environ, {"OIN_RC1_PROPAGATE": "0"}):
            oin = XYZToSMILES().convert(TULTAX)
        eta = dict(_eta_fragments(oin))
        cp = [f for f in eta if "1[cH]" in f or "[cH]" in f and f.endswith("1")]
        self.assertTrue(cp, oin)
        # The Cp ring carries only the allyl's three markers ...
        self.assertEqual(eta[cp[0]], 3, oin)
        # ... and two of its five ring carbons are written bare (unmarked, no H bracket).
        self.assertIn("cc1", cp[0], oin)

    def test_lever_gives_every_eta_carbon_its_marker(self):
        with mock.patch.dict(os.environ, {"OIN_RC1_PROPAGATE": "1"}):
            oin = XYZToSMILES().convert(TULTAX)
        eta = sorted(_eta_fragments(oin), key=lambda e: e[1])
        self.assertEqual([n for _f, n in eta], [3, 5], oin)
        allyl, cp = eta[0][0], eta[1][0]
        self.assertEqual(len(re.findall(r"\[CH2?\]", allyl)), 3, oin)  # [CH2]{n>}[CH]{n}=[CH2]{n}
        self.assertEqual(cp.count("[cH]"), 5, oin)  # all five ring carbons marked, none bare
        # Everything else is untouched by the permutation: same fragments, same picolinate.
        with mock.patch.dict(os.environ, {"OIN_RC1_PROPAGATE": "0"}):
            shipped = XYZToSMILES().convert(TULTAX)
        pic = [f for f in oin.split(".") if "C(=O)" in f]
        self.assertEqual(pic, [f for f in shipped.split(".") if "C(=O)" in f])

    def test_lever_is_a_no_op_without_two_eta_groups_of_one_element(self):
        # Ferrocene: two IDENTICAL Cp rings -- RC1 may swap them, the strings cannot differ.
        ferrocene = os.path.join(os.path.dirname(__file__), "..", "fixtures", "Ferrocene.xyz")
        if not os.path.exists(ferrocene):
            self.skipTest("no ferrocene fixture")
        with mock.patch.dict(os.environ, {"OIN_RC1_PROPAGATE": "1"}):
            on = XYZToSMILES().convert(ferrocene)
        with mock.patch.dict(os.environ, {"OIN_RC1_PROPAGATE": "0"}):
            off = XYZToSMILES().convert(ferrocene)
        self.assertEqual(on, off)


if __name__ == "__main__":
    unittest.main()
