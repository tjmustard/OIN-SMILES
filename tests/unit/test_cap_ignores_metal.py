"""OIN_CAP_IGNORES_METAL (v0.4.19, H4 of the serializer lane map): a ligand atom's valence cap
must not cut a ligand bond to make room for the metal contact.

`perception_core`'s per-atom cap counts every AC neighbour, the metal included, and when an atom
is over its maximum valence `remove_weakest_bond` deletes the neighbour with the largest excess
`d - r_i - r_j`. A metal contact's excess is the most negative of the set (Pd-Se 2.38 A against
covalent radii is -0.21; Se-C 1.94 A is -0.02), so the ligand bond goes: KICSUM's PhSe-CH2 is
split into two fragments and the string names a different ligand (census E1_GRAPH/LIGAND_DIFF).
With the lever a heavy non-metal atom's cap counts and cuts ligand bonds only.

Hydrogen is exempt on purpose: an H between a carbon and the metal must lose one of the two
(max valence 1), and the shipped excess rule is right there -- with H included, 8 of 296 probe
rows died on "Explicit valence for atom H, 2", one of them a verified pass. UMENEG is the
smallest of the eight. Promoted to default-ON in v0.4.19; both arms are pinned.
"""

from __future__ import annotations

import os
import sys
import unittest
from unittest import mock

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.oin.levers import _DEFAULT_ON, _HELD_OFF  # noqa: E402

FIX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
KICSUM = os.path.join(FIX, "KICSUM_comp_0.xyz")  # Pd, PhSe-CH2-pyrazole, 29 atoms
UMENEG = os.path.join(FIX, "UMENEG_comp_0.xyz")  # Ru, an H that touches C and Ru, 61 atoms


def _fragments(oin: str) -> list[str]:
    return oin.split(".")[1:]


class TestRegistry(unittest.TestCase):
    def test_held_off_with_a_reason(self):
        # Promoted in v0.4.19 (owner, 2026-09-24, conditional on the candidate sweep landing on
        # its prediction -- 87.78% / 79.12% against 87.56 / 79.02).
        self.assertIn("OIN_CAP_IGNORES_METAL", _DEFAULT_ON)
        self.assertNotIn("OIN_CAP_IGNORES_METAL", _HELD_OFF)


class TestSelenoetherIsOneLigand(unittest.TestCase):
    def test_shipped_default_splits_the_ligand(self):
        # Written as the defect, so a default flip is visible here and not only in a sweep.
        with mock.patch.dict(os.environ, {"OIN_CAP_IGNORES_METAL": "0"}):
            oin = XYZToSMILES().convert(KICSUM)
        frags = _fragments(oin)
        self.assertEqual(len(frags), 4, oin)  # PhSe and the CH2-pyrazole are two fragments
        self.assertTrue(any(f.startswith("[Se]") for f in frags), oin)
        self.assertTrue(any("[CH2]" in f for f in frags), oin)  # the orphaned methylene

    def test_lever_keeps_the_se_c_bond(self):
        with mock.patch.dict(os.environ, {"OIN_CAP_IGNORES_METAL": "1"}):
            oin = XYZToSMILES().convert(KICSUM)
        frags = _fragments(oin)
        self.assertEqual(len(frags), 3, oin)  # one bidentate N,Se ligand + two chlorides
        big = max(frags, key=len)
        self.assertIn("[Se]{", big, oin)  # selenium still a donor ...
        self.assertIn("n{", big, oin)  # ... alongside the pyrazole N, in ONE fragment
        self.assertNotIn("[CH2]", oin)  # no orphaned methylene anywhere
        self.assertEqual(sum(f.count("[Cl]") for f in frags), 2, oin)


class TestHydrogenKeepsTheShippedRule(unittest.TestCase):
    def test_encodes_identically_both_ways(self):
        with mock.patch.dict(os.environ, {"OIN_CAP_IGNORES_METAL": "1"}):
            on = XYZToSMILES().convert(UMENEG)
        with mock.patch.dict(os.environ, {"OIN_CAP_IGNORES_METAL": "0"}):
            off = XYZToSMILES().convert(UMENEG)
        self.assertEqual(on, off)
        self.assertNotIn("[H]", on)  # the H stayed on its carbon, not a hydride


if __name__ == "__main__":
    unittest.main()
