"""Structural stability test for BINAP-containing Pd complex.

Per MiniPRD_ChiralTests (US-001 / Negative Space): the pipeline must handle
axial-chiral biaryl ligands without crashing.  No @/@@  assertion on P atoms
is made because BINAP chirality is axial (atropisomeric), not P-centred.

Fixture: tests/fixtures/PdCl2-R-BINAP.xyz  (81 atoms, K:19,38;)
Candidate OIN artifact: tests/candidate_outputs/binap_oin.txt (2026-03-04)
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES

_BINAP_XYZ = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../fixtures/PdCl2-R-BINAP.xyz")
)


class TestBinapStability(unittest.TestCase):

    def test_binap_does_not_crash(self):
        """Pipeline returns a non-None OIN string for BINAP complex."""
        result = XYZToSMILES().convert(_BINAP_XYZ)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, str)
        self.assertGreater(len(result), 0)

    def test_binap_contains_metal_fragment(self):
        """OIN string contains a metal fragment marker."""
        result = XYZToSMILES().convert(_BINAP_XYZ)
        # Metal fragment identified by _GEO suffix (e.g. [Pd@SP1_SPL])
        self.assertIn("_SPL", result)

    def test_binap_contains_phosphorus_slots(self):
        """Both P donor atoms appear as slot-tagged fragments."""
        result = XYZToSMILES().convert(_BINAP_XYZ)
        # BINAP is bidentate — two P{N} slot markers expected
        self.assertIn("P{0}", result)
        self.assertIn("P{1}", result)


if __name__ == "__main__":
    unittest.main()
