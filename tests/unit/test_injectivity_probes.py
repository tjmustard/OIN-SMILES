"""Injectivity-audit guards (Y1): the encoder's confirmed chirality blind spots.

These lock in the Wave-1 findings from ``tools/injectivity`` and the independent oracle.
Two flavours, following the repo convention (see
``tests/integration/test_isomer_divergence.py::test_metal_stereo_raw_only``):

* **assert-current-behavior** -- document the blindness as a fact so the suite stays green
  and a regression that *changes* it is caught.
* **aspirational ``@expectedFailure``** -- encode the *desired* divergence. It fails today
  (kept green as an expected failure); when a future encoder fix lands it flips to an
  unexpected success and turns the suite red, flagging that the guard should be promoted.

Rigid, curated fixtures only -- the geometric oracle is valid for rigid species.
"""

import sys
import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from tools.injectivity.oracle import is_distinct_enantiomer  # noqa: E402
from tools.injectivity.twin_collision import (  # noqa: E402
    VERDICT_ENCODER_BLIND,
    VERDICT_KEY_BLIND,
    probe_mirror,
)

FIX = _ROOT / "tests" / "fixtures"
CISPLATIN = FIX / "CisPlatin.xyz"
FAC_IRPPY3 = FIX / "fac-Ir(ppy)3.xyz"
BINAP = FIX / "PdCl2-R-BINAP.xyz"
POJJOP = FIX / "POJJOP.xyz"  # Pd complex; sole stereocentre is a metal-bound 2° amine (P3)


class TestIndependentOracle(unittest.TestCase):
    """The geometric enantiomer oracle, independent of the OIN encoder."""

    def test_achiral_control_not_distinct(self):
        v = is_distinct_enantiomer(CISPLATIN)
        self.assertFalse(v.distinct, f"cisplatin mirror should superimpose (rmsd={v.rmsd})")

    def test_metal_chirality_is_distinct(self):
        v = is_distinct_enantiomer(FAC_IRPPY3)
        self.assertTrue(v.distinct, "fac-Ir(ppy)3 Δ/Λ mirror must be a distinct isomer")

    def test_axial_chirality_is_distinct(self):
        v = is_distinct_enantiomer(BINAP)
        self.assertTrue(v.distinct, "R/S-BINAP mirror must be a distinct isomer")

    def test_metal_bound_amine_is_distinct(self):
        v = is_distinct_enantiomer(POJJOP)
        self.assertTrue(v.distinct, "POJJOP metal-bound-amine mirror must be a distinct isomer")


class TestConfirmedBlindSpots(unittest.TestCase):
    """Assert-current-behavior: the two Wave-1 collisions, documented as facts."""

    def test_metal_delta_lambda_is_key_blind(self):
        """fac-Ir(ppy)3 Δ/Λ enantiomers collapse at the round-trip KEY (batch false-positive)."""
        o = probe_mirror(FAC_IRPPY3)
        self.assertTrue(o.oracle_distinct)
        self.assertTrue(o.key_equal, "metal Δ/Λ key regression: enantiomers now key-diverge")
        self.assertEqual(o.verdict, VERDICT_KEY_BLIND)

    def test_axial_is_encoder_blind(self):
        """R-BINAP and S-BINAP encode to byte-identical OIN strings (total blindness)."""
        o = probe_mirror(BINAP)
        self.assertTrue(o.oracle_distinct)
        self.assertTrue(o.raw_equal, "axial regression: BINAP enantiomers now raw-diverge")
        self.assertEqual(o.verdict, VERDICT_ENCODER_BLIND)

    def test_metal_bound_amine_is_encoder_blind(self):
        """POJJOP's sole stereocentre is a metal-bound 2° amine; its two enantiomers
        encode to byte-identical OIN strings (the Zone-A total_degree<4 clear -- P3)."""
        o = probe_mirror(POJJOP)
        self.assertTrue(o.oracle_distinct)
        self.assertTrue(o.raw_equal, "P3 regression: POJJOP enantiomers now raw-diverge")
        self.assertEqual(o.verdict, VERDICT_ENCODER_BLIND)

    def test_achiral_control_invariant(self):
        o = probe_mirror(CISPLATIN)
        self.assertFalse(o.oracle_distinct)
        self.assertTrue(o.raw_equal)


class TestAspirationalDivergence(unittest.TestCase):
    """The behavior a lossless encoder SHOULD have. Fails today -> expected failures."""

    @unittest.expectedFailure
    def test_metal_chirality_should_diverge_at_key(self):
        o = probe_mirror(FAC_IRPPY3)
        self.assertFalse(o.key_equal, "metal Δ/Λ should diverge at the round-trip key")

    @unittest.expectedFailure
    def test_axial_should_diverge_in_raw_string(self):
        o = probe_mirror(BINAP)
        self.assertFalse(o.raw_equal, "R/S-BINAP should not encode to identical strings")

    @unittest.expectedFailure
    def test_metal_bound_amine_should_diverge_in_raw_string(self):
        """Still an expected failure, and that is now a statement about the DEFAULT.

        The v0.4.5 Lane 6 plan expected this xfail to flip to passing. It does not, and
        should not: the fix is gated behind ``OIN_EMIT_LOCKED_DONOR`` (default OFF, like
        every other v0.4.5 lever), and "levers-OFF output is byte-identical" is a harder
        acceptance requirement than flipping this marker. The capability itself is guarded
        permanently in ``tests/unit/test_locked_donor.py``, which asserts divergence with
        the lever on plus the full three-property test. Flip this marker if and when the
        lever is promoted to default-ON.
        """
        o = probe_mirror(POJJOP)
        self.assertFalse(o.raw_equal, "metal-bound 2° amine enantiomers should not encode alike")


if __name__ == "__main__":
    unittest.main()
