"""Tests for string-exact acceptance (v0.4.15 Lane 2, ``OIN_ACCEPT_STRING_EXACT``).

The defect: ``compare._parse_vertex_colors`` folds reflection deliberately -- it colours every
donor of a ligand with that ligand's whole canonical body, so a transposition of two
same-coloured donors (an ODD permutation, i.e. a reflection) is invisible to the vertex
signature. ``accept_fn`` decides by that key, so the generator can return the MIRROR IMAGE and
the harness files it as a benign string difference.

**A lossy key must never be reused as an acceptance predicate for an axis it folds.** Third
instance in this project, after v0.4.8 and v0.4.11.

These tests pin three things:

1. The key really does fold what the normalized string keeps -- the premise, not an assumption.
2. Lever off is byte-identical: the predicate is never consulted and ``accept_fn`` keeps
   returning plain booleans.
3. The ``ACCEPT_INCUMBENT`` contract, which is what makes the lever cost latency instead of
   accuracy: a bool-returning ``accept_fn`` is unaffected, and the sentinel is distinguishable
   from ``True`` by identity, never by truthiness.
"""

import os
import unittest

from oinsmiles.generation.metallogen_adapter import _string_exact_match
from oinsmiles.generator3d import ACCEPT_INCUMBENT
from oinsmiles.oin.compare import canonical_roundtrip_key, normalize_oin_for_comparison
from oinsmiles.oin.levers import held_off, lever_enabled

# A real enantiomer pair from the v0.4.14 baseline sweep: AGAVIQ_comp_0's input and the structure
# the round trip built. The two P donors are transposed between slots 3 and 4 -- nothing else
# differs. `veto_residue_chirality.py` classifies this molecule MIRROR_MATCH.
AGAVIQ_IN = "[Co_SPY].c{0}1ccccc1.CC(C)P{3}(Oc1c{1}c(OP{4}(C(C)C)C(C)C)ccc1)C(C)C.S{2}c1ccccc1"
AGAVIQ_GEN = "[Co_SPY].c{0}1ccccc1.CC(C)P{4}(Oc1c{1}c(OP{3}(C(C)C)C(C)C)ccc1)C(C)C.S{2}c1ccccc1"


class _Parsed:
    """The one attribute the predicate reads off ``ParsedOIN``."""

    def __init__(self, oin):
        self.original_oin = oin


class _Lever:
    def __init__(self, **levers):
        self.levers = levers
        self.prior = {}

    def __enter__(self):
        for k, v in self.levers.items():
            self.prior[k] = os.environ.get(k)
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        return self

    def __exit__(self, *exc):
        for k, v in self.prior.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        return False


class TestThePremise(unittest.TestCase):
    """🔴 The whole lane rests on this: the key folds what the string keeps.

    If these fail, the predicate cannot separate an enantiomer pair and the lane is pointless --
    so they are asserted on a real measured pair rather than a constructed one.
    """

    def test_the_two_strings_really_are_different(self):
        self.assertNotEqual(AGAVIQ_IN, AGAVIQ_GEN)

    def test_the_round_trip_key_CANNOT_see_the_difference(self):
        self.assertEqual(
            canonical_roundtrip_key(AGAVIQ_IN),
            canonical_roundtrip_key(AGAVIQ_GEN),
            "premise: the key folds reflection, which is why acceptance accepted the mirror",
        )

    def test_the_normalized_string_CAN_see_the_difference(self):
        self.assertNotEqual(
            normalize_oin_for_comparison(AGAVIQ_IN),
            normalize_oin_for_comparison(AGAVIQ_GEN),
            "the normalization keeps absolute slots, so the donor transposition survives it",
        )

    def test_the_difference_is_a_transposition_of_two_same_element_donors(self):
        """Names the mechanism: P{3}<->P{4}, an odd permutation of two same-coloured vertices."""
        self.assertEqual(
            AGAVIQ_IN.replace("P{3}", "P{X}").replace("P{4}", "P{3}").replace("P{X}", "P{4}"),
            AGAVIQ_GEN,
        )


class TestPredicate(unittest.TestCase):
    def test_the_same_string_matches(self):
        self.assertTrue(_string_exact_match(_Parsed(AGAVIQ_IN), AGAVIQ_IN))

    def test_the_mirror_does_not_match(self):
        self.assertFalse(
            _string_exact_match(_Parsed(AGAVIQ_IN), AGAVIQ_GEN),
            "this is the 201-molecule class the lane exists for",
        )

    def test_a_slot_renumbering_also_does_not_match(self):
        """Deliberately in scope: slot drift fails ``byte_exact`` too, so it is the same defect."""
        self.assertFalse(
            _string_exact_match(
                _Parsed("[Pt_SPL].N{0}.[Cl]{1}"),
                "[Pt_SPL].N{1}.[Cl]{0}",
            )
        )

    def test_bound_water_notation_is_still_folded(self):
        """``[OH2]``/``O`` equivalence is a benign encoder difference the normalization owns."""
        self.assertTrue(
            _string_exact_match(_Parsed("[Pt_SPL].[OH2]{0}.[Cl]{1}"), "[Pt_SPL].O{0}.[Cl]{1}")
        )

    def test_no_target_abstains_rather_than_rejecting(self):
        self.assertTrue(_string_exact_match(_Parsed(None), AGAVIQ_GEN))

    def test_no_reencode_abstains_rather_than_rejecting(self):
        self.assertTrue(_string_exact_match(_Parsed(AGAVIQ_IN), None))

    def test_a_raising_target_abstains_rather_than_rejecting(self):
        """Errors stay PERMISSIVE: a conformer is demoted on evidence, never on ignorance."""

        class Boom:
            @property
            def original_oin(self):
                raise RuntimeError("perception exploded")

        self.assertTrue(_string_exact_match(Boom(), AGAVIQ_GEN))


class TestIncumbentSentinel(unittest.TestCase):
    """The contract that keeps a stricter predicate from costing accuracy."""

    def test_the_sentinel_is_not_true_and_not_false(self):
        self.assertIsNot(ACCEPT_INCUMBENT, True)
        self.assertIsNot(ACCEPT_INCUMBENT, False)

    def test_the_sentinel_is_truthy_which_is_exactly_why_identity_checks_are_required(self):
        """Pins the trap: every consumer must test ``is ACCEPT_INCUMBENT``, not truthiness.

        A truthiness read would accept the very conformer the lever exists to step over -- the
        same shape of defect as reading ``os.environ.get(lever)`` as a boolean.
        """
        self.assertTrue(bool(ACCEPT_INCUMBENT))

    def test_a_plain_bool_accept_fn_is_unaffected_by_the_sentinel_protocol(self):
        for verdict in (True, False):
            self.assertIsNot(verdict, ACCEPT_INCUMBENT)


class TestLeverRegistration(unittest.TestCase):
    def test_lever_is_registered_as_held_off_with_a_reason(self):
        self.assertIn("OIN_ACCEPT_STRING_EXACT", held_off())

    def test_the_rationale_states_the_measurement_and_the_refuted_alternative(self):
        why = held_off()["OIN_ACCEPT_STRING_EXACT"]
        self.assertIn("183", why, "the measured 183/183 separation belongs in the justification")
        self.assertIn("metal_config", why, "why the chartered approach was NOT used belongs there")
        self.assertIn(
            "ACCEPT_INCUMBENT", why, "the non-regression mechanism must be stated, not implied"
        )

    def test_lever_defaults_off(self):
        with _Lever(OIN_ACCEPT_STRING_EXACT=None):
            self.assertFalse(lever_enabled("OIN_ACCEPT_STRING_EXACT"))

    def test_a_zero_string_disables_rather_than_enables(self):
        with _Lever(OIN_ACCEPT_STRING_EXACT="0"):
            self.assertFalse(lever_enabled("OIN_ACCEPT_STRING_EXACT"))


class TestWiring(unittest.TestCase):
    """The predicate must be READ through ``lever_enabled`` at the acceptance seam."""

    def test_the_acceptance_seam_gates_on_the_lever(self):
        import inspect

        from oinsmiles.generation import metallogen_adapter as MA

        src = inspect.getsource(MA._reencode_key_matches)
        self.assertIn('lever_enabled("OIN_ACCEPT_STRING_EXACT")', src)
        self.assertIn("ACCEPT_INCUMBENT", src)

    def test_the_early_exit_rescan_uses_an_identity_check_not_truthiness(self):
        """Regression guard for the trap in ``TestIncumbentSentinel`` above."""
        import inspect

        from oinsmiles.generation import metallogen_adapter as MA

        src = inspect.getsource(MA._select_by_geometry_impl)
        self.assertIn("_verdict is True", src, "truthiness here would accept the mirror")


if __name__ == "__main__":
    unittest.main()
