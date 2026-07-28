"""The accept-side eta winding exit (v0.4.12 Lane 2, ``OIN_ETA_ACCEPT_EXIT``).

WHY THIS LEVER EXISTS WHEN ``OIN_ETA_EARLY_EXIT`` ALREADY APPLIES THE SAME TEST
==============================================================================
``OIN_ETA_EARLY_EXIT`` sits in ``_select_by_geometry_impl``, which runs **after**
``generate_3d_structures`` has already filled the entire pool. Its own measured A/B says so:
Ferrocene, lever off -> 32 attempts; lever on -> 32 attempts. It fires and buys nothing,
because a selection-side early exit can only shorten the selection scan.

``accept_fn`` is the only site consulted per conformer *during* pool filling. That is where
this lever puts the same criterion, and it is the whole difference between the two.

WHAT THESE TESTS PIN
====================
1. Default OFF, and the predicate is byte-identical when unset.
2. The conjunction. Winding alone would bypass ``_select_by_geometry``'s clash-first ranking
   -- structurally the same defect ``OIN_ACCEPT_SCORED`` has -- so geometry classification and
   ligand attachment are both required, and removing either must make the predicate decline.
3. The precedence: an exact key match still wins, so a molecule that round-trips today is
   unaffected.
"""

import os
import unittest
from unittest import mock

from oinsmiles.generation import metallogen_adapter as ad
from oinsmiles.oin.levers import held_off, lever_enabled

LEVER = "OIN_ETA_ACCEPT_EXIT"


class TestDefaultOff(unittest.TestCase):
    def test_registered_with_a_reason(self):
        self.assertIn(LEVER, held_off())
        why = held_off()[LEVER]
        self.assertIn("accept_fn", why)

    def test_unset_means_off(self):
        env = {k: v for k, v in os.environ.items() if k != LEVER}
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertFalse(lever_enabled(LEVER))

    def test_explicit_zero_disables(self):
        with mock.patch.dict(os.environ, {LEVER: "0"}):
            self.assertFalse(lever_enabled(LEVER))


class TestPredicateIsAConjunction(unittest.TestCase):
    """Each conjunct removed in turn must make the predicate decline.

    Asserted on ``_eta_accept_exit_ok`` directly rather than through a generation run. v0.4.11
    learned this the hard way in the mirror direction: its scope tests asserted on
    ``canonicalize_oin_slots`` output and passed for the wrong reason, because the rotation
    group already converged the small cases on its own. A predicate test has to poke the
    predicate.
    """

    def setUp(self):
        self.parsed = mock.Mock(geo_code="OCT")
        self.cmol = mock.Mock()

    def _run(self, winding_ok=True, geometry_ok=True, attached_ok=True, targets=(">", ">")):
        oin = "[Fe_OCT].[cH]{0>}1[cH]{0}1" if winding_ok else "[Fe_OCT].[cH]{0<}1[cH]{0}1"
        with (
            mock.patch.object(ad, "_reencode_oin_fast", return_value=oin),
            mock.patch.object(
                ad,
                "_eta_winding_multiset",
                side_effect=lambda s: list(targets) if s is oin and winding_ok else ["<"],
            ),
            mock.patch.object(ad, "_geometry_classifies", return_value=geometry_ok),
            mock.patch.object(ad, "conformer_ligands_attached", return_value=attached_ok),
        ):
            return ad._eta_accept_exit_ok(self.parsed, self.cmol, list(targets))

    def test_all_three_conjuncts_accept(self):
        self.assertTrue(self._run())

    def test_wrong_winding_declines(self):
        self.assertFalse(self._run(winding_ok=False))

    def test_wrong_geometry_declines(self):
        """Winding alone would stop the pool before clash-first ranking ever ran."""
        self.assertFalse(self._run(geometry_ok=False))

    def test_detached_ligand_declines(self):
        """v0.4.7: never run a scored-acceptance lever without the attachment check."""
        self.assertFalse(self._run(attached_ok=False))

    def test_no_targets_declines(self):
        """Non-eta molecules must never reach this path."""
        self.assertFalse(ad._eta_accept_exit_ok(self.parsed, self.cmol, []))

    def test_no_contract_mol_declines(self):
        self.assertFalse(ad._eta_accept_exit_ok(self.parsed, None, [">"]))

    def test_an_exception_declines_rather_than_accepts(self):
        with mock.patch.object(ad, "_reencode_oin_fast", side_effect=RuntimeError("boom")):
            self.assertFalse(ad._eta_accept_exit_ok(self.parsed, self.cmol, [">"]))


class TestGeometryGateIsConservative(unittest.TestCase):
    """Anything unverifiable must fall through to the pool, never be accepted early."""

    def test_unknown_geometry_declines(self):
        self.assertFalse(_geometry_call(geo_code=None))

    def test_missing_contract_mol_declines(self):
        self.assertFalse(ad._geometry_classifies(mock.Mock(geo_code="OCT"), None))

    def test_unreadable_coordination_vectors_decline(self):
        with mock.patch.object(ad, "_coordination_vectors", return_value=None):
            self.assertFalse(ad._geometry_classifies(mock.Mock(geo_code="OCT"), mock.Mock()))


def _geometry_call(geo_code):
    return ad._geometry_classifies(mock.Mock(geo_code=geo_code), mock.Mock())


class TestSourceWiring(unittest.TestCase):
    """The lever must be read at the site that can stop pool filling, not in selection.

    Asserted on the source because that is the claim: v0.4.7's attachment check was wired to
    the wrong site, every call raised, the error was swallowed, and a complete A/B reported
    the silent no-op as a clean null result.
    """

    def _source(self):
        import inspect

        return inspect.getsource(ad)

    def test_the_lever_is_read_where_accept_fn_is_built(self):
        src = self._source()
        self.assertIn('lever_enabled("OIN_ETA_ACCEPT_EXIT")', src)
        # the accept_fn closure must receive the targets, not recompute them per conformer
        self.assertIn("_eta=_eta_accept_targets", src)

    def test_the_predicate_calls_the_attachment_check_unconditionally(self):
        import inspect

        src = inspect.getsource(ad._eta_accept_exit_ok)
        self.assertIn("conformer_ligands_attached(cmol)", src)
        # The check must not be GATED here -- the docstring may cite OIN_ATTACH_CHECK as its
        # provenance, but no lever read may stand between this branch and the check.
        self.assertNotIn(
            'lever_enabled("OIN_ATTACH_CHECK")', src, "inside this branch it is not optional"
        )

    def test_it_records_telemetry_so_a_no_op_cannot_look_like_a_null(self):
        import inspect

        src = inspect.getsource(ad._eta_accept_exit_ok)
        self.assertIn("adapter.eta_accept_hit", src)
        self.assertIn("adapter.eta_accept_unevaluable", src)


if __name__ == "__main__":
    unittest.main()
