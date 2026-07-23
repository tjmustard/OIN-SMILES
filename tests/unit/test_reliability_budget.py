"""Guards for v0.4.4 SL4 (reliability): the two levers that shrink the hard-fail bucket.

* RMSD is demoted from a pass/fail gate to a reported diagnostic
  (``tools/test_dataset_roundtrip.py::_attempt_generation``). For a lossless-hash the
  contract is the canonical-key string match; coordination-sphere RMSD is only
  ~0.22-correlated with geometric quality, so a string-exact round-trip that only
  exceeds the tightness threshold is now a success carrying an ``rmsd_over_gate``
  diagnostic, not a failure. A whole-complex ``clash_count`` is recorded alongside it.

* A gated, OFF-by-default no-acceptance-progress cutoff
  (``generator3d/__init__.py::generate_3d_structures``) stops the attempt loop once the
  vdW acceptance gate has rejected every embed for N consecutive attempts (the OSIHUU
  pattern), so a stuck molecule fails fast instead of burning the full ``max_attempts``.
  Unset, the loop count is byte-identical to pristine.

Each test fails against the pre-SL4 code (RMSD >= gate returned False; the cutoff and the
clash diagnostic did not exist).
"""

import os
import sys
import types
import unittest
from unittest import mock

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../tools")))

from oinsmiles.generation.metallogen_adapter import convert_parsed_to_msmiles
from oinsmiles.generation.oin_parser import OINParser
from oinsmiles.generator3d import embed, generate_3d_structures

# A trivial, always-parseable complex (as in test_embed_budget).
CISPLATIN_OIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"


def _msmiles(oin):
    return convert_parsed_to_msmiles(OINParser().parse(oin))


class TestRmsdDemotedToDiagnostic(unittest.TestCase):
    """A string-exact round-trip must not fail on RMSD alone; clash is recorded."""

    def _run_attempt(self, rmsd_value, clash=(2, 0, 0.5)):
        import test_dataset_roundtrip as harness

        gen = mock.MagicMock()
        gen.generate.return_value = types.SimpleNamespace(xyz="1\n\nPt 0.0 0.0 0.0\n", mol=None)
        xyz_conv = mock.MagicMock()
        xyz_conv.convert.return_value = "OIN"
        fake_mol = mock.MagicMock()
        fake_mol.GetAtoms.return_value = []
        report = {"metrics": {}, "smiles_2": None, "error": None, "status": "pending"}
        detail = (rmsd_value, None if rmsd_value is not None else "no_metal")
        with (
            mock.patch.object(harness, "XYZToSMILES", return_value=xyz_conv),
            mock.patch.object(harness, "canonical_roundtrip_key", side_effect=lambda s: "K"),
            mock.patch.object(harness, "normalize_oin_for_comparison", side_effect=lambda s: s),
            mock.patch.object(harness.Chem, "MolFromXYZFile", return_value=fake_mol),
            mock.patch.object(harness, "vdw_clash_count", return_value=clash),
            mock.patch.object(harness, "calculate_tmc_rmsd_detailed", return_value=detail),
            mock.patch.object(harness, "read_atom_count", return_value=5),
            mock.patch.dict(harness.RUN_ENV, {"optimizer_effective": "ff"}, clear=True),
        ):
            ok, _ = harness._attempt_generation("UFF_1", gen, "OIN", "/nonexistent/in.xyz", report)
        return ok, report

    def test_high_rmsd_string_exact_is_success(self):
        ok, report = self._run_attempt(1.5)  # >= RMSD_GATE (1.0)
        self.assertTrue(ok, "a string-exact round-trip must not fail on RMSD alone")
        self.assertEqual(report["status"], "success")
        self.assertEqual(report["metrics"]["rmsd"], 1.5)
        self.assertTrue(report["metrics"]["rmsd_over_gate"])
        self.assertTrue(report.get("ff_floor"), "FF-floor provenance flag is still set")

    def test_low_rmsd_carries_no_over_gate_flag(self):
        ok, report = self._run_attempt(0.4)
        self.assertTrue(ok)
        self.assertEqual(report["status"], "success")
        self.assertNotIn("rmsd_over_gate", report["metrics"])

    def test_clash_count_recorded(self):
        _ok, report = self._run_attempt(0.4, clash=(3, 1, 0.4))
        self.assertEqual(report["metrics"]["clash_count"], 3)

    def test_rmsd_mapping_failure_still_fails(self):
        # rmsd is None (the metric could not map the sphere) stays a failure -- a
        # metric error, not a demonstrated geometry pass (approved SL4 scope).
        ok, report = self._run_attempt(None)
        self.assertFalse(ok)
        self.assertNotEqual(report["status"], "success")
        self.assertIn("RMSD mapping failed", report["error"])


class TestNoProgressCutoff(unittest.TestCase):
    """The gated no-acceptance-progress cutoff; unset, the loop count is unchanged."""

    def setUp(self):
        self.msmiles = _msmiles(CISPLATIN_OIN)

    def _count_calls_returning_none(self, ff_params=None, env=None):
        """Run generation with get_embedding stubbed to return None; count the calls.

        Returning None (not raising) marks had_nonstructural_embed True and files
        nothing, i.e. the pool never grows -- exactly the gate-rejection pattern the
        cutoff targets, without needing a real embed.
        """
        calls = {"n": 0}

        def none_embed(*a, **k):
            calls["n"] += 1
            return None

        with mock.patch.dict(os.environ):
            os.environ.pop("OIN_EMBED_NO_PROGRESS", None)
            if env:
                os.environ.update(env)
            with mock.patch.object(embed, "get_embedding", none_embed):
                res = generate_3d_structures(self.msmiles, ff_params=ff_params)
        return calls["n"], res

    def test_cutoff_stops_early_via_ff_params(self):
        n, res = self._count_calls_returning_none(
            ff_params={"max_attempts": 200, "embed_no_progress_attempts": 5}
        )
        self.assertEqual(res, [], "no embed validated, so the pool must be empty")
        self.assertEqual(n, 5, "the cutoff must stop after 5 no-progress attempts")

    def test_cutoff_stops_early_via_env(self):
        n, res = self._count_calls_returning_none(
            ff_params={"max_attempts": 200}, env={"OIN_EMBED_NO_PROGRESS": "5"}
        )
        self.assertEqual(res, [])
        self.assertEqual(n, 5)

    def test_unset_runs_full_attempt_budget(self):
        # Default (knob unset): the added bookkeeping must not cut the loop short.
        n, res = self._count_calls_returning_none(ff_params={"max_attempts": 12})
        self.assertEqual(res, [])
        self.assertEqual(n, 12, "an unset cutoff must run every attempt (byte-identical)")

    def test_never_firing_cutoff_matches_default_pool(self):
        # A cutoff so high it never fires must produce the identical conformer pool as
        # the default path -- the byte-identity guard for the success path (seed=42).
        def coords(mols):
            return [
                tuple(tuple(round(c, 6) for c in a.get_coordinate()) for a in m.atom_list)
                for m in mols
            ]

        with mock.patch.dict(os.environ):
            os.environ.pop("OIN_EMBED_NO_PROGRESS", None)
            default_pool = generate_3d_structures(self.msmiles, uff_pool_size=2)
            gated_pool = generate_3d_structures(
                self.msmiles, uff_pool_size=2, ff_params={"embed_no_progress_attempts": 100000}
            )
        self.assertTrue(default_pool, "cisplatin must generate at least one conformer")
        self.assertEqual(coords(default_pool), coords(gated_pool))


if __name__ == "__main__":
    unittest.main()
