"""Guard tests for S7 harness honesty (tools/test_dataset_roundtrip.py).

Two labeling-only additions, each pinned here (each fails against pre-S7 code):

* ``_build_run_env`` stamps ``mol_timeout`` and ``rmsd_gate`` (alongside the existing
  ``quick`` / ``optimizer_effective``) into every report, so a failed row is
  self-describing about the budget and gate it was judged under.
* ``_honesty_breakdown`` splits failed rows into FF-floor high_rmsd and quick-timeout
  artifacts vs real failures, so the backlog stops conflating them. It changes no
  status and no gate threshold -- it is purely descriptive.
"""

import os
import sys
import types
import unittest
from unittest import mock

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../tools")))

from test_dataset_roundtrip import RMSD_GATE, _build_run_env, _honesty_breakdown


def _args(quick=True, mol_timeout=30):
    return types.SimpleNamespace(quick=quick, mol_timeout=mol_timeout)


class TestRunEnvProvenance(unittest.TestCase):
    def test_stamps_budget_and_gate(self):
        env = _build_run_env(_args(quick=True, mol_timeout=30))
        # New fields (absent pre-S7):
        self.assertEqual(env["mol_timeout"], 30)
        self.assertEqual(env["rmsd_gate"], RMSD_GATE)
        # Existing fields still present:
        self.assertTrue(env["quick"])
        self.assertIn("optimizer_effective", env)

    def test_optimizer_effective_ff_without_xtb(self):
        with mock.patch("test_dataset_roundtrip.resolve_xtb_binary", return_value=None):
            env = _build_run_env(_args())
        self.assertFalse(env["xtb_available"])
        self.assertEqual(env["optimizer_effective"], "ff")

    def test_optimizer_effective_gxtb_with_xtb(self):
        with mock.patch("test_dataset_roundtrip.resolve_xtb_binary", return_value="/usr/bin/xtb"):
            env = _build_run_env(_args())
        self.assertTrue(env["xtb_available"])
        self.assertEqual(env["optimizer_effective"], "g-xtb")


class TestHonestyBreakdown(unittest.TestCase):
    def test_classifies_each_bucket(self):
        reports = [
            {"status": "success", "error": None},
            # FF-floor high_rmsd: string matched, geometry only, under FF
            {
                "status": "failed",
                "error": "High RMSD at FF_reroll_5: 1.2345",
                "optimizer_effective": "ff",
            },
            # FF-floor via explicit flag (optimizer field absent)
            {"status": "failed", "error": "High RMSD at UFF_1: 1.10", "ff_floor": True},
            # quick timeout: budget artifact
            {
                "status": "failed",
                "error": "TimeoutException at UFF_1: exceeded 30s ...",
                "quick": True,
            },
            # full-budget timeout: kept separate from quick
            {
                "status": "failed",
                "error": "TimeoutException at UFF_1: exceeded 1800s ...",
                "quick": False,
            },
            # real failure: string mismatch is a genuine accuracy defect
            {"status": "failed", "error": "String mismatch at FF_reroll_5. Exp: ..., Got: ..."},
        ]
        bd = _honesty_breakdown(reports)
        self.assertEqual(bd["ff_floor_high_rmsd"], 2)
        self.assertEqual(bd["quick_timeout"], 1)
        self.assertEqual(bd["timeout_full_budget"], 1)
        self.assertEqual(bd["real_failure"], 1)

    def test_high_rmsd_under_gxtb_is_not_ff_floor(self):
        # A high_rmsd produced with a real optimizer is not an FF-floor artifact.
        reports = [
            {
                "status": "failed",
                "error": "High RMSD at g-xTB_5: 1.5",
                "optimizer_effective": "g-xtb",
            }
        ]
        bd = _honesty_breakdown(reports)
        self.assertEqual(bd["ff_floor_high_rmsd"], 0)
        self.assertEqual(bd["real_failure"], 1)


if __name__ == "__main__":
    unittest.main()
