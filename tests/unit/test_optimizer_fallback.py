"""Guards for P6 (harness honesty): the g-xTB optimizer degrades loudly, once.

Three behaviors this pins:

* ``ASEOptimizer.optimize`` with no ``xtb`` binary returns the *input* molecule
  object unchanged, before any deepcopy / ASE work, and warns exactly once per
  process (not once per conformer). This is the deliberate non-fatal fallback
  (contrast MACE, which raises ``ImportError``).
* ``test_dataset_roundtrip._pass2_config`` yields ``optimizer=None`` -- so PASS 2
  builds no ``ASEOptimizer`` and spawns no subprocess -- and the honest
  ``FF_reroll_*`` tier names when ``xtb`` is absent, keeping the ``g-xTB_*`` names
  when it is present.
* ``classify_failures.classify`` is independent of ``tier_passed``, so renaming
  the tier does not reclassify any historical report.
"""

import os
import sys
import unittest
from unittest import mock

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(os.path.join(_ROOT, "src"))
sys.path.append(os.path.join(_ROOT, "tools"))

import classify_failures  # noqa: E402
from test_dataset_roundtrip import _pass2_config  # noqa: E402

import oinsmiles.generator3d.ml_optimizer as ml_optimizer  # noqa: E402
from oinsmiles.generator3d.ml_optimizer import ASEOptimizer  # noqa: E402


class TestLoudXtbFallback(unittest.TestCase):
    def setUp(self):
        # The once-per-process guard is module state; reset it so each test
        # observes the first-warning behavior independently.
        ml_optimizer._XTB_FALLBACK_WARNED = False

    @mock.patch.object(ml_optimizer.shutil, "which", return_value=None)
    def test_returns_input_mol_unchanged_before_any_work(self, _which):
        # A bare object has no .atom_list / .get_coordinate: if optimize() touched
        # it (deepcopy, ASE Atoms, subprocess) this would raise. A clean return of
        # the same object proves the early bail happened before any of that work.
        sentinel = object()
        success, energy, out = ASEOptimizer(method="g-xtb").optimize(sentinel)
        self.assertFalse(success)
        self.assertEqual(energy, 0.0)
        self.assertIs(out, sentinel)

    @mock.patch.object(ml_optimizer.shutil, "which", return_value=None)
    def test_warns_exactly_once_per_process(self, _which):
        opt = ASEOptimizer(method="g-xtb")
        with self.assertLogs(ml_optimizer.logger, level="WARNING") as cm:
            opt.optimize(object())
            opt.optimize(object())  # second call must NOT warn again
        self.assertEqual(len(cm.records), 1, cm.output)
        self.assertIn("xtb", cm.output[0].lower())

    @mock.patch.object(ml_optimizer.shutil, "which", return_value=None)
    def test_xtb_alias_also_falls_back(self, _which):
        # OIN3DGenerator's public default is optimizer="xtb"; it must hit the same
        # fallback as the explicit "g-xtb".
        sentinel = object()
        success, _energy, out = ASEOptimizer(method="xtb").optimize(sentinel)
        self.assertFalse(success)
        self.assertIs(out, sentinel)


class TestPass2Config(unittest.TestCase):
    def test_ff_reroll_when_xtb_absent(self):
        # optimizer=None => PASS 2 constructs no ASEOptimizer and spawns no
        # subprocess, while the ensemble 1->5 recovery is preserved by the caller.
        self.assertEqual(_pass2_config(False), (None, "FF_reroll_1", "FF_reroll_5"))

    def test_gxtb_when_xtb_present(self):
        self.assertEqual(_pass2_config(True), ("g-xtb", "g-xTB_1", "g-xTB_5"))


class TestTierVocabularyIndependentOfClassifier(unittest.TestCase):
    """Renaming tiers must not reclassify historical rows: classify() ignores tier_passed."""

    def _report(self, tier, **over):
        rep = {
            "molecule": "TESTMOL_comp_0",
            "status": "success",
            "tier_passed": tier,
            "error": None,
            "smiles_1": "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
            "smiles_2": "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
        }
        rep.update(over)
        return rep

    def test_success_row_same_class_across_tier_names(self):
        old = classify_failures.classify(self._report("g-xTB_5"))
        new = classify_failures.classify(self._report("FF_reroll_5"))
        self.assertEqual(old, new)
        self.assertEqual(old[0], "success")

    def test_failed_row_same_class_across_tier_names(self):
        # A soft high-RMSD failure: classification keys on the error, not the tier.
        err = "High RMSD at TESTMOL_comp_0: 1.42"
        old = classify_failures.classify(self._report("g-xTB_5", status="failed", error=err))
        new = classify_failures.classify(self._report("FF_reroll_5", status="failed", error=err))
        self.assertEqual(old, new)
        self.assertEqual(old[0], "high_rmsd")


if __name__ == "__main__":
    unittest.main()
