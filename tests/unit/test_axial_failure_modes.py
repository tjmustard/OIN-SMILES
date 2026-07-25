"""Guards against the axial path failing SILENTLY (Y2 P2 hardening).

Three defects were found while wiring the axial fix. Each is repaired, but the property worth
locking in is not the individual bug -- it is that this path must never again produce a
plausible wrong answer without saying so. A wrong atropisomer still satisfies the round-trip
key (the key deliberately folds the axial token), so nothing downstream catches it: exactly
the "self-consistent but wrong" failure the Y1 audit exists to hunt.

1. **Trailing metadata broke the generation parser.** ``oin_parser`` decided inline-vs-sidecar
   by the ABSENCE of ``|``, so a trailing `` |ax:-|`` misrouted an inline OIN to the sidecar
   branch and lost the geometry code. (This one at least raised.)
2. **A wrong atropisomer was returned silently.** The axial branches fall through on a miss,
   so a mismatch must be reported at the single exit.
3. **Perception failed on every candidate.** Unsanitized conformers raise in
   ``CanonicalRankAtoms``; the ``try/except`` turned that into ``None`` and the filter
   compared ``None`` against ``'-'``. A perception failure must be distinguishable from a
   genuine miss.
"""

import logging
import sys
import unittest
import unittest.mock
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from rdkit import Chem  # noqa: E402

from oinsmiles.generation.oin_parser import OINParser  # noqa: E402
from oinsmiles.oin.axial import mol_axial_token, parse_axial_token  # noqa: E402

PLAIN_OIN = "[Pd_SPL].c1ccccc1P{0}.[Cl]{1}.[Cl]{2}"


class TestParserToleratesTrailingMetadata(unittest.TestCase):
    """Bug 1: the axial suffix must be transparent to the generation parser."""

    def test_token_does_not_change_the_parse(self):
        plain = OINParser().parse(PLAIN_OIN)
        tokened = OINParser().parse(PLAIN_OIN + " |ax:-|")
        self.assertEqual(tokened.geo_code, plain.geo_code)
        self.assertEqual(tokened.smiles, plain.smiles)
        self.assertEqual(tokened.fragments, plain.fragments)

    def test_geometry_code_survives(self):
        # the original failure signature was geo_code == "" -> "Geometry code '' not supported"
        self.assertEqual(OINParser().parse(PLAIN_OIN + " |ax:+|").geo_code, "SPL")

    def test_original_oin_keeps_the_token_for_the_adapter(self):
        parsed = OINParser().parse(PLAIN_OIN + " |ax:-|")
        self.assertEqual(parse_axial_token(parsed.original_oin), "-")


class TestPerceptionFailureIsRecoverableAndVisible(unittest.TestCase):
    """Bug 3: unsanitized mols must not silently read as 'no axial token'."""

    @staticmethod
    def _unsanitized_biaryl():
        """A hindered biaryl mol with 3D coords but NO sanitization -- as pool mols arrive."""
        from rdkit.Chem import AllChem

        m = Chem.AddHs(Chem.MolFromSmiles("COc1cccc(C)c1-c1c(C)cccc1[N+](=O)[O-]"))
        AllChem.EmbedMolecule(m, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(m)
        raw = Chem.Mol(m)
        raw.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(raw)  # ring info without full sanitization
        return raw

    def test_sanitize_fallback_recovers_the_token(self):
        # the fix: mol_axial_token retries on a sanitized copy rather than giving up
        token = mol_axial_token(self._unsanitized_biaryl())
        self.assertIn(token, ("+", "-"), "must recover a token from an unsanitized mol")

    def test_total_failure_is_logged_not_swallowed(self):
        class _Exploding:
            def GetNumConformers(self):
                raise RuntimeError("boom")

        with self.assertLogs("oinsmiles.oin.axial", level=logging.DEBUG) as cm:
            self.assertIsNone(mol_axial_token(_Exploding()))
        self.assertTrue(any("axial perception failed" in m for m in cm.output))


class TestNarrowingSeparatesMissFromBlindness(unittest.TestCase):
    """Bug 3, at the selection layer: a miss and a blind spot must not look identical."""

    def setUp(self):
        import oinsmiles.generation.metallogen_adapter as MA

        self.MA = MA

    def test_genuine_miss_reports_zero_blind(self):
        # every candidate perceived, none matching the request
        with unittest.mock.patch.object(self.MA, "mol_axial_token", return_value="+"):
            kept, blind = self.MA._axial_narrow(["a", "b"], lambda c: c, "-")
        self.assertEqual(kept, [])
        self.assertEqual(blind, 0, "perceived-but-different is NOT blindness")

    def test_perception_failure_is_counted(self):
        with unittest.mock.patch.object(self.MA, "mol_axial_token", return_value=None):
            kept, blind = self.MA._axial_narrow(["a", "b"], lambda c: c, "-")
        self.assertEqual(kept, [])
        self.assertEqual(blind, 2, "unperceived candidates must be counted separately")

    def test_match_is_kept(self):
        with unittest.mock.patch.object(self.MA, "mol_axial_token", side_effect=["-", "+"]):
            kept, blind = self.MA._axial_narrow(["hit", "miss"], lambda c: c, "-")
        self.assertEqual(kept, ["hit"])
        self.assertEqual(blind, 0)

    def test_blind_pool_warns(self):
        with self.assertLogs("oinsmiles.generation.metallogen_adapter", level=logging.WARNING):
            self.MA._axial_report_miss(3, 3, "-", key="n_scored")

    def test_ordinary_miss_does_not_warn(self):
        logger = logging.getLogger("oinsmiles.generation.metallogen_adapter")
        with unittest.mock.patch.object(logger, "warning") as warn:
            self.MA._axial_report_miss(3, 0, "-", key="n_scored")
        warn.assert_not_called()


class TestUnhonoredAxialIsReported(unittest.TestCase):
    """Bug 2: returning the wrong atropisomer must be announced, not just tolerated."""

    def setUp(self):
        import oinsmiles.generation.metallogen_adapter as MA

        self.MA = MA

    class _Parsed:
        def __init__(self, oin):
            self.original_oin = oin

    def test_mismatch_warns(self):
        with unittest.mock.patch.object(self.MA, "mol_axial_token", return_value="+"):
            with self.assertLogs(
                "oinsmiles.generation.metallogen_adapter", level=logging.WARNING
            ) as cm:
                self.MA._verify_axial_honored(self._Parsed("X |ax:-|"), object(), object())
        self.assertTrue(any("axial NOT honored" in m for m in cm.output))

    def test_match_is_quiet(self):
        logger = logging.getLogger("oinsmiles.generation.metallogen_adapter")
        with unittest.mock.patch.object(self.MA, "mol_axial_token", return_value="-"):
            with unittest.mock.patch.object(logger, "warning") as warn:
                self.MA._verify_axial_honored(self._Parsed("X |ax:-|"), object(), object())
        warn.assert_not_called()

    def test_no_token_requested_is_a_noop(self):
        # the default path: nothing requested, nothing perceived, nothing said
        logger = logging.getLogger("oinsmiles.generation.metallogen_adapter")
        with unittest.mock.patch.object(logger, "warning") as warn:
            self.MA._verify_axial_honored(self._Parsed(PLAIN_OIN), object(), object())
        warn.assert_not_called()


if __name__ == "__main__":
    unittest.main()
