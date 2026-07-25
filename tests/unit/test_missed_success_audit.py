"""Guards for the Wave-3 missed-success attribution (Y3).

The audit's whole value is that it does NOT score generator failures against the notation.
Its headline claim -- "77.8 % of round-trip failures never test the OIN" -- rests entirely on
the attribution rules below, so they are pinned here. In particular the audit must refuse to
call a divergent isomer evidence against the encoder: that case is equally consistent with
the generator having built the wrong structure, and conflating them is exactly the reasoning
error the Y1 audit exists to prevent.
"""

import sys
import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from tools.injectivity.missed_success_audit import (  # noqa: E402
    TIMEOUT_S,
    UNINFORMATIVE,
    attribute,
)

OIN_A = "[Pd_SPL].c1ccccc1P{0}.[Cl]{1}.[Cl]{2}"
OIN_B = "[Pd_SPL].c1ccccc1N{0}.[Cl]{1}.[Cl]{2}"  # a different ligand -> different isomer


def _row(**kw):
    base = {
        "bucket": "hard_fail",
        "error": None,
        "elapsed_s": 5.0,
        "smiles_1": OIN_A,
        "smiles_2": None,
    }
    base.update(kw)
    return base


class TestGeneratorCauses(unittest.TestCase):
    def test_explicit_timeout_error(self):
        cause, _ = attribute(_row(error="generation timeout after 300s", elapsed_s=300.1))
        self.assertEqual(cause, "generator_timeout")

    def test_at_the_wall_without_output_is_a_timeout(self):
        # no error text, but it burned the whole budget and produced nothing
        cause, _ = attribute(_row(error=None, elapsed_s=TIMEOUT_S, smiles_2=None))
        self.assertEqual(cause, "generator_timeout")

    def test_fast_death_is_not_a_timeout(self):
        cause, _ = attribute(_row(error="MetalloGen failed", elapsed_s=3.0))
        self.assertEqual(cause, "generator_no_output")

    def test_generator_causes_are_uninformative(self):
        for cause in ("generator_timeout", "generator_no_output"):
            self.assertIn(cause, UNINFORMATIVE)


class TestOutputComparison(unittest.TestCase):
    def test_same_isomer_is_canonicalization_noise(self):
        # identical strings trivially share a key -> the test was too strict, not the OIN wrong
        cause, _ = attribute(_row(smiles_2=OIN_A, elapsed_s=12.0))
        self.assertEqual(cause, "canonicalization_noise")

    def test_different_isomer_is_ambiguous_not_notation_evidence(self):
        cause, _ = attribute(_row(smiles_2=OIN_B, elapsed_s=12.0))
        self.assertEqual(cause, "divergent_isomer_ambiguous")
        self.assertNotIn(cause, UNINFORMATIVE, "it is not *uninformative* -- it is unattributable")

    def test_divergence_is_never_scored_against_the_notation(self):
        # the load-bearing guard: no attribution may claim a lossy OIN, because the sweep
        # cannot distinguish generator error from notation ambiguity.
        causes = {attribute(_row(smiles_2=s, elapsed_s=12.0))[0] for s in (OIN_A, OIN_B)}
        self.assertNotIn("lossy_oin", causes)
        self.assertNotIn("genuinely_different", causes)


class TestEncodeFail(unittest.TestCase):
    def test_encode_fail_bucket_wins(self):
        # an encoder refusal is a coverage limit, not a round-trip result -- even if the row
        # also looks like it timed out.
        cause, _ = attribute(_row(bucket="encode_fail", error="timeout", elapsed_s=300.0))
        self.assertEqual(cause, "encode_fail")
        self.assertNotIn(cause, UNINFORMATIVE)


if __name__ == "__main__":
    unittest.main()
