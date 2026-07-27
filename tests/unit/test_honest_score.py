"""The honest round-trip verdict, pinned on the three fixtures that define it.

WHAT THIS GUARDS
================
A round-trip verdict may not consult any artifact the generator produced except the
coordinates. Bonds, stereo and coordination must be re-derived from the XYZ alone.

The harness used to score with ``get_oin_string(gen_result.mol, coords)`` -- the generator's
own bond graph. That is not merely inaccurate, it is **circular**: ``gen_result.mol`` is
exactly the artifact that would have to be wrong for the test to fail. One shortcut produced
errors in both directions at once, and all three are pinned here.

Fixtures are vendored under ``tests/fixtures/honest_score/`` rather than read from a sweep
directory. The three molecules live in ``results-v0.4.5-rebaseline`` and **none of them is in
the 5000-molecule v0.4.6 corpus** -- a test pointed at the sweep would pass vacuously on a
machine that has it and fail on one that does not. The generated geometry is the exact string
the generator emitted, so no generator run is needed and the test is deterministic.
"""

import os
import unittest

from oinsmiles import XYZToSMILES
from oinsmiles.oin.compare import canonical_roundtrip_key
from oinsmiles.oin.coordination import coordination_report
from oinsmiles.oin.levers import default_on

_FIX = os.path.join(os.path.dirname(__file__), "..", "fixtures", "honest_score")

#: ``smiles_1`` exactly as the v0.4.5 rebaseline sweep recorded it. Hard-coded rather than
#: recomputed: the point of a golden is to fail when the encoder moves, and an expectation
#: derived from the encoder under test cannot do that.
_INPUT_OIN = {
    "FIYHUT_comp_0": (
        "[Fe_LIN].CCCCCN(C)(C)Cc{0}1[cH]{0}[cH]{0>}[cH]{0}[cH]{0}1."
        "[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
    ),
    "YOSXIP_comp_0": "[Ru_OCT].CS{0}(C)=O.CS{1}CC[S@]{5}(C)=O.CS{2}(C)=O.[Cl]{3}.[Cl]{4}",
    "OGARAP_comp_0": (
        "[Pd_TPL].[CH]{0<}([CH]{0}=C{0}(c1ccccc1)c1ccccc1)c1ccccc1."
        "c1ccc([C@@H]2COC(c3ccccc3-p{2}3c4ccccc4c4ccccc34)=N{1}2)cc1"
    ),
}

#: ``smiles_2`` as the sweep scored it -- ``get_oin_string(gen_result.mol, coords)``.
_SCORED_OIN = {
    "FIYHUT_comp_0": _INPUT_OIN["FIYHUT_comp_0"],  # scored a byte-exact PASS
    "YOSXIP_comp_0": "[Ru_OCT].CS{0}(C)=O.CS{1}CCS{5}(C)=O.CS{2}(C)=O.[Cl]{3}.[Cl]{4}",
    "OGARAP_comp_0": _INPUT_OIN["OGARAP_comp_0"],  # scored a byte-exact PASS
}


def _generated(name):
    return os.path.join(_FIX, f"{name}_generated.xyz")


def _honest(name):
    """The honest round-trip string: coordinates in, nothing else."""
    return XYZToSMILES().convert(_generated(name))


class TestHonestScoreFixtures(unittest.TestCase):
    def test_fiyhut_scored_pass_is_a_honest_failure(self):
        """§1 -- the generator's graph ASSERTS bonds the coordinates do not support.

        Both Cp rings sit 0.85 A off the iron: Fe-C goes 2.02-2.05 A in the input to
        2.84-2.96 A out, so all 10 bonded carbons are gone. ``gen_result.mol`` still calls
        them bonded, the string matches byte-for-byte, and the molecule scores a pass.
        """
        name = "FIYHUT_comp_0"
        self.assertEqual(
            _SCORED_OIN[name],
            _INPUT_OIN[name],
            "precondition: the scored arm called this byte-exact",
        )
        honest = _honest(name)
        self.assertNotEqual(_INPUT_OIN[name], honest, "honest arm must NOT call this byte-exact")
        self.assertNotEqual(
            canonical_roundtrip_key(_INPUT_OIN[name]),
            canonical_roundtrip_key(honest),
            "the key must not rescue it either -- this is a real coordination loss",
        )

    def test_fiyhut_coordination_agrees(self):
        """The cheap coordinate-only tripwire must corroborate the expensive judge."""
        with open(os.path.join(_FIX, "FIYHUT_comp_0_input.xyz")) as a:
            with open(_generated("FIYHUT_comp_0")) as b:
                report = coordination_report(a.read(), b.read())
        self.assertIs(report["intact"], False)
        self.assertIn("lost", report["reason"])

    def test_yosxip_scored_failure_is_a_honest_pass(self):
        """§6 -- the generator's graph LACKS stereo the coordinates do support.

        ``[S@]{5}`` is a genuine sulfoxide stereocentre present in the generated geometry.
        ``gen_result.mol`` flattens it to ``S{5}``, so the harness scored a ``String
        mismatch`` on a round trip that in fact succeeded. Re-perceiving from coordinates
        recovers the tag.
        """
        name = "YOSXIP_comp_0"
        self.assertNotEqual(
            _SCORED_OIN[name],
            _INPUT_OIN[name],
            "precondition: the scored arm called this a mismatch",
        )
        self.assertIn("[S@]{5}", _INPUT_OIN[name])
        self.assertIn("S{5}", _SCORED_OIN[name])
        self.assertEqual(
            _INPUT_OIN[name], _honest(name), "honest arm must recover the stereocentre exactly"
        )

    def test_ogarap_is_caught_only_by_the_honest_arm(self):
        """The blind spot that justifies keeping BOTH instruments.

        ``coordination.py`` is loss-based: it compares per-element contact COUNTS, so an
        eta-3 -> eta-2 rearrangement that keeps the count is invisible to it. The honest
        re-perception sees it, because the winding token changes. Neither instrument
        dominates the other -- that is the whole argument for running both.
        """
        name = "OGARAP_comp_0"
        with open(os.path.join(_FIX, f"{name}_input.xyz")) as a:
            with open(_generated(name)) as b:
                report = coordination_report(a.read(), b.read())
        self.assertIsNot(report["intact"], False, "coordination is expected to be BLIND here")
        self.assertNotEqual(_INPUT_OIN[name], _honest(name), "the honest arm must still catch it")


class TestHonestScoreIsTheDefault(unittest.TestCase):
    def test_indep_score_lever_ships_enabled(self):
        """v0.4.8 promoted the honest verdict from a diagnostic to the reported number.

        If this ever flips back, every downstream accuracy figure silently reverts to the
        over-stated one, and nothing else in the suite would notice.
        """
        self.assertIn("OIN_INDEP_SCORE", default_on())


if __name__ == "__main__":
    unittest.main()
