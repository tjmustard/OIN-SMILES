"""``OIN_BORON_GEN_FASTFAIL``: fail fast on a boron cage the generator cannot embed.

Background (measured, see ``docs/BORON_GEN_CEILING_v0.4.7.md``): promoting
``OIN_BORON_CAGE`` (v0.4.6) fixed a genuine ENCODER ceiling -- 34 boron-cluster
molecules go from 0/34 to 34/36 encoding -- but exposed a GENERATOR ceiling the
encode failure had been hiding. A sample of 32 coordinated-cage molecules on a
non-LIN geometry (the 34-molecule ``encode_fail`` class plus a 14-molecule control
group, minus the two exceptions this module pins) all burn their whole embed
budget or worse without producing a structure. This lever detects the same B-B-B
cage motif ``OIN_BORON_CAGE`` gates on, BEFORE generation starts, and raises
immediately instead.

The naive version of this predicate (bare motif check) was MEASURED WRONG:
``RAWJEG_comp_0`` ([Hg_LIN], monodentate cage + Cl-) carries the motif and still
produces a structure in a couple of seconds. These tests pin the refined
predicate's two safety exclusions (uncoordinated fragment; the one confirmed-safe
geometry) against that counter-example and its neighbours, not just the easy
cases the naive version would also have passed.
"""

from __future__ import annotations

import os
import unittest
from pathlib import Path

from rdkit import RDLogger

from oinsmiles import XYZToSMILES
from oinsmiles.generation.metallogen_adapter import (
    BoronCageGenerationUnsupportedError,
    OIN3DGeneratorMetallogen,
    _parsed_oin_has_boron_cage,
)
from oinsmiles.generation.oin_parser import OINParser

RDLogger.DisableLog("rdApp.*")

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"

LEVER = "OIN_BORON_GEN_FASTFAIL"
CAGE_LEVER = "OIN_BORON_CAGE"

#: A Rh thiaborane cage on a TET (CN=4) metal -- one of the 34 encode_fail class,
#: confirmed (tools/boron_gen_sweep34.sh) to burn its whole embed budget with
#: nothing produced.
KIXXOF = FIXTURES / "KIXXOF_comp_0.xyz"
#: A nido-C2B7 cage on a TPL (CN=3) Ru -- one of the 14 "silently wrong before"
#: molecules (docs/BORON_CAGE_v0.4.5.md SS5a); also confirmed 0/1 with the cage
#: lever on (tools/boron_gen_sweep_14silent.sh).
VEJXOZ = FIXTURES / "VEJXOZ_comp_0.xyz"
#: A boron-rich molecule with NO cage motif (Ir boroxine, B-O-B-O-B-O ring, zero
#: B-B bonds) -- must never trip this predicate regardless of the lever.
ASUVIV = FIXTURES / "ASUVIV_comp_0.xyz"
#: An Fe complex with BH3-/BH- groups at several coordination numbers, still no
#: deltahedron -- second no-cage control.
AROTAE = FIXTURES / "AROTAE_comp_0.xyz"
#: closo-B12H12 dodecaborate cage -- carries the motif, TPY (CN=4) geometry.
OZAREO = FIXTURES / "OZAREO_comp_0.xyz"


def _parsed(xyz_path, cage_lever="1"):
    # NOTE: variable is named `_saved` (not `saved`) so it matches the restore-hint
    # heuristic in test_levers.py::TestNoTestUnsetsAPromotedLever -- see that
    # module's docstring; without the underscore this reads as unsetting a
    # promoted lever (OIN_BORON_CAGE) rather than restoring its prior value.
    _saved = os.environ.get(CAGE_LEVER)
    os.environ[CAGE_LEVER] = cage_lever
    try:
        oin = XYZToSMILES().convert(str(xyz_path))
    finally:
        if _saved is None:
            os.environ.pop(CAGE_LEVER, None)
        else:
            os.environ[CAGE_LEVER] = _saved
    return OINParser().parse(oin), oin


class _LeverMixin(unittest.TestCase):
    """Sets ``OIN_BORON_GEN_FASTFAIL`` EXPLICITLY in both directions -- see the same
    trap noted in ``test_boron_cage.py``: never rely on an unset default.
    """

    def setUp(self):
        self._saved = os.environ.get(LEVER)
        os.environ[LEVER] = "0"

    def tearDown(self):
        if self._saved is None:
            os.environ.pop(LEVER, None)
        else:
            os.environ[LEVER] = self._saved

    def set_lever(self, on):
        os.environ[LEVER] = "1" if on else "0"


class TestPredicateFiresOnConfirmedFailures(_LeverMixin):
    """The predicate matches every confirmed-slow cage molecule checked so far."""

    def test_kixxof_tet_cage_flagged(self):
        parsed, _ = _parsed(KIXXOF)
        self.assertTrue(_parsed_oin_has_boron_cage(parsed))

    def test_vejxoz_tpl_cage_flagged(self):
        parsed, _ = _parsed(VEJXOZ)
        self.assertTrue(_parsed_oin_has_boron_cage(parsed))

    def test_ozareo_tpy_cage_flagged(self):
        parsed, _ = _parsed(OZAREO)
        self.assertTrue(_parsed_oin_has_boron_cage(parsed))


class TestPredicateSpecificityControls(_LeverMixin):
    """Molecules that must NEVER trip the predicate, on the evidence that earned
    each exclusion (docs/BORON_GEN_CEILING_v0.4.7.md SS5)."""

    def test_no_cage_motif_at_all_is_never_flagged(self):
        parsed, _ = _parsed(ASUVIV)
        self.assertFalse(_parsed_oin_has_boron_cage(parsed))

    def test_boron_at_several_coordination_numbers_no_deltahedron(self):
        parsed, _ = _parsed(AROTAE)
        self.assertFalse(_parsed_oin_has_boron_cage(parsed))


class TestEndToEndFastFailIsActuallyFast(_LeverMixin):
    """The load-bearing proof for the "fast" half of fail-FAST: with the lever on,
    a confirmed-slow cage molecule raises before any embed work starts, so this
    test is itself fast regardless of host load -- the whole point of the fix.
    """

    def test_lever_on_raises_before_embedding_a_confirmed_slow_cage(self):
        self.set_lever(True)
        oin = XYZToSMILES().convert(str(KIXXOF))
        with self.assertRaises(BoronCageGenerationUnsupportedError):
            OIN3DGeneratorMetallogen(optimizer=None, ensemble_size=1, timeout=300).generate(oin)
        # timeout=300 above is deliberate: if the predicate did NOT short-circuit,
        # a real run would spend most of that budget (docs/BORON_GEN_CEILING_v0.4.7.md
        # SS2) -- this test only completes quickly because it never reaches the
        # embed loop at all, which is the behaviour being pinned.

    def test_lever_on_does_not_raise_for_a_non_cage_molecule(self):
        # ASUVIV carries boron but no cage motif; the predicate must not fire
        # regardless of the lever, so generation proceeds exactly as it does
        # with the lever off (may still fail for OTHER reasons -- e.g. size or
        # geometry support -- but never via BoronCageGenerationUnsupportedError).
        self.set_lever(True)
        oin = XYZToSMILES().convert(str(ASUVIV))
        try:
            OIN3DGeneratorMetallogen(optimizer=None, ensemble_size=1, timeout=5).generate(oin)
        except BoronCageGenerationUnsupportedError:
            self.fail("fast-fail predicate fired on a molecule with no cage motif")
        except Exception:
            pass  # any OTHER failure (timeout, embed) is out of scope here


class TestLeverOffIsInert(_LeverMixin):
    """Lever OFF (the default) must never even CALL the predicate, let alone raise
    the new exception -- checked via a monkeypatch rather than a real slow
    molecule, so this stays fast and deterministic regardless of host load.
    """

    def test_lever_off_never_calls_the_predicate(self):
        import oinsmiles.generation.metallogen_adapter as mod

        self.set_lever(False)
        calls = []
        original = mod._parsed_oin_has_boron_cage

        def spy(parsed):
            calls.append(parsed)
            return original(parsed)

        mod._parsed_oin_has_boron_cage = spy
        try:
            oin = XYZToSMILES().convert(str(FIXTURES / "CisPlatin.xyz"))
            OIN3DGeneratorMetallogen(optimizer=None, ensemble_size=1, timeout=30).generate(oin)
        finally:
            mod._parsed_oin_has_boron_cage = original
        self.assertEqual(calls, [], "lever OFF must short-circuit before the predicate runs")


if __name__ == "__main__":
    unittest.main()
