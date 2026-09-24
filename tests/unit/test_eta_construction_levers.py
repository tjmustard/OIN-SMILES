"""Guards for v0.4.18 L2: an eta ring must be BUILT inside the metal's bonding range.

The census filed 343 molecules / 6.86 pts under ``G_CONSTRUCTION / DETACHED``, 86% eta-bound, and
the defect was a distance: ``ff_clean`` aimed an eta5 ring at 1.2 x the covalent-radius sum (the
edge of the encoder's contact cutoff), and for a small metal the FF scan's vdW guard read the
ring's own approach as a clash and reverted it. FERROCENE ITSELF was generated detached -- Fe-C
2.85 A against a real 2.05 A, coordination not intact, honest round trip failing.

``OIN_ETA_COVALENT_TARGET`` + ``OIN_VDW_EXEMPT_BINDING`` were promoted to default-ON by the owner
(2026-09-20) on a harness A/B over all 1,146 eta-bound molecules: +195/-29 self-consistent,
+189/-39 VERIFIED, noise floor zero. ``OIN_ETA_TARGET_UNSCALED`` was measured as a full third arm
and is dominated -- it stays held off.

Levers are WRITTEN ("0" / "1"), never unset: ``test_levers::TestNoTestUnsetsAPromotedLever``.
"""

import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

import numpy as np
from rdkit import RDLogger

from oinsmiles import XYZToSMILES
from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen
from oinsmiles.oin import levers
from oinsmiles.oin.coordination import coordination_report

RDLogger.DisableLog("rdApp.*")

FERROCENE = os.path.join(os.path.dirname(__file__), "../fixtures/Ferrocene.xyz")
PROMOTED = ("OIN_ETA_COVALENT_TARGET", "OIN_VDW_EXEMPT_BINDING")
UPSTREAM = {name: "0" for name in PROMOTED}


class TestRegistry(unittest.TestCase):
    def test_the_pair_is_promoted_and_not_held_off(self):
        for name in PROMOTED:
            self.assertIn(name, levers._DEFAULT_ON)
            self.assertNotIn(name, levers._HELD_OFF)

    def test_zero_still_selects_upstream(self):
        with mock.patch.dict(os.environ, UPSTREAM):
            for name in PROMOTED:
                self.assertFalse(levers.lever_enabled(name))

    def test_the_dominated_third_lever_stays_held_off_with_its_verdict(self):
        self.assertNotIn("OIN_ETA_TARGET_UNSCALED", levers._DEFAULT_ON)
        self.assertIn("DOMINATED", levers._HELD_OFF["OIN_ETA_TARGET_UNSCALED"])


def _build_ferrocene():
    """(mean Fe-C of the ten ring carbons, coordination intact, honest round trip)."""
    oin = XYZToSMILES().convert(FERROCENE)
    res = OIN3DGeneratorMetallogen(optimizer=None, ensemble_size=1, timeout=120).generate(oin)
    rows = [ln.split() for ln in res.xyz.strip().splitlines()[2:]]
    xyz = np.array([[float(v) for v in r[1:4]] for r in rows])
    fe = next(i for i, r in enumerate(rows) if r[0] == "Fe")
    d_c = sorted(np.linalg.norm(xyz[i] - xyz[fe]) for i, r in enumerate(rows) if r[0] == "C")
    with tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False) as fh:
        fh.write(res.xyz)
    try:
        # the HONEST round trip: an independent encode of the coordinates, never the
        # generator's own bond graph
        back = XYZToSMILES().convert(fh.name)
    finally:
        os.unlink(fh.name)
    with open(FERROCENE) as fh:
        intact = coordination_report(fh.read(), res.xyz)["intact"]
    return float(np.mean(d_c[:10])), bool(intact), back == oin


class TestFerroceneIsBuiltAttached(unittest.TestCase):
    """Fails against pre-L2 code and against either lever alone being dropped from the default:
    Fe is a dead-zone metal, so the covalent target without the exemption leaves the ring at the
    embed's ~2.87 A."""

    def test_shipped_defaults_build_a_bonded_ring_that_round_trips(self):
        mean_fe_c, intact, round_trip = _build_ferrocene()
        self.assertLess(mean_fe_c, 2.30, "real ferrocene is 2.05 A; the contact cutoff is 2.53 A")
        self.assertTrue(intact)
        self.assertTrue(round_trip)

    def test_upstream_behaviour_is_the_defect_and_zero_brings_it_back(self):
        # pins WHAT WAS WRONG, so the A/B's OFF arm stays reproducible and nobody "fixes" the
        # "0" path into the shipped one
        with mock.patch.dict(os.environ, UPSTREAM):
            mean_fe_c, intact, round_trip = _build_ferrocene()
        self.assertGreater(mean_fe_c, 2.60)
        self.assertFalse(intact)
        self.assertFalse(round_trip)


if __name__ == "__main__":
    unittest.main()
