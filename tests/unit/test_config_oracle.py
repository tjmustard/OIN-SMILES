"""Configurational-oracle guards (Y2): the three blind-spot axes are RECOVERABLE from 3D.

Wave 1 proved the encoder is *blind* to metal Δ/Λ, axial atropisomerism, and metal-bound
amine chirality. Wave 2's feasibility claim is the complement: the information the encoder
discards is still present in the input coordinates. ``tools/injectivity/config_oracle``
recovers a per-axis descriptor that FLIPS between enantiomers; these guards lock that in.

If any of these ever stops flipping, the recovery route regressed and the "recoverable, not
permanent" disposition in ``docs/INJECTIVITY_Y2_FEASIBILITY.md`` no longer holds.
"""

import sys
import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from tools.injectivity.config_oracle import (  # noqa: E402
    configurational_signature,
    load_mol,
    mirror_flip_report,
)

FIX = _ROOT / "tests" / "fixtures"
CISPLATIN = FIX / "CisPlatin.xyz"
FAC_IRPPY3 = FIX / "fac-Ir(ppy)3.xyz"
BINAP = FIX / "PdCl2-R-BINAP.xyz"
POJJOP = FIX / "POJJOP.xyz"


class TestMetalStereoRecovered(unittest.TestCase):
    """P1 metal Δ/Λ: RDKit perceives the octahedral permutation from 3D and it flips."""

    def test_fac_irppy3_metal_permutation_flips(self):
        rep = mirror_flip_report(FAC_IRPPY3)
        self.assertTrue(rep.metal_flips, "octahedral Δ/Λ permutation must flip for the mirror")
        self.assertEqual(rep.base.metal[0].shape, "OH")

    def test_cisplatin_metal_does_not_flip(self):
        # square-planar cisplatin is achiral: the metal descriptor must NOT flip.
        rep = mirror_flip_report(CISPLATIN)
        self.assertFalse(rep.metal_flips)


class TestAxialRecovered(unittest.TestCase):
    """P2 axial: the signed biaryl dihedral is the recoverable configuration."""

    def test_binap_axis_detected_and_flips(self):
        rep = mirror_flip_report(BINAP)
        self.assertTrue(rep.axial_flips, "the BINAP biaryl axis sign must flip for the mirror")
        self.assertEqual(len(rep.base.axial), 1, "exactly one hindered biaryl axis in BINAP")

    def test_ppy_planar_biaryl_not_flagged(self):
        # the chelated phenyl-pyridyl bonds of fac-Ir(ppy)3 are planar, not atropisomeric.
        sig = configurational_signature(load_mol(FAC_IRPPY3))
        self.assertEqual(sig.axial, [], "planar chelated biaryls are not atropisomer axes")


class TestBoundAmineRecovered(unittest.TestCase):
    """P3 metal-bound amine: signed tetrahedral volume at the locked N is recoverable."""

    def test_pojjop_bound_amine_flips(self):
        rep = mirror_flip_report(POJJOP)
        self.assertTrue(rep.amine_flips, "metal-bound amine configuration must flip for the mirror")

    def test_cisplatin_ammine_not_a_stereocentre(self):
        # NH3 ligands (M,H,H,H) have equivalent H's -- not stereogenic, must be excluded.
        sig = configurational_signature(load_mol(CISPLATIN))
        self.assertEqual(sig.bound_amine, [])


class TestControlIsInert(unittest.TestCase):
    def test_cisplatin_recovers_nothing(self):
        rep = mirror_flip_report(CISPLATIN)
        self.assertFalse(rep.any_recovered, "the achiral control must recover no configuration")


if __name__ == "__main__":
    unittest.main()
