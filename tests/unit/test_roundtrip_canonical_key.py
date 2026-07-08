"""Unit tests for the round-trip canonical comparison key.

`canonical_roundtrip_key` decides whether two OIN strings describe the *same*
coordination compound. It must collapse chemically-meaningless notation drift
(the dominant tmCAT/tmPHOTO round-trip failure mode) while still keeping
genuinely-different structures, metals, geometries, and eta windings distinct.

Every Exp/Got pair below is a real case taken from a tmCAT/tmPHOTO
`summary_roundtrip.json` run (molecule id in the test name), so the assertions
double as a regression guard for the Track A2 comparator relaxation.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.oin.compare import canonical_roundtrip_key


class TestCanonicalRoundtripKey(unittest.TestCase):
    def assertEquivalent(self, exp, got):
        self.assertEqual(
            canonical_roundtrip_key(exp),
            canonical_roundtrip_key(got),
            f"expected these to be the SAME structure:\n  Exp: {exp}\n  Got: {got}",
        )

    def assertDistinct(self, exp, got):
        self.assertNotEqual(
            canonical_roundtrip_key(exp),
            canonical_roundtrip_key(got),
            f"expected these to be DIFFERENT structures (comparator must not mask it):"
            f"\n  Exp: {exp}\n  Got: {got}",
        )

    # -- Notation-only drift: MUST collapse (same molecule, different string) --

    def test_carboxylate_symmetric_o_swap(self):
        """ABAZIO: which of the two equivalent carboxylate O carries the slot."""
        exp = (
            "[Pd_SPL].CN(C)c1ccc(P{0}(C(C)(C)C)C(C)(C)C)cc1.CN{1}(C)Cc1c{2}cccc1.O{3}C(=O)C(F)(F)F"
        )
        got = (
            "[Pd_SPL].CN(C)c1ccc(P{0}(C(C)(C)C)C(C)(C)C)cc1.CN{1}(C)Cc1c{2}cccc1.OC(=O{3})C(F)(F)F"
        )
        self.assertEquivalent(exp, got)

    def test_binder_explicit_vs_implicit_h(self):
        """ABESAD: metal-binding amine written [NH] (explicit) vs N (implicit)."""
        exp = "[Ti_OCT].[Cl]{0}.[Cl]{1}.CC(C)c1cc(C(C)C)c(-c2cccc(N{2}c3c(C(C)C)cccc3C(C)C)n{3}2)c(C(C)C)c1.C[NH]{4}C.[Cl]{5}"
        got = "[Ti_OCT].[Cl]{0}.[Cl]{1}.CC(C)c1cc(C(C)C)c(-c2cccc(N{2}c3c(C(C)C)cccc3C(C)C)n{3}2)c(C(C)C)c1.CN{4}C.[Cl]{5}"
        self.assertEquivalent(exp, got)

    def test_nhc_carbene_bare_c_vs_ch2(self):
        """ACAWOR: NHC carbene carbon serialized as bare C vs [CH2]."""
        exp = "[Pd_SPL].C=CCN1C{0}N(CC=C)c2ccccc21.[Br]{1}.C=CCN1C{2}N(CC=C)c2ccccc21.[Br]{3}"
        got = (
            "[Pd_SPL].C=CCN1[CH2]{0}N(CC=C)c2ccccc21.[Br]{1}.C=CCN1[CH2]{2}N(CC=C)c2ccccc21.[Br]{3}"
        )
        self.assertEquivalent(exp, got)

    def test_fragment_order_and_slot_permutation(self):
        """ACEPUT: CO and Br ligands swap position/slot but the compound is identical."""
        exp = "[Re_OCT].C{0}#O.[Br]{1}.Clc1ccn{2}c2c1ccc1ccc3c(Cl)ccn{3}c3c21.C{4}#O.C{5}#O"
        got = "[Re_OCT].[Br]{0}.C{1}#O.Clc1ccn{2}c2c1ccc1ccc3c(Cl)ccn{3}c3c21.C{4}#O.C{5}#O"
        self.assertEquivalent(exp, got)

    def test_bound_water_notation(self):
        """[OH2] and O are the same bound-water ligand."""
        exp = "[Zn_SPL].[OH2]{0}.[Cl]{1}"
        got = "[Zn_SPL].O{0}.[Cl]{1}"
        self.assertEquivalent(exp, got)

    # -- Genuinely different structures: MUST stay distinct (no masking) --

    def test_allyl_double_bond_loss(self):
        """ABAZEK: eta-allyl double bonds lost + pendant phenyl dearomatized."""
        exp = "[Pd_TPL].[CH2]{0>}[CH]{0}=[CH]{0}c1ccccc1.CN(C)c1ccc(P{1}(C(C)(C)C)C(C)(C)C)cc1.[Cl]{2}"
        got = "[Pd_TPL].[CH2]{0>}[CH]{0}[CH]{0}C1[CH][CH][CH][CH][CH]1.CN(C)c1ccc(P{1}(C(C)(C)C)C(C)(C)C)cc1.[Cl]{2}"
        self.assertDistinct(exp, got)

    def test_p_stereocenter_loss(self):
        """ABOPOY: chiral phosphorus [P@] regenerated as achiral P."""
        exp = "[Zn_TPL].c{0}1ccccc1.CC(C)c1cccc(C(C)C)c1N{1}=[P@](C)(c1ccccc1)c1cccc2c1oc1ccccc21.c{2}1ccccc1"
        got = "[Zn_TPL].c{0}1ccccc1.CC(C)c1cccc(C(C)C)c1N{1}=P(C)(c1ccccc1)c1cccc2c1oc1ccccc21.c{2}1ccccc1"
        self.assertDistinct(exp, got)

    def test_eta_winding_flip(self):
        """ADEYIR: eta-arene winding flips (>/<) -- a real face/stereo difference."""
        exp = "[Ta_SPY].Cc{0}1c{0>}(C)c{0}(C)c{0}(C)c{0}1C.CC(C)(C)P(=N{1})(C(C)(C)C)C(C)(C)C.[Cl]{2}.[Cl]{3}.[Cl]{4}"
        got = "[Ta_SPY].Cc{0}1c{0<}(C)c{0}(C)c{0}(C)c{0}1C.CC(C)(C)P(=[NH]{1})(C(C)(C)C)C(C)(C)C.[Cl]{2}.[Cl]{3}.[Cl]{4}"
        self.assertDistinct(exp, got)

    def test_geometry_spy_vs_tbp(self):
        """AFEJAX: five-coordinate geometry reassigned SPY -> TBP."""
        exp = "[Zn_SPY].O{0}.Cc1cccc(C=N{1}c2ccccc2N{2}=Cc2cccc(C)c2O{3})c1O{4}"
        got = "[Zn_TBP].Cc1cccc(C=N{0}c2ccccc2N{1}=Cc2cccc(C)c2O{2})c1O{3}.O{4}"
        self.assertDistinct(exp, got)

    def test_different_metal(self):
        """A different metal element is never the same compound."""
        exp = "[Pd_SPL].[Cl]{0}.[Cl]{1}"
        got = "[Pt_SPL].[Cl]{0}.[Cl]{1}"
        self.assertDistinct(exp, got)


if __name__ == "__main__":
    unittest.main()
