"""Regression tests for donor-atom hydrogen handling in the MetalloGen m-SMILES.

`convert_parsed_to_msmiles` decides, per binding atom, whether the metal bond
replaces a hydrogen. The rule: honor an EXPLICIT bracket-H count from the OIN
(a neutral dative L-type donor keeps its H) and only reinterpret a BARE binding
atom (its implicit H is a phantom the metal bond replaces).

These guard the two generation defects the tmCAT/tmPHOTO atom-count check
surfaced: an NHC carbene carbon built as CH2 (+2 H), and a dative secondary
amine stripped of its N-H (-1 H).
"""

import os
import re
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.metallogen_adapter import convert_oin_to_msmiles


class TestMSmilesDonorHydrogens(unittest.TestCase):
    def test_nhc_carbene_is_zero_h(self):
        """NHC carbene carbon must be a 0-H [C:n] donor, never a [CH2:n] methylene."""
        oin = "[Pd_SPL].C=CCN1C{0}N(CC=C)c2ccccc21.[Br]{1}.C=CCN1C{2}N(CC=C)c2ccccc21.[Br]{3}"
        msmiles = convert_oin_to_msmiles(oin)
        self.assertNotIn("[CH2:", msmiles, f"carbene built as CH2 in: {msmiles}")
        # Two mapped bare carbene carbons expected.
        self.assertGreaterEqual(len(re.findall(r"\[C:\d+\]", msmiles)), 2, msmiles)

    def test_dative_secondary_amine_keeps_h(self):
        """A dative amine written [NH] keeps its N-H; the amido N (bare) is stripped."""
        oin = (
            "[Ti_OCT].[Cl]{0}.[Cl]{1}"
            ".CC(C)c1cc(C(C)C)c(-c2cccc(N{2}c3c(C(C)C)cccc3C(C)C)n{4}2)c(C(C)C)c1"
            ".C[NH]{3}C.[Cl]{5}"
        )
        msmiles = convert_oin_to_msmiles(oin)
        # dative dimethylamine keeps its H
        self.assertRegex(msmiles, r"C\[NH:\d+\]C", f"dative amine lost its H in: {msmiles}")
        # anionic amido N (bare in the OIN) is stripped to 0 H
        self.assertRegex(msmiles, r"\[N:\d+\]", f"amido N not stripped in: {msmiles}")

    def test_carbon_monoxide_carbon_has_no_h(self):
        """CO carbon (triple bond) stays a 0-H donor (pre-existing behavior guard)."""
        oin = "[Re_OCT].C{0}#O.[Br]{1}.[Br]{2}.C{3}#O.C{4}#O.C{5}#O"
        msmiles = convert_oin_to_msmiles(oin)
        self.assertNotIn("[CH", msmiles, msmiles)

    def test_bare_alkoxide_oxygen_is_stripped(self):
        """A bare chalcogen donor (alkoxide) is an anionic 0-H donor (pre-existing guard)."""
        oin = "[Ti_OCT].[Cl]{0}.[Cl]{1}.[Cl]{2}.[Cl]{3}.CCO{4}.CCO{5}"
        msmiles = convert_oin_to_msmiles(oin)
        # ethoxide O bound to one C + metal -> [O:n], not [OH:n]
        self.assertNotIn("[OH:", msmiles, msmiles)
        self.assertRegex(msmiles, r"\[O:\d+\]", msmiles)


if __name__ == "__main__":
    unittest.main()
