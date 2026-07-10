"""Regression tests for BARE anionic donor atoms in the MetalloGen m-SMILES.

`convert_parsed_to_msmiles` decides, per binding atom, whether the metal bond
replaces a hydrogen. The bracket/bare split is made upstream by ``replace_map``
in ``oin/inline.py``, which de-brackets a binding atom only when its bracket
content is a bare organic-subset symbol. So an N donor that keeps a hydrogen
always serializes as ``[NH2]``/``[NH]``; the sole bare H-bearing N is ammine
NH3, which has no heavy neighbour.

A bare ``N{n}`` with at least one heavy neighbour therefore carries exactly zero
hydrogens -- amido, anilide, silylamide, azide, phosphinimide. These tests pin
that rule from both sides: the X-type donors must be stripped to 0 H, and the
L-type donors that legitimately keep H must not be touched.

Companion file: ``test_msmiles_donor_hydrogens.py`` (NHC carbene, dative amine,
CO, alkoxide). Kept separate so the pre-existing guards stay readable.
"""

import os
import re
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.metallogen_adapter import (  # noqa: E402
    UncoordinatedFragmentError,
    convert_oin_to_msmiles,
)


class TestBareDonorStripped(unittest.TestCase):
    """A bare N donor with >=1 heavy neighbour is a 0-H X-type donor."""

    def assertZeroHNitrogen(self, msmiles):
        self.assertNotRegex(msmiles, r"\[NH\d?:", f"donor kept a hydrogen in: {msmiles}")
        self.assertRegex(msmiles, r"\[N:\d+\]", f"no stripped N donor in: {msmiles}")

    def test_silylamide_nitrogen_is_zero_h(self):
        """C[Si](C)(C)N{1} -- one heavy neighbour (Si). XADYAC_comp_0."""
        oin = "[V_TBP].S{3}CCN{0}(CCS{2})CCS{4}.C[Si](C)(C)N{1}"
        self.assertZeroHNitrogen(convert_oin_to_msmiles(oin))

    def test_anilide_nitrogen_is_zero_h(self):
        """N{n}c1ccccc1 -- one heavy neighbour (aryl C). UDIVUY_comp_0."""
        oin = "[Pd_SPL].[Cl]{0}.[Cl]{1}.[Cl]{2}.N{3}c1ccccc1"
        self.assertZeroHNitrogen(convert_oin_to_msmiles(oin))

    def test_azide_binding_nitrogen_is_zero_h(self):
        """N{n}N#N -- one heavy neighbour, single bond. FENMIX_comp_0."""
        oin = "[Mn_OCT].N{0}N#N.N{1}N#N.[Cl]{2}.[Cl]{3}.[Cl]{4}.[Cl]{5}"
        self.assertZeroHNitrogen(convert_oin_to_msmiles(oin))

    def test_phosphinimide_nitrogen_is_zero_h(self):
        """P(=N{n}) -- one heavy neighbour through a double bond. ADEYIR_comp_0."""
        oin = "[Ti_TET].[Cl]{0}.[Cl]{1}.[Cl]{2}.CC(C)(C)P(=N{3})(C(C)(C)C)C(C)(C)C"
        self.assertZeroHNitrogen(convert_oin_to_msmiles(oin))

    def test_amido_nitrogen_is_zero_h(self):
        """CC(C)N{n}C(C)C -- two heavy neighbours (pre-existing behaviour)."""
        oin = "[Ti_OCT].[Cl]{0}.[Cl]{1}.[Cl]{2}.[Cl]{3}.CC(C)N{4}C(C)C.[Cl]{5}"
        self.assertZeroHNitrogen(convert_oin_to_msmiles(oin))


class TestDativeDonorsKeepHydrogens(unittest.TestCase):
    """L-type donors that legitimately carry H must not be stripped."""

    def test_ammine_keeps_three_hydrogens(self):
        """A lone bare N{n} is ammine NH3, not a nitride: it has no heavy neighbour.

        Regression guard for AFAVIO/RIZVAY/OQIHUT/XILBIF, which pass today only
        because the strip rule stops at heavy >= 1.
        """
        oin = "[Ni_OCT].N{0}.N{1}.c1c[nH]cn{2}1.c1c[nH]cn{3}1.[Cl]{4}.[Cl]{5}"
        msmiles = convert_oin_to_msmiles(oin)
        self.assertEqual(
            len(re.findall(r"\[NH3:\d+\]", msmiles)), 2, f"ammine lost its H in: {msmiles}"
        )

    def test_dative_secondary_amine_keeps_its_hydrogen(self):
        """An explicit [NH]{n} bracket is authoritative."""
        oin = "[Ti_OCT].[Cl]{0}.[Cl]{1}.[Cl]{2}.C[NH]{3}C.[Cl]{4}.[Cl]{5}"
        msmiles = convert_oin_to_msmiles(oin)
        self.assertRegex(msmiles, r"C\[NH:\d+\]C", f"dative amine lost its H in: {msmiles}")

    def test_aqua_keeps_two_hydrogens(self):
        """An explicit [OH2]{n} bracket survives the unconditional bare-O strip."""
        oin = "[Pt_SPL].[Cl]{0}.[Cl]{1}.[OH2]{2}.[OH2]{3}"
        msmiles = convert_oin_to_msmiles(oin)
        self.assertEqual(
            len(re.findall(r"\[OH2:\d+\]", msmiles)), 2, f"aqua lost its H in: {msmiles}"
        )

    def test_nhc_carbene_and_alkoxide_stay_zero_h(self):
        """The C and O rules are unchanged by the N gate widening."""
        carbene = convert_oin_to_msmiles(
            "[Pd_SPL].C=CCN1C{0}N(CC=C)c2ccccc21.[Br]{1}.C=CCN1C{2}N(CC=C)c2ccccc21.[Br]{3}"
        )
        self.assertNotIn("[CH2:", carbene, carbene)
        alkoxide = convert_oin_to_msmiles("[Ti_OCT].[Cl]{0}.[Cl]{1}.[Cl]{2}.[Cl]{3}.CCO{4}.CCO{5}")
        self.assertNotIn("[OH:", alkoxide, alkoxide)
        self.assertRegex(alkoxide, r"\[O:\d+\]", alkoxide)


class TestUncoordinatedFragment(unittest.TestCase):
    """A fragment with no binding slot cannot be expressed in m-SMILES."""

    def test_uncoordinated_solvent_raises(self):
        """A free water carries no {slot}; previously an IndexError at frag_vectors[0]."""
        oin = "[Zn_TET].[Cl]{0}.[Cl]{1}.[Cl]{2}.[Cl]{3}.O"
        with self.assertRaises(UncoordinatedFragmentError) as ctx:
            convert_oin_to_msmiles(oin)
        self.assertIn("has no binding slot", str(ctx.exception))

    def test_uncoordinated_counterion_raises(self):
        """An outer-sphere borate anion. BEYHEU_comp_0 class."""
        oin = "[Sc_LIN].[Cl]{0}.[Cl]{1}.CB(c1ccccc1)(c1ccccc1)c1ccccc1"
        with self.assertRaises(UncoordinatedFragmentError):
            convert_oin_to_msmiles(oin)

    def test_error_is_a_value_error(self):
        """Callers catching ValueError keep working."""
        self.assertTrue(issubclass(UncoordinatedFragmentError, ValueError))


if __name__ == "__main__":
    unittest.main()
