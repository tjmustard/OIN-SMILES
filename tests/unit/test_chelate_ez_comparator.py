"""The round-trip comparator must treat metal-chelate-locked C=C/C=N E/Z as
notation, not structure — but still catch a genuinely free (pendant) E/Z flip.

Regression guard for the v0.3.6 chelate-locked E/Z asymmetry: S6's encoder
(`_clear_chelate_locked_bond_stereo`) drops the directional marker on a double
bond a metal ring holds rigid, but only where the *input* structure's donor
bonds make that ring perceivable. A generated structure whose donor is bonded
differently keeps the marker, so `smiles_1` (no slash) and `smiles_2` (slash)
described the same chelate yet failed the string gate — 62 dataset rows. The
comparator now reconstructs the chelate ring (dummy metal bonded to each slot
atom) and clears E/Z on the locked bonds on both sides, so they key equal.
"""

import unittest

from oinsmiles.oin.compare import canonical_roundtrip_key as key


class TestChelateLockedEZComparator(unittest.TestCase):
    def test_chelate_locked_cn_slash_asymmetry_keys_equal(self):
        """C=N in a ring closed through the metal: input cleared, generated
        slashed — same structure, must compare equal."""
        # O{0} and N{1} both bind the metal, closing a 5-ring over the C=N.
        cleared = "[Ni_SPL].O{0}C(=NN{1}=Cc1ccccc1)c1ccccc1"
        slashed = "[Ni_SPL].O{0}/C(=N\\N{1}=Cc1ccccc1)c1ccccc1"
        self.assertEqual(key(cleared), key(slashed))

    def test_pendant_ez_flip_still_distinguished(self):
        """A freely-rotatable C=C in no metal ring keeps its E/Z: an E vs Z
        flip is a real diastereomer and must NOT compare equal."""
        e = "[Pt_SPL].C/C=C/CN{0}.N{1}.[Cl]{2}.[Cl]{3}"
        z = "[Pt_SPL].C/C=C\\CN{0}.N{1}.[Cl]{2}.[Cl]{3}"
        self.assertNotEqual(key(e), key(z))

    def test_chelate_normalization_does_not_collapse_different_skeletons(self):
        """Clearing chelate E/Z must not make two different ligands compare
        equal — only the locked double bond's direction is normalized."""
        salicyl = "[Ni_SPL].O{0}C(=NN{1}=Cc1ccccc1)c1ccccc1"
        pyridyl = "[Ni_SPL].O{0}C(=NN{1}=Cc1ccncc1)c1ccccc1"
        self.assertNotEqual(key(salicyl), key(pyridyl))

    def test_non_chelate_fragments_unchanged(self):
        """A ligand with no metal-locked double bond keys exactly as before
        (the fix is a no-op outside chelate rings)."""
        a = "[Pt_SPL].N{0}.N{1}.[Cl]{2}.[Cl]{3}"
        b = "[Pt_SPL].N{0}.N{1}.[Cl]{2}.[Cl]{3}"
        self.assertEqual(key(a), key(b))

    # -- The chelate-lock clear must reach fragments that only PARTIALLY sanitize.
    #    A slot-stripped eta-Cp/Cp* ring or a bare-`n` azole/pyridyl donor ring raises
    #    KekulizeException on a full sanitize; before the fix the comparator bailed to
    #    the `RAW:` fallback, which preserves the raw `/`\ verbatim, so exactly the
    #    fragments carrying the ring-locked slash never got cleared (62-row v0.3.6
    #    regression; these 10 were the residual after the first hotfix). --

    def test_unkekulizable_dipyrrin_chelate_keys_equal(self):
        """CUHSEE: the meso C=C of a dipyrrin sits in the N{0}...metal...n{1}
        macrocycle. The fragment won't kekulize (bare-`n` azole ring), but the
        ring-locked slash must still normalize away."""
        cleared = (
            "[Rh_TET].N#Cc1ccc(C(=C2C=CC=N{0}2)c2cccn{1}2)cc1."
            "Cc{2}1c{2>}(C)c{2}(C)c{2}(C)c{2}1C.[Cl]{3}"
        )
        slashed = (
            "[Rh_TET].N#Cc1ccc(/C(=C2\\C=CC=N{0}2)c2cccn{1}2)cc1."
            "Cc{2}1c{2>}(C)c{2}(C)c{2}(C)c{2}1C.[Cl]{3}"
        )
        self.assertEqual(key(cleared), key(slashed))

    def test_unkekulizable_azole_donor_ring_keys_equal(self):
        """KAQLIZ: azomethine C=N in a Ti chelate whose pyridyl donor ring
        (bare `n`) won't kekulize once the metal is stripped."""
        cleared = (
            "[Ti_OCT].C1CCO{0}C1.[Cl]{1}.CC(C)(C)c1cc(C2=N{2}C(=C(c3ccccc3)"
            "c3ccc(-c4cc(C(C)(C)C)cc(C(C)(C)C)c4O{3})n{4}3)C=C2)c(O{5})c(C(C)(C)C)c1"
        )
        slashed = (
            "[Ti_OCT].C1CCO{0}C1.[Cl]{1}.CC(C)(C)c1cc(C2=N{2}/C(=C(/c3ccccc3)"
            "c3ccc(-c4cc(C(C)(C)C)cc(C(C)(C)C)c4O{3})n{4}3)C=C2)c(O{5})c(C(C)(C)C)c1"
        )
        self.assertEqual(key(cleared), key(slashed))

    def test_genuine_ez_flip_direction_stays_distinct(self):
        """RIQFON: a real diastereomer flip (same slash COUNT, opposite
        DIRECTION `/N=C/` vs `/N=C\\`) on an EXOCYCLIC imine -- the ring-locked
        predicate must leave it alone, so the two forms stay distinct."""
        one = (
            "[Ni_SPL].CC(=N{0}c1c(C(C)C)cccc1C(C)C)/C(O{1})=N/c1c(C(C)C)cccc1C(C)C."
            "c1ccn{2}cc1.O=C{3}c1ccccc1"
        )
        two = (
            "[Ni_SPL].CC(=N{0}c1c(C(C)C)cccc1C(C)C)/C(O{1})=N\\c1c(C(C)C)cccc1C(C)C."
            "c1ccn{2}cc1.O=C{3}c1ccccc1"
        )
        self.assertNotEqual(key(one), key(two))


if __name__ == "__main__":
    unittest.main()
