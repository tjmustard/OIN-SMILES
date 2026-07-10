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


if __name__ == "__main__":
    unittest.main()
