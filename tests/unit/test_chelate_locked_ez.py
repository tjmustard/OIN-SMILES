"""A double bond held rigid by a chelate ring must not carry an E/Z marker.

A chelating ligand closes its ring *through the metal*. ``get_oin_string``
fragments the complex and strips the metal, which opens that ring, so a C=C or
C=N that is physically ring-locked looks acyclic to RDKit and gets a directional
marker whose sign falls out of the SMILES traversal order. Re-encoding the
generated 3D structure traverses differently and flips it, failing the round trip
on a bond that never had a free E/Z.

Encoder and generator must agree on exactly which bonds are locked, so both use
the same virtual-metal ring test. These tests pin that agreement.
"""

import os
import unittest

from rdkit import Chem

from oinsmiles import XYZToSMILES
from oinsmiles.core.translator import _clear_chelate_locked_bond_stereo
from oinsmiles.generator3d.ligand import get_ligand_from_smiles

_FIXTURES = os.path.join(os.path.dirname(__file__), "..", "fixtures")


class TestRDKitIgnoresDativeBondsForRings(unittest.TestCase):
    """The gotcha the whole fix rests on: SSSR does not traverse DATIVE bonds.

    This is why a chelate ring is invisible to RDKit even with the metal present,
    and hence why DetectBondStereochemistry marks a ring-locked alkene at all.
    Both ring tests work around it by upgrading DATIVE to SINGLE on a scratch copy.
    """

    @staticmethod
    def _dative_ring():
        rw = Chem.RWMol()
        for symbol in "NCCN":
            rw.AddAtom(Chem.Atom(symbol))
        rw.AddAtom(Chem.Atom("Pt"))
        rw.AddBond(0, 1, Chem.BondType.SINGLE)
        rw.AddBond(1, 2, Chem.BondType.SINGLE)
        rw.AddBond(2, 3, Chem.BondType.SINGLE)
        rw.AddBond(0, 4, Chem.BondType.DATIVE)
        rw.AddBond(3, 4, Chem.BondType.DATIVE)
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        return mol

    def test_dative_ring_is_invisible_to_sssr(self):
        self.assertEqual(len(Chem.GetSymmSSSR(self._dative_ring())), 0)

    def test_upgrading_dative_to_single_reveals_the_ring(self):
        rw = Chem.RWMol(self._dative_ring())
        for bond in rw.GetBonds():
            if bond.GetBondType() == Chem.BondType.DATIVE:
                bond.SetBondType(Chem.BondType.SINGLE)
        mol = rw.GetMol()
        Chem.FastFindRings(mol)
        self.assertEqual(len(Chem.GetSymmSSSR(mol)), 1)


class TestGeneratorFilterMatchesEncoderChelateRing(unittest.TestCase):
    """The generator's double-bond filter now uses the encoder's chelate-ring test.

    ``ligand.py`` used to drop a broader set than the encoder (any bond touching a
    donor or a donor's neighbour), so a free C=N hanging off a monodentate donor
    was suppressed and its genuine E/Z never round-tripped. It now mirrors
    ``core/translator._clear_chelate_locked_bond_stereo`` (S6a): only a bond both of
    whose atoms lie in a ring closing through the metal is dropped. A metal-ring
    imine (AGULIX) stays dropped; a free monodentate arm (AFECIZ) and a pendant
    alkene are now enforced. Enforcing the free arm relies on the charge-aware
    ``_apply_double_bond_stereo`` (see ``TestChargeAwareDoubleBondPromotion``): the
    PuLP-demoted C=N is promoted back to DOUBLE with a +1 formal charge on the N,
    not a 4-valent neutral N that would crash every ff_clean.
    """

    def test_agulix_chelate_imine_is_not_enforced(self):
        # The C=N is ring-locked through the metal (C-[S:5]...metal...[N:3]-N=C);
        # constraining it cut AGULIX's conformer yield 9/9 -> 3/9, so it stays
        # dropped by the chelate-ring test.
        lig = get_ligand_from_smiles(r"CS/C(=N/[N:3]=Cc1ccccc1[O:6])[S:5]")
        self.assertEqual(lig.molecule.stereo_bonds, [])

    def test_pendant_alkene_is_enforced(self):
        lig = get_ligand_from_smiles("C/C=C/CC[NH2:1]")
        self.assertTrue(lig.molecule.stereo_bonds)

    def test_dangling_imine_is_enforced(self):
        # AFECIZ's monodentate salicylaldiminate arm: one donor ([O:2]) closes no
        # metal ring, so the free C=N is now enforced and round-trips. This is the
        # bond the "pending charge fix" guard used to hold back; the charge-aware
        # promotion is what makes enforcing it safe.
        lig = get_ligand_from_smiles(r"Cc1cccc(C)c1/N=C(\[O:2])c1ccccc1")
        self.assertTrue(
            lig.molecule.stereo_bonds,
            "the free monodentate C=N must now be enforced (chelate-ring narrowing)",
        )


class TestChargeAwareDoubleBondPromotion(unittest.TestCase):
    """``_apply_double_bond_stereo`` promotes a demoted C=N with a +1 charge.

    PuLP re-perception can hand the embed mol a carried C=N as a SINGLE bond with a
    neutral 3-valent N (the electrons paid for elsewhere). Restoring the double bond
    over-fills the N; the fix bumps its formal charge to +1 (the encoder's charged
    Lewis form) so the mol stays valence-valid instead of raising in every ff_clean.
    Pre-fix the promotion was declined and the bond left SINGLE, so this fails.
    """

    @staticmethod
    def _demoted_nitrone():
        # C0=N1 written as SINGLE (PuLP's demotion); N1 also bonded to O2 and a
        # methyl C4, so promoting C0-N1 to DOUBLE makes N1 4-valent. C3 / C4 are
        # the stereo reference neighbours.
        rw = Chem.RWMol()
        for symbol in "CNOCC":
            rw.AddAtom(Chem.Atom(symbol))
        rw.AddBond(0, 1, Chem.BondType.SINGLE)  # C=N carried as single
        rw.AddBond(1, 2, Chem.BondType.SINGLE)  # N-O
        rw.AddBond(0, 3, Chem.BondType.SINGLE)  # C-CH3 (ref for C0)
        rw.AddBond(1, 4, Chem.BondType.SINGLE)  # N-CH3 (ref for N1)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        return mol

    def test_promotion_bumps_nitrogen_charge_and_stays_valid(self):
        from oinsmiles.generator3d.embed import _apply_double_bond_stereo

        mol = self._demoted_nitrone()
        stereo_bonds = [(0, 1, Chem.BondStereo.STEREOCIS, 3, 4)]
        _apply_double_bond_stereo(mol, stereo_bonds)

        bond = mol.GetBondBetweenAtoms(0, 1)
        self.assertEqual(bond.GetBondType(), Chem.BondType.DOUBLE)
        self.assertEqual(mol.GetAtomWithIdx(1).GetFormalCharge(), 1)
        # The whole mol now passes a full sanitize (the 4-valent N is now N+).
        Chem.SanitizeMol(Chem.Mol(mol))

    def test_doubly_overfull_bond_is_declined(self):
        # A bond whose promotion over-fills BOTH endpoints (no single +1 bump
        # helps, e.g. FIXYER's C#6 -> valence 5) must degrade to random: the bond
        # stays as re-perceived and no charge is invented.
        from oinsmiles.generator3d.embed import _apply_double_bond_stereo

        mol = Chem.RWMol(Chem.MolFromSmiles("C(C)(C)(C)C"))  # neopentane, C0 valence 4
        mol = mol.GetMol()
        stereo_bonds = [(0, 1, Chem.BondStereo.STEREOCIS, 2, 3)]
        _apply_double_bond_stereo(mol, stereo_bonds)
        self.assertEqual(mol.GetBondBetweenAtoms(0, 1).GetBondType(), Chem.BondType.SINGLE)
        self.assertEqual(mol.GetAtomWithIdx(0).GetFormalCharge(), 0)


class TestEncoderSuppression(unittest.TestCase):
    """End-to-end on the fixture that has always expected a bare C=C."""

    def test_voacac2_emits_no_direction_markers(self):
        oin = XYZToSMILES().convert(os.path.join(_FIXTURES, "VOacac2.xyz"))
        self.assertEqual(oin, "[V_SPY].O{0}.CC(=O{1})C=C(C)O{4}.CC(=O{2})C=C(C)O{3}")
        self.assertNotIn("/", oin)
        self.assertNotIn("\\", oin)

    def test_no_metal_is_a_no_op(self):
        mol = Chem.MolFromSmiles(r"C/C=C/C")
        Chem.DetectBondStereochemistry(mol, -1)
        Chem.AssignStereochemistry(mol, force=True)
        _clear_chelate_locked_bond_stereo(mol)
        stereo = [b.GetStereo() for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]
        self.assertNotEqual(stereo, [Chem.BondStereo.STEREONONE])


if __name__ == "__main__":
    unittest.main()
