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


class TestGeneratorFilterIsBroaderThanTheEncoder(unittest.TestCase):
    """The generator's near-donor proxy drops MORE than the encoder suppresses.

    The encoder's chelate-ring test is the physically correct predicate, and
    narrowing ligand.py to match it is the obvious next step -- but it detonates a
    latent bug. ``_apply_double_bond_stereo`` force-sets a carried bond back to
    DOUBLE without touching formal charges; on AFECIZ the PuLP charge assignment
    wants that C=N single with a charged N, so ``SanitizeMol`` then rejects a
    4-valent N and *every* ff_clean raises. AFECIZ goes from 553s (clean succeeds on
    attempt 1) to exhausting the whole attempt budget at 1.6s per failed clean.

    So the two halves deliberately disagree today: the encoder emits an E/Z on a
    free C=N hanging off a monodentate donor, and the generator does not enforce it.
    AFECIZ and XIZXAG still fail their round trip on that bond, exactly as before
    this change. Fix _apply_double_bond_stereo's charge handling first, then narrow.
    """

    def test_agulix_chelate_imine_is_not_enforced(self):
        # The C=N sits next to the metal-binding N/S donors; constraining it cut
        # AGULIX's conformer yield 9/9 -> 3/9.
        lig = get_ligand_from_smiles(r"CS/C(=N/[N:3]=Cc1ccccc1[O:6])[S:5]")
        self.assertEqual(lig.molecule.stereo_bonds, [])

    def test_pendant_alkene_is_enforced(self):
        lig = get_ligand_from_smiles("C/C=C/CC[NH2:1]")
        self.assertTrue(lig.molecule.stereo_bonds)

    def test_dangling_imine_is_not_enforced_pending_charge_fix(self):
        # AFECIZ's monodentate salicylaldiminate arm. Its E/Z is genuine and the
        # encoder emits it, but enforcing it here breaks ff_clean (see class doc).
        lig = get_ligand_from_smiles(r"Cc1cccc(C)c1/N=C(\[O:2])c1ccccc1")
        self.assertEqual(
            lig.molecule.stereo_bonds,
            [],
            "enforcing this bond breaks ff_clean until _apply_double_bond_stereo "
            "adjusts formal charges when it restores the double bond",
        )


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
