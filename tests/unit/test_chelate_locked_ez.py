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
from oinsmiles.generator3d.ligand import _chelate_locked_atoms, get_ligand_from_smiles

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


class TestChelateLockedPredicate(unittest.TestCase):
    """The generator's virtual-metal ring test, on the ligands that motivated it."""

    @staticmethod
    def _locked_stereo_bonds(smiles):
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        Chem.SetBondStereoFromDirections(mol)
        locked = _chelate_locked_atoms(mol)
        return [
            (b.GetBeginAtomIdx() in locked and b.GetEndAtomIdx() in locked)
            for b in mol.GetBonds()
            if b.GetBondType() == Chem.BondType.DOUBLE
            and b.GetStereo() != Chem.BondStereo.STEREONONE
        ]

    def test_acac_alkene_is_locked(self):
        self.assertEqual(self._locked_stereo_bonds(r"CC(=[O:1])/C=C(/C)[O:2]"), [True])

    def test_agulix_chelate_imine_is_locked(self):
        # The C=N sits on the S->N chelate path; constraining it cut AGULIX's
        # conformer yield 9/9 -> 3/9, which is what the filter exists to prevent.
        self.assertEqual(self._locked_stereo_bonds(r"CS/C(=N/[N:3]=Cc1ccccc1[O:6])[S:5]"), [True])

    def test_pendant_alkene_is_free(self):
        self.assertEqual(self._locked_stereo_bonds("C/C=C/CC[NH2:1]"), [False])

    def test_dangling_imine_of_monodentate_arm_is_free(self):
        # AFECIZ's third salicylaldiminate arm binds through O only; its imine N
        # is a leaf, so the C=N lies on no chelate ring and its E/Z is genuine.
        # The old "touches a donor or a donor's neighbour" proxy wrongly dropped it.
        self.assertEqual(self._locked_stereo_bonds(r"Cc1cccc(C)c1/N=C(\[O:2])c1ccccc1"), [False])

    def test_monodentate_ligand_has_no_locked_atoms(self):
        mol = Chem.MolFromSmiles("C/C=C/CC[NH2:1]", sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        self.assertEqual(_chelate_locked_atoms(mol), set())


class TestGeneratorFilterAgrees(unittest.TestCase):
    """ligand.py must enforce exactly the bonds the encoder still emits."""

    def test_agulix_imine_is_not_enforced(self):
        lig = get_ligand_from_smiles(r"CS/C(=N/[N:3]=Cc1ccccc1[O:6])[S:5]")
        self.assertEqual(lig.molecule.stereo_bonds, [])

    def test_pendant_alkene_is_enforced(self):
        lig = get_ligand_from_smiles("C/C=C/CC[NH2:1]")
        self.assertTrue(lig.molecule.stereo_bonds)

    def test_dangling_imine_is_enforced(self):
        lig = get_ligand_from_smiles(r"Cc1cccc(C)c1/N=C(\[O:2])c1ccccc1")
        self.assertTrue(
            lig.molecule.stereo_bonds,
            "a free C=N hanging off a monodentate donor must keep its constraint",
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
