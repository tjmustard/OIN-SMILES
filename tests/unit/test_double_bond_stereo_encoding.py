"""Guard that the OIN encoder preserves C=C (cis/trans, E/Z) double-bond stereo.

`get_oin_string` rebuilds each ligand fragment atom-by-atom; that copied bond
TYPE but not the STEREOE/Z reference atoms, so cis/trans alkene geometry was
silently dropped (both E and Z re-encoded as a bare ``C=C``). These tests feed a
metal complex whose alkenyl ligand has a KNOWN E or Z 3D geometry (a clean RDKit
embed, independent of the 3D generator) and assert the encoder emits the matching
directional markers.
"""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from oinsmiles.utils.xyz2mol import get_oin_string


def _encode_alkenyl_amine_tmc(lig_smiles):
    """Embed a pent-3-enylamine ligand (defined E/Z), bolt it onto PtCl3, and
    return the OIN string produced by the encoder for that 3D structure."""
    lig = Chem.AddHs(Chem.MolFromSmiles(lig_smiles))
    AllChem.EmbedMolecule(lig, randomSeed=3)
    AllChem.MMFFOptimizeMolecule(lig)
    lig = Chem.RemoveHs(lig)

    rw = Chem.RWMol(lig)
    conf = rw.GetConformer()
    n_idx = next(a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 7)
    n_pos = np.array(conf.GetAtomPosition(n_idx))

    pt = rw.AddAtom(Chem.Atom(78))
    conf.SetAtomPosition(pt, Point3D(*(n_pos + np.array([2.1, 0.0, 0.0]))))
    rw.AddBond(n_idx, pt, Chem.BondType.DATIVE)
    for off in ([2.1, 2.0, 0.0], [2.1, -2.0, 0.0], [4.2, 0.0, 0.0]):
        cl = rw.AddAtom(Chem.Atom(17))
        conf.SetAtomPosition(cl, Point3D(*(n_pos + np.array(off))))
        rw.AddBond(pt, cl, Chem.BondType.DATIVE)

    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.AssignStereochemistryFrom3D(mol)
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    return get_oin_string(mol, coords)


def _alkene_stereo(oin_string):
    """Perceived C=C stereo of the alkenyl ligand fragment in an OIN string."""
    import re

    for frag in oin_string.split("."):
        clean = re.sub(r"\{[^}]*\}", "", frag)
        if "C=C" not in clean:
            continue
        m = Chem.MolFromSmiles(clean)
        if m is None:
            continue
        Chem.AssignStereochemistry(m, cleanIt=True, force=True)
        for b in m.GetBonds():
            if (
                b.GetBondType() == Chem.BondType.DOUBLE
                and b.GetStereo() != Chem.BondStereo.STEREONONE
            ):
                return b.GetStereo()
    return Chem.BondStereo.STEREONONE


class TestDoubleBondStereoEncoding(unittest.TestCase):
    def test_e_alkene_encoded_with_direction(self):
        oin = _encode_alkenyl_amine_tmc("C/C=C/CCN")
        self.assertIn("/", oin, f"no directional bond marker emitted: {oin}")
        self.assertNotEqual(_alkene_stereo(oin), Chem.BondStereo.STEREONONE, oin)

    def test_z_alkene_encoded_with_direction(self):
        oin = _encode_alkenyl_amine_tmc(r"C/C=C\CCN")
        self.assertIn("/", oin, f"no directional bond marker emitted: {oin}")
        self.assertNotEqual(_alkene_stereo(oin), Chem.BondStereo.STEREONONE, oin)

    def test_e_and_z_encode_differently(self):
        """The whole point: E and Z geometry must not collapse to the same OIN."""
        e = _alkene_stereo(_encode_alkenyl_amine_tmc("C/C=C/CCN"))
        z = _alkene_stereo(_encode_alkenyl_amine_tmc(r"C/C=C\CCN"))
        self.assertNotEqual(e, Chem.BondStereo.STEREONONE)
        self.assertNotEqual(z, Chem.BondStereo.STEREONONE)
        self.assertNotEqual(e, z, "E and Z alkene geometry encoded to the same stereo")


if __name__ == "__main__":
    unittest.main()
