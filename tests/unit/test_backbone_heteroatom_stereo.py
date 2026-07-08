"""Guards for backbone P/S/Si stereocentre recovery in the OIN re-encode.

When a generated 3D structure is re-encoded, ``ChiralityRecoveryUtility.recover``
runs on each ligand fragment. It is P/N-aware: a 4-coordinate backbone P with a
valid chiral tag but no encoded CIP is treated as a stray and cleared. So the
MetalloGen adapter stamps the template CIP as ``_OIN_CIPCode`` on a backbone P
(see build_contract_mol), which makes recover() KEEP and orient the tag.

Si and S stereocentres are not P/N, so recover() leaves them alone -- they only
need the adapter's perceive-then-flip loop to match the encoded handedness.

These tests pin that contract at the recover() level (no 3D generation needed).
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem

from oinsmiles.core.chirality import ChiralityRecoveryUtility


def _recovered_smiles(mol):
    return Chem.MolToSmiles(ChiralityRecoveryUtility().recover(mol))


class TestBackboneHeteroatomStereo(unittest.TestCase):
    def _iminophosphorane(self, stamp_cip=None):
        # 4-coordinate P stereocentre: =N-Me, Me, Et, Ph (four distinct groups).
        mol = Chem.MolFromSmiles("CN=[P@](C)(CC)c1ccccc1")
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        p = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 15)
        self.assertEqual(p.GetTotalDegree(), 4)
        if stamp_cip is not None:
            p.SetProp("_OIN_CIPCode", stamp_cip)
        return mol

    def test_backbone_p_kept_when_oin_cipcode_stamped(self):
        """A backbone P survives recover() when the encoded CIP is stamped."""
        out = _recovered_smiles(self._iminophosphorane(stamp_cip="R"))
        self.assertIn("@", out, f"backbone P stereo dropped: {out}")

    def test_backbone_p_cleared_without_oin_cipcode(self):
        """Without the stamp, recover() clears the P tag (the behavior the adapter
        works around by stamping _OIN_CIPCode)."""
        out = _recovered_smiles(self._iminophosphorane(stamp_cip=None))
        self.assertNotIn("@", out, f"expected P stereo cleared without stamp: {out}")

    def test_silicon_stereo_survives_recover(self):
        """Si is not P/N, so recover() must leave its tag untouched."""
        mol = Chem.MolFromSmiles("C[Si@](CC)(CCC)c1ccccc1")
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        self.assertIn("@", _recovered_smiles(mol))

    def test_sulfur_stereo_survives_recover(self):
        """A chiral sulfoxide S is not P/N, so recover() must leave its tag."""
        mol = Chem.MolFromSmiles("C[S@](=O)CC")
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        self.assertIn("@", _recovered_smiles(mol))


if __name__ == "__main__":
    unittest.main()
