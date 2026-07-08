"""Guards for Zone-A P *donor* lone-pair stereocentre recovery in the re-encode.

A Zone-A P is a phosphorus that binds the metal directly; after the OINSanitizer
strips the metal, the ligand fragment carries a 3-coordinate P whose stereogenic
lone pair is the 4th "substituent". ``ChiralityRecoveryUtility.recover`` has a
dedicated lone-pair branch that KEEPS and orients such a P *iff* it carries the
``_OIN_CIPCode_LP`` property (rdCIPLabeler convention) plus a specified chiral
tag; without the property it falls through to the degree-keyed branch and is
cleared (``total_degree < 4``).

The MetalloGen embed is stereo-blind for lone-pair chirality, so the generated
geometry has random handedness at the donor P. ``build_contract_mol`` therefore
re-asserts the encoded handedness the same way the forward CIPAssigner does:
stamp ``_OIN_CIPCode_LP`` (taken from rdCIPLabeler on the OIN template, NOT the
legacy ``_CIPCode`` -- the two disagree for a 3-coordinate P) and seed a tag, so
recover() flips to match the stored label regardless of the geometry.

These tests pin that contract at the recover() level (no 3D generation needed),
mirroring tests/unit/test_backbone_heteroatom_stereo.py.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from oinsmiles.core.chirality import _LP_CIP_PROP, ChiralityRecoveryUtility


def _p_rdcip_label(mol):
    """rdCIPLabeler 'R'/'S' for the (single) phosphorus in mol, or None."""
    m = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(m)
    p = next(a for a in m.GetAtoms() if a.GetAtomicNum() == 15)
    return p.GetProp("_CIPCode") if p.HasProp("_CIPCode") else None


class TestZoneAPDonorStereo(unittest.TestCase):
    def _chiral_phosphine(self, lp_cip=None, seed_tag=None):
        """A metal-free 3-coordinate chiral phosphine (Me, Et, Ph, lone pair).

        Stands in for a Zone-A P donor fragment after the metal is stripped.
        Optionally stamp ``_OIN_CIPCode_LP`` and seed a specified chiral tag,
        exactly as ``build_contract_mol`` does for the generated contract mol.
        """
        mol = Chem.MolFromSmiles("[P@](C)(CC)c1ccccc1")
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        p = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 15)
        self.assertEqual(p.GetTotalDegree(), 3)  # 3 real substituents + lone pair
        if lp_cip is not None:
            p.SetProp(_LP_CIP_PROP, lp_cip)
        if seed_tag is not None:
            p.SetChiralTag(seed_tag)
        return mol

    def test_zone_a_p_kept_when_lp_cipcode_stamped(self):
        """A Zone-A P donor survives recover() when _OIN_CIPCode_LP is stamped."""
        base = self._chiral_phosphine()
        label = _p_rdcip_label(base)
        self.assertIn(label, ("R", "S"))
        mol = self._chiral_phosphine(lp_cip=label, seed_tag=Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        out = ChiralityRecoveryUtility().recover(mol)
        self.assertIn("@", Chem.MolToSmiles(out), "Zone-A P donor stereo dropped")
        self.assertEqual(_p_rdcip_label(out), label, "recovered enantiomer != stored label")

    def test_zone_a_p_cleared_without_lp_cipcode(self):
        """Without the stamp, recover() clears the 3-coordinate P tag -- the exact
        bug the adapter fix works around by stamping _OIN_CIPCode_LP."""
        mol = self._chiral_phosphine(lp_cip=None, seed_tag=Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        out = Chem.MolToSmiles(ChiralityRecoveryUtility().recover(mol))
        self.assertNotIn("@", out, f"expected Zone-A P stereo cleared without stamp: {out}")

    def test_zone_a_p_reoriented_to_stored_label(self):
        """recover() orients the tag to the STORED label, not the seeded handedness:
        seed one tag but stamp the opposite label -> output matches the label. This
        is what lets build_contract_mol seed an arbitrary tag over the stereo-blind
        embed geometry and still round-trip the encoded enantiomer."""
        base_label = _p_rdcip_label(self._chiral_phosphine())
        opposite = "S" if base_label == "R" else "R"
        # Seed the tag that (bare) yields base_label, but stamp the opposite label.
        mol = self._chiral_phosphine(lp_cip=opposite, seed_tag=Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        out = ChiralityRecoveryUtility().recover(mol)
        self.assertIn("@", Chem.MolToSmiles(out))
        self.assertEqual(
            _p_rdcip_label(out), opposite, "recover() did not flip to the stored label"
        )


if __name__ == "__main__":
    unittest.main()
