"""sp3 atom chirality must survive the MetalloGen embed.

The ace_mol that MetalloGen embeds carries only atomic number, formal charge and
bond order, so the '@'/'@@' parsed from the m-SMILES used to be discarded before
the embed ever ran. With no chiral constraint and an unseeded random seed, each
sp3 stereocentre landed on a random enantiomer -- independently, per run. These
tests pin the three pieces that fix it: the capture in ligand.py, the
parity-correct re-application in embed.py, and the signed-volume convention the
post-embed verifier relies on.
"""

import tempfile
import unittest

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from oinsmiles import XYZToSMILES
from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generator3d.embed import _apply_atom_chirality, _permutation_is_odd
from oinsmiles.generator3d.ligand import get_ligand_from_smiles

# A square-planar Pt with one chiral carbon on a pendant aminoalcohol arm. The
# stereocentre is three bonds from the donor, so it is genuinely free -- nothing
# but an enforced chiral constraint can make the embed reproduce it.
OIN_R = "[Pt_SPL].C[C@H](O)CN{0}.[Cl]{1}.[Cl]{2}.N{3}"
OIN_S = "[Pt_SPL].C[C@@H](O)CN{0}.[Cl]{1}.[Cl]{2}.N{3}"


def _generated_ligand_fragment(oin):
    """Generate 3D from ``oin``, re-encode it, and return the ligand fragment."""
    result = OIN3DGenerator(engine="metallogen", optimizer="ff").generate(oin)
    with tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False) as handle:
        handle.write(result.xyz)
        path = handle.name
    return XYZToSMILES().convert(path).split(".")[1]


class TestSignedVolumeConvention(unittest.TestCase):
    """Pin RDKit's tag/volume convention that _chiral_targets_satisfied assumes."""

    def test_cw_is_positive_volume(self):
        for smiles in ("C[C@H](N)O", "C[C@@H](N)O", "C[C@H](N)C(=O)O"):
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            params = AllChem.ETKDGv3()
            params.randomSeed = 7
            AllChem.EmbedMolecule(mol, params)
            positions = mol.GetConformer().GetPositions()
            for atom in mol.GetAtoms():
                tag = atom.GetChiralTag()
                if tag == Chem.ChiralType.CHI_UNSPECIFIED:
                    continue
                nbrs = [b.GetOtherAtomIdx(atom.GetIdx()) for b in atom.GetBonds()]
                if len(nbrs) != 4:
                    continue
                p0, p1, p2, p3 = (np.array(positions[i]) for i in nbrs)
                volume = float(np.dot(p1 - p0, np.cross(p2 - p0, p3 - p0)))
                expected_positive = tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW
                self.assertEqual(
                    expected_positive,
                    volume > 0,
                    f"{smiles}: {tag} should have {'positive' if expected_positive else 'negative'}"
                    f" signed volume, got {volume:+.3f}",
                )


class TestChiralCentreCapture(unittest.TestCase):
    """ligand.py must record the stereocentre, anchored to a neighbour ordering."""

    def test_pendant_stereocentre_is_captured(self):
        lig = get_ligand_from_smiles("C[C@H](O)C[NH2:1]")
        centers = lig.molecule.chiral_centers
        self.assertEqual(len(centers), 1, "the single sp3 stereocentre must be captured")
        _center, nbrs, tag = centers[0]
        self.assertEqual(len(nbrs), 4, "the neighbour tuple anchors the tag's parity")
        self.assertIn(
            tag, (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        )

    def test_achiral_ligand_captures_nothing(self):
        lig = get_ligand_from_smiles("CC(O)C[NH2:1]")
        self.assertEqual(lig.molecule.chiral_centers, [])

    def test_metal_binding_centre_is_skipped(self):
        # A donor atom gains a metal bond in the complex, so its neighbour set --
        # and hence the meaning of its tag -- changes. Zone-A donor stereo is
        # re-asserted from the template by build_contract_mol instead.
        lig = get_ligand_from_smiles("C[C@H](O)[NH:1]C")
        self.assertEqual(
            [c for c in lig.molecule.chiral_centers if c[0] == 3],
            [],
            "a metal-binding stereocentre must not be carried into the embed",
        )


class TestApplyAtomChirality(unittest.TestCase):
    """Re-applying a tag must respect the neighbour ordering it was read against."""

    def test_odd_permutation_flips_the_tag(self):
        self.assertFalse(_permutation_is_odd((0, 1, 2, 3), (0, 1, 2, 3)))
        self.assertTrue(_permutation_is_odd((0, 1, 2, 3), (1, 0, 2, 3)))
        self.assertFalse(_permutation_is_odd((0, 1, 2, 3), (1, 0, 3, 2)))

    def test_reapplied_tag_preserves_configuration(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](N)O"))
        center = 1
        atom = mol.GetAtomWithIdx(center)
        nbrs = tuple(b.GetOtherAtomIdx(center) for b in atom.GetBonds())
        tag = atom.GetChiralTag()
        original = Chem.MolToSmiles(mol)

        # Wipe the tag, then restore it through the same path the embed uses but
        # with the neighbour tuple deliberately given in a swapped (odd) order.
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        swapped = (nbrs[1], nbrs[0], nbrs[2], nbrs[3])
        flipped = (
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW
            if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW
            else Chem.ChiralType.CHI_TETRAHEDRAL_CW
        )
        _apply_atom_chirality(mol, [(center, swapped, flipped)])
        self.assertEqual(
            Chem.MolToSmiles(mol),
            original,
            "an odd neighbour permutation with an inverted tag denotes the same 3D centre",
        )

    def test_changed_neighbour_set_is_skipped(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CC(N)O"))
        # Neighbour indices that do not match the mol: must be ignored, not raise.
        _apply_atom_chirality(mol, [(1, (0, 2, 3, 99), Chem.ChiralType.CHI_TETRAHEDRAL_CW)])
        self.assertEqual(
            mol.GetAtomWithIdx(1).GetChiralTag(),
            Chem.ChiralType.CHI_UNSPECIFIED,
        )


class TestGeneratorRoundTripsChirality(unittest.TestCase):
    """The real generator must return the requested enantiomer, reproducibly.

    On pristine main both of these failed: the '@' case flipped on one run in
    three, and the '@@' case came back as its mirror image every time.
    """

    def test_r_enantiomer_round_trips(self):
        self.assertEqual(_generated_ligand_fragment(OIN_R), "C[C@H](O)CN{0}")

    def test_s_enantiomer_round_trips(self):
        self.assertEqual(_generated_ligand_fragment(OIN_S), "C[C@@H](O)CN{0}")

    def test_generation_is_deterministic(self):
        runs = {_generated_ligand_fragment(OIN_R) for _ in range(3)}
        self.assertEqual(len(runs), 1, f"the seeded embed must be reproducible, got {sorted(runs)}")


if __name__ == "__main__":
    unittest.main()
