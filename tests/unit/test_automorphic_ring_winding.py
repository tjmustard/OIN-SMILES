"""Guards for eta-ring winding canonicalization (v0.3.6, S4).

An eta ring's winding marker ({n>} / {n<}) records which face of the ring the
metal sees. That is only structural information when the ring cannot be turned
over onto itself. Turning a ring over is a 180 deg rotation about an in-plane
axis: a PROPER rotation that leaves the metal and every other ligand alone and
relabels the ring's atoms in reverse cyclic order. So the winding is meaningless
exactly when some automorphism of the ligand fragment reverses the ring's cyclic
order -- and it is load-bearing otherwise (an ansa-bis(indenyl)'s rac vs meso).

These tests pin both halves. The bearing half is the regression floor: if it
goes green-to-red, the encoder has started throwing away real stereochemistry.
"""

import os
import unittest

from rdkit import Chem

from oinsmiles import XYZToSMILES
from oinsmiles.utils.oin_aligner import (
    OINSanitizer,
    _eta_traversal_order,
    _orientation_symmetry_graph,
    _winding_is_orientation_free,
)

_FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "..", "fixtures")


def _ring_atoms(smiles, size, which=0):
    """Indices of the `which`-th SSSR ring of `size` atoms, in SMILES atom order."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    rings = [tuple(sorted(r)) for r in Chem.GetSymmSSSR(mol) if len(r) == size]
    return rings[which]


# Fragment SMILES exactly as the encoder emits them (marker-stripped).
_CP = "[cH]1[cH][cH][cH][cH]1"
_CP_STAR = "Cc1c(C)c(C)c(C)c1C"
_BENZENE = "[cH]1[cH][cH][cH][cH][cH]1"
_MESITYLENE = "Cc1[cH]c(C)[cH]c(C)[cH]1"
_BPH4 = "c1ccc(B(c2[cH][cH][cH][cH][cH]2)(c2ccccc2)c2ccccc2)cc1"
_ANSA_BIS_CP = "C[Si](C)(c1cccc1)c1cccc1"
_CGC_CP = "CC(C)(C)N[Si](C)(C)c1cccc1"
_ANSA_BIS_INDENYL = "C[Si](C)(c1[cH][cH]c2ccccc12)c1[cH][cH]c2ccccc12"
_ARM_CP_STAR = "Cc1c(C)c(C)c(CB(c2cc(F)c(F)c(F)c2F)(c2c(F)c(F)c(F)c(F)c2F)c2c(F)c(F)c(F)c(F)c2F)c1C"
_QOFTOU = "Cc1[cH]c2c(-c3ccccc3)c(C)sc2c1[Si](C)(C)c1c(C)[cH]c2c(-c3ccccc3)c(C)sc12"
_HALIDE_FACE_CP = "Oc1[cH]c(Cl)c(Br)c1I"
_ASYMMETRIC_ALLYL = "COC(=O)C(=[CH][CH]c1ccccc1)c1ccc(Br)cc1"


class TestOrientationFreeRings(unittest.TestCase):
    """Rings that CAN be turned over onto themselves: winding is notation."""

    def test_unsubstituted_rings_are_free(self):
        self.assertTrue(_winding_is_orientation_free(_CP, (0, 1, 2, 3, 4)))
        self.assertTrue(_winding_is_orientation_free(_BENZENE, (0, 1, 2, 3, 4, 5)))

    def test_fully_substituted_cp_star_is_free(self):
        """CAHKEE. Every ring atom bears a methyl -- one automorphism class."""
        self.assertTrue(_winding_is_orientation_free(_CP_STAR, (1, 2, 4, 6, 8)))

    def test_mesitylene_is_free_despite_two_symmetry_classes(self):
        """GELBAD. Its ring atoms are NOT all one symmetry class, yet a reflection
        through any methyl-bearing carbon reverses the cycle. A criterion based on
        `CanonicalRankAtoms(breakTies=False)` transitivity would wrongly call this
        load-bearing and leave the round trip broken."""
        self.assertTrue(_winding_is_orientation_free(_MESITYLENE, (1, 2, 3, 5, 6, 8)))

    def test_arm_substituted_cp_star_is_free(self):
        """NEFNER slot 2. Four methyls plus one -CH2B(ArF)3 arm: four symmetry
        classes, but the reflection through the arm-bearing carbon reverses it."""
        self.assertTrue(_winding_is_orientation_free(_ARM_CP_STAR, _ring_atoms(_ARM_CP_STAR, 5)))

    def test_borate_phenyl_is_free(self):
        """SOJMIQ. The ipso carbon is distinct from ortho/meta/para, yet the C2
        through ipso and para reverses the ring. Also exercises the valence-tolerant
        graph: this fragment's 4-coordinate neutral boron defeats SanitizeMol."""
        self.assertTrue(_winding_is_orientation_free(_BPH4, (5, 6, 7, 8, 9, 10)))

    def test_ansa_bis_cp_rings_are_free(self):
        """TiCat1. Each unsubstituted Cp turns over about the axis through its own
        bridgehead carbon; the bridge and the sibling ring never move. Hence both
        rings read '>' and there is no rac/meso to lose."""
        for which in (0, 1):
            atoms = _ring_atoms(_ANSA_BIS_CP, 5, which)
            self.assertTrue(_winding_is_orientation_free(_ANSA_BIS_CP, atoms), atoms)

    def test_constrained_geometry_cp_is_free(self):
        """TiCat2. The amide donor sits off the ring and does not move with it."""
        self.assertTrue(_winding_is_orientation_free(_CGC_CP, _ring_atoms(_CGC_CP, 5)))

    def test_eta2_is_degenerate(self):
        """Two atoms have no circulation; signed_circulation already returns '>'."""
        self.assertTrue(_winding_is_orientation_free("[CH2]=[CH2]", (0, 1)))


class TestOrientationBearingRings(unittest.TestCase):
    """Rings that CANNOT be turned over: winding is real stereochemistry.

    This class is the regression floor for the 2026-07-06 eta-winding fix.
    """

    def test_ansa_bis_indenyl_rings_are_bearing(self):
        """TiCat3/TiCat4. Reversing the five-ring through its bridgehead would have
        to swap a CH for a ring-fusion carbon. No automorphism does that, so the
        winding distinguishes rac from meso and must survive."""
        for which in (0, 1):
            atoms = _ring_atoms(_ANSA_BIS_INDENYL, 5, which)
            self.assertFalse(_winding_is_orientation_free(_ANSA_BIS_INDENYL, atoms), atoms)

    def test_qoftou_ansa_zirconocene_is_bearing(self):
        self.assertFalse(_winding_is_orientation_free(_QOFTOU, _ring_atoms(_QOFTOU, 5)))

    def test_four_different_substituents_is_bearing(self):
        """The Ferrocene-halide-face / ChiralityWitnessRing fixtures rely on this."""
        atoms = _ring_atoms(_HALIDE_FACE_CP, 5)
        self.assertFalse(_winding_is_orientation_free(_HALIDE_FACE_CP, atoms))

    def test_asymmetric_eta3_allyl_path_is_bearing(self):
        """XIVJEU slot 1. An open chain, not a ring: its two termini differ, so the
        coordinated face is a genuine diastereomeric choice."""
        self.assertFalse(_winding_is_orientation_free(_ASYMMETRIC_ALLYL, (4, 5, 6)))


class TestSymmetryGraph(unittest.TestCase):
    """The reduced graph must forget bond orders but never hydrogen counts."""

    def test_h_count_survives_bond_flattening(self):
        """Cp- and cyclopentadiene have the same skeleton and the same flattened
        bonds. Only the sp3 CH2's extra hydrogen separates them, so the graph must
        carry H counts or the two ligands would be conflated."""
        cp = _orientation_symmetry_graph(Chem.MolFromSmiles(_CP, sanitize=False))
        cpd = _orientation_symmetry_graph(
            Chem.MolFromSmiles("[CH2]1[CH]=[CH][CH]=[CH]1", sanitize=False)
        )
        self.assertNotEqual(Chem.MolToSmiles(cp), Chem.MolToSmiles(cpd))

    def test_kekule_alternation_does_not_break_ring_symmetry(self):
        """A Kekule benzene's alternating bonds must not fake an asymmetry."""
        kekule = "C1=CC=CC=C1"
        self.assertTrue(_winding_is_orientation_free(kekule, (0, 1, 2, 3, 4, 5)))

    def test_traversal_rejects_non_orientable_groups(self):
        """Two disjoint eta rings on one fragment, and a ring plus a pendant donor,
        are not single orientable groups. Both must be refused rather than guessed."""
        mol = Chem.MolFromSmiles(_ANSA_BIS_CP, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        both = sorted(_ring_atoms(_ANSA_BIS_CP, 5, 0) + _ring_atoms(_ANSA_BIS_CP, 5, 1))
        self.assertEqual(_eta_traversal_order(mol, both), (None, None))

        cgc = Chem.MolFromSmiles(_CGC_CP, sanitize=False)
        cgc.UpdatePropertyCache(strict=False)
        ring_plus_donor = sorted(_ring_atoms(_CGC_CP, 5)) + [4]  # +N
        self.assertEqual(_eta_traversal_order(cgc, ring_plus_donor), (None, None))

    def test_undecidable_fragment_falls_back(self):
        """An unparseable fragment returns None so the caller keeps geometric winding."""
        self.assertIsNone(_winding_is_orientation_free("not a smiles", (0, 1, 2)))
        self.assertIsNone(_winding_is_orientation_free(_CP, (0, 1, 99)))


class TestCanonicalEtaSetRepresentative(unittest.TestCase):
    """Which of several equivalent rings carries the eta slot (SOJMIQ)."""

    @staticmethod
    def _mol(smiles):
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        return mol

    def test_all_four_borate_phenyls_converge_on_one_ring(self):
        """Whichever phenyl physically faces the metal, the marker lands on the same
        canonical ring -- so an input structure and its regenerated twin agree."""
        mol = self._mol("c1ccc(B(c2ccccc2)(c2ccccc2)c2ccccc2)cc1")
        targets = set()
        for ring in Chem.GetSymmSSSR(mol):
            ring = sorted(ring)
            mapping = OINSanitizer.canonical_eta_set_representative(mol, ring)
            targets.add(tuple(sorted(mapping.values())) if mapping else tuple(ring))
        self.assertEqual(len(targets), 1, targets)

    def test_mapping_is_a_bijection_onto_the_chosen_ring(self):
        mol = self._mol("c1ccc(B(c2ccccc2)(c2ccccc2)c2ccccc2)cc1")
        ring = sorted([r for r in Chem.GetSymmSSSR(mol) if 5 in r][0])
        mapping = OINSanitizer.canonical_eta_set_representative(mol, ring)
        self.assertTrue(mapping)
        self.assertEqual(sorted(mapping.keys()), ring)
        self.assertEqual(len(set(mapping.values())), len(ring))

    def test_lone_ring_is_untouched(self):
        """No sibling to canonicalize against -- must stay byte-identical."""
        mol = self._mol(_CP)
        self.assertEqual(OINSanitizer.canonical_eta_set_representative(mol, [0, 1, 2, 3, 4]), {})

    def test_ansa_metallocene_is_untouched(self):
        """Its two eta rings are one fragment but not one ring; relabeling either
        onto the other would destroy the per-ring slot assignment."""
        mol = self._mol(_ANSA_BIS_CP)
        both = sorted(_ring_atoms(_ANSA_BIS_CP, 5, 0) + _ring_atoms(_ANSA_BIS_CP, 5, 1))
        self.assertEqual(OINSanitizer.canonical_eta_set_representative(mol, both), {})

    def test_bearing_ring_is_never_relabeled(self):
        """Guard 2: relabeling a load-bearing ring could hide a real face difference."""
        mol = self._mol(_ANSA_BIS_INDENYL)
        atoms = sorted(_ring_atoms(_ANSA_BIS_INDENYL, 5, 0))
        self.assertEqual(OINSanitizer.canonical_eta_set_representative(mol, atoms), {})

    def test_failure_returns_empty_mapping(self):
        mol = self._mol(_CP)
        self.assertEqual(OINSanitizer.canonical_eta_set_representative(mol, [0, 1]), {})
        self.assertEqual(OINSanitizer.canonical_eta_set_representative(mol, [0, 1, 99]), {})


class TestEncoderEndToEnd(unittest.TestCase):
    """The encoder itself, on real fixtures."""

    def _encode(self, name):
        return XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, name))

    def test_free_rings_always_emit_forward(self):
        """TiCat1's two orientation-free Cp rings both read '>', whatever face the
        embedding happened to present."""
        oin = self._encode("TiCat1.xyz")
        self.assertIn("{0>}", oin)
        self.assertIn("{1>}", oin)
        self.assertNotIn("<", oin)

    def test_rac_and_meso_stay_distinct(self):
        """THE regression floor. TiCat3 (rac) and TiCat4 (meso) are the same
        constitution and differ only in one indenyl's coordinated face. Their
        winding multisets must not collapse."""
        rac = sorted(c for c in self._encode("TiCat3.xyz") if c in "<>")
        meso = sorted(c for c in self._encode("TiCat4.xyz") if c in "<>")
        self.assertEqual(rac, [">", ">"])
        self.assertEqual(meso, ["<", ">"])
        self.assertNotEqual(rac, meso)

    def test_enantiomeric_ansa_stays_distinct_from_rac(self):
        """TiCat5 is TiCat3 with Z inverted: both rings flip, so the multiset does."""
        rac = sorted(c for c in self._encode("TiCat3.xyz") if c in "<>")
        ent = sorted(c for c in self._encode("TiCat5.xyz") if c in "<>")
        self.assertEqual(ent, ["<", "<"])
        self.assertNotEqual(rac, ent)


if __name__ == "__main__":
    unittest.main()
