"""Guard that the OIN encoder canonicalizes which symmetric-donor atom carries
the binding slot (Track A1, gap 1).

A monodentate carboxylate binds through one of two resonance-equivalent oxygens.
In any single Kekule structure they are distinguishable (=O vs -O), so which O
the 3D bond perception picked -- and thus which carries the ``{slot}`` marker --
drifts between the two round-trip directions (ABAZIO: ``O{n}C(=O)`` vs
``OC(=O{n})``). ``OINSanitizer.canonical_donor_representative`` remaps the binder
to a canonical member of its symmetry class before serialization, so the two
Kekule forms of one physical structure encode to a byte-identical OIN string --
without merging genuinely different (asymmetric) donors.
"""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem
from rdkit.Geometry import Point3D

from oinsmiles.utils.oin_aligner import OINSanitizer
from oinsmiles.utils.xyz2mol import get_oin_string


def _attach_conformer(rw, coords):
    conf = Chem.Conformer(rw.GetNumAtoms())
    for i in range(rw.GetNumAtoms()):
        conf.SetAtomPosition(i, Point3D(*[float(x) for x in coords[i]]))
    rw.AddConformer(conf)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return mol


def _acetate_pt(double_on_bound_o):
    """CH3-CO2 bound to PtCl3 through O1, identical coords for both variants;
    only the C=O double bond moves (bound O vs far O) -- the two Kekule forms
    of one physical acetate."""
    coords = np.array(
        [
            [0.00, 0.00, 0.00],  # 0 C methyl
            [1.50, 0.00, 0.00],  # 1 C carboxyl
            [2.20, 1.10, 0.00],  # 2 O1 (binds Pt)
            [2.20, -1.10, 0.00],  # 3 O2
            [2.20, 3.00, 0.00],  # 4 Pt
            [4.00, 3.00, 0.00],  # 5 Cl
            [0.40, 3.00, 0.00],  # 6 Cl
            [2.20, 4.80, 0.00],  # 7 Cl
            [-0.60, 0.80, 0.30],  # 8 H
            [-0.60, -0.80, 0.30],  # 9 H
            [-0.40, 0.00, -0.90],  # 10 H
        ]
    )
    rw = Chem.RWMol()
    for z in (6, 6, 8, 8, 78, 17, 17, 17, 1, 1, 1):
        rw.AddAtom(Chem.Atom(z))
    rw.AddBond(0, 1, Chem.BondType.SINGLE)
    if double_on_bound_o:
        rw.AddBond(1, 2, Chem.BondType.DOUBLE)
        rw.AddBond(1, 3, Chem.BondType.SINGLE)
    else:
        rw.AddBond(1, 2, Chem.BondType.SINGLE)
        rw.AddBond(1, 3, Chem.BondType.DOUBLE)
    rw.AddBond(2, 4, Chem.BondType.DATIVE)
    for c in (5, 6, 7):
        rw.AddBond(4, c, Chem.BondType.DATIVE)
    for h in (8, 9, 10):
        rw.AddBond(0, h, Chem.BondType.SINGLE)
    return _attach_conformer(rw, coords), coords


def _sulfonate_pt(double_on_bound_o):
    """CH3-SO3 bound to PtCl3 through O1, donor set size 3. Same coords for both
    variants; the bound O toggles single<->double (the far O compensates), so S
    valence stays 6 while the metal-facing O changes bond order."""
    coords = np.array(
        [
            [0.00, 0.00, 0.00],  # 0 C methyl
            [1.60, 0.00, 0.00],  # 1 S
            [2.30, 1.20, 0.00],  # 2 O1 (binds Pt)
            [2.30, -1.20, 0.00],  # 3 O2
            [1.60, 0.00, 1.40],  # 4 O3
            [2.30, 3.10, 0.00],  # 5 Pt
            [4.10, 3.10, 0.00],  # 6 Cl
            [0.50, 3.10, 0.00],  # 7 Cl
            [2.30, 4.90, 0.00],  # 8 Cl
            [-0.60, 0.80, 0.30],  # 9 H
            [-0.60, -0.80, 0.30],  # 10 H
            [-0.40, 0.00, -0.90],  # 11 H
        ]
    )
    rw = Chem.RWMol()
    for z in (6, 16, 8, 8, 8, 78, 17, 17, 17, 1, 1, 1):
        rw.AddAtom(Chem.Atom(z))
    rw.AddBond(0, 1, Chem.BondType.SINGLE)  # C-S
    if double_on_bound_o:
        rw.AddBond(1, 2, Chem.BondType.DOUBLE)  # S=O1 (bound)
        rw.AddBond(1, 3, Chem.BondType.SINGLE)  # S-O2
        rw.AddBond(1, 4, Chem.BondType.DOUBLE)  # S=O3
    else:
        rw.AddBond(1, 2, Chem.BondType.SINGLE)  # S-O1 (bound)
        rw.AddBond(1, 3, Chem.BondType.DOUBLE)  # S=O2
        rw.AddBond(1, 4, Chem.BondType.DOUBLE)  # S=O3
    rw.AddBond(2, 5, Chem.BondType.DATIVE)
    for c in (6, 7, 8):
        rw.AddBond(5, c, Chem.BondType.DATIVE)
    for h in (9, 10, 11):
        rw.AddBond(0, h, Chem.BondType.SINGLE)
    return _attach_conformer(rw, coords), coords


def _neutral_carboxylate_fragment():
    """CH3-C(=O)-O with BOTH oxygens neutral and 0-H, mimicking the encoder's
    post-neutralization frag_mol (SetFormalCharge(0) + NoImplicit)."""
    rw = Chem.RWMol()
    for z in (6, 6, 8, 8):
        rw.AddAtom(Chem.Atom(z))
    rw.AddBond(0, 1, Chem.BondType.SINGLE)
    rw.AddBond(1, 2, Chem.BondType.DOUBLE)
    rw.AddBond(1, 3, Chem.BondType.SINGLE)
    for a in rw.GetAtoms():
        a.SetNoImplicit(True)
    rw.GetAtomWithIdx(0).SetNumExplicitHs(3)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return mol, 2, 3  # (mol, o1_idx, o2_idx)


class TestCanonicalDonorBinding(unittest.TestCase):
    # -- Integration: the two Kekule forms MUST encode byte-identically --

    def test_carboxylate_kekule_forms_collapse(self):
        """Bound-O-single vs bound-O-double acetate -> one canonical OIN."""
        a, ca = _acetate_pt(double_on_bound_o=False)
        b, cb = _acetate_pt(double_on_bound_o=True)
        oin_a = get_oin_string(a, ca.copy())
        oin_b = get_oin_string(b, cb.copy())
        self.assertEqual(oin_a, oin_b, f"\n  A: {oin_a}\n  B: {oin_b}")

    def test_sulfonate_kekule_forms_collapse(self):
        """Donor set size 3: sulfonate O's are interchangeable -> one OIN."""
        a, ca = _sulfonate_pt(double_on_bound_o=False)
        b, cb = _sulfonate_pt(double_on_bound_o=True)
        oin_a = get_oin_string(a, ca.copy())
        oin_b = get_oin_string(b, cb.copy())
        self.assertEqual(oin_a, oin_b, f"\n  A: {oin_a}\n  B: {oin_b}")

    # -- The helper: symmetric donors collapse, asymmetric ones do not --

    def test_carboxylate_both_o_map_to_same_rep(self):
        """Either carboxyl O canonicalizes to the same representative atom."""
        mol, o1, o2 = _neutral_carboxylate_fragment()
        r1 = OINSanitizer.canonical_donor_representative(mol, o1)
        r2 = OINSanitizer.canonical_donor_representative(mol, o2)
        self.assertEqual(r1, r2)
        self.assertIn(r1, (o1, o2))

    def test_ester_asymmetric_donor_not_remapped(self):
        """Methyl acetate's two O's are inequivalent -> no remap (fail-safe)."""
        mol = Chem.MolFromSmiles("CC(=O)OC")
        carbonyl_o = next(
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetAtomicNum() == 8
            and any(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds())
        )
        self.assertEqual(OINSanitizer.canonical_donor_representative(mol, carbonyl_o), carbonyl_o)

    def test_single_donor_is_noop(self):
        """A donor with no symmetric partner (pyridine N) is left untouched."""
        mol = Chem.MolFromSmiles("c1ccncc1")
        n_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
        self.assertEqual(OINSanitizer.canonical_donor_representative(mol, n_idx), n_idx)

    # -- Anti-over-collapse: genuinely different donors stay distinct --

    def test_distinct_donors_stay_distinct(self):
        """Acetate and formate are different molecules -> different OIN."""
        acetate, ca = _acetate_pt(double_on_bound_o=False)
        oin_acetate = get_oin_string(acetate, ca.copy())

        # Formate: swap the methyl carbon (+3 H) for a lone H on the carboxyl C.
        rw = Chem.RWMol()
        for z in (1, 6, 8, 8, 78, 17, 17, 17):  # H, C, O1, O2, Pt, 3 Cl
            rw.AddAtom(Chem.Atom(z))
        rw.AddBond(0, 1, Chem.BondType.SINGLE)  # H-C
        rw.AddBond(1, 2, Chem.BondType.SINGLE)  # C-O1 (bound)
        rw.AddBond(1, 3, Chem.BondType.DOUBLE)  # C=O2
        rw.AddBond(2, 4, Chem.BondType.DATIVE)
        for c in (5, 6, 7):
            rw.AddBond(4, c, Chem.BondType.DATIVE)
        fcoords = np.array(
            [
                [-0.60, 0.00, 0.00],  # 0 H
                [1.50, 0.00, 0.00],  # 1 C
                [2.20, 1.10, 0.00],  # 2 O1 (binds Pt)
                [2.20, -1.10, 0.00],  # 3 O2
                [2.20, 3.00, 0.00],  # 4 Pt
                [4.00, 3.00, 0.00],  # 5 Cl
                [0.40, 3.00, 0.00],  # 6 Cl
                [2.20, 4.80, 0.00],  # 7 Cl
            ]
        )
        formate = _attach_conformer(rw, fcoords)
        oin_formate = get_oin_string(formate, fcoords.copy())
        self.assertNotEqual(oin_acetate, oin_formate)


if __name__ == "__main__":
    unittest.main()
