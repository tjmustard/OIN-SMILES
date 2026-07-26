"""Guards for the torsion-aware configurational oracle (Lane 7 / Task A).

The oracle answers "is the mirror reachable from this structure by rotating bonds?", which is
what separates *conformational* chirality (collapsing it is correct) from *configurational*
chirality (collapsing it is a losslessness failure). Three layers of guard:

1. **The torsion model itself** -- rotations must be rigid (bond lengths untouched), the
   dihedral-negating seed must actually negate dihedrals, and the batched Kabsch must agree
   with the scalar one in ``oracle.py`` including the reflection ban.
2. **Two bugs the cross-validation caught**, each with a regression test, because both made
   the tool confidently wrong in a way no verdict-level assertion would have localised:
   - cutting the *atom* instead of the *bond* when testing whether a bond is rotatable. A
     metal is a cut vertex, so blocking it detaches every ligand from every other and each
     metal-donor bond of a **chelate** looked rotatable -- which let fac-Ir(ppy)3's
     Delta/Lambda mirror be "reached" by swinging whole ligands off their own chelate rings.
   - enumerating automorphisms on the **H-explicit** graph, where methyls consume the whole
     ``maxMatches`` budget on permutations that leave every heavy atom fixed. Starving the
     automorphism set inflates every RMSD, and a freely rotating biphenyl read as
     configurational because of it.
3. **Verdicts on fixtures whose answer is known** by construction.
"""

import sys
import unittest
from pathlib import Path

import numpy as np
from rdkit import Chem

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from tools.injectivity.oracle import _kabsch_proper_rmsd, load_mol  # noqa: E402
from tools.injectivity.torsion_oracle import (  # noqa: E402
    _automorphism_perms,
    _batch_proper_rmsd,
    _heavy_indices,
    apply_torsions,
    configurational_verdict,
    dihedral_negating_theta,
    rotatable_torsions,
)

FIX = _ROOT / "tests" / "fixtures"
CISPLATIN = FIX / "CisPlatin.xyz"
FAC_IRPPY3 = FIX / "fac-Ir(ppy)3.xyz"
BINAP = FIX / "PdCl2-R-BINAP.xyz"
ZUMNEC = FIX / "ZUMNEC.xyz"
JEGKOW = FIX / "JEGKOW.xyz"
PHENPHOS = FIX / "PdCl2PhenPhosMe.xyz"


def _root(mol):
    return max(range(mol.GetNumAtoms()), key=lambda i: mol.GetAtomWithIdx(i).GetAtomicNum())


class TestTorsionModel(unittest.TestCase):
    def test_rotation_preserves_every_bond_length(self):
        """A torsion rotation is rigid: if it changed a bond it would not be a conformer move."""
        mol, coords = load_mol(BINAP)
        tors = rotatable_torsions(mol, _root(mol))
        self.assertGreater(len(tors), 0)
        rng = np.random.default_rng(7)
        moved = apply_torsions(coords, tors, rng.uniform(0, 360, len(tors)))
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            before = np.linalg.norm(coords[i] - coords[j])
            after = np.linalg.norm(moved[i] - moved[j])
            self.assertAlmostEqual(before, after, places=6)

    def test_dihedral_negating_seed_negates_dihedrals(self):
        from rdkit.Chem import rdMolTransforms as rmt

        mol, coords = load_mol(BINAP)
        tors = rotatable_torsions(mol, _root(mol))
        theta = dihedral_negating_theta(mol, coords, tors)
        moved = apply_torsions(coords, tors, theta)

        conf_a, conf_b = Chem.Conformer(mol.GetNumAtoms()), Chem.Conformer(mol.GetNumAtoms())
        for k in range(mol.GetNumAtoms()):
            conf_a.SetAtomPosition(k, coords[k].tolist())
            conf_b.SetAtomPosition(k, moved[k].tolist())
        heavy = set(_heavy_indices(mol).tolist())
        for a, b, _m, _d in tors:
            ra = sorted(
                n.GetIdx()
                for n in mol.GetAtomWithIdx(a).GetNeighbors()
                if n.GetIdx() != b and n.GetIdx() in heavy
            )
            rb = sorted(
                n.GetIdx()
                for n in mol.GetAtomWithIdx(b).GetNeighbors()
                if n.GetIdx() != a and n.GetIdx() in heavy
            )
            if not ra or not rb:
                continue
            before = rmt.GetDihedralDeg(conf_a, ra[0], a, b, rb[0])
            after = rmt.GetDihedralDeg(conf_b, ra[0], a, b, rb[0])
            self.assertAlmostEqual(after, -before, places=3)

    def test_batched_kabsch_matches_the_scalar_one(self):
        """Same numbers as ``oracle._kabsch_proper_rmsd``, reflection ban included."""
        mol, coords = load_mol(FAC_IRPPY3)
        heavy = _heavy_indices(mol)
        perms = _automorphism_perms(mol, heavy, 4000)
        mirror = coords.copy()
        mirror[:, 2] *= -1.0
        batched = _batch_proper_rmsd(coords[heavy], mirror[heavy], perms)
        for k in range(len(perms)):
            scalar = _kabsch_proper_rmsd(coords[heavy], mirror[heavy][perms[k]])
            self.assertAlmostEqual(batched[k], scalar, places=8)

    def test_reflection_is_forbidden(self):
        """A pure reflection must NOT score zero -- that is the whole point of the proper fit."""
        mol, coords = load_mol(FAC_IRPPY3)
        heavy = _heavy_indices(mol)
        mirror = coords.copy()
        mirror[:, 2] *= -1.0
        identity = np.arange(len(heavy))[None, :]
        self.assertGreater(_batch_proper_rmsd(coords[heavy], mirror[heavy], identity)[0], 1.0)


class TestChelateBondsAreNotRotatable(unittest.TestCase):
    """Regression: cut the BOND, not the atom (see the module docstring)."""

    def test_chelate_locked_metal_bonds_are_excluded(self):
        for path in (FAC_IRPPY3, ZUMNEC):
            mol, _coords = load_mol(path)
            tors = rotatable_torsions(mol, _root(mol))
            self.assertEqual(tors, [], f"{path.name}: a chelate-locked bond was called rotatable")

    def test_a_monodentate_metal_bond_is_still_rotatable(self):
        """The exclusion must not over-fire: a lone donor bond IS conformational freedom."""
        mol, _coords = load_mol(BINAP)
        tors = rotatable_torsions(mol, _root(mol))
        self.assertGreater(len(tors), 0)


class TestAutomorphismCompleteness(unittest.TestCase):
    """Regression: enumerate on the heavy skeleton, not the H-explicit graph."""

    def _starved_count(self, mol, heavy, cap=4000):
        n = mol.GetNumAtoms()
        matches = mol.GetSubstructMatches(mol, uniquify=False, useChirality=False, maxMatches=cap)
        pos = np.full(n, -1, dtype=int)
        pos[heavy] = np.arange(len(heavy))
        seen = set()
        for m in matches:
            if len(m) != n:
                continue
            img = np.fromiter((m[i] for i in heavy), dtype=int, count=len(heavy))
            if np.any(pos[img] < 0):
                continue
            seen.add(pos[img].tobytes())
        return len(seen)

    def test_heavy_skeleton_finds_more_than_the_starved_path(self):
        mol, _coords = load_mol(PHENPHOS)
        heavy = _heavy_indices(mol)
        complete = len(_automorphism_perms(mol, heavy, 4000))
        starved = self._starved_count(mol, heavy)
        self.assertGreater(
            complete, starved, "H-explicit starvation is no longer reproduced -- retune the guard"
        )

    def test_every_permutation_is_a_real_automorphism(self):
        """Preserves adjacency, element, H count and formal charge."""
        mol, _coords = load_mol(PHENPHOS)
        heavy = _heavy_indices(mol)
        perms = _automorphism_perms(mol, heavy, 4000)
        pos = {int(a): k for k, a in enumerate(heavy)}
        edges = {
            frozenset((pos[b.GetBeginAtomIdx()], pos[b.GetEndAtomIdx()]))
            for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in pos and b.GetEndAtomIdx() in pos
        }
        z = [mol.GetAtomWithIdx(int(i)).GetAtomicNum() for i in heavy]
        nh = [
            sum(1 for x in mol.GetAtomWithIdx(int(i)).GetNeighbors() if x.GetAtomicNum() == 1)
            for i in heavy
        ]
        self.assertTrue(any(np.array_equal(p, np.arange(len(heavy))) for p in perms))
        for p in perms:
            for u, v in (tuple(e) for e in edges):
                self.assertIn(frozenset((int(p[u]), int(p[v]))), edges)
            for k in range(len(heavy)):
                self.assertEqual(z[p[k]], z[k])
                self.assertEqual(nh[p[k]], nh[k])


class TestVerdictsOnKnownAnswers(unittest.TestCase):
    """Fixtures whose answer is known by construction. Budget trimmed for suite runtime."""

    def test_achiral_control_is_rigid_achiral(self):
        v = configurational_verdict(CISPLATIN, restarts=2)
        self.assertEqual(v.verdict, "rigid_achiral")

    def test_square_planar_four_donors_is_rigid_achiral(self):
        """A square-planar complex is planar, hence achiral however different its donors."""
        v = configurational_verdict(JEGKOW, restarts=2)
        self.assertEqual(v.verdict, "rigid_achiral")

    def test_metal_delta_lambda_is_configurational(self):
        for path in (FAC_IRPPY3, ZUMNEC):
            v = configurational_verdict(path, restarts=2)
            self.assertEqual(v.verdict, "configurational", path.name)
            self.assertEqual(v.n_torsions, 0)

    def test_hindered_biaryl_is_configurational_with_a_passing_control(self):
        """The strong form: the mirror is unreachable while a reachable target is recovered."""
        v = configurational_verdict(BINAP, restarts=4)
        self.assertEqual(v.verdict, "configurational")
        self.assertGreater(v.d_mirror, 1.0)
        self.assertLess(v.d_control, v.threshold)

    def test_verdict_carries_its_own_search_size(self):
        """A conclusion must never be quotable without the budget that produced it."""
        d = configurational_verdict(BINAP, restarts=2).to_dict()
        for key in ("budget_evaluations", "threshold", "d_control", "n_autos_total"):
            self.assertIn(key, d)


if __name__ == "__main__":
    unittest.main()
