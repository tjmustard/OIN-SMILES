"""Guards for the structural-edit twin operators (Lane 7 / Task C).

``mirror_z`` is a coordinate transform and cannot produce an invalid structure. These two are
**structural edits** and can, so the guards fall into three groups:

1. **The edit is rigid where it claims to be.** ``swap_donor`` moves whole ligands about the
   metal, so every intra-ligand bond length must survive untouched; ``invert_axial`` must move
   exactly one dihedral.
2. **The edit is filtered, not trusted.** The vdW clash gate has to actually reject a swap that
   drives ligands together, and ``geometry_ok`` has to be false when it does.
3. **The distinctness certification is honest in both directions.** Swapping two *trans*
   donors is a 180 deg rotation of the whole complex -- the same isomer -- and must be reported
   as NOT distinct. An instrument that called those distinct would manufacture blind spots,
   which is exactly how the hand-written donor-swap probe went wrong. Swapping two *cis*
   donors, or exchanging symmetry-equivalent donors, are the matching positive and negative
   controls.
"""

import os
import sys
import unittest
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from tools.injectivity.oracle import load_mol  # noqa: E402
from tools.injectivity.twin_operators import (  # noqa: E402
    donor_groups,
    enumerate_donor_swaps,
    invert_axial,
    invert_tetrahedral,
    probe_operator,
    swap_donor,
)

FIX = _ROOT / "tests" / "fixtures"
PTMIXED = FIX / "PtMeNH3ClBr-Cis.xyz"  # four different donors, square planar
CISPLATIN = FIX / "CisPlatin.xyz"
YESKOZ = FIX / "YESKOZ.xyz"  # two symmetry-equivalent hindered axes, opposite sign
JEGKOW = FIX / "JEGKOW.xyz"
POJJOP = FIX / "POJJOP.xyz"  # sole stereocentre is a metal-bound 2 deg amine (P3)


def _amine_centre(mol):
    from tools.injectivity.config_oracle import bound_amine_centers

    centres = bound_amine_centers(mol)
    assert centres, "POJJOP must expose a metal-bound amine stereocentre"
    return centres[0].atom_idx


class TestSwapDonorIsRigid(unittest.TestCase):
    def test_intra_ligand_bond_lengths_are_untouched(self):
        """Only which site each donor occupies may change -- nothing inside a ligand."""
        mol, coords = load_mol(PTMIXED)
        groups = donor_groups(mol)
        a, b = enumerate_donor_swaps(mol)[0]
        twin = swap_donor(mol, coords, a, b)
        self.assertEqual(twin.error, "")
        member = {}
        for donor, atoms in groups.items():
            for at in atoms:
                member[at] = donor
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if member.get(i) is None or member.get(i) != member.get(j):
                continue  # a metal-donor bond, whose length is deliberately re-seated
            self.assertAlmostEqual(
                float(np.linalg.norm(coords[i] - coords[j])),
                float(np.linalg.norm(twin.coords[i] - twin.coords[j])),
                places=6,
            )

    def test_chelate_donors_are_refused(self):
        """One arm of a chelate cannot be moved rigidly, so the operator must decline."""
        mol, coords = load_mol(JEGKOW)
        donors = sorted(
            n.GetIdx()
            for n in next(a for a in mol.GetAtoms() if a.GetSymbol() == "Rh").GetNeighbors()
        )
        chelate = [d for d in donors if d not in donor_groups(mol)]
        self.assertTrue(chelate, "JEGKOW has a bidentate; its donors must be non-detachable")
        other = next(d for d in donor_groups(mol))
        self.assertIn("chelate", swap_donor(mol, coords, chelate[0], other).error)

    def test_symmetry_equivalent_donors_are_not_enumerated(self):
        """Exchanging equivalent donors is the identity on the isomer, not a twin."""
        mol, _coords = load_mol(CISPLATIN)
        for a, b in enumerate_donor_swaps(mol):
            self.assertNotEqual(
                mol.GetAtomWithIdx(a).GetSymbol(),
                mol.GetAtomWithIdx(b).GetSymbol(),
                "two equivalent donors were offered as a swap",
            )


class TestClashGateFilters(unittest.TestCase):
    def test_some_swap_is_rejected_by_the_clash_gate(self):
        """The gate must bite: cisplatin has swaps that drive ligands into each other."""
        mol, coords = load_mol(CISPLATIN)
        twins = [swap_donor(mol, coords, a, b) for a, b in enumerate_donor_swaps(mol)]
        self.assertTrue(
            any(not t.geometry_ok for t in twins), "clash gate never fired on cisplatin"
        )
        self.assertTrue(any(t.geometry_ok for t in twins), "clash gate rejected everything")

    def test_a_rejected_twin_is_never_scored(self):
        mol, coords = load_mol(CISPLATIN)
        bad = next(
            t
            for t in (swap_donor(mol, coords, a, b) for a, b in enumerate_donor_swaps(mol))
            if not t.geometry_ok
        )
        o = probe_operator(CISPLATIN, bad)
        self.assertIsNone(o.oracle_distinct)
        self.assertIsNone(o.key_equal)


class TestDistinctnessIsHonestBothWays(unittest.TestCase):
    """The negative control matters as much as the positive one."""

    def test_trans_swap_is_the_same_isomer(self):
        """Exchanging two donors that sit trans is a 180 deg rotation -- NOT a new isomer."""
        mol, coords = load_mol(PTMIXED)
        results = {}
        for a, b in enumerate_donor_swaps(mol):
            o = probe_operator(PTMIXED, swap_donor(mol, coords, a, b), restarts=2)
            if o.geometry_ok:
                results[(a, b)] = o
        same = [o for o in results.values() if o.oracle_distinct is False]
        self.assertTrue(same, "no trans swap was recognised as the same isomer")
        for o in same:
            self.assertTrue(o.raw_equal, "the encoder must give one string for one isomer")

    def test_cis_swap_is_a_distinct_isomer_the_encoder_separates(self):
        mol, coords = load_mol(PTMIXED)
        distinct = []
        for a, b in enumerate_donor_swaps(mol):
            o = probe_operator(PTMIXED, swap_donor(mol, coords, a, b), restarts=2)
            if o.geometry_ok and o.oracle_distinct:
                distinct.append(o)
        self.assertTrue(distinct, "no cis swap produced a distinct isomer")
        for o in distinct:
            self.assertFalse(o.key_equal, f"key folds a distinct positional isomer: {o.detail}")


class TestInvertAxial(unittest.TestCase):
    def test_flips_exactly_one_axis(self):
        from oinsmiles.oin.axial import detect_axial_axes

        mol, coords = load_mol(YESKOZ)
        base = [ax for ax in detect_axial_axes(mol) if ax.hindered]
        self.assertGreaterEqual(len(base), 2, "YESKOZ must carry two hindered axes")
        twin = invert_axial(mol, coords, 0)
        self.assertEqual(twin.error, "")
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, twin.coords[i].tolist())
        after = {(a.a1, a.a2): a.sign for a in detect_axial_axes(mol) if a.hindered}
        flipped = [
            k for k, v in after.items() if v != dict(((a.a1, a.a2), a.sign) for a in base)[k]
        ]
        self.assertEqual(len(flipped), 1, f"expected one axis to flip, got {flipped}")

    def test_the_flipped_axis_is_frozen_during_certification(self):
        """Otherwise the torsion orbit simply rotates the edit back and nothing is ever distinct."""
        mol, coords = load_mol(YESKOZ)
        twin = invert_axial(mol, coords, 0)
        self.assertEqual(twin.freeze_bonds, ((twin_axis(twin)[0], twin_axis(twin)[1]),))
        o = probe_operator(YESKOZ, twin, restarts=2)
        self.assertTrue(o.oracle_distinct, "single-axis flip must be a distinct diastereomer")

    @unittest.skipIf(os.environ.get("OIN_EMIT_AXIAL"), "measures the default (token OFF) build")
    def test_default_build_collapses_the_single_axis_flip(self):
        """Assert-current-behavior: with the axial token off, the OIN is byte-identical."""
        mol, coords = load_mol(YESKOZ)
        o = probe_operator(YESKOZ, invert_axial(mol, coords, 0), restarts=2)
        self.assertTrue(o.oracle_distinct)
        self.assertTrue(o.raw_equal, "P2 multi-axis regression: the flip now raw-diverges")


def twin_axis(twin):
    """The (a1, a2) the axial operator recorded as frozen."""
    return twin.freeze_bonds[0]


class TestInvertTetrahedral(unittest.TestCase):
    def test_non_stereocentre_is_refused(self):
        """A metal-bound ammine has three equivalent H: exchanging one is the identity."""
        mol, coords = load_mol(PTMIXED)
        ammine_n = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "N")
        self.assertIn("stereocentre", invert_tetrahedral(mol, coords, ammine_n).error)

    def test_ring_locked_centre_is_refused(self):
        mol, coords = load_mol(JEGKOW)
        ring_atom = next(
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.IsInRing() and a.GetSymbol() == "C" and a.GetDegree() >= 3
        )
        self.assertNotEqual(invert_tetrahedral(mol, coords, ring_atom).error, "")

    def test_branch_internal_geometry_survives(self):
        """The exchange must be an isometry: a rotate-and-scale would shrink the moved branch."""
        mol, coords = load_mol(POJJOP)
        centre = _amine_centre(mol)
        twin = invert_tetrahedral(mol, coords, centre)
        self.assertEqual(twin.error, "")
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if centre in (i, j):
                continue
            self.assertAlmostEqual(
                float(np.linalg.norm(coords[i] - coords[j])),
                float(np.linalg.norm(twin.coords[i] - twin.coords[j])),
                places=6,
            )

    def test_it_actually_inverts_the_centre(self):
        """The signed tetrahedral volume the P3 oracle reads must change sign."""
        from tools.injectivity.config_oracle import bound_amine_centers

        mol, coords = load_mol(POJJOP)
        centre = _amine_centre(mol)
        before = {c.atom_idx: c.sign for c in bound_amine_centers(mol)}
        twin = invert_tetrahedral(mol, coords, centre)
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, twin.coords[i].tolist())
        after = {c.atom_idx: c.sign for c in bound_amine_centers(mol)}
        self.assertEqual(after[centre], -before[centre])


if __name__ == "__main__":
    unittest.main()
