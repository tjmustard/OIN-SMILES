"""Tests for the coordinate-only ligand-attachment check (v0.4.7, `OIN_ATTACH_CHECK`).

The tests that matter here are the ones that would fail if somebody "simplified" the check
into the trap it was written to avoid -- reading attachment off the generator's own bond
graph. `test_detached_ligand_is_caught_even_though_its_bond_survives` is that guard, and it
is the whole point of the module.
"""

import unittest

import numpy as np
from rdkit import Chem

from oinsmiles.generation.attach_check import (
    conformer_ligands_attached,
    encoder_donor_set,
    group_sites,
    ligands_attached,
)
from oinsmiles.oin.levers import held_off, lever_enabled


def _square_planar_ptcl2n2():
    """Pt with 4 donors at 2.0 A along +-x / +-y. Atom 0 is the metal."""
    z = [78, 17, 17, 7, 7]
    c = np.array(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [-2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, -2.0, 0.0]]
    )
    return z, c


class TestEncoderDonorSet(unittest.TestCase):
    def test_finds_all_four_donors_at_bonding_distance(self):
        z, c = _square_planar_ptcl2n2()
        midx, donors = encoder_donor_set(z, c)
        self.assertEqual(midx, 0)
        self.assertEqual(donors, {1, 2, 3, 4})

    def test_a_donor_moved_out_of_range_is_not_reported(self):
        z, c = _square_planar_ptcl2n2()
        c[3] = [0.0, 6.0, 0.0]  # push one N well beyond any covalent tolerance
        _midx, donors = encoder_donor_set(z, c)
        self.assertNotIn(3, donors)
        self.assertEqual(donors, {1, 2, 4})

    def test_no_transition_metal_yields_no_donor_set(self):
        midx, donors = encoder_donor_set([6, 1, 1, 1, 1], np.zeros((5, 3)))
        self.assertIsNone(midx)
        self.assertEqual(donors, set())


class TestGroupSites(unittest.TestCase):
    def test_a_cp_ring_collapses_to_one_site(self):
        # 5 carbons on a 1.2 A-radius circle: adjacent separation ~1.41 A < 1.6 cutoff.
        ang = np.linspace(0, 2 * np.pi, 5, endpoint=False)
        ring = np.stack([1.2 * np.cos(ang), 1.2 * np.sin(ang), np.full(5, 2.0)], axis=1)
        coords = np.vstack([np.zeros((1, 3)), ring])
        self.assertEqual(len(group_sites([1, 2, 3, 4, 5], coords)), 1)

    def test_well_separated_donors_stay_separate_sites(self):
        _z, c = _square_planar_ptcl2n2()
        self.assertEqual(len(group_sites([1, 2, 3, 4], c)), 4)


class TestLigandsAttached(unittest.TestCase):
    def test_intact_complex_is_accepted(self):
        z, c = _square_planar_ptcl2n2()
        ok, detail = ligands_attached([1, 2, 3, 4], z, c)
        self.assertTrue(ok)
        self.assertEqual(detail["sites_lost"], 0)

    def test_a_site_that_has_left_the_sphere_is_rejected(self):
        z, c = _square_planar_ptcl2n2()
        c[3] = [0.0, 6.0, 0.0]
        ok, detail = ligands_attached([1, 2, 3, 4], z, c)
        self.assertFalse(ok)
        self.assertEqual(detail["sites_lost"], 1)
        self.assertEqual(detail["lost_atom_indices"], [3])

    def test_ring_slip_keeps_the_site_and_is_accepted(self):
        """eta5 -> eta2 must NOT be rejected. This is the failure mode that killed both
        predicates the promote lane proposed: MEDZUR's ring slips, the raw donor count drops
        10 -> 7, and the molecule still round-trips. A site keeps its meaning as long as ONE
        of its atoms is still bonded.

        The ring is translated RIGIDLY -- every C-C stays 1.41 A -- so this is a slip and not
        a shattered ligand. An earlier version of this fixture moved two carbons out of the
        ring plane, which detaches them from the ring rather than the ring from the metal, and
        the check rejected it. Correctly: those two carbons then form a site of their own with
        nothing bonded in it.
        """
        ang = np.linspace(0, 2 * np.pi, 5, endpoint=False)
        ring = np.stack([1.9 + 1.2 * np.cos(ang), 1.2 * np.sin(ang), np.full(5, 2.0)], axis=1)
        coords = np.vstack([np.zeros((1, 3)), ring])
        z = [40] + [6] * 5  # Zr + Cp

        self.assertEqual(len(group_sites([1, 2, 3, 4, 5], coords)), 1, "the ring is still a ring")
        _m, donors = encoder_donor_set(z, coords)
        self.assertTrue(donors, "some ring atoms must still be bonded")
        self.assertLess(len(donors), 5, "premise: the ring has slipped, not stayed eta5")

        ok, detail = ligands_attached([1, 2, 3, 4, 5], z, coords)
        self.assertTrue(ok, "a slipped ring is still an attached ligand")
        self.assertEqual(detail["sites_lost"], 0)

    def test_the_same_ring_fully_departed_is_rejected(self):
        """The other side of the slip test: move the identical rigid ring right out of the
        coordination sphere and the site must go empty."""
        ang = np.linspace(0, 2 * np.pi, 5, endpoint=False)
        ring = np.stack([8.0 + 1.2 * np.cos(ang), 1.2 * np.sin(ang), np.full(5, 2.0)], axis=1)
        coords = np.vstack([np.zeros((1, 3)), ring])
        z = [40] + [6] * 5
        ok, detail = ligands_attached([1, 2, 3, 4, 5], z, coords)
        self.assertFalse(ok)
        self.assertEqual(detail["sites_lost"], 1)

    def test_no_claimed_donors_abstains_rather_than_rejecting(self):
        z, c = _square_planar_ptcl2n2()
        ok, _ = ligands_attached([], z, c)
        self.assertTrue(ok)


class TestTrapAvoidance(unittest.TestCase):
    """The §6.1 trap, as an executable guard.

    `_coordination_vectors` reads donors from `metal.GetBonds()`. A ligand that has left the
    coordination sphere KEEPS its bond object, so a check built that way certifies exactly the
    structures it exists to catch. This test builds precisely that situation -- an RDKit mol
    whose metal-N bond is intact while the N sits 6 A away -- and requires the check to reject.
    """

    def _mol_with_surviving_bond_to_a_departed_ligand(self):
        rw = Chem.RWMol()
        for zn in (78, 17, 17, 7, 7):
            rw.AddAtom(Chem.Atom(zn))
        for j in (1, 2, 3, 4):
            rw.AddBond(0, j, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, p in enumerate(
            [(0, 0, 0), (2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, 6, 0)]  # atom 4 has departed
        ):
            conf.SetAtomPosition(i, tuple(float(v) for v in p))
        mol.AddConformer(conf)
        return mol

    def test_detached_ligand_is_caught_even_though_its_bond_survives(self):
        mol = self._mol_with_surviving_bond_to_a_departed_ligand()
        metal = mol.GetAtomWithIdx(0)
        self.assertEqual(len(metal.GetBonds()), 4, "premise: the bond object still exists")
        self.assertFalse(
            conformer_ligands_attached(mol),
            "a check that trusted GetBonds() would pass this; it must not",
        )

    def test_the_same_mol_with_the_ligand_in_place_is_accepted(self):
        mol = self._mol_with_surviving_bond_to_a_departed_ligand()
        mol.GetConformer().SetAtomPosition(4, (0.0, -2.0, 0.0))
        self.assertTrue(conformer_ligands_attached(mol))

    def test_unevaluable_input_abstains_instead_of_rejecting(self):
        """No conformer -> the predicate has not judged the structure, so it must not reject
        one. Failing closed here would silently disable the lever it guards."""
        mol = Chem.MolFromSmiles("CCO")
        self.assertTrue(conformer_ligands_attached(mol))


class TestRealCrystalInputsAreNotRejected(unittest.TestCase):
    """The check must not reject genuine structures. Ferrocene is the sharpest available
    fixture: two eta5 rings, i.e. exactly the coordination the check reasons about."""

    def test_ferrocene_like_sandwich_keeps_both_ring_sites(self):
        ang = np.linspace(0, 2 * np.pi, 5, endpoint=False)
        top = np.stack([1.2 * np.cos(ang), 1.2 * np.sin(ang), np.full(5, 1.65)], axis=1)
        bot = np.stack([1.2 * np.cos(ang), 1.2 * np.sin(ang), np.full(5, -1.65)], axis=1)
        coords = np.vstack([np.zeros((1, 3)), top, bot])
        z = [26] + [6] * 10
        ok, detail = ligands_attached(list(range(1, 11)), z, coords)
        self.assertTrue(ok)
        self.assertEqual(detail["n_sites_claimed"], 2, "the two rings are two sites, not ten")


class TestLeverRegistration(unittest.TestCase):
    def test_lever_is_registered_as_held_off_with_a_reason(self):
        reasons = held_off()
        self.assertIn("OIN_ATTACH_CHECK", reasons)
        why = reasons["OIN_ATTACH_CHECK"]
        self.assertIn("7", why, "the measured separation belongs in the justification")
        self.assertIn("POVPIA", why, "the known residual must be stated, not hidden")

    def test_lever_defaults_off(self):
        import os

        prior = os.environ.pop("OIN_ATTACH_CHECK", None)
        try:
            self.assertFalse(lever_enabled("OIN_ATTACH_CHECK"))
        finally:
            if prior is not None:
                os.environ["OIN_ATTACH_CHECK"] = prior


class TestAcceptancePredicateWiring(unittest.TestCase):
    def test_check_is_unreachable_when_accept_scored_is_off(self):
        """With `independent_confirm=True` (the default path) the lever's branch is never
        entered, so the check cannot alter the default output. Asserted on the source rather
        than by generating, because 'byte-identical by construction' is the claim being made."""
        import inspect

        from oinsmiles.generation import metallogen_adapter

        src = inspect.getsource(metallogen_adapter._reencode_key_matches)
        guard = "if not independent_confirm and fast is not None:"
        self.assertIn(guard, src)
        after = src.split(guard, 1)[1]
        before = src.split(guard, 1)[0]
        self.assertIn("conformer_ligands_attached", after)
        self.assertNotIn("conformer_ligands_attached", before)


if __name__ == "__main__":
    unittest.main()
