"""Unit guards for SL3 greedy (difficulty-ordered, collision-aware) placement.

Greedy placement is an opt-in *variant* of the ``option=3`` kabsch embed
(``get_embedding(..., greedy=True)`` / ``OIN_GREEDY_PLACEMENT`` / ``ff_params["greedy"]``).
Default OFF keeps the embed pool byte-identical to pristine (guarded by the existing
``test_kabsch_placement`` end-to-end tests + an external vs-main diff in the squash body).

These tests pin the greedy path's contracts, and each fails against pre-SL3 code
(``get_embedding`` had no ``greedy`` kwarg; the helper functions did not exist):

* the difficulty key orders most-constrained-first (chelate/pincer -> haptic -> bulky
  monodentate -> lone halide last);
* greedy placement keeps every donor on its slot vector, the metal at the origin, no
  NaNs, and is deterministic for a fixed seed;
* the collision-aware spin and the prefix settle are STRICTLY guarded -- neither ever
  ships a clashier structure than it started with;
* the gate is OFF by default.
"""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.metallogen_adapter import convert_oin_to_msmiles
from oinsmiles.generator3d import _greedy_enabled, om
from oinsmiles.generator3d.embed import (
    _interligand_clash_cost,
    _ligand_difficulty_key,
    _pack_ligands,
    _place_ligand_collision_aware,
    get_embedding,
)


class _FakeMol:
    def __init__(self, n):
        self.atom_list = [None] * n


class _FakeLigand:
    """Minimal stand-in exposing just what ``_ligand_difficulty_key`` reads."""

    def __init__(self, binding_infos, n_atoms):
        self.binding_infos = binding_infos
        self.molecule = _FakeMol(n_atoms)


class TestDifficultyOrdering(unittest.TestCase):
    """The placement-order key sorts the most-constrained ligand first."""

    def test_chelate_before_haptic_before_mono_before_halide(self):
        # denticity-2 chelate (6 atoms), haptic 5-ring, bulky monodentate (5 atoms),
        # lone halide (1 atom). Free coords: give each a plausible spread.
        chelate = _FakeLigand([([0], 1), ([5], 2)], 6)
        haptic = _FakeLigand([([0, 1, 2, 3, 4], 1)], 5)
        bulky_mono = _FakeLigand([([0], 1)], 5)
        halide = _FakeLigand([([0], 1)], 1)

        coords = {
            id(chelate): np.random.RandomState(0).randn(6, 3) * 1.5,
            id(haptic): np.random.RandomState(1).randn(5, 3) * 1.2,
            id(bulky_mono): np.random.RandomState(2).randn(5, 3) * 0.9,
            id(halide): np.zeros((1, 3)),
        }
        ligs = [halide, bulky_mono, haptic, chelate]  # deliberately scrambled
        ordered = sorted(ligs, key=lambda lg: _ligand_difficulty_key(lg, coords[id(lg)]))

        self.assertIs(ordered[0], chelate, "chelate (highest denticity) must place first")
        self.assertIs(ordered[-1], halide, "lone halide must place last")
        # Among the denticity-1 ligands, the haptic ring outranks the monodentates.
        self.assertLess(ordered.index(haptic), ordered.index(bulky_mono))

    def test_bulkier_monodentate_sorts_before_smaller(self):
        big = _FakeLigand([([0], 1)], 6)
        small = _FakeLigand([([0], 1)], 2)
        big_coords = np.array(
            [[0, 0, 0], [2, 0, 0], [0, 2, 0], [0, 0, 2], [2, 2, 0], [1, 1, 2]], float
        )
        small_coords = np.array([[0, 0, 0], [0.3, 0, 0]], float)
        order = sorted(
            [small, big],
            key=lambda lg: _ligand_difficulty_key(lg, big_coords if lg is big else small_coords),
        )
        self.assertIs(order[0], big)


class TestGreedyPlacement(unittest.TestCase):
    """End-to-end greedy placement holds the geometry invariants and is deterministic."""

    CISPLATIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
    FERROCENE = (
        "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
    )
    EN_CHELATE = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}CCN{3}"

    def _place(self, oin, greedy):
        mc = om.get_om_from_modified_smiles(convert_oin_to_msmiles(oin))
        pos = get_embedding(mc, scale=1.0, option=3, align=True, seed=42, greedy=greedy)
        self.assertIsNotNone(pos, f"option=3 greedy={greedy} returned None")
        return mc, np.asarray(pos)

    def _assert_donors_on_slots(self, mc, pos, tol=0.95):
        metal = pos[mc.metal_index]
        dvec = np.asarray(mc.geometry_type.direction_vector, dtype=float)
        aidx = mc.get_atom_indices_for_each_ligand()
        for i, lig in enumerate(mc.ligands):
            for atom_locals, slot in lig.binding_infos:
                donor = pos[[aidx[i][a] for a in atom_locals]].mean(axis=0)
                rel = donor - metal
                v = dvec[slot - 1] / np.linalg.norm(dvec[slot - 1])
                cos = float(np.dot(rel, v) / (np.linalg.norm(rel) + 1e-9))
                self.assertGreater(cos, tol, f"lig{i} slot{slot} off-vector (cos={cos:.3f})")

    def test_monodentate_greedy_donors_on_slots(self):
        mc, pos = self._place(self.CISPLATIN, greedy=True)
        self.assertEqual(pos.shape, (mc.num_atom, 3))
        self.assertFalse(bool(np.isnan(pos).any()))
        np.testing.assert_allclose(pos[mc.metal_index], 0.0, atol=1e-6)
        self._assert_donors_on_slots(mc, pos)

    def test_haptic_greedy_donors_on_slots(self):
        mc, pos = self._place(self.FERROCENE, greedy=True)
        self.assertFalse(bool(np.isnan(pos).any()))
        np.testing.assert_allclose(pos[mc.metal_index], 0.0, atol=1e-6)
        self._assert_donors_on_slots(mc, pos)

    def test_chelate_greedy_donors_on_slots(self):
        mc, pos = self._place(self.EN_CHELATE, greedy=True)
        self.assertFalse(bool(np.isnan(pos).any()))
        self._assert_donors_on_slots(mc, pos)

    def test_greedy_is_deterministic(self):
        _, a = self._place(self.CISPLATIN, greedy=True)
        _, b = self._place(self.CISPLATIN, greedy=True)
        np.testing.assert_array_equal(a, b)


class TestGreedyGate(unittest.TestCase):
    """The greedy gate is OFF by default and honours env + ff_params."""

    def setUp(self):
        self._saved = os.environ.pop("OIN_GREEDY_PLACEMENT", None)

    def tearDown(self):
        os.environ.pop("OIN_GREEDY_PLACEMENT", None)
        if self._saved is not None:
            os.environ["OIN_GREEDY_PLACEMENT"] = self._saved

    def test_off_by_default(self):
        self.assertFalse(_greedy_enabled(None))
        self.assertFalse(_greedy_enabled({}))
        self.assertFalse(_greedy_enabled({"greedy": False}))

    def test_ff_params_enables(self):
        self.assertTrue(_greedy_enabled({"greedy": True}))

    def test_env_enables(self):
        os.environ["OIN_GREEDY_PLACEMENT"] = "1"
        self.assertTrue(_greedy_enabled(None))


class TestGuardedMonotonicity(unittest.TestCase):
    """The collision-aware spin and the prefix settle never worsen clash."""

    def test_collision_aware_never_increases_cross_clash(self):
        # A ligand of 3 atoms sitting on an axis through the origin, and one committed
        # atom placed to clash with the ligand's initial orientation. The spin search
        # must not return a pose with a higher cross-clash cost than the base pose.
        positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.2, 0.0]])  # metal + 1 committed
        committed = [0, 1]
        placed = np.array([[0.0, 0.0, 1.5], [0.9, 0.1, 1.5], [-0.9, -0.1, 1.5]])
        axis = np.array([0.0, 0.0, 1.0])
        pivot = np.array([0.0, 0.0, 1.5])
        base_cost = _interligand_clash_cost(placed, positions[committed])
        out = _place_ligand_collision_aware(placed, (pivot, axis), positions, committed)
        out_cost = _interligand_clash_cost(out, positions[committed])
        self.assertLessEqual(out_cost, base_cost + 1e-9)

    def test_collision_aware_no_axis_is_identity(self):
        placed = np.random.RandomState(3).randn(4, 3)
        out = _place_ligand_collision_aware(placed, None, np.zeros((1, 3)), [0])
        np.testing.assert_array_equal(out, placed)

    def test_prefix_settle_scores_only_committed_and_never_worsens(self):
        # Two ligands' worth of atoms; only the first (indices 1..3) is committed, the
        # rest sit at the origin. Settling with score_indices must ignore the origin
        # pile and never increase the committed-subset repulsive potential.
        from oinsmiles.generator3d.embed import get_repulsive_potential

        positions = np.zeros((7, 3))
        positions[0] = [0, 0, 0]  # metal
        positions[1:4] = [[1.0, 0.0, 0.0], [1.5, 0.5, 0.0], [1.5, -0.5, 0.0]]
        committed = [0, 1, 2, 3]
        axes = [([1, 2, 3], (np.array([1.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])))]
        before = get_repulsive_potential(positions[committed])
        out = _pack_ligands(positions, axes, 0, n_passes=1, score_indices=committed)
        after = get_repulsive_potential(out[committed])
        self.assertLessEqual(after, before + 1e-9)
        # The not-yet-placed atoms (still at the origin) are untouched by the settle.
        np.testing.assert_array_equal(out[4:], positions[4:])


if __name__ == "__main__":
    unittest.main()
