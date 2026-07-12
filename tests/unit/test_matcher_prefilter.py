"""Guard tests for the geometry-matcher's batched-numpy candidate prefilter.

``OINDiscreteAligner._map_to_template`` no longer scores all ``n!`` slot
assignments with scipy; a vectorised batched-numpy Kabsch nominates a small
candidate set (``_candidate_permutations``) and scipy scores only those. These
tests pin the two load-bearing properties:

  * **Byte-identity** -- the prefiltered result (winning permutation, best_rmsd
    to the bit, best rotation matrix to the bit) equals the full exhaustive
    scipy sweep, across coordination numbers 2-6, ideal + noisy + the worst
    symmetric-degeneracy cases. If this ever breaks, generated geometry changes.
  * **Pruning is real** -- for a symmetric sphere the candidate set is a strict
    subset of all permutations (otherwise the refactor bought nothing).
"""

import itertools
import unittest

import numpy as np
from scipy.spatial.transform import Rotation

from oinsmiles.utils.oin_aligner import TEMPLATES, OINDiscreteAligner


def _virtual_atoms(coords):
    return [{"coords": np.asarray(c, dtype=float)} for c in coords]


def _exhaustive(virtual_atoms, template_vectors):
    """Verbatim pre-refactor _map_to_template: score every permutation with scipy."""
    n_atoms = len(virtual_atoms)
    n_slots = len(template_vectors)
    input_vecs = np.array([a["coords"] for a in virtual_atoms])
    input_norms = input_vecs / (np.linalg.norm(input_vecs, axis=1)[:, None] + 1e-9)
    best_rmsd = float("inf")
    best_perm = None
    best_R = np.eye(3)
    for slot_indices in itertools.permutations(range(n_slots), n_atoms):
        target_vecs = template_vectors[list(slot_indices)]
        try:
            R, rmsd = Rotation.align_vectors(target_vecs, input_norms)
        except Exception:  # noqa: BLE001
            continue
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_perm = slot_indices
            best_R = R.as_matrix()
    return best_perm, best_rmsd, best_R


def _perm_of(mapping, virtual_atoms):
    """Reconstruct the winning slot assignment from _map_to_template's mapping."""
    id_to_slot = {id(mapping[s]): s for s in range(len(mapping)) if mapping[s] is not None}
    return tuple(id_to_slot[id(a)] for a in virtual_atoms)


class TestMatcherPrefilterByteIdentity(unittest.TestCase):
    def _assert_identical(self, template, coords):
        va = _virtual_atoms(coords)
        aligner = OINDiscreteAligner(0, [])
        mapping, rmsd, R = aligner._map_to_template(va, template)
        exp_perm, exp_rmsd, exp_R = _exhaustive(va, template)
        got_perm = _perm_of(mapping, va)
        ctx = f"template n_slots={len(template)} n_atoms={len(coords)}"
        self.assertEqual(got_perm, exp_perm, f"winning permutation differs ({ctx})")
        self.assertEqual(rmsd, exp_rmsd, f"best_rmsd not bit-identical ({ctx})")
        self.assertTrue(
            np.array_equal(R.as_matrix(), exp_R), f"best rotation not bit-identical ({ctx})"
        )

    def test_battery_ideal_and_noisy(self):
        rng = np.random.default_rng(7)
        # CN 2-6 is the production range (CN7/8 exhaustive is slow and ~absent in data).
        for geo in ("LIN", "TPL", "SPL", "TET", "TPY", "TBP", "SPY", "OCT"):
            template = np.asarray(TEMPLATES[geo], dtype=float)
            m = len(template)
            for n in range(2, m + 1):
                for noise in (0.0, 0.05, 0.2):
                    slots = rng.permutation(m)[:n]
                    v = Rotation.random(random_state=rng).apply(template[slots])
                    v = v + noise * rng.standard_normal(v.shape)
                    self._assert_identical(template, v)

    def test_ideal_symmetric_ties(self):
        # Pure template slots -> many exactly-degenerate assignments; the lex-first
        # tie-break must survive the prefilter. Hardest case for byte-identity.
        for geo in ("SPL", "TET", "OCT"):
            template = np.asarray(TEMPLATES[geo], dtype=float)
            self._assert_identical(template, template.copy())


class TestMatcherPrefilterPrunes(unittest.TestCase):
    def test_candidate_set_is_a_strict_subset(self):
        template = np.asarray(TEMPLATES["OCT"], dtype=float)  # 720 permutations
        va = _virtual_atoms(template.copy())
        input_vecs = np.array([a["coords"] for a in va])
        input_norms = input_vecs / (np.linalg.norm(input_vecs, axis=1)[:, None] + 1e-9)
        perms = list(itertools.permutations(range(6), 6))
        cand = OINDiscreteAligner._candidate_permutations(perms, template, input_norms)
        self.assertLess(len(cand), len(perms), "prefilter did not prune anything")
        self.assertGreaterEqual(len(cand), 1)

    def test_candidate_set_contains_exhaustive_argmin(self):
        rng = np.random.default_rng(11)
        template = np.asarray(TEMPLATES["OCT"], dtype=float)
        for _ in range(20):
            v = Rotation.random(random_state=rng).apply(template) + 0.15 * rng.standard_normal(
                (6, 3)
            )
            va = _virtual_atoms(v)
            input_vecs = np.array([a["coords"] for a in va])
            input_norms = input_vecs / (np.linalg.norm(input_vecs, axis=1)[:, None] + 1e-9)
            perms = list(itertools.permutations(range(6), 6))
            cand = set(
                int(k)
                for k in OINDiscreteAligner._candidate_permutations(perms, template, input_norms)
            )
            exp_perm, _, _ = _exhaustive(va, template)
            exp_idx = perms.index(exp_perm)
            self.assertIn(exp_idx, cand, "exhaustive argmin fell outside the candidate set")


if __name__ == "__main__":
    unittest.main()
