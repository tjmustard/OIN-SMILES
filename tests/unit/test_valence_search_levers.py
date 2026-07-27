"""Guards for the valence-search instrumentation and its two default-OFF levers.

The point of these tests is not that the levers work -- it is that **the default path is
unchanged**. Both levers alter perceived bond orders when engaged, so the thing that must
be pinned is that an unset environment reproduces the historical constants and the
historical matcher, and that ``possible_valences`` (extracted out of ``AC2BO`` so the
over-cap population could be scanned without running the search) is a faithful extraction
rather than a re-implementation that drifted.
"""

import os
import unittest
from unittest import mock

import networkx as nx
import numpy as np

from oinsmiles.utils import perception_core as xl


class TestFallbackTriesLever(unittest.TestCase):
    def test_unset_env_is_the_historical_constant(self):
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop(xl._FALLBACK_TRIES_ENV, None)
            self.assertEqual(xl._fallback_tries(), 20_000)
            self.assertEqual(xl._VALENCE_FALLBACK_TRIES, 20_000)

    def test_empty_string_is_treated_as_unset(self):
        with mock.patch.dict(os.environ, {xl._FALLBACK_TRIES_ENV: ""}):
            self.assertEqual(xl._fallback_tries(), 20_000)

    def test_env_overrides(self):
        with mock.patch.dict(os.environ, {xl._FALLBACK_TRIES_ENV: "250"}):
            self.assertEqual(xl._fallback_tries(), 250)

    def test_garbage_and_nonpositive_fall_back_to_the_default(self):
        # A typo in an env var must not silently make the encoder try one candidate --
        # that would change every over-cap perception without anyone noticing.
        for bad in ("banana", "0", "-5", "1e4"):
            with mock.patch.dict(os.environ, {xl._FALLBACK_TRIES_ENV: bad}):
                self.assertEqual(xl._fallback_tries(), 20_000, bad)


class TestMatcherLever(unittest.TestCase):
    def _graph(self):
        g = nx.Graph()
        g.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
        return g

    def test_unset_env_is_nx_max_weight_matching(self):
        g = self._graph()
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop(xl._MATCHER_ENV, None)
            self.assertEqual(xl._maximum_matching(g), nx.max_weight_matching(g))

    def test_unknown_matcher_falls_back_to_nx(self):
        g = self._graph()
        with mock.patch.dict(os.environ, {xl._MATCHER_ENV: "not-a-matcher"}):
            self.assertEqual(xl._maximum_matching(g), nx.max_weight_matching(g))

    def test_named_arms_are_reachable_and_return_matchings(self):
        g = self._graph()
        for arm in ("nx", "maxcard", "greedy"):
            with mock.patch.dict(os.environ, {xl._MATCHER_ENV: arm}):
                matching = xl._maximum_matching(g)
            seen = set()
            for u, v in matching:
                self.assertNotIn(u, seen, arm)
                self.assertNotIn(v, seen, arm)
                seen.update((u, v))

    def test_the_graph_get_ua_pairs_builds_carries_no_weights(self):
        """The premise of Q3: this is maximum *cardinality* matching in disguise."""
        captured = {}

        def spy(g):
            captured["weights"] = [d for _, _, d in g.edges(data=True)]
            return nx.max_weight_matching(g)

        # A 4-atom chain, every atom unsaturated with one degree of unsaturation.
        AC = np.array(
            [[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]],
            dtype=int,
        )
        with mock.patch.object(xl, "_maximum_matching", spy):
            xl.get_UA_pairs([0, 1, 2, 3], AC, [1, 1, 1, 1], use_graph=True)
        self.assertTrue(captured["weights"])
        self.assertTrue(all(d == {} for d in captured["weights"]), captured["weights"])


class TestPossibleValencesExtraction(unittest.TestCase):
    """``possible_valences`` was lifted out of ``AC2BO``; pin the behaviour it must keep."""

    def test_matches_an_inline_replica_of_the_original_loop(self):
        atoms = [6, 6, 7, 8, 16, 15, 6, 8]
        AC_valence = [4, 2, 2, 1, 1, 3, 1, 2]
        got = xl.possible_valences(AC_valence, atoms, allow_carbenes=True)

        expected = []
        for atomicNum, valence in zip(atoms, AC_valence):
            pv = [x for x in xl.atomic_valence[atomicNum] if x >= valence]
            if atomicNum == 6 and valence == 1:
                pv.remove(2)
            if atomicNum == 6 and valence == 2:
                pv.append(3)
            if atomicNum == 16 and valence == 1:
                pv = [1, 2]
            expected.append(pv)
        self.assertEqual(got, expected)

    def test_allow_carbenes_false_drops_the_divalent_carbon_option(self):
        # A 2-coordinate carbon: carbenes allowed keeps 2, disallowed drops it. Both
        # append 3, which is what the original loop does.
        allowed = xl.possible_valences([2], [6], allow_carbenes=True)[0]
        disallowed = xl.possible_valences([2], [6], allow_carbenes=False)[0]
        self.assertIn(2, allowed)
        self.assertNotIn(2, disallowed)

    def test_combo_size_short_circuits_at_the_cap_like_ac2bo_does(self):
        # 30 atoms with 2 options each = 2**30, far over the cap. The reported number is
        # only required to be > cap, matching AC2BO's own early break.
        vll = [[1, 2]] * 30
        self.assertGreater(xl.valence_combo_size(vll), xl._VALENCE_COMBO_CAP)
        # Under the cap it is the exact product.
        self.assertEqual(xl.valence_combo_size([[1, 2], [1, 2, 3], [1]]), 6)


class TestStatsCounters(unittest.TestCase):
    def test_reset_zeroes_every_counter(self):
        xl.AC2BO_STATS["candidates"] = 99
        xl.reset_ac2bo_stats()
        self.assertTrue(all(v == 0 for v in xl.AC2BO_STATS.values()))

    def test_a_sub_cap_perception_counts_candidates_and_finds_a_valid_structure(self):
        # Ethene: C=C with 4 H. AC2BO should find a real Lewis structure, so the
        # over-cap counters must stay at zero.
        atoms = [6, 6, 1, 1, 1, 1]
        AC = np.zeros((6, 6), dtype=int)
        for i, j in [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]:
            AC[i, j] = AC[j, i] = 1
        xl.reset_ac2bo_stats()
        BO, _ = xl.AC2BO(AC, atoms, 0)
        self.assertEqual(xl.AC2BO_STATS["ac2bo_calls"], 1)
        self.assertEqual(xl.AC2BO_STATS["over_cap_calls"], 0)
        self.assertEqual(xl.AC2BO_STATS["over_cap_exhausted"], 0)
        self.assertEqual(xl.AC2BO_STATS["found_valid"], 1)
        self.assertEqual(BO[0, 1], 2)


if __name__ == "__main__":
    unittest.main()
