"""Guards for the over-cap valence ENUMERATION levers, and for the claim each one rests on.

Two levers, both default OFF, both read only inside ``AC2BO``'s ``if over_cap:`` branch --
which is why 99.8% of the corpus is byte-identical *by construction* rather than by sampling:

* ``OIN_VALENCE_ORDERED_FALLBACK`` -- take the bounded prefix in ``_ordered_valences``' order
  instead of the raw product order. Its correctness claim is an **equality**: the lazy
  generator must produce the sorted path's order element for element, or "sub-cap is
  untouched" stops being provable. ``test_lazy_order_is_identical_to_the_sorted_path`` is
  that proof.
* ``OIN_VALENCE_CHARGE_FILTER`` -- skip candidates that provably cannot satisfy ``AC2BO``'s
  own predicate. Its correctness claim is that the skip condition is **necessary**. That is a
  derivation about ``charge_is_OK``, and a derivation is exactly the kind of thing that is
  wrong in one unnoticed case, so it is brute-forced here instead of argued:
  ``test_every_candidate_the_real_predicate_accepts_survives_the_filter`` enumerates whole
  small spaces and asserts no accepted candidate is ever dropped.

The third thing pinned here is the **sense** of both levers. ``os.environ.get("X")`` is truthy
for the string ``"0"``, so the obvious way to opt out of a bare-truthiness lever turns it on.
``OIN_VALENCE_ORDERED_FALLBACK=0`` must mean off.
"""

import contextlib
import itertools
import os
import random
import unittest
from unittest import mock

import numpy as np

from oinsmiles.utils import xyz2mol_local as xl


def chain_ac(n, atomic_num=6):
    """Adjacency matrix of a linear ``n``-atom chain, and its atom list."""
    AC = np.zeros((n, n), dtype=int)
    for i in range(n - 1):
        AC[i, i + 1] = AC[i + 1, i] = 1
    return AC, [atomic_num] * n


def real_predicate_valid(AC, atoms, charge, valences, allow_carbenes=True):
    """Exactly what ``AC2BO``'s inner loop accepts, lifted out so a test can enumerate it."""
    AC_valence = list(AC.sum(axis=1))
    UA, DU = xl.get_UA(list(valences), AC_valence)
    if not UA:
        return bool(
            xl.BO_is_OK(
                AC,
                AC,
                charge,
                DU,
                xl.atomic_valence_electrons,
                atoms,
                list(valences),
                allow_charged_fragments=True,
                allow_carbenes=allow_carbenes,
            )
        )
    for UA_pairs in xl.get_UA_pairs(UA, AC, DU, use_graph=True):
        BO = xl.get_BO(AC, UA, DU, list(valences), UA_pairs, use_graph=True)
        if xl.BO_is_OK(
            BO,
            AC,
            charge,
            DU,
            xl.atomic_valence_electrons,
            atoms,
            list(valences),
            allow_charged_fragments=True,
            allow_carbenes=allow_carbenes,
        ):
            return True
    return False


class TestLazyOrderedValences(unittest.TestCase):
    """``iter_ordered_valences`` must reproduce ``_ordered_valences`` exactly, not roughly."""

    def test_lazy_order_is_identical_to_the_sorted_path(self):
        rng = random.Random(11)
        elements = [1, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 33, 34]
        compared = 0
        for _ in range(400):
            n = rng.randint(1, 7)
            atoms = [rng.choice(elements) for _ in range(n)]
            vll = []
            for z in atoms:
                base = list(xl.atomic_valence[z])
                lst = base[: rng.randint(1, len(base))]
                if z == 6 and rng.random() < 0.4:
                    lst = lst + [3]  # the `C valence == 2` append: [4, 2, 3]
                if rng.random() < 0.15:
                    lst = list(reversed(lst))  # exercise non-ascending lists (As is [5, 3])
                vll.append(lst)
            size = 1
            for lst in vll:
                size *= len(lst)
            if size > 5000:
                continue
            compared += 1
            self.assertEqual(
                xl._ordered_valences(vll, atoms),
                list(xl.iter_ordered_valences(vll, atoms)),
                f"atoms={atoms} vll={vll}",
            )
        # A vacuous pass is the failure mode here: if the size filter rejected everything the
        # assertion above would never run.
        self.assertGreater(compared, 150, "too few configurations actually compared")

    def test_covers_the_element_groups_the_heuristic_names(self):
        """A hand-built case touching every group plus a non-grouped, non-ascending atom."""
        atoms = [8, 7, 6, 15, 16, 33]  # O, N, C, P, S, As -- As is the tie-broken one
        vll = [[2, 1], [3, 4], [4, 2], [3, 5], [2, 4], [5, 3]]
        self.assertEqual(
            xl._ordered_valences(vll, atoms), list(xl.iter_ordered_valences(vll, atoms))
        )

    def test_first_candidate_comes_from_an_unmaterialisable_space(self):
        """The whole point: no exponential materialisation before the first yield."""
        atoms = [6] * 200
        vll = [[4, 2, 3]] * 200
        first = next(iter(xl.iter_ordered_valences(vll, atoms)))
        self.assertEqual(first, tuple([4] * 200))

    def test_lazy_order_matches_on_CORPUS_SHAPED_configurations(self):
        """The 400 random configs above are 1-7 atoms; real ligand fragments are not.

        ``iter_ordered_valences`` is now the **sub-cap** enumeration (99.8% of the corpus),
        so the equality it rests on has to hold on the population the corpus actually emits,
        not only on the one the generator above happens to produce. That population was
        harvested from 15 slow / high-``combo_size`` molecules (50 ligand fragments across
        both ``allow_carbenes`` arms, incl. ``KEMTED`` 168 atoms and ``XIRMER`` 108 atoms):

        * fragments are 2-168 atoms, but only **1-12** atoms have more than one candidate
          valence -- which is why a 168-atom fragment can have ``combo_size`` 16;
        * every per-atom option list has length **1, 2 or 3**, never more;
        * in every fragment measured, the multi-option atoms were C/N/O/P/S, i.e. exactly
          the elements the heuristic groups.

        That last point is why this test does not simply reuse the generator above: when
        every multi-option atom is grouped, ``order_idx`` pins the whole candidate and
        ``sorted()``'s lexicographic tie-break is never exercised. So the corpus shape and
        the tie-break shape are *different* risks, and both are generated here -- the
        ``forced_tiebreak`` arm puts the options on a NON-grouped element (As/Cl/Se) at
        large ``n``, which is the case the corpus does not currently show but the code
        still has to get right.
        """
        rng = random.Random(4711)
        grouped = [6, 7, 8, 15, 16]
        ungrouped = [33, 17, 34, 9, 35]
        filler = [1, 1, 1, 6, 6, 9, 17]
        compared = 0
        for trial in range(60):
            forced_tiebreak = trial % 3 == 2
            n = rng.randint(20, 200)
            atoms = [rng.choice(filler) for _ in range(n)]
            vll = [[xl.atomic_valence[z][0]] for z in atoms]
            k = rng.randint(1, 12)
            size = 1
            for idx in rng.sample(range(n), k):
                z = rng.choice(ungrouped if forced_tiebreak else grouped)
                base = list(xl.atomic_valence[z])
                if len(base) < 2:
                    base = base + [base[0] + 2]
                opts = base[: rng.randint(2, min(3, len(base)))]
                if rng.random() < 0.2:
                    opts = list(reversed(opts))  # As is [5, 3]: non-ascending is real
                if size * len(opts) > 20000:
                    continue
                size *= len(opts)
                atoms[idx] = z
                vll[idx] = opts
            if size == 1:
                continue  # a single candidate proves nothing about ORDER
            compared += 1
            self.assertEqual(
                xl._ordered_valences(vll, atoms),
                list(xl.iter_ordered_valences(vll, atoms)),
                f"n={n} combo={size} forced_tiebreak={forced_tiebreak}",
            )
        self.assertGreater(compared, 30, "too few corpus-shaped configurations compared")

    def test_subcap_AC2BO_does_not_materialise_the_product(self):
        """Structural pin: the sub-cap branch must take the LAZY path.

        ``_ordered_valences`` builds the full Cartesian product twice plus a dict over the
        five-group product, and the candidate loop that consumes it measurably exits after
        a handful of iterations -- on ``NOCGAN_comp_0`` one ``AC2BO`` call consumed **1**
        candidate out of a materialised 20 736 (99.2% of that call's wall). This asserts
        the dead work is gone rather than trusting a code reading: if anything reintroduces
        the eager call on the sub-cap branch, ``_ordered_valences`` raises and this fails.
        """
        AC, atoms = chain_ac(6)
        atoms = [6, 7, 6, 8, 16, 6]

        def boom(*_a, **_kw):
            raise AssertionError("sub-cap AC2BO materialised the full product")

        with mock.patch.object(xl, "_ordered_valences", boom):
            BO, _ = xl._AC2BO_core(AC, atoms, 0, allow_charged_fragments=True, use_graph=True)
        self.assertEqual(BO.shape, AC.shape)


class TestChargeFilterIsNecessary(unittest.TestCase):
    """The filter may drop candidates only if none of them could ever have been accepted."""

    def _spaces(self):
        """Small (AC, atoms, charge) cases whose whole candidate space can be enumerated."""
        cases = []
        # Ethene: a valid Lewis structure exists, so these are non-vacuous.
        AC = np.zeros((6, 6), dtype=int)
        for i, j in [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]:
            AC[i, j] = AC[j, i] = 1
        cases.append((AC, [6, 6, 1, 1, 1, 1], 0))
        # Acetate-like: C-C(-O)(-O), charges 0 and -1 both probed.
        AC = np.zeros((4, 4), dtype=int)
        for i, j in [(0, 1), (1, 2), (1, 3)]:
            AC[i, j] = AC[j, i] = 1
        for charge in (0, -1, -2):
            cases.append((AC, [6, 6, 8, 8], charge))
        # A carbon chain with heteroatoms, several charges.
        AC, atoms = chain_ac(6)
        atoms = [6, 7, 6, 8, 16, 6]
        for charge in (0, -1, 1, -3):
            cases.append((AC, atoms, charge))
        # A 5-membered ring of N/C.
        AC = np.zeros((5, 5), dtype=int)
        for i in range(5):
            AC[i, (i + 1) % 5] = AC[(i + 1) % 5, i] = 1
        for charge in (0, -1, -2):
            cases.append((AC, [6, 7, 6, 7, 6], charge))
        return cases

    def test_every_candidate_the_real_predicate_accepts_survives_the_filter(self):
        accepted_total = 0
        for AC, atoms, charge in self._spaces():
            AC_valence = list(AC.sum(axis=1))
            vll = xl.possible_valences(AC_valence, atoms, allow_carbenes=True)
            kept = set(xl.iter_charge_feasible_valences(vll, atoms, charge, AC_valence))
            for cand in itertools.product(*vll):
                if real_predicate_valid(AC, atoms, charge, cand):
                    accepted_total += 1
                    self.assertIn(
                        cand,
                        kept,
                        f"filter dropped an ACCEPTED candidate: atoms={atoms} "
                        f"charge={charge} valences={cand}",
                    )
        # If no case in the corpus above had a valid structure, the loop asserts nothing.
        self.assertGreater(accepted_total, 0, "no accepted candidate anywhere -- test vacuous")

    def test_output_is_a_subsequence_of_the_raw_product_in_the_same_order(self):
        """It skips candidates; it must never reorder them."""
        for AC, atoms, charge in self._spaces():
            AC_valence = list(AC.sum(axis=1))
            vll = xl.possible_valences(AC_valence, atoms, allow_carbenes=True)
            kept = list(xl.iter_charge_feasible_valences(vll, atoms, charge, AC_valence))
            raw = list(itertools.product(*vll))
            positions = [raw.index(c) for c in kept]
            self.assertEqual(positions, sorted(positions), f"reordered: atoms={atoms}")
            self.assertEqual(len(set(kept)), len(kept), "duplicate candidate emitted")

    def test_an_unreachable_charge_yields_nothing_at_all(self):
        """LIYFAA's shape: no candidate can meet the target, so the search is hopeless."""
        AC, atoms = chain_ac(4)
        AC_valence = list(AC.sum(axis=1))
        vll = xl.possible_valences(AC_valence, atoms, allow_carbenes=True)
        kept = list(xl.iter_charge_feasible_valences(vll, atoms, +99, AC_valence))
        self.assertEqual(kept, [])

    def test_it_actually_prunes(self):
        """A filter that keeps everything would pass every test above and buy nothing."""
        AC, atoms = chain_ac(8)
        AC_valence = list(AC.sum(axis=1))
        vll = xl.possible_valences(AC_valence, atoms, allow_carbenes=True)
        space = 1
        for lst in vll:
            space *= len(lst)
        kept = list(xl.iter_charge_feasible_valences(vll, atoms, 0, AC_valence))
        self.assertLess(len(kept), space)
        self.assertGreater(len(kept), 0)


class TestInfeasibleChargeFallsBackToTheHistoricalPath(unittest.TestCase):
    """The lever must never turn one guess into a *different* guess.

    When no candidate can satisfy the charge, there is no valid structure to find -- but
    ``best_BO`` is still what downstream judges, and it is assembled from candidates the
    filter drops. So the over-cap branch falls back to the raw product in that case, which
    is what makes the lever's blast radius exactly "a guess becomes a real Lewis structure".
    LIYFAA_comp_0 is the real molecule this covers: 0 feasible candidates at charge -10.
    """

    def _call(self, charge, env):
        AC, atoms = chain_ac(17)
        env = dict(env)
        env.setdefault(xl._FALLBACK_TRIES_ENV, "5")
        with mock.patch.dict(os.environ, env, clear=False):
            for name in (xl._ORDERED_FALLBACK_ENV, xl._CHARGE_FILTER_ENV):
                if name not in env:
                    os.environ.pop(name, None)
            xl.reset_ac2bo_stats()
            BO, _ = xl.AC2BO(AC, atoms, charge)
        return BO, dict(xl.AC2BO_STATS)

    def test_an_infeasible_charge_reproduces_the_default_answer_exactly(self):
        # +99 on a 17-carbon chain is unreachable, so the feasible set is empty.
        base, base_stats = self._call(99, {})
        BO, stats = self._call(99, {xl._CHARGE_FILTER_ENV: "1"})
        self.assertEqual(stats["over_cap_infeasible"], 1)
        self.assertEqual(base_stats["over_cap_infeasible"], 0)
        self.assertTrue(np.array_equal(BO, base), "infeasible fallback changed best_BO")
        self.assertEqual(stats["candidates"], base_stats["candidates"])


class TestChargeFilterDeclinesWhatItCannotReasonAbout(unittest.TestCase):
    """A fragment containing a transition metal must not be made *worse* by the lever.

    All 30 transition metals have an ``atomic_valence`` entry (``[20]``) and **no**
    ``atomic_valence_electrons`` entry, so ``get_atomic_charge``'s call signature is a
    ``KeyError`` for them. The filter would hit that before examining a single candidate.
    """

    def test_the_gap_this_guard_exists_for_is_real(self):
        missing = [z for z in xl.atomic_valence if z not in xl.atomic_valence_electrons]
        self.assertTrue(missing, "no element lacks a valence-electron count any more")
        self.assertFalse(xl.charge_filter_supported([missing[0], 6]))
        self.assertTrue(xl.charge_filter_supported([6, 7, 8, 1]))

    def test_a_metal_bearing_fragment_takes_the_historical_path(self):
        AC, atoms = chain_ac(17)
        atoms = list(atoms)
        atoms[0] = next(z for z in xl.atomic_valence if z not in xl.atomic_valence_electrons)
        env = {xl._CHARGE_FILTER_ENV: "1", xl._FALLBACK_TRIES_ENV: "3"}
        with mock.patch.dict(os.environ, env, clear=False):
            os.environ.pop(xl._ORDERED_FALLBACK_ENV, None)
            xl.reset_ac2bo_stats()
            with contextlib.suppress(KeyError):
                # The default path raises this too, from charge_is_OK. What is asserted is
                # only that the filter declined rather than crashing earlier than the default.
                xl.AC2BO(AC, atoms, 0)
            stats = dict(xl.AC2BO_STATS)
        self.assertEqual(stats["over_cap_filter_unsupported"], 1)
        self.assertEqual(stats["over_cap_filtered_calls"], 0)


class TestOverCapLeverDefaults(unittest.TestCase):
    """Both levers OFF by default, and ``"0"`` means off."""

    def _over_cap_call(self, env):
        # 17-atom carbon chain: the 15 interior atoms have 3 candidate valences each,
        # 3**15 = 14 348 907, comfortably over _VALENCE_COMBO_CAP.
        AC, atoms = chain_ac(17)
        env = dict(env)
        env.setdefault(xl._FALLBACK_TRIES_ENV, "5")  # keep the test fast
        with mock.patch.dict(os.environ, env, clear=False):
            for name in (xl._ORDERED_FALLBACK_ENV, xl._CHARGE_FILTER_ENV):
                if name not in env:
                    os.environ.pop(name, None)
            xl.reset_ac2bo_stats()
            BO, _ = xl.AC2BO(AC, atoms, 0)
        return BO, dict(xl.AC2BO_STATS)

    def test_the_branch_under_test_is_actually_the_over_cap_one(self):
        _, stats = self._over_cap_call({})
        self.assertEqual(stats["over_cap_calls"], 1, "test case is not over cap any more")

    def test_unset_environment_takes_neither_lever(self):
        _, stats = self._over_cap_call({})
        self.assertEqual(stats["over_cap_ordered_calls"], 0)
        self.assertEqual(stats["over_cap_filtered_calls"], 0)

    def test_zero_disables_rather_than_enabling(self):
        for name in (xl._ORDERED_FALLBACK_ENV, xl._CHARGE_FILTER_ENV):
            for falsey in ("0", "", "false", "no", "off"):
                _, stats = self._over_cap_call({name: falsey})
                self.assertEqual(stats["over_cap_ordered_calls"], 0, f"{name}={falsey!r}")
                self.assertEqual(stats["over_cap_filtered_calls"], 0, f"{name}={falsey!r}")

    def test_each_lever_is_reachable_when_enabled(self):
        _, stats = self._over_cap_call({xl._ORDERED_FALLBACK_ENV: "1"})
        self.assertEqual(stats["over_cap_ordered_calls"], 1)
        _, stats = self._over_cap_call({xl._CHARGE_FILTER_ENV: "1"})
        self.assertEqual(stats["over_cap_filtered_calls"], 1)

    def test_the_filter_wins_when_both_are_set(self):
        _, stats = self._over_cap_call({xl._ORDERED_FALLBACK_ENV: "1", xl._CHARGE_FILTER_ENV: "1"})
        self.assertEqual(stats["over_cap_filtered_calls"], 1)
        self.assertEqual(stats["over_cap_ordered_calls"], 0)


class TestSubCapIsUntouched(unittest.TestCase):
    """The 99.8% claim, checked functionally as well as structurally."""

    def _sub_cap(self, env):
        # Ethene -- sub cap by a wide margin.
        atoms = [6, 6, 1, 1, 1, 1]
        AC = np.zeros((6, 6), dtype=int)
        for i, j in [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]:
            AC[i, j] = AC[j, i] = 1
        with mock.patch.dict(os.environ, env, clear=False):
            for name in (xl._ORDERED_FALLBACK_ENV, xl._CHARGE_FILTER_ENV):
                if name not in env:
                    os.environ.pop(name, None)
            xl.reset_ac2bo_stats()
            BO, _ = xl.AC2BO(AC, atoms, 0)
        return BO, dict(xl.AC2BO_STATS)

    def test_neither_lever_can_change_a_sub_cap_perception(self):
        base, base_stats = self._sub_cap({})
        self.assertEqual(base_stats["over_cap_calls"], 0)
        for env in (
            {xl._ORDERED_FALLBACK_ENV: "1"},
            {xl._CHARGE_FILTER_ENV: "1"},
            {xl._ORDERED_FALLBACK_ENV: "1", xl._CHARGE_FILTER_ENV: "1"},
        ):
            BO, stats = self._sub_cap(env)
            self.assertTrue(np.array_equal(BO, base), env)
            self.assertEqual(stats["candidates"], base_stats["candidates"], env)
            self.assertEqual(stats["over_cap_ordered_calls"], 0, env)
            self.assertEqual(stats["over_cap_filtered_calls"], 0, env)


if __name__ == "__main__":
    unittest.main()
