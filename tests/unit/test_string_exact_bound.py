"""Tests for the incumbent search bound (v0.4.16 Lane 1, ``OIN_STRING_EXACT_BOUND``).

``OIN_ACCEPT_STRING_EXACT`` (v0.4.15) buys +48 molecules / +0.96 pts with zero losses, and pays
**4.00x** runtime for them: the predicate declines to STOP the pool, so the pool fills to budget.
Measured on the frozen arm, the 317 molecules that never gain consume 16149 s of the lever's
17191 s -- **93.9% of the entire bill buys nothing.** This bound caps that tail.

What these tests pin, in the order the risk actually sits:

1. **bound=0 is byte-identical to the pre-lever answer.** This is the wiring gate, and it is a
   test rather than only a measurement because it is the one property that makes every derived
   number trustworthy. A broken bound and a working bound print the same recovered count at large
   N; they differ here.
2. **The off-by-one.** ``bound=N`` must allow exactly N ``accept_fn`` evaluations *after* the one
   that recorded the incumbent. The first draft of this loop counted the recording evaluation
   itself, which silently shifts every point on the knee curve by one.
3. **The bound cannot fire before an incumbent exists.** Until then there is nothing to fall back
   on, and stopping early would hand the caller the energy-sorted pool instead -- a real
   regression, and exactly the case the derivation must exclude rather than assume away.
4. **Unset is unbounded**, so every pre-v0.4.16 caller is untouched.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.metallogen_adapter import convert_parsed_to_msmiles
from oinsmiles.generation.oin_parser import OINParser
from oinsmiles.generator3d import ACCEPT_INCUMBENT, generate_3d_structures
from oinsmiles.oin.levers import held_off, lever_int

SEED = 42
CISPLATIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"


def _msmiles(oin):
    return convert_parsed_to_msmiles(OINParser().parse(oin))


class _CountingAcceptFn:
    """An ``accept_fn`` that records how many times it was consulted.

    ``verdicts`` is consumed in order; the last value repeats once exhausted. That lets a test
    say "incumbent, then incumbent forever" or "incumbent, then a string-exact hit on the 3rd
    call" without reaching into the loop.
    """

    def __init__(self, *verdicts):
        self.verdicts = list(verdicts)
        self.calls = 0

    def __call__(self, _mol):
        i = min(self.calls, len(self.verdicts) - 1)
        self.calls += 1
        return self.verdicts[i]


def _positions(mols):
    """A comparable fingerprint of what the generator handed back.

    ``get_coordinate_list`` and not ``get_position``: these are MetalloGen ``chem.Molecule``
    objects, the same ones the pipeline's standing note warns are not ``rdkit.Chem.Mol``.
    """
    return [tuple(tuple(round(c, 9) for c in xyz) for xyz in m.get_coordinate_list()) for m in mols]


class TestBoundZeroIsThePreLeverAnswer(unittest.TestCase):
    """🔴 THE WIRING GATE. bound=0 must reproduce the lever-OFF answer byte-for-byte.

    Lever OFF, the predicate returns ``True`` on the key match and the pool stops there. Lever ON
    with bound=0, the same conformer comes back as the recorded incumbent. Same geometry, one
    conformer, same number of embeds -- the difference is only which return path carried it.

    If this test fails, every recovered-vs-bound number derived from the instrumented run is
    void, because the zero point of the curve does not sit where the OFF arm sits.
    """

    def test_bound_zero_matches_accept_true(self):
        msmiles = _msmiles(CISPLATIN)

        off = _CountingAcceptFn(True)
        mols_off = generate_3d_structures(
            msmiles, num_conformers=5, uff_pool_size=10, seed=SEED, accept_fn=off
        )
        on = _CountingAcceptFn(ACCEPT_INCUMBENT)
        mols_on = generate_3d_structures(
            msmiles,
            num_conformers=5,
            uff_pool_size=10,
            seed=SEED,
            accept_fn=on,
            incumbent_bound=0,
        )

        self.assertEqual(len(mols_off), 1, "lever-OFF stops on the first accepted conformer")
        self.assertEqual(len(mols_on), 1, "bound=0 returns the incumbent as the sole pool member")
        self.assertEqual(
            _positions(mols_off),
            _positions(mols_on),
            "bound=0 must be BYTE-IDENTICAL to the pre-lever answer -- this is the wiring gate",
        )
        self.assertEqual(
            off.calls, on.calls, "bound=0 consults the predicate exactly as often as lever-OFF"
        )


class TestTheOffByOne(unittest.TestCase):
    """``bound=N`` allows exactly N evaluations AFTER the incumbent was recorded."""

    def _calls_at_bound(self, bound):
        fn = _CountingAcceptFn(ACCEPT_INCUMBENT)
        generate_3d_structures(
            _msmiles(CISPLATIN),
            num_conformers=5,
            uff_pool_size=10,
            seed=SEED,
            accept_fn=fn,
            incumbent_bound=bound,
        )
        return fn.calls

    def test_each_bound_adds_exactly_one_evaluation(self):
        # The incumbent is recorded on call 1 (cisplatin's first conformer is always accepted by
        # this stub), so bound=N costs N+1 calls. Asserted as a sequence rather than a single
        # point: an off-by-one is invisible at any one bound and obvious across three.
        self.assertEqual(
            [self._calls_at_bound(n) for n in (0, 1, 2)],
            [1, 2, 3],
            "bound=N must permit N evaluations BEYOND the one that recorded the incumbent",
        )

    def test_unset_bound_is_unbounded(self):
        # None is the v0.4.15 behaviour and must stay the default: the pool fills past the
        # incumbent rather than stopping at it.
        fn = _CountingAcceptFn(ACCEPT_INCUMBENT)
        generate_3d_structures(
            _msmiles(CISPLATIN),
            num_conformers=5,
            uff_pool_size=10,
            seed=SEED,
            accept_fn=fn,
            incumbent_bound=None,
        )
        self.assertGreater(fn.calls, 3, "unset must not bound the search")


class TestTheBoundCannotFireWithoutAnIncumbent(unittest.TestCase):
    """Until an incumbent exists there is nothing to fall back on, so the bound must not stop."""

    def test_all_false_verdicts_fill_the_pool_even_at_bound_zero(self):
        # An accept_fn that never returns ACCEPT_INCUMBENT leaves the generator on its ordinary
        # energy-sorted-pool path. bound=0 must be inert here -- stopping would return a SHORTER
        # pool to `_select_by_geometry`, which can then select differently. This is the one case
        # where truncation is not answer-neutral, and it is why the derivation excludes
        # no-incumbent molecules instead of assuming they behave.
        fn = _CountingAcceptFn(False)
        mols = generate_3d_structures(
            _msmiles(CISPLATIN),
            num_conformers=5,
            uff_pool_size=10,
            seed=SEED,
            accept_fn=fn,
            incumbent_bound=0,
        )
        self.assertGreater(len(mols), 1, "a bound with no incumbent recorded must not truncate")


class TestLeverIntSemantics(unittest.TestCase):
    """``lever_int`` exists because 0 is a meaningful bound and must survive as 0."""

    def setUp(self):
        self._prior = os.environ.get("OIN_STRING_EXACT_BOUND")

    def tearDown(self):
        if self._prior is None:
            os.environ.pop("OIN_STRING_EXACT_BOUND", None)
        else:
            os.environ["OIN_STRING_EXACT_BOUND"] = self._prior

    def test_zero_is_zero_not_unset(self):
        # 🔴 The whole reason this reader exists. `os.environ.get(name) or default` would turn
        # "0" into the default, i.e. turn the tightest bound into NO bound -- the loudest
        # possible wrong answer, and the same shape as the truthiness trap `lever_enabled` fixed.
        os.environ["OIN_STRING_EXACT_BOUND"] = "0"
        self.assertEqual(lever_int("OIN_STRING_EXACT_BOUND"), 0)
        self.assertIsNotNone(lever_int("OIN_STRING_EXACT_BOUND"))

    def test_unset_is_none(self):
        os.environ.pop("OIN_STRING_EXACT_BOUND", None)
        self.assertIsNone(lever_int("OIN_STRING_EXACT_BOUND"))

    def test_garbage_and_negatives_fall_back_rather_than_raise(self):
        for raw in ("", "  ", "abc", "1.5", "-1"):
            with self.subTest(raw=raw):
                os.environ["OIN_STRING_EXACT_BOUND"] = raw
                self.assertIsNone(lever_int("OIN_STRING_EXACT_BOUND"))

    def test_whitespace_is_tolerated(self):
        os.environ["OIN_STRING_EXACT_BOUND"] = " 7 "
        self.assertEqual(lever_int("OIN_STRING_EXACT_BOUND"), 7)


class TestTheBoundIsCoupledToTheLever(unittest.TestCase):
    """Modelled on ``test_levers::TestDonorFoldAndParityVetoAreCoupled``.

    Promoting ``OIN_ACCEPT_STRING_EXACT`` without the bound reinstates the full 4.00x, which is
    the cost that kept it default-OFF in the first place. The coupling has to be stated where a
    future promoter will read it, not only in a release note.
    """

    BOUND = "OIN_STRING_EXACT_BOUND"
    LEVER = "OIN_ACCEPT_STRING_EXACT"

    def test_both_are_held_off_together(self):
        off = held_off()
        self.assertIn(self.BOUND, off)
        self.assertIn(self.LEVER, off)

    def test_the_bound_states_the_coupling(self):
        why = held_off()[self.BOUND]
        self.assertIn(self.LEVER, why, "the bound must name the lever it is coupled to")
        self.assertIn("COUPLED", why)

    def test_the_bound_states_why_lever_enabled_is_wrong_for_it(self):
        why = held_off()[self.BOUND]
        self.assertIn("lever_int", why)
        self.assertIn("0 is a meaningful", why)


if __name__ == "__main__":
    unittest.main()
