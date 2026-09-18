"""``OIN_EXACT_DONOR_FOLD`` -- fold over true fragment automorphisms, not symmetry buckets.

Built and promoted to default-ON in v0.4.17.

WHAT IS BEING PINNED, AND WHY EACH TEST EXISTS
==============================================
The bucket fold (``_donor_swap_permutations``) permutes every symmetry class of a fragment
independently. A ligand's C2 axis moves all of its classes at once, so a one-class swap is not a
symmetry of the ligand: it relabels the complex as a DIFFERENT arrangement. The exact fold
(``_donor_automorphism_permutations``) returns only the slot permutations a real automorphism
induces. Four properties, and only the first two are what the lever is *for*:

1. an achiral complex and its mirror image get ONE string, with no mirror encode;
2. a complex that is chiral only through how a ligand WRAPS (cis-alpha) keeps two strings --
   and the bucket fold demonstrably does not, which is the defect in one fixture;
3. unset means ON, and ``=0`` still returns the bucket fold and its veto;
4. the emitted string is invariant to which member of its own orbit it was handed.

The fixtures are rotation-only strings of real cohort molecules (census ``veto_probe.jsonl``,
``fold_off`` arm), x and its z-mirror, so each test is the corpus measurement in miniature.

⚠ Coverage is five ligand motifs. The corpus instrument is ``tools/v0417/autofold_audit.py``.
"""

import os
import unittest
from unittest import mock

from oinsmiles.oin import fold_parity
from oinsmiles.oin.canonical_slots import (
    _donor_automorphism_permutations,
    _donor_swap_permutations,
    canonicalize_oin_slots,
)
from oinsmiles.oin.compare import _parse_vertex_colors, normalize_oin_for_comparison
from oinsmiles.oin.levers import default_on, held_off, lever_enabled

LEVER = "OIN_EXACT_DONOR_FOLD"

#: ``(name, rotation-only E(x), rotation-only E(mirror x))``
ACHIRAL = [
    (
        "ACEPUT 4,7-dichlorophen",
        "[Re_OCT].Clc1ccn{0}c2c1ccc1ccc3c(Cl)ccn{2}c3c21.C{1}#O.C{3}#O.C{4}#O.[Br]{5}",
        "[Re_OCT].Clc1ccn{2}c2c1ccc1ccc3c(Cl)ccn{0}c3c21.C{1}#O.C{3}#O.C{4}#O.[Br]{5}",
    ),
    (
        "ADABIT terpy outer N",
        "[Ru_OCT].O{1}C(=O)c1ccc2ccc3cccn{0}c3c2n{2}1.c1ccc(-c2cccc(-c3ccccn{5}3)n{3}2)n{4}c1",
        "[Ru_OCT].O{1}C(=O)c1ccc2ccc3cccn{0}c3c2n{2}1.c1ccc(-c2cccc(-c3ccccn{4}3)n{3}2)n{5}c1",
    ),
    (
        # One pyrrole perceived as a diradical. CanonicalRankAtoms ignores radicals, the
        # substructure matcher does not: without zeroing them the automorphism is never found.
        "LAMTAX A2B corrole, [CH][CH] radicals",
        "[Co_SPY].c1ccc(P{0}(c2ccccc2)c2ccccc2)cc1.O=[N+]([O-])c1cccc(C2=C3C=CC(=N{1}3)C(c3c(Cl)cccc3Cl)=C3C=CC(=N{4}3)C(c3cccc([N+](=O)[O-])c3)=C3C=CC(=N{2}3)C3[CH][CH]C2N{3}3)c1",
        "[Co_SPY].c1ccc(P{0}(c2ccccc2)c2ccccc2)cc1.O=[N+]([O-])c1cccc(C2=C3C=CC(=N{1}3)C(c3c(Cl)cccc3Cl)=C3C=CC(=N{3}3)C(c3cccc([N+](=O)[O-])c3)=C3C=CC(=N{2}3)C3[CH][CH]C2N{4}3)c1",
    ),
    (
        # The two strings are tmeda written forwards and backwards: one molecule by any reading.
        # The shipped veto "protects" this pair; there is nothing in the strings to protect.
        "PEBVEZ tmeda",
        "[Zn_TET].CN{0}(C)CCN{1}(C)C.CN(C)Cc1ccccc1[CH]{2}[Si](C)(C)C.[Cl]{3}",
        "[Zn_TET].CN{1}(C)CCN{0}(C)C.CN(C)Cc1ccccc1[CH]{2}[Si](C)(C)C.[Cl]{3}",
    ),
]

#: Chiral by WRAP: phenolates trans, ethers cis -- a C2-symmetric cis-alpha tetradentate. The
#: mirror differs by exchanging the ether pair ALONE, which is not a ligand automorphism.
WRAP_CHIRAL = [
    (
        "VOLGOV Hf bis(phenolate)-bis(ether)",
        "[Hf_OCT].O{0}c1c(-c2ccccc2)cccc1-c1ccccc1O{4}CCCO{2}c1ccccc1-c1cccc(-c2ccccc2)c1O{1}.[CH3]{3}.[CH3]{5}",
        "[Hf_OCT].O{0}c1c(-c2ccccc2)cccc1-c1ccccc1O{2}CCCO{4}c1ccccc1-c1cccc(-c2ccccc2)c1O{1}.[CH3]{3}.[CH3]{5}",
    ),
    (
        "GEFGUW helical PNNNP",
        "[Mn_OCT].CC(=N{0}CCCP{4}(c1ccccc1)c1ccccc1)c1cccc(C(C)=N{1}CCCP{5}(c2ccccc2)c2ccccc2)n{2}1.[H]{3}",
        "[Mn_OCT].CC(=N{0}CCCP{5}(c1ccccc1)c1ccccc1)c1cccc(C(C)=N{1}CCCP{4}(c2ccccc2)c2ccccc2)n{2}1.[H]{3}",
    ),
]


def _env(exact):
    # Written as "1"/"0", never deleted: `lever_enabled` reads "0" as off and unset as default.
    return mock.patch.dict(os.environ, {LEVER: "1" if exact else "0"})


def _perms(fn, oin):
    frags = [f for f in oin.split(".") if f]
    _metal, _geo, vcolor = _parse_vertex_colors(normalize_oin_for_comparison(oin))
    return fn(frags, vcolor)


class TestLeverIsPromoted(unittest.TestCase):
    """Default-ON since v0.4.17 (owner decision 2026-09-18). Was ``TestLeverIsHeldOff``.

    Inverted rather than deleted, for the reason ``test_fold_parity.TestDefaultOn`` gives: the
    property worth pinning is that UNSET resolves to the shipped answer and ``=0`` still disables.
    """

    def test_default_on_and_no_longer_held_off(self):
        self.assertIn(LEVER, default_on())
        self.assertNotIn(LEVER, held_off())

    def test_unset_means_on(self):
        env = {k: v for k, v in os.environ.items() if k != LEVER}
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertTrue(lever_enabled(LEVER))

    def test_explicit_zero_disables(self):
        """Load-bearing: ``=0`` is the one spelling that brings the bucket fold AND its veto back."""
        with mock.patch.dict(os.environ, {LEVER: "0"}):
            self.assertFalse(lever_enabled(LEVER))


class TestAchiralPairsGetOneString(unittest.TestCase):
    def test_exact_fold_unifies_x_and_mirror(self):
        with _env(True):
            for name, x, m in ACHIRAL:
                with self.subTest(name):
                    self.assertEqual(canonicalize_oin_slots(x), canonicalize_oin_slots(m))


class TestWrapChiralityIsNotFolded(unittest.TestCase):
    def test_exact_fold_keeps_the_pair_apart(self):
        with _env(True):
            for name, x, m in WRAP_CHIRAL:
                with self.subTest(name):
                    self.assertNotEqual(canonicalize_oin_slots(x), canonicalize_oin_slots(m))

    def test_bucket_fold_collapses_the_same_pair(self):
        """The defect, pinned: without this the test above could pass on a fold that never fires."""
        with _env(False):
            for name, x, m in WRAP_CHIRAL:
                with self.subTest(name):
                    self.assertEqual(canonicalize_oin_slots(x), canonicalize_oin_slots(m))

    def test_a_one_class_swap_is_not_an_automorphism(self):
        _name, x, _m = WRAP_CHIRAL[0]
        bucket = _perms(_donor_swap_permutations, x)
        exact = _perms(_donor_automorphism_permutations, x)
        self.assertEqual(len(bucket), 4)  # {phenolates} x {ethers}, each swapped or not
        self.assertEqual(len(exact), 2)  # identity, and the C2 axis moving BOTH pairs
        moved = [p for p in exact if any(k != v for k, v in p.items())]
        self.assertEqual(len(moved), 1)
        self.assertEqual({k for k, v in moved[0].items() if k != v}, {0, 1, 2, 4})


class TestUnsetIsTheExactFold(unittest.TestCase):
    """Was ``TestOffIsByteIdentical`` (unset == "0") while the lever was held off.

    After promotion that equality is FALSE by design, and a test still asserting it would be
    asserting the old default. What must hold now: unset is byte-identical to "1", and "0" is
    still reachable and still the bucket fold -- it differs on exactly the fixtures the fold
    fires on, so this cannot pass on a lever that stopped reaching the post-pass.
    """

    def test_unset_agrees_with_one_on_every_fixture(self):
        env = {k: v for k, v in os.environ.items() if k != LEVER}
        for _name, x, m in ACHIRAL + WRAP_CHIRAL:
            for s in (x, m):
                with _env(True):
                    on = canonicalize_oin_slots(s)
                with mock.patch.dict(os.environ, env, clear=True):
                    self.assertEqual(on, canonicalize_oin_slots(s))

    def test_zero_still_selects_the_bucket_fold(self):
        _name, x, m = WRAP_CHIRAL[0]
        with _env(False):
            self.assertEqual(canonicalize_oin_slots(x), canonicalize_oin_slots(m))
        with _env(True):
            self.assertNotEqual(canonicalize_oin_slots(x), canonicalize_oin_slots(m))


class TestOrbitInvariance(unittest.TestCase):
    def test_idempotent_and_blind_to_the_bucket_relabeling(self):
        """The shipped string is another member of the input's orbit -- when it is one at all."""
        with _env(True):
            for name, x, m in ACHIRAL:
                with self.subTest(name):
                    out = canonicalize_oin_slots(x)
                    self.assertEqual(canonicalize_oin_slots(out), out)
                    self.assertEqual(canonicalize_oin_slots(m), out)


class TestResolveSkipsTheVeto(unittest.TestCase):
    def test_no_mirror_is_built_and_the_branch_is_observable(self):
        _name, x, _m = ACHIRAL[0]
        with _env(True), mock.patch.object(fold_parity, "_encode_pairs") as enc:
            fold_parity._state.outcome = None  # last_outcome() is sticky per thread
            out = fold_parity.resolve(x, None, None)
            self.assertEqual(out, canonicalize_oin_slots(x))
            self.assertEqual(fold_parity.last_outcome(), "exact_fold")
            enc.assert_not_called()

    def test_with_the_lever_off_the_veto_path_is_untouched(self):
        _name, x, _m = ACHIRAL[1]  # the fold fires on this presentation
        with _env(False):
            fold_parity._state.outcome = None
            fold_parity.resolve(x, None, None)
            self.assertEqual(fold_parity.last_outcome(), "declined_no_conformer")


if __name__ == "__main__":
    unittest.main()
