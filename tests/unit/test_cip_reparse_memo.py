"""Guards for the CIP re-parse memo (v0.4.10, ``OIN_MEMO_CIP_REPARSE``).

``chirality._reparse_cip_label_once`` is 2.43 s a call and 99% of ``VAFMIA_comp_0``'s
generation, and the in-flight one of these is what sets v0.4.9's budget epsilon. v0.4.10
memoises it, and v0.4.10 is a **byte-identical by construction** release -- so the memo's
whole licence to exist is the claim that a hit returns exactly what a miss would compute.

These tests pin the four properties that claim rests on:

1. **Warm matches cold.** For every key shape (plain sp3, aromatic-adjacent, open-valence
   donor, non-stereocentre), the memoised answer equals the un-memoised one.
2. **The key is complete.** ``fill_deficit`` genuinely changes the answer for at least one
   input, so a memo that dropped it from the key would be caught here rather than in a
   corpus sweep. Same for ``probe``.
3. **The lever's sense is right.** Default OFF means *no memo traffic at all*, and
   ``OIN_MEMO_CIP_REPARSE=0`` disables rather than enables -- the ``"0"``-is-truthy trap
   that cost this project 23 test failures across two promotions.
4. **It actually hits.** A memo that never hits is dead weight plus a risk. The repeat
   traffic this lever targets is structural (``accept_fn`` runs per conformer while the key
   is derived from the conformer-INVARIANT OIN template), so a hit on a repeated key is
   the behaviour under test, not an implementation detail.
"""

import os
import unittest

from oinsmiles.core import chirality as ch
from oinsmiles.oin.levers import lever_enabled

LEVER = "OIN_MEMO_CIP_REPARSE"

#: (smiles, probe, fill_deficit) triples spanning every branch of the worker: a plain sp3
#: centre both ways round, an aromatic-adjacent centre (the rdCIPLabeler-flips-on-kekulisation
#: case the helper exists for), a fused-ring variant, a metal-stripped open-valence donor
#: (the ``fill_deficit`` path), and two inputs that legitimately yield ``None``.
CASES = [
    ("[C@H:99](F)(Cl)Br", 99, True),
    ("[C@@H:99](F)(Cl)Br", 99, True),
    ("C[C@H:99](N)C(=O)O", 99, True),
    ("c1ccccc1[C@H:99](C)CC", 99, True),
    ("c1ccc2c(c1)[C@H:99](C)CC2", 99, True),
    ("[C@H:99](C)(CC)[O]", 99, True),
    ("[C@:99]([O])(C)(CC)CCC", 99, False),
    ("[S@@:99](=O)([N])c1ccccc1", 99, True),
    ("CC(C)C", 99, True),
    ("[S@:99](=O)(=N)C", 99, False),
]

#: The one input in this file whose label DEPENDS on ``fill_deficit`` (None vs "R").
#: If a refactor drops ``fill_deficit`` from the memo key, this is what fails.
FILL_SENSITIVE = "O[C@:99]([O])(C)CC"


class _LeverCase(unittest.TestCase):
    """Save/restore the lever by VALUE, never by deleting it.

    Deleting would spell "off" as "take the default", which stops meaning "off" the moment
    the lever is promoted -- see ``test_levers::TestNoTestUnsetsAPromotedLever``.
    """

    def setUp(self):
        self._saved = os.environ.get(LEVER)
        ch._reparse_cip_memo_clear()

    def tearDown(self):
        if self._saved is None:
            os.environ.pop(LEVER, None)
        else:
            os.environ[LEVER] = self._saved
        ch._reparse_cip_memo_clear()


class TestMemoIsTransparent(_LeverCase):
    def test_warm_matches_cold_on_every_case(self):
        """The memo returns exactly what the worker would have computed."""
        os.environ[LEVER] = "1"
        for smiles, probe, fill in CASES:
            with self.subTest(smiles=smiles, fill=fill):
                cold = ch._reparse_cip_label_once_uncached(smiles, probe, fill)
                warm_miss = ch._reparse_cip_label_once(smiles, probe, fill)
                warm_hit = ch._reparse_cip_label_once(smiles, probe, fill)
                self.assertEqual(cold, warm_miss)
                self.assertEqual(cold, warm_hit)

    def test_lever_off_matches_lever_on(self):
        """Turning the lever on changes no answer -- the release's whole rule."""
        os.environ[LEVER] = "0"
        off = [ch._reparse_cip_label_once(s, p, f) for s, p, f in CASES]
        os.environ[LEVER] = "1"
        on = [ch._reparse_cip_label_once(s, p, f) for s, p, f in CASES]
        self.assertEqual(off, on)


class TestMemoKeyIsComplete(_LeverCase):
    def test_fill_deficit_is_part_of_the_key(self):
        """A real input whose answer depends on ``fill_deficit``, both orders.

        Run fill=True first and then fill=False, and again in the opposite order on a
        cleared cache: a key missing ``fill_deficit`` would return the first-computed
        answer for both, so one of the two orders would fail whichever way the bug went.
        """
        expected_fill = ch._reparse_cip_label_once_uncached(FILL_SENSITIVE, 99, True)
        expected_nofill = ch._reparse_cip_label_once_uncached(FILL_SENSITIVE, 99, False)
        self.assertNotEqual(
            expected_fill,
            expected_nofill,
            "fixture no longer discriminates fill_deficit -- pick another input, do not "
            "delete the test: it is the only guard on that key component",
        )

        os.environ[LEVER] = "1"
        for first, second in ((True, False), (False, True)):
            with self.subTest(order=(first, second)):
                ch._reparse_cip_memo_clear()
                got_first = ch._reparse_cip_label_once(FILL_SENSITIVE, 99, first)
                got_second = ch._reparse_cip_label_once(FILL_SENSITIVE, 99, second)
                self.assertEqual(got_first, expected_fill if first else expected_nofill)
                self.assertEqual(got_second, expected_fill if second else expected_nofill)

    def test_probe_is_part_of_the_key(self):
        """A probe that matches no atom yields None, even after a hit on a live probe."""
        os.environ[LEVER] = "1"
        smiles = "[C@H:99](F)(Cl)Br"
        self.assertIsNotNone(ch._reparse_cip_label_once(smiles, 99, True))
        self.assertIsNone(ch._reparse_cip_label_once(smiles, 7, True))

    def test_smiles_is_part_of_the_key(self):
        """Enantiomers differing only in the SMILES must not collide."""
        os.environ[LEVER] = "1"
        left = ch._reparse_cip_label_once("[C@H:99](F)(Cl)Br", 99, True)
        right = ch._reparse_cip_label_once("[C@@H:99](F)(Cl)Br", 99, True)
        self.assertEqual({left, right}, {"R", "S"})


class TestLeverSense(_LeverCase):
    def test_default_is_off_and_generates_no_traffic(self):
        os.environ.pop(LEVER, None)  # lever-lint: intentional-unset
        self.assertFalse(lever_enabled(LEVER), "OIN_MEMO_CIP_REPARSE must ship OFF")
        for smiles, probe, fill in CASES[:3]:
            ch._reparse_cip_label_once(smiles, probe, fill)
        info = ch._reparse_cip_memo_info()
        self.assertEqual((info.hits, info.misses, info.currsize), (0, 0, 0))

    def test_zero_disables_rather_than_enables(self):
        """``OIN_MEMO_CIP_REPARSE=0`` is OFF. ``os.environ.get`` would make it ON."""
        os.environ[LEVER] = "0"
        self.assertFalse(lever_enabled(LEVER))
        ch._reparse_cip_label_once(*CASES[0])
        self.assertEqual(ch._reparse_cip_memo_info().misses, 0)


class TestMemoActuallyHits(_LeverCase):
    def test_repeated_key_is_a_hit(self):
        """The traffic this lever targets: the same key asked for again.

        ``accept_fn`` runs per conformer, but both call sites key on a SMILES with no
        coordinates in it, so every conformer of a molecule asks the same question.
        """
        os.environ[LEVER] = "1"
        smiles, probe, fill = CASES[0]
        for _ in range(5):
            ch._reparse_cip_label_once(smiles, probe, fill)
        info = ch._reparse_cip_memo_info()
        self.assertEqual(info.misses, 1)
        self.assertEqual(info.hits, 4)

    def test_clear_drops_everything(self):
        os.environ[LEVER] = "1"
        ch._reparse_cip_label_once(*CASES[0])
        self.assertEqual(ch._reparse_cip_memo_info().currsize, 1)
        ch._reparse_cip_memo_clear()
        self.assertEqual(ch._reparse_cip_memo_info().currsize, 0)

    def test_cache_is_bounded(self):
        """A long single-interpreter sweep must not grow the memo without limit."""
        self.assertIsNotNone(ch._reparse_cip_memo_info().maxsize)
        self.assertEqual(ch._reparse_cip_memo_info().maxsize, ch._CIP_REPARSE_MEMO_MAX)


class TestPublicHelperUnchanged(_LeverCase):
    def test_aromatic_helper_agrees_across_the_lever(self):
        """``_reparse_aromatic_cip_label`` is the real caller -- pin it, not just the worker."""
        from rdkit import Chem

        mol = Chem.MolFromSmiles("c1ccc2c(c1)[C@H](C)CC2")
        self.assertIsNotNone(mol)
        idx = next(
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        )

        os.environ[LEVER] = "0"
        off = ch._reparse_aromatic_cip_label(mol, idx)
        os.environ[LEVER] = "1"
        on = ch._reparse_aromatic_cip_label(mol, idx)
        self.assertEqual(off, on)
        self.assertIsNotNone(off)


if __name__ == "__main__":
    unittest.main()
