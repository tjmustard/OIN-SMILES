"""`OIN_ACCEPT_SCORED` / `independent_confirm` — the acceptance predicate's two modes.

The lever exists because `accept_fn` was measured to be STRICTER than the predicate that scores
a round trip: it adds an independent `XYZToSMILES().convert` re-perception the score never asks
for. On HIDCIH_comp_1 the cheap test matched at pool conformer 0 (1.66 s) while the strict one
first matched at conformer 25 (49.4 s) — 44 conformers were scored-successes that acceptance
discarded. See `docs/agentic-notes/v0.4.5-retrospective/LANE-eta-runtime-30s.md`.

These tests pin the three behaviours that matter and are cheap to get wrong:

1. default (`independent_confirm=True`) still demands the strict confirm — so the shipped path is
   unchanged and a cheap-only match is NOT accepted;
2. with the lever on, a cheap match accepts without the strict call ever being made;
3. `fast is None` (perception failed) does NOT accept even with the lever on — an unperceivable
   conformer is not a scored success, and must fall through to the strict path.

Mocked at the seams rather than generating a molecule: the question here is predicate logic, and
a real pool fill costs minutes.
"""

import unittest
from unittest import mock

from oinsmiles.generation import metallogen_adapter as MA

TARGET = "TARGET-KEY"


def _key(s):
    """Stand-in for canonical_roundtrip_key: the string IS the key."""
    return s


class TestIndependentConfirm(unittest.TestCase):
    def _run(self, *, fast, full, confirm):
        """Evaluate the predicate with the two re-encode paths stubbed.

        Returns ``(verdict, strict_was_called)``.
        """
        calls = {"strict": 0}

        def fake_full(_m):
            calls["strict"] += 1
            return full

        with (
            mock.patch.object(MA, "build_contract_mol", return_value=object()),
            mock.patch.object(MA, "_reencode_oin_fast", return_value=fast),
            mock.patch.object(MA, "_reencode_oin", side_effect=fake_full),
            mock.patch.object(MA, "canonical_roundtrip_key", side_effect=_key),
        ):
            verdict = MA._reencode_key_matches(
                object(), object(), TARGET, independent_confirm=confirm
            )
        return verdict, calls["strict"] > 0

    def test_default_still_requires_the_strict_confirm(self):
        """A cheap match with a strict MISMATCH is rejected on the shipped default path."""
        verdict, strict_called = self._run(fast=TARGET, full="SOMETHING-ELSE", confirm=True)
        self.assertFalse(verdict, "default acceptance must not trust the cheap test alone")
        self.assertTrue(strict_called, "the default path must actually run the strict re-encode")

    def test_lever_accepts_a_cheap_match_without_the_strict_call(self):
        """This is the entire runtime win: skip the expensive independent re-perception."""
        verdict, strict_called = self._run(fast=TARGET, full="SOMETHING-ELSE", confirm=False)
        self.assertTrue(verdict, "a scored-success conformer must be accepted")
        self.assertFalse(strict_called, "the strict re-encode is what the lever exists to skip")

    def test_cheap_mismatch_is_rejected_in_both_modes(self):
        """The prefilter's REJECT power is unchanged -- it is not a shortcut to acceptance."""
        for confirm in (True, False):
            with self.subTest(independent_confirm=confirm):
                verdict, _ = self._run(fast="WRONG", full=TARGET, confirm=confirm)
                self.assertFalse(verdict)

    def test_unperceivable_conformer_falls_through_rather_than_accepting(self):
        """``fast is None`` means perception FAILED, which is not a scored success.

        With the lever on it must still consult the strict path, and accept only if that agrees.
        Accepting blind here would hand back a conformer nothing has verified.
        """
        verdict, strict_called = self._run(fast=None, full=TARGET, confirm=False)
        self.assertTrue(strict_called, "a perception failure must not short-circuit to accept")
        self.assertTrue(verdict, "...but a strict match on that conformer is still acceptable")

        verdict, strict_called = self._run(fast=None, full="SOMETHING-ELSE", confirm=False)
        self.assertFalse(verdict)
        self.assertTrue(strict_called)


class TestLeverIsOffByDefault(unittest.TestCase):
    """Guard the shipped configuration, matching the Lane 4 `OIN_EMIT_AXIAL` pattern.

    Promoting this on the 22-molecule cohort would be the mistake the lane exists to correct;
    the recorded gate is a corpus A/B.
    """

    def test_default_off(self):
        from oinsmiles.oin.levers import held_off, lever_enabled

        self.assertIn("OIN_ACCEPT_SCORED", held_off())
        import os

        env = {k: v for k, v in os.environ.items() if k != "OIN_ACCEPT_SCORED"}
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertFalse(lever_enabled("OIN_ACCEPT_SCORED"))


if __name__ == "__main__":
    unittest.main()
