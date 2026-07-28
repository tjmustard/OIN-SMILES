"""Guards for the lever registry, including the sense-inversion trap it was built to fix."""

import os
import unittest
from unittest import mock

from oinsmiles.oin.levers import default_on, held_off, lever_enabled


class TestLeverSense(unittest.TestCase):
    def test_zero_disables(self):
        """The trap this module exists to close.

        The older bare read was ``os.environ.get("OIN_EMIT_AXIAL")``, which is truthy for the
        string ``"0"`` -- so ``OIN_EMIT_AXIAL=0`` ENABLED the lever. Anyone opting out the
        obvious way got the opposite of what they asked for.
        """
        for val in ("0", "false", "FALSE", "no", "off", "", "  ", " Off "):
            with self.subTest(value=val), mock.patch.dict(os.environ, {"OIN_CANONICAL_BODY": val}):
                self.assertFalse(lever_enabled("OIN_CANONICAL_BODY"), f"{val!r} should disable")

    def test_other_values_enable(self):
        for val in ("1", "true", "yes", "on", "anything"):
            with self.subTest(value=val), mock.patch.dict(os.environ, {"OIN_EMIT_AXIAL": val}):
                self.assertTrue(lever_enabled("OIN_EMIT_AXIAL"), f"{val!r} should enable")


class TestDefaults(unittest.TestCase):
    def test_promoted_levers_default_on(self):
        for name in default_on():
            with self.subTest(lever=name):
                env = {k: v for k, v in os.environ.items() if k != name}
                with mock.patch.dict(os.environ, env, clear=True):
                    self.assertTrue(lever_enabled(name), f"{name} should be on by default")

    def test_held_levers_default_off(self):
        for name in held_off():
            with self.subTest(lever=name):
                env = {k: v for k, v in os.environ.items() if k != name}
                with mock.patch.dict(os.environ, env, clear=True):
                    self.assertFalse(lever_enabled(name), f"{name} must stay opt-in")

    def test_promoted_and_held_are_disjoint(self):
        self.assertFalse(default_on() & set(held_off()))

    def test_every_held_lever_states_a_reason(self):
        """A lever held back without a recorded reason becomes folklore."""
        for name, reason in held_off().items():
            with self.subTest(lever=name):
                self.assertGreater(len(reason), 40, f"{name}'s reason is too thin to be useful")

    def test_unknown_lever_defaults_off(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertFalse(lever_enabled("OIN_NOT_A_REAL_LEVER"))


class TestDonorFoldAndParityVetoAreCoupled(unittest.TestCase):
    """The donor fold and its reflection-parity veto ship together or not at all (v0.4.13).

    ``OIN_CANONICAL_DONOR_FOLD`` ON with ``OIN_FOLD_PARITY_VETO`` OFF is the exact configuration
    v0.4.11 measured and refuted: it buys +7.86 ``byte_exact`` points and COLLAPSES ENANTIOMERS in
    221 of those same 393 gains (56.2%), plus 19 of a uniform 250-molecule draw.

    What makes that worth a dedicated guard rather than a comment: **neither the headline metric
    nor the round-trip comparison key can see the damage.** ``compare._parse_vertex_colors``
    colours every donor of a ligand with that ligand's whole body, so the axis this fold can
    destroy is deliberately folded by the key as well. A future session that demotes the veto
    alone -- to "isolate a regression", or because the veto costs an encode -- would see a green
    suite, green gate arms, and a HIGHER pass rate, with the loss invisible in every instrument
    the project routinely reads. Only ``tools/mirror_audit_donor_fold.py`` can detect it.

    So the invariant is pinned here, at the registry, where it cannot be missed.
    """

    FOLD = "OIN_CANONICAL_DONOR_FOLD"
    VETO = "OIN_FOLD_PARITY_VETO"

    def test_both_are_promoted_together(self):
        on = default_on()
        self.assertEqual(
            self.FOLD in on,
            self.VETO in on,
            f"{self.FOLD} and {self.VETO} must share a default. The fold without the veto "
            "collapses enantiomers and no gate arm, golden, or comparison key can see it -- "
            "mirror_audit_donor_fold.py is the only instrument that can. Demote both or neither.",
        )

    def test_the_shipped_default_is_both_on(self):
        """Pins the v0.4.13 promotion itself, so a silent revert is a test failure."""
        env = {k: v for k, v in os.environ.items() if k not in (self.FOLD, self.VETO)}
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertTrue(lever_enabled(self.FOLD), f"{self.FOLD} ships ON since v0.4.13")
            self.assertTrue(lever_enabled(self.VETO), f"{self.VETO} ships ON since v0.4.13")

    def test_neither_is_listed_as_held_off(self):
        held = held_off()
        self.assertNotIn(self.FOLD, held)
        self.assertNotIn(self.VETO, held)


class TestOverride(unittest.TestCase):
    def test_explicit_override_beats_environment(self):
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_BODY": "1"}):
            self.assertFalse(lever_enabled("OIN_CANONICAL_BODY", override=False))
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_BODY": "0"}):
            self.assertTrue(lever_enabled("OIN_CANONICAL_BODY", override=True))

    def test_none_override_falls_through_to_environment(self):
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_BODY": "0"}):
            self.assertFalse(lever_enabled("OIN_CANONICAL_BODY", override=None))


class TestNoTestUnsetsAPromotedLever(unittest.TestCase):
    """A lint, not a unit test: no test may spell "lever off" by DELETING the variable.

    Deleting a lever's env var means "take the default". While every lever defaulted OFF that was
    indistinguishable from disabling it. The moment a lever is promoted, every
    ``env.pop(LEVER)`` or ``{k: v for k, v in os.environ.items() if k != LEVER}`` silently becomes
    a second lever-ON test -- and a test asserting the OLD behaviour then fails for a reason
    unrelated to what it tests, while a test asserting *stability* goes quietly vacuous instead.

    This trap has cost **23 test failures across two promotions** (17 in v0.4.5, 6 more when
    ``OIN_BORON_CAGE`` was promoted in v0.4.6) and was diagnosed from scratch each time. The fix
    is always the same one line -- write ``"0"`` instead of deleting -- so it becomes mechanical
    here. Prose in a docstring did not prevent occurrences three and four; a failing test will.
    """

    #: `pop` on a lever is legitimate when RESTORING a saved pre-test value. The guard is usually
    #: on a PRECEDING line (``if self._saved is None:`` / ``if cls._prior is None:``), so the scan
    #: looks back a couple of lines -- checking only the pop's own line produced two false
    #: positives against correct teardown code on first run.
    _RESTORE_HINTS = ("_saved", "_prev", "_prior", "restore", "tearDown")
    _LOOKBACK = 3

    #: Explicit opt-out for the rare test whose SUBJECT is the unset state itself -- e.g. "an
    #: unset variable must take the promoted default". Put this marker on the pop line. An
    #: auditable marker is deliberately preferred over a looser regex: widening the heuristic
    #: to accommodate one legitimate case would blind the lint to the illegitimate ones.
    _ALLOW = "lever-lint: intentional-unset"

    def test_no_test_file_unsets_a_default_on_lever(self):
        import pathlib
        import re

        promoted = sorted(default_on())
        here = pathlib.Path(__file__).resolve()
        pop_re = re.compile(r"(?:environ|env)\.pop\(\s*[\"']?([A-Za-z_]+)")
        strip_re = re.compile(r"if\s+k\s*!=\s*[\"']([A-Z_]+)[\"']")

        offenders = []
        for path in sorted(here.parent.glob("test_*.py")):
            if path.name == here.name:
                continue
            text = path.read_text()
            # resolve a module-level `LEVER = "OIN_..."` alias so `env.pop(LEVER)` is caught too
            aliases = dict(re.findall(r"^([A-Z_]+)\s*=\s*[\"'](OIN_[A-Z_]+)[\"']", text, re.M))
            lines = text.splitlines()
            for lineno, line in enumerate(lines, 1):
                if self._ALLOW in line:
                    continue
                window = lines[max(0, lineno - 1 - self._LOOKBACK) : lineno]
                if any(h in w for w in window for h in self._RESTORE_HINTS):
                    continue
                for m in list(pop_re.finditer(line)) + list(strip_re.finditer(line)):
                    name = aliases.get(m.group(1), m.group(1))
                    if name in promoted:
                        offenders.append(f'{path.name}:{lineno} unsets {name} -- write "0"')

        self.assertEqual(
            offenders,
            [],
            "these express a promoted lever's OFF state by DELETING the variable, which now "
            "means ON:\n  " + "\n  ".join(offenders),
        )


if __name__ == "__main__":
    unittest.main()
