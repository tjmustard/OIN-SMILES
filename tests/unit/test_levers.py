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


class TestOverride(unittest.TestCase):
    def test_explicit_override_beats_environment(self):
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_BODY": "1"}):
            self.assertFalse(lever_enabled("OIN_CANONICAL_BODY", override=False))
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_BODY": "0"}):
            self.assertTrue(lever_enabled("OIN_CANONICAL_BODY", override=True))

    def test_none_override_falls_through_to_environment(self):
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_BODY": "0"}):
            self.assertFalse(lever_enabled("OIN_CANONICAL_BODY", override=None))


if __name__ == "__main__":
    unittest.main()
