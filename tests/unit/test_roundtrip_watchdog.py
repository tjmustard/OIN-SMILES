"""Guard tests for the round-trip harness watchdog (tools/test_dataset_roundtrip.py).

``--mol-timeout`` used to be a ``signal.alarm``, which cannot interrupt a hang inside
native code or a tight loop that never returns to the interpreter. UGUHAH_comp_0 wedged
a Phase-0 shard for 35+ minutes despite a 420 s cap -- and it hangs inside
``XYZToSMILES.convert``, i.e. the *encode*, not the generator. ``_supervise`` replaces
the alarm with a real subprocess kill.

These tests drive ``_supervise`` with fake process/queue doubles. Spawning real children
from under ``python -m unittest`` would make each child re-run test discovery.
"""

import os
import queue as queue_mod
import sys
import time
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../tools")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../integration")))

from test_dataset_roundtrip import _ENCODED, _supervise

TIMEOUT = 0.4


class FakeQueue:
    """Yields queued messages, then behaves like a queue that never fills."""

    def __init__(self, messages=()):
        self._messages = list(messages)

    def get(self, timeout=None):
        if self._messages:
            return self._messages.pop(0)
        time.sleep(0.02)  # keep the supervisor's poll loop off a hot spin
        raise queue_mod.Empty


class FakeProc:
    """A stand-in child process.

    Args:
        alive_for: Number of ``is_alive`` polls before it reports dead. None = never
            exits on its own (a wedged child).
        dies_on_join: Whether ``join`` reaps it. True models a child that exits once it
            has delivered its payload; False models one that lingers afterwards.
    """

    def __init__(self, alive_for=None, exitcode=0, dies_on_join=True):
        self.alive_for = alive_for
        self.exitcode = exitcode
        self.dies_on_join = dies_on_join
        self.started = False
        self.killed = False
        self._dead = False
        self._polls = 0

    def start(self):
        self.started = True

    def is_alive(self):
        if self.killed or self._dead:
            return False
        self._polls += 1
        return self.alive_for is None or self._polls <= self.alive_for

    def kill(self):
        self.killed = True

    def join(self, timeout=None):
        if self.dies_on_join:
            self._dead = True


class TestSupervisorHappyPath(unittest.TestCase):
    def test_payload_is_returned_and_child_not_killed(self):
        payload = (True, {"status": "success"}, "xyz")
        proc, report = FakeProc(), {}
        result = _supervise(proc, FakeQueue([payload]), TIMEOUT, "UFF_1", report)
        self.assertEqual(result, payload)
        self.assertTrue(proc.started)
        self.assertFalse(proc.killed, "a child that exits on its own must not be killed")
        self.assertEqual(report, {})

    def test_child_lingering_after_its_payload_is_still_reaped(self):
        payload = (True, {"status": "success"}, "xyz")
        proc = FakeProc(dies_on_join=False)
        result = _supervise(proc, FakeQueue([payload]), TIMEOUT, "UFF_1", {})
        self.assertEqual(result, payload)
        self.assertTrue(proc.killed, "a child still alive after join must be killed")

    def test_encoded_marker_is_not_mistaken_for_the_result(self):
        payload = (True, {"status": "success"}, "xyz", "[Pd_SPL].[Cl]{0}")
        progress = {}
        result = _supervise(
            FakeProc(),
            FakeQueue([(_ENCODED, "[Pd_SPL].[Cl]{0}"), payload]),
            TIMEOUT,
            "UFF_1",
            {},
            progress,
            stage="encoding",
        )
        self.assertEqual(result, payload)
        self.assertEqual(progress["oin1"], "[Pd_SPL].[Cl]{0}")


class TestSupervisorTimeout(unittest.TestCase):
    def test_hang_during_encode_is_killed_and_named(self):
        proc, report, progress = FakeProc(), {}, {}
        result = _supervise(proc, FakeQueue(), TIMEOUT, "UFF_1", report, progress, stage="encoding")
        self.assertIsNone(result)
        self.assertTrue(proc.killed, "a wedged child must be SIGKILLed, not waited on")
        self.assertEqual(report["status"], "failed")
        self.assertIn("TimeoutException at UFF_1", report["error"])
        self.assertIn("while encoding", report["error"])
        self.assertNotIn("oin1", progress)

    def test_hang_during_generate_keeps_the_encode(self):
        # The child got as far as reporting smiles_1 before wedging; losing that would
        # blank out the report of every timed-out molecule.
        proc, report, progress = FakeProc(), {}, {}
        result = _supervise(
            proc,
            FakeQueue([(_ENCODED, "[Ni_SPL].[Cl]{0}")]),
            TIMEOUT,
            "UFF_1",
            report,
            progress,
            stage="encoding",
        )
        self.assertIsNone(result)
        self.assertTrue(proc.killed)
        self.assertIn("while generating", report["error"])
        self.assertEqual(progress["oin1"], "[Ni_SPL].[Cl]{0}")

    def test_pass_two_timeout_never_claims_to_be_encoding(self):
        # The pass-2 child is handed an OIN string; it has no encode stage at all.
        report = {}
        _supervise(FakeProc(), FakeQueue(), TIMEOUT, "g-xTB_5", report)
        self.assertIn("while generating", report["error"])
        self.assertNotIn("encoding", report["error"])

    def test_timeout_error_string_stays_in_the_classifier_timeout_bin(self):
        from classify_failures import classify

        report = {}
        _supervise(FakeProc(), FakeQueue(), TIMEOUT, "UFF_1", report, stage="encoding")
        cls, _evidence = classify(
            {"status": "failed", "error": report["error"], "smiles_1": "x", "smiles_2": "y"}
        )
        self.assertEqual(cls, "timeout")


class TestSupervisorCrash(unittest.TestCase):
    def test_dead_child_is_reported_as_a_hard_failure_without_waiting(self):
        # A segfaulted child must not burn the whole timeout budget before being noticed.
        proc, report = FakeProc(alive_for=0, exitcode=-11), {}
        started = time.monotonic()
        result = _supervise(proc, FakeQueue(), TIMEOUT, "UFF_1", report)
        self.assertIsNone(result)
        self.assertLess(time.monotonic() - started, TIMEOUT)
        self.assertFalse(proc.killed, "an already-dead child must not be killed again")
        self.assertEqual(report["status"], "failed")
        self.assertIn("child process died with exit code -11", report["error"])

    def test_crash_error_string_is_a_hard_failure_for_the_pass_one_gate(self):
        # PASS 1 skips the g-xTB tiers only for errors with these prefixes.
        report = {}
        _supervise(FakeProc(alive_for=0, exitcode=1), FakeQueue(), TIMEOUT, "UFF_1", report)
        self.assertTrue(
            report["error"].startswith(("Generation/Verification failed at", "TimeoutException"))
        )


if __name__ == "__main__":
    unittest.main()
