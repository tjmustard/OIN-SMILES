"""``OIN_ENFORCE_BUDGET`` -- the requested budget is a BOUND, not a hint.

WHAT IS ACTUALLY BEING PINNED
=============================
``OIN3DGenerator(timeout=)`` is passed down as ``embed_time_budget`` and becomes a
deadline that, unset, is checked ONLY at the top of the embed attempt loop
(``generator3d/__init__.py``). An in-flight attempt therefore always runs to
completion, and two independent probes measured the consequence directly: 60 s asked,
60.7-137.9 s spent on an eta sample and 60.0-172.8 s on the boron set.

With the lever ON the same deadline is threaded into ``embed.get_embedding`` and
checked inside its two nested Python loops, and an empty pool at the deadline raises
``BudgetExhaustedError`` instead of returning ``[]``.

⚠ NOT PINNED HERE: the 759.9 s corpus figure the release was chartered on. That number
is a SUM over up to three separately SIGKILLed harness attempts, and all 4658
single-attempt rows in the 5k sweep finish within 0.2 s of their 300 s cap. See
``docs/agentic-notes/v0.4.9/ELAPSED_S_IS_A_SUM_v0.4.9.md``. These tests pin the
mechanism the code actually has, not the one the charter inferred.

WHY THE EMBED IS MONKEYPATCHED
==============================
A test that reached a real distance-geometry embed would take minutes and its timing
would depend on the box -- and a wall-clock assertion on a contended box is exactly the
kind of flake this project has already paid for. ``embed.get_embedding`` is replaced by
a stub that burns a controlled amount of time, so the bound is tested against a clock we
own. ``test_deadline_is_threaded_into_get_embedding`` separately pins that the real
function accepts and honours the parameter, which is the part a stub cannot prove.
"""

from __future__ import annotations

import os
import sys
import time
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "src"))

from oinsmiles.generation.metallogen_adapter import (  # noqa: E402
    convert_parsed_to_msmiles,
)
from oinsmiles.generation.oin_parser import OINParser  # noqa: E402
from oinsmiles.generator3d import (  # noqa: E402
    BudgetExhaustedError,
    StructuralAssemblyError,
    embed,
    generate_3d_structures,
)
from oinsmiles.oin.levers import lever_enabled  # noqa: E402

CISPLATIN_OIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"


def _msmiles(oin):
    return convert_parsed_to_msmiles(OINParser().parse(oin))


class TestBudgetIsABound(unittest.TestCase):
    def setUp(self):
        self.msmiles = _msmiles(CISPLATIN_OIN)
        self._real_embed = embed.get_embedding
        self.addCleanup(setattr, embed, "get_embedding", self._real_embed)

    def _slow_never_embeds(self, sleep_s):
        """An embed that burns ``sleep_s`` and never produces positions."""

        def stub(*a, **k):
            time.sleep(sleep_s)
            return None

        return stub

    def test_empty_pool_at_the_deadline_raises_typed(self):
        """Budget exhaustion is its own error, not a generic empty pool.

        Without this distinction v0.4.10 cannot tell its own regressions from v0.4.9's
        intended behaviour: both would arrive as ``MetalloGen failed to generate any
        conformers``, which ``tools/classify_failures.py`` buckets as ``no_conformers``.
        """
        embed.get_embedding = self._slow_never_embeds(0.05)
        with self.assertRaises(BudgetExhaustedError) as ctx:
            generate_3d_structures(
                self.msmiles, num_conformers=1, embed_time_budget=0.2, enforce_budget=True
            )
        msg = str(ctx.exception)
        # The budget AND the spend, because a bound is only interpretable next to the
        # number it was given -- ULODUU assembles at 60 s and not at 30 s.
        self.assertIn("0.2s requested", msg)
        self.assertIn("spent", msg)

    def test_lever_off_returns_empty_not_raises(self):
        """Default path is byte-identical: the same run returns ``[]`` as before."""
        embed.get_embedding = self._slow_never_embeds(0.05)
        out = generate_3d_structures(
            self.msmiles, num_conformers=1, embed_time_budget=0.2, enforce_budget=False
        )
        self.assertEqual(out, [])

    def test_no_budget_means_no_bound(self):
        """``embed_time_budget=None`` is unbounded even with the lever on.

        Direct callers that never asked for a budget must not acquire one. The attempt
        loop still terminates via ``max_attempts``.
        """
        embed.get_embedding = self._slow_never_embeds(0.0)
        out = generate_3d_structures(
            self.msmiles, num_conformers=1, embed_time_budget=None, enforce_budget=True
        )
        self.assertEqual(out, [])

    def test_structural_failure_outranks_budget_exhaustion(self):
        """A uniformly structural failure is the better diagnosis.

        It would have happened at any budget, so ``StructuralAssemblyError`` must win --
        otherwise a real, permanent assembly defect gets relabelled "out of time" and
        the next release wastes a lane chasing a compute problem that does not exist.
        """

        def raises_structural(*a, **k):
            time.sleep(0.05)
            raise IndexError("no binding vector at slot 3")

        embed.get_embedding = raises_structural
        with self.assertRaises(StructuralAssemblyError):
            generate_3d_structures(
                self.msmiles, num_conformers=1, embed_time_budget=0.2, enforce_budget=True
            )

    def test_deadline_is_threaded_into_get_embedding(self):
        """The bound reaches the cost sink, and only when enforcing.

        Profiled on ``FOSNEI_comp_0`` at a 300 s budget, ``get_embedding`` is **61.5 s of
        self time** out of an 82.4 s generation, against 1.74 s for the CBC solve (2.1%)
        and 0.63 s for the ``accept_fn`` re-encode (0.8%). The charter named the latter
        two as the likely sinks; bounding either would have measured as "no change".
        A deadline that does not arrive here is not a bound.
        """
        seen = []

        def capture(*a, **k):
            seen.append(k.get("deadline"))
            return None

        embed.get_embedding = capture
        generate_3d_structures(
            self.msmiles, num_conformers=1, embed_time_budget=5.0, enforce_budget=True
        )
        self.assertTrue(seen, "get_embedding was never called")
        self.assertTrue(
            all(d is not None for d in seen),
            "enforcing, but get_embedding received deadline=None -- the bound never reaches "
            "the loop where the time is actually spent",
        )

        seen.clear()
        generate_3d_structures(
            self.msmiles, num_conformers=1, embed_time_budget=5.0, enforce_budget=False
        )
        self.assertTrue(seen)
        self.assertTrue(
            all(d is None for d in seen),
            "lever OFF must leave get_embedding unbounded and byte-identical to pristine",
        )

    def test_real_get_embedding_accepts_and_honours_deadline(self):
        """Against the REAL function, not the stub.

        The stub tests above would all pass against a ``get_embedding`` that accepted
        ``deadline`` and ignored it. An already-expired deadline must make the real
        function return ``None`` without doing the work.
        """
        parsed = OINParser().parse(CISPLATIN_OIN)
        from oinsmiles.generation.metallogen_adapter import (  # noqa: PLC0415
            _prepare_ligand_fragments,
        )
        from oinsmiles.generator3d import om

        metal_frag, ligand_specs, geo = _prepare_ligand_fragments(parsed)
        complex_ = om.get_om_from_parsed(metal_frag, ligand_specs, geo)

        t0 = time.monotonic()
        out = self._real_embed(complex_, 1.0, 0, align=True, seed=42, deadline=time.monotonic() - 1)
        self.assertIsNone(out, "an expired deadline must abandon the embed")
        self.assertLess(
            time.monotonic() - t0,
            2.0,
            "returned None but still spent the time -- the check is in the wrong place",
        )


class TestLeverRegistration(unittest.TestCase):
    def test_default_off(self):
        """v0.4.9 ships this OFF. Promotion is an accuracy decision, not a runtime one.

        93.1% of honest passes already finish under 30 s, so a 30 s bound recovers
        37.8 CPU-h per 5000-molecule sweep but costs 251 passes = 5.02 points of
        ``byte_exact``, against a headline goal of 100%.
        """
        os.environ.pop("OIN_ENFORCE_BUDGET", None)
        self.assertFalse(lever_enabled("OIN_ENFORCE_BUDGET"))

    def test_zero_disables(self):
        """``"0"`` must DISABLE. ``os.environ.get`` would return a truthy non-empty string.

        That trap cost 23 test failures across two promotions, which is why
        ``lever_enabled`` exists and why reading the env directly is linted against.
        """
        os.environ["OIN_ENFORCE_BUDGET"] = "0"
        self.addCleanup(os.environ.pop, "OIN_ENFORCE_BUDGET", None)
        self.assertFalse(lever_enabled("OIN_ENFORCE_BUDGET"))

    def test_one_enables(self):
        os.environ["OIN_ENFORCE_BUDGET"] = "1"
        self.addCleanup(os.environ.pop, "OIN_ENFORCE_BUDGET", None)
        self.assertTrue(lever_enabled("OIN_ENFORCE_BUDGET"))

    def test_has_a_rationale(self):
        """Every lever carries why it exists and what would promote it."""
        from oinsmiles.oin import levers

        rationale = levers._HELD_OFF["OIN_ENFORCE_BUDGET"]
        self.assertIn("BOUND", rationale)
        self.assertIn("ELAPSED_S_IS_A_SUM", rationale)


if __name__ == "__main__":
    unittest.main()
