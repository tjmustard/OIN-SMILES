"""Guard tests for the per-stage benchmark harness (tools/benchmark_generation.py).

The harness measures generator wall-clock, so its own logic -- stage accounting, the
median/IQR/UNSTABLE math, CN-stratified sampling, the serial guard, and the streamed-run
watchdog -- must be correct independently of any real (slow, env-heavy) generation. These
tests drive the pure helpers directly and the watchdog with injected fakes, exactly as
``test_roundtrip_watchdog.py`` does for the round-trip supervisor. No real generation runs.
"""

import os
import sys
import tempfile
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../tools")))

import benchmark_generation as B  # noqa: E402


class TestGeoToCN(unittest.TestCase):
    def test_known_tokens(self):
        self.assertEqual(B.geo_to_cn("Pd_SPL"), 4)
        self.assertEqual(B.geo_to_cn("Ir_OCT"), 6)
        self.assertEqual(B.geo_to_cn("Au_LIN"), 2)
        self.assertEqual(B.geo_to_cn("Fe_TBP"), 5)
        self.assertEqual(B.geo_to_cn("X_PBP"), 7)
        self.assertEqual(B.geo_to_cn("Y_SQA"), 8)
        self.assertEqual(B.geo_to_cn("Mn_TPL"), 3)

    def test_bare_code_and_unknowns(self):
        self.assertEqual(B.geo_to_cn("SPL"), 4)  # no metal prefix
        self.assertIsNone(B.geo_to_cn("Foo_ZZZ"))
        self.assertIsNone(B.geo_to_cn(""))
        self.assertIsNone(B.geo_to_cn(None))


class TestSummarize(unittest.TestCase):
    def test_tight_series_is_stable(self):
        s = B.summarize([5.05, 5.47, 5.48, 5.50, 5.51])
        self.assertAlmostEqual(s["median"], 5.48, places=2)
        self.assertFalse(s["unstable"])
        self.assertLess(s["iqr_pct"], 20.0)
        self.assertEqual(s["n"], 5)

    def test_wide_series_is_unstable(self):
        s = B.summarize([1.0, 5.0, 5.0, 9.0, 10.0])
        self.assertTrue(s["unstable"])
        self.assertGreater(s["iqr_pct"], 20.0)

    def test_single_value_has_zero_iqr(self):
        s = B.summarize([3.3])
        self.assertEqual(s["iqr"], 0.0)
        self.assertFalse(s["unstable"])

    def test_submilli_stage_not_flagged_unstable(self):
        # A jittery but truly negligible stage (sub-ms median) is reported count-first,
        # not flagged UNSTABLE.
        s = B.summarize([0.0001, 0.0002, 0.0009, 0.0003, 0.0005])
        self.assertTrue(s["subms"])
        self.assertFalse(s["unstable"])

    def test_empty(self):
        s = B.summarize([])
        self.assertEqual(s["median"], 0.0)
        self.assertFalse(s["unstable"])


class TestStageRecorderAndTiming(unittest.TestCase):
    def test_add_accumulates(self):
        rec = B.StageRecorder()
        rec.add("s", 1.0, fail=0)
        rec.add("s", 2.0, fail=1)
        self.assertEqual(rec.stages["s"]["wall"], 3.0)
        self.assertEqual(rec.stages["s"]["count"], 2)
        self.assertEqual(rec.stages["s"]["fail"], 1)
        rec.reset()
        self.assertEqual(rec.stages, {})
        self.assertEqual(rec.depth, 0)

    def test_recursion_timed_and_counted_once(self):
        # Model get_valid_molecule's self-recursion: the wrapper must attribute exactly
        # one call even though the target calls itself.
        rec = B.StageRecorder()
        calls = {"n": 0}

        def recursive(self_obj, depth):
            calls["n"] += 1
            if depth > 0:
                # nested invocation goes back through the same reentrancy-guarded timer
                B._timed_outermost(
                    rec, "perception.PuLP", recursive, (self_obj, depth - 1), {}, lambda r: False
                )
            return None

        B._timed_outermost(rec, "perception.PuLP", recursive, (object(), 3), {}, lambda r: False)
        self.assertEqual(calls["n"], 4)  # target really recursed 4 deep
        self.assertEqual(rec.stages["perception.PuLP"]["count"], 1)  # ...but timed once

    def test_fail_predicate_marks_failure(self):
        rec = B.StageRecorder()
        B._timed_outermost(rec, "perception.PuLP", lambda s: None, (None,), {}, lambda r: r is None)
        self.assertEqual(rec.stages["perception.PuLP"]["fail"], 1)

    def test_embed_failed_semantics(self):
        self.assertEqual(B.embed_failed(-1), 1)  # RDKit failure sentinel
        self.assertEqual(B.embed_failed(0), 0)  # conformer id 0 is success
        self.assertEqual(B.embed_failed(3), 0)


class TestAggregateRuns(unittest.TestCase):
    def _rec(self, total, pulp_wall, sha):
        return {
            "total": total,
            "atoms": 11,
            "xyz_sha": sha,
            "stages": {"perception.PuLP": {"wall": pulp_wall, "count": 198, "fail": 158}},
        }

    def test_folds_and_detects_byte_identity(self):
        runs = [self._rec(5.0, 4.5, "abc"), self._rec(5.2, 4.6, "abc"), self._rec(5.1, 4.55, "abc")]
        agg = B.aggregate_runs(runs)
        self.assertEqual(agg["runs"], 3)
        self.assertTrue(agg["byte_identical"])
        self.assertAlmostEqual(agg["total"]["median"], 5.1, places=2)
        self.assertEqual(agg["stages"]["perception.PuLP"]["count"], 198)
        self.assertEqual(agg["stages"]["perception.PuLP"]["fail"], 158)

    def test_detects_xyz_variation(self):
        runs = [self._rec(5.0, 4.5, "abc"), self._rec(5.2, 4.6, "DIFFERENT")]
        self.assertFalse(B.aggregate_runs(runs)["byte_identical"])

    def test_empty(self):
        self.assertEqual(B.aggregate_runs([])["runs"], 0)


class TestFindConflicts(unittest.TestCase):
    def test_excludes_self_pids(self):
        out = {"benchmark_generation": "999\n1000\n", "test_dataset_roundtrip": ""}
        self.assertEqual(B.find_conflicts(out, {999, 1000}), [])

    def test_detects_competing_sweep(self):
        out = {"test_dataset_roundtrip": "4242 4243\n", "benchmark_generation": "999\n"}
        conflicts = B.find_conflicts(out, {999})
        self.assertEqual(len(conflicts), 1)
        self.assertIn("test_dataset_roundtrip", conflicts[0])
        self.assertIn("4242", conflicts[0])

    def test_ancestor_pids_includes_self_and_parent(self):
        # The guard excludes the whole ancestry (python + its uv/bash launchers), so a
        # `bash -c "uv run python ...benchmark_generation..."` launch is not misread as a
        # competing run. At minimum this process and its parent are present.
        anc = B._ancestor_pids()
        self.assertIn(os.getpid(), anc)
        self.assertIn(os.getppid(), anc)


class TestStratifiedSample(unittest.TestCase):
    def _registry(self):
        reg = []
        # 10 CN4 (SPL), 6 CN6 (OCT), 3 CN2 (LIN), 1 CN5 (TBP)
        for i in range(10):
            reg.append({"molecule": f"spl{i}", "metal_geo": "Pd_SPL"})
        for i in range(6):
            reg.append({"molecule": f"oct{i}", "metal_geo": "Ir_OCT"})
        for i in range(3):
            reg.append({"molecule": f"lin{i}", "metal_geo": "Au_LIN"})
        reg.append({"molecule": "tbp0", "metal_geo": "Fe_TBP"})
        return reg

    def _index(self, reg):
        return {r["molecule"]: f"/data/{r['molecule']}.xyz" for r in reg}

    def test_deterministic_and_in_range(self):
        reg = self._registry()
        idx = self._index(reg)
        sel1, observed = B.stratified_sample(reg, idx, 8, 2, 7, seed=42)
        sel2, _ = B.stratified_sample(reg, idx, 8, 2, 7, seed=42)
        self.assertEqual([d["molecule"] for d in sel1], [d["molecule"] for d in sel2])
        self.assertEqual(len(sel1), 8)
        self.assertEqual(observed, {2: 3, 4: 10, 5: 1, 6: 6})
        for d in sel1:
            self.assertIn(d["cn"], (2, 4, 5, 6))

    def test_small_buckets_represented(self):
        reg = self._registry()
        idx = self._index(reg)
        sel, _ = B.stratified_sample(reg, idx, 8, 2, 7, seed=1)
        cns = {d["cn"] for d in sel}
        # every non-empty in-range CN gets at least one slot
        self.assertEqual(cns, {2, 4, 5, 6})

    def test_only_eligible_molecules(self):
        reg = self._registry()
        idx = self._index(reg)
        del idx["spl0"]  # remove one from the path index -> ineligible
        sel, _ = B.stratified_sample(reg, idx, 20, 2, 7, seed=1)
        self.assertNotIn("spl0", [d["molecule"] for d in sel])

    def test_cn_range_filter(self):
        reg = self._registry()
        idx = self._index(reg)
        sel, observed = B.stratified_sample(reg, idx, 20, 4, 4, seed=1)  # CN4 only
        self.assertTrue(all(d["cn"] == 4 for d in sel))


class TestBuildPathIndex(unittest.TestCase):
    def test_scans_cat_and_photo_excluding_generated(self):
        with tempfile.TemporaryDirectory() as d:
            for sub, names in (("cat", ["A_comp_0", "B_comp_0"]), ("photo", ["C_comp_0"])):
                os.makedirs(os.path.join(d, sub))
                for n in names:
                    open(os.path.join(d, sub, f"{n}.xyz"), "w").close()
            # a generated output must be excluded
            open(os.path.join(d, "cat", "A_comp_0_generated.xyz"), "w").close()
            idx = B.build_path_index(d)
            self.assertEqual(set(idx), {"A_comp_0", "B_comp_0", "C_comp_0"})
            self.assertTrue(idx["A_comp_0"].endswith("cat/A_comp_0.xyz"))


class FakeClock:
    """Monotonic clock that advances a fixed step on each read."""

    def __init__(self, step=1.0):
        self.t = 0.0
        self.step = step

    def __call__(self):
        self.t += self.step
        return self.t


class TestCollectRuns(unittest.TestCase):
    def _msgs(self, seq):
        seq = list(seq)

        def get_msg():
            return seq.pop(0) if seq else B._EMPTY

        return get_msg

    def test_happy_path_collects_all_runs(self):
        runs = [
            (B.MSG_RUN, {"total": 1.0, "stages": {}, "atoms": 5, "xyz_sha": "x"}) for _ in range(3)
        ]
        out = B._collect_runs(
            self._msgs([*runs, (B.MSG_DONE,)]),
            lambda: True,
            3,
            240.0,
            FakeClock(step=0.0),
            initial_oin="OIN",
        )
        self.assertEqual(len(out["records"]), 3)
        self.assertFalse(out["timed_out"])
        self.assertIsNone(out["error"])
        self.assertEqual(out["oin"], "OIN")

    def test_encoded_updates_oin(self):
        seq = [
            (B.MSG_ENCODED, "ENCODED-OIN"),
            (B.MSG_RUN, {"total": 1.0, "stages": {}}),
            (B.MSG_DONE,),
        ]
        out = B._collect_runs(self._msgs(seq), lambda: True, 1, 240.0, FakeClock(0.0), None)
        self.assertEqual(out["oin"], "ENCODED-OIN")
        self.assertEqual(len(out["records"]), 1)

    def test_timeout_is_recorded_not_raised(self):
        # get_msg always empty, child alive, clock marches past the deadline.
        out = B._collect_runs(
            self._msgs([]),
            lambda: True,
            5,
            mol_timeout=3.0,
            monotonic=FakeClock(step=1.0),
            initial_oin="O",
        )
        self.assertTrue(out["timed_out"])
        self.assertEqual(out["records"], [])
        self.assertIsNone(out["error"])

    def test_child_death_recorded_as_error(self):
        out = B._collect_runs(
            self._msgs([]), lambda: False, 5, 240.0, FakeClock(step=0.0), initial_oin="O"
        )
        self.assertIn("died", out["error"])
        self.assertFalse(out["timed_out"])

    def test_child_error_message_propagated(self):
        seq = [(B.MSG_ERROR, {"type": "ValueError", "msg": "no conformers"})]
        out = B._collect_runs(self._msgs(seq), lambda: True, 5, 240.0, FakeClock(0.0), "O")
        self.assertIn("ValueError", out["error"])
        self.assertIn("no conformers", out["error"])
        self.assertEqual(out["records"], [])


if __name__ == "__main__":
    unittest.main()
