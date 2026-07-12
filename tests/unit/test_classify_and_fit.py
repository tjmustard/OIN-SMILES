"""Unit tests for the single-pass ``classify_and_fit`` geometry helper (P5).

``classify_and_fit(vectors, target)`` folds the two-call
``classify_coordination_geometry`` + ``coordination_geometry_fit`` sequence the
geometry-aware conformer selector used into ONE permutation/Kabsch match. These
tests pin two properties:

  * **Byte-identity** -- ``classify_and_fit`` returns exactly the
    ``(classify_coordination_geometry(v), coordination_geometry_fit(v, target))``
    pair the old sequence produced, across coordination numbers 2-7, for in-CN
    targets, cross-CN targets (fallback branch), and an unknown code. It is the
    same computation, so equality is exact, not approximate.
  * **Single pass** -- for a conformer whose classified label equals the target,
    ``_map_to_template`` runs once per candidate template, not twice for the
    target (the redundancy: 20 -> 10 calls on Ir(ppy)3's CN-6 sphere).
"""

import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.utils.oin_aligner import (
    TEMPLATES,
    OINDiscreteAligner,
    classify_and_fit,
    classify_coordination_geometry,
    coordination_geometry_fit,
)

# Candidate template lists, keyed by coordination number -- must mirror
# ``_match_geometry_candidates``.
CN_CANDIDATES = {
    2: ["LIN"],
    3: ["TPL"],
    4: ["SPL", "TET", "TPY"],
    5: ["TBP", "SPY"],
    6: ["OCT"],
    7: ["PBP"],
}


def _perturb(vectors, scale):
    """Deterministic, RNG-free perturbation so fits are non-trivial but reproducible."""
    v = np.asarray(vectors, dtype=float)
    off = np.array(
        [[((i * 7 + j * 3) % 5) - 2 for j in range(3)] for i in range(len(v))],
        dtype=float,
    )
    return v + scale * off


def _assert_float_identical(test, a, b, msg):
    """Exact equality, treating inf/inf and nan/nan as equal."""
    if math.isinf(a) or math.isinf(b):
        test.assertEqual(a, b, msg)
    elif math.isnan(a) or math.isnan(b):
        test.assertTrue(math.isnan(a) and math.isnan(b), msg)
    else:
        test.assertEqual(a, b, msg)  # same permutation/SVD path -> bit-identical


class TestByteIdentity(unittest.TestCase):
    def test_matches_old_two_call_sequence(self):
        for cn, cands in CN_CANDIDATES.items():
            for base_name in cands:
                base = TEMPLATES[base_name]
                for scale in (0.0, 0.15, 0.4):
                    v = _perturb(base, scale)
                    old_label = classify_coordination_geometry(v)
                    # in-CN targets, a cross-CN target (fallback branch), unknown code
                    for tgt in cands + ["OCT", "ZZZ"]:
                        old_fit = coordination_geometry_fit(v, tgt)
                        new_label, new_fit = classify_and_fit(v, tgt)
                        ctx = f"cn={cn} base={base_name} scale={scale} tgt={tgt}"
                        self.assertEqual(new_label, old_label, f"label: {ctx}")
                        _assert_float_identical(self, new_fit, old_fit, f"fit: {ctx}")

    def test_unknown_code_is_infinite(self):
        v = TEMPLATES["SPL"]
        label, fit = classify_and_fit(v, "ZZZ")
        self.assertEqual(label, "SPL")
        self.assertEqual(fit, float("inf"))

    def test_ideal_target_fit_is_near_zero(self):
        label, fit = classify_and_fit(TEMPLATES["OCT"], "OCT")
        self.assertEqual(label, "OCT")
        self.assertLess(fit, 1e-6)


class _MapSpy:
    """Context manager counting ``_map_to_template`` invocations."""

    def __init__(self):
        self.n = 0

    def __enter__(self):
        self._orig = OINDiscreteAligner._map_to_template
        spy = self

        def wrapped(inner_self, va, tv):
            spy.n += 1
            return spy._orig(inner_self, va, tv)

        OINDiscreteAligner._map_to_template = wrapped
        return self

    def __exit__(self, *exc):
        OINDiscreteAligner._map_to_template = self._orig
        return False


class TestSinglePass(unittest.TestCase):
    def test_cn6_target_matched_once_not_twice(self):
        # The Ir(ppy)3 case: one CN-6 candidate (OCT). classify_and_fit matches it
        # once; the old classify + fit sequence matched it twice.
        v = TEMPLATES["OCT"]
        with _MapSpy() as combined:
            classify_and_fit(v, "OCT")
        with _MapSpy() as two_call:
            classify_coordination_geometry(v)
            coordination_geometry_fit(v, "OCT")
        self.assertEqual(combined.n, 1, "single-pass should match OCT once")
        self.assertEqual(two_call.n, 2, "old sequence matched OCT twice")

    def test_cn4_target_not_matched_twice(self):
        # CN-4 evaluates 3 candidates. classify_and_fit -> 3 matches; the old
        # sequence -> 3 (classify) + 1 (fit of the target) = 4, i.e. the target
        # template SPL was matched a redundant second time.
        v = TEMPLATES["SPL"]
        with _MapSpy() as combined:
            classify_and_fit(v, "SPL")
        with _MapSpy() as two_call:
            classify_coordination_geometry(v)
            coordination_geometry_fit(v, "SPL")
        self.assertEqual(combined.n, 3)
        self.assertEqual(two_call.n, 4)


if __name__ == "__main__":
    unittest.main()
