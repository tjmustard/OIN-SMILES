"""Guards for the AC2BO candidate-generation memo (v0.4.5 encspeed lane).

The memo is what makes a slow encode ~2x faster, and it is also the one piece of new
machinery that could silently return a *wrong* cached answer. These tests pin the three
properties the safety argument rests on:

1. ``get_UA_pairs`` / ``get_bonds`` return the same values with the memo warm as cold.
2. A cache hit hands back a **fresh mutable** object -- callers mutate what they get
   (``get_UA_pairs`` appends virtual-node edges to ``get_bonds``' list), so a shared
   object would poison every later hit.
3. The memo is only read for an array the cache holds a **live reference** to, so a
   different matrix with a recycled ``id()`` can never alias a stale entry.

Plus an end-to-end check that a repeated encode in one process is byte-identical -- the
cross-encode reuse is the whole point of the LRU, so it must not drift.
"""

import unittest
from pathlib import Path

import numpy as np

from oinsmiles.utils import perception_core as loc

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"


def _ring_ac(n=8):
    """Adjacency matrix of an n-membered ring: enough unsaturated atoms to matter."""
    AC = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        AC[i, (i + 1) % n] = 1
        AC[(i + 1) % n, i] = 1
    return AC


class TestAC2BOMemo(unittest.TestCase):
    def setUp(self):
        loc._ac2bo_memo_clear()

    def tearDown(self):
        loc._ac2bo_memo_clear()

    def test_warm_matches_cold_for_get_UA_pairs(self):
        AC = _ring_ac(8)
        UA = list(range(8))
        DU = [1] * 8

        cold = loc.get_UA_pairs(UA, AC, DU)
        loc._ac2bo_memo_anchor(AC)
        warm_miss = loc.get_UA_pairs(UA, AC, DU)  # populates
        warm_hit = loc.get_UA_pairs(UA, AC, DU)  # reads back

        self.assertEqual(sorted(map(sorted, cold[0])), sorted(map(sorted, warm_miss[0])))
        self.assertEqual(sorted(map(sorted, cold[0])), sorted(map(sorted, warm_hit[0])))

    def test_warm_matches_cold_for_get_bonds(self):
        AC = _ring_ac(6)
        UA = [0, 1, 2, 3, 4, 5]

        cold = loc.get_bonds(UA, AC)
        loc._ac2bo_memo_anchor(AC)
        loc.get_bonds(UA, AC)
        warm_hit = loc.get_bonds(UA, AC)
        self.assertEqual(cold, warm_hit)

    def test_hit_returns_a_fresh_mutable_object(self):
        """Callers mutate the returned lists; a shared object would poison later hits."""
        AC = _ring_ac(6)
        UA = [0, 1, 2, 3, 4, 5]
        loc._ac2bo_memo_anchor(AC)

        first = loc.get_bonds(UA, AC)
        first.append(("poison", "poison"))
        second = loc.get_bonds(UA, AC)
        self.assertNotIn(("poison", "poison"), second)
        self.assertIsNot(first, second)

        pairs_a = loc.get_UA_pairs(UA, AC, [1] * 6)
        pairs_a[0].append(("poison", "poison"))
        pairs_b = loc.get_UA_pairs(UA, AC, [1] * 6)
        self.assertNotIn(("poison", "poison"), pairs_b[0])
        self.assertIsNot(pairs_a[0], pairs_b[0])

    def test_different_matrix_is_not_served_from_the_slot(self):
        """A second AC with different contents must miss, not alias the first slot."""
        AC1 = _ring_ac(6)
        AC2 = _ring_ac(6)
        AC2[0, 3] = 1  # add a cross-ring bond -> genuinely different connectivity
        AC2[3, 0] = 1
        UA = [0, 1, 2, 3, 4, 5]

        loc._ac2bo_memo_anchor(AC1)
        bonds1 = loc.get_bonds(UA, AC1)
        loc._ac2bo_memo_anchor(AC2)
        bonds2 = loc.get_bonds(UA, AC2)

        self.assertNotEqual(sorted(bonds1), sorted(bonds2))
        self.assertIn((0, 3), sorted(bonds2))

    def test_unregistered_array_bypasses_the_memo(self):
        """An array no slot holds must compute directly, never read another slot."""
        AC1 = _ring_ac(6)
        loc._ac2bo_memo_anchor(AC1)
        loc.get_bonds([0, 1, 2, 3, 4, 5], AC1)

        AC2 = _ring_ac(6)  # identical contents, different object, NOT anchored
        self.assertEqual((None, None), loc._ac2bo_memo_for(AC2))

    def test_identical_contents_adopt_the_existing_entries(self):
        """A re-perceived conformer with the same connectivity should reuse the cache."""
        AC1 = _ring_ac(7)
        loc._ac2bo_memo_anchor(AC1)
        loc.get_bonds([0, 1, 2, 3], AC1)
        entries_before = loc._ac2bo_memo_entries()
        self.assertGreater(entries_before, 0)

        AC2 = _ring_ac(7)  # same bytes, different array object
        loc._ac2bo_memo_anchor(AC2)
        self.assertEqual(entries_before, loc._ac2bo_memo_entries())
        bonds, uap = loc._ac2bo_memo_for(AC2)
        self.assertIsNotNone(bonds)
        self.assertIsNotNone(uap)

    def test_slot_count_is_bounded(self):
        for n in range(4, 4 + loc._AC2BO_MEMO_SLOTS + 3):
            loc._ac2bo_memo_anchor(_ring_ac(n))
        self.assertLessEqual(len(loc._AC2BO_SLOTS), loc._AC2BO_MEMO_SLOTS)


class TestRepeatEncodeIsStable(unittest.TestCase):
    """A second encode in the same process reads the warm memo -- it must not drift."""

    def test_repeat_encode_byte_identical(self):
        from oinsmiles import XYZToSMILES

        loc._ac2bo_memo_clear()
        path = str(FIXTURES / "Ferrocene.xyz")
        first = XYZToSMILES().convert(path)
        second = XYZToSMILES().convert(path)
        third = XYZToSMILES().convert(path)
        self.assertEqual(first, second)
        self.assertEqual(first, third)


if __name__ == "__main__":
    unittest.main()
