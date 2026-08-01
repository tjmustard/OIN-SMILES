"""Tests for Bug 5 fix: normalize_oin_for_comparison must fold '^' → '>'.

Before the fix, normalize_oin_for_comparison left '^' winding markers
untouched, so two OIN strings that differed only in '^' vs '>' compared
as unequal after normalization even though they describe the same structure.
_parse_vertex_colors (inside the key functions) folded it internally,
but the normalized *string* was not canonical.

After the fix, re.sub replaces {n^} with {n>} during normalization.
"""

import unittest

from oinsmiles.oin.compare import (
    canonical_roundtrip_key,
    normalize_oin_for_comparison,
    winding_canonical_key,
)

# [Fe_OCT] is a valid 3-letter geo code so the metal regex strips it correctly.
_S_GT = "[Fe_OCT].C1=CC=CC=C1{0>}.C1=CC=CC=C1{1<}"
_S_CARET = "[Fe_OCT].C1=CC=CC=C1{0^}.C1=CC=CC=C1{1<}"


class TestNormalizeCaretWinding(unittest.TestCase):
    def test_normalized_strings_are_equal(self):
        """After normalization, '^' and '>' variants must produce the same string."""
        n1 = normalize_oin_for_comparison(_S_GT)
        n2 = normalize_oin_for_comparison(_S_CARET)
        self.assertEqual(n1, n2)

    def test_caret_replaced_by_gt_in_normalized_form(self):
        """The normalized string must not contain '^' in any slot marker."""
        n2 = normalize_oin_for_comparison(_S_CARET)
        self.assertNotIn("^", n2)

    def test_keys_still_equal(self):
        """The canonical keys must remain equal (pre-existing behaviour)."""
        n1 = normalize_oin_for_comparison(_S_GT)
        n2 = normalize_oin_for_comparison(_S_CARET)
        self.assertEqual(winding_canonical_key(n1), winding_canonical_key(n2))

    def test_canonical_roundtrip_key_equal(self):
        self.assertEqual(canonical_roundtrip_key(_S_GT), canonical_roundtrip_key(_S_CARET))

    def test_non_winding_caret_unaffected(self):
        """A '^' that is NOT inside a slot marker must not be touched."""
        s = "[Fe_OCT].C^C{0>}"  # '^' outside a slot
        result = normalize_oin_for_comparison(s)
        self.assertIn("C^C", result)

    def test_idempotent(self):
        """Normalizing an already-normalized string must be a no-op."""
        n1 = normalize_oin_for_comparison(_S_GT)
        self.assertEqual(n1, normalize_oin_for_comparison(n1))


if __name__ == "__main__":
    unittest.main()

