"""Tests for Bug 2 fix: Ligand.__str__() returned None.

Before the fix __str__ was a bare ``pass``, causing Python to raise
TypeError: __str__ returned non-string (type NoneType) on str(lig).
After the fix __str__ returns a proper string.
"""

import unittest

from oinsmiles.generator3d.ligand import Ligand


class TestLigandStr(unittest.TestCase):
    def test_str_returns_string(self):
        lig = Ligand(molecule=None, binding_infos=[])
        # Pre-fix: TypeError: __str__ returned non-string (type NoneType)
        result = str(lig)
        self.assertIsInstance(result, str)

    def test_str_non_empty(self):
        lig = Ligand(molecule=None, binding_infos=[])
        self.assertTrue(len(str(lig)) > 0)

    def test_str_includes_denticity(self):
        """Denticity (number of binding infos) appears in the string."""
        lig_mono = Ligand(molecule=None, binding_infos=[(0, 1)])
        lig_bi = Ligand(molecule=None, binding_infos=[(0, 1), (1, 2)])
        self.assertIn("1", str(lig_mono))
        self.assertIn("2", str(lig_bi))

    def test_repr_usable_in_fstring(self):
        lig = Ligand(molecule=None, binding_infos=[])
        # Should not raise
        msg = f"created: {lig}"
        self.assertIsInstance(msg, str)


if __name__ == "__main__":
    unittest.main()

