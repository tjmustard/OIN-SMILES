"""Tests for Bug 4 fix: SDF bond-parsing no longer crashes after a warning.

Before the fix, get_molecule_info_from_sdf() called logger.debug() for a
malformed bond record but then unconditionally called int(s) - 1 on the next
line, crashing with ValueError on the same invalid token.

After the fix, the malformed bond record is skipped with `continue` so the
rest of the file parses successfully.
"""

import logging
import os
import tempfile
import unittest

from oinsmiles.generator3d.process import get_molecule_info_from_sdf


def _write_sdf(content: str) -> str:
    """Write *content* to a temp SDF file and return its path."""
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".sdf", delete=False)
    f.write(content)
    f.close()
    return f.name


# Minimal valid SDF: 2 atoms, 1 good bond + 1 bad bond.
_SDF_ONE_BAD_BOND = (
    "test\n  test\n\n"
    "  2  2  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0\n"   # valid
    "  1  X  1  0  0  0\n"   # malformed — non-numeric end atom
    "M  END\n$$$$\n"
)

# SDF with only a bad bond (1 atom, 1 malformed bond).
_SDF_ONLY_BAD_BOND = (
    "test\n  test\n\n"
    "  1  1  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  X  1  0  0  0\n"
    "M  END\n$$$$\n"
)


class TestSdfBondParsingSkipsBadRecord(unittest.TestCase):
    def setUp(self):
        logging.getLogger("oinsmiles.generator3d.process").setLevel(logging.DEBUG)

    def test_no_value_error_on_malformed_bond(self):
        """A malformed bond record must not raise ValueError."""
        path = _write_sdf(_SDF_ONLY_BAD_BOND)
        try:
            # Pre-fix: ValueError: invalid literal for int() with base 10: 'X'
            try:
                get_molecule_info_from_sdf(path)
            except ValueError as e:
                self.fail(f"get_molecule_info_from_sdf() raised ValueError: {e}")
        finally:
            os.unlink(path)

    def test_valid_bond_still_recorded(self):
        """After skipping the bad record the valid bond must still be in adj_matrix."""
        path = _write_sdf(_SDF_ONE_BAD_BOND)
        try:
            _, _, adj_matrix, _, _ = get_molecule_info_from_sdf(path)
            # atoms 0 and 1 are bonded by the valid record
            self.assertEqual(adj_matrix[0][1], 1)
            self.assertEqual(adj_matrix[1][0], 1)
        finally:
            os.unlink(path)

    def test_warning_logged_for_bad_bond(self):
        """The library must still log a warning for the skipped record."""
        path = _write_sdf(_SDF_ONLY_BAD_BOND)
        try:
            with self.assertLogs("oinsmiles.generator3d.process", level=logging.DEBUG) as cm:
                get_molecule_info_from_sdf(path)
            combined = " ".join(cm.output)
            self.assertIn("WRONG SDF", combined)
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()

