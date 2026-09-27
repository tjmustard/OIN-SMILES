"""get_molecule_info_from_sdf() rejects a malformed bond record with a clear ValueError.

Before the fix, a non-numeric end index logged a debug warning and then crashed on
int() of the same token, and a non-numeric start index crashed on int() before the
check was reached. A malformed record now raises ValueError naming the record, rather
than being skipped: a dropped bond would return a silently wrong adjacency matrix.
"""

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


def _sdf(n_atoms: int, bond_lines: list) -> str:
    """A V2000 SDF of *n_atoms* carbons on a line, with the given bond-block lines."""
    atoms = "".join(
        f"{1.5 * k:10.4f}    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        for k in range(n_atoms)
    )
    counts = f"{n_atoms:3d}{len(bond_lines):3d}  0  0  0  0  0  0  0  0999 V2000\n"
    return "test\n  test\n\n" + counts + atoms + "".join(bond_lines) + "M  END\n$$$$\n"


class TestSdfBondParsing(unittest.TestCase):
    def _parse(self, content):
        path = _write_sdf(content)
        try:
            return get_molecule_info_from_sdf(path)
        finally:
            os.unlink(path)

    def test_well_formed_bond_is_recorded(self):
        _, _, adj_matrix, _, _ = self._parse(_sdf(2, ["  1  2  1  0  0  0\n"]))
        self.assertEqual(adj_matrix[0][1], 1)
        self.assertEqual(adj_matrix[1][0], 1)

    def test_bad_end_index_raises(self):
        content = _sdf(2, ["  1  2  1  0  0  0\n", "  1  X  1  0  0  0\n"])
        with self.assertRaisesRegex(ValueError, "malformed SDF bond record 2"):
            self._parse(content)

    def test_bad_start_index_raises_the_same_clear_error(self):
        # Previously int(s.strip()) crashed here with "invalid literal for int()".
        with self.assertRaisesRegex(ValueError, "malformed SDF bond record 1"):
            self._parse(_sdf(2, ["  X  2  1  0  0  0\n"]))

    def test_out_of_range_index_raises(self):
        # Index 0 used to become -1 and silently bond to the last atom.
        with self.assertRaisesRegex(ValueError, "malformed SDF bond record 1"):
            self._parse(_sdf(2, ["  0  2  1  0  0  0\n"]))

    def test_run_together_indices_past_99_still_parse(self):
        # V2000 fields are 3 wide, so atoms 100 and 101 are written "100101".
        _, _, adj_matrix, _, _ = self._parse(_sdf(101, ["100101  1  0  0  0\n"]))
        self.assertEqual(adj_matrix[99][100], 1)
        self.assertEqual(adj_matrix[100][99], 1)


if __name__ == "__main__":
    unittest.main()
