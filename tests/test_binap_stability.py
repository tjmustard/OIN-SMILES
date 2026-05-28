"""BINAP stability test: graceful handling of axial-chiral ligands.

BINAP (2,2'-bis(diphenylphosphino)-1,1'-binaphthyl) exhibits axial chirality
(atropisomerism), not P-centered chirality. This test verifies that the pipeline
handles axial-chiral compounds without crashing and produces non-None output.

**Constraint:** No @/@@ assertion on BINAP P atoms — the pipeline does not
interpret axial chirality as P-centered stereochemistry.
"""

import unittest

from oinsmiles import XYZToSMILES
from .test_helpers import get_fixture_path


class TestBINAPStability(unittest.TestCase):
    """Verify graceful handling of axial-chiral BINAP ligands."""

    def test_binap_does_not_crash(self):
        """PdCl₂-R-BINAP should encode without crashing, returning non-None OIN."""
        xyz_path = get_fixture_path("PdCl2-R-BINAP.xyz")
        result = XYZToSMILES().convert(str(xyz_path))
        self.assertIsNotNone(result)
        self.assertIsInstance(result, str)
        self.assertGreater(len(result), 0)

    def test_binap_oin_is_valid(self):
        """Encoded BINAP OIN should be a valid OIN string (contains metal + ligands)."""
        xyz_path = get_fixture_path("PdCl2-R-BINAP.xyz")
        result = XYZToSMILES().convert(str(xyz_path))
        # Basic validation: should contain metal marker and slot markers
        self.assertIn("[Pd", result)
        self.assertIn("}", result)  # Slot markers present


if __name__ == "__main__":
    unittest.main()
