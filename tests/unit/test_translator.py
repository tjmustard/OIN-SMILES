import unittest
from unittest.mock import MagicMock

from oinsmiles.core.translator import SMILESToXYZ
from oinsmiles.generation.structure import GeneratedStructure

# Cisplatin PtCl2(NH3)2 -- one of the four golden fixtures.
CISPLATIN_OIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"


class TestSMILESToXYZDelegation(unittest.TestCase):
    """SMILESToXYZ is a thin public wrapper over OIN3DGenerator."""

    def test_convert_returns_generator_xyz(self):
        """convert() delegates to the generator and returns its .xyz block."""
        converter = SMILESToXYZ(optimizer="ff")
        converter.generator = MagicMock()
        converter.generator.generate.return_value = GeneratedStructure(
            xyz="5\n\nPt 0.0 0.0 0.0\n...", mol=None
        )

        out = converter.convert("fake_oin")

        converter.generator.generate.assert_called_once_with("fake_oin")
        self.assertEqual(out, "5\n\nPt 0.0 0.0 0.0\n...")

    def test_generate_returns_full_structure(self):
        """generate() returns the full GeneratedStructure (xyz + bonded mol)."""
        converter = SMILESToXYZ(optimizer="ff")
        converter.generator = MagicMock()
        sentinel = GeneratedStructure(xyz="x", mol=None)
        converter.generator.generate.return_value = sentinel

        self.assertIs(converter.generate("fake_oin"), sentinel)


class TestSMILESToXYZEndToEnd(unittest.TestCase):
    """End-to-end: the public reverse API produces a real 3D structure.

    Regression guard against the former stub, which emitted dummy ``X`` atoms at
    the origin instead of generating a conformer.
    """

    def test_convert_produces_real_xyz_block(self):
        # FF-only optimizer keeps this hermetic (no g-xTB/MACE binaries); the
        # fixed seed makes generation deterministic.
        converter = SMILESToXYZ(optimizer="ff", seed=42)
        xyz = converter.convert(CISPLATIN_OIN)

        self.assertIsInstance(xyz, str)
        lines = xyz.strip().splitlines()
        declared = int(lines[0])  # XYZ header: atom count
        atom_lines = lines[2:]  # skip count + comment lines
        self.assertEqual(len(atom_lines), declared)
        self.assertGreaterEqual(declared, 5)  # Pt + 2 Cl + 2 N at minimum

        elements = [ln.split()[0] for ln in atom_lines]
        self.assertIn("Pt", elements)
        self.assertIn("Cl", elements)
        self.assertNotIn("X", elements)  # the old stub emitted dummy 'X' atoms


if __name__ == "__main__":
    unittest.main()
