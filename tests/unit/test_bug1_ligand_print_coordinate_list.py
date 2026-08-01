"""Tests for Bug 1 fix: Ligand.print_coordinate_list() AttributeError.

Before the fix the method called self.coordinate_list() and
self.get_coordinate_list(), neither of which exists, raising AttributeError.
After the fix it delegates to self.molecule.atom_list and
self.molecule.get_coordinate_list().
"""

import logging
import unittest
from unittest.mock import MagicMock

from oinsmiles.generator3d.ligand import Ligand


def _make_ligand_with_atoms(atoms):
    """Return a Ligand whose molecule is a mock with the given atom list.

    Each entry in *atoms* is an ``(element, x, y, z)`` tuple.
    """
    mock_atoms = []
    for element, x, y, z in atoms:
        a = MagicMock()
        a.get_element.return_value = element
        mock_atoms.append(a)

    import numpy as np

    coords = [[x, y, z] for _, x, y, z in atoms]

    mol = MagicMock()
    mol.atom_list = mock_atoms
    mol.get_coordinate_list.return_value = coords

    return Ligand(molecule=mol, binding_infos=[])


class TestPrintCoordinateListNoError(unittest.TestCase):
    """print_coordinate_list() must not raise after the fix."""

    def test_does_not_raise_attribute_error(self):
        lig = _make_ligand_with_atoms([("Cl", 2.0, 0.0, 0.0)])
        # Pre-fix: AttributeError: 'Ligand' object has no attribute 'coordinate_list'
        try:
            lig.print_coordinate_list()
        except AttributeError as e:
            self.fail(f"print_coordinate_list() raised AttributeError: {e}")

    def test_logs_each_atom(self):
        atoms = [("Pt", 0.0, 0.0, 0.0), ("Cl", 2.0, 0.0, 0.0)]
        lig = _make_ligand_with_atoms(atoms)

        with self.assertLogs("oinsmiles.generator3d.ligand", level=logging.DEBUG) as cm:
            lig.print_coordinate_list()

        # Both element symbols must appear somewhere in the debug output.
        combined = " ".join(cm.output)
        self.assertIn("Pt", combined)
        self.assertIn("Cl", combined)

    def test_delegates_to_molecule(self):
        """get_coordinate_list() on the molecule must be called exactly once."""
        lig = _make_ligand_with_atoms([("N", 1.0, 1.0, 1.0)])
        lig.print_coordinate_list()
        lig.molecule.get_coordinate_list.assert_called_once()


if __name__ == "__main__":
    unittest.main()

