"""Unit tests for the vendored MetalloGen 3D engine (``oinsmiles.generator3d``).

These characterize the pure/near-pure helpers in ``chem``, ``process``,
``ml_optimizer`` and the package ``__init__`` so refactors of the vendored engine
(dead-code removal, the frag->PuLP charge/bond-order reroute, docstring/style
cleanup) stay behavior-preserving. They are fast and deterministic; the slow
end-to-end ``generate_3d_structures`` path is smoke-tested by
``tests/test_generator3d.py``.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

import numpy as np

from oinsmiles.generator3d import calculate_heavy_atom_rmsd, chem, get_xyz_string, process


class TestAtom(unittest.TestCase):
    def test_element_to_atomic_number(self):
        self.assertEqual(chem.Atom("C").get_atomic_number(), 6)
        self.assertEqual(chem.Atom("H").get_atomic_number(), 1)
        self.assertEqual(chem.Atom("Fe").get_atomic_number(), 26)

    def test_atomic_number_to_element(self):
        self.assertEqual(chem.Atom(6).get_element(), "C")
        self.assertEqual(chem.Atom(26).get_element(), "Fe")

    def test_set_atomic_number_sets_element(self):
        a = chem.Atom()
        a.set_atomic_number(8)
        self.assertEqual(a.element, "O")

    def test_set_element_sets_atomic_number(self):
        a = chem.Atom()
        a.set_element("Cl")
        self.assertEqual(a.atomic_number, 17)

    def test_get_atomic_number_normalizes_case(self):
        # Upper-cased multi-char symbols (e.g. "CL") normalize to "Cl" -> Z=17.
        a = chem.Atom("CL")
        self.assertEqual(a.get_atomic_number(), 17)
        self.assertEqual(a.element, "Cl")

    def test_is_same_atom_and_eq(self):
        self.assertTrue(chem.Atom("C").is_same_atom(chem.Atom("C")))
        self.assertFalse(chem.Atom("C").is_same_atom(chem.Atom("N")))
        self.assertTrue(chem.Atom("O") == chem.Atom("O"))
        self.assertFalse(chem.Atom("O") == chem.Atom("S"))

    def test_copy_is_independent(self):
        a = chem.Atom("C")
        a.set_coordinate([1.0, 2.0, 3.0])
        b = a.copy()
        self.assertEqual(b.get_element(), "C")
        b.x = 99.0
        self.assertEqual(a.x, 1.0)

    def test_set_coordinate_3d_and_2d(self):
        a = chem.Atom("C")
        a.set_coordinate([1.0, 2.0, 3.0])
        np.testing.assert_allclose(a.get_coordinate(), [1.0, 2.0, 3.0])
        a.set_coordinate([4.0, 5.0])
        np.testing.assert_allclose(a.get_coordinate(), [4.0, 5.0, 0.0])

    def test_get_content_element_and_z(self):
        a = chem.Atom("C")
        a.set_coordinate([0.0, 0.0, 0.0])
        self.assertTrue(a.get_content("element").startswith("C "))
        self.assertTrue(a.get_content("z").startswith("6 "))


class TestMolecule(unittest.TestCase):
    def _water(self):
        mol = chem.Molecule()
        mol.atom_list = [chem.Atom("O"), chem.Atom("H"), chem.Atom("H")]
        mol.adj_matrix = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]])
        mol.chg = 0
        return mol

    def test_get_z_list(self):
        self.assertEqual(list(self._water().get_z_list()), [8, 1, 1])

    def test_get_adj_matrix_roundtrip(self):
        mol = self._water()
        np.testing.assert_array_equal(mol.get_adj_matrix(), mol.adj_matrix)

    def test_get_chg(self):
        self.assertEqual(self._water().get_chg(), 0)


class TestProcessHelpers(unittest.TestCase):
    def test_z_list_roundtrip(self):
        z = [6, 1, 8]
        atoms = process.get_atom_list_from_z_list(z)
        self.assertEqual(process.get_z_list_from_atom_list(atoms), z)

    def test_element_list_roundtrip(self):
        els = ["C", "H", "O"]
        atoms = process.get_atom_list_from_element_list(els)
        self.assertEqual(process.get_element_list_from_atom_list(atoms), els)

    def test_copy_atom_list_is_independent(self):
        atoms = process.get_atom_list_from_element_list(["C", "O"])
        atoms[0].x = 1.0
        copied = process.copy_atom_list(atoms)
        copied[0].x = 42.0
        self.assertEqual(atoms[0].x, 1.0)

    def test_block_diagonal_adj(self):
        f1 = np.array([[0, 1], [1, 0]])
        f2 = np.array([[0, 1], [1, 0]])
        block = process.get_block_diagonal_adj_from_fragments([f1, f2])
        expected = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=float)
        np.testing.assert_array_equal(block, expected)

    def test_group_molecules_two_components(self):
        # Two disjoint 2-atom molecules -> two groups of indices.
        adj = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
        groups = process.group_molecules(adj)
        self.assertEqual(len(groups), 2)
        self.assertEqual(sorted(sorted(g) for g in groups), [[0, 1], [2, 3]])

    def test_get_molecule_group_connected_component(self):
        adj = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
        self.assertEqual(process.get_molecule_group(adj, 0), {0, 1})
        self.assertEqual(process.get_molecule_group(adj, 2), {2})

    def test_check_geometry(self):
        well_separated = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        self.assertTrue(process.check_geometry(well_separated))
        overlapping = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]])
        self.assertFalse(process.check_geometry(overlapping))

    def test_check_atom_validity(self):
        self.assertTrue(process.check_atom_validity(group=4, bo=4, chg=0))
        # group - bo - chg < 0 -> invalid (negative lone pairs)
        self.assertFalse(process.check_atom_validity(group=1, bo=3, chg=0))

    def test_get_rmsd_identity_is_zero(self):
        mol = chem.Molecule()
        mol.atom_list = [chem.Atom("C"), chem.Atom("O")]
        mol.atom_list[0].set_coordinate([0.0, 0.0, 0.0])
        mol.atom_list[1].set_coordinate([1.2, 0.0, 0.0])
        self.assertAlmostEqual(process.get_rmsd(mol, mol), 0.0)


class TestChargeBondOrderReroute(unittest.TestCase):
    """The frag-based perceiver was rerouted to the PuLP solver; verify results."""

    def _water(self):
        mol = chem.Molecule()
        mol.atom_list = [chem.Atom("O"), chem.Atom("H"), chem.Atom("H")]
        mol.adj_matrix = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]])
        mol.chg = 0
        return mol

    def test_get_chg_and_bo_water(self):
        chg_list, bo = process.get_chg_and_bo(self._water(), 0)
        bo = np.asarray(bo)
        self.assertEqual(bo.shape, (3, 3))
        # Neutral water: all formal charges zero, two O-H single bonds.
        self.assertEqual([int(c) for c in np.asarray(chg_list)], [0, 0, 0])
        self.assertEqual(int(bo[0, 1]), 1)
        self.assertEqual(int(bo[0, 2]), 1)

    def test_get_bo_matrix_from_adj_matrix_shape(self):
        bo = process.get_bo_matrix_from_adj_matrix(self._water(), 0)
        self.assertEqual(np.asarray(bo).shape, (3, 3))

    def test_get_bo_matrix_resonance_returns_list(self):
        result = process.get_bo_matrix_from_adj_matrix(self._water(), 0, obtain_all_resonance=True)
        self.assertIsInstance(result, list)


class TestPackageHelpers(unittest.TestCase):
    def _mol(self, coords):
        mol = chem.Molecule()
        mol.atom_list = []
        for el, xyz in coords:
            a = chem.Atom(el)
            a.set_coordinate(xyz)
            mol.atom_list.append(a)
        mol.chg = 0
        mol.multiplicity = 1
        return mol

    def test_heavy_atom_rmsd_identical_is_zero(self):
        mol = self._mol([("C", [0.0, 0.0, 0.0]), ("O", [1.2, 0.0, 0.0])])
        self.assertAlmostEqual(calculate_heavy_atom_rmsd(mol, mol), 0.0, places=5)

    def test_heavy_atom_rmsd_mismatched_count_is_inf(self):
        m1 = self._mol([("C", [0.0, 0.0, 0.0]), ("O", [1.2, 0.0, 0.0])])
        m2 = self._mol([("C", [0.0, 0.0, 0.0])])
        self.assertEqual(calculate_heavy_atom_rmsd(m1, m2), float("inf"))

    def test_get_xyz_string_format(self):
        mol = self._mol([("O", [0.0, 0.0, 0.0]), ("H", [0.96, 0.0, 0.0])])
        xyz = get_xyz_string(mol)
        lines = xyz.splitlines()
        self.assertEqual(lines[0], "2")  # atom count
        self.assertTrue(lines[2].startswith("O"))
        self.assertTrue(lines[3].startswith("H"))


class TestASEOptimizer(unittest.TestCase):
    def test_xtb_constructs(self):
        from oinsmiles.generator3d.ml_optimizer import ASEOptimizer

        opt = ASEOptimizer(method="xtb")
        self.assertEqual(opt.method, "xtb")

    def test_method_is_case_insensitive(self):
        from oinsmiles.generator3d.ml_optimizer import ASEOptimizer

        self.assertEqual(ASEOptimizer(method="XTB").method, "xtb")

    def test_invalid_method_raises(self):
        from oinsmiles.generator3d.ml_optimizer import ASEOptimizer

        with self.assertRaises(ValueError):
            ASEOptimizer(method="not-a-real-optimizer")

    def test_mace_constructs_or_skips(self):
        from oinsmiles.generator3d.ml_optimizer import ASEOptimizer

        try:
            opt = ASEOptimizer(method="mace-omol-0-extra-large-1024")
        except ImportError:
            self.skipTest("mace-torch not installed (optional extra).")
        self.assertEqual(opt.method, "mace-omol-0-extra-large-1024")

    def test_xtb_optimize_returns_energy_without_calculator(self):
        # Regression: the g-xTB (subprocess) path must NOT call
        # atoms.get_potential_energy() -- that raises "Atoms object has no
        # calculator" because only the MACE path attaches an ASE calculator.
        import shutil

        from oinsmiles.generator3d.ml_optimizer import ASEOptimizer

        if shutil.which("xtb") is None:
            self.skipTest("xtb binary not on PATH.")
        mol = chem.Molecule()
        mol.atom_list = []
        for el, xyz in [
            ("O", [0.0, 0.0, 0.117]),
            ("H", [0.0, 0.757, -0.467]),
            ("H", [0.0, -0.757, -0.467]),
        ]:
            a = chem.Atom(el)
            a.set_coordinate(xyz)
            mol.atom_list.append(a)
        mol.chg = 0
        mol.multiplicity = 1
        success, energy, opt_mol = ASEOptimizer(method="xtb").optimize(mol)
        self.assertTrue(success)
        self.assertEqual(len(opt_mol.atom_list), 3)


if __name__ == "__main__":
    unittest.main()
