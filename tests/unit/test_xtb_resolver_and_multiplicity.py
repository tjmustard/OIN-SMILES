"""Guard tests for A1 (v0.4.3): xtb binary resolution + weak-field multiplicity.

Two independent fixes that must land together (see
``spec/handoffs/v0.4.3/A1-xtb-multiplicity.md``):

* ``resolve_xtb_binary`` makes the project's ``xtb`` discoverable under the harness
  interpreter (env override -> PATH -> interpreter-adjacent) without breaking the
  fail-soft path when ``xtb`` is genuinely absent.
* ``_spin_multiplicity`` replaces the blanket low-spin forcing at
  ``get_om_from_modified_smiles`` with a bounded high-spin assignment for
  unambiguously weak-field first-row complexes, while leaving strong-field /
  non-cohort complexes (and therefore every golden) unchanged.
"""

import os
import stat
import tempfile
import unittest
from unittest import mock

import numpy as np

from oinsmiles.generator3d import ml_optimizer, om
from oinsmiles.generator3d.ml_optimizer import resolve_xtb_binary
from oinsmiles.generator3d.om import _spin_multiplicity


def _make_executable(dirpath, name="xtb"):
    """Create an executable stub file and return its path."""
    path = os.path.join(dirpath, name)
    with open(path, "w") as fh:
        fh.write("#!/bin/sh\n")
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


class TestXtbResolver(unittest.TestCase):
    """resolve_xtb_binary() resolution order and fail-soft contract."""

    def test_returns_none_when_absent_everywhere(self):
        # No override, nothing on PATH, nothing beside the interpreter -> None, so
        # the caller degrades to the force-field geometry (fail-soft preserved).
        with (
            tempfile.TemporaryDirectory() as empty,
            mock.patch.dict(os.environ, {}, clear=False),
            mock.patch.object(ml_optimizer.shutil, "which", return_value=None),
            mock.patch.object(ml_optimizer.sys, "executable", os.path.join(empty, "python")),
        ):
            os.environ.pop("OIN_XTB_BIN", None)
            self.assertIsNone(resolve_xtb_binary())

    def test_interpreter_adjacent_found_when_not_on_path(self):
        # The venv-split case: xtb sits next to the interpreter but its bin/ is not
        # on PATH. shutil.which misses it; the resolver must still find it.
        with (
            tempfile.TemporaryDirectory() as venv_bin,
            mock.patch.dict(os.environ, {}, clear=False),
            mock.patch.object(ml_optimizer.shutil, "which", return_value=None),
            mock.patch.object(ml_optimizer.sys, "executable", os.path.join(venv_bin, "python")),
        ):
            os.environ.pop("OIN_XTB_BIN", None)
            exe = _make_executable(venv_bin)
            self.assertEqual(resolve_xtb_binary(), exe)

    def test_env_override_wins(self):
        with tempfile.TemporaryDirectory() as d:
            exe = _make_executable(d, "custom-xtb")
            with mock.patch.dict(os.environ, {"OIN_XTB_BIN": exe}):
                self.assertEqual(resolve_xtb_binary(), exe)

    def test_bad_override_falls_through(self):
        # A non-executable / missing override must not shadow a real xtb on PATH.
        with (
            mock.patch.dict(os.environ, {"OIN_XTB_BIN": "/does/not/exist/xtb"}),
            mock.patch.object(ml_optimizer.shutil, "which", return_value="/usr/bin/xtb"),
        ):
            self.assertEqual(resolve_xtb_binary(), "/usr/bin/xtb")


class TestWeakFieldMultiplicity(unittest.TestCase):
    """Bounded weak-field high-spin assignment via get_om_from_modified_smiles."""

    def _complex(self, msmiles):
        return om.get_om_from_modified_smiles(msmiles)

    def _total_electrons(self, metal_complex):
        mol = metal_complex.get_molecule()
        return int(np.sum(mol.get_z_list())) - metal_complex.chg

    def test_mn_ii_hexafluoride_is_sextet(self):
        # Mn(II) d5 octahedral high-spin -> sextet (5 unpaired). The pre-fix formula
        # forced the parity minimum (doublet, mult=2).
        mc = self._complex("[Mn]|[F:1]|[F:2]|[F:3]|[F:4]|[F:5]|[F:6]|6_octahedral")
        self.assertEqual(mc.multiplicity, 6)

    def test_ni_ii_hexaaqua_is_triplet(self):
        # Ni(II) d8 octahedral -> triplet (2 unpaired). Pre-fix: singlet (mult=1).
        mc = self._complex("[Ni]|[OH2:1]|[OH2:2]|[OH2:3]|[OH2:4]|[OH2:5]|[OH2:6]|6_octahedral")
        self.assertEqual(mc.multiplicity, 3)

    def test_co_ii_tetrachloride_is_quartet(self):
        # Co(II) d7 tetrahedral high-spin -> quartet (3 unpaired). Pre-fix: doublet.
        mc = self._complex("[Co]|[Cl:1]|[Cl:2]|[Cl:3]|[Cl:4]|4_tetrahedral")
        self.assertEqual(mc.multiplicity, 4)

    def test_assigned_uhf_matches_electron_parity(self):
        # The .UHF an optimizer consumes (multiplicity - 1) must have the same parity
        # as the drawn electron count, or the spin state is impossible.
        for msmiles in (
            "[Mn]|[F:1]|[F:2]|[F:3]|[F:4]|[F:5]|[F:6]|6_octahedral",
            "[Ni]|[OH2:1]|[OH2:2]|[OH2:3]|[OH2:4]|[OH2:5]|[OH2:6]|6_octahedral",
            "[Co]|[Cl:1]|[Cl:2]|[Cl:3]|[Cl:4]|4_tetrahedral",
        ):
            mc = self._complex(msmiles)
            uhf = mc.multiplicity - 1
            self.assertEqual(uhf % 2, self._total_electrons(mc) % 2, msmiles)

    def test_strong_field_carbon_donor_unchanged(self):
        # CO donors are strong-field: FeCO5 stays the low-spin singlet it always was.
        mc = self._complex("[Fe]|O#[C:1]|O#[C:2]|O#[C:3]|O#[C:4]|O#[C:5]|5_trigonal_bipyramidal")
        self.assertEqual(mc.multiplicity, 1)

    def test_non_cohort_metal_unchanged(self):
        # Pt is a third-row metal outside the first-row cohort -> low-spin default.
        mc = self._complex("[Pt]|[Cl:1]|[NH3:2]|[Cl:3]|[NH3:4]|4_square_planar")
        self.assertEqual(mc.multiplicity, 1)

    def test_mixed_donor_set_not_bumped(self):
        # Halide + N is not 'unambiguously weak-field' (N is excluded) -> low-spin.
        mc = self._complex("[Fe]|[Cl:1]|[Cl:2]|[NH3:3]|[NH3:4]|4_tetrahedral")
        self.assertEqual(mc.multiplicity, 1)


class TestSpinMultiplicityHelper(unittest.TestCase):
    """Unit-level properties of the _spin_multiplicity helper."""

    def test_bare_metal_falls_back_to_low_spin(self):
        # No donors -> low-spin parity default (Fe z=26, chg 0 -> singlet).
        self.assertEqual(_spin_multiplicity("Fe", 26, 0, []), 1)

    def test_never_lowers_spin(self):
        # The high-spin branch may only raise (or match) the low-spin multiplicity.
        for msmiles in (
            "[Mn]|[F:1]|[F:2]|[F:3]|[F:4]|[F:5]|[F:6]|6_octahedral",
            "[Ni]|[OH2:1]|[OH2:2]|[OH2:3]|[OH2:4]|[OH2:5]|[OH2:6]|6_octahedral",
            "[Cu]|[Cl:1]|[Cl:2]|[Cl:3]|[Cl:4]|4_tetrahedral",
        ):
            mc = om.get_om_from_modified_smiles(msmiles)
            total = int(np.sum(mc.get_molecule().get_z_list())) - mc.chg
            self.assertGreaterEqual(mc.multiplicity, total % 2 + 1, msmiles)


if __name__ == "__main__":
    unittest.main()
