"""E2_P_FRAGILE (v0.4.19 E2 lane): three renumbering leaks found by bisecting ONE molecule.

XIVMEX_comp_0 (Zn, 5,15-ditolyl-porphyrin dianion, 63 atoms) is a census `E2_P_FRAGILE` row:
its string changes under every renumbering the canonicality audit applies. Stage-by-stage
bisection of `get_tmc_mol` on two numberings found the order-dependence entering at three
places, each behind its own held-off lever:

* `OIN_N_VALENCE_2` -- `atomic_valence[7] = [3, 4]`, so the valence search cannot write an
  N(-) and every candidate at charge -2 fails; `AC2BO` returns `best_BO` = the AC itself (zero
  double bonds) and the string is whatever RDKit's resonance enumeration makes of the
  zwitterion `set_atomic_charges` builds from it.
* `OIN_CANONICAL_CHARGES` -- that zwitterion is numbering-dependent, because the charge walk
  reads a running total in INPUT order while `charge_is_OK` accepted the BO in canonical order.
* `OIN_CANONICAL_RESONANCE` -- `ResonanceMolSupplier` returns a numbering-dependent SUBSET of
  forms (167 vs 131 on the same molecule), so the sorted-first selection still moves.

The shipped defect is written as a test on purpose, so a default flip has to change this file.
The renumbering is the audit's own (`tools/census/e_selfconsistency.py`, seed 42, tag renum0).
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
import zlib
from unittest import mock

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem  # noqa: E402

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.oin.levers import _DEFAULT_ON, _HELD_OFF  # noqa: E402
from oinsmiles.utils import perception_core as pc  # noqa: E402
from oinsmiles.utils import perception_tmc as pt  # noqa: E402

XIVMEX = os.path.join(os.path.dirname(__file__), "..", "fixtures", "XIVMEX_comp_0.xyz")
LEVERS = ("OIN_N_VALENCE_2", "OIN_CANONICAL_CHARGES", "OIN_CANONICAL_RESONANCE")
OFF = {k: "0" for k in LEVERS}


def _read_xyz(path):
    with open(path) as fh:
        n = int(fh.readline())
        comment = fh.readline().rstrip("\n")
        syms, xyz = [], []
        for _ in range(n):
            s, x, y, z = fh.readline().split()[:4]
            syms.append(s)
            xyz.append((float(x), float(y), float(z)))
    return syms, np.array(xyz), comment


def _renumbered_copy(tmpdir, tag="renum0"):
    """The audit's renum0 variant of XIVMEX: same seeds, same file name."""
    syms, xyz, comment = _read_xyz(XIVMEX)
    rng = np.random.default_rng([42, zlib.crc32(b"XIVMEX_comp_0"), zlib.crc32(tag.encode())])
    perm = rng.permutation(len(syms))
    while np.array_equal(perm, np.arange(len(syms))):
        perm = rng.permutation(len(syms))
    path = os.path.join(tmpdir, "XIVMEX_comp_0.xyz")
    with open(path, "w") as fh:
        fh.write(f"{len(syms)}\n{comment}\n")
        for i in perm:
            fh.write(f"{syms[i]:<3} {xyz[i][0]:>14.8f} {xyz[i][1]:>14.8f} {xyz[i][2]:>14.8f}\n")
    return path


class TestRegistry(unittest.TestCase):
    def test_all_three_held_off_with_the_mechanism(self):
        for k in LEVERS:
            self.assertIn(k, _HELD_OFF)
            self.assertNotIn(k, _DEFAULT_ON)
        self.assertIn("atomic_valence[7] = [3, 4]", _HELD_OFF["OIN_N_VALENCE_2"])
        self.assertIn("running total", _HELD_OFF["OIN_CANONICAL_CHARGES"])
        self.assertIn("ResonanceMolSupplier", _HELD_OFF["OIN_CANONICAL_RESONANCE"])


class TestShippedDefect(unittest.TestCase):
    def test_renumbering_moves_the_string_and_writes_radicals(self):
        with tempfile.TemporaryDirectory() as td, mock.patch.dict(os.environ, OFF):
            base = XYZToSMILES().convert(XIVMEX)
            moved = XYZToSMILES().convert(_renumbered_copy(td))
        self.assertNotEqual(base, moved)
        self.assertNotIn("[CH]", base)  # the file's numbering happens to get the dianion
        self.assertIn("[CH]", moved)  # the renumbering gets radicals

    def test_nitrogen_cannot_be_valence_2_in_the_search(self):
        with mock.patch.dict(os.environ, OFF):
            self.assertEqual(pc.possible_valences([2], [7]), [[3, 4]])
        with mock.patch.dict(os.environ, {"OIN_N_VALENCE_2": "1"}):
            self.assertEqual(pc.possible_valences([2], [7]), [[3, 4, 2]])  # appended LAST
            self.assertEqual(pc.possible_valences([3], [7]), [[3, 4]])  # three-coordinate: no


class TestLeversTogether(unittest.TestCase):
    def test_four_numberings_one_string(self):
        on = {"OIN_N_VALENCE_2": "1", "OIN_CANONICAL_RESONANCE": "1", "OIN_CANONICAL_CHARGES": "1"}
        with tempfile.TemporaryDirectory() as td, mock.patch.dict(os.environ, on):
            base = XYZToSMILES().convert(XIVMEX)
            others = []
            for t in ("renum0", "renum1", "renum2"):
                os.makedirs(os.path.join(td, t))
                others.append(XYZToSMILES().convert(_renumbered_copy(os.path.join(td, t), t)))
        for s in others:
            self.assertEqual(base, s)
        self.assertNotIn("[CH]", base)
        self.assertIn("n{0}", base)  # four N donors, the dianion
        self.assertEqual(base.count("n{"), 4, base)

    def test_n_valence_alone_fixes_the_chemistry_not_the_stability(self):
        on = {"OIN_N_VALENCE_2": "1", "OIN_CANONICAL_RESONANCE": "0", "OIN_CANONICAL_CHARGES": "0"}
        with tempfile.TemporaryDirectory() as td, mock.patch.dict(os.environ, on):
            base = XYZToSMILES().convert(XIVMEX)
            moved = XYZToSMILES().convert(_renumbered_copy(td, "renum1"))
        self.assertNotIn("[CH]", base)
        self.assertNotIn("[CH]", moved)  # a real dianion from both numberings ...
        self.assertNotEqual(base, moved)  # ... but a different Kekule form: the supplier leak


class TestCanonicalResonanceFrame(unittest.TestCase):
    def test_forms_come_back_in_the_callers_numbering(self):
        lig = Chem.AddHs(Chem.MolFromSmiles("c1cc[n-]c1"))
        canon, back = pt._canonical_resonance_frame(lig)
        self.assertIsNotNone(back)
        self.assertEqual(canon.GetNumAtoms(), lig.GetNumAtoms())
        restored = Chem.RenumberAtoms(canon, back)
        for a, b in zip(restored.GetAtoms(), lig.GetAtoms()):
            self.assertEqual(a.GetAtomicNum(), b.GetAtomicNum())
            self.assertEqual(a.GetDegree(), b.GetDegree())

    def test_charge_walk_order_reproduces_the_acceptance(self):
        # A BO accepted by charge_is_OK in canonical order is charged to exactly `charge`
        # when set_atomic_charges walks the same order (AC2mol's formal-charge check).
        with mock.patch.dict(os.environ, {"OIN_N_VALENCE_2": "0", "OIN_CANONICAL_CHARGES": "1"}):
            mol, _ = pt.get_basic_mol(XIVMEX, 0)
            tm = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 30)
            rw = Chem.RWMol(mol)
            for nb in list(rw.GetAtomWithIdx(tm).GetNeighbors()):
                rw.RemoveBond(tm, nb.GetIdx())
            rw.RemoveAtom(tm)
            lig = rw.GetMol()
            atoms = [a.GetAtomicNum() for a in lig.GetAtoms()]
            AC = Chem.GetAdjacencyMatrix(lig)
            got = pc.AC2mol(lig, AC, atoms, -2, allow_charged_fragments=True, use_atom_maps=False)
        self.assertIsNotNone(got)
        self.assertEqual(Chem.GetFormalCharge(got), -2)


if __name__ == "__main__":
    unittest.main()
