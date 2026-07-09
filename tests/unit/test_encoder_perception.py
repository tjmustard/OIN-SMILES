"""Unit tests for bond-order / charge perception crashes in the XYZ->OIN encoder.

Three round-trip failure classes traced to ``xyz2mol``'s input-side perception:

* ``xyz2mol_none_crash`` -- ``lig_checks`` called ``len(res_mols)`` on a
  ``ResonanceMolSupplier`` and then iterated it. The supplier is a stateful
  iterator: ``len()`` runs the enumeration and leaves the cursor at the end, so
  the loop yielded ``None`` and ``res_mol.GetAtoms()`` raised ``AttributeError``.
  Indexing the supplier never returns ``None`` for the same ligand.

* ``garbled_aromatic`` (the nitro sub-bucket) -- ``get_oin_string`` zeroed every
  formal charge, so a nitro group's ``[N+](=O)[O-]`` became a four-bonded neutral
  nitrogen and serialized as the unparseable ``N(O)=O``.
"""

import os
import re
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem, RDLogger

from oinsmiles.core.translator import XYZToSMILES
from oinsmiles.utils.xyz2mol import lig_checks

RDLogger.DisableLog("rdApp.*")

_FIXTURES = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))


class TestResonanceSupplierGuard(unittest.TestCase):
    """`lig_checks` must survive a ResonanceMolSupplier that yields None."""

    def test_dixhok_encodes(self):
        # DIXHOK's bis(phosphine) ligand drove the supplier past its cursor:
        # XYZToSMILES.convert raised "'NoneType' object has no attribute
        # 'GetAtoms'" for this and 7 sibling molecules.
        oin = XYZToSMILES().convert(os.path.join(_FIXTURES, "DIXHOK.xyz"))
        self.assertTrue(oin.startswith("[Ni_SPL]"), oin)
        self.assertIn("P{0}", oin)

    def test_lig_checks_returns_a_candidate_for_every_ligand(self):
        # Every entry must be a real mol, never None, and the aromatic count must
        # be an int -- i.e. res_mol.GetAtoms() was reachable.
        mol = Chem.MolFromSmiles("COc1ccccc1P(CN(C)CP(c1ccccc1OC)c1ccccc1OC)c1ccccc1OC")
        self.assertIsNotNone(mol)
        results = lig_checks(mol, [])
        self.assertTrue(results)
        for res_mol, n_pos, n_neg, n_aromatic in results:
            self.assertIsNotNone(res_mol)
            self.assertGreater(res_mol.GetNumAtoms(), 0)
            self.assertIsInstance(n_aromatic, int)
            self.assertIsInstance(n_pos, int)
            self.assertIsInstance(n_neg, int)

    def test_lig_checks_falls_back_when_supplier_is_empty(self):
        # A ligand the supplier cannot enumerate must still yield one candidate
        # (the un-resonated mol) rather than an empty list, which would IndexError
        # in get_lig_mol's `possible_res_mols[0]`.
        mol = Chem.MolFromSmiles("[Cl-]")
        results = lig_checks(mol, [])
        self.assertEqual(len(results), 1)
        self.assertIsNotNone(results[0][0])


_METAL_FRAGMENT = re.compile(r"^\[[A-Z][a-z]?_[A-Z]{3}\]$")


def _ligand_fragments(oin):
    """Slot-stripped ligand fragments, as oin/compare.py sees them."""
    for fragment in oin.split("."):
        clean = re.sub(r"\{\d+[<>^]?\}", "", fragment)
        if clean and not _METAL_FRAGMENT.match(clean):
            yield clean


class TestChargeSeparatedGroups(unittest.TestCase):
    """Groups with no neutral Lewis structure keep their charges."""

    def test_nitro_survives_neutralization(self):
        oin = XYZToSMILES().convert(os.path.join(_FIXTURES, "ABUVUP.xyz"))
        # Before the fix this emitted `N(O)=O` -- a four-bonded neutral nitrogen.
        self.assertNotIn("N(O)=O", oin)
        self.assertIn("[N+]", oin)
        self.assertIn("[O-]", oin)

    def test_abuvup_fragments_all_reparse(self):
        # The property oin/compare.py needs: an unparseable fragment degrades to a
        # RAW: token and the round-trip key can never match.
        oin = XYZToSMILES().convert(os.path.join(_FIXTURES, "ABUVUP.xyz"))
        for fragment in _ligand_fragments(oin):
            self.assertIsNotNone(Chem.MolFromSmiles(fragment), f"unparseable: {fragment}")

    def test_bound_carbonyl_is_not_rewritten(self):
        # [C-]#[O+]'s oxygen is hypervalent when neutral, but its balancing
        # carbanion is the Zone A donor. Restoring the [O+] alone would turn the
        # well-formed `C{0}#O` into `C{0}#[O+]`. Charges are restored in pairs.
        oin = XYZToSMILES().convert(os.path.join(_FIXTURES, "FeCO5.xyz"))
        self.assertNotIn("#[O+]", oin)
        self.assertIn("#O", oin)


if __name__ == "__main__":
    unittest.main()
