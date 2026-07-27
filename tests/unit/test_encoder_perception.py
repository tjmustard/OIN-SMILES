"""Unit tests for bond-order / charge perception crashes in the XYZ->OIN encoder.

Three round-trip failure classes traced to ``perception_tmc``'s input-side perception:

* ``xyz2mol_none_crash`` -- ``lig_checks`` called ``len(res_mols)`` on a
  ``ResonanceMolSupplier`` and then iterated it. The supplier is a stateful
  iterator: ``len()`` runs the enumeration and leaves the cursor at the end, so
  the loop yielded ``None`` and ``res_mol.GetAtoms()`` raised ``AttributeError``.
  Indexing the supplier never returns ``None`` for the same ligand.

* ``garbled_aromatic`` (the nitro sub-bucket) -- ``get_oin_string`` zeroed every
  formal charge, so a nitro group's ``[N+](=O)[O-]`` became a four-bonded neutral
  nitrogen and serialized as the unparseable ``N(O)=O``.

* ``kekulize_encode_crash`` -- ``AC2mol`` draws a double bond from phosphorus onto
  an aromatic ring carbon (a PPN+ counter-cation, a phosphonium ylide). The ring
  is unkekulizable and the first full ``SanitizeMol`` raised. All 17 tracebacks
  pointed at ``fix_equivalent_Os``, which rewrites nothing on these molecules --
  it is simply whichever sanitize ran first. The underlying cause is a wrong
  ligand charge from extended Huckel, so an unusable perception is re-perceived
  over a charge sweep.
"""

import os
import re
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem, RDLogger

from oinsmiles.core.translator import XYZToSMILES
from oinsmiles.utils.aromaticity import (
    OINEncodeError,
    kekulize_safe_sanitize,
    stuck_ring_atoms,
)
from oinsmiles.utils.perception_tmc import _perception_is_usable, lig_checks

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


class TestStuckRingDetection(unittest.TestCase):
    def test_detects_exocyclic_double_bond_on_aromatic_ring(self):
        # A quinoid ring: aromatic flags kept, exocyclic C=N hanging off the ring.
        mol = Chem.MolFromSmiles("Cc1cccc(=[N]c2ccccc2)[n]1", sanitize=False)
        self.assertIsNotNone(mol)
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        self.assertTrue(stuck_ring_atoms(mol))

    def test_plain_benzene_is_not_stuck(self):
        self.assertEqual(stuck_ring_atoms(Chem.MolFromSmiles("Cc1ccccc1")), set())

    def test_pyridone_exocyclic_oxygen_is_reported(self):
        # Pyridone genuinely has an exocyclic C=O on an aromatic ring; RDKit cannot
        # kekulize it either, so reporting it as stuck is correct, not a false hit.
        mol = Chem.MolFromSmiles("O=c1cccc[nH]1")
        self.assertTrue(stuck_ring_atoms(mol))


class TestKekulizeSafeSanitize(unittest.TestCase):
    def test_noop_on_a_clean_molecule(self):
        mol = Chem.MolFromSmiles("Cc1ccccc1")
        self.assertIs(kekulize_safe_sanitize(mol), mol)

    def test_quinoid_ylide_input_encodes(self):
        # NAXDOI: a [P+]=c phosphonium ylide. Before the fix, XYZToSMILES raised
        # "Can't kekulize mol. Unkekulized atoms: 14 26 34 42 46".
        oin = XYZToSMILES().convert(os.path.join(_FIXTURES, "NAXDOI.xyz"))
        self.assertTrue(oin.startswith("[Mn_TET]"), oin)

    def test_oin_encode_error_is_a_value_error(self):
        # core.translator catches ValueError; the new type must not slip past it.
        self.assertTrue(issubclass(OINEncodeError, ValueError))


class TestChargeRescue(unittest.TestCase):
    """A ligand whose perceived bond orders are unusable is re-perceived."""

    def test_ppn_cation_is_perceived_correctly(self):
        # AGUFEN carries a PPN+ counter-cation, [Ph3P=N=PPh3]+. Extended Huckel
        # proposes -1; AC2mol succeeds there and at -3, both times by drawing double
        # bonds from phosphorus onto phenyl carbons -- leaving the ipso carbon
        # pentavalent, which no ring relaxation can rescue. A charge sweep finds the
        # chemically correct q=+1, where the ligand is stuck-ring- and radical-free.
        oin = XYZToSMILES().convert(os.path.join(_FIXTURES, "AGUFEN.xyz"))
        self.assertTrue(oin.startswith("[Y_TPL]"), oin)
        # The P-N-P backbone, with intact phenyls rather than P=c quinoid rings.
        self.assertIn("P(N", oin)
        self.assertIn("c1ccc(P", oin)

    def test_rescue_does_not_fire_on_a_usable_perception(self):
        # The rescue is additive: a ligand whose perception already sanitizes never
        # enters it, which is what keeps the 1944 passing molecules still.
        mol = Chem.MolFromSmiles("c1ccc(P(c2ccccc2)c2ccccc2)cc1")
        self.assertTrue(_perception_is_usable(mol))

    def test_genuine_quinoid_ligand_is_not_rescued_away(self):
        # A 2-iminopyridine really is quinoid. Its stuck ring is chemically correct,
        # kekulize_safe_sanitize relaxes it, and so it must stay usable -- otherwise
        # the rescue would re-perceive ABERIK's ligand at some other charge.
        quinoid = Chem.MolFromSmiles(
            "Cc1cc(C)c(-c2cccc(=[N]c3c(C)cc(C)cc3C)[n]2)c(C)c1", sanitize=False
        )
        self.assertIsNotNone(quinoid)
        quinoid.UpdatePropertyCache(strict=False)
        self.assertTrue(_perception_is_usable(quinoid))


if __name__ == "__main__":
    unittest.main()
