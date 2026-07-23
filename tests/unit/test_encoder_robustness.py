"""SL5 (v0.4.4) encoder-robustness guards for the XYZ->OIN encode-fail cohort.

The v0.4.2 capstone had 48 molecules whose crystal XYZ produced no OIN string at all
(``smiles_1 is None``). Sub-triage (``tools/sl5_triage.py``) split them into:

* ~34 electron-deficient boron clusters (carboranes / closo-nido boranes): RDKit's
  2-center-2-electron valence model has no Lewis structure for a 3-center-2-electron cage,
  so these are an irreducible ceiling. The encoder now *classifies* the failure with a
  typed ``OINEncodeError`` instead of a bare ``ValueError`` (W1).
* A timeout cohort that hung inside xyz2mol perception on large conjugated ligands --
  ``AC2BO`` materialising an exponential valence-order product, and ``ResonanceMolSupplier``
  building conjugation groups. Both are now bounded so oversized ligands complete (encode
  or fail fast) instead of hanging; small/medium ligands are byte-identical (W3).

Each test fails against pre-SL5 code (which raised ``ValueError`` untyped, or hung).
Recovery fixtures are distilled from the tmCAT/tmPHOTO sweep. Assertions stay
rdkit-version robust (no exact bond-direction strings): a non-empty encode led by the
correct metal token, plus forward-encode stability.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem, RDLogger

from oinsmiles.core.translator import XYZToSMILES
from oinsmiles.utils.aromaticity import OINEncodeError
from oinsmiles.utils.xyz2mol import (
    _is_electron_deficient_cluster,
    get_tmc_mol,
)

RDLogger.DisableLog("rdApp.*")

_FIXTURES = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))


def _fixture(name):
    return os.path.join(_FIXTURES, name)


class TestBoronClusterTypedCeiling(unittest.TestCase):
    """W1: an electron-deficient boron cluster fails with a classified OINEncodeError."""

    FIXTURE = _fixture("RAWJEG_comp_0.xyz")

    @unittest.skipUnless(os.path.exists(FIXTURE), "fixture missing")
    def test_get_tmc_mol_raises_typed_error(self):
        with self.assertRaises(OINEncodeError) as ctx:
            get_tmc_mol(self.FIXTURE, 0, with_stereo=False)
        self.assertIn("boron cluster", str(ctx.exception))

    @unittest.skipUnless(os.path.exists(FIXTURE), "fixture missing")
    def test_convert_propagates_typed_error(self):
        # convert() re-raises the typed limitation instead of flattening it to a bare
        # ValueError, so callers can tell a known ceiling from an unexpected failure.
        with self.assertRaises(OINEncodeError):
            XYZToSMILES().convert(self.FIXTURE)

    def test_typed_error_is_a_value_error(self):
        # Back-compat: existing ``except ValueError`` handlers still catch it.
        self.assertTrue(issubclass(OINEncodeError, ValueError))


class TestElectronDeficientClusterDetector(unittest.TestCase):
    """W1: the detector fires on cage boron (B-B bonds) but not on borate/organic."""

    @staticmethod
    def _mol_from_bonds(elements, bonds):
        rw = Chem.RWMol()
        for z in elements:
            rw.AddAtom(Chem.Atom(z))
        for i, j in bonds:
            rw.AddBond(i, j, Chem.BondType.SINGLE)
        return rw.GetMol()

    def test_cage_boron_is_detected(self):
        # 3 borons in a B-B-B chain: >=3 B and a B-B bond.
        mol = self._mol_from_bonds([5, 5, 5], [(0, 1), (1, 2)])
        self.assertTrue(_is_electron_deficient_cluster(mol))

    def test_bph4_borate_is_not_detected(self):
        # One boron, four B-C bonds, no B-B bond (BPh4- perceives fine).
        mol = self._mol_from_bonds([5, 6, 6, 6, 6], [(0, 1), (0, 2), (0, 3), (0, 4)])
        self.assertFalse(_is_electron_deficient_cluster(mol))

    def test_scattered_borons_without_bb_bond_not_detected(self):
        # 3 borons but each only bonds carbon -- not a cluster.
        mol = self._mol_from_bonds([5, 6, 5, 6, 5], [(0, 1), (2, 3), (4, 1)])
        self.assertFalse(_is_electron_deficient_cluster(mol))

    def test_organic_ligand_is_not_detected(self):
        self.assertFalse(_is_electron_deficient_cluster(Chem.MolFromSmiles("c1ccncc1")))


class TestAc2boCapIsByteIdentical(unittest.TestCase):
    """W3: the AC2BO valence-combo cap is inert on ordinary (sub-cap) ligands.

    The cap only bypasses the valence-ordering sort when the per-atom valence product
    exceeds _VALENCE_COMBO_CAP; every ordinary ligand stays under it and takes the
    original sorted path, so encodes are byte-identical. (Corpus/golden byte-identity is
    covered by test_regression_stability; this pins the constant is present and sane so a
    future edit cannot silently lower it into the range that would perturb real ligands.)
    """

    def test_cap_is_large(self):
        from oinsmiles.utils.xyz2mol_local import _VALENCE_COMBO_CAP

        # A cisplatin/ferrocene-scale ligand's valence product is a few dozen at most;
        # the cap must stay far above that so ordinary ligands never take the fallback.
        self.assertGreaterEqual(_VALENCE_COMBO_CAP, 100_000)

    def test_ordered_valences_matches_unsorted_content(self):
        # The extracted sorter is a permutation of the raw product -- same members.
        import itertools

        from oinsmiles.utils.xyz2mol_local import _ordered_valences

        vlls = [[4], [3, 4], [2, 1]]
        atoms = [6, 7, 8]
        got = sorted(map(tuple, _ordered_valences(vlls, atoms)))
        expect = sorted(itertools.product(*vlls))
        self.assertEqual(got, expect)


class TestForkedResonanceRecovery(unittest.TestCase):
    """W3: a molecule whose ResonanceMolSupplier hangs is recovered by the forked, time-
    bounded resonance -- the child is killed at the budget and perception degrades to the
    single form. A completing large ligand instead stays byte-identical (covered by the
    branch-vs-main firing-set gate and test_regression_stability)."""

    FIXTURE = _fixture("BENVOG_comp_0.xyz")

    @unittest.skipUnless(os.path.exists(FIXTURE) and hasattr(os, "fork"), "fixture/fork missing")
    def test_benvog_recovers_via_cpu_budget_fallback(self):
        from oinsmiles.utils import xyz2mol

        # BENVOG's macrocycle burns CPU in ResonanceMolSupplier without ever finishing;
        # shrink the CPU budget so the kernel (SIGXCPU) kills the child quickly and perception
        # falls back to the single form.
        orig = xyz2mol._RESONANCE_CPU_BUDGET_S
        xyz2mol._RESONANCE_CPU_BUDGET_S = 2
        try:
            oin = XYZToSMILES().convert(self.FIXTURE)
            self.assertTrue(oin, "expected BENVOG to encode via the resonance-budget fallback")
            self.assertIn("[Ni", oin)
            # Deterministic across a repeat (same fallback form).
            self.assertEqual(oin, XYZToSMILES().convert(self.FIXTURE))
        finally:
            xyz2mol._RESONANCE_CPU_BUDGET_S = orig


if __name__ == "__main__":
    unittest.main()
