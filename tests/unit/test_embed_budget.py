"""Guards for R2 (v0.3.7 no-conformer wave): the generation-internal wall-clock
budget and the double-bond-stereo over-valence guard.

Both target the ``no_conformers`` failure class:

* ``embed_time_budget`` bounds the FF-only attempt loop. Before it, ``timeout`` was
  consumed only by the ASE optimizer, so a molecule whose embed never validated ran
  the full ``max_attempts`` budget (ZIHGEE_comp_0: ~1696 s) before returning nothing.

* ``_apply_double_bond_stereo`` used to force a carried stereo bond back to DOUBLE
  unconditionally. When PuLP had relocated the double bond, that made an endpoint
  over-valent (FIXYER_comp_0: a 5-valent carbon), so every downstream
  ``SanitizeMol`` / ``MolToSmiles`` raised and generation produced no conformer. The
  guard promotes only when it keeps the molecule valence-valid.

Each test fails against the pre-fix code (the budget param did not exist; the
promotion was unconditional).
"""

import os
import sys
import time
import unittest
from unittest import mock

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem

from oinsmiles.generation.metallogen_adapter import convert_parsed_to_msmiles
from oinsmiles.generation.oin_parser import OINParser
from oinsmiles.generator3d import embed, generate_3d_structures

# A trivial, always-parseable complex; the loop body is monkeypatched, so the only
# requirement is that om.get_om_from_modified_smiles() succeeds on it.
CISPLATIN_OIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"


def _msmiles(oin):
    return convert_parsed_to_msmiles(OINParser().parse(oin))


class TestEmbedTimeBudget(unittest.TestCase):
    def setUp(self):
        self.msmiles = _msmiles(CISPLATIN_OIN)

    def test_budget_stops_before_max_attempts(self):
        """A tiny wall-clock budget cuts the attempt loop far short of max_attempts."""
        calls = {"n": 0}

        def slow_fail(*a, **k):
            calls["n"] += 1
            time.sleep(0.05)
            return None

        with mock.patch.object(embed, "get_embedding", slow_fail):
            res = generate_3d_structures(
                self.msmiles, ff_params={"max_attempts": 200}, embed_time_budget=0.3
            )
        self.assertEqual(res, [], "no embed validated, so the pool must be empty")
        self.assertLess(
            calls["n"], 60, "the 0.3 s budget must stop the loop well short of 200 attempts"
        )

    def test_none_budget_runs_full_attempt_budget(self):
        """embed_time_budget=None preserves the prior unbounded behavior."""
        calls = {"n": 0}

        def fast_fail(*a, **k):
            calls["n"] += 1
            return None

        with mock.patch.object(embed, "get_embedding", fast_fail):
            res = generate_3d_structures(
                self.msmiles, ff_params={"max_attempts": 15}, embed_time_budget=None
            )
        self.assertEqual(res, [])
        self.assertEqual(calls["n"], 15, "an unbounded budget must run every attempt")


class TestDoubleBondStereoValenceGuard(unittest.TestCase):
    """``_apply_double_bond_stereo`` must not promote a bond that over-valences an atom."""

    def test_overvalent_promotion_is_skipped(self):
        # Neopentane: the central carbon (idx 1) already has four single bonds.
        # Promoting bond 0-1 to DOUBLE would make it 5-valent -- the FIXYER failure.
        mol = Chem.MolFromSmiles("CC(C)(C)C", sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        stereo_bonds = [(0, 1, Chem.BondStereo.STEREOE, 2, 3)]
        embed._apply_double_bond_stereo(mol, stereo_bonds)  # must not raise
        bond = mol.GetBondBetweenAtoms(0, 1)
        self.assertEqual(
            bond.GetBondType(),
            Chem.BondType.SINGLE,
            "guard must leave an over-valent promotion as a single bond",
        )

    def test_valid_promotion_still_happens(self):
        # Butane: promoting the central bond 1-2 to DOUBLE yields but-2-ene --
        # valence-valid, so the E/Z enforcement the S6 pendant-alkene case relies on
        # must still apply. Refs 0/3 are genuine neighbours, so the stereo set is clean.
        mol = Chem.MolFromSmiles("CCCC", sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        stereo_bonds = [(1, 2, Chem.BondStereo.STEREOZ, 0, 3)]
        embed._apply_double_bond_stereo(mol, stereo_bonds)
        bond = mol.GetBondBetweenAtoms(1, 2)
        self.assertEqual(
            bond.GetBondType(),
            Chem.BondType.DOUBLE,
            "a valence-valid promotion must still restore the double bond",
        )
        self.assertEqual(bond.GetStereo(), Chem.BondStereo.STEREOZ)

    def test_promotion_keeps_valence_helper(self):
        mol = Chem.MolFromSmiles("CC(C)(C)C", sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        mol.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.DOUBLE)
        self.assertFalse(embed._promotion_keeps_valence(mol), "5-valent carbon must be rejected")
        ok = Chem.MolFromSmiles("CCC", sanitize=False)
        ok.UpdatePropertyCache(strict=False)
        ok.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.DOUBLE)
        self.assertTrue(embed._promotion_keeps_valence(ok), "propene is valence-valid")


if __name__ == "__main__":
    unittest.main()
