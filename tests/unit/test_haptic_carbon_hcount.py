"""Regression tests for haptic (eta) ring-carbon hydrogen handling in the
MetalloGen m-SMILES.

`convert_parsed_to_msmiles` keeps the implicit hydrogen of a haptic ring carbon
(a Cp/arene carbon that coordinates the metal as part of an eta face). For a
GENUINE C-H haptic carbon that is correct -- it re-materializes as 1 H. But a
BARE 0-H haptic carbon (an ipso / ring-fusion / substituted carbon written
`c{n}` with no explicit H) is over-protonated by the anionic-aromatic re-parse
in `get_ligand_from_smiles`: `AddHs(explicitOnly=False)` puts a phantom H on a
carbon that already has three heavy neighbours. The full-dataset atom-count
check surfaced this as +H round-trip failures on fused-ring eta ligands
(ARONEA +4, BOXJUU +6). The fix locks a bare 0-H haptic carbon at 0 H
(NoImplicit) so it survives the re-parse without gaining a hydrogen, while a
genuine C-H haptic carbon (its `[cH]` bracket already sets NoImplicit) is left
untouched -- so Cp/arene passers stay byte-identical.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem

from oinsmiles.generation.metallogen_adapter import convert_oin_to_msmiles


def _eta_ligand_frag(msmiles: str) -> str:
    """Return the ligand fragment bound at slot 1 (the eta face under test)."""
    return next(p for p in msmiles.split("|") if ":1]" in p)


def _addhs_h_count(frag: str) -> int:
    mol = Chem.MolFromSmiles(frag, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    return sum(1 for a in Chem.AddHs(mol, explicitOnly=False).GetAtoms() if a.GetAtomicNum() == 1)


class TestHapticCarbonHCount(unittest.TestCase):
    def test_fused_eta_ring_fusion_carbons_stay_zero_h(self):
        """An eta indenyl's ring-fusion / ipso carbons must not gain a phantom H.

        Pre-fix the anionic-aromatic re-parse protonated all five 5-ring carbons
        (eta ligand -> 9 H); the two fused and one ipso carbons are 0-H, so the
        correct count is 6.
        """
        oin = "[Fe_LIN].c{0}1[cH]{0}c{0}2ccccc{0}2[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
        frag = _eta_ligand_frag(convert_oin_to_msmiles(oin))
        self.assertEqual(_addhs_h_count(frag), 6, f"phantom H on eta fusion carbon in: {frag}")

        # structural guard: no aromatic carbon with 3 heavy neighbours carries an H
        mol = Chem.MolFromSmiles(frag, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        for a in mol.GetAtoms():
            heavy = sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1)
            if a.GetSymbol() == "C" and a.GetIsAromatic() and heavy >= 3:
                self.assertEqual(
                    a.GetTotalNumHs(), 0, f"deg-{heavy} aromatic C over-protonated in {frag}"
                )

    def test_genuine_ch_haptic_carbons_keep_their_h(self):
        """A ferrocene Cp ring (all `[cH]`) keeps 1 H per carbon -- the fix must
        not strip a genuine haptic C-H (negative control / byte-identity)."""
        oin = (
            "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
        )
        frag = _eta_ligand_frag(convert_oin_to_msmiles(oin))
        self.assertEqual(_addhs_h_count(frag), 5, f"ferrocene Cp lost/gained H in: {frag}")

    def test_substituted_eta_ipso_is_not_protonated(self):
        """A methyl-Cp's substituted ipso carbon (`Cc{0}`) stays 0-H; the four
        ring CH keep their H (methyl 3 + 4 ring CH = 7)."""
        oin = "[Fe_LIN].Cc{0}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
        frag = _eta_ligand_frag(convert_oin_to_msmiles(oin))
        self.assertEqual(_addhs_h_count(frag), 7, f"methyl-Cp ipso mis-protonated in: {frag}")

    def test_sigma_aryl_carbanion_still_strips_its_h(self):
        """A non-haptic sigma-aryl carbon (single coordinating atom in the ring)
        is still stripped to 0 H -- the fix must not disturb the sigma branch."""
        oin = "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1"
        msmiles = convert_oin_to_msmiles(oin)
        # the cyclometalated aryl carbon binds as a bare [c:n] (0 H), not [cH:n]
        self.assertNotIn("[cH:", msmiles, f"sigma-aryl carbon protonated in: {msmiles}")


if __name__ == "__main__":
    unittest.main()
