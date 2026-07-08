"""Regression: eta3-allyl (and any metal-bound, under-valent ligand) bond-order
transfer in ``build_contract_mol``.

An OIN ligand atom that binds the metal is under-valent once the metal is
stripped, so RDKit gives it a radical electron -- e.g. the three metal-bound
carbons of an eta3-allyl, ``[CH2][CH]=[CH]...`` (and after the template is
flattened to all-single bonds for connectivity matching, every one of them is a
radical). ``build_contract_mol`` transfers bond orders/aromaticity by a
connectivity-only substructure match of ``_flatten_template(template)`` into the
MetalloGen-generated fragment. RDKit's substructure matcher treats
radical-electron count as a match constraint, so a radical-bearing query
silently fails to match the H-saturated generated fragment: no bond orders or
aromaticity transfer, and the ligand is emitted all-single / dearomatized on the
round trip (the eta3-allyl "double-bond loss" failure, e.g. ABAZEK).

``_flatten_template`` must therefore normalize the query to a plain connectivity
graph (0 radicals, implicit H) so the match succeeds. These tests pin that.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem

from oinsmiles.generation.metallogen_adapter import _flatten_template, _oin_fragment_templates
from oinsmiles.generation.oin_parser import OINParser

# ABAZEK: Pd eta3-(1-phenylallyl) + aminophenyl-phosphine + Cl. The allyl-phenyl
# ligand is the one whose bond-order transfer used to fail.
_ABAZEK_OIN = (
    "[Pd_TPL].[CH2]{0>}[CH]{0}=[CH]{0}c1ccccc1.CN(C)c1ccc(P{1}(C(C)(C)C)C(C)(C)C)cc1.[Cl]{2}"
)

# The MetalloGen-generated allyl-phenyl fragment is H-saturated with all-single
# bonds: a 3-carbon chain on a 6-membered carbocycle (propyl-on-ring skeleton).
_SATURATED_ALLYL_PHENYL = "CCCC1CCCCC1"


def _allyl_phenyl_template():
    """Build the allyl-phenyl heavy-atom template the way _oin_fragment_templates does."""
    t = Chem.MolFromSmiles("[CH2][CH]=[CH]c1ccccc1", sanitize=False)
    Chem.SanitizeMol(t)
    return Chem.RemoveHs(t, sanitize=False)


class TestAllylTemplateFlatten(unittest.TestCase):
    def test_precondition_template_has_a_radical(self):
        """Documents WHY the fix is needed: the metal-bound template is a radical."""
        t = _allyl_phenyl_template()
        self.assertTrue(
            any(a.GetNumRadicalElectrons() for a in t.GetAtoms()),
            "eta3-allyl template should carry a radical (metal-bound, under-valent)",
        )

    def test_flatten_clears_radicals(self):
        fl = _flatten_template(_allyl_phenyl_template())
        self.assertTrue(
            all(a.GetNumRadicalElectrons() == 0 for a in fl.GetAtoms()),
            "flattened connectivity query must have no radical electrons",
        )

    def test_flatten_matches_saturated_generated_fragment(self):
        """The core guard: the flattened template must substructure-match the
        H-saturated generated fragment (this is the match build_contract_mol runs;
        before the fix it returned empty and no bond orders transferred)."""
        fl = _flatten_template(_allyl_phenyl_template())
        gen = Chem.MolFromSmiles(_SATURATED_ALLYL_PHENYL)
        self.assertTrue(
            gen.GetSubstructMatch(fl),
            "radical-bearing eta3-allyl template must match the saturated fragment",
        )

    def test_parsed_oin_allyl_template_matches(self):
        """End-to-end template path: parse a real OIN, and its 9-heavy-atom
        allyl-phenyl template must match the saturated skeleton."""
        parsed = OINParser().parse(_ABAZEK_OIN)
        flats = [_flatten_template(t) for t in _oin_fragment_templates(parsed)]
        gen = Chem.MolFromSmiles(_SATURATED_ALLYL_PHENYL)
        allyl = [f for f in flats if f.GetNumAtoms() == gen.GetNumAtoms()]
        self.assertTrue(allyl, "expected a 9-heavy-atom allyl-phenyl template")
        self.assertTrue(
            any(gen.GetSubstructMatch(f) for f in allyl),
            "allyl-phenyl template from the parsed OIN must match the saturated fragment",
        )


if __name__ == "__main__":
    unittest.main()
