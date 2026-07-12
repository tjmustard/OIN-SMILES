"""Stereo Phase 4 (Zone-A P Stereocenter Encoding) — parser-side passthrough test.

The generation-side enforcement tests in this module previously exercised the
legacy Molassembler/stitch backend (``OIN3DGenerator(engine="legacy")`` plus its
``_stitch_fragment`` / ``_verify_zone_a_p`` DG-fallback machinery). That backend
was removed, and MetalloGen's own Zone-A P / atom-stereo behaviour is covered by
the round-trip diagnostics and the ``verify_xyz_to_oin`` integration checks. The
one engine-agnostic assertion — that bracketed P chiral tags survive parsing into
``ParsedOIN`` — is retained below.
"""

import os
import sys
import unittest

from rdkit import Chem

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.oin_parser import OINParser


class TestZoneAPParsedOINPassthrough(unittest.TestCase):
    """[P@]{0} / [P@@]{1>} survive intact at the ParsedOIN level -- asserts the
    tag-preservation property at OINParser.parse()'s output, which is what the
    generation backend consumes.
    """

    def test_bracket_p_chiral_tags_survive_to_parsed_oin(self):
        inline = "[Rh_SPL].c1ccccc1[P@]{0}(C)c1ccccc1.[P@@]{1}.[Cl]{2}.[Cl]{3}"
        parsed = OINParser().parse(inline)

        self.assertIn("[P@]", parsed.fragments[1])
        self.assertEqual(parsed.fragments[2], "[P@@]")

        # atom_in_fragment_idx for the tagged P in fragment 1 must point at
        # the P atom itself (index 6 in "c1ccccc1[P@](C)c1ccccc1": 6 ring
        # atoms at indices 0-5, P at index 6).
        p_vector = next(v for v in parsed.vectors if v.fragment_idx == 1)
        self.assertEqual(p_vector.atom_in_fragment_idx, 6)
        frag_mol = Chem.MolFromSmiles(parsed.fragments[1], sanitize=False)
        self.assertEqual(frag_mol.GetAtomWithIdx(6).GetAtomicNum(), 15)
        self.assertNotEqual(
            frag_mol.GetAtomWithIdx(6).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED
        )


if __name__ == "__main__":
    unittest.main()
