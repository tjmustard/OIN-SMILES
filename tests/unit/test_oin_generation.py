import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.oin_parser import OINParser


class TestOINParser(unittest.TestCase):
    def test_parse_simple(self):
        oin = "[Pt].[Cl] |g:SPL|w:1.0:0|"
        parser = OINParser()
        parsed = parser.parse(oin)
        self.assertEqual(parsed.smiles, "[Pt].[Cl]")
        self.assertEqual(len(parsed.vectors), 1)
        self.assertEqual(
            parsed.vectors[0].atom_idx, -1
        )  # Parser sets this to -1 for vector/template based
        self.assertEqual(parsed.vectors[0].fragment_idx, 1)  # [Pt] is 0, [Cl] is 1
        self.assertEqual(parsed.vectors[0].atom_in_fragment_idx, 0)

    def test_parse_complex(self):
        oin = "[Pt].[NH2]CC[NH2] |g:SPL|w:1.0:0;1.3:1|"
        # [Pt] -> 0
        # [NH2]CC[NH2] -> 1. Atoms: N(0), C(1), C(2), N(3). Global: 1, 2, 3, 4.
        parser = OINParser()
        parsed = parser.parse(oin)

        self.assertEqual(len(parsed.vectors), 2)

        v1 = parsed.vectors[0]
        self.assertEqual(v1.fragment_idx, 1)
        self.assertEqual(v1.atom_in_fragment_idx, 0)

        v2 = parsed.vectors[1]
        self.assertEqual(v2.fragment_idx, 1)
        self.assertEqual(v2.atom_in_fragment_idx, 3)

    def test_parse_inline_winding_threads_to_template_vectors(self):
        # Ferrocene: LIN geometry (template present), heading atom on each
        # ring carries '>', all other ring atoms carry None.
        oin = (
            "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
        )
        parser = OINParser()
        parsed = parser.parse(oin)

        windings = [v.winding for v in parsed.vectors]
        self.assertEqual(windings.count(">"), 2)  # one heading atom per ring
        self.assertTrue(all(w in (">", None) for w in windings))

    def test_parse_inline_winding_survives_template_less_geometry(self):
        # NON geometry has no TEMPLATES entry, so no OINVector is emitted,
        # but winding_by_slot must still capture the winding (Red Team gap).
        oin = "[Fe_NON].[Cl]{0>}"
        parser = OINParser()
        parsed = parser.parse(oin)

        self.assertEqual(parsed.vectors, [])
        self.assertEqual(parsed.winding_by_slot[0], ">")

    def test_winding_by_slot_survives_multi_atom_eta_slots(self):
        # Multi-atom eta slots: ferrocene rings with heading atoms carrying
        # winding ('>' or '<') and trailing atoms carrying None. Non-heading atoms
        # must not clobber the heading atom's winding.
        oin = (
            "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1<}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
        )
        parser = OINParser()
        parsed = parser.parse(oin)

        self.assertEqual(parsed.winding_by_slot, {0: ">", 1: "<"})


if __name__ == "__main__":
    unittest.main()
