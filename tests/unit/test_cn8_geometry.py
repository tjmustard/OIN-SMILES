"""Unit tests for CN-8 (square-antiprismatic, ``SQA``) geometry support.

Coordination number 8 (e.g. Hf/Zr with four bidentate salicylaldimine ligands)
previously had no encoder template, so the aligner emitted the undetermined
``NON`` code and the generator raised
``Geometry code 'NON' not supported by MetalloGen mapping``. This adds ``SQA``
to both the encoder matcher (``_find_best_geometry_match`` / ``TEMPLATES``) and
the generator mapping (``OIN_TO_METALLOGEN_GEO`` -> ``8_squre_antiprismatic``),
so a CN-8 OIN now round-trips through a concrete geometry with slot markers.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.metallogen_adapter import (
    OIN_TO_METALLOGEN_GEO,
    convert_parsed_to_msmiles,
)
from oinsmiles.generation.oin_parser import TEMPLATES as PARSER_TEMPLATES
from oinsmiles.generation.oin_parser import OINParser
from oinsmiles.utils.oin_aligner import (
    TEMPLATES,
    classify_coordination_geometry,
    coordination_geometry_fit,
)

# Ideal square-antiprismatic donor directions (metal at origin), the same two
# staggered squares MetalloGen's ``8_squre_antiprismatic`` uses, left
# un-normalized to confirm the matcher is scale-invariant.
_SQA = [
    [-1, 1, 1],
    [1, 1, 1],
    [-1, -1, 1],
    [1, -1, 1],
    [-1.41, 0, -1],
    [0, -1.41, -1],
    [1.41, 0, -1],
    [0, 1.41, -1],
]

# A real CN-8 OIN (AFEPIM: Hf with four salicylaldimine N,O bidentate ligands),
# as emitted by the encoder after the SQA extension. Used as a fixture so the
# parse -> m-SMILES path is exercised without the heavy XYZ encoder / scine stack.
_AFEPIM_OIN = (
    "[Hf_SQA]."
    "Cc1cccc(C)c1N{0}=C(O{1})c1ccccc1."
    "Cc1cccc(C)c1N{3}=C(O{2})c1ccccc1."
    "Cc1cccc(C)c1N{4}=C(O{5})c1ccccc1."
    "Cc1cccc(C)c1N{6}=C(O{7})c1ccccc1"
)


class TestSqaTemplate(unittest.TestCase):
    def test_encoder_template_has_eight_slots(self):
        self.assertIn("SQA", TEMPLATES)
        self.assertEqual(len(TEMPLATES["SQA"]), 8)

    def test_parser_template_has_eight_slots(self):
        # Second (duplicated) TEMPLATES dict lives in oin_parser; both must know SQA
        # or the parser drops every CN-8 slot vector.
        self.assertIn("SQA", PARSER_TEMPLATES)
        self.assertEqual(len(PARSER_TEMPLATES["SQA"]), 8)

    def test_generator_mapping(self):
        self.assertEqual(OIN_TO_METALLOGEN_GEO.get("SQA"), "8_squre_antiprismatic")


class TestClassifySqa(unittest.TestCase):
    def test_square_antiprismatic(self):
        self.assertEqual(classify_coordination_geometry(_SQA), "SQA")

    def test_ideal_fit_is_small(self):
        # Truncated 1.41 (vs sqrt2) makes it a near-, not exactly-ideal antiprism.
        self.assertLess(coordination_geometry_fit(_SQA, "SQA"), 1e-2)


class TestSqaMsmilesRoundTrip(unittest.TestCase):
    def test_parse_emits_eight_vectors(self):
        parsed = OINParser().parse(_AFEPIM_OIN)
        self.assertEqual(parsed.geo_code, "SQA")
        self.assertEqual(len(parsed.vectors), 8)

    def test_convert_to_msmiles_maps_to_square_antiprism(self):
        parsed = OINParser().parse(_AFEPIM_OIN)
        # This raised "Geometry code 'NON' not supported" before the fix.
        msmiles = convert_parsed_to_msmiles(parsed)
        self.assertIn("8_squre_antiprismatic", msmiles)
        # Eight binding atoms -> eight 1-based MetalloGen map numbers :1..:8.
        for n in range(1, 9):
            self.assertIn(f":{n}]", msmiles)


if __name__ == "__main__":
    unittest.main()
