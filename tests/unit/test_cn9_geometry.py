"""Unit tests for CN-9 (tricapped trigonal prismatic, ``TCT``) geometry support.

Coordination number 9 (e.g. Y/lanthanide with three bidentate + one tridentate
donor, XERTUK_comp_3) previously had no encoder template, so the aligner emitted
the undetermined ``NON`` code and the generator raised
``Geometry code 'NON' not supported by MetalloGen mapping``. This adds ``TCT``
to both the encoder matcher (``_find_best_geometry_match`` / ``TEMPLATES``) and
the generator mapping (``OIN_TO_METALLOGEN_GEO`` -> ``9_tricapped_trigonal_prismatic``),
so a CN-9 OIN now round-trips through a concrete geometry with slot markers.
Mirrors ``test_cn8_geometry.py``.
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

# Ideal tricapped-trigonal-prismatic donor directions (metal at origin), the same
# unit-normalized vectors MetalloGen's ``9_tricapped_trigonal_prismatic`` uses.
# Left un-normalized here (scaled by 2) to confirm the matcher is scale-invariant.
_TCT = [
    [1.5111472, 0.0000000, 1.3101276],
    [-0.7561420, 1.3092458, 1.3092458],
    [-0.7561420, -1.3092458, 1.3092458],
    [1.5111472, 0.0000000, -1.3101276],
    [-0.7561420, 1.3092458, -1.3092458],
    [-0.7561420, -1.3092458, -1.3092458],
    [-2.0000000, 0.0000000, 0.0000000],
    [1.1643746, -1.6261094, 0.0000000],
    [1.1643746, 1.6261094, 0.0000000],
]

# A real CN-9 OIN (XERTUK_comp_3: Y with three bidentate nitrates + one tridentate
# N,N,N ligand), as emitted by the encoder after the TCT extension. Used as a
# fixture so the parse -> m-SMILES path is exercised without the heavy XYZ encoder
# stack.
_XERTUK_OIN = (
    "[Y_TCT]."
    "O{0}N(O{7})=O."
    "CCCCN(CCCC)c1ccc(-c2c3c(n{4}c4c2CCc2c(C)c5ccccc5n{5}c42)"
    "c2n{1}c4ccccc4c(C)c2CC3)cc1."
    "O{6}N(O{2})=O."
    "O{8}N(O{3})=O"
)


class TestTctTemplate(unittest.TestCase):
    def test_encoder_template_has_nine_slots(self):
        self.assertIn("TCT", TEMPLATES)
        self.assertEqual(len(TEMPLATES["TCT"]), 9)

    def test_parser_template_has_nine_slots(self):
        # Second (duplicated) TEMPLATES dict lives in oin_parser; both must know TCT
        # or the parser drops every CN-9 slot vector (-> UncoordinatedFragmentError).
        self.assertIn("TCT", PARSER_TEMPLATES)
        self.assertEqual(len(PARSER_TEMPLATES["TCT"]), 9)

    def test_generator_mapping(self):
        self.assertEqual(OIN_TO_METALLOGEN_GEO.get("TCT"), "9_tricapped_trigonal_prismatic")


class TestClassifyTct(unittest.TestCase):
    def test_tricapped_trigonal_prismatic(self):
        self.assertEqual(classify_coordination_geometry(_TCT), "TCT")

    def test_ideal_fit_is_small(self):
        self.assertLess(coordination_geometry_fit(_TCT, "TCT"), 1e-2)


class TestTctMsmilesRoundTrip(unittest.TestCase):
    def test_parse_emits_nine_vectors(self):
        parsed = OINParser().parse(_XERTUK_OIN)
        self.assertEqual(parsed.geo_code, "TCT")
        self.assertEqual(len(parsed.vectors), 9)

    def test_convert_to_msmiles_maps_to_tricapped_prism(self):
        parsed = OINParser().parse(_XERTUK_OIN)
        # This raised "Geometry code 'NON' not supported" before the fix.
        msmiles = convert_parsed_to_msmiles(parsed)
        self.assertIn("9_tricapped_trigonal_prismatic", msmiles)
        # Nine binding atoms -> nine 1-based MetalloGen map numbers :1..:9.
        for n in range(1, 10):
            self.assertIn(f":{n}]", msmiles)


if __name__ == "__main__":
    unittest.main()
