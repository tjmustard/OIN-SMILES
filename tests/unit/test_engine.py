import os
import sys
import unittest
from unittest.mock import MagicMock

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.engine import OIN3DGenerator


class TestOIN3DGenerator(unittest.TestCase):
    def setUp(self):
        self.generator = OIN3DGenerator()
        # Mock internal dependencies
        self.generator.parser = MagicMock()
        self.generator.adapter = MagicMock()

    def test_generate_flow(self):
        """generate() should parse the OIN string then delegate to adapter.generate()."""
        oin_string = "fake_oin"
        mock_parsed = MagicMock()
        expected_xyz = "5\n\nPt 0.0 0.0 0.0\n..."

        self.generator.parser.parse.return_value = mock_parsed
        self.generator.adapter.generate.return_value = expected_xyz

        result = self.generator.generate(oin_string)

        self.generator.parser.parse.assert_called_once_with(oin_string)
        self.generator.adapter.generate.assert_called_once_with(mock_parsed)
        self.assertEqual(result, expected_xyz)


if __name__ == "__main__":
    unittest.main()
