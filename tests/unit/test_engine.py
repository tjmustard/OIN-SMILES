import unittest
from unittest.mock import MagicMock, patch
import sys
import os

from oinsmiles.generation.engine import OIN3DGenerator

class TestOIN3DGenerator(unittest.TestCase):
    def setUp(self):
        self.generator = OIN3DGenerator()
        # Mock dependencies
        self.generator.parser = MagicMock()
        self.generator.adapter = MagicMock()

    def test_generate_flow(self):
        oin_string = "fake_oin"
        extra_params = {"optimize": True}
        
        # Setup mocks
        mock_parsed = MagicMock()
        self.generator.parser.parse.return_value = mock_parsed
        
        mock_structure = "XYZ CONTENT"
        self.generator.adapter.generate.return_value = mock_structure

        # Execute
        result = self.generator.generate(oin_string, extra_params=extra_params)

        # Verify interactions
        self.generator.parser.parse.assert_called_with(oin_string)
        self.generator.adapter.generate.assert_called_with(mock_parsed)
        
        self.assertEqual(result, mock_structure)

if __name__ == '__main__':
    unittest.main()
