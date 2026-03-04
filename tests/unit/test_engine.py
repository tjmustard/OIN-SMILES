import unittest
from unittest.mock import MagicMock, patch
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.generation.engine import OIN3DGenerator

class TestOIN3DGenerator(unittest.TestCase):
    def setUp(self):
        self.generator = OIN3DGenerator()
        # Mock dependencies
        self.generator.parser = MagicMock()
        self.generator.adapter = MagicMock()
        self.generator.wrapper = MagicMock()

    def test_generate_flow(self):
        oin_string = "fake_oin"
        extra_params = {"optimize": True}
        
        # Setup mocks
        mock_parsed = MagicMock()
        self.generator.parser.parse.return_value = mock_parsed
        
        mock_arch_args = {"metal": "Pt", "ligands": []}
        self.generator.adapter.convert.return_value = mock_arch_args
        
        mock_structure = MagicMock()
        self.generator.wrapper.run.return_value = mock_structure

        # Execute
        result = self.generator.generate(oin_string, extra_params=extra_params)

        # Verify interactions
        self.generator.parser.parse.assert_called_with(oin_string)
        self.generator.adapter.convert.assert_called_with(mock_parsed)
        self.generator.wrapper.run.assert_called()
        
        # Verify params merged
        args, kwargs = self.generator.wrapper.run.call_args
        # Kwargs should be mock_arch_args + extra_params
        self.assertEqual(kwargs['metal'], "Pt")
        self.assertEqual(kwargs['parameters']['optimize'], True)
        self.assertEqual(result, mock_structure)

if __name__ == '__main__':
    unittest.main()
