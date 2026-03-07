import unittest
import subprocess
import os

class TestCLI(unittest.TestCase):
    def test_cli_help(self):
        result = subprocess.run(["uv", "run", "python", "-m", "oinsmiles.cli", "--help"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn("xyz2oin", result.stdout)
        self.assertIn("oin2xyz", result.stdout)

    def test_cli_xyz2oin_success(self):
        # Use existing cisplatin fixture
        fixture_path = "tests/fixtures/cisplatin.xyz"
        if not os.path.exists(fixture_path):
            self.skipTest("Fixture not found")
            
        result = subprocess.run(["uv", "run", "python", "-m", "oinsmiles.cli", "xyz2oin", fixture_path], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn("[Pt_SPL]", result.stdout)

    def test_cli_xyz2oin_file_not_found(self):
        result = subprocess.run(["uv", "run", "python", "-m", "oinsmiles.cli", "xyz2oin", "/nonexistent_path.xyz"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn("FileNotFoundError", result.stderr)

    def test_cli_oin2xyz_invalid(self):
        result = subprocess.run(["uv", "run", "python", "-m", "oinsmiles.cli", "oin2xyz", "invalid_oin_string"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn("Error", result.stderr)

    def test_cli_oin2xyz_valid(self):
        # US-002, Novel Test 5
        oin = "[Pt_SPL].N{0}.N{1}.[Cl]{2}.[Cl]{3}"
        result = subprocess.run(["uv", "run", "python", "-m", "oinsmiles.cli", "oin2xyz", oin], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn("Pt", result.stdout)
        self.assertIn("Cl", result.stdout)
        
        # Save candidate artifact
        os.makedirs("tests/candidate_outputs", exist_ok=True)
        with open("tests/candidate_outputs/cli_cisplatin.xyz", "w") as f:
            f.write(result.stdout)

if __name__ == '__main__':
    unittest.main()
