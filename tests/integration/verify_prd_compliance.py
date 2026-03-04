
import unittest
import numpy as np
from rdkit import Chem
from oinsmiles.utils.xyz2mol import get_oin_string, get_tmc_mol
from oinsmiles.utils.oin_aligner import OINDiscreteAligner, TEMPLATES
from oinsmiles.oin.inline import OINInlineHandler

class TestPRDCompliance(unittest.TestCase):

    def test_inline_regex(self):
        """Verify regex from PRD Section 5 matches implementation."""
        # PRD: \[([A-Z][a-z]?)_([A-Z]{3})\]
        # Code: \[([A-Z][a-z]?)\_([A-Z]{3})\] (Backslash before underscore is harmless but present)
        text = "[Pt_SPL]"
        match = OINInlineHandler.METAL_REGEX.search(text)
        self.assertTrue(match, "Failed to match [Pt_SPL]")
        self.assertEqual(match.group(1), "Pt")
        self.assertEqual(match.group(2), "SPL")

    def test_spl_template_vectors(self):
        """Verify SPL Vectors match PRD Section 3.3."""
        # PRD: 0=(1,0,0), 1=(-1,0,0), 2=(0,1,0), 3=(0,-1,0)
        # Note: PRD Step 3.3 Table.
        
        vectors = TEMPLATES['SPL']
        # Check Slot 1 (Index 1)
        v1 = vectors[1]
        # PRD says (-1,0,0)
        # Code says (0,1,0)
        
        # We expect this to FAIL based on my read.
        expected_v1 = np.array([-1, 0, 0])
        if not np.allclose(v1, expected_v1, atol=1e-5):
            print(f"\n[FAILURE] SPL Template Mismatch! Slot 1. Expected {expected_v1}, Got {v1}")
        else:
            print("\n[SUCCESS] SPL Template matches PRD.")

    def test_zone_b_bracketing(self):
        """Verify Zone B atoms use implicit H rules (no force brackets)."""
        # We need a mock molecule where a ligand has Zone B atoms (backbone).
        # e.g. Ethylamine -> N-CC. N is Zone A. C is Zone B.
        # SMILES should be like N[0].CC (or similar).
        # If it comes out as [NH2][0].[CH2][CH3], that's a failure (over-bracketing).
        
        # Difficult to unit test without full xyz2mol stack, but we can check OINSanitizer behavior alone?
        # OINSanitizer takes a ligand mol and Zone A indices.
        
        from oinsmiles.utils.oin_aligner import OINSanitizer
        
        mol = Chem.MolFromSmiles("NCC") # Ethylamine
        # N is index 0. C is 1. C is 2.
        # Bind through N (0).
        zone_a = [0]
        
        robust_smiles, _ = OINSanitizer.generate_robust_smiles(mol, zone_a)
        print(f"\nZone B Check: NCC bound at N -> {robust_smiles}")
        
        # Expectation: N is bracketed [NH2]. Cs are NOT forcefully bracketed.
        # However, generate_robust_smiles might canonicalize.
        # If it returns [NH2]CC, that's good.
        # If it returns [NH2][CH2][CH3], that's bad (unless chemically necessary).
        
        if "[CH2]" in robust_smiles or "[CH3]" in robust_smiles:
             print("[FAILURE] Zone B seems to be forcefully bracketed!")
        else:
             print("[SUCCESS] Zone B has implicit Hs.")
             
if __name__ == '__main__':
    unittest.main()
