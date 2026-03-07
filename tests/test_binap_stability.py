import unittest
import os
from oinsmiles.core.translator import XYZToSMILES

class TestBinapStability(unittest.TestCase):
    def setUp(self):
        self.translator = XYZToSMILES()
        self.integration_dir = "tests/integration"
        
    def test_binap_stability(self):
        """
        US-001: Verifies BINAP stability without @ assertions on P atoms.
        (Using PdCl2-R-BINAP.xyz)
        """
        xyz_filename = "PdCl2-R-BINAP.xyz"
        path = os.path.join(self.integration_dir, xyz_filename)
        if not os.path.exists(path):
            self.skipTest(f"Fixture {xyz_filename} not found")
        
        oin = self.translator.convert(path)
        
        # We expect a valid OIN string. 
        # For BINAP, it should fragment into Pd, 2 Cl, and 1 BINAP ligand.
        # Pd-BINAP complex is usually [Pd_SPL].(BINAP).Cl.Cl
        
        self.assertIn("[Pd_SPL]", oin)
        self.assertIn("[Cl]{2}", oin)
        self.assertIn("[Cl]{3}", oin)
        
        # Check that BINAP SMILES is present (contains P{0} and P{1})
        self.assertIn("P{0}", oin)
        self.assertIn("P{1}", oin)
        
        # US-001 negative constraint: No @ on P atoms in this baseline test 
        # (Since we are testing the complex stability, not necessarily the ligand chirality recovery yet)
        # Note: If P-chiral atoms were present and perceived, they would have @. 
        # BINAP P atoms are not stereocenters themselves (they are symmetrical P-Ph2).
        # The axial chirality is the key.
        self.assertNotIn("P@", oin, "P atoms should not have chiral descriptors in BINAP")

if __name__ == "__main__":
    unittest.main()
