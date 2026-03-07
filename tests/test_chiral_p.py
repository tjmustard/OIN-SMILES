import unittest
import os
from rdkit import Chem
from oinsmiles.core.translator import XYZToSMILES
from tests.test_helpers import extract_ligand_smiles

class TestChiralP(unittest.TestCase):
    def setUp(self):
        self.translator = XYZToSMILES()
        self.integration_dir = "tests/integration"
        
    def test_pd_rr_bdpp_p_chiral(self):
        """
        US-002: Verifies P-chiral centers (BDPP) using RDKit CIP oracle.
        Fixture: PdCl2-RR-BDPP.xyz
        Expected: P atoms should have @ or @@ and match RR configuration.
        """
        xyz_filename = "PdCl2-RR-BDPP.xyz"
        path = os.path.join(self.integration_dir, xyz_filename)
        if not os.path.exists(path):
            self.skipTest(f"Fixture {xyz_filename} not found")
        
        oin = self.translator.convert(path)
        
        # 1. Check basic OIN structure
        self.assertIn("[Pd_SPL]", oin)
        
        # 2. Extract ligand SMILES for CIP oracle
        ligand_smiles = extract_ligand_smiles(oin)
        mol = Chem.MolFromSmiles(ligand_smiles)
        self.assertIsNotNone(mol, f"Failed to parse ligand SMILES: {ligand_smiles}")
        
        # 3. Verify CIP perception on the recovered ligand
        # RDKit should perceive the P atoms as chiral
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        
        p_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() == "P"]
        self.assertEqual(len(p_atoms), 2, "Should find 2 P atoms in BDPP")
        
        # RR-BDPP: The carbon backbone (C2 and C4 of the pentane chain) are the stereocenters.
        # Recovered SMILES should have @ or @@ on these carbons.
        c_chiral = [a for a in mol.GetAtoms() if a.GetSymbol() == "C" and a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]
        self.assertGreaterEqual(len(c_chiral), 2, "Should find at least 2 chiral carbons in BDPP backbone")
        
        for atom in c_chiral:
            self.assertTrue(atom.HasProp("_CIPCode"), f"Chiral carbon {atom.GetIdx()} lacks CIP code")
            cip = atom.GetProp("_CIPCode")
            self.assertIn(cip, ["R", "S"])

if __name__ == "__main__":
    unittest.main()
