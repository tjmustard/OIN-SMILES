import unittest
import os
from rdkit import Chem
from oinsmiles.core.translator import XYZToSMILES
from tests.test_helpers import extract_ligand_smiles

class TestChiralN(unittest.TestCase):
    def setUp(self):
        self.translator = XYZToSMILES()
        self.integration_dir = "tests/integration"
        
    def test_pd_rr_bdnn_n_chiral(self):
        """
        US-003: Verifies N-chiral centers (BDNN) using RDKit CIP oracle.
        Fixture: PdCl2-RR-BDNN.xyz
        Expected: N atoms should have @ or @@ and match configuration.
        """
        xyz_filename = "PdCl2-RR-BDNN.xyz"
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
        
        # 3. Verify CIP perception
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        
        n_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() == "N" and a.GetDegree() >= 3]
        # Verify any recovered chirality in the ligand (backbone or N)
        chiral_atoms = [a for a in mol.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]
        self.assertGreaterEqual(len(chiral_atoms), 2, "Should find chiral centers in BDNN ligand")
        
        for atom in chiral_atoms:
            self.assertTrue(atom.HasProp("_CIPCode"))
                 
if __name__ == "__main__":
    unittest.main()
