import pytest
from rdkit import Chem
from oinsmiles.core.chirality import CIPAssigner, ChiralityRecoveryUtility

def test_cip_assignment_and_recovery():
    # 1. Create a chiral molecule
    # Cl[C@H](F)I is (S)
    mol = Chem.MolFromSmiles("Cl[C@H](F)I")
    
    # 2. Assign CIP
    CIPAssigner.assign_all(mol)
    
    # Check property exists
    atom = mol.GetAtomWithIdx(1)
    assert atom.HasProp("OIN_CIP")
    assert atom.GetProp("OIN_CIP") == "S"
    
    # 3. Simulate fragmentation loss (clear ChiralTag)
    test_mol = Chem.Mol(mol)
    for a in test_mol.GetAtoms():
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    
    # Verify it's lost
    assert test_mol.GetAtomWithIdx(1).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
    
    # 4. Recover
    ChiralityRecoveryUtility.recover(test_mol)
    
    # Verify it's back to CCW (2)
    assert test_mol.GetAtomWithIdx(1).GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    
def test_p_chiral_assignment():
    # Cis-Platin style P-chiral simulation
    # Simple P stereocenter
    mol = Chem.MolFromSmiles("F[P@@H](Cl)Br")
    CIPAssigner.assign_all(mol)
    
    atom = mol.GetAtomWithIdx(1)
    assert atom.HasProp("OIN_CIP")
    cip = atom.GetProp("OIN_CIP")
    assert cip in ["S", "R"]
    
    # Recover and re-gen SMILES
    ChiralityRecoveryUtility.recover(mol)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    assert "@" in smiles
