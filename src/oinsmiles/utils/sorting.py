"""
Deterministic sorting logic for OIN ligands.
Waterfall: Fragment MW (Desc) -> Binding Atom Mass (Desc) -> SMILES (Alpha)
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from typing import List, Tuple, Any

def get_binding_atom_mass(mol: Chem.Mol, binding_indices: List[int]) -> float:
    """Returns the maximum mass of any binding atom."""
    max_mass = 0.0
    for idx in binding_indices:
        if idx >= mol.GetNumAtoms():
            continue
        atom = mol.GetAtomWithIdx(idx)
        # atomic mass, not number, although number is proxy for mass usually.
        # GetMass returns float.
        m = atom.GetMass()
        if m > max_mass:
            max_mass = m
    return max_mass

def get_ligand_sort_key(item: Tuple[Chem.Mol, List[int]]) -> Tuple[float, float, str]:
    """
    Generates a sort key for a ligand tuple (Mol, BindingIndices).
    
    Returns:
        tuple: (-MW, -BindingMass, SMILES)
        Negated for Descending sort on first two keys.
    """
    mol, binding_indices = item
    
    # 1. Fragment Molecular Weight
    mw = Descriptors.MolWt(mol)
    
    # 2. Binding Atom Mass
    b_mass = get_binding_atom_mass(mol, binding_indices)
    
    # 3. Canonical SMILES
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    
    return (-mw, -b_mass, smiles)

def sort_ligands(ligands: List[Tuple[Chem.Mol, List[int]]]) -> List[Tuple[Chem.Mol, List[int]]]:
    """
    Sorts a list of ligands deterministically.
    
    Args:
        ligands: List of (RDKit Mol, List of Binding Atom Indices)
        
    Returns:
        Sorted list.
    """
    return sorted(ligands, key=get_ligand_sort_key)
