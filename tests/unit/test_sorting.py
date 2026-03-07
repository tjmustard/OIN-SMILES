"""
Tests for deterministic sorting.
"""
import pytest
from rdkit import Chem
from oinsmiles.utils.sorting import sort_ligands

def test_sort_waterfall():
    """
    Test: 1. MW, 2. Binding Mass, 3. SMILES
    """
    # Ligands:
    # A: [Br] (MW ~80)
    # B: [Cl] (MW ~35.5)
    # C: [Pt] (Heavy, but let's assume it's a ligand for sorting test)
    # D: Isomer of A?
    
    # 1. MW Check
    # Br vs Cl
    mol_br = Chem.MolFromSmiles("[Br]")
    mol_cl = Chem.MolFromSmiles("[Cl]")
    
    # List of (Mol, BindingIndices)
    items = [
        (mol_cl, [0]),
        (mol_br, [0])
    ]
    
    sorted_items = sort_ligands(items)
    
    # Expect Br first (MW 80 > 35)
    assert sorted_items[0][0] == mol_br
    assert sorted_items[1][0] == mol_cl
    
def test_sort_binding_atom_mass():
    """
    Ligands with same MW (approx) or override.
    Construct artificial scenario or find real isomers.
    Let's rely on Binding Atom Mass diff.
    Ligand X: C-O (Bind C) vs Ligand Y: C-O (Bind O)
    Ligands identical MW (assume same molecule).
    """
    # Methanol C-O
    mol1 = Chem.MolFromSmiles("CO")
    mol2 = Chem.MolFromSmiles("CO")
    
    # Item 1: Binds C (idx 0). Mass 12.
    # Item 2: Binds O (idx 1). Mass 16.
    
    items = [
        (mol1, [0]), # Binds C
        (mol2, [1])  # Binds O
    ]
    
    sorted_items = sort_ligands(items)
    
    # Expect O-bound first (Mass 16 > 12)
    # MW is equal.
    assert sorted_items[0][1] == [1]
    
def test_sort_smiles_alphabetical():
    """
    Equal MW, Equal Binding Mass.
    Sort by SMILES.
    """
    # Isomers C3H6: Cyclopropane vs Propene?
    # Cyclopropane: C1CC1
    # Propene: CC=C
    # MW is identical (42.08)
    # Binding C (12)
    
    mol_cyc = Chem.MolFromSmiles("C1CC1")
    mol_prop = Chem.MolFromSmiles("CC=C")
    
    items = [
        (mol_prop, [0]),
        (mol_cyc, [0])
    ]
    
    # C1CC1 vs CC=C
    # "C1CC1" < "CC=C" (alphabetically)
    # Waterfall: (Desc, Desc, Asc).
    # Expected: Cyclopropane first.
    
    sorted_items = sort_ligands(items)
    
    assert sorted_items[0][0] == mol_cyc
