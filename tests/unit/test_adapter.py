"""
Tests for ArchitectorAdapter generation.
"""
import pytest
import numpy as np
from oinsmiles.generation.architector_adapter import ArchitectorAdapter
from oinsmiles.generation.oin_parser import ParsedOIN, OINVector

def test_adapter_haptic_expansion():
    """
    Test that a Cp (N=5) ligand expands to 5 vectors in core.coordList.
    """
    # Create Mock ParsedOIN matches [Pt_LIN].[cH]{0>}1[cH]...
    # We simulate parsed output manually.
    
    # Metal Fragment
    frags = ["Pt", "c1cccc1"] # 0: Metal, 1: Cp
    
    # Vectors for Cp
    # 5 Atoms. All Haptic. Slot 0.
    # Atom indices in fragment: 0,1,2,3,4.
    vectors = []
    
    # Using LIN geometry (Slot 0 is Top). 
    # NOTE: LIN has 2 slots. 0 and 1.
    # Cp at Slot 0.
    
    for i in range(5):
        # We dummy the 'vector' field as it comes from template resolution in parser.
        # Adapter re-resolves using geometry templates internally if needed?
        # Actually Adapter uses `group_vecs` and looks up TEMPLATES again to get `pos` and `ref`.
        # So we just need valid Slot Idx.
        
        vectors.append(OINVector(
            atom_idx=-1,
            vector=(0,0,1), # Dummy
            fragment_idx=1,
            atom_in_fragment_idx=i,
            haptic_heading=(i==0), # Atom 0 is heading
            haptic_direction=1,
            slot_idx=0
        ))
        
    parsed = ParsedOIN(
        smiles="[Pt_LIN].c1cccc1",
        fragments=frags,
        metal_fragment_idx=0,
        vectors=vectors,
        original_oin="...",
        geometry="LIN"
    )
    
    adapter = ArchitectorAdapter()
    output = adapter.convert(parsed)
    
    # VERIFY SCHEMA
    assert "core" in output
    assert "ligands" in output
    assert "parameters" in output
    
    core = output["core"]
    assert core["metal"] == "Pt"
    assert core["coreType"] == "user_core"
    
    # CHECK EXPANSION
    # Expect 5 vectors in coordList (for the 1 ligand).
    assert len(core["coordList"]) == 5
    
    ligand = output["ligands"][0]
    assert ligand["smiles"] == "c1cccc1"
    assert len(ligand["coordinating_atoms"]) == 5
    assert len(ligand["coordList"]) == 5 
    
    # coordList in ligand is [[atom_idx, global_vec_idx], ...]
    # Atoms 0..4 should map to Vectors 0..4
    atom_indices = sorted([param[0] for param in ligand["coordList"]])
    vec_indices = sorted([param[1] for param in ligand["coordList"]])
    
    assert atom_indices == [0,1,2,3,4]
    assert vec_indices == [0,1,2,3,4]

def test_adapter_mixed_ligands():
    """Test Monodentate + Haptic."""
    # [Pt_SPL].[Cl]{0}.[Cp]{1} (Assume SPL supports Haptic at slot 1 - technically SPL is 4 monodentate slots but logic holds)
    
    frags = ["Pt", "Cl", "c1cccc1"]
    vectors = []
    
    # Cl at Slot 0
    vectors.append(OINVector(-1, (1,0,0), 1, 0, False, 1, 0))
    
    # Cp at Slot 1 (N=5)
    for i in range(5):
        vectors.append(OINVector(-1, (0,1,0), 2, i, (i==0), 1, 1))
        
    parsed = ParsedOIN(
        smiles="...",
        fragments=frags,
        metal_fragment_idx=0,
        vectors=vectors,
        original_oin="...",
        geometry="SPL"
    )
    
    adapter = ArchitectorAdapter()
    output = adapter.convert(parsed)
    
    # Cl (1 vector) + Cp (5 vectors) = 6 total vectors
    assert len(output["core"]["coordList"]) == 6
    assert len(output["ligands"]) == 2
