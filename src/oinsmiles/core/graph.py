"""
Graph representation and wrappers for molecule generation.
"""
from enum import Enum, auto
from typing import List, Tuple, Optional, Dict
from rdkit import Chem
import numpy as np

# We assume xyz2mol is available in utils
try:
    from ..utils import xyz2mol
except ImportError:
    xyz2mol = None
    
class BondType(Enum):
    SINGLE = auto()
    DOUBLE = auto()
    TRIPLE = auto()
    AROMATIC = auto()
    DATIVE = auto() # Important for Metals

class AtomNode:
    def __init__(self, element: str, coords: Tuple[float, float, float]):
        self.element = element
        self.coords = coords
        self.id = -1 # Assigned when adding to graph

class TMCGraph:
    """
    Transition Metal Complex Graph.
    Lightweight intermediate representation.
    """
    def __init__(self):
        self.atoms: List[AtomNode] = []
        self.bonds: List[Tuple[int, int, BondType]] = []
        
    def add_atom(self, element: str, coords: Tuple[float, float, float]) -> int:
        node = AtomNode(element, coords)
        node.id = len(self.atoms)
        self.atoms.append(node)
        return node.id
        
    def add_bond(self, src: int, tgt: int, bond_type: BondType):
        self.bonds.append((src, tgt, bond_type))
        
    def to_rdkit(self) -> Chem.Mol:
        """Converts graph to RDKit Mol (RW)."""
        mol = Chem.RWMol()
        node_idx_map = {}
        
        conf = Chem.Conformer()
        
        for i, node in enumerate(self.atoms):
            idx = mol.AddAtom(Chem.Atom(node.element))
            node_idx_map[i] = idx
            
            # Set coords
            conf.SetAtomPosition(idx, node.coords)
            
        mol.AddConformer(conf)
            
        for src, tgt, btype in self.bonds:
            if src in node_idx_map and tgt in node_idx_map:
                # Map BondType to RDKit
                r_type = Chem.BondType.SINGLE
                if btype == BondType.DOUBLE: r_type = Chem.BondType.DOUBLE
                elif btype == BondType.TRIPLE: r_type = Chem.BondType.TRIPLE
                elif btype == BondType.AROMATIC: r_type = Chem.BondType.AROMATIC
                elif btype == BondType.DATIVE: r_type = Chem.BondType.DATIVE
                
                mol.AddBond(node_idx_map[src], node_idx_map[tgt], r_type)
                
        return mol.GetMol()

def mol_from_xyz_file(path: str, charge: int = 0) -> Chem.Mol:
    """
    Wraps xyz2mol to read XYZ file and return RDKit Mol.
    """
    if xyz2mol is None:
        raise ImportError("xyz2mol module not found in utils.")
        
    atoms, charge, xyz_coordinates = xyz2mol.read_xyz_file(path)
    # xyz2mol logic:
    # mols = xyz2mol(atoms, xyz_coordinates, charge=charge, ...)
    mols = xyz2mol.xyz2mol(atoms, xyz_coordinates, charge=charge, use_graph=True)
    if not mols:
        raise ValueError("xyz2mol failed to generate molecule.")
    return mols[0]
