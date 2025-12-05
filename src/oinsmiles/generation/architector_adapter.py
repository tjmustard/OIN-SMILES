from typing import Dict, List, Any, Optional
from .oin_parser import ParsedOIN, OINVector
import numpy as np
from rdkit import Chem
try:
    from mendeleev import element
except ImportError:
    element = None

class ArchitectorAdapter:
    def __init__(self, scaling_factor: float = 1.0):
        self.scaling_factor = scaling_factor

    def convert(self, parsed_oin: ParsedOIN) -> Dict[str, Any]:
        """
        Converts ParsedOIN to arguments for architector.build_complex
        """
        # 1. Extract Metal
        metal_frag = parsed_oin.fragments[parsed_oin.metal_fragment_idx]
        metal_symbol = self._get_metal_symbol(metal_frag)
        
        # 2. Prepare Ligands and Geometry
        ligands = []
        coordination_geometry = {}
        
        # Group vectors by fragment
        frag_vectors = {}
        for vec in parsed_oin.vectors:
            if vec.fragment_idx not in frag_vectors:
                frag_vectors[vec.fragment_idx] = []
            frag_vectors[vec.fragment_idx].append(vec)
            
        all_site_coords = []
        global_site_idx = 0
        
        # Iterate fragments (skip metal)
        for i, frag_smiles in enumerate(parsed_oin.fragments):
            if i == parsed_oin.metal_fragment_idx:
                continue
            
            vectors = frag_vectors.get(i, [])
            
            # If no vectors, skip? Or add as outer sphere?
            # PRD implies we are building the coordination sphere.
            if not vectors:
                continue
                
            # Sort vectors by atom_in_fragment_idx to ensure consistent order if needed
            vectors.sort(key=lambda v: v.atom_in_fragment_idx)
            
            connect_atoms = []
            lig_map = []
            
            for v in vectors:
                connect_atoms.append(v.atom_in_fragment_idx)
                
                # Calculate geometry position
                ligand_atom_symbol = self._get_atom_symbol(frag_smiles, v.atom_in_fragment_idx)
                dist = self._get_bond_distance(metal_symbol, ligand_atom_symbol) * self.scaling_factor
                
                # Scale vector
                # Vector is x,y,z.
                pos = np.array(v.vector) * dist
                
                all_site_coords.append(pos.tolist())
                lig_map.append([v.atom_in_fragment_idx, global_site_idx])
                global_site_idx += 1
            
            # Create Ligand object representation
            # We will pass this to Wrapper, which will instantiate Architector objects
            ligands.append({
                "smiles": frag_smiles,
                "coordinating_atoms": connect_atoms,
                "coordList": lig_map # List of [atom_idx, site_idx]
            })
            
        return {
            "metal": metal_symbol,
            "ligands": ligands,
            "parameters": {
                "site_coords": all_site_coords
            }
        }

    def _get_metal_symbol(self, smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.GetNumAtoms() == 1:
            return mol.GetAtomWithIdx(0).GetSymbol()
        return smiles.replace("[", "").replace("]", "")

    def _get_atom_symbol(self, smiles: str, atom_idx: int) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            if atom_idx < mol.GetNumAtoms():
                return mol.GetAtomWithIdx(atom_idx).GetSymbol()
        return "C" # Fallback

    def _get_bond_distance(self, metal_symbol: str, ligand_symbol: str) -> float:
        if element is None:
            return 2.0 # Fallback
        
        try:
            r_metal = element(metal_symbol).covalent_radius_pyykko
            r_ligand = element(ligand_symbol).covalent_radius_pyykko
            
            # Convert pm to Angstrom
            if r_metal is None: r_metal = 130.0
            if r_ligand is None: r_ligand = 70.0
            
            return (r_metal + r_ligand) / 100.0
        except Exception:
            return 2.0
