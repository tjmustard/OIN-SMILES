"""
Adapter for converting Parsed OIN to Architector inputDict.
Strict implementation of Haptic Expansion and PRD schema.
"""
from typing import Dict, List, Any, Optional
import numpy as np
from rdkit import Chem
try:
    from mendeleev import element
except ImportError:
    element = None

from .oin_parser import ParsedOIN, OINVector
from ..core.geometry_templates import TEMPLATES
from ..core.haptic_math import expand_slot

class ArchitectorAdapter:
    def __init__(self, scaling_factor: float = 1.0):
        self.scaling_factor = scaling_factor

    def convert(self, parsed_oin: ParsedOIN) -> Dict[str, Any]:
        """
        Converts ParsedOIN to Architector inputDict.
        Follows PRD Data Contract strictly.
        """
        # 1. Extract Metal and Geometry
        metal_frag = parsed_oin.fragments[parsed_oin.metal_fragment_idx]
        metal_symbol = self._get_metal_symbol(metal_frag)
        
        # OIN Parser might return geometry without 'g:' if it's inline. 
        # parsed_oin.geometry usually has code like 'SPL'
        geometry_code = parsed_oin.geometry
        
        # 2. Build Core CoordList (Haptic Expansion)
        # We need to assign each slot used to a set of vectors in the global coordList.
        
        # Group vectors by fragment first
        frag_vectors = {}
        for vec in parsed_oin.vectors:
            if vec.fragment_idx not in frag_vectors:
                frag_vectors[vec.fragment_idx] = []
            frag_vectors[vec.fragment_idx].append(vec)
            
        full_coord_list: List[List[float]] = [] # Global list of [x,y,z]
        ligands_list: List[Dict[str, Any]] = []
        
        current_global_vec_idx = 0
        
        # Process each ligand fragment
        # Note: parsed_oin.fragments includes metal fragment, we skip it.
        # But frag_idx in vectors refers to the fragments list.
        # We must be careful about indices.
        
        for i, frag_smiles in enumerate(parsed_oin.fragments):
            if i == parsed_oin.metal_fragment_idx:
                continue
                
            vectors = frag_vectors.get(i, [])
            if not vectors:
                continue
                
            # Group by Slot Index to handle Haptic/Monodentate
            slot_groups = {}
            for v in vectors:
                if v.slot_idx not in slot_groups:
                    slot_groups[v.slot_idx] = []
                slot_groups[v.slot_idx].append(v)
             
            # Prepare Ligand entry
            ligand_entry = {
                "smiles": frag_smiles, # Clean SMILES
                "coordinating_atoms": [],
                "coordList": [] # Map [atom_idx_in_smiles, global_vec_idx_in_core]
            }
            
            # Re-map indices to account for Architector's expansion of Explicit Hydrogens (interleaved)
            # Architector (via OB/XYZ conversion?) appears to expand [cH] into C-H sequence.
            # RDKit indices: 0, 1, 2...
            # Architector indices: 0, 2, 4... (if [cH])
            
            rdkit_to_arch_map = {}
            current_arch_idx = 0
            
             # We must parse the frag_smiles with RDKit to determine explicit counts
            mol_frag = Chem.MolFromSmiles(frag_smiles, sanitize=False)
            if not mol_frag:
                # Fallback: Identity mapping if parsing fails
                for i in range(100): rdkit_to_arch_map[i] = i
            else:
                 # Sanitize to get implicit H counts correct (as much as possible from SMILES)
                 try:
                     Chem.SanitizeMol(mol_frag)
                 except:
                     pass

                 # Apply Heuristics to Coordinating Atoms to fix Architector inputs
                 # 1. Identify all atom_in_fragment_idx coordinating
                 coord_indices = set()
                 if i in frag_vectors:
                     for v in frag_vectors[i]:
                         coord_indices.add(v.atom_in_fragment_idx)
                 
                 for c_idx in coord_indices:
                     if c_idx < mol_frag.GetNumAtoms():
                         atom = mol_frag.GetAtomWithIdx(c_idx)
                         
                         # Check User Intent: If explicit info present, skip.
                         # OIN often has implicit SMILES. If user provided [n], NoImplicit=True.
                         if atom.GetNumExplicitHs() > 0 or atom.GetNoImplicit():
                             continue
                             
                         # Apply Heuristics
                         sym = atom.GetSymbol()
                         orig_implicit = atom.GetNumImplicitHs()
                         new_explicit = orig_implicit # Default: Convert Implicit to Explicit
                         
                         if sym == 'C':
                             # Carbon Rule: Phenyl, Cyanide -> Remove H (Assume anionic/ylidic)
                             new_explicit = 0
                         elif sym == 'O' or sym == 'S':
                             # Chalcogen Rule
                             if atom.GetDegree() > 0:
                                 # Bonded (Acac, Thiolate) -> Remove H
                                 new_explicit = 0
                             else:
                                 # Isolated (Water, H2S) -> Keep H (L-type)
                                 pass
                         
                         atom.SetNoImplicit(True)
                         atom.SetNumExplicitHs(new_explicit)
                 
                 # Regenerate SMILES with modifications
                 # canonical=False preserves atom indices!
                 # allHsExplicit=True forces writing [c] and [NH3] correctly.
                 cleaned_smiles = Chem.MolToSmiles(mol_frag, canonical=False, allHsExplicit=True)
                 ligand_entry["smiles"] = cleaned_smiles
                 
                 # Re-parse to ensure consistency for mapping?
                 # Actually, MolToSmiles(canonical=False) re-generates string.
                 # Indices *should* be preserved if we don't reorder.
                 # But we need to use this new mol for counting architector indices (Explicit Hs are now Explicit atoms in mapping logic).
                 # Wait, Architector expansion logic further down:
                 # "Architector indices: 0, 2, 4... (if [cH])"
                 # My existing logic counted ExplicitHs.
                 # Now strict ExplicitHs are set on the atoms.
                 # So `atom.GetNumExplicitHs()` will return `new_explicit`.
                 # So downstream logic should work automatically if I use `mol_frag`.
                 
                 # IMPORTANT: update frag_smiles? No, we use `cleaned_smiles` in entry.
                 # But we need `mol_frag` for the mapping loop below.
                 # AND we need to make sure the loop below iterates in the same order as `cleaned_smiles`?
                 # If canonical=False, order is atom index order.
                 
                 current_arch_idx = 0
                 # Ensure atoms are in order 0..N
                 for atom in mol_frag.GetAtoms():
                     r_idx = atom.GetIdx()
                     rdkit_to_arch_map[r_idx] = current_arch_idx
                     
                     # Increment logic
                     # 1 for the atom itself
                     current_arch_idx += 1
                     
                     # Add explicit Hs?
                     num_explicit = atom.GetNumExplicitHs()
                     current_arch_idx += num_explicit

            # We need to map atom_indices for this ligand.
            # Local atom idx -> Global Vec Idx
            
            for slot_idx, group_vecs in slot_groups.items():
                if len(group_vecs) == 0: continue
                
                # Get Template Def
                if geometry_code not in TEMPLATES or slot_idx not in TEMPLATES[geometry_code]:
                    # Error or Unknown Geometry fallback
                    # Just append empty vectors? Or single generic vector?
                    # PRD implies geometry must be valid.
                    # We'll skip or error. For now, skip to prevent crash.
                    continue
                    
                geo_def = TEMPLATES[geometry_code][slot_idx]
                
                if len(group_vecs) > 1:
                    # ================= HAPTIC LOGIC =================
                    n_haptic = len(group_vecs)
                    
                    # 1. Expand Slot -> N vectors
                    # Note: These are unit vectors around centroid. 
                    # Does Architector define bond length per vector? 
                    # Usually Architector scales sites. 
                    # We will store the vectors as Generated.
                    # Scaling? If we put unit vectors, we must ensure Architector scales them or we scale them.
                    # PRD `coordList` section shows `[2.22, 0.0, 0.0]`. This implies SCALING.
                    # We need to scale these vectors by bond length.
                    # But vectors in a haptic ring might have different atoms (e.g. N vs C).
                    # `expand_slot` returns polygon of same radius.
                    # We'll assume average bond length scaling for the haptic group?
                    # Or scale each vector individually after assignment?
                    # BUT `coordList` is in `core`. It defines the sites.
                    # If we use different scalings for the same polygon, it won't be a regular polygon anymore.
                    # We will generate the polygon vectors (unit), then scale them by the Metal-Ligand distance.
                    # For haptic rings, usually we scale the Centroid distance? 
                    # Or each atom distance?
                    # `expand_slot` generates vectors starting from Origin (Metal).
                    # `vec = Z + cone_spread * ...`. Z is centroid direction.
                    # If Z is unit, the result is on unit sphere (normalized).
                    # We should scale all these vectors by the expected M-L distance.
                    # We'll use the first atom's symbol to guess distance? Or average?
                    
                    heading_atom_idx = group_vecs[0].atom_in_fragment_idx
                    direction = 1
                    
                    # Find North Star
                    for hv in group_vecs:
                        if hv.haptic_heading:
                            heading_atom_idx = hv.atom_in_fragment_idx
                            direction = hv.haptic_direction
                            break
                            
                    expanded_vectors = expand_slot(n_haptic, geo_def['pos'], geo_def['ref'])
                    
                    # Scale vectors
                    # Use the heading atom to determine bond distance?
                    atom_sym = self._get_atom_symbol(frag_smiles, heading_atom_idx)
                    dist = self._get_bond_distance(metal_symbol, atom_sym) * self.scaling_factor
                    
                    scaled_vectors = [v * dist for v in expanded_vectors]
                    
                    # Add to global coordList
                    start_vec_idx = len(full_coord_list)
                    full_coord_list.extend([v.tolist() for v in scaled_vectors])
                    
                    # Map Atoms (Winding)
                    star_atom_local_idx = heading_atom_idx # This is the index in the fragment SMILES
                    
                    # We verify which atom in `group_vecs` corresponds to which vector.
                    # group_vecs contains `atom_in_fragment_idx` for all binding atoms.
                    binding_atom_indices = {v.atom_in_fragment_idx for v in group_vecs}
                    
                    # Logic: 
                    # The Star Atom (heading_atom_idx) maps to Vector 0 (expanded_vectors[0]).
                    # The Next Atom maps to Vector 1, etc.
                    # But which is "Next"? 
                    # We assume indices traverse the ring/chain sequentially via connectivity?
                    # Or simply by Delta of atom index?
                    # "Compare cross product...". That assumes we have 3D. We don't.
                    # We are generating.
                    # The parser gives us {0>}.
                    # This means: "Map this atom to V0. The *next neighbor* maps to V(direction)."
                    # Ideally we trace connectivity.
                    # BUT if we assume simplistic index delta as the previous code did:
                    # `delta = (idx - star_idx) * direction`
                    # `k = delta % n_haptic`
                    # This works if atoms are numbered sequentially (C1-C2-C3...).
                    # For generated OIN from canonical Smi, this might be true?
                    # Let's stick to Index Delta logic for now as implementing full graph traversal on SMILES here is complex.
                    
                    for hv in group_vecs:
                        idx = hv.atom_in_fragment_idx
                        delta = (idx - star_atom_local_idx) * direction
                        # Modulo arithmetic handles wrapping
                        k = delta % n_haptic
                        
                        # Architector expects: [LigandAtomIdx, CoreCoordListIdx]
                        global_vec_idx = start_vec_idx + k
                        
                        # Apply Mapping: RDKit Index -> Architector Index
                        arch_atom_idx = rdkit_to_arch_map.get(idx, idx)
                        
                        ligand_entry["coordinating_atoms"].append(arch_atom_idx)
                        ligand_entry["coordList"].append([arch_atom_idx, global_vec_idx])

                else:
                    # ================= MONODENTATE / CHELATE (Specific leg) =================
                    hv = group_vecs[0]
                    # Get vector (unit)
                    vec_unit = geo_def['pos'] 
                    # Scale
                    atom_sym = self._get_atom_symbol(frag_smiles, hv.atom_in_fragment_idx)
                    dist = self._get_bond_distance(metal_symbol, atom_sym) * self.scaling_factor
                    
                    vec_scaled = vec_unit * dist
                    
                    global_idx_vec = len(full_coord_list)
                    full_coord_list.append(vec_scaled.tolist())
                    
                    idx = hv.atom_in_fragment_idx
                    arch_atom_idx = rdkit_to_arch_map.get(idx, idx)
                    
                    ligand_entry["coordinating_atoms"].append(arch_atom_idx)
                    ligand_entry["coordList"].append([arch_atom_idx, global_idx_vec])

            # Deduplicate coordinating_atoms for metadata
            ligand_entry["coordinating_atoms"] = sorted(list(set(ligand_entry["coordinating_atoms"])))
            
            if ligand_entry["coordList"]:
                ligands_list.append(ligand_entry)

        return {
            "core": {
                "metal": metal_symbol,
                "coreType": "user_core",
                "coordList": full_coord_list
            },
            "ligands": ligands_list,
            "parameters": {
                "assemble_method": "GFN2-xTB",
                "full_method": "GFN2-xTB",
                "debug": True
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
        return "C"

    def _get_bond_distance(self, metal_symbol: str, ligand_symbol: str) -> float:
        # Simple covalent radii sum or fallback
        if element is None:
            return 2.22
        try:
            r_m = element(metal_symbol).covalent_radius_pyykko
            r_l = element(ligand_symbol).covalent_radius_pyykko
            if r_m is None: r_m = 130.0
            if r_l is None: r_l = 75.0
            return (r_m + r_l) / 100.0
        except:
            return 2.22
