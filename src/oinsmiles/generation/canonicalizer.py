"""
Canonicalization Pipeline (XYZ -> OIN).
Orchestrates: PAI Alignment -> Fragmentation -> Mass-First Sort -> Geometry Detection -> OIN Generation.
Migrated from utils/xyz2mol.py.
"""
import numpy as np
import logging
from rdkit import Chem
import re
from typing import List, Tuple, Dict, Any

from ..utils.oin_aligner import _align_to_pai, OINDiscreteAligner
from ..oin.sanitizer import OINSanitizer
# We need OINInlineHandler for the final formatting (V3.0/V4.0 requirement)
try:
    from ..oin.inline import OINInlineHandler
except ImportError:
    OINInlineHandler = None

logger = logging.getLogger(__name__)

class Canonicalizer:
    @staticmethod
    def canonicalize(tmc_mol: Chem.Mol, xyz_coords: List[List[float]]) -> str:
        """
        Generates the Open Isomer Notation (OIN) string for the molecule.
        
        Args:
            tmc_mol: RDKit Molecule object (can contain Metal coordinates, etc.)
            xyz_coords: List of [x,y,z] coordinates corresponding to atoms in tmc_mol.
        
        Returns:
            OIN String (e.g. [Pt_SPL].Cl{0}.Cl{1}...)
        """
        # 1. Identify Metal
        metal_idx = -1
        # Simple TM check logic or assume single metal?
        # We search for Transition Metals (AtomicNum 21-30, 39-48, 57-80)
        TRANSITION_METALS_NUM = [21,22,23,24,25,26,27,57,28,29,30,39,40,41,
                                 42,43,44,45,46,47,48,71,72,73,74,75,76,77,78,79,80]
        
        for atom in tmc_mol.GetAtoms():
            if atom.GetAtomicNum() in TRANSITION_METALS_NUM:
                metal_idx = atom.GetIdx()
                break
                
        if metal_idx == -1:
            raise ValueError("No transition metal found in molecule!")

        # 2. CANONICALIZE ORIENTATION (Translation + PAI Alignment)
        # Updates coordinates to be canonically aligned
        canonical_coords = _align_to_pai(tmc_mol, xyz_coords, metal_idx)
        
        # 3. Fragment Molecule into Ligands
        mol = Chem.RWMol(tmc_mol)
        metal_atom = tmc_mol.GetAtomWithIdx(metal_idx)
        
        # Identify bonds to remove
        metal_bonds = mol.GetAtomWithIdx(metal_idx).GetBonds()
        coordinating_atoms = []
        bonds_to_remove = []
        
        for bond in metal_bonds:
            other_atom = bond.GetOtherAtomIdx(metal_idx)
            coordinating_atoms.append(other_atom)
            bonds_to_remove.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            
        for u, v in bonds_to_remove:
            mol.RemoveBond(u, v)
        
        # Neutralize Charges for clean SMILES generation
        for atom in mol.GetAtoms():
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True) # Relies on Sanitizer later

        # Get Fragments
        frags_indices = Chem.GetMolFrags(mol, asMols=False)
        fragments_data = []

        # Map Mol Idx -> Original XYZ Index
        atom_map_to_xyz = {}
        for atom in mol.GetAtoms():
            if atom.HasProp("__origIdx"):
                atom_map_to_xyz[atom.GetIdx()] = atom.GetIntProp("__origIdx")
            else:
                atom_map_to_xyz[atom.GetIdx()] = atom.GetIdx()
        
        # Process Fragments Data gathering
        for i, indices in enumerate(frags_indices):
            is_metal = (metal_idx in indices)
            
            mass = 0.0
            binding_mass = 0.0
            frag_binding_atoms = [] # (global_idx, mass, coords)
            
            for idx in indices:
                atom = mol.GetAtomWithIdx(idx)
                m = atom.GetMass()
                mass += m
                if idx in coordinating_atoms:
                    binding_mass = max(binding_mass, m)
                    orig_i = atom_map_to_xyz[idx]
                    coords = np.array(canonical_coords[orig_i])
                    frag_binding_atoms.append((idx, m, coords))
            
            metal_coords = None
            if is_metal and indices:
                 m_idx = indices[0]
                 orig_m_idx = atom_map_to_xyz.get(m_idx, m_idx)
                 metal_coords = np.array(canonical_coords[orig_m_idx]) # Origin usually

            fragments_data.append({
                'indices': indices,
                'is_metal': is_metal,
                'mass': mass,
                'binding_mass': binding_mass,
                'metal_coords': metal_coords,
                'smiles': "",
                'binding_atoms': frag_binding_atoms
            })

        # Process each fragment to generate SMILES (Sanitizer)
        for item in fragments_data:
            indices = item['indices']
            frag_binding_atoms = item['binding_atoms']
            
            # Logic to extract clean fragment mol
            # (Simplified from xyz2mol details: extract submol)
            # Re-implementation of extraction logic:
            
            # Identify atoms to keep vs merge
            # RATIONALE: We merge non-binding explicit Hydrogens into the 'NumExplicitHs' count
            # of their heavy atom neighbor. 
            # If we keep them as separate atoms in the RWMol, RDKit will generate SMILES like
            # [C]([H])([H])[H] instead of [CH3]. 
            # Both are chemically valid, but the verbose [H] nodes clutter the OIN string 
            # and break some downstream assumptions.
            # DO NOT REVERT THIS 'MERGE' LOGIC unless you want verbose OIN strings.
            atoms_to_keep = []
            hs_to_merge = []
            
            # Global Set of binding indices for fast lookup
            # frag_binding_atoms is list of (global_idx, mass, coords)
            binding_set = {x[0] for x in frag_binding_atoms}
            
            # Create fresh RWMol
            mw = Chem.RWMol()
            old_to_new = {}

            for old_idx in indices:
                atom = mol.GetAtomWithIdx(old_idx)
                if atom.GetAtomicNum() == 1 and old_idx not in binding_set:
                    hs_to_merge.append(old_idx)
                else:
                    atoms_to_keep.append(old_idx)
            
            # 1. Add Keep Atoms
            for old_idx in atoms_to_keep:
                atom = mol.GetAtomWithIdx(old_idx)
                new_atom = Chem.Atom(atom.GetAtomicNum())
                new_atom.SetFormalCharge(atom.GetFormalCharge())
                new_atom.SetIsotope(atom.GetIsotope())
                new_atom.SetChiralTag(atom.GetChiralTag())
                new_atom.SetNoImplicit(True) # We control H count explicitly
                new_atom.SetNumExplicitHs(0) # Start with 0, increment later
                
                new_idx = mw.AddAtom(new_atom)
                old_to_new[old_idx] = new_idx

            # 2. Merge H Atoms
            for h_idx in hs_to_merge:
                h_atom = mol.GetAtomWithIdx(h_idx)
                # Find neighbor to attach count to
                # In fragments, H should effectively have 1 neighbor
                neighbors = h_atom.GetNeighbors()
                for nbr in neighbors:
                    nbr_idx = nbr.GetIdx()
                    # Check if neighbor is in this fragment AND is a kept atom
                    if nbr_idx in old_to_new:
                        # Increment H count
                        target_atom = mw.GetAtomWithIdx(old_to_new[nbr_idx])
                        target_atom.SetNumExplicitHs(target_atom.GetNumExplicitHs() + 1)
            
            # 3. Add Bonds (Between Kept Atoms)
            for old_idx in atoms_to_keep:
                atom = mol.GetAtomWithIdx(old_idx)
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    # Only add if neighbor is also kept and index > old_idx (avoid duplicates)
                    if nbr_idx in old_to_new and nbr_idx > old_idx:
                        bond = mol.GetBondBetweenAtoms(old_idx, nbr_idx)
                        if bond:
                             mw.AddBond(old_to_new[old_idx], old_to_new[nbr_idx], bond.GetBondType())
            
            frag_mol = mw.GetMol()
            
            # Indentify binding atoms in local frame
            local_binding_indices = []
            for b_item in frag_binding_atoms:
                g_idx = b_item[0]
                if g_idx in old_to_new:
                    local_binding_indices.append(old_to_new[g_idx])
            
            sanitized_smiles = ""
            if not item['is_metal']:
                sanitized_smiles, _ = OINSanitizer.generate_robust_smiles(frag_mol, local_binding_indices)
            else:
                 sanitized_smiles = f"[{mol.GetAtomWithIdx(metal_idx).GetSymbol()}]"

            # Determine Map Number mapping to get `local_idx` matching SMILES order
            # (Need to run MolFromSmiles on sanitized string)
            smiles_mol = Chem.MolFromSmiles(sanitized_smiles, sanitize=False)
            frag_to_smiles_idx = {}
            if smiles_mol:
                for s_atom in smiles_mol.GetAtoms():
                    map_num = s_atom.GetAtomMapNum()
                    if map_num > 0:
                        l_idx_in_frag = map_num - 1
                        frag_to_smiles_idx[l_idx_in_frag] = s_atom.GetIdx()
            
            # Update binding_atoms to include local_idx (SMILES idx)
            final_binding_atoms = []
            for g_idx, m, coords in frag_binding_atoms:
                if g_idx in old_to_new:
                    l_idx = old_to_new[g_idx]
                    s_idx = frag_to_smiles_idx.get(l_idx, 0)
                    final_binding_atoms.append((g_idx, m, coords, s_idx))
            
            # Cleanup map numbers in SMILES
            clean_smiles = re.sub(r':\d+\]', ']', sanitized_smiles)
            item['smiles'] = clean_smiles
            item['binding_atoms'] = final_binding_atoms

        # 4. Canonical Sorting (Mass-First)
        def get_input_order_key(item):
            if item['is_metal']: return -1
            valid_indices = [atom_map_to_xyz.get(idx, idx) for idx in item['indices']]
            return min(valid_indices) if valid_indices else float('inf')

        def get_canonical_sort_key(item):
            if item['is_metal']: return (-float('inf'),)
            mw = item['mass']
            b_mass = item['binding_mass']
            smiles = item['smiles']
            # Tie-break with input order
            return (-mw, -b_mass, smiles, get_input_order_key(item))
            
        fragments_data.sort(key=get_canonical_sort_key)
        
        # 5. Geometry Detection & Alignment
        aligner = OINDiscreteAligner(0, fragments_data) # Metal is at 0 after sort
        geometry_string_raw = aligner.generate_canonical_vectors()
        
        # 6. Re-Sort and Serialize
        if "w:NON" in geometry_string_raw or "error" in geometry_string_raw:
             # Fallback
             full_smiles = ".".join([f['smiles'] for f in fragments_data])
             sidecar_oin = f"{full_smiles} |{geometry_string_raw}|"
        else:
             # Parse result to re-sort by Slot
             parts = geometry_string_raw.split("|")
             geo_tag = parts[0]
             w_content = parts[1][2:] # Remove w:
             
             rank_to_slots = {}
             pair_data = [] # (Rank, LocalIdx, Slot, HeadingChar)
             
             if w_content and w_content != "NON":
                 for entry in w_content.split(";"):
                     if ":" not in entry: continue
                     left, slot_str = entry.split(":")
                     
                     heading_char = ""
                     for char in ['^', '>', '<']:
                         if char in slot_str:
                             heading_char = char
                             break
                     slot_str_clean = slot_str.replace('^', '').replace('>', '').replace('<', '')
                     slot = int(slot_str_clean)
                     
                     if "." in left:
                         r_str, l_str = left.split(".")
                         rank = int(float(r_str))
                         l_idx = int(l_str)
                     else:
                         rank = int(float(left))
                         l_idx = 0
                     
                     if rank not in rank_to_slots: rank_to_slots[rank] = []
                     rank_to_slots[rank].append(slot)
                     pair_data.append((rank, l_idx, slot, heading_char))

             for r, frag in enumerate(fragments_data):
                 if frag['is_metal']: frag['_sort_slot'] = -1
                 elif r in rank_to_slots:
                     frag['_sort_slot'] = min(rank_to_slots[r])
                 else:
                     frag['_sort_slot'] = float('inf')
                 frag['_orig_rank'] = r
            
             # Re-Sort
             fragments_data.sort(key=lambda x: (x['_sort_slot'], get_input_order_key(x)))
             
             # Map Old -> New Ranks
             old_to_new_rank = {frag['_orig_rank']: new_r for new_r, frag in enumerate(fragments_data)}
             
             # Reconstruct Tag
             new_pair_data = []
             for old_r, l_idx, slot, hd in pair_data:
                 if old_r in old_to_new_rank:
                     new_r = old_to_new_rank[old_r]
                     new_pair_data.append((new_r, l_idx, slot, hd))
             
             new_pair_data.sort(key=lambda x: (x[0], x[1]))
             new_w_parts = [f"{nr}.{li}:{sl}{hd}" for nr, li, sl, hd in new_pair_data]
             new_geometry_string = f"{geo_tag}|w:{';'.join(new_w_parts)}"
             
             full_smiles = ".".join([f['smiles'] for f in fragments_data])
             sidecar_oin = f"{full_smiles} |{new_geometry_string}|"

        # Final Inline Conversion
        if OINInlineHandler:
            return OINInlineHandler.generate_inline_string(sidecar_oin)
        return sidecar_oin
