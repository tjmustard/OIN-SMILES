import numpy as np
from scipy.spatial.transform import Rotation
from collections import defaultdict
import itertools
import logging
import sys
import warnings

# Re-export sanitizer for compatibility
from ..oin.sanitizer import OINSanitizer

logger = logging.getLogger(__name__)

def _align_to_pai(tmc_mol, xyz_coords, metal_idx):
    """
    Canonicalizes the orientation of the molecule:
    1. Translates so the metal is at (0,0,0).
    2. Rotates so the Principal Axes of Inertia (PAI) align with the Cartesian axes.
       - Highest Moment of Inertia -> Z
       - Lowest Moment of Inertia -> X
       - Enforce Right-Handed System
    """
    coords = np.array(xyz_coords)
    masses = np.array([a.GetMass() for a in tmc_mol.GetAtoms()])
    
    # 1. Translate Metal to Origin
    metal_pos = coords[metal_idx]
    coords -= metal_pos
    
    # 2. Calculate Inertia Tensor relative to Origin (Metal)
    I = np.zeros((3, 3))
    for i in range(len(coords)):
        m = masses[i]
        pos = coords[i]
        sq_norm = np.dot(pos, pos)
        
        # Diagonal elements
        I[0, 0] += m * (sq_norm - pos[0]*pos[0])
        I[1, 1] += m * (sq_norm - pos[1]*pos[1])
        I[2, 2] += m * (sq_norm - pos[2]*pos[2])
        
        # Off-diagonal elements (symmetric, negative product)
        I[0, 1] -= m * (pos[0]*pos[1])
        I[0, 2] -= m * (pos[0]*pos[2])
        I[1, 2] -= m * (pos[1]*pos[2])
        
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]
    
    # 3. Diagonalize
    evals, evecs = np.linalg.eigh(I)
    
    x_axis = evecs[:, 0] # v1 (Lowest)
    z_axis = evecs[:, 2] # v3 (Highest)
    
    # Enforce Right-Handed System: Y = Z x X
    y_axis = np.cross(z_axis, x_axis)
    
    # Construct Rotation Matrix (Rows are new basis vectors)
    R = np.vstack([x_axis, y_axis, z_axis])
    
    # Apply Rotation
    new_coords = coords @ R.T
    
    # 4. Handle Degeneracy (Axial Rotation) & Sign Ambiguity
    # Determine Pivot Atom to fix X-axis rotation
    dists_sq = np.sum(new_coords**2, axis=1)
    max_dist_sq = np.max(dists_sq)
    tolerance = 1e-5
    candidates = np.where(dists_sq >= max_dist_sq - tolerance)[0]
    pivot_idx = np.min(candidates)
    
    pivot_pos = new_coords[pivot_idx]
    angle = np.arctan2(pivot_pos[1], pivot_pos[0])
    
    c = np.cos(-angle)
    s = np.sin(-angle)
    
    R_pivot = np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])
    
    canonical_coords = new_coords @ R_pivot.T
    
    # 5. Handle Z-Axis Sign Ambiguity
    z_moment_idx = 0.0
    for i in range(len(canonical_coords)):
        z_moment_idx += canonical_coords[i][2] * (i + 1)**3
        
    if z_moment_idx < 0:
        canonical_coords[:, 1] *= -1
        canonical_coords[:, 2] *= -1
    
    return canonical_coords.tolist()

def normalize_template(arr):
    return arr / np.linalg.norm(arr, axis=1)[:, None]

def normalize_cols(arr):
    return arr / np.linalg.norm(arr, axis=1)[:, None]

TEMPLATE_SPECS = {
    'LIN': {
        0: {'pos': [0, 0, 1],  'ref': [1, 0, 0]},
        1: {'pos': [0, 0, -1], 'ref': [1, 0, 0]}
    },
    'TPL': {
        0: {'pos': [0, 1, 0],            'ref': [0, 0, 1]},
        1: {'pos': [0.8660254, -0.5, 0], 'ref': [0, 0, 1]},
        2: {'pos': [-0.8660254, -0.5, 0],'ref': [0, 0, 1]}
    },
    'SPL': {
        0: {'pos': [1, 0, 0],  'ref': [0, 0, 1]},
        1: {'pos': [0, 1, 0],  'ref': [0, 0, 1]},
        2: {'pos': [-1, 0, 0], 'ref': [0, 0, 1]},
        3: {'pos': [0, -1, 0], 'ref': [0, 0, 1]}
    },
    'TET': {
        0: {'pos': [1, 1, 1],    'ref': [-1, 1, 0]},
        1: {'pos': [1, -1, -1],  'ref': [0, 1, -1]},
        2: {'pos': [-1, 1, -1],  'ref': [1, 1, 0]},
        3: {'pos': [-1, -1, 1],  'ref': [1, 0, 1]}
    },
    'TPY': {
        0: {'pos': [0, 0, 1],             'ref': [1, 0, 0]},
        1: {'pos': [0, 1, 0],             'ref': [0, 0, 1]},
        2: {'pos': [0.8660254, -0.5, 0],  'ref': [0, 0, 1]},
        3: {'pos': [-0.8660254, -0.5, 0], 'ref': [0, 0, 1]}
    },
    'TBP': {
        0: {'pos': [0, 0, 1],             'ref': [1, 0, 0]},
        1: {'pos': [0, 0, -1],            'ref': [1, 0, 0]},
        2: {'pos': [0, 1, 0],             'ref': [0, 0, 1]},
        3: {'pos': [0.8660254, -0.5, 0],  'ref': [0, 0, 1]},
        4: {'pos': [-0.8660254, -0.5, 0], 'ref': [0, 0, 1]}
    },
    'SPY': {
        0: {'pos': [0, 0, 1],  'ref': [1, 0, 0]},
        1: {'pos': [1, 0, 0],  'ref': [0, 0, 1]},
        2: {'pos': [-1, 0, 0], 'ref': [0, 0, 1]},
        3: {'pos': [0, 1, 0],  'ref': [0, 0, 1]},
        4: {'pos': [0, -1, 0], 'ref': [0, 0, 1]}
    },
    'OCT': {
        0: {'pos': [0, 0, 1],  'ref': [1, 0, 0]},
        1: {'pos': [0, 0, -1], 'ref': [1, 0, 0]},
        2: {'pos': [1, 0, 0],  'ref': [0, 0, 1]},
        3: {'pos': [-1, 0, 0], 'ref': [0, 0, 1]},
        4: {'pos': [0, 1, 0],  'ref': [0, 0, 1]},
        5: {'pos': [0, -1, 0], 'ref': [0, 0, 1]}
    },
    'PBP': {
        0: {'pos': [0, 0, 1],             'ref': [1, 0, 0]},
        1: {'pos': [0, 0, -1],            'ref': [1, 0, 0]},
        2: {'pos': [1, 0, 0],             'ref': [0, 0, 1]},
        3: {'pos': [0.309017, 0.951057, 0], 'ref': [0, 0, 1]},
        4: {'pos': [-0.809017, 0.587785, 0], 'ref': [0, 0, 1]},
        5: {'pos': [-0.809017, -0.587785, 0], 'ref': [0, 0, 1]},
        6: {'pos': [0.309017, -0.951057, 0], 'ref': [0, 0, 1]}
    }
}

SYMMETRIC_LIGANDS = {
    "C=C",
    "[CH2]=[CH2]",
    "c1cccc1",
    "C1=C-C=C-[CH-]1",
    "[cH]1[cH][cH][cH][cH]1",
    "c1ccccc1",
    "[cH]1[cH][cH][cH][cH][cH]1",
    "c1cccccc1",
    "[cH]1[cH][cH][cH][cH][cH][cH]1",
}

# Generate Legacy TEMPLATES
TEMPLATES = {}
for geo, specs in TEMPLATE_SPECS.items():
    sorted_slots = sorted(specs.keys())
    pos_vecs = np.array([specs[idx]['pos'] for idx in sorted_slots], dtype=float)
    if geo == 'TET':
        pos_vecs = normalize_cols(pos_vecs)
        for i, vec in enumerate(pos_vecs):
            TEMPLATE_SPECS[geo][sorted_slots[i]]['pos'] = vec.tolist()
    TEMPLATES[geo] = pos_vecs

class OINDiscreteAligner:
    def __init__(self, metal_idx, ligands):
        """
        ligands: List of dicts.
        Each dict MUST contain 'smiles' generated by OINSanitizer.
        And 'binding_atoms' list of tuples/lists: [(global_idx, mass, coords, local_idx)]
        """
        self.metal_idx = metal_idx
        self.ligands = ligands
        
    def generate_canonical_vectors(self):
        # 1. Haptic Reduction
        virtual_atoms = self._reduce_hapticity()

        # 2. Competitive Geometry Detection
        n_eff = len(virtual_atoms)
        if n_eff < 2: 
             tmpl_name = 'LIN' 
             tmpl_vectors = TEMPLATES['LIN']
             mapping = [None, None] 
             R_mat = None
        else:
            best_res = self._find_best_geometry_match(n_eff, virtual_atoms)
            if best_res:
                tmpl_name, tmpl_vectors, mapping, R_mat = best_res
            else:
                 return "g:NON|w:NON"

        # 3. Canonicalize
        canonical_str = self._permute_and_serialize(mapping, tmpl_vectors, geometry_name=tmpl_name, alignment_rotation=R_mat)
        
        return f"g:{tmpl_name}|w:{canonical_str}"

    def _reduce_hapticity(self):
        virtual_atoms = []
        for i, lig in enumerate(self.ligands):
            if i == self.metal_idx:
                continue
            if not lig.get('binding_atoms'):
                continue
            
            first_binding_atom_mass = lig['binding_atoms'][0][1]
            base_sort_key = (i, first_binding_atom_mass, lig['smiles']) 
            
            binding_coords = np.array([ba[2] for ba in lig['binding_atoms']])
            
            metal_frag = self.ligands[self.metal_idx]
            metal_origin = metal_frag.get('metal_coords')
            if metal_origin is None:
                metal_origin = np.array([0.0, 0.0, 0.0])
                
            zone_a_info = lig['binding_atoms'] 
            n_b = len(binding_coords)
            
            groups = []
            visited = set()
            for j in range(n_b):
                if j in visited: continue
                stack = [j]
                component = []
                while stack:
                    curr = stack.pop()
                    if curr in visited: continue
                    visited.add(curr)
                    component.append(curr)
                    for k in range(n_b):
                        if k in visited: continue
                        if np.linalg.norm(binding_coords[curr] - binding_coords[k]) < 1.6:
                            stack.append(k)
                groups.append(component)
            
            for grp in groups:
                grp.sort(key=lambda k: zone_a_info[k][3])
                grp_coords = binding_coords[grp]
                if len(grp_coords) == 0:
                    continue
                centroid = np.mean(grp_coords, axis=0)
                rep_idx = zone_a_info[grp[0]][3]
                constituent_indices = sorted([zone_a_info[k][3] for k in grp])
                
                virtual_atoms.append({
                    'rank': i, 
                    'local_idx': rep_idx,
                    'constituent_indices': constituent_indices,
                    'key': base_sort_key,
                    'coords': centroid - metal_origin,
                    'group_coords': grp_coords - metal_origin,
                    'chem_id': (first_binding_atom_mass, lig['smiles'])
                })
        
        return virtual_atoms

    def _find_best_geometry_match(self, n, virtual_atoms):
        candidates = []
        if n == 2: candidates = ['LIN']
        elif n == 3: candidates = ['TPL']
        elif n == 4: candidates = ['SPL', 'TET', 'TPY']
        elif n == 5: candidates = ['TBP', 'SPY']
        elif n == 6: candidates = ['OCT']
        elif n == 7: candidates = ['PBP']
        else: 
            if n > 7: candidates = ['OCT'] 
            else: candidates = ['LIN']
        
        min_rmsd = float('inf')
        best_result = None
        
        for name in candidates:
            vectors = TEMPLATES.get(name)
            if vectors is None: continue
            if n > len(vectors): continue

            mapping, rmsd, R_mat = self._map_to_template(virtual_atoms, vectors)
            
            if mapping is not None and rmsd < min_rmsd:
                min_rmsd = rmsd
                best_result = (name, vectors, mapping, R_mat)
        
        return best_result
        
    def _map_to_template(self, virtual_atoms, template_vectors):
        n_atoms = len(virtual_atoms)
        n_slots = len(template_vectors)
        
        if n_atoms == 0: 
            return [None]*n_slots, 0.0, Rotation.from_matrix(np.eye(3))

        input_vecs = np.array([a['coords'] for a in virtual_atoms])
        input_norms = input_vecs / (np.linalg.norm(input_vecs, axis=1)[:,None] + 1e-9)
        
        best_rmsd = float('inf')
        best_mapping = None
        best_R = Rotation.from_matrix(np.eye(3))

        perm_iterator = itertools.permutations(range(n_slots), n_atoms)
        
        for slot_indices in perm_iterator:
             target_vecs = template_vectors[list(slot_indices)]
             
             try:
                 with warnings.catch_warnings():
                     warnings.filterwarnings("ignore", message="Optimal rotation is not uniquely or poorly defined")
                     R, rmsd = Rotation.align_vectors(target_vecs, input_norms)
             except:
                 continue
                 
             if rmsd < best_rmsd:
                 best_rmsd = rmsd
                 current_mapping = [None] * n_slots
                 for atom_idx, slot_idx in enumerate(slot_indices):
                     current_mapping[slot_idx] = virtual_atoms[atom_idx]
                     
                 best_mapping = current_mapping
                 best_R = R
        
        return best_mapping, best_rmsd, best_R


    def _permute_and_serialize(self, slot_assignment, tmpl_vectors, geometry_name=None, alignment_rotation=None):
        symmetries = self._brute_force_symmetries(tmpl_vectors)
        best_sequence = None
        best_final_map = None 
        
        if not slot_assignment: return ""

        for perm in symmetries:
            current_view_map = []
            
            for old_slot_idx, atom in enumerate(slot_assignment):
                if atom is None: continue
                new_slot_idx = perm[old_slot_idx]
                ideal_vec = tmpl_vectors[new_slot_idx]
                v_clean = []
                for val in ideal_vec:
                    if abs(val) < 1e-9: val = 0.0
                    v_clean.append(val)
                ideal_vec = tuple(v_clean)
                
                current_view_map.append({
                    'rank': atom['rank'],
                    'local_idx': atom['local_idx'], 
                    'constituent_indices': atom.get('constituent_indices', [atom['local_idx']]), 
                    'group_coords': atom.get('group_coords'), 
                    'chem_id': atom['chem_id'],
                    'vec': ideal_vec,
                    'slot': new_slot_idx
                })
            
            if not current_view_map: continue

            grouped = defaultdict(list)
            for item in current_view_map:
                grouped[item['chem_id']].append(item)
            
            final_sorted_view = []
            for chem_id, items in grouped.items():
                frag_groups = defaultdict(list)
                for it in items:
                    frag_groups[it['rank']].append(it)
                
                target_ranks = sorted(list(frag_groups.keys()))
                
                available_sets = []
                for rank in target_ranks:
                    f_items = frag_groups[rank]
                    f_items.sort(key=lambda x: x['local_idx'])
                    vec_set = tuple([x['vec'] for x in f_items])
                    available_sets.append({
                        'vec_set': vec_set,
                        'items': f_items
                    })
                
                available_sets.sort(key=lambda x: x['vec_set'], reverse=True)
                
                for i, target_rank in enumerate(target_ranks):
                    assigned_set = available_sets[i]
                    for atom_data in assigned_set['items']:
                        final_sorted_view.append({
                            'rank': target_rank,
                            'local_idx': atom_data['local_idx'],
                            'constituent_indices': atom_data['constituent_indices'],
                            'group_coords': atom_data['group_coords'],
                            'vec': atom_data['vec'],
                            'slot': atom_data['slot'],
                            'chem_id': chem_id
                        })

            final_sorted_view.sort(key=lambda x: (x['rank'], x['local_idx']))
            current_sequence = [x['vec'] for x in final_sorted_view]
            
            if best_sequence is None or current_sequence > best_sequence:
                best_sequence = current_sequence
                best_final_map = final_sorted_view
        
        if not best_final_map: return "error"

        heading_local_indices = set()
        
        if geometry_name and alignment_rotation is not None and geometry_name in TEMPLATE_SPECS:
            by_rank = defaultdict(list)
            for x in best_final_map:
                by_rank[x['rank']].append(x)
            
            template_spec = TEMPLATE_SPECS[geometry_name]
            
            for rank, items in by_rank.items():
                first_item = items[0]
                slot_idx = first_item['slot']
                
                if slot_idx not in template_spec or 'ref' not in template_spec[slot_idx]:
                    continue
                
                ref_vec = np.array(template_spec[slot_idx]['ref'])
                grp_coords = first_item.get('group_coords') 
                if grp_coords is None: continue
                centroid = np.mean(grp_coords, axis=0)
                
                if len(grp_coords) < 2: continue
                
                best_dot = -float('inf')
                best_idx = -1
                ordered_indices = first_item['constituent_indices']
                
                for k, coord in enumerate(grp_coords):
                    v_mol = coord - centroid
                    v_tmpl = alignment_rotation.apply(v_mol)
                    norm = np.linalg.norm(v_tmpl)
                    if norm > 1e-6:
                        v_tmpl_n = v_tmpl / norm
                        dot = np.dot(v_tmpl_n, ref_vec)
                        if dot > best_dot:
                            best_dot = dot
                            best_idx = ordered_indices[k]
                
                if best_idx != -1:
                    heading_local_indices.add((rank, best_idx))
            
            for rank, items in by_rank.items():
                first_item = items[0]
                smiles = first_item['chem_id'][1]
                if smiles in SYMMETRIC_LIGANDS:
                    ordered_indices = sorted(first_item['constituent_indices'])
                    forced_idx = ordered_indices[0]
                    to_remove = [idx for r, idx in heading_local_indices if r == rank]
                    for idx in to_remove:
                        heading_local_indices.remove((rank, idx))
                    heading_local_indices.add((rank, forced_idx))

        parts = []
        for x in best_final_map:
            slot = x['slot']
            rank = x['rank']
            indices = x.get('constituent_indices', [x['local_idx']])
            
            for idx in indices:
                tag = f"{rank}.{idx}:{slot}"
                if (rank, idx) in heading_local_indices:
                     direction_char = ">"
                     if geometry_name in TEMPLATE_SPECS and slot in TEMPLATE_SPECS[geometry_name]:
                         slot_def = TEMPLATE_SPECS[geometry_name][slot]
                         c_indices = sorted(x.get('constituent_indices', [idx]))
                         direction_char = self._determine_winding(
                             grp_coords=x.get('group_coords'),
                             star_idx=idx,
                             constituent_indices=c_indices,
                             slot_z=np.array(slot_def['pos']),
                             slot_x_ref=np.array(slot_def['ref']),
                             alignment_rotation=alignment_rotation
                         )
                     tag += direction_char
                parts.append(tag)
        
        return ";".join(parts)

    def _determine_winding(self, grp_coords, star_idx, constituent_indices, slot_z, slot_x_ref, alignment_rotation=None):
        n = len(constituent_indices)
        if n < 3: return ">"
        
        try:
            list_idx = constituent_indices.index(star_idx)
        except ValueError:
            return ">"
            
        coord_star = grp_coords[list_idx]
        next_list_idx = (list_idx + 1) % n
        coord_next = grp_coords[next_list_idx]
        centroid = np.mean(grp_coords, axis=0)
        
        v_star_mol = coord_star - centroid
        v_next_mol = coord_next - centroid
        
        if alignment_rotation:
            v_star = alignment_rotation.apply(v_star_mol)
            v_next = alignment_rotation.apply(v_next_mol)
        else:
            v_star = v_star_mol
            v_next = v_next_mol
            
        winding_normal = np.cross(v_star, v_next)
        dot = np.dot(winding_normal, slot_z)
        
        return ">" if dot >= 0 else "<"

    def _brute_force_symmetries(self, vectors):
        n = len(vectors)
        valid = set()
        steps = [0, 90, 120, 180, 240, 270]
        for rx, ry, rz in itertools.product(steps, repeat=3):
            R = Rotation.from_euler('xyz', [rx, ry, rz], degrees=True)
            rot = R.apply(vectors)
            perm = [-1]*n
            matches = 0
            for i in range(n):
                dists = np.linalg.norm(vectors - rot[i], axis=1)
                best = np.argmin(dists)
                if dists[best] < 0.1:
                    perm[i] = best
                    matches += 1
            if matches == n: 
                valid.add(tuple(perm))
        return sorted(list(valid))
