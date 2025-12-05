
import numpy as np

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm

def get_rotation_matrix_to_align_z(v_target):
    """
    Returns a rotation matrix that aligns v_target to the global Z-axis (0,0,1).
    """
    v_target = normalize(v_target)
    z_axis = np.array([0.0, 0.0, 1.0])

    # Check for already aligned
    if np.allclose(v_target, z_axis):
        return np.eye(3)

    # Check for perfectly opposite (180 deg)
    if np.allclose(v_target, -z_axis):
        return np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])

    # Axis of rotation is cross product of target and Z
    axis = np.cross(v_target, z_axis)
    axis = normalize(axis)
    # Clamp dot product to avoid acos domain errors
    angle = np.arccos(np.clip(np.dot(v_target, z_axis), -1.0, 1.0))

    # Rodrigues' rotation formula
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])

    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)
    return R

def get_rotation_to_fix_x_plane(coords, p2_index):
    """
    Rotates around Z-axis so that the atom at p2_index has y=0 and x>0.
    """
    p2 = coords[p2_index]
    # Project p2 onto XY plane (z is invariant)
    x, y = p2[0], p2[1]

    current_angle = np.arctan2(y, x)
    rotation_angle = -current_angle

    c, s = np.cos(rotation_angle), np.sin(rotation_angle)
    R_z = np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])

    return R_z

class OINCanonicalAligner:
    def __init__(self, ligands):
        """
        ligands: List of objects containing:
                 - 'smiles': Canonical SMILES of fragment
                 - 'mass': Molecular weight
                 - 'binding_atoms': list of (original_index, atomic_mass, coords_xyz)
        """
        self.ligands = ligands

    def get_best_alignment(self):
        """
        V1.6 Standard Entry Point: Double Exhaustion Strategy.
        Returns:
            R_final: The rotation matrix (3x3) for the canonical alignment.
        """
        # --- Step 1: Mass-First Sort ---
        # 1. Fragment Mass (Desc)
        # 2. Heaviest Binding Atom Mass (Desc)
        # 3. SMILES String (Asc)
        ranked_ligands = sorted(
            self.ligands,
            key=lambda x: (-x['mass'], -x['binding_atoms'][0][1], x['smiles'])
        )

        # --- Step 2: Identify Z-Axis Candidates (P1) ---
        # Returns list of (LigandObject, Coords) for all atoms tying for Rank 0
        p1_candidates = self._get_z_anchor_candidates(ranked_ligands)
        
        candidates_results = []

        # --- Step 3: Outer Loop (Iterate through all valid P1s) ---
        for p1_ligand, p1_atom_coords in p1_candidates:
            
            # --- Step 4: Identify X-Axis Candidates (P2) ---
            # V1.7 CRITICAL CHANGE: Scans ALL ranks for valid anchors.
            p2_candidates = self._get_global_p2_candidates(p1_atom_coords, ranked_ligands)
            
            if not p2_candidates:
                # Linear molecule case: P2 is effectively None
                p2_candidates = [None]

            # --- Step 5: Inner Loop (Iterate through all valid P2s) ---
            for p2_atom_coords in p2_candidates:
                
                # Phase C: Alignment
                aligned_coords_map, R_final = self._align_structure(
                    p1_atom_coords, 
                    p2_atom_coords, 
                    ranked_ligands
                )
                
                # Phase D: Serialization (Sanitized)
                oin_string = self._serialize_vectors(aligned_coords_map)
                candidates_results.append((oin_string, R_final))

        # --- Step 6: Final Selection ---
        if not candidates_results:
             raise ValueError("No valid binding atoms found.")
               
        # Lexicographical sort picks the canonical orientation
        candidates_results.sort(key=lambda x: x[0])
        return candidates_results[0][1]

    def _get_z_anchor_candidates(self, ranked_ligands):
        """
        Identifies all atoms that tie for the highest priority rank.
        """
        candidates = []
        best_ligand = ranked_ligands[0]
        
        # Max binding mass within the best ligand
        best_atom_mass = -1
        for batom in best_ligand['binding_atoms']:
            if batom[1] > best_atom_mass:
                best_atom_mass = batom[1]
                
        # Search all ligands for matches
        for lig in ranked_ligands:
            is_tied_ligand = (
                lig['mass'] == best_ligand['mass'] and 
                lig['smiles'] == best_ligand['smiles']
            )
            
            if not is_tied_ligand:
                break 
            
            for batom in lig['binding_atoms']:
                if batom[1] == best_atom_mass:
                    candidates.append((lig, batom[2]))
                    
        return candidates

    def _get_global_p2_candidates(self, p1_coords, ranked_ligands):
        """
        V1.7 LOGIC: Returns ALL non-linear atoms from ALL ranks.
        Does not stop at the first valid rank.
        """
        p1_norm = normalize(p1_coords)
        candidates = []

        # Iterate through EVERY ligand in the complex
        for lig in ranked_ligands:
            for batom in lig['binding_atoms']:
                p2_coords = batom[2]
                
                # Identity Check
                if np.array_equal(p2_coords, p1_coords): continue
                
                # Linearity Check
                p2_norm = normalize(p2_coords)
                dot = np.dot(p1_norm, p2_norm)
                
                # If not collinear (0 or 180 deg), it is a valid candidate
                if abs(dot) < 0.99:
                     candidates.append(p2_coords)

        return candidates

    def _align_structure(self, p1_coords, p2_coords, all_ligands):
        """
        Performs the matrix rotations.
        """
        # 1. Rotate P1 to +Z
        R1 = get_rotation_matrix_to_align_z(p1_coords)
        
        if p2_coords is not None:
            # 2. Rotate P2 intermediate to +X plane
            p2_intermediate = np.dot(R1, p2_coords)
            R2 = get_rotation_to_fix_x_plane([p2_intermediate], 0)
            R_final = np.dot(R2, R1)
        else:
            R_final = R1

        # Apply to all atoms
        aligned_map = []
        for i, lig in enumerate(all_ligands):
            for batom in lig['binding_atoms']:
                vec = batom[2]
                new_vec = np.dot(R_final, vec)
                aligned_map.append({
                    'ligand_rank': i,
                    'vec': new_vec
                })
            
        return aligned_map, R_final

    def _serialize_vectors(self, aligned_map):
        """
        V1.6: Zero Sanitation Included.
        """
        parts = []
        aligned_map.sort(key=lambda x: (x['ligand_rank'], x['vec'].tolist()))
        
        for item in aligned_map:
            v = item['vec']
            
            # --- SANITATION ---
            v_clean = []
            for val in v:
                if abs(val) < 1e-9:
                    val = 0.0
                v_clean.append(val)
            v = np.array(v_clean)
            # ------------------
            
            s = f"{item['ligand_rank']}:{v[0]:.3f},{v[1]:.3f},{v[2]:.3f}"
            parts.append(s)
            
        return ";".join(parts)
