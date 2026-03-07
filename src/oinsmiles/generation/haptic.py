import numpy as np
from typing import List, Dict, Any, Optional

def normalize(v):
    norm = np.linalg.norm(v)
    return v / norm if norm > 0 else v

# V3.7 Templates with Reference Vectors
# Structure: GEO -> Slot Idx -> {pos, ref}
TEMPLATES = {
    'LIN': {
        0: {'pos': np.array([0,0,1]), 'ref': np.array([1,0,0])},
        1: {'pos': np.array([0,0,-1]), 'ref': np.array([1,0,0])}
    },
    'TPL': { # Trigonal Planar (xy plane)
        0: {'pos': np.array([0,1,0]), 'ref': np.array([0,0,1])},
        1: {'pos': np.array([0.8660254,-0.5,0]), 'ref': np.array([0,0,1])},
        2: {'pos': np.array([-0.8660254,-0.5,0]), 'ref': np.array([0,0,1])}
    },
    'SPL': { # Square Planar (xy plane)
        0: {'pos': np.array([1,0,0]), 'ref': np.array([0,0,1])},
        1: {'pos': np.array([0,1,0]), 'ref': np.array([0,0,1])},
        2: {'pos': np.array([-1,0,0]), 'ref': np.array([0,0,1])},
        3: {'pos': np.array([0,-1,0]), 'ref': np.array([0,0,1])}
    },
    'OCT': {
        0: {'pos': np.array([0,0,1]), 'ref': np.array([1,0,0])}, # Axial Top
        1: {'pos': np.array([0,0,-1]), 'ref': np.array([1,0,0])}, # Axial Bottom
        2: {'pos': np.array([1,0,0]), 'ref': np.array([0,0,1])}, # Eq
        3: {'pos': np.array([-1,0,0]), 'ref': np.array([0,0,1])}, # Eq
        4: {'pos': np.array([0,1,0]), 'ref': np.array([0,0,1])}, # Eq
        5: {'pos': np.array([0,-1,0]), 'ref': np.array([0,0,1])}  # Eq
    },
    'SPY': {
        0: {'pos': np.array([0,0,1]), 'ref': np.array([1,0,0])}, # Axial
        1: {'pos': np.array([1,0,0]), 'ref': np.array([0,0,1])}, # Eq
        2: {'pos': np.array([-1,0,0]), 'ref': np.array([0,0,1])}, # Eq
        3: {'pos': np.array([0,1,0]), 'ref': np.array([0,0,1])}, # Eq
        4: {'pos': np.array([0,-1,0]), 'ref': np.array([0,0,1])}  # Eq
    },
    'TET': {
        0: {'pos': np.array([1,1,1]), 'ref': np.array([0,0,1])},
        1: {'pos': np.array([1,-1,-1]), 'ref': np.array([0,0,1])},
        2: {'pos': np.array([-1,1,-1]), 'ref': np.array([0,0,1])},
        3: {'pos': np.array([-1,-1,1]), 'ref': np.array([0,0,1])}
    },
    'TPY': {
        0: {'pos': np.array([0,0,1]), 'ref': np.array([1,0,0])},
        1: {'pos': np.array([0,1,0]), 'ref': np.array([0,0,1])},
        2: {'pos': np.array([0.8660254,-0.5,0]), 'ref': np.array([0,0,1])},
        3: {'pos': np.array([-0.8660254,-0.5,0]), 'ref': np.array([0,0,1])}
    },
    'TBP': {
        0: {'pos': np.array([0,0,1]), 'ref': np.array([1,0,0])},
        1: {'pos': np.array([0,0,-1]), 'ref': np.array([1,0,0])},
        2: {'pos': np.array([0,1,0]), 'ref': np.array([0,0,1])},
        3: {'pos': np.array([0.8660254,-0.5,0]), 'ref': np.array([0,0,1])},
        4: {'pos': np.array([-0.8660254,-0.5,0]), 'ref': np.array([0,0,1])}
    },
    'PBP': {
        0: {'pos': np.array([0,0,1]), 'ref': np.array([1,0,0])}, # Axial
        1: {'pos': np.array([0,0,-1]), 'ref': np.array([1,0,0])}, # Axial
        2: {'pos': np.array([1,0,0]), 'ref': np.array([0,0,1])},
        3: {'pos': np.array([0.30901699, 0.95105652, 0]), 'ref': np.array([0,0,1])},
        4: {'pos': np.array([-0.80901699, 0.58778525, 0]), 'ref': np.array([0,0,1])},
        5: {'pos': np.array([-0.80901699, -0.58778525, 0]), 'ref': np.array([0,0,1])},
        6: {'pos': np.array([0.30901699, -0.95105652, 0]), 'ref': np.array([0,0,1])}
    }
}

# Fallback for other geometries:
# If not distinct, assume POS from oin_parser definitions and REF = [0,0,1] or [1,0,0]
# We can inject them on demand if needed or expand this table later.

class HapticTransformer:
    
    @staticmethod
    def get_step_angle(n: int) -> float:
        if n == 2: return 180.0
        if n == 3: return 72.0 # Subset of Pentagon
        if n == 4: return 60.0 # Subset of Hexagon
        if n == 5: return 72.0
        if n == 6: return 60.0
        if n == 7: return 360.0 / 7.0
        if n >= 8: return 360.0 / float(n) # General Ring
        return 0.0

    @staticmethod
    def transform_vectors(n_haptic: int, slot_z: np.ndarray, slot_x_ref: np.ndarray, cone_spread: float = 0.2) -> List[np.ndarray]:
        """
        Generates N vectors. 
        Vector 0 is always aligned with slot_x_ref (projected).
        """
        if n_haptic <= 1:
            return [normalize(slot_z)]

        # 1. Basis Frame
        Z = normalize(slot_z)
        X = normalize(slot_x_ref)
        # Orthogonalize X w.r.t Z
        # Remove component of X parallel to Z
        X = X - (np.dot(X, Z) * Z)
        
        norm_x = np.linalg.norm(X)
        if norm_x < 1e-6:
            # Degenerate X (X was parallel to Z). Pick arbitrary Y.
            # Try Global X
            X = np.array([1,0,0], dtype=float) - (np.dot(np.array([1,0,0]), Z) * Z)
            norm_x = np.linalg.norm(X)
            if norm_x < 1e-6:
                # Try Global Y
                X = np.array([0,1,0], dtype=float) - (np.dot(np.array([0,1,0]), Z) * Z)
                
        X = normalize(X)
        Y = np.cross(Z, X)
        
        vectors = []
        step = HapticTransformer.get_step_angle(n_haptic)
        
        for k in range(n_haptic):
            theta_deg = k * step
            theta_rad = np.radians(theta_deg)
            
            # Planar components
            # cos(0) = 1 -> Aligns with X (North Star)
            vx = np.cos(theta_rad)
            vy = np.sin(theta_rad)
            
            # Construct Vector: Z + spread*(x*X + y*Y)
            vec = Z + cone_spread * (vx * X + vy * Y)
            vectors.append(normalize(vec))
            
        return vectors

class HapticResolver:
    """
    Helper to resolve vectors for a fragment based on Heading/Direction
    """
    @staticmethod
    def resolve_fragment_vectors(
        geometry: str,
        slot_idx: int,
        n_atoms: int,
        heading_atom_in_frag_idx: int, # The atom marked with {0} (or first if implicit)
        direction: int, # 1 or -1
        cone_spread: float = 0.2
    ) -> List[Dict[str, Any]]:
        """
        Returns list of dicts: {'atom_idx_in_frag': int, 'vector': np.array}
        """
        
        # 1. Get Template Slot
        geo_def = TEMPLATES.get(geometry)
        if not geo_def or slot_idx not in geo_def:
             # Fallback: Just return Z for all atoms (collapse)
             # Or try to construct on fly
             return []
             
        slot_def = geo_def[slot_idx]
        slot_z = slot_def['pos']
        slot_ref = slot_def['ref']
        
        # 2. Generate Master Pool [V0, V1... Vn-1]
        master_pool = HapticTransformer.transform_vectors(n_atoms, slot_z, slot_ref, cone_spread)
        
        results = []
        
        # 3. Determine Stride/Mapping
        # Star Index is the atom index in the ligand fragment that matches V0 (North Star)
        star_idx = heading_atom_in_frag_idx
        
        for i in range(n_atoms):
            # i is atom index in fragment (0..N-1)
            
            # Logic: Delta from Star
            delta = (i - star_idx) * direction
            
            # Determine Vector Index k
            if n_atoms >= 5:
                # Rings wrap
                k = delta % n_atoms
            else:
                # Chains/Subsets
                # Logic from PRD: "Assuming subsets effectively behave like cyclic definitions for assignment"
                # "k = delta % n"
                k = delta % n_atoms
            
            assigned_vec = master_pool[k]
            results.append({
                'atom_idx_in_frag': i,
                'vector': assigned_vec
            })
            
        return results
