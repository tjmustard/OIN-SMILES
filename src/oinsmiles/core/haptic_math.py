"""
Core mathematical logic for haptic ligand expansion and North Star calculations.
"""
import numpy as np
from typing import List, Tuple, Optional

def normalize(v: np.ndarray) -> np.ndarray:
    """Returns the normalized vector v."""
    norm = np.linalg.norm(v)
    return v / norm if norm > 1e-9 else v

def get_step_angle(n: int) -> float:
    """Returns the angular step in degrees for an N-membered ring/chain."""
    if n <= 1: return 0.0
    if n == 2: return 180.0
    if n == 3: return 72.0 # Subset of Pentagon (common for Allyl/Cp?) Or 120? 
                           # PRD doesn't specify rigid angles for 3, but Alpha used 72.0.
                           # Let's stick to regular polygon logic unless specified.
                           # Regular triangle = 120. But for Cp subset (Allyl), it is 72.
                           # The existing code used 72.0 for n=3. I will keep it for now.
    if n == 4: return 60.0 # Subset of Hexagon? Or 90?
    if n == 5: return 72.0
    if n == 6: return 60.0
    if n == 7: return 360.0 / 7.0
    if n >= 8: return 360.0 / float(n)
    return 360.0 / float(n)

def expand_slot(n_points: int, slot_z: np.ndarray, slot_x_ref: np.ndarray, cone_spread: float = 1.0) -> List[np.ndarray]:
    """
    Expands a single haptic slot into N vectors arranged in a polygon/cone.
    
    Args:
        n_points: Number of binding atoms (vertices of the polygon).
        slot_z: The primary slot vector (Centroid location).
        slot_x_ref: The reference vector defining the 'North Star' (0 deg) direction.
        cone_spread: The spread of the cone (radial distance from centroid vector).

    Returns:
        List of N normalized vectors. Vector 0 is aligned with slot_x_ref.
    """
    if n_points <= 1:
        return [normalize(slot_z)]

    # 1. Establish Basis Frame at the Slot
    Z = normalize(slot_z)
    
    # Project X_ref onto the plane perpendicular to Z
    X_ref = normalize(slot_x_ref)
    proj_X = X_ref - (np.dot(X_ref, Z) * Z)
    
    norm_x = np.linalg.norm(proj_X)
    if norm_x < 1e-6:
        # X_ref is parallel to Z (Degenerate). Pick arbitrary X.
        # Try Global X
        proj_X = np.array([1.0, 0.0, 0.0]) - (np.dot(np.array([1.0, 0.0, 0.0]), Z) * Z)
        if np.linalg.norm(proj_X) < 1e-6:
             # Try Global Y
            proj_X = np.array([0.0, 1.0, 0.0]) - (np.dot(np.array([0.0, 1.0, 0.0]), Z) * Z)
            
    X = normalize(proj_X)
    Y = np.cross(Z, X) # Right-handed system
    
    vectors = []
    
    # Logic: For n=3 (Allyl) and n=4, we might want specific subsets of 5/6 rings.
    # But for general "Regular Polygon" expansion:
    step_deg = get_step_angle(n_points)
    
    # If n=3 is allyl (subset of 5), the arc is 144 degrees total?
    # Or is it a regular triangle?
    # PRD says "Regular polygon/cone". 
    # But existing code had specific overrides.
    # "Union of OIN v3.6 Logic": Haptic expansion usually assumes a regular polygon for full rings.
    # For open chains, it might differ. However, to ensure invertibility, we must be consistent.
    # If we generate a specific set of vectors, the canonicalizer will snap to them.
    # The crucial part is that OIN string {0} maps to Vector 0 (Aligned with X).
    
    for k in range(n_points):
        # We distribute them symmetrically?
        # Standard: Start at 0, step by step_deg.
        # For full rings (n=5, 6), sum(step) = 360.
        # For chains (n=3), if step=72, total arc = 144. It's an open arc.
        # This matches the "subset of ring" logic for Allyl.
        
        theta_rad = np.radians(k * step_deg)
        
        # Planar components on the cone base
        vx = np.cos(theta_rad)
        vy = np.sin(theta_rad)
        
        # Construct Vector: Z + spread * (vx*X + vy*Y)
        # Note: This creates a cone.
        vec = Z + cone_spread * (vx * X + vy * Y)
        vectors.append(normalize(vec))
        
    return vectors

def identify_north_star_index(ligand_local_vectors: List[np.ndarray], ref_vector: np.ndarray) -> int:
    """
    Identifies the index of the 'North Star' atom in the ligand.
    The North Star is the atom with the Maximum Dot Product against X_ref.
    
    Args:
        ligand_local_vectors: List of vectors for each atom relative to the centroid (or origin).
        ref_vector: The reference vector (X_ref) to align against.
        
    Returns:
        Index of the atom in the list that is the North Star.
    """
    best_idx = -1
    max_dot = -float('inf')
    
    ref_norm = normalize(ref_vector)
    
    for i, vec in enumerate(ligand_local_vectors):
        dot = np.dot(normalize(vec), ref_norm)
        if dot > max_dot:
            max_dot = dot
            best_idx = i
            
    return best_idx
