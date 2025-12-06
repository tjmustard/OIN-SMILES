import dataclasses
from dataclasses import dataclass
from typing import List, Tuple
import numpy as np

# --- V2.3 Templates for Parser Resolution ---
def normalize_template(arr):
    return arr / np.linalg.norm(arr, axis=1)[:, None]

TEMPLATES = {
    'LIN': np.array([[0,0,1], [0,0,-1]]),
    'TPL': np.array([[0,1,0], [0.8660254,-0.5,0], [-0.8660254,-0.5,0]]),
    'SPL': np.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0]]),
    'TET': np.array([
        [ 1,  1,  1], [ 1, -1, -1], [-1,  1, -1], [-1, -1,  1]
    ]),
    'TPY': np.array([
        [0,0,1], 
        [0,1,0], [0.8660254,-0.5,0], [-0.8660254,-0.5,0] 
    ]),
    'TBP': np.array([
        [0,0,1], [0,0,-1],
        [0,1,0], [0.8660254,-0.5,0], [-0.8660254,-0.5,0]
    ]),
    'OCT': np.array([
        [0,0,1], [0,0,-1],
        [1,0,0], [-1,0,0], [0,1,0], [0,-1,0]
    ])
}

@dataclass
class OINVector:
    atom_idx: int
    vector: Tuple[float, float, float]
    fragment_idx: int
    atom_in_fragment_idx: int

@dataclass
class ParsedOIN:
    smiles: str
    fragments: List[str]
    metal_fragment_idx: int
    vectors: List[OINVector]
    original_oin: str

class OINParser:
    def parse(self, oin_string: str) -> ParsedOIN:
        parts = oin_string.split("|")
        smiles = parts[0].strip()
        metadata = parts[1:] if len(parts) > 1 else []
        
        fragments = smiles.split(".")
        
        # Identify metal fragment (usually 0, but could check symbol if needed)
        metal_fragment_idx = 0
        
        # 1. Identify Geometry Template First
        tmpl_vectors = None
        for meta in metadata:
            if meta.startswith("g:"):
                # g:SPL
                geo_code = meta[2:]
                tmpl_vectors = TEMPLATES.get(geo_code)
                # Don't break, continue processing
                
        vectors = []
        
        for meta in metadata:
            if meta.startswith("w:"):
                if tmpl_vectors is None:
                    # Cannot resolve vectors without geometry
                    continue
                    
                # Format: w:Rank.Idx:Slot;...
                content = meta[2:]
                items = content.split(";")
                for item in items:
                    if not item: continue
                    try:
                        # item format "Rank.Idx:Slot"
                        if ":" not in item: continue
                        
                        indices_str, slot_str = item.split(":", 1)
                        slot_idx = int(slot_str)
                        
                        if "." not in indices_str: continue
                             
                        frag_idx_str, atom_idx_str = indices_str.split(".")
                        frag_idx = int(frag_idx_str)
                        atom_in_frag_idx = int(atom_idx_str)
                        
                        # Resolve Vector
                        if slot_idx >= len(tmpl_vectors): continue # Safety
                        resolved_vec = tmpl_vectors[slot_idx]
                        
                        vectors.append(OINVector(
                            atom_idx= -1, 
                            vector=tuple(resolved_vec.tolist()),
                            fragment_idx=frag_idx,
                            atom_in_fragment_idx=atom_in_frag_idx
                        ))
                    except ValueError:
                        continue 
                        
        return ParsedOIN(
            smiles=smiles,
            fragments=fragments,
            metal_fragment_idx=metal_fragment_idx,
            vectors=vectors,
            original_oin=oin_string
        )
