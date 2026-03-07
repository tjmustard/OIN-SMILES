import dataclasses
from dataclasses import dataclass
from typing import List, Tuple
import numpy as np

from ..core.geometry_templates import TEMPLATES

@dataclass
class OINVector:
    atom_idx: int
    vector: Tuple[float, float, float]
    fragment_idx: int
    atom_in_fragment_idx: int
    haptic_heading: bool = False
    haptic_direction: int = 1
    slot_idx: int = -1

@dataclass
class ParsedOIN:
    smiles: str
    fragments: List[str]
    metal_fragment_idx: int
    vectors: List[OINVector]
    original_oin: str
    geometry: str = "UNK" # NEW: Architecture Template Code

class OINParser:
    def parse(self, oin_string: str) -> ParsedOIN:
        # Check for V3.0 Inline Topology
        # Heuristic: No "|" separator AND contains Metal tag like [Pt_SPL]
        from ..oin.inline import OINInlineHandler
        
        is_inline = False
        parts = oin_string.split("|")
        smiles = parts[0].strip()
        metadata = parts[1:] if len(parts) > 1 else []
        
        if len(parts) == 1 and OINInlineHandler.METAL_REGEX.search(oin_string):
             is_inline = True
             
        if is_inline:
             # Convert Inline -> Standard (Sidecar) components
             # Note: OINInlineHandler.parse_inline_string returns (smiles, geo, vector_list)
             # vector_list struct: (Rank, AtomIdx, Slot, IsHeading, Direction)
             
             smiles, geo_code, vector_data = OINInlineHandler.parse_inline_string(oin_string)
             fragments = smiles.split(".")
             metal_fragment_idx = 0 # Assumption
             
             tmpl_vectors = TEMPLATES.get(geo_code)
             vectors = []
             
             if tmpl_vectors is not None:
                for item in vector_data:
                    # Unpack
                    if len(item) == 5:
                        lig_rank, atom_in_fragment_idx, slot_idx, is_heading, direction = item
                    else:
                        # Fallback for old signature if needed (Rank, Atom, Slot)
                        lig_rank, atom_in_fragment_idx, slot_idx = item
                        is_heading = False
                        direction = 1

                    if slot_idx in tmpl_vectors:
                        resolved_vec = tmpl_vectors[slot_idx]['pos']
                        
                        vectors.append(OINVector(
                            atom_idx=-1,
                            vector=tuple(resolved_vec.tolist()),
                            fragment_idx=lig_rank,
                            atom_in_fragment_idx=atom_in_fragment_idx,
                            haptic_heading=is_heading,
                            haptic_direction=direction,
                            slot_idx=slot_idx
                        ))
             
             return ParsedOIN(
                smiles=smiles,
                fragments=fragments,
                metal_fragment_idx=metal_fragment_idx,
                vectors=vectors,
                original_oin=oin_string,
                geometry=geo_code
             )

        # Standard / Legacy Parsing
        
        fragments = smiles.split(".")
        
        # Identify metal fragment (usually 0, but could check symbol if needed)
        metal_fragment_idx = 0
        
        # 1. Identify Geometry Template First
        tmpl_vectors = None
        geo_code = "UNK"
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
                        # TODO: V2.4 Legacy does not support advanced Heading/Direction in w-tag explicitly??
                        # The regex in inline.py suggested it might be there.
                        # For now, assume default.
                        
                        if ":" not in item: continue
                        
                        indices_str, slot_str = item.split(":", 1)
                        # Check for < > ^ in slot_str
                        
                        is_heading = False
                        direction = 1
                        
                        if '<' in slot_str:
                            is_heading = True
                            direction = -1
                            slot_str = slot_str.replace('<', '')
                        elif '>' in slot_str:
                            is_heading = True
                            direction = 1
                            slot_str = slot_str.replace('>', '')
                        elif '^' in slot_str:
                            is_heading = True
                            direction = 1
                            slot_str = slot_str.replace('^', '')
                            
                        slot_idx = int(slot_str)
                        
                        if "." not in indices_str: continue
                             
                        frag_idx_str, atom_idx_str = indices_str.split(".")
                        frag_idx = int(frag_idx_str)
                        atom_in_frag_idx = int(atom_idx_str)
                        
                        # Resolve Vector
                        if slot_idx not in tmpl_vectors: continue # Safety
                        resolved_vec = tmpl_vectors[slot_idx]['pos']
                        
                        vectors.append(OINVector(
                            atom_idx= -1, 
                            vector=tuple(resolved_vec.tolist()),
                            fragment_idx=frag_idx,
                            atom_in_fragment_idx=atom_in_frag_idx,
                            haptic_heading=is_heading,
                            haptic_direction=direction,
                            slot_idx=slot_idx
                        ))
                    except ValueError:
                        continue 
                        
        return ParsedOIN(
            smiles=smiles,
            fragments=fragments,
            metal_fragment_idx=metal_fragment_idx,
            vectors=vectors,
            original_oin=oin_string,
            geometry=geo_code
        )
