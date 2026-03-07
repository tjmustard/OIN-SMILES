from typing import Dict, List, Tuple
import re

class OINParser:
    def __init__(self):
        pass

    def parse(self, oin_string: str) -> Tuple[str, Dict[str, str]]:
        """
        Parses an OIN string into the canonical SMILES and a dictionary of tags.
        Format: [Canonical SMILES] | [OIN Block]
        Example: [Pt]... |w:0:0,0,0|d:1.0|
        """
        # Check for V3.0 Inline Topology First
        # Heuristic: No "|" separator AND contains Metal tag like [Pt_SPL]
        from .inline import OINInlineHandler
        
        is_inline = False
        parts = oin_string.split("|")
        
        if len(parts) == 1 and OINInlineHandler.METAL_REGEX.search(oin_string):
             is_inline = True
             
        if is_inline:
             smiles, geo_str, vector_data = OINInlineHandler.parse_inline_string(oin_string)
             tags = {}
             if geo_str:
                 tags['g'] = geo_str
                 
             # Reconstruct 'w' tag from vector_data list of (rank, atom_idx, slot, is_heading, direction)
             w_entries = []
             for rank, a_idx, slot, is_h, direction in vector_data:
                 dir_char = ""
                 if is_h:
                     dir_char = ">" if direction > 0 else "<"
                 w_entries.append(f"{rank}.{a_idx}:{slot}{dir_char}")
                 
             if w_entries:
                 tags['w'] = ";".join(w_entries)
                 
             return smiles, tags

        if len(parts) < 2:
            # No OIN block found, return just the SMILES
            return parts[0].strip(), {}

        smiles = parts[0].strip()
        # The rest are tags, possibly with empty strings if there are trailing pipes
        raw_tags = [p.strip() for p in parts[1:] if p.strip()]
        
        tags = {}
        for tag_str in raw_tags:
            # Each tag_str should be like "w:..." or "d:..."
            if ":" in tag_str:
                key, value = tag_str.split(":", 1)
                tags[key] = value
            
        return smiles, tags

    def _parse_tags(self, oin_block: str) -> Dict[str, str]:
        """
        Deprecated internal method. Parsing is now done in parse().
        Kept for compatibility if needed, but logic moved to parse.
        """
        return {}

    def parse_coordinates(self, tag_content: str) -> List[Tuple[int, float, float, float]]:
        """
        Parses the 'v' tag into a list of (idx, x, y, z) tuples.
        Supports V1.4 (v:Metal.Ligand:x,y,z) and legacy V1.3 (w:Ligand:x,y,z).
        Returns list of (LigandIdx, x, y, z).
        """
        coords = []
        if not tag_content:
            return coords
            
        entries = tag_content.split(";")
        for entry in entries:
            try:
                if ":" not in entry: continue
                parts = entry.split(":")
                
                # Check for V1.4 format: Metal.Ligand
                if "." in parts[0]:
                    metal_str, ligand_str = parts[0].split(".")
                    idx = int(ligand_str)
                else:
                    # Legacy or simple format
                    idx = int(parts[0])
                    
                x, y, z = map(float, parts[1].split(","))
                coords.append((idx, x, y, z))
            except ValueError:
                continue 
        return coords

    def parse_connectivity(self, tag_content: str) -> List[Tuple[int, int]]:
        """
        Parses the 'v' tag to extract connectivity.
        Returns list of (LigandIdx, MetalIdx) tuples.
        """
        bonds = []
        if not tag_content:
            return bonds
            
        entries = tag_content.split(";")
        for entry in entries:
            try:
                if ":" not in entry: continue
                parts = entry.split(":")
                if "." in parts[0]:
                     metal_str, ligand_str = parts[0].split(".")
                     bonds.append((int(ligand_str), int(metal_str)))
            except ValueError:
                continue
        return bonds
    def parse_bond_stereo(self, b_tag: str) -> List[Tuple[int, int, str]]:
        """
        Parses OIN '@b' tag: 'u1-v1:S1;u2-v2:S2'
        Returns: List of (u, v, stereo_str)
        """
        results = []
        if not b_tag:
            return results
        for entry in b_tag.split(";"):
            if not entry: continue
            try:
                # Format: u-v:stereo
                indices_part, stereo = entry.split(":")
                u_str, v_str = indices_part.split("-")
                results.append((int(u_str), int(v_str), stereo))
            except ValueError:
                continue
        return results
