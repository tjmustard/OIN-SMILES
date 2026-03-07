import re
from typing import Tuple, List, Optional

class OINInlineHandler:  
    """  
    Experimental Handler for V3.0 OIN-Inline strings.  
    """  
      
    # Regex to find [Element_Geo] e.g. [Pt_SPL]  
    METAL_REGEX = re.compile(r"\[([A-Z][a-z]?)\_([A-Z]{3})\]")  
      
    # Regex to find [SlotIndex] e.g. {0}, {12}  
    # We look for braces containing ONLY digits  
    SLOT_REGEX = re.compile(r"\{(\d+)\}")

    @staticmethod  
    def generate_inline_string(oin_v2_string: str) -> str:
        """  
        Converts a V2.4 string (|w:..|) to V3.0 Inline format.  
        Input: "[Pt].[Cl]... |g:SPL|w:1.0:0;2.0:1..."  
        """  
        if "|" not in oin_v2_string:  
            return oin_v2_string # Not an OIN string  
              
        # 1. Parse V2.4 components  
        smiles_part, metadata = oin_v2_string.split(" |", 1)  
          
        # Parse Geometry
        geo_match = re.search(r"g:([A-Z]{3})", metadata)
        
        # Parse Bond Stereo
        b_match = re.search(r"b:([^|]+)", metadata)
        b_tag = f"|b:{b_match.group(1)}" if b_match else ""

        if not geo_match:
            # It might be a partial OIN string or missing geometry, just return as is
            return oin_v2_string
            
        geometry = geo_match.group(1)  
          
        # Parse Vector Map: {LigandRank: SlotIndex}  
        w_match = re.search(r"w:([^|]+)", metadata)
        if not w_match:
             return oin_v2_string

        w_tag = w_match.group(1)  
        # w:1.0:0;2.0:1 -> Map {1: 0, 2: 1}  
        # Note: Ligand Rank is 0-indexed based on the order in the SMILES string.  
        # But V2.4 OIN SMILES is dot-separated.  
          
        slot_map = {} # Key: Ligand Index (0-based), Value: Slot Index  
        for item in w_tag.split(";"):  
            if not item: continue
            if ":" not in item: continue
            rank_str, slot_str = item.split(":")  
            # Rank might be "1.0" or "1". Convert to int.  
            rank = int(float(rank_str))  
            slot = int(slot_str.replace('^', '').replace('>', '').replace('<', '')) # Sanitize heading marker
            slot_map[rank] = slot  
              
        # 2. Tokenize SMILES to inject tags  
        fragments = smiles_part.split(".")  
          
        # Metal is usually first (Rank 0? No, Metal is separate).  
        # We assume standard OIN order: [Metal].L1.L2...  
          
        new_fragments = []  
          
        # Handle Metal (Index 0)  
        metal_frag = fragments[0]  
        # Inject Geometry: [Pt] -> [Pt_SPL]  
        # Remove trailing ']' and append '_GEO]'  
        if metal_frag.endswith("]"):  
            new_metal = metal_frag[:-1] + f"_{geometry}]"  
        else:  
            # Metal might be unbracketed? Unlikely in OIN.  
            new_metal = f"[{metal_frag}_{geometry}]"  
        new_fragments.append(new_metal)  
          
        # Handle Ligands (Indices 1..N)  
        # Re-parse w-tag properly to handle atom indices
        # Map: Rank -> List of (AtomIdx, Slot)
        detailed_map = {}
        
        for item in w_tag.split(";"):
            if not item: continue
            # Format: Rank.AtomIdx:Slot
            if ":" not in item: continue
            left, slot_str = item.split(":")
            
            # Need to strip heading chars ^, >, <
            heading_char = ""
            for c in ['^', '>', '<']:
                if c in slot_str:
                    heading_char = c
                    break
            
            slot_str = slot_str.replace('^', '').replace('>', '').replace('<', '')
            slot = int(slot_str)
            
            if "." in left:
                rank_str, atom_idx_str = left.split(".")
                rank = int(rank_str)
                atom_idx = int(atom_idx_str)
            else:
                 # Fallback just rank
                 rank = int(float(left))
                 atom_idx = 0 # Assume 0
            
            if rank not in detailed_map:
                detailed_map[rank] = []
            detailed_map[rank].append((atom_idx, slot, heading_char))
            
        # Re-build fragments
        # Skip metal 0
        import rdkit.Chem as Chem

        for i in range(1, len(fragments)):
            frag_smiles = fragments[i]
            lig_rank = i
            
            if lig_rank in detailed_map:
                binders = detailed_map[lig_rank]
                # binders is list of (atom_idx, slot)
                
                try:
                    mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
                    if not mol:
                        # Fallback to simple replace
                        raise ValueError("Invalid SMILES")
                        
                    # Apply Map Numbers = Slot + 1000 (normal) or + 2000 (heading >) or + 3000 (heading <)
                    for atom_idx, slot, heading_char in binders:
                        if atom_idx < mol.GetNumAtoms():
                            # Validate Slot
                            offset = 1000
                            if heading_char == '<': offset = 3000
                            elif heading_char in ['>', '^']: offset = 2000
                            
                            mol.GetAtomWithIdx(atom_idx).SetAtomMapNum(slot + offset)
                            
                    # Restore ring info for aromaticity in SMILES
                    Chem.FastFindRings(mol)
                    # Generate SMILES with maps
                    mapped_smiles = Chem.MolToSmiles(mol, canonical=False) 
                    
                    # Regex replacement: [Symbol:Map] -> Symbol[Slot]
                    def replace_map(match):
                        content = match.group(1)
                        map_num = int(match.group(2))
                        
                        is_heading = False
                        heading_char = ""
                        
                        if map_num >= 3000:
                            slot = map_num - 3000
                            is_heading = True
                            heading_char = "<"
                        elif map_num >= 2000:
                            slot = map_num - 2000
                            is_heading = True
                            heading_char = ">" # Normalize ^ to >
                        else:
                            slot = map_num - 1000
                        
                        suffix = heading_char if is_heading else ""
                        
                        if content == 'NH3':
                            return f"N{{{slot}{suffix}}}"

                        is_pure_organic = re.fullmatch(r"^(N|P|S|c|n|o|p|s)$", content)
                        
                        if is_pure_organic:
                            return f"{content}{{{slot}{suffix}}}"
                        else:
                            return f"[{content}]{{{slot}{suffix}}}"

                    tagged_frag = re.sub(r"\[([^:\]]+):(\d+)\]", replace_map, mapped_smiles)
                    new_fragments.append(tagged_frag)
                    
                except Exception:
                     # Fallback to old REPLACE logic if RDKit fails or something
                    slots = [s for _, s in binders]
                    unique_slots = list(set(slots))
                    if len(unique_slots) == 1:
                        slot = unique_slots[0]
                        if "]" in frag_smiles:
                            tagged_frag = frag_smiles.replace("]", f"]{{{slot}}}")
                        else:
                            tagged_frag = f"{frag_smiles}{{{slot}}}"
                        new_fragments.append(tagged_frag)
                    else:
                        # Fallback for chelator: tag first slot only (Incorrect but safe-ish)
                        slot = unique_slots[0]
                        if "]" in frag_smiles:
                            tagged_frag = frag_smiles.replace("]", f"]{{{slot}}}")
                        else:
                            tagged_frag = f"{frag_smiles}{{{slot}}}"
                        new_fragments.append(tagged_frag)
            else:
                new_fragments.append(frag_smiles)
                  
        return ".".join(new_fragments) + b_tag

    @staticmethod  
    def parse_inline_string(inline_string: str) -> Tuple[str, str, List[Tuple[int, int, int, bool, int]]]:  
        """  
        Parses OIN-Inline back to (SMILES, Geometry, VectorData).  
        Returns:  
            clean_smiles: Standard SMILES  
            geometry: String code  
            vector_data: List of (LigandRank, AtomInFragIdx, SlotIndex, IsHeading, Direction)
        """  
        import rdkit.Chem as Chem

        # 1. Extract Geometry from Metal  
        # Look for [Pt_SPL]  
        metal_match = OINInlineHandler.METAL_REGEX.search(inline_string)  
        if not metal_match:  
            return inline_string, "", []
              
        element_sym = metal_match.group(1)  
        geometry = metal_match.group(2)  
          
        # Revert metal to standard [Pt]  
        clean_string = OINInlineHandler.METAL_REGEX.sub(f"[{element_sym}]", inline_string)  
        
        # 2. Extract Slots via RDKit Map Numbers
        fragments = clean_string.split(".")  
        clean_fragments = []  
        vector_data = [] # List of (LigandRank, AtomIdx, SlotIndex, IsHeading, Direction)
          
        # Metal is frag 0  
        clean_fragments.append(fragments[0])  
          
        for i in range(1, len(fragments)):  
            raw_frag = fragments[i]  
            lig_rank = i 
            
            # Identify slots and inject map numbers for RDKit parsing
            # Pattern 1: [Element][Slot] -> [Element:1xxx]
            # Pattern 2: Element[Slot] -> [Element:1xxx] (Unbracketed)
            
            # We need to act on the string.
            # Convert [0] -> :1000] inside brackets?
            # Or append map number?
            
            # Strategy: Replace `{Slot}` with `:<1000+slot>`.
            # But where?
            # If `[NH]{0}`, we want `[NH:1000]`.
            # If `N{0}`, we want `[N:1000]`.
            
            # NOTE: We need to handle < and > and ^ inside keys too if they exist in the input string?
            # The regex for keys was `\{(\d+)[\>\<]?\}`?
            # Current SLOT_REGEX is `\{(\d+)\}` which misses `<`.
            # We need to update regex handling here locally since we are parsing raw string again.
            
            processed_frag = raw_frag
            
            # 2a. Handle Directionality in Tag
            # Format: {0}, {0<}, {0>}, {0^}
            
            def get_map_num(slot_str, direction_char):
                slot = int(slot_str)
                if direction_char == '<':
                    return 3000 + slot
                elif direction_char in ['>', '^']:
                    return 2000 + slot
                else:
                    return 1000 + slot

            # Regex for Bracketed Atom + Tag: `(\[[^\]]+\])\{(\d+)([\>\<]|\^)?\}`
            
            def sub_bracketed(m):
                # m.group(1) is `[NH]` or `[c]`
                # m.group(2) is `0`
                # m.group(3) is `<` or None
                atom_block = m.group(1)
                slot_code = m.group(2)
                dir_char = m.group(3) if m.group(3) else ""
                
                map_num = get_map_num(slot_code, dir_char)
                
                # Insert :map_num before last ]
                if ":" in atom_block:
                    # Already mapped? Replace it?
                    # `[NH:1]` -> `[NH:1000]`
                    return re.sub(r":\d+\]", f":{map_num}]", atom_block)
                else:
                    return atom_block[:-1] + f":{map_num}]"

            processed_frag = re.sub(r"(\[[^\]]+\])\{(\d+)([\>\<]|\^)?\}", sub_bracketed, processed_frag)
            
            # Regex for Unbracketed Atom + Tag: `(?<!\])([A-Za-z][a-z]?)\{(\d+)([\>\<]|\^)?\}`
            
            def sub_unbracketed(match): 
                sym = match.group(1)
                slot_code = match.group(2)
                dir_char = match.group(3) if match.group(3) else ""
                
                map_num = get_map_num(slot_code, dir_char)
                
                # print(f"DEBUG: Match {sym}{{{slot}}} -> [{sym}:{map_num}]")
                return f"[{sym}:{map_num}]"
                
            # RATIONALE: We removed the negative lookbehind '(?<!\])' that was here previously.
            # This allows us to parse tags on atoms that immediately follow a bracketed group,
            # e.g. [H]N{2} or [CH3]C{0}.
            # Without this, such atoms were skipped and resulted in parse errors.
            processed_frag = re.sub(r"(se|as|[A-Z][a-z]?|[bcnops])\{(\d+)([\>\<]|\^)?\}", sub_unbracketed, processed_frag)
            # print(f"DEBUG: Processed frag: {processed_frag}")
            
            # Now parse with RDKit to get canonical index
            mol = Chem.MolFromSmiles(processed_frag)
            if not mol:
                mol = Chem.MolFromSmiles(processed_frag, sanitize=False)
                
            if mol:
                # Iterate atoms to find maps
                # final_atoms_found = []
                for atom in mol.GetAtoms():
                    map_num = atom.GetAtomMapNum()
                    
                    slot = -1
                    is_heading = False
                    direction = 1
                    
                    if map_num >= 3000:
                        slot = map_num - 3000
                        is_heading = True
                        direction = -1
                    elif map_num >= 2000:
                        slot = map_num - 2000
                        is_heading = True
                        direction = 1
                    elif map_num >= 1000:
                        slot = map_num - 1000
                        is_heading = False
                        direction = 1
                    
                    if slot >= 0:
                        atom_idx = atom.GetIdx()
                        vector_data.append( (lig_rank, atom_idx, slot, is_heading, direction) )
                        # Clear map for clean SMILES
                        atom.SetAtomMapNum(0)
                        
                        # Fix for [N] vs N (and similar):
                        pass
                
                # Generate clean SMILES
                # Use regex stripping on raw_frag to strictly preserve atom order/indices 
                # AND chiral markers (@/@@) which RDKit might otherwise normalize away 
                # during round-trips if not careful.
                clean_frag = re.sub(r"\{\d+([\>\<]|\^)?\}", "", raw_frag)
                
                # Manual fix for [N] -> N if it came from N{0} and looks bracketed?
                # If input was N{0}, raw_frag is N{0}, clean is N.
                # If input was [N]{0}, raw_frag is [N]{0}, clean is [N].
                clean_fragments.append(clean_frag)
            else:
                # Fallback: Tag extraction via regex failed or RDKit didn't find maps.
                # This happens for strings like `c1cccc1[0]` where tag is appended to ring number.
                # In this case, we assume the tag applies to the Representative Atom (Idx 0).
                
                # Check if we have slots we found earlier
                # Update Regex to support direction
                RAW_SLOT_REGEX = re.compile(r"\{(\d+)([\>\<]|\^)?\}")
                slots_found = RAW_SLOT_REGEX.findall(raw_frag)
                if slots_found:
                    # We need to clean the string manually
                    clean_frag = re.sub(r"\{\d+([\>\<]|\^)?\}", "", raw_frag)
                    clean_fragments.append(clean_frag)
                    
                    # Extract unique slots
                    # slots_found is list of tuples (slot, dir)
                    # We assume 1 slot for the whole fragment if regex fallback is used?
                    # Or multiple?
                    
                    # Let's deduplicate by slot ID, but keep first direction found?
                    # Usually `c1ccccc1{0}` means whole thing is 0.
                    
                    seen_slots = set()
                    for slot_str, dir_char in slots_found:
                        slot = int(slot_str)
                        if slot in seen_slots: continue
                        seen_slots.add(slot)
                        
                        is_heading = False
                        direction = 1
                        if dir_char == '<':
                             is_heading = True
                             direction = -1
                        elif dir_char in ['>', '^']:
                             is_heading = True
                             direction = 1
                             
                        # Assume Atom 0 is the anchor/representative
                        vector_data.append( (lig_rank, 0, slot, is_heading, direction) )
                else:
                    clean_frag = re.sub(r"\{\d+([\>\<]|\^)?\}", "", raw_frag)
                    clean_fragments.append(clean_frag)

        final_smiles = ".".join(clean_fragments)  
        return final_smiles, geometry, vector_data
