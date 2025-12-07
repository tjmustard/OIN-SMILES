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
        # ... (rest of method unchanged until loop)
        # We need to preserve the surrounding code, so I will provide the chunks.
        pass # Placeholder for replace logic

    # ... (skipping to parse_inline_string modification)
  
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
            slot = int(slot_str)  
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
        # Note: The 'w' tag ranks correspond to Ligand Indices (0..N-1)  
        # BUT the 'w' tag usually ignores the metal.  
        # Let's assume w-tag rank 0 is the first LIGAND (fragment 1).  
        # Actually in V2.4, OIN aligner (oin_aligner.py line 720+) usually uses rank based on sorted fragments.
        # But fragments in SMILES are also sorted.
        # So fragments[1] corresponds to rank 1 (if we count metal as 0)?
        # Let's check xyz2mol/oin_aligner.
        # In oin_aligner.py, fragments are processed.
        # The w tag is `rank.atom_idx:slot`.
        
        #Wait, w-tag format in V2.4 (OINDiscreteAligner):
        # `{rank_in_smiles}.{atom_idx}: {slot_idx}`
        # rank_in_smiles is the index of the fragment in the final sorted SMILES string.
        # So fragment 0 is metal. fragments[1] is rank 1.
        
        for i in range(1, len(fragments)):  
            lig_rank = i # 0-based index in SMILES string
            frag_smiles = fragments[i]  
              
            if lig_rank in slot_map:  
                slot = slot_map[lig_rank]  
                  
                # INJECTION LOGIC:  
                # We need to append [slot] to the Binding Atom(s).  
                # OIN-Inline syntax: Atom[Slot].  
                # [cH] -> [cH][0]  
                # [NH3] -> [NH3][2]  
                
                # We need to append [Slot] to the atoms that are actuall binding.
                # But the w-tag from aligner doesn't tell us WHICH atom is binding in this simple map logic, 
                # wait, w-tag DOES tell us: `rank.atom_idx:slot` 
                # But here we simplified `slot_map[rank] = slot` in the code above.
                # This assumes monodentate per fragment? Or 1 slot per fragment?
                # The PRD code example assumed: `w:1.0:0;2.0:1` -> Rank 1 (atom 0) binds to slot 0.
                
                # If we have `w:1.0:0`, it means Fragment 1, Atom 0 binds to Slot 0.
                # If we have chelators, we might have `w:1.0:0;1.5:1`.
                
                # So we should parse w-tag fully.
                pass
        
        # Re-parse w-tag properly to handle atom indices
        # Map: Rank -> List of (AtomIdx, Slot)
        detailed_map = {}
        
        for item in w_tag.split(";"):
            if not item: continue
            # Format: Rank.AtomIdx:Slot
            if ":" not in item: continue
            left, slot_str = item.split(":")
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
            detailed_map[rank].append((atom_idx, slot))
            
        # Re-build fragments
        # Skip metal 0
        
        
        # Mapping helpers
        import rdkit.Chem as Chem

        for i in range(1, len(fragments)):
            frag_smiles = fragments[i]
            lig_rank = i
            
            if lig_rank in detailed_map:
                binders = detailed_map[lig_rank]
                # binders is list of (atom_idx, slot)
                
                # Check consistency: if all atoms map to same slot, use simple replacement if brackets exist?
                # No, mixed strategy is bad. Let's use RDKit mapping for robustness if possible.
                
                try:
                    mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
                    if not mol:
                        # Fallback to simple replace
                        raise ValueError("Invalid SMILES")
                        
                    # Apply Map Numbers = Slot + 1 (to avoid 0)
                    for atom_idx, slot in binders:
                        if atom_idx < mol.GetNumAtoms():
                            # We use offset 1000 to avoid conflict with real maps (unlikely in OIN)
                            mol.GetAtomWithIdx(atom_idx).SetAtomMapNum(slot + 1000)
                            
                    # Generate SMILES with maps
                    mapped_smiles = Chem.MolToSmiles(mol, canonical=False) 
                    # canonical=False to hopefully preserve atom order if input was canonical?
                    # Ideally we want output to match input structure but with tags.
                    # RDKit might reorder.
                    # BUT xyz2mol generated the input Frag SMILES using RDKit canonical.
                    # So Chem.MolToSmiles(mol) should likely produce the same string + maps.
                    # Let's trust RDKit.
                    
                    # Regex replacement: [Symbol:Map] -> Symbol[Slot]
                    # Pattern: \[([^:\]]+):(\d+)\]
                    def replace_map(match):
                        content = match.group(1)
                        map_num = int(match.group(2))
                        slot = map_num - 1000
                        
                        # Heuristic to decide on brackets:
                        # If content is simple organic subset (c, n, C, N, O, Cl, F, Br, I, etc)
                        # AND contains no other characters (like H, +, -)... 
                        # We might strip brackets.
                        # BUT we need to be careful about implicit Hs.
                        # [c] (Zone A) -> [c][0] (Keep brackets to imply no implicit H if that was intention?)
                        # But [c:1000] -> [c][0]. Content is 'c'.
                        # If we output c[0], we lose the "bracket-ness".
                        # However, for [n:1000] -> content 'n'.
                        # User wants n[0] for Pyridine N.
                        # Pyridine N is 'n'. [n] is also valid but implies 0 H.
                        # 
                        # Safe Strategy: ALWAYS output [Content][Slot]?
                        # User output shows: N[3]. (Unbracketed N).
                        # So we prefer unbracketed if possible.
                        #
                        # Check if 'content' is valid unbracketed SMILES atom?
                        # B, C, N, O, P, S, F, Cl, Br, I + aromatic b, c, n, o, p, s.
                        

                        if content == 'NH3':
                            return f"N{{{slot}}}"

                        # Use explicit list for pure organic atoms that can be unbracketed
                        # B, C, N, O, P, S and aromatic versions.
                        # Exclude Halogens (Cl, Br, I, F) so they remain bracketed [Cl].
                        is_pure_organic = re.fullmatch(r"^(B|C|N|O|P|S|c|n|o|p|s)$", content)
                        
                        if is_pure_organic:
                            return f"{content}{{{slot}}}"
                        else:
                            return f"[{content}]{{{slot}}}"

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
                  
        return ".".join(new_fragments)

    @staticmethod  
    def parse_inline_string(inline_string: str) -> Tuple[str, str, List[Tuple[int, int, int]]]:  
        """  
        Parses OIN-Inline back to (SMILES, Geometry, VectorData).  
        Returns:  
            clean_smiles: Standard SMILES  
            geometry: String code  
            vector_data: List of (LigandRank, AtomInFragIdx, SlotIndex)
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
        vector_data = [] # List of (LigandRank, AtomIdx, SlotIndex)
          
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
            
            processed_frag = raw_frag
            
            # Helper to replace {Slot} with appropriate RDKit map syntax
            def inject_map(match):
                # match.g(0) is the whole thing e.g. "N{0}" or "[NH]{0}"
                pass
            
            # Regex for Bracketed Atom + Tag: `(\[[^\]]+\])\{(\d+)\}` -> `[NH]{0}`
            # We convert to `[NH:1000]`. Note: `[NH]` becomes `[NH:1000]`.
            
            def sub_bracketed(m):
                # m.group(1) is `[NH]` or `[c]`
                # m.group(2) is `0`
                atom_block = m.group(1)
                slot = int(m.group(2))
                map_num = 1000 + slot
                # Insert :map_num before last ]
                if ":" in atom_block:
                    # Already mapped? Replace it?
                    # `[NH:1]` -> `[NH:1000]`
                    return re.sub(r":\d+\]", f":{map_num}]", atom_block)
                else:
                    return atom_block[:-1] + f":{map_num}]"

            processed_frag = re.sub(r"(\[[^\]]+\])\{(\d+)\}", sub_bracketed, processed_frag)
            
            # Regex for Unbracketed Atom + Tag: `(?<!\])([A-Za-z][a-z]?)\{(\d+)\}` -> `N{0}`, `c{0}`
            # (?<!\]) ensures we don't match the closing bracket of previous match (redundant if ordered?)
            # Match `N{0}`, `Cl{1}`, `c{0}`
            
            def sub_unbracketed(match): 
                sym = match.group(1)
                slot = int(match.group(2))
                map_num = 1000 + slot
                # print(f"DEBUG: Match {sym}{{{slot}}} -> [{sym}:{map_num}]")
                return f"[{sym}:{map_num}]"
                
            processed_frag = re.sub(r"(?<!\])(se|as|[A-Z][a-z]?|[bcnops])\{(\d+)\}", sub_unbracketed, processed_frag)
            # print(f"DEBUG: Processed frag: {processed_frag}")
            
            # Now parse with RDKit to get canonical index
            mol = Chem.MolFromSmiles(processed_frag)
            if mol:
                # Iterate atoms to find maps
                final_atoms_found = []
                for atom in mol.GetAtoms():
                    map_num = atom.GetAtomMapNum()
                    if map_num >= 1000:
                        slot = map_num - 1000
                        atom_idx = atom.GetIdx()
                        vector_data.append( (lig_rank, atom_idx, slot) )
                        # Clear map for clean SMILES
                        atom.SetAtomMapNum(0)
                        
                        # Fix for [N] vs N (and similar):
                        # If the atom was N, and we forced [N:1000], RDKit MolToSmiles might output [N].
                        # But standard N (NH3) is 'N'. [N] is distinct.
                        # If we have [N] with 0 explicit Hs, it's a radical or N-? 
                        # Wait, [N] is nitrogen with no Hydrogens? 
                        # RDKit interpretation of [N] vs N:
                        # [N] has H count 0.
                        # N has valence derived H count (NH3).
                        # If we force [N:1000], we force current H count.
                        # If mol was built from 'N', it has implicit Hs. 
                        # When we build from '[N:1000]', RDKit assumes H count is 0?
                        # No, MolFromSmiles('[N:1000]') -> Atom N, Map 1000. H count?
                        # If we assume input OIN used 'N[0]', we want 'N'.
                        # If input OIN used '[NH2][0]', we want '[NH2]'.
                        
                        # So, simply clearing the map and regenerating SMILES should work IF RDKit respects the valence.
                        # However, MolFromSmiles('[N:1000]') creates an N with explicit degree 0?
                        # Let's check if we can remove explicit H setting or check radical.
                        pass
                
                # Generate clean SMILES from Mol to ensure canonical consistency
                # isomericSmiles=True, canonical=True
                # Ideally this matches output SMILES structure.
                clean_frag = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                
                # Manual fix for [N] -> N if it came from N[0]
                # If clean_frag is '[N]' but we expect 'N'?
                # The issue is RDKit might default to brackets if it feels like it.
                # Specifically for Ammonia (N) in OIN, we want 'N'. 
                if clean_frag == '[N]':
                     clean_frag = 'N'
                elif clean_frag == '[n]': # Pyridine N?
                     clean_frag = 'n'
                
                clean_fragments.append(clean_frag)
            else:
                # Fallback: Tag extraction via regex failed or RDKit didn't find maps.
                # This happens for strings like `c1cccc1[0]` where tag is appended to ring number.
                # In this case, we assume the tag applies to the Representative Atom (Idx 0).
                
                # Check if we have slots we found earlier
                slots_found = OINInlineHandler.SLOT_REGEX.findall(raw_frag)
                if slots_found:
                    # We need to clean the string manually
                    clean_frag = re.sub(r"\{\d+\}", "", raw_frag)
                    clean_fragments.append(clean_frag)
                    
                    # Extract unique slots
                    unique_slots = sorted(list(set([int(s) for s in slots_found])))
                    for slot in unique_slots:
                        # Assume Atom 0 is the anchor/representative
                        vector_data.append( (lig_rank, 0, slot) )
                else:
                    clean_frag = re.sub(r"\{\d+\}", "", raw_frag)
                    clean_fragments.append(clean_frag)

        final_smiles = ".".join(clean_fragments)  
        return final_smiles, geometry, vector_data
