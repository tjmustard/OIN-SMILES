import re
from typing import List, NamedTuple, Optional, Tuple


def _count_smiles_atoms_before(smiles: str, pos: int) -> int:
    """Return the 0-based index of the last SMILES atom seen before *pos*.

    Scans *smiles* character by character, counting atoms (bracket atoms and
    organic-subset atoms).  Stops just before *pos*, which is the start of a
    ``{slot}`` marker.  The returned value is the atom index of the atom
    immediately preceding that marker.

    Used by ``parse_inline_string`` to convert slot-marker positions into
    RDKit atom indices without a full SMILES round-trip.
    """
    atom_count = -1
    i = 0
    while i < pos:
        ch = smiles[i]
        if ch == "[":
            end = smiles.index("]", i)
            atom_count += 1
            i = end + 1
        elif smiles[i : i + 2] in ("Cl", "Br"):
            atom_count += 1
            i += 2
        elif ch in "BCNOPSFIcnops":
            atom_count += 1
            i += 1
        else:
            i += 1
    return max(atom_count, 0)


class SlotAssignment(NamedTuple):
    """A single ligand-atom-to-slot binding, with optional winding direction.

    Only the producer in ``parse_inline_string`` constructs this; a positional
    4th value present means winding was captured (`'>'`/`'<'`/`None`).
    """

    lig_rank: int
    atom_idx: int
    slot: int
    winding: Optional[str] = None


class OINInlineHandler:
    """Experimental handler for V3.0 OIN-Inline strings."""

    # Regex to find [Element_Geo] e.g. [Pt_SPL] or [Pt@SP1_SPL] / [Ir@OH10_OCT].
    # The optional (?:@[A-Z0-9]+)? group captures V3.4+ absolute-config markers.
    METAL_REGEX = re.compile(r"\[([A-Z][a-z]?)(?:@[A-Z0-9]+)?\_([A-Z]{3})\]")

    # Regex to find slot markers e.g. {0}, {12}, {0>}, {1<}, {0^} (winding/heading markers)
    SLOT_REGEX = re.compile(r"\{(\d+)([><^])?\}")

    @staticmethod
    def generate_inline_string(oin_v2_string: str) -> str:
        """Convert a V2.4 sidecar OIN string into the V3.0 inline format."""
        # ... (rest of method unchanged until loop)
        # We need to preserve the surrounding code, so I will provide the chunks.
        pass  # Placeholder for replace logic

        # ... (skipping to parse_inline_string modification)

        """
        Converts a V2.4 string (|w:..|) to V3.0 Inline format.
        Input: "[Pt].[Cl]... |g:SPL|w:1.0:0;2.0:1..."
        """
        if "|" not in oin_v2_string:
            return oin_v2_string  # Not an OIN string

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

        slot_map = {}  # Key: Ligand Index (0-based), Value: Slot Index
        for item in w_tag.split(";"):
            if not item:
                continue
            if ":" not in item:
                continue
            rank_str, slot_str = item.split(":")
            # Rank might be "1.0" or "1". Convert to int.
            rank = int(float(rank_str))
            slot = int(
                slot_str.replace("^", "").replace(">", "").replace("<", "")
            )  # Sanitize heading marker
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
        # Actually in V2.4, OIN aligner (oin_aligner.py line 720+) usually uses
        # rank based on sorted fragments.
        # But fragments in SMILES are also sorted.
        # So fragments[1] corresponds to rank 1 (if we count metal as 0)?
        # Let's check xyz2mol/oin_aligner.
        # In oin_aligner.py, fragments are processed.
        # The w tag is `rank.atom_idx:slot`.

        # Wait, w-tag format in V2.4 (OINDiscreteAligner):
        # `{rank_in_smiles}.{atom_idx}: {slot_idx}`
        # rank_in_smiles is the index of the fragment in the final sorted SMILES string.
        # So fragment 0 is metal. fragments[1] is rank 1.

        for i in range(1, len(fragments)):
            lig_rank = i  # 0-based index in SMILES string
            frag_smiles = fragments[i]

            if lig_rank in slot_map:
                slot = slot_map[lig_rank]

                # INJECTION LOGIC:
                # We need to append [slot] to the Binding Atom(s).
                # OIN-Inline syntax: Atom[Slot].
                # [cH] -> [cH][0]
                # [NH3] -> [NH3][2]

                # We need to append [Slot] to the atoms that are actuall binding.
                # But the w-tag from aligner doesn't tell us WHICH atom is
                # binding in this simple map logic,
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
        detailed_map: dict = {}

        for item in w_tag.split(";"):
            if not item:
                continue
            # Format: Rank.AtomIdx:Slot
            if ":" not in item:
                continue
            left, slot_str = item.split(":")

            # Need to strip heading chars ^, >, <
            heading_char = ""
            for c in ["^", ">", "<"]:
                if c in slot_str:
                    heading_char = c
                    break

            slot_str = slot_str.replace("^", "").replace(">", "").replace("<", "")
            slot = int(slot_str)

            if "." in left:
                rank_str, atom_idx_str = left.split(".")
                rank = int(rank_str)
                atom_idx = int(atom_idx_str)
            else:
                # Fallback just rank
                rank = int(float(left))
                atom_idx = 0  # Assume 0

            if rank not in detailed_map:
                detailed_map[rank] = []
            detailed_map[rank].append((atom_idx, slot, heading_char))

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

                # Check consistency: if all atoms map to same slot, use simple
                # replacement if brackets exist?
                # No, mixed strategy is bad. Let's use RDKit mapping for robustness if possible.

                try:
                    mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
                    if not mol:
                        # Fallback to simple replace
                        raise ValueError("Invalid SMILES")

                    # Apply Map Numbers = Slot + 1000 (normal) or + 2000
                    # (heading >) or + 3000 (heading <)
                    for atom_idx, slot, heading_char in binders:
                        if atom_idx < mol.GetNumAtoms():
                            # Validate Slot
                            offset = 1000
                            if heading_char == "<":
                                offset = 3000
                            elif heading_char in [">", "^"]:
                                offset = 2000

                            mol.GetAtomWithIdx(atom_idx).SetAtomMapNum(slot + offset)

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

                        is_heading = False
                        heading_char = ""

                        if map_num >= 3000:
                            slot = map_num - 3000
                            is_heading = True
                            heading_char = "<"
                        elif map_num >= 2000:
                            slot = map_num - 2000
                            is_heading = True
                            heading_char = ">"  # Normalize ^ to >
                        else:
                            slot = map_num - 1000

                        suffix = heading_char if is_heading else ""

                        # Heuristic to decide on brackets:
                        # If content is simple organic subset (c, n, C, N, O, Cl, F, Br, I, etc)
                        # AND contains no other characters (like H, +, -)...
                        # We might strip brackets.
                        # BUT we need to be careful about implicit Hs.
                        # [c] (Zone A) -> [c][0] (Keep brackets to imply no
                        # implicit H if that was intention?)
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

                        if content == "NH3":
                            return f"N{{{slot}{suffix}}}"

                        # Terminal (heavy==0) nitride 0-H marker. A bare `N{n}` donor
                        # with no heavy neighbour is undecidable between a nitride
                        # `[N]` (0 H) and an ammine `[NH3]` (3 H) -- both otherwise
                        # serialize identically. Ammine is emitted as bare `N{n}` by
                        # the `content == "NH3"` branch above; for the 0-H nitride we
                        # keep the bracket (`[N]{n}`) so the round trip stays 0-H: the
                        # generator's `NoImplicit` guard treats a bracket atom as
                        # authoritative and does not re-protonate it, and the
                        # comparator keys `[N]` distinct from ammine. A bare `N` with a
                        # heavy neighbour is 0-H amido/imido by S1's exact convention
                        # and MUST stay de-bracketed (falls through below), so this is
                        # gated strictly on GetDegree() == 0.
                        if content == "N":
                            n_atom = next(
                                (a for a in mol.GetAtoms() if a.GetAtomMapNum() == map_num),
                                None,
                            )
                            if n_atom is not None and n_atom.GetDegree() == 0:
                                return f"[N]{{{slot}{suffix}}}"

                        # Use explicit list for pure organic atoms that can be unbracketed
                        # B, C, N, O, P, S and aromatic versions.
                        # Exclude Halogens (Cl, Br, I, F) so they remain bracketed [Cl].
                        # Also handle NH3 special case
                        is_pure_organic = re.fullmatch(r"^(B|C|N|O|P|S|c|n|o|p|s)$", content)

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

        return ".".join(new_fragments)

    @staticmethod
    def parse_inline_string(inline_string: str) -> Tuple[str, str, List[SlotAssignment]]:
        """Parse OIN-Inline back to (SMILES, Geometry, VectorData).

        Returns non-canonical SMILES with @/@@ markers preserved.
        Uses regex-only slot stripping — no MolFromSmiles/MolToSmiles round-trip.

        Returns:
        -------
        clean_smiles : str
            SMILES with slot markers stripped.  @/@@ descriptors are preserved.
        geometry : str
            Three-letter geometry code (e.g. 'SPL', 'OC6').
        vector_data : list of SlotAssignment
            Slot assignments, including winding direction (`'>'`/`'<'`/`None`,
            with `'^'` normalized to `'>'` on capture).
        """
        # 1. Extract geometry from [Metal_GEO] token
        metal_match = OINInlineHandler.METAL_REGEX.search(inline_string)
        if not metal_match:
            return inline_string, "", []

        element_sym = metal_match.group(1)
        geometry = metal_match.group(2)

        # Revert metal token to plain [Metal]
        clean_string = OINInlineHandler.METAL_REGEX.sub(f"[{element_sym}]", inline_string)

        # 2. Split on fragment separator and process each ligand fragment
        fragments = clean_string.split(".")
        clean_fragments: List[str] = [fragments[0]]  # metal fragment unchanged
        vector_data: List[SlotAssignment] = []

        for lig_rank in range(1, len(fragments)):
            raw_frag = fragments[lig_rank]

            # Extract slot assignments with actual atom indices derived from
            # the SMILES string position immediately before each {slot} marker.
            for slot_match in OINInlineHandler.SLOT_REGEX.finditer(raw_frag):
                slot = int(slot_match.group(1))
                atom_idx = _count_smiles_atoms_before(raw_frag, slot_match.start())
                winding = slot_match.group(2)
                if winding == "^":
                    winding = ">"
                vector_data.append(SlotAssignment(lig_rank, atom_idx, slot, winding))

            # Strip {slot} markers while preserving all other content (@/@@, brackets, etc.)
            clean_frag = OINInlineHandler.SLOT_REGEX.sub("", raw_frag)
            clean_fragments.append(clean_frag)

        final_smiles = ".".join(clean_fragments)
        return final_smiles, geometry, vector_data
