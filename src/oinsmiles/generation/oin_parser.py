from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

# --- V2.3 Templates for Parser Resolution ---
TEMPLATES = {
    "LIN": np.array([[0, 0, 1], [0, 0, -1]]),
    "TPL": np.array([[0, 1, 0], [0.8660254, -0.5, 0], [-0.8660254, -0.5, 0]]),
    "SPL": np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]),
    "SPY": np.array([[0, 0, 1], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]]),
    "TET": np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]),
    "TPY": np.array([[0, 0, 1], [0, 1, 0], [0.8660254, -0.5, 0], [-0.8660254, -0.5, 0]]),
    "TBP": np.array(
        [[0, 0, 1], [0, 0, -1], [0, 1, 0], [0.8660254, -0.5, 0], [-0.8660254, -0.5, 0]]
    ),
    "OCT": np.array([[0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]]),
    "PBP": np.array(
        [
            [0, 0, 1],
            [0, 0, -1],  # Axial
            [1, 0, 0],  # Eq 1 (0 deg)
            [0.30901699, 0.95105652, 0],  # Eq 2 (72 deg)
            [-0.80901699, 0.58778525, 0],  # Eq 3 (144 deg)
            [-0.80901699, -0.58778525, 0],  # Eq 4 (216 deg)
            [0.30901699, -0.95105652, 0],  # Eq 5 (288 deg)
        ]
    ),
    # CN 8 -- square antiprismatic. Same /sqrt(3) vectors as MetalloGen's
    # `8_squre_antiprismatic` and the encoder's TEMPLATE_SPECS["SQA"], so
    # slot s here maps to generator slot s (identity nearest-vector match).
    "SQA": np.array(
        [
            [-0.5773503, 0.5773503, 0.5773503],
            [0.5773503, 0.5773503, 0.5773503],
            [-0.5773503, -0.5773503, 0.5773503],
            [0.5773503, -0.5773503, 0.5773503],
            [-0.8141210, 0.0, -0.5773503],
            [0.0, -0.8141210, -0.5773503],
            [0.8141210, 0.0, -0.5773503],
            [0.0, 0.8141210, -0.5773503],
        ]
    ),
}


@dataclass
class OINVector:
    """A coordination vector for one binding atom, with its fragment location."""

    atom_idx: int
    vector: Tuple[float, float, float]
    fragment_idx: int
    atom_in_fragment_idx: int
    slot: int = -1
    winding: Optional[str] = None


@dataclass
class ParsedOIN:
    """Structured result of parsing an OIN string for 3D generation."""

    smiles: str
    fragments: List[str]
    metal_fragment_idx: int
    vectors: List[OINVector]
    original_oin: str
    geo_code: str = ""
    winding_by_slot: Dict[int, Optional[str]] = field(default_factory=dict)


class OINParser:
    """Parse an OIN string into a ParsedOIN for the 3D generation pipeline."""

    def parse(self, oin_string: str) -> ParsedOIN:
        """Parse an OIN string into a ParsedOIN structure."""
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
            # but vector_list is just (Rank, Slot). We need to map to vectors.

            smiles, geo_code, vector_data = OINInlineHandler.parse_inline_string(oin_string)
            fragments = smiles.split(".")
            metal_fragment_idx = 0  # Assumption

            tmpl_vectors = TEMPLATES.get(geo_code)
            vectors = []
            winding_by_slot = {}

            for sa in vector_data:
                # Universal channel: populated for every slot assignment,
                # regardless of whether a geometry template exists, so
                # template-less (NON/eta) paths never lose winding.
                # Guard prevents non-heading atoms (with winding=None) from clobbering
                # the heading atom's winding in multi-atom slots (e.g., eta rings).
                if sa.slot not in winding_by_slot or sa.winding is not None:
                    winding_by_slot[sa.slot] = sa.winding

                if tmpl_vectors is not None and sa.slot < len(tmpl_vectors):
                    resolved_vec = tmpl_vectors[sa.slot]

                    vectors.append(
                        OINVector(
                            atom_idx=-1,
                            vector=tuple(resolved_vec.tolist()),
                            fragment_idx=sa.lig_rank,
                            atom_in_fragment_idx=sa.atom_idx,
                            slot=sa.slot,
                            winding=sa.winding,
                        )
                    )

            return ParsedOIN(
                smiles=smiles,
                fragments=fragments,
                metal_fragment_idx=metal_fragment_idx,
                vectors=vectors,
                original_oin=oin_string,
                geo_code=geo_code,
                winding_by_slot=winding_by_slot,
            )

        # Standard / Legacy Parsing

        fragments = smiles.split(".")

        # Identify metal fragment (usually 0, but could check symbol if needed)
        metal_fragment_idx = 0

        # 1. Identify Geometry Template First
        geo_code = ""
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
                    if not item:
                        continue
                    try:
                        # item format "Rank.Idx:Slot"
                        if ":" not in item:
                            continue

                        indices_str, slot_str = item.split(":", 1)
                        slot_idx = int(slot_str)

                        if "." not in indices_str:
                            continue

                        frag_idx_str, atom_idx_str = indices_str.split(".")
                        frag_idx = int(frag_idx_str)
                        atom_in_frag_idx = int(atom_idx_str)

                        # Resolve Vector
                        if slot_idx >= len(tmpl_vectors):
                            continue  # Safety
                        resolved_vec = tmpl_vectors[slot_idx]

                        vectors.append(
                            OINVector(
                                atom_idx=-1,
                                vector=tuple(resolved_vec.tolist()),
                                fragment_idx=frag_idx,
                                atom_in_fragment_idx=atom_in_frag_idx,
                                slot=slot_idx,
                            )
                        )
                    except ValueError:
                        continue

        return ParsedOIN(
            smiles=smiles,
            fragments=fragments,
            metal_fragment_idx=metal_fragment_idx,
            vectors=vectors,
            original_oin=oin_string,
            geo_code=geo_code,
        )
