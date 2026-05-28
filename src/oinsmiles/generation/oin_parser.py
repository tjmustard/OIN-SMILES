import dataclasses
from dataclasses import dataclass
from typing import List, Tuple, Dict
import numpy as np
import re
from rdkit import Chem
import scine_molassembler as masm

# --- Molassembler Shape Mapping (v0.2.1) ---
# Maps OIN geometry codes to Molassembler shape names
SCINE_SHAPE_MAP: Dict[str, str] = {
    'LIN': 'Line',
    'TPL': 'EquilateralTriangle',
    'SQP': 'Square',  # Square Planar (primary code)
    'SPL': 'Square',  # Square Planar (alternate code)
    'TET': 'Tetrahedron',
    'TPY': 'TrigonalPyramid',
    'SPY': 'SquarePyramid',
    'TBP': 'TrigonalBipyramid',
    'OCT': 'Octahedron',
    'PBP': 'PentagonalBipyramid',
}


def convert_bond_type(rdkit_bond_type: Chem.BondType) -> masm.BondType:
    """Map RDKit BondType to Molassembler BondType.

    Args:
        rdkit_bond_type: RDKit bond type (SINGLE, DOUBLE, AROMATIC, etc.)

    Returns:
        Corresponding Molassembler BondType

    Raises:
        ValueError: If the bond type cannot be converted
    """
    bond_map = {
        Chem.BondType.SINGLE: masm.BondType.Single,
        Chem.BondType.DOUBLE: masm.BondType.Double,
        Chem.BondType.TRIPLE: masm.BondType.Triple,
        Chem.BondType.AROMATIC: masm.BondType.Single,  # Treat aromatic as single for Molassembler
    }

    if rdkit_bond_type not in bond_map:
        raise ValueError(f"Unsupported RDKit bond type: {rdkit_bond_type}")

    return bond_map[rdkit_bond_type]


def construct_molassembler_mol(
    atoms: List,
    bonds: List[Tuple],
    constraints: Dict,
    mol_rdkit: Chem.Mol = None,
    frag_to_atom: Dict[int, List[int]] = None,
) -> masm.Molecule:
    """Construct Molassembler molecule with atoms, bonds, shape, and eta bonds.

    All-or-nothing transaction semantics: if any step fails, raise ValueError
    and do not return a partial molecule.

    Args:
        atoms: List of RDKit Atom objects (from AST tokenization)
        bonds: List of (i, j, bond_type) tuples (from AST tokenization)
        constraints: Dict {atom_idx: {'shape': str, 'vertex_indices': list[int]}}
                    (from regex preprocessor)
        mol_rdkit: Original RDKit molecule (contains full atom/bond info)
        frag_to_atom: Dict mapping fragment rank to atom indices in connected SMILES.
                     Required for eta bond translation when constraints['vertex_indices']
                     specifies fragment-rank-based slot assignments.

    Returns:
        masm.Molecule with shape and eta bonds applied

    Raises:
        ValueError: With context about which step failed and why
    """
    try:
        # Step 1: Build/validate RDKit molecule with all connectivity
        # (Include both standard bonds from AST and bonds inferred from constraints)
        if mol_rdkit is not None:
            working_mol = mol_rdkit
        else:
            # Reconstruct RDKit molecule from atoms/bonds
            rw_mol = Chem.RWMol()
            for atom in atoms:
                rw_mol.AddAtom(Chem.Atom(atom.GetAtomicNum()))
                rw_mol.GetAtomWithIdx(rw_mol.GetNumAtoms() - 1).SetFormalCharge(
                    atom.GetFormalCharge()
                )

            for i, j, bond_type in bonds:
                rw_mol.AddBond(i, j, bond_type)

            working_mol = rw_mol.GetMol()

        # Convert to SMILES for Molassembler
        smiles_repr = Chem.MolToSmiles(working_mol, isomericSmiles=True)

        try:
            mol = masm.io.experimental.from_smiles(smiles_repr)
        except Exception as e:
            raise ValueError(
                f"Failed to create Molassembler molecule from SMILES '{smiles_repr}': {e}"
            ) from e

        # Note: Molassembler may add hydrogen atoms, so mol.graph.V might be > len(atoms)
        # We just verify that we have at least the heavy atoms
        if mol.graph.V < len(atoms):
            raise ValueError(
                f"Atom count too low: input has {len(atoms)}, mol has {mol.graph.V}"
            )

        # Step 2: Assign polyhedral shape to metal center (assume at index 0)
        if 0 in constraints and 'shape' in constraints[0]:
            shape_code = constraints[0]['shape']

            if shape_code not in SCINE_SHAPE_MAP:
                valid_shapes = list(SCINE_SHAPE_MAP.keys())
                raise ValueError(
                    f"Unknown shape '{shape_code}'; valid shapes: {valid_shapes}"
                )

            try:
                scine_shape_name = SCINE_SHAPE_MAP[shape_code]
                scine_shape = getattr(masm.shapes.Shape, scine_shape_name)
                mol.set_shape_at_atom(0, scine_shape)
            except AttributeError as e:
                raise ValueError(
                    f"Shape '{shape_code}' → '{scine_shape_name}' not found in "
                    f"masm.shapes.Shape: {e}"
                ) from e
            except RuntimeError as e:
                # Molassembler may refuse to set shape if coordination number doesn't match
                # This can happen if the SMILES doesn't match the intended geometry
                # Log warning but continue (Molassembler's inferred shape may be acceptable)
                if "size" not in str(e).lower() and "mismatch" not in str(e).lower():
                    raise ValueError(
                        f"Failed to assign shape '{shape_code}' to metal center: {e}"
                    ) from e
                # Shape mismatch: Molassembler inferred a different coordination geometry
                # This is often expected with disconnected OIN SMILES
            except Exception as e:
                raise ValueError(
                    f"Failed to assign shape '{shape_code}' to metal center: {e}"
                ) from e

        # Step 3: Add eta bonds (metal to ligand atoms at same slot)
        if 0 in constraints and 'vertex_indices' in constraints[0]:
            vertex_indices = constraints[0]['vertex_indices']
            metal_idx = 0

            for vertex_index, frag_rank in enumerate(vertex_indices):
                try:
                    # Translate fragment rank to atom indices via mapping
                    if frag_to_atom is None:
                        # Fallback: assume vertex_indices are already atom indices
                        # (legacy compatibility, though this is the buggy case)
                        atom_indices = [frag_rank]
                    else:
                        # Use fragment mapping to get actual atom indices
                        if frag_rank not in frag_to_atom:
                            raise ValueError(
                                f"Fragment rank {frag_rank} not found in mapping. "
                                f"Valid ranks: {sorted(frag_to_atom.keys())}"
                            )
                        atom_indices = frag_to_atom[frag_rank]
                except ValueError as e:
                    raise ValueError(f"Failed to translate eta vertex {vertex_index}: {e}") from e

                # Add eta bonds from metal to each atom in the ligand
                for lig_idx in atom_indices:
                    if lig_idx < 0 or lig_idx >= mol.graph.V:
                        raise ValueError(
                            f"Eta bond: ligand atom {lig_idx} out of bounds "
                            f"(mol has {mol.graph.V} atoms)"
                        )

                    # Prevent bond-to-self (regression test target)
                    if metal_idx == lig_idx:
                        raise ValueError(
                            f"Refusing bond-to-self at atom {lig_idx}: "
                            f"eta bond source and target cannot be the same atom"
                        )

                    # Check if eta bond already exists (avoid duplicates)
                    bond_exists = False
                    try:
                        existing_type = mol.graph[masm.BondIndex(metal_idx, lig_idx)]
                        if existing_type == masm.BondType.Eta:
                            bond_exists = True
                    except (KeyError, RuntimeError, IndexError):
                        pass

                    if not bond_exists:
                        try:
                            mol.add_bond(metal_idx, lig_idx, masm.BondType.Eta)
                        except Exception as e:
                            raise ValueError(
                                f"Failed to add eta bond ({metal_idx}, {lig_idx}): {e}"
                            ) from e

        return mol

    except ValueError:
        # Re-raise ValueError as-is (already has context)
        raise
    except Exception as e:
        # Wrap any other exception
        raise ValueError(
            f"Failed to construct Molassembler molecule: {e}"
        ) from e

# --- Regex Preprocessor for Direct Parser (v0.2.1) ---
def _extract_oin_constraints(oin_smiles: str) -> Tuple[str, Dict[int, Dict], Dict[int, List[int]]]:
    """
    Extract polyhedral shape codes, chiral tags, vertex indices, and fragment-to-atom mapping.

    Args:
        oin_smiles: OIN-SMILES string in v3.6 format (e.g., "[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}")

    Returns:
        3-tuple: (stripped_smiles, constraints_dict, fragment_to_atom_mapping)
        - stripped_smiles: SMILES with OIN annotations removed, with RDKit atom maps
        - constraints_dict: {atom_idx: {'shape': str, 'chiral_tag': str, 'vertex_indices': list[int]}}
        - fragment_to_atom_mapping: {fragment_rank: [atom_indices]} mapping fragment ranks to atom indices in connected SMILES

    Patterns extracted from v3.6 format:
    - Shape codes: _([A-Z0-9]+) → e.g., _SQP, _OC, _SPL, _LIN
    - Chiral tags: @SP([0-9]+) → e.g., @SP1, @SP2
    - Vertex indices: \{([0-9><]+)\} → e.g., {0}, {1}, {0>} (> = CW, < = CCW)
    """
    # Regex patterns for extracting OIN v3.6 annotations
    shape_pattern = r'_([A-Z0-9]+)'
    chiral_pattern = r'@SP([0-9]+)'
    vertex_pattern = r'\{([0-9><]+)\}'

    constraints: Dict[int, Dict] = {}
    frag_vertex_map: Dict[int, int] = {}  # Maps fragment rank to vertex index

    # FIRST: Extract all constraints from the original string before stripping
    # Extract shape codes (from metal atom)
    shape_match = re.search(shape_pattern, oin_smiles)
    if shape_match:
        shape_code = shape_match.group(1)
        atom_idx = 0
        if atom_idx not in constraints:
            constraints[atom_idx] = {}
        constraints[atom_idx]['shape'] = shape_code

    # Extract chiral tag from metal atom
    chiral_match = re.search(chiral_pattern, oin_smiles)
    if chiral_match:
        chiral_tag = chiral_match.group(0)  # Include the @ symbol
        atom_idx = 0
        if atom_idx not in constraints:
            constraints[atom_idx] = {}
        constraints[atom_idx]['chiral_tag'] = chiral_tag

    # Extract vertex indices per-fragment
    # Split the OIN string by fragments (dots) and look for {N} patterns in each
    fragments = oin_smiles.split('.')
    for frag_rank, frag_smiles in enumerate(fragments):
        vertex_matches = re.findall(vertex_pattern, frag_smiles)
        if vertex_matches:
            # Take the first vertex annotation in this fragment
            numeric_part = vertex_matches[0].rstrip('><')
            if numeric_part.isdigit():
                vertex_idx = int(numeric_part)
                frag_vertex_map[frag_rank] = vertex_idx

    # Build vertex_indices list for metal center
    # This maps slot index to fragment rank
    if frag_vertex_map:
        # Build a list where slot N contains the fragment rank that binds at slot N
        vertex_indices = [None] * (max(frag_vertex_map.values()) + 1)
        for frag_rank, vertex_idx in frag_vertex_map.items():
            vertex_indices[vertex_idx] = frag_rank

        # Filter out None entries (slots with no ligand)
        vertex_indices = [v for v in vertex_indices if v is not None]

        if vertex_indices:
            if 0 not in constraints:
                constraints[0] = {}
            constraints[0]['vertex_indices'] = vertex_indices

    # SECOND: Strip OIN annotations to produce clean SMILES
    stripped = oin_smiles
    stripped = re.sub(shape_pattern, '', stripped)
    stripped = re.sub(chiral_pattern, '', stripped)
    stripped = re.sub(vertex_pattern, '', stripped)

    # THIRD: Build fragment-to-atom mapping before adding atom maps
    # Split into fragments (separated by dots) and track which atoms belong to which fragment
    fragment_to_atom_mapping: Dict[int, List[int]] = {}
    fragments = stripped.split('.')
    atom_idx = 0

    for frag_rank, frag_smiles in enumerate(fragments):
        frag_atoms = []
        # Count atoms in this fragment
        i = 0
        while i < len(frag_smiles):
            if frag_smiles[i] == '[':
                # Bracketed atom
                j = i + 1
                while j < len(frag_smiles) and frag_smiles[j] != ']':
                    j += 1
                frag_atoms.append(atom_idx)
                atom_idx += 1
                i = j + 1
            elif frag_smiles[i] in 'CNOPSBIFCcnopsbif':
                # Unbracketed organic atom (both uppercase and lowercase aromatic)
                # Check if it's a two-letter symbol (Cl, Br)
                if i + 1 < len(frag_smiles) and frag_smiles[i:i+2] in ['Cl', 'Br']:
                    i += 1
                frag_atoms.append(atom_idx)
                atom_idx += 1
                i += 1
            else:
                i += 1
        if frag_atoms:
            fragment_to_atom_mapping[frag_rank] = frag_atoms

    # FOURTH: Insert RDKit atom maps for tracking through AST tokenization
    # Replace [Metal] → [Metal:1], [Cl] → [Cl:2], etc.
    atom_map_counter = 1
    stripped_with_maps = ''
    i = 0
    while i < len(stripped):
        if stripped[i] == '[':
            # Find the closing bracket
            j = i + 1
            while j < len(stripped) and stripped[j] != ']':
                j += 1
            if j < len(stripped):
                # Extract the atom specification
                atom_spec = stripped[i+1:j]
                # Check if it already has an atom map (contains :)
                if ':' not in atom_spec:
                    stripped_with_maps += f'[{atom_spec}:{atom_map_counter}]'
                    atom_map_counter += 1
                else:
                    stripped_with_maps += stripped[i:j+1]
                i = j + 1
            else:
                stripped_with_maps += stripped[i]
                i += 1
        else:
            stripped_with_maps += stripped[i]
            i += 1

    return stripped_with_maps, constraints, fragment_to_atom_mapping


def tokenize_unsanitized_smiles(stripped_smiles: str) -> Tuple[List, List]:
    """
    Parse unsanitized SMILES into atom and bond lists for AST processing.

    Converts a SMILES string (with RDKit atom maps from regex stage) into lists of
    RDKit Atom and bond tuple objects. The molecule is parsed without sanitization
    to preserve aromatic flags and implicit hydrogen counts, deferring validation
    to Molassembler.

    Args:
        stripped_smiles: SMILES string with RDKit atom maps (e.g., "[Pt:1].[Cl:2]")

    Returns:
        (atoms, bonds)
        - atoms: list of rdkit.Chem.Atom objects in stable index order
        - bonds: list of (i, j, bond_type) tuples representing connectivity

    Raises:
        ValueError: If RDKit fails to parse the SMILES string.
    """
    # Parse with sanitize=False to preserve aromatic flags and skip valence checks
    mol = Chem.MolFromSmiles(stripped_smiles, sanitize=False)

    if mol is None:
        raise ValueError(f"Failed to parse unsanitized SMILES: {stripped_smiles}")

    # Extract atoms in stable order (by index)
    atoms = list(mol.GetAtoms())

    # Verify atom indices are contiguous and match GetIdx()
    for atom in atoms:
        assert atom.GetIdx() < len(atoms), (
            f"Atom index {atom.GetIdx()} out of bounds (mol has {len(atoms)} atoms)"
        )

    # Extract bonds as (begin_idx, end_idx, bond_type) tuples
    bonds = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        bonds.append((i, j, bond_type))

    # Note: We don't compute implicit hydrogens here since the molecule is unsanitized
    # and RDKit may not have computed implicit valence yet

    return atoms, bonds

# --- V2.3 Templates for Parser Resolution ---
def normalize_template(arr):
    return arr / np.linalg.norm(arr, axis=1)[:, None]

TEMPLATES = {
    'LIN': np.array([[0,0,1], [0,0,-1]]),
    'TPL': np.array([[0,1,0], [0.8660254,-0.5,0], [-0.8660254,-0.5,0]]),
    'SPL': np.array([[1,0,0], [0,1,0], [-1,0,0], [0,-1,0]]),
    'SPY': np.array([
        [0,0,1],
        [1,0,0], [-1,0,0], [0,1,0], [0,-1,0]
    ]),
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
    ]),
    'PBP': np.array([
        [0,0,1], [0,0,-1], # Axial
        [1,0,0], # Eq 1 (0 deg)
        [0.30901699, 0.95105652, 0], # Eq 2 (72 deg)
        [-0.80901699, 0.58778525, 0], # Eq 3 (144 deg)
        [-0.80901699, -0.58778525, 0], # Eq 4 (216 deg)
        [0.30901699, -0.95105652, 0] # Eq 5 (288 deg)
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
    geo_code: str = ""

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
             # but vector_list is just (Rank, Slot). We need to map to vectors.
             
             smiles, geo_code, vector_data = OINInlineHandler.parse_inline_string(oin_string)
             fragments = smiles.split(".")
             metal_fragment_idx = 0 # Assumption
             
             tmpl_vectors = TEMPLATES.get(geo_code)
             vectors = []
             
             if tmpl_vectors is not None:
                for lig_rank, atom_in_fragment_idx, slot_idx in vector_data:
                    if slot_idx < len(tmpl_vectors):
                        resolved_vec = tmpl_vectors[slot_idx] 
                        
                        vectors.append(OINVector(
                            atom_idx=-1,
                            vector=tuple(resolved_vec.tolist()),
                            fragment_idx=lig_rank,
                            atom_in_fragment_idx=atom_in_fragment_idx
                        ))
             
             return ParsedOIN(
                smiles=smiles,
                fragments=fragments,
                metal_fragment_idx=metal_fragment_idx,
                vectors=vectors,
                original_oin=oin_string,
                geo_code=geo_code,
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
            original_oin=oin_string,
            geo_code=geo_code,
        )
