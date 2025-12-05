from dataclasses import dataclass
from typing import List, Tuple
from rdkit import Chem

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
        
        # Identify metal fragment (Fragment 0 usually)
        metal_fragment_idx = 0
        
        vectors = []
        
        for meta in metadata:
            if meta.startswith("v:"):
                # Format: v:MetalIdx.LigandIdx:x,y,z;...
                content = meta[2:]
                items = content.split(";")
                for item in items:
                    if not item: continue
                    try:
                        # item might be "0.1:0.774,0.633,-0.000"
                        if ":" not in item: continue
                        
                        indices_str, vec_str = item.split(":", 1)
                        # indices_str is "Metal.Ligand"
                        if "." in indices_str:
                            metal_idx_str, ligand_idx_str = indices_str.split(".")
                            ligand_idx = int(ligand_idx_str)
                            # We use ligand_idx as the primary target for the vector
                            atom_idx = ligand_idx
                        else:
                            # Fallback or invalid format
                            continue
                            
                        x, y, z = map(float, vec_str.split(","))
                        
                        # Map global atom_idx to fragment
                        frag_idx, atom_in_frag_idx = self._map_index(atom_idx, fragments)
                        
                        vectors.append(OINVector(
                            atom_idx=atom_idx,
                            vector=(x, y, z),
                            fragment_idx=frag_idx,
                            atom_in_fragment_idx=atom_in_frag_idx
                        ))
                    except ValueError:
                        continue # Handle malformed items
                        
        return ParsedOIN(
            smiles=smiles,
            fragments=fragments,
            metal_fragment_idx=metal_fragment_idx,
            vectors=vectors,
            original_oin=oin_string
        )

    def _map_index(self, global_idx: int, fragments: List[str]) -> Tuple[int, int]:
        current_idx = 0
        for i, frag in enumerate(fragments):
            mol = Chem.MolFromSmiles(frag)
            if not mol:
                # Try to sanitize? Or just assume it works?
                # If it fails, we can't count atoms.
                # For now, raise error.
                raise ValueError(f"Invalid SMILES fragment: {frag}")
            
            num_atoms = mol.GetNumAtoms()
            if global_idx < current_idx + num_atoms:
                return i, global_idx - current_idx
            current_idx += num_atoms
            
        raise ValueError(f"Atom index {global_idx} out of bounds (total atoms: {current_idx})")
