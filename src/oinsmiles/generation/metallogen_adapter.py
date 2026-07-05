from rdkit import Chem
from .oin_parser import OINParser
from ..generator3d import generate_3d_structures
from .molassembler_adapter import GeneratedStructure

OIN_TO_METALLOGEN_GEO = {
    "LIN": "2_linear",
    "TPL": "3_trigonal_planar",
    "SQP": "4_square_planar",
    "SPL": "4_square_planar",
    "TET": "4_tetrahedral",
    "TPY": "4_trigonal_pyramidal", 
    "SPY": "5_square_pyramidal",
    "TBP": "5_trigonal_bipyramidal",
    "OCT": "6_octahedral",
    "PBP": "7_pentagonal_bipyramidal",
}

def convert_oin_to_msmiles(oin_string: str) -> str:
    parser = OINParser()
    parsed = parser.parse(oin_string)
    
    geo = OIN_TO_METALLOGEN_GEO.get(parsed.geo_code, "")
    if not geo:
        raise ValueError(f"Geometry code '{parsed.geo_code}' not supported by MetalloGen mapping.")
        
    from ..generator3d import globalvars
    import numpy as np
    metallogen_vectors = globalvars.known_geometries_vector_dict[geo]
    
    # Pre-allocate array for fragments based on MetalloGen slots
    num_slots = len(metallogen_vectors)
    ligand_parts = [None] * num_slots
    
    # Metal fragment (strip OIN annotations like _OCT)
    metal_frag = parsed.fragments[parsed.metal_fragment_idx]
    import re
    metal_frag = re.sub(r"_[A-Z0-9]+", "", metal_frag)
    metal_frag = re.sub(r"@SP[0-9]+", "", metal_frag)
    
    # Map ligands
    for i, frag_smiles in enumerate(parsed.fragments):
        if i == parsed.metal_fragment_idx:
            continue
            
        mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
        if mol is None:
            raise ValueError(f"Failed to parse fragment {i}: {frag_smiles}")
            
        frag_vectors = [v for v in parsed.vectors if v.fragment_idx == i]
        
        for v in frag_vectors:
            # Find the matching MetalloGen slot index
            target_vec = np.array([v.vector[0], v.vector[1], v.vector[2]])
            dists = np.linalg.norm(metallogen_vectors - target_vec, axis=1)
            mg_slot_idx = np.argmin(dists)
            
            # Map number starts at 1, so slot index + 1
            atom = mol.GetAtomWithIdx(v.atom_in_fragment_idx)
            atom.SetAtomMapNum(int(mg_slot_idx + 1))
            
        mapped_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        # Store in the list using the FIRST slot index matched for this fragment (for monodentate, it's just 1 slot)
        # For multidentate, it will be placed at the index of its first attachment point
        first_mg_slot_idx = np.argmin(np.linalg.norm(metallogen_vectors - np.array([frag_vectors[0].vector[0], frag_vectors[0].vector[1], frag_vectors[0].vector[2]]), axis=1))
        ligand_parts[first_mg_slot_idx] = mapped_smiles
        
    # Filter out None and join
    msmiles_parts = [metal_frag] + [p for p in ligand_parts if p is not None]
    m_smiles = "|".join(msmiles_parts) + f"|{geo}"
    return m_smiles

class OIN3DGeneratorMetallogen:
    def __init__(self, timeout: int = 60, ensemble_size: int = 1, dg_strategy: str = "single"):
        self.ensemble_size = ensemble_size

    def generate(self, oin_string: str) -> GeneratedStructure:
        msmiles = convert_oin_to_msmiles(oin_string)
        print(f"DEBUG: Converted OIN '{oin_string}' -> m-SMILES '{msmiles}'")
        
        mols = generate_3d_structures(msmiles, num_conformers=self.ensemble_size)
        if not mols:
            raise ValueError("MetalloGen failed to generate any conformers")
            
        mol_obj = mols[0]
        
        from ..generator3d import get_xyz_string
        xyz_str = get_xyz_string(mol_obj)
        
        return GeneratedStructure(xyz=xyz_str, mol=None)
