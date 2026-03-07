import re

def extract_ligand_smiles(oin_string: str) -> str:
    """
    Strips the metal fragment and slot markers from an OIN string.
    Returns the ligand SMILES for CIP oracle testing.
    Test-only utility.
    """
    # 1. Remove geometry tag and vectors if present (Zone C)
    if "|" in oin_string:
        oin_string = oin_string.split("|")[0]
        
    # 2. Split by dots
    frags = oin_string.split(".")
    
    # 3. Identify and remove fragments that look like metals [Pt_SPL] or [Pd]
    ligand_frags = []
    for f in frags:
        if "[" in f and "_" in f and "]" in f:
            # Likely [M_GEO]
            continue
        if f.startswith("[") and f.endswith("]") and len(f) <= 4:
            # Likely single atom metal [Pd]
            continue
        
        # 4. Strip slot markers like {0}, {1}
        f_clean = re.sub(r'\{\d+\}', '', f)
        ligand_frags.append(f_clean)
        
    return ".".join(ligand_frags)
