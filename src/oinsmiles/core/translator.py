from .graph import TMCGraph, BondType
from ..oin.parser import OINParser
from ..oin.writer import OINWriter
# from .utils.xyz2mol import xyz2mol # Placeholder

class XYZToSMILES:
    def __init__(self):
        self.writer = OINWriter()

    def convert(self, xyz_file_path: str) -> str:
        """
        Converts an XYZ file to an OIN-SMILES string.
        """
        from ..utils.xyz2mol import read_xyz_file, xyz2mol
        from rdkit import Chem
        import numpy as np

    def convert(self, xyz_file_path: str) -> str:
        """
        Converts an XYZ file to an OIN-SMILES string.
        """
        from ..utils.xyz2mol import get_tmc_mol, get_oin_string
        
        # 1. Get TMC Mol and Coords (using updated xyz2mol)
        # We need to know the charge. 
        # The current signature of convert(xyz_file_path) doesn't accept charge.
        # We might need to guess it or default to 0, or update the signature.
        # For now, let's assume 0 or try to infer? 
        # xyz2mol requires charge.
        # Let's default to 0 and maybe allow passing it?
        # But convert signature is fixed?
        # Let's check if we can parse charge from file or just use 0.
        
        charge = 0 # Default
        
        # Check if xyz_file_path is a path object or string
        from pathlib import Path
        path = Path(xyz_file_path)
        
        try:
            tmc_mol, xyz_coords = get_tmc_mol(path, charge, with_stereo=False)
        except Exception as e:
            # Maybe try different charges? 
            # For now, just raise.
            raise ValueError(f"xyz2mol failed: {e}")

        # 2. Generate OIN
        oin_string = get_oin_string(tmc_mol, xyz_coords)
        
        return oin_string

class SMILESToXYZ:
    def __init__(self):
        self.parser = OINParser()

    def convert(self, oin_string: str) -> TMCGraph:
        """
        Converts an OIN-SMILES string to a TMCGraph (which can be written to XYZ).
        """
        smiles, tags = self.parser.parse(oin_string)
        
        graph = TMCGraph()
        
        # 1. Parse SMILES to get basic connectivity (using RDKit in real app)
        # For now, we just create a dummy graph if we can't parse SMILES without RDKit
        # In a real implementation, we would parse the SMILES to get all atoms.
        # Here we only populate atoms that have coordinates in the w tag.
        
        # 2. Apply Coordinates
        # Check for 'v' tag first, fallback to 'w'
        vector_tag = tags.get('v', tags.get('w', ''))
        
        if vector_tag:
            coords_with_idx = self.parser.parse_coordinates(vector_tag)
            # We need to add atoms. 
            # Since we don't have the full molecule from SMILES here (placeholder),
            # we'll just add atoms found in the tag.
            
            max_idx = 0
            if coords_with_idx:
                max_idx = max(idx for idx, x, y, z in coords_with_idx)
                
            # Initialize with dummy atoms up to max_idx
            for _ in range(max_idx + 1):
                graph.add_atom("X", (0.0, 0.0, 0.0))
                
            for idx, x, y, z in coords_with_idx:
                if idx < len(graph.atoms):
                    graph.atoms[idx].coords = (x, y, z)
                    # We could try to infer element from SMILES if we parsed it
        
        # 3. Apply Connectivity
        if 'v' in tags:
            # V1.4: v tag defines connectivity
            bonds = self.parser.parse_connectivity(tags['v'])
            for src, tgt in bonds:
                if src < len(graph.atoms) and tgt < len(graph.atoms):
                    graph.add_bond(src, tgt, BondType.DATIVE)
        elif 'd' in tags:
            # Legacy support if parser still supports it (it doesn't, so this branch is dead unless I revert parser)
            # We assume input is V1.4 or we fail gracefully
            pass

        return graph
