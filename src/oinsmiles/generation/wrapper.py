from typing import Any, Dict, List
try:
    from architector.complex_construction import build_complex
    from architector.io_molecule import Molecule
except ImportError:
    build_complex = None
    Molecule = None

import tempfile
import os
from rdkit import Chem
from rdkit.Chem import AllChem

class ArchitectorWrapper:
    def run(self, metal: str, ligands: List[Dict[str, Any]], parameters: Dict[str, Any]) -> Any:
        """
        Runs architector.build_complex
        """
        if build_complex is None:
            raise ImportError("Architector is not installed or functional.")
            
        # Construct inputDict for Architector
        inputDict = {
            "core": {
                "metal": metal,
                "coreType": "user_core",
                "coordList": parameters.get("site_coords", [])
            },
            "ligands": ligands, # List of dicts with 'smiles' and 'coordList'
            "parameters": {
                "force_field": "UFF",
                "debug": True # Helpful for debugging
            }
        }
        
        # Update parameters (excluding site_coords from parameters dict if needed, but Architector might ignore extra keys)
        # We should remove site_coords from parameters to avoid confusion or errors if Architector checks keys
        params_copy = parameters.copy()
        if "site_coords" in params_copy:
            del params_copy["site_coords"]
            
        inputDict["parameters"].update(params_copy)
        
        # Check for debug dump
        if "_debug_dump_path" in inputDict["parameters"]:
            dump_path = inputDict["parameters"].pop("_debug_dump_path")
            try:
                import pprint
                with open(dump_path, 'w') as f:
                    f.write(pprint.pformat(inputDict))
            except Exception as e:
                print(f"Failed to dump Architector inputDict: {e}")
        
        # Call build_complex
        # It returns an ordered dictionary of conformers
        complex_out = build_complex(inputDict)
        
        # The output is a dict of conformers. We usually want the best one.
        # Keys are usually names like '1_1_1'.
        # We can return the first one or the whole dict.
        # Let's return the first conformer object (which should be a Molecule-like object or dict)
        
        if not complex_out:
            raise ValueError("Architector failed to generate any complexes.")
            
        # Get the first key
        first_key = list(complex_out.keys())[0]
        return complex_out[first_key]
