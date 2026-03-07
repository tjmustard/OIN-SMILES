from typing import Any, Dict, List
try:
    from architector.complex_construction import build_complex
    from architector.io_molecule import Molecule
except ImportError as e:
    print(f"DEBUG: Architector Import Failed in Wrapper: {e}")
    build_complex = None
    Molecule = None

import tempfile
import os
from rdkit import Chem
from rdkit.Chem import AllChem

class ArchitectorWrapper:
    def run(self, input_dict: Dict[str, Any]) -> Any:
        """
        Runs architector.build_complex with the provided dictionary.
        """
        if build_complex is None:
            raise ImportError("Architector build_complex has no information. Install Architector or check configuration.")



        # Check for debug dump
        if "parameters" in input_dict and "_debug_dump_path" in input_dict["parameters"]:
            dump_path = input_dict["parameters"].pop("_debug_dump_path")
            try:
                import pprint
                with open(dump_path, 'w') as f:
                    f.write(pprint.pformat(input_dict))
            except Exception as e:
                print(f"Failed to dump Architector inputDict: {e}")

        # Call build_complex
        # It returns an ordered dictionary of conformers
        complex_out = build_complex(input_dict)
        
        if not complex_out:
            raise ValueError("Architector failed to generate any complexes.")
            
        # Get the first key
        first_key = list(complex_out.keys())[0]
        return complex_out[first_key]
