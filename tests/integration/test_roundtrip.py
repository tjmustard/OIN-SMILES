"""
Integration test: XYZ -> OIN -> Architector Dict.
"""
import pytest
import os
from pathlib import Path
from oinsmiles.core.translator import XYZToSMILES
from oinsmiles.generation.oin_parser import OINParser as GenParser
from oinsmiles.generation.architector_adapter import ArchitectorAdapter

TEST_DIR = Path(__file__).parent

def test_ticp2me2_roundtrip():
    xyz_path = TEST_DIR / "TiCp2Me2.xyz"
    if not xyz_path.exists():
        pytest.skip("TiCp2Me2.xyz not found")
        
    print(f"Testing with {xyz_path}")
    
    # 1. XYZ -> OIN
    # Note: Requires xyz2mol functionality in utils
    try:
        translator = XYZToSMILES()
        oin_string = translator.convert(str(xyz_path))
        print(f"OIN: {oin_string}")
        
    except ImportError:
        pytest.skip("xyz2mol dependencies missing")
    except Exception as e:
        pytest.fail(f"Translation failed: {e}")

    # Check for tags
    # Inline format: [Ti_TET] ...
    assert "_TET" in oin_string or "g:TET" in oin_string
    
    # 2. OIN -> Dict
    gen_parser = GenParser()
    parsed_oin = gen_parser.parse(oin_string)
    
    adapter = ArchitectorAdapter()
    input_dict = adapter.convert(parsed_oin)
    
    print("Generated InputDict Keys:", input_dict.keys())
    
    assert "core" in input_dict
    core = input_dict["core"]
    assert core["coreType"] == "user_core"
    
    # TiCp2Me2 has 2 Cp rings (5+5 atoms) and 2 Me (1+1 atoms).
    # Total vectors in core.coordList should be 5+5+1+1 = 12 vectors?
    # Or if Me is defined as Monodentate C, yes.
    # Provided xyz2mol detects them correctly.
    
    coord_list = core["coordList"]
    num_vecs = len(coord_list)
    print(f"Number of vectors: {num_vecs}")
    
    # Assert at least 10 vectors (for 2 Cp)
    assert num_vecs >= 10
    
    # Verify Ligand count
    # Ligands: Cp, Cp, Me, Me -> 4 ligands.
    ligands = input_dict["ligands"]
    assert len(ligands) >= 3 # At least Cp and Me
