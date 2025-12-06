import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.generation.engine import OIN3DGenerator
import numpy as np
from rdkit import Chem
from oinsmiles import XYZToSMILES
import tempfile
import argparse
import shutil

def calculate_rmsd_mols(mol1, mol2):
    try:
        rmsd = Chem.rdMolAlign.GetBestRMS(mol2, mol1)
        return rmsd
    except Exception as e:
        print(f"RMSD calculation failed: {e}")
        return 999.0

from verify_xyz_to_oin import get_examples, Example

def main():
    parser = argparse.ArgumentParser(description="Verify OIN Round-Trip")
    parser.add_argument("--output-dir", type=str, help="Directory to save verification artifacts")
    parser.add_argument("--limit", type=int, help="Limit number of examples to run (for fast testing)")
    args = parser.parse_args()

    output_dir = args.output_dir
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Saving artifacts to: {output_dir}")

    print("Verifying Round-Trip Tests (V2.3 String Identity)")
    
    xyz_to_smiles = XYZToSMILES()
    generator = OIN3DGenerator()
    
    examples = get_examples()
    if args.limit:
        print(f"Limiting to first {args.limit} examples.")
        examples = examples[:args.limit]

    print(f"Loaded {len(examples)} examples.")
    
    for i, example in enumerate(examples, 1):
        print(f"\n==================================================")
        print(f"Running Round-Trip for Example {i}: {example.name}")
        print(f"==================================================")
        
        safe_name = "".join(x for x in example.name if x.isalnum() or x in (' ', '-', '_')).replace(' ', '_')
        
        # --- Test A: Geometric Round-Trip (XYZ -> OIN -> XYZ) ---
        if example.xyz_content:
            print(f"\n--- Test A: Geometric Round-Trip (XYZ -> OIN -> XYZ) ---")
            try:
                # Setup Paths
                if output_dir:
                    xyz_path = os.path.join(output_dir, f"Ex{i}_{safe_name}_TestA_input.xyz")
                    gen_xyz_path = os.path.join(output_dir, f"Ex{i}_{safe_name}_TestA_generated.xyz")
                    with open(xyz_path, 'w') as f: f.write(example.xyz_content)
                else:
                    tmp_in = tempfile.NamedTemporaryFile(suffix='.xyz', delete=False, mode='w')
                    tmp_in.write(example.xyz_content)
                    tmp_in.close()
                    xyz_path = tmp_in.name
                    tmp_out = tempfile.NamedTemporaryFile(suffix='.xyz', delete=False, mode='w')
                    tmp_out.close()
                    gen_xyz_path = tmp_out.name

                # 1. XYZ -> OIN
                oin_string = xyz_to_smiles.convert(xyz_path)
                print(f"Generated OIN: {oin_string}")
                
                if output_dir:
                    with open(os.path.join(output_dir, f"Ex{i}_{safe_name}_TestA_generated.oin"), 'w') as f:
                        f.write(oin_string)
                
                # 2. OIN -> XYZ
                extra_params = {}
                if output_dir:
                    extra_params["_debug_dump_path"] = os.path.join(output_dir, f"Ex{i}_{safe_name}_TestA_inputDict.txt")
                
                structure = generator.generate(oin_string, extra_params=extra_params)
                
                # Save XYZ
                if isinstance(structure, dict) and 'ase_atoms' in structure:
                    from ase.io import write
                    write(gen_xyz_path, structure['ase_atoms'], format='xyz')
                elif hasattr(structure, 'write_xyz'):
                    structure.write_xyz(gen_xyz_path)
                elif hasattr(structure, 'write_file'):
                    structure.write_file(gen_xyz_path)
                else:
                    with open(gen_xyz_path, 'w') as f: f.write(str(structure))
                
                # 3. RMSD check
                mol_orig = Chem.MolFromXYZFile(xyz_path)
                mol_gen = Chem.MolFromXYZFile(gen_xyz_path)
                if mol_orig and mol_gen:
                    rmsd = calculate_rmsd_mols(mol_orig, mol_gen)
                    print(f"Geometric Round-Trip execution successful. RMSD: {rmsd:.4f}")
                else:
                    print("Geometric Round-Trip execution successful (Result Check Skipped).")
                
                if not output_dir:
                    if os.path.exists(xyz_path): os.remove(xyz_path)
                    if os.path.exists(gen_xyz_path): os.remove(gen_xyz_path)
                
            except Exception as e:
                print(f"Test A FAILED: {e}")

        # --- Test B: Stability Round-Trip (OIN -> XYZ -> OIN) ---
        test_oin = example.oin_string or example.expected_oin_string
        if test_oin:
            print(f"\n--- Test B: Stability Round-Trip (OIN -> XYZ -> OIN) ---")
            try:
                if output_dir:
                    gen_xyz_path_b = os.path.join(output_dir, f"Ex{i}_{safe_name}_TestB_generated.xyz")
                else:
                    tmp = tempfile.NamedTemporaryFile(suffix='.xyz', delete=False, mode='w')
                    tmp.close()
                    gen_xyz_path_b = tmp.name

                print(f"Input OIN: {test_oin}")
                
                # 1. OIN -> XYZ
                extra_params = {}
                if output_dir:
                    extra_params["_debug_dump_path"] = os.path.join(output_dir, f"Ex{i}_{safe_name}_TestB_inputDict.txt")
                
                structure = generator.generate(test_oin, extra_params=extra_params)
                
                if isinstance(structure, dict) and 'ase_atoms' in structure:
                    from ase.io import write
                    write(gen_xyz_path_b, structure['ase_atoms'], format='xyz')
                elif hasattr(structure, 'write_xyz'): structure.write_xyz(gen_xyz_path_b)
                elif hasattr(structure, 'write_file'): structure.write_file(gen_xyz_path_b)
                else:
                    with open(gen_xyz_path_b, 'w') as f: f.write(str(structure))
                    
                # 2. XYZ -> OIN
                new_oin = xyz_to_smiles.convert(gen_xyz_path_b)
                print(f"Output OIN: {new_oin}")
                
                if output_dir:
                    with open(os.path.join(output_dir, f"Ex{i}_{safe_name}_TestB_output.oin"), 'w') as f:
                        f.write(new_oin)
                
                # 3. VERIFICATION (String Identity)
                print(f"Expected OIN: {test_oin}")
                print(f"Got OIN:      {new_oin}")
                
                if test_oin.strip() == new_oin.strip():
                    print("  [PASS] OIN Strings are identical.")
                    print("Stability Round-Trip execution successful.")
                else:
                    print("  [FAIL] OIN Strings match FAILED.")
                    failed = True
                    print("Stability Round-Trip FAILED verification.")

                if not output_dir:
                    if os.path.exists(gen_xyz_path_b): os.remove(gen_xyz_path_b)
            except Exception as e:
                print(f"Test B FAILED: {e}")

if __name__ == "__main__":
    main()
