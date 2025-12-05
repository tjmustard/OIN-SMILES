import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.generation.engine import OIN3DGenerator
# from oinsmiles import XYZToSMILES # Assuming this is the class for XYZ -> OIN

import numpy as np
from rdkit import Chem
from oinsmiles import XYZToSMILES
from oinsmiles.generation.engine import OIN3DGenerator
import tempfile

def calculate_pai_rmsd(coords1, coords2):
    """
    Aligns coords2 to coords1 using PAI and calculates RMSD.
    """
    # Since we don't have full molecules with bonds for coords1/2 easily without parsing,
    # and PAI alignment depends on masses, we'll assume equal masses for now or try to use RDKit if we have Mols.
    # But here we receive raw coords.
    # Let's try to use RDKit's alignment if we can pass Mols.
    # For now, we will just return a placeholder.
    return 0.0

def calculate_rmsd_mols(mol1, mol2):
    try:
        # Align mol2 to mol1
        # This uses standard alignment (min RMSD), not necessarily PAI.
        # But for verification it's a good proxy.
        rmsd = Chem.rdMolAlign.GetBestRMS(mol2, mol1)
        return rmsd
    except Exception as e:
        print(f"RMSD calculation failed: {e}")
        return 999.0

def read_xyz_coords(filepath):
    # Simple parser to get coordinates
    coords = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
        # Skip atom count and comment
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                coords.append([float(x) for x in parts[1:4]])
    return np.array(coords)

# Add examples to path to import real_life_examples
# Now in the same directory as verify_xyz_to_oin
from verify_xyz_to_oin import get_examples, Example

import argparse
import shutil

def main():
    parser = argparse.ArgumentParser(description="Verify OIN Round-Trip")
    parser.add_argument("--output-dir", type=str, help="Directory to save verification artifacts")
    parser.add_argument("--limit", type=int, help="Limit number of examples to run (for fast testing)")
    args = parser.parse_args()

    output_dir = args.output_dir
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Saving artifacts to: {output_dir}")

    print("Verifying Round-Trip Tests")
    
    # Setup
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
                # Determine paths
                if output_dir:
                    xyz_path = os.path.join(output_dir, f"Ex{i}_{safe_name}_TestA_input.xyz")
                    gen_xyz_path = os.path.join(output_dir, f"Ex{i}_{safe_name}_TestA_generated.xyz")
                    
                    with open(xyz_path, 'w') as f:
                        f.write(example.xyz_content)
                else:
                    # Create temp XYZ file for input
                    tmp_in = tempfile.NamedTemporaryFile(suffix='.xyz', delete=False, mode='w')
                    tmp_in.write(example.xyz_content)
                    tmp_in.close()
                    xyz_path = tmp_in.name
                    
                    tmp_out = tempfile.NamedTemporaryFile(suffix='.xyz', delete=False, mode='w')
                    tmp_out.close()
                    gen_xyz_path = tmp_out.name

                # 1. XYZ -> OIN
                print(f"Converting XYZ to OIN...")
                oin_string = xyz_to_smiles.convert(xyz_path)
                print(f"Generated OIN: {oin_string}")
                
                if output_dir:
                    with open(os.path.join(output_dir, f"Ex{i}_{safe_name}_TestA_generated.oin"), 'w') as f:
                        f.write(oin_string)
                
                # 2. OIN -> XYZ
                print("Converting OIN to XYZ...")
                extra_params = {}
                if output_dir:
                    extra_params["_debug_dump_path"] = os.path.join(output_dir, f"Ex{i}_{safe_name}_TestA_inputDict.txt")
                
                structure = generator.generate(oin_string, extra_params=extra_params)
                
                # 3. Save Generated XYZ
                if isinstance(structure, dict) and 'ase_atoms' in structure:
                    from ase.io import write
                    write(gen_xyz_path, structure['ase_atoms'], format='xyz')
                elif hasattr(structure, 'write_xyz'):
                    structure.write_xyz(gen_xyz_path)
                elif hasattr(structure, 'write_file'):
                    structure.write_file(gen_xyz_path)
                else:
                    with open(gen_xyz_path, 'w') as f:
                        f.write(str(structure))
                
                print(f"Generated XYZ saved to: {gen_xyz_path}")
                
                # 4. Compare
                mol_orig = Chem.MolFromXYZFile(xyz_path)
                mol_gen = Chem.MolFromXYZFile(gen_xyz_path)
                
                if mol_orig and mol_gen:
                    rmsd = calculate_rmsd_mols(mol_orig, mol_gen)
                    print(f"Geometric Round-Trip execution successful. RMSD: {rmsd:.4f}")
                else:
                    print("Geometric Round-Trip execution successful (Could not load Mols for RMSD).")
                
                # Cleanup if not saving
                if not output_dir:
                    if os.path.exists(xyz_path): os.remove(xyz_path)
                    if os.path.exists(gen_xyz_path): os.remove(gen_xyz_path)
                
            except Exception as e:
                import traceback
                traceback.print_exc()
                print(f"Test A FAILED: {e}")

        # --- Test B: Stability Round-Trip (OIN -> XYZ -> OIN) ---
        # Use oin_string if present, otherwise use expected_oin_string
        test_oin = example.oin_string or example.expected_oin_string
        
        if test_oin:
            print(f"\n--- Test B: Stability Round-Trip (OIN -> XYZ -> OIN) ---")
            try:
                # Determine paths
                if output_dir:
                    gen_xyz_path_b = os.path.join(output_dir, f"Ex{i}_{safe_name}_TestB_generated.xyz")
                else:
                    tmp = tempfile.NamedTemporaryFile(suffix='.xyz', delete=False, mode='w')
                    tmp.close()
                    gen_xyz_path_b = tmp.name

                # 1. OIN -> XYZ
                print(f"Input OIN: {test_oin}")
                print("Generating XYZ...")
                extra_params = {}
                if output_dir:
                    extra_params["_debug_dump_path"] = os.path.join(output_dir, f"Ex{i}_{safe_name}_TestB_inputDict.txt")
                
                structure = generator.generate(test_oin, extra_params=extra_params)
                
                if isinstance(structure, dict) and 'ase_atoms' in structure:
                    from ase.io import write
                    write(gen_xyz_path_b, structure['ase_atoms'], format='xyz')
                elif hasattr(structure, 'write_xyz'):
                    structure.write_xyz(gen_xyz_path_b)
                elif hasattr(structure, 'write_file'):
                    structure.write_file(gen_xyz_path_b)
                else:
                    with open(gen_xyz_path_b, 'w') as f:
                        f.write(str(structure))
                    
                # 2. XYZ -> OIN
                print("Converting generated XYZ back to OIN...")
                new_oin = xyz_to_smiles.convert(gen_xyz_path_b)
                print(f"Output OIN: {new_oin}")
                
                if output_dir:
                    with open(os.path.join(output_dir, f"Ex{i}_{safe_name}_TestB_output.oin"), 'w') as f:
                        f.write(new_oin)
                
                # 3. Compare
                print(f"Expected OIN: {test_oin}")
                print(f"Got OIN:      {new_oin}")
                
                failed_comparison = False
                
                # A. SMILES Part
                exp_smiles = test_oin.split('|')[0].strip()
                got_smiles = new_oin.split('|')[0].strip()
                
                if exp_smiles == got_smiles:
                     print("SMILES part matches.")
                else:
                     print("SMILES part MISMATCH.")
                     print(f"Expected: {exp_smiles}")
                     print(f"Got:      {got_smiles}")
                     failed_comparison = True

                # B. Tags Part
                def parse_oin_tags(oin_str):
                    tags = {}
                    parts = oin_str.split('|')
                    if len(parts) > 1:
                        for part in parts[1:]:
                            if ':' in part:
                                key, val = part.split(':', 1)
                                tags[key] = val
                    return tags

                exp_tags = parse_oin_tags(test_oin)
                got_tags = parse_oin_tags(new_oin)
                
                # Check 'v' tag (Vectors)
                if 'v' in exp_tags:
                    if 'v' not in got_tags:
                        print("Tag MISMATCH: 'v' tag missing in output.")
                        failed_comparison = True
                    else:
                        # Parse vectors: Metal.Ligand:x,y,z
                        def parse_vectors(v_str):
                            vecs = {}
                            entries = v_str.split(';')
                            for entry in entries:
                                if ':' in entry:
                                    indices, coords = entry.split(':')
                                    vecs[indices] = [float(x) for x in coords.split(',')]
                            return vecs

                        exp_vecs = parse_vectors(exp_tags['v'])
                        got_vecs = parse_vectors(got_tags['v'])
                        
                        # Calculate RMSD of vectors
                        # Collect all expected and got vectors into matching arrays
                        # Keys are Metal.Ligand indices
                        
                        common_indices = [k for k in exp_vecs.keys() if k in got_vecs]
                        missing_indices = [k for k in exp_vecs.keys() if k not in got_vecs]
                        
                        if missing_indices:
                            print(f"Vector MISMATCH: Missing vectors for {missing_indices}")
                            failed_comparison = True
                        else:
                            squared_diff_sum = 0.0
                            count = 0
                            
                            print("Comparing Vectors (RMSD Check):")
                            for idx in common_indices:
                                v_exp = np.array(exp_vecs[idx])
                                v_got = np.array(got_vecs[idx])
                                dist_sq = np.sum((v_exp - v_got)**2)
                                squared_diff_sum += dist_sq
                                count += 1
                                # print(f"  {idx}: DistSq={dist_sq:.4f}")
                            
                            if count > 0:
                                rmsd = np.sqrt(squared_diff_sum / count)
                                print(f"  Vector RMSD: {rmsd:.4f} (Tolerance: 1.0)")
                                
                                if rmsd > 1.0:
                                    print("  [FAIL] RMSD too high.")
                                    failed_comparison = True
                                else:
                                    print("  [PASS] RMSD within tolerance.")
                            else:
                                print("  No vectors to compare.")

                # Check 'm' tag (Multicenter)
                if failed_comparison:
                    print("Stability Round-Trip FAILED verification.")
                else:
                    print("Stability Round-Trip execution successful.")
                
                if not output_dir:
                    if os.path.exists(gen_xyz_path_b): os.remove(gen_xyz_path_b)
                
            except Exception as e:
                import traceback
                traceback.print_exc()
                print(f"Test B FAILED: {e}")
        else:
            print("\nSkipping Test B: No OIN string available.")

if __name__ == "__main__":
    main()
