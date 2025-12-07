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
from reporting import VerificationReporter

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
    
    reporter = VerificationReporter("Round-Trip Verification Report")
    
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
        
        # --- Unified Round-Trip Test ---
        # Flow: XYZ -> OIN(1) -> XYZ(Gen) -> OIN(2)
        # Checks: RMSD(XYZ, XYZ_Gen) < 1.0  AND  OIN(1) == OIN(2)
        
        test_name = f"{example.name} (Unified Round-Trip)"
        print(f"\n--- Unified Round-Trip (XYZ -> OIN -> XYZ -> OIN) ---")
        
        try:
            # 0. Setup Paths
            base_name = f"Ex{i}_{safe_name}"
            
            if output_dir:
                input_xyz_path = os.path.join(output_dir, f"{base_name}_input.xyz")
                gen_xyz_path = os.path.join(output_dir, f"{base_name}_generated.xyz")
                oin1_path = os.path.join(output_dir, f"{base_name}_step1.oin")
                oin2_path = os.path.join(output_dir, f"{base_name}_step2.oin")
                debug_dump_path = os.path.join(output_dir, f"{base_name}_inputDict.txt")
            else:
                tmp_dir = tempfile.mkdtemp()
                input_xyz_path = os.path.join(tmp_dir, "input.xyz")
                gen_xyz_path = os.path.join(tmp_dir, "gen.xyz")
                oin1_path = os.path.join(tmp_dir, "step1.oin")
                oin2_path = os.path.join(tmp_dir, "step2.oin")
                debug_dump_path = None # Or temp if needed

            # -------------------------------------------------------------
            # Step 1: START (Determine Input)
            # -------------------------------------------------------------
            start_oin = None
            
            if example.xyz_content:
                # Flow A: Start from XYZ
                with open(input_xyz_path, 'w') as f: f.write(example.xyz_content)
                
                # XYZ -> OIN(1)
                print("Step 1: Convert Input XYZ -> OIN(1)")
                oin1_string = xyz_to_smiles.convert(input_xyz_path)
                print(f"  OIN(1): {oin1_string}")
                
                if output_dir:
                    with open(oin1_path, 'w') as f: f.write(oin1_string)
                
                start_oin = oin1_string
            else:
                # Flow B: Start from OIN (if no XYZ provided)
                # Use example.oin_string or expected
                start_oin = example.oin_string or example.expected_oin_string
                if not start_oin:
                    print("Skipping: No XYZ content and no OIN string provided.")
                    reporter.log_failure(test_name, "No data provided")
                    continue
                
                print(f"Step 1: Start from provided OIN: {start_oin}")
                oin1_string = start_oin
                # No input XYZ path to compare RMSD against later (unless we generate one?)
                # We'll skip RMSD check if started from OIN.
                input_xyz_path = None

            # -------------------------------------------------------------
            # Step 2: OIN(1) -> XYZ(Gen) (Architector)
            # -------------------------------------------------------------
            print("Step 2: Generate Structure OIN(1) -> XYZ(Gen)")
            extra_params = {}
            if debug_dump_path:
                extra_params["_debug_dump_path"] = debug_dump_path
            
            structure = generator.generate(oin1_string, extra_params=extra_params)
            
            # Save generated XYZ
            if isinstance(structure, dict) and 'ase_atoms' in structure:
                from ase.io import write
                write(gen_xyz_path, structure['ase_atoms'], format='xyz')
            elif hasattr(structure, 'write_xyz'):
                structure.write_xyz(gen_xyz_path)
            elif hasattr(structure, 'write_file'):
                structure.write_file(gen_xyz_path)
            else:
                with open(gen_xyz_path, 'w') as f: f.write(str(structure))
                
            # -------------------------------------------------------------
            # Step 3: XYZ(Gen) -> OIN(2)
            # -------------------------------------------------------------
            print("Step 3: Convert XYZ(Gen) -> OIN(2)")
            oin2_string = xyz_to_smiles.convert(gen_xyz_path)
            print(f"  OIN(2): {oin2_string}")
            
            if output_dir:
                with open(oin2_path, 'w') as f: f.write(oin2_string)

            # -------------------------------------------------------------
            # Step 4: Verification
            # -------------------------------------------------------------
            print("\n--- Verification Results ---")
            passed = True
            details = []
            
            # Check 1: String Identity (OIN 1 vs OIN 2)
            # Normalize
            s1 = oin1_string.strip()
            s2 = oin2_string.strip()
            
            if s1 == s2:
                msg = "[PASS] OIN Stability: Strings Identical"
                print(msg)
                details.append(msg)
            else:
                passed = False
                msg = f"[FAIL] OIN Stability: Mismatch"
                print(msg)
                print(f"  Expected: {s1}")
                print(f"  Got:      {s2}")
                details.append(f"{msg}<br>Exp: `{s1}`<br>Got: `{s2}`")

            # Check 2: Geometric Fidelity (RMSD) - Only if we started from XYZ
            if input_xyz_path and os.path.exists(input_xyz_path):
                mol_orig = Chem.MolFromXYZFile(input_xyz_path)
                mol_gen = Chem.MolFromXYZFile(gen_xyz_path)
                
                if mol_orig and mol_gen:
                    rmsd = calculate_rmsd_mols(mol_orig, mol_gen)
                    print(f"  RMSD Input vs Generated: {rmsd:.4f}")
                    
                    if rmsd < 1.0:
                        msg = f"[PASS] Geometry: RMSD {rmsd:.4f} < 1.0"
                        print(msg)
                        details.append(f"RMSD: {rmsd:.4f}")
                    else:
                        passed = False
                        msg = f"[FAIL] Geometry: High RMSD {rmsd:.4f}"
                        print(msg)
                        details.append(f"<b>High RMSD: {rmsd:.4f}</b>")
                else:
                    msg = "[WARN] RDKit failed to support RMSD calc (atom mismatch or parsing error)"
                    print(msg)
                    details.append(msg)
            else:
                details.append("(Skipped RMSD - No Input XYZ)")

            if passed:
                reporter.log_success(test_name, " | ".join(details))
            else:
                reporter.log_failure(test_name, "Validations Failed", got="<br>".join(details))

            # Cleanup
            if not output_dir:
                if 'tmp_dir' in locals() and os.path.exists(tmp_dir):
                    shutil.rmtree(tmp_dir)

        except Exception as e:
            print(f"Unified Test FAILED: {e}")
            import traceback
            traceback.print_exc()
            reporter.log_failure(test_name, f"Exception: {str(e)}")

    # Final Summary
    reporter.print_summary()

if __name__ == "__main__":
    main()
