import sys
import os
import tempfile
import shutil
import argparse
from rdkit import Chem

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles import XYZToSMILES
from verify_xyz_to_oin import get_examples
from reporting import VerificationReporter
from rmsd_utils import calculate_tmc_rmsd

# Expected OIN strings (updated for V3.0 curly brace syntax)
FAC_EXPECTED = r"[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1"
MER_EXPECTED = r"[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{1}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{4}1"

def main():
    print("Verifying Ir(ppy)3 Complex Fixes (Unified Round-Trip)")
    
    # 1. Setup
    xyz_to_smiles = XYZToSMILES()
    generator = OIN3DGenerator()
    reporter = VerificationReporter("Ir(ppy)3 Verification Report")
    
    # 2. Get Examples and Filter
    all_examples = get_examples()
    target_names = ["fac-Ir(ppy)3", "mer-Ir(ppy)3"]
    examples = [ex for ex in all_examples if any(t in ex.name for t in target_names)]
    
    print(f"Found {len(examples)} target examples: {[ex.name for ex in examples]}")
    
    output_dir = "verification_artifacts_IR_TEST"
    if os.path.exists(output_dir): shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    
    for i, example in enumerate(examples, 1):
        test_name = f"{example.name} (Unified Round-Trip)"
        safe_name = "".join(x for x in example.name if x.isalnum() or x in (' ', '-', '_')).replace(' ', '_')
        print(f"\n==================================================")
        print(f"Running Example: {example.name}")
        print(f"==================================================")
        
        try:
            # Paths
            base_path = os.path.join(output_dir, f"{safe_name}")
            input_xyz_path = f"{base_path}_input.xyz"
            gen_xyz_path = f"{base_path}_generated.xyz"
            oin1_path = f"{base_path}_step1.oin"
            oin2_path = f"{base_path}_step2.oin"
            
            # Step 1: XYZ -> OIN(1)
            with open(input_xyz_path, 'w') as f: f.write(example.xyz_content)
            print("Step 1: Convert Input XYZ -> OIN(1)")
            oin1_string = xyz_to_smiles.convert(input_xyz_path)
            print(f"  OIN(1): {oin1_string}")
            with open(oin1_path, 'w') as f: f.write(oin1_string)
            
            # Step 2: OIN(1) -> XYZ(Gen)
            print("Step 2: Generate Structure OIN(1) -> XYZ(Gen)")
            structure = generator.generate(oin1_string)
            with open(gen_xyz_path, 'w') as f:
                f.write(structure.xyz)
                
            # Step 3: XYZ(Gen) -> OIN(2)
            print("Step 3: Convert XYZ(Gen) -> OIN(2)")
            oin2_string = xyz_to_smiles.convert(gen_xyz_path)
            print(f"  OIN(2): {oin2_string}")
            with open(oin2_path, 'w') as f: f.write(oin2_string)
            
            # Verification
            passed = True
            details = []
            
            # String Identity
            if oin1_string.strip() == oin2_string.strip():
                print("[PASS] Strings Identical")
                details.append("Strings Match")
            else:
                passed = False
                print("[FAIL] String Mismatch")
                details.append(f"Mismatch: <br>Exp: {oin1_string}<br>Got: {oin2_string}")
                
            # RMSD (coordination sphere only)
            mol_orig = Chem.MolFromXYZFile(input_xyz_path)
            mol_gen = Chem.MolFromXYZFile(gen_xyz_path)
            mol_gen_bonded = structure.mol if structure.mol is not None else mol_gen
            if mol_orig and mol_gen:
                rmsd = calculate_tmc_rmsd(mol_orig, mol_gen, mol2_bonded=mol_gen_bonded)
                print(f"  RMSD (coord sphere): {rmsd:.4f}")
                if rmsd < 1.0:
                    details.append(f"RMSD: {rmsd:.4f}")
                else:
                    passed = False
                    details.append(f"High RMSD: {rmsd:.4f}")
            else:
                details.append("RMSD Skipped (RDKit Error)")
                
            if passed:
                reporter.log_success(test_name, " | ".join(details))
            else:
                reporter.log_failure(test_name, "Validations Failed", got=" | ".join(details))
                
        except Exception as e:
            print(f"Failed: {e}")
            reporter.log_failure(test_name, f"Exception: {e}")
            
    reporter.print_summary()

if __name__ == "__main__":
    main()
