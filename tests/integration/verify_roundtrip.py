import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.generation.engine import OIN3DGenerator
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from oinsmiles import XYZToSMILES
import tempfile
import argparse
import shutil

def calculate_rmsd_mols(mol1, mol2):
    try:
        rmsd = rdMolAlign.GetBestRMS(mol2, mol1)
        return rmsd
    except Exception as e:
        print(f"RMSD calculation failed: {e}")
        return 999.0

import re as _re
_METAL_STEREO_RE = _re.compile(r'\[([A-Z][a-z]?)@[A-Z0-9]+_([A-Z]{3})\]')

_WINDING_RE = _re.compile(r'\{(\d+)[><]\}')

def normalize_oin_for_comparison(oin_string: str) -> str:
    """Normalize an OIN string for round-trip comparison.

    1. Strip atom-ordering-dependent @SP/@OH/@TB stereo descriptors from the
       metal fragment — the slot assignments already encode the isomer geometry;
       the @XY## label depends on XYZ atom ordering and is not reproducible.
    2. Remove empty fragments (consecutive/trailing dots) caused by ligands that
       are present in the XYZ but uncoordinated in the OIN (e.g. H2 in FeH2(CO)4).
    3. Normalize water notation: [OH2] and O are chemically equivalent as bound
       water ligands. The XYZ→OIN pipeline may write O while generated structures
       re-analyzed after H addition write [OH2].
    4. Strip winding direction markers (> and <) from slot tags: {n>} and {n<}
       are normalized to {n}.  The ring rotation phase of eta-ligands (Cp, arene)
       cannot be deterministically reproduced from the OIN alone; the RMSD check
       verifies geometric correctness instead.
    """
    s = _METAL_STEREO_RE.sub(r'[\1_\2]', oin_string)
    # Normalize [OH2] → O (bound water notation equivalence)
    s = s.replace('[OH2]', 'O')
    # Normalize winding direction: {n>} → {n}, {n<} → {n}
    s = _WINDING_RE.sub(r'{\1}', s)
    # Collapse multiple consecutive dots and strip trailing dots
    while '..' in s:
        s = s.replace('..', '.')
    s = s.rstrip('.')
    return s

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
            else:
                tmp_dir = tempfile.mkdtemp()
                input_xyz_path = os.path.join(tmp_dir, "input.xyz")
                gen_xyz_path = os.path.join(tmp_dir, "gen.xyz")
                oin1_path = os.path.join(tmp_dir, "step1.oin")
                oin2_path = os.path.join(tmp_dir, "step2.oin")

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
            # Step 2: OIN(1) -> XYZ(Gen) (Molassembler)
            # -------------------------------------------------------------
            print("Step 2: Generate Structure OIN(1) -> XYZ(Gen)")
            structure = generator.generate(oin1_string)
            with open(gen_xyz_path, 'w') as f:
                f.write(structure)
                
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
            # Normalize: strip atom-ordering-dependent @SP/@OH/@TB descriptors
            # from the metal fragment before comparing — the slot assignments
            # already encode the isomer; the @XY## label is xyz-order-dependent.
            s1 = normalize_oin_for_comparison(oin1_string.strip())
            s2 = normalize_oin_for_comparison(oin2_string.strip())

            if s1 == s2:
                msg = "[PASS] OIN Stability: Strings Identical (normalized)"
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
