import os
import glob
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES
from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen as OIN3DGenerator

def run_stress_test():
    tmqm_dir = os.path.join(os.path.dirname(__file__), "tmQM")
    files = sorted(glob.glob(os.path.join(tmqm_dir, "*.xyz")))[:3]
    
    print(f"Found {len(files)} tmQM examples to stress test.")
    
    for xyz_file in files:
        name = os.path.basename(xyz_file)
        with open(xyz_file, "r") as f:
            xyz_content = f.read()
            
        print(f"\n======================================")
        print(f"--- Processing {name} ---")
        try:
            converter = XYZToSMILES()
            oin_string = converter.convert(xyz_file)
            print(f"OIN-SMILES: {oin_string}")
            
            # Generate 3D with MACE optimizer enabled
            generator = OIN3DGenerator()
            
            print(f"Attempting 3D generation + MLIP optimization...")
            result_xyzs = generator.generate(oin_string)
            
            if result_xyzs:
                print(f"Successfully generated {len(result_xyzs)} conformer(s) for {name}.")
            else:
                print(f"Failed to generate 3D for {name}.")
                
        except Exception as e:
            print(f"Error processing {name}: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    run_stress_test()
