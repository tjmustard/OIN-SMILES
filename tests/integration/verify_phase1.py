import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.generation.engine import OIN3DGenerator

def main():
    print("Verifying Phase 1: Rigid MVP")
    # Example 1: Cisplatin
    oin = "[Pt].[Cl].[Cl].[NH3].[NH3] |v:0.1:0.774,0.633,-0.000;0.2:-0.689,0.724,-0.000;0.3:-0.797,-0.604,-0.001;0.4:0.937,-0.350,-0.002|"
    
    generator = OIN3DGenerator()
    try:
        structure = generator.generate(oin)
        print("Structure generated successfully.")
        if structure:
            # Assuming structure has a write method or similar
            # structure.write_xyz("cisplatin_gen.xyz")
            print(f"Structure type: {type(structure)}")
    except Exception as e:
        print(f"Generation failed: {e}")

if __name__ == "__main__":
    main()
