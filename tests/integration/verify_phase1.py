import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.engine import OIN3DGenerator


def main():
    print("Verifying Phase 1: Rigid MVP")
    # OIN for Cisplatin (Square Planar like) - V2.3/V2.4 format
    oin = "[Pt].[Cl].[Cl].N.N |g:SPL|w:1.0:0;2.0:2;3.0:3;4.0:1|"

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
