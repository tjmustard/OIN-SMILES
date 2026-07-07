import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.engine import OIN3DGenerator


import argparse

def main():
    parser = argparse.ArgumentParser(description="Verify Phase 1: Rigid MVP")
    parser.add_argument("--optimizer", type=str, default=None, help="Optimizer to use")
    parser.add_argument("--ff-preset", type=str, default=None, help="FF preset")
    parser.add_argument("--ensemble-size", type=int, default=None, help="Number of conformers")
    parser.add_argument("--only", type=str, help="Ignored in phase1, added for bash compatibility")
    parser.add_argument("--cpu", action="store_true", help="Force CPU execution")
    args = parser.parse_args()

    if args.cpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    print("Verifying Phase 1: Rigid MVP")
    # OIN for Cisplatin (Square Planar like) - V2.3/V2.4 format
    oin = "[Pt].[Cl].[Cl].N.N |g:SPL|w:1.0:0;2.0:2;3.0:3;4.0:1|"

    generator = OIN3DGenerator(
        optimizer=args.optimizer, 
        ff_preset=args.ff_preset,
        ensemble_size=args.ensemble_size
    )
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
