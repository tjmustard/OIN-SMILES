import sys
import argparse
import pathlib

from oinsmiles.core.translator import XYZToSMILES
from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generation.molassembler_adapter import MolassemblerTimeoutError

def main():
    parser = argparse.ArgumentParser(description="OIN-SMILES Command Line Interface")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Command: xyz2oin
    parser_xyz2oin = subparsers.add_parser("xyz2oin", help="Convert XYZ file to OIN-SMILES")
    parser_xyz2oin.add_argument("path", type=str, help="Path to the XYZ file")

    # Command: oin2xyz
    parser_oin2xyz = subparsers.add_parser("oin2xyz", help="Convert OIN-SMILES to XYZ block")
    parser_oin2xyz.add_argument("oin", type=str, help="OIN-SMILES string to convert")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(0)

    if args.command == "xyz2oin":
        if not pathlib.Path(args.path).exists():
            print(f"FileNotFoundError: {args.path}", file=sys.stderr)
            sys.exit(1)
        
        try:
            translator = XYZToSMILES()
            oin_string = translator.convert(args.path)
            print(oin_string)
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

    elif args.command == "oin2xyz":
        try:
            generator = OIN3DGenerator()
            xyz_block = generator.generate(args.oin)
            print(xyz_block)
        except MolassemblerTimeoutError:
            print("Error: Molassembler timed out", file=sys.stderr)
            sys.exit(2)
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

if __name__ == "__main__":
    main()
