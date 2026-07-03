"""Command-line interface for OIN-SMILES.

Entry point registered as ``oin-smiles`` in pyproject.toml.

Subcommands
-----------
xyz2oin <path>     Convert an XYZ file to an OIN-SMILES string (stdout).
oin2xyz <oin>      Generate a 3D XYZ block from an OIN-SMILES string (stdout).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _cmd_xyz2oin(args: argparse.Namespace) -> None:
    path = Path(args.path)
    if not path.exists():
        print(f"FileNotFoundError: {args.path}", file=sys.stderr)
        sys.exit(1)

    from oinsmiles import XYZToSMILES  # noqa: PLC0415

    try:
        oin_string = XYZToSMILES().convert(str(path))
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    print(oin_string)


def _cmd_oin2xyz(args: argparse.Namespace) -> None:
    from oinsmiles.generation.engine import (  # noqa: PLC0415
        MolassemblerTimeoutError,
        OIN3DGenerator,
    )

    try:
        result = OIN3DGenerator().generate(args.oin)
    except MolassemblerTimeoutError:
        print("Error: Molassembler timed out", file=sys.stderr)
        sys.exit(2)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    print(result.xyz)


def main() -> None:
    """Entry point for the ``oin-smiles`` command-line interface."""
    parser = argparse.ArgumentParser(
        prog="oin-smiles",
        description="OIN-SMILES toolkit: convert between XYZ structures and OIN-SMILES strings.",
    )
    subparsers = parser.add_subparsers(dest="command", metavar="subcommand")
    subparsers.required = True

    # xyz2oin subcommand
    p_xyz2oin = subparsers.add_parser(
        "xyz2oin",
        help="Convert an XYZ file to an OIN-SMILES string.",
        description="Read an XYZ file and print the corresponding OIN-SMILES string to stdout.",
    )
    p_xyz2oin.add_argument("path", type=str, help="Path to the input XYZ file.")
    p_xyz2oin.set_defaults(func=_cmd_xyz2oin)

    # oin2xyz subcommand
    p_oin2xyz = subparsers.add_parser(
        "oin2xyz",
        help="Generate a 3D XYZ block from an OIN-SMILES string.",
        description=(
            "Generate a 3D conformer from an OIN-SMILES string and print the XYZ block to stdout."
        ),
    )
    p_oin2xyz.add_argument("oin", type=str, help="OIN-SMILES string.")
    p_oin2xyz.set_defaults(func=_cmd_oin2xyz)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
