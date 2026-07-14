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
    from oinsmiles.generation.engine import OIN3DGenerator  # noqa: PLC0415

    # Opt-in knobs -> ff_params. --embed-threads != 1 engages the batched (not
    # byte-identical) embed; --optimize-workers overrides the parallel g-xTB worker
    # count (default: a safe number below the core total).
    ff_params: dict[str, int] = {}
    if args.embed_threads != 1:
        ff_params["embed_num_threads"] = args.embed_threads
    if args.optimize_workers > 0:
        ff_params["optimize_num_workers"] = args.optimize_workers
    try:
        result = OIN3DGenerator(
            optimizer=args.optimizer, seed=args.seed, ff_params=ff_params or None
        ).generate(args.oin)
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
    p_oin2xyz.add_argument(
        "--optimizer",
        type=str,
        default="g-xtb",
        help=(
            "Geometry optimizer (default: xtb, standard g-xTB). Use "
            "'mace-omol-0-extra-large-1024' or 'mace-omol25' for higher accuracy, "
            "or 'ff' for the fast FF-only path."
        ),
    )
    p_oin2xyz.add_argument(
        "--seed",
        type=int,
        default=42,
        help=(
            "Random seed for the metallogen ETKDG embed (default: 42, "
            "deterministic). Change it to sample a different reproducible "
            "conformer."
        ),
    )
    p_oin2xyz.add_argument(
        "--embed-threads",
        type=int,
        default=1,
        help=(
            "Opt-in parallel conformer embedding across N cores (0 = all cores). "
            "Default 1 keeps the serial, byte-identical embed. Any other value uses "
            "RDKit's batched EmbedMultipleConfs -- faster for large complexes but "
            "NOT byte-identical (it samples conformers differently)."
        ),
    )
    p_oin2xyz.add_argument(
        "--optimize-workers",
        type=int,
        default=0,
        help=(
            "Number of conformers to optimize concurrently with g-xTB/MACE. Each worker "
            "runs a single-threaded xtb, so this is a real speedup. 0 (default) uses a "
            "safe count below the machine's core total; set a positive integer to override."
        ),
    )
    p_oin2xyz.set_defaults(func=_cmd_oin2xyz)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
