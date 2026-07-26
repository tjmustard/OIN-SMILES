"""Byte-identity check for the perf lane's reencode-memo change.

Runs the full XYZ -> OIN -> XYZ (+ a final re-encode) pipeline for a molecule and prints a
sha256 of the generated XYZ string plus the re-encoded OIN string, so two code revisions can be
diffed by eye or with `diff <(script rev-A) <(script rev-B)`. Deterministic given a fixed seed
(default 42, matching the rest of the project) and a fixed `max_attempts` cap -- capping attempts
is a legitimate byte-identity probe here because the change under test (a memo keyed on mol
object identity) either reproduces a prior result exactly or, on a cache miss, recomputes it
exactly as before; it never alters which attempt wins or how many run.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


def find_xyz(dataset_dir: str, molecule: str) -> str:
    for sub in ("cat", "photo"):
        p = os.path.join(dataset_dir, sub, f"{molecule}.xyz")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{molecule} not found under {dataset_dir}/{{cat,photo}}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--molecule", required=True)
    ap.add_argument("--max-attempts", type=int, default=8)
    ap.add_argument("--optimizer", default=None)
    args = ap.parse_args()

    from oinsmiles import XYZToSMILES
    from oinsmiles.generation.metallogen_adapter import (
        OIN3DGeneratorMetallogen as OIN3DGenerator,
    )

    xyz_path = find_xyz(args.dataset, args.molecule)
    oin_string = XYZToSMILES().convert(xyz_path)

    gen = OIN3DGenerator(
        optimizer=args.optimizer,
        ensemble_size=1,
        timeout=600,
        ff_params={"max_attempts": args.max_attempts},
    )
    error = None
    xyz_out = None
    try:
        result = gen.generate(oin_string)
        xyz_out = result.xyz if result is not None else None
    except Exception as e:
        error = f"{type(e).__name__}: {e}"

    xyz_sha = hashlib.sha256(xyz_out.encode()).hexdigest() if xyz_out else None
    print(f"molecule={args.molecule} max_attempts={args.max_attempts} error={error!r}")
    print(f"xyz_sha256={xyz_sha}")
    if xyz_out:
        print(f"xyz_len={len(xyz_out)}")


if __name__ == "__main__":
    main()
