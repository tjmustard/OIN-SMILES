#!/usr/bin/env python
"""Pinpoint where a slow XYZ->OIN encode is stuck.

A Python ``signal.alarm`` cannot interrupt a long C-level RDKit call, so it also
cannot tell us *which* C call is hanging. ``faulthandler.dump_traceback_later``
runs in a watchdog thread and dumps every thread's Python stack after N seconds --
the main thread's top Python frame is the one that entered the stuck C call.

Usage:
    PYTHONPATH=src python tools/sl5_profile_hang.py <MOL> --dataset-dir <abs> [--after 25]
"""

from __future__ import annotations

import argparse
import faulthandler
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")


def find_xyz(mol: str, dataset_dir: str) -> str | None:
    for sub in ("cat", "photo"):
        p = os.path.join(dataset_dir, sub, f"{mol}_comp_0.xyz")
        if os.path.exists(p):
            return p
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mol")
    ap.add_argument("--dataset-dir", required=True)
    ap.add_argument("--after", type=float, default=25.0)
    args = ap.parse_args()

    from oinsmiles import XYZToSMILES

    path = find_xyz(args.mol, args.dataset_dir)
    if path is None:
        print(f"{args.mol}: NO_FILE", flush=True)
        return

    # Dump all thread stacks after `after` seconds, repeat every `after` seconds,
    # so a genuine hang prints its stuck frame while still running.
    faulthandler.dump_traceback_later(args.after, repeat=True, file=sys.stderr)
    t0 = time.monotonic()
    try:
        oin = XYZToSMILES().convert(path)
        print(f"{args.mol}: OK in {time.monotonic() - t0:.1f}s -> {oin[:60]}", flush=True)
    except Exception as e:  # noqa: BLE001
        print(
            f"{args.mol}: {type(e).__name__} in {time.monotonic() - t0:.1f}s -> {str(e)[:120]}",
            flush=True,
        )
    finally:
        faulthandler.cancel_dump_traceback_later()


if __name__ == "__main__":
    main()
