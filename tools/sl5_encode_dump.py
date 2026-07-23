#!/usr/bin/env python
"""Encode a molecule list with whatever oinsmiles is on PYTHONPATH; dump {mol: oin} JSON.

Run once with the main checkout's src and once with the branch src, then diff the two
JSONs: that is the correct SL5 byte-identity gate (branch vs main), independent of the
stale capstone `smiles_1`. In-process (one import) so it is fast; a per-mol signal alarm
guards against any residual slow case.
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import warnings

warnings.filterwarnings("ignore")


class _TO(Exception):
    pass


def _h(s, f):
    raise _TO()


def find_xyz(mol, dataset_dir):
    for sub in ("cat", "photo"):
        p = os.path.join(dataset_dir, sub, f"{mol}.xyz")
        if os.path.exists(p):
            return p
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mols-file", required=True, help="one molecule name per line")
    ap.add_argument("--dataset-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-mol-timeout", type=int, default=120)
    args = ap.parse_args()

    from oinsmiles import XYZToSMILES

    signal.signal(signal.SIGALRM, _h)
    mols = [ln.strip() for ln in open(args.mols_file) if ln.strip()]
    out = {}
    for i, mol in enumerate(mols):
        xyz = find_xyz(mol, args.dataset_dir)
        if not xyz:
            out[mol] = "__NO_FILE__"
            continue
        signal.alarm(args.per_mol_timeout)
        try:
            out[mol] = XYZToSMILES().convert(xyz)
        except _TO:
            out[mol] = "__TIMEOUT__"
        except Exception as e:  # noqa: BLE001
            out[mol] = f"__ERR__:{type(e).__name__}"
        finally:
            signal.alarm(0)
        # Write incrementally so a kill under load does not lose the whole pass.
        if (i + 1) % 10 == 0 or i + 1 == len(mols):
            json.dump(out, open(args.out, "w"), indent=0)
            print(f"  ...{i + 1}/{len(mols)}", flush=True)
    json.dump(out, open(args.out, "w"), indent=0)
    print(f"wrote {len(out)} -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
