"""Measure how much of AC2BO's cost is *redundant* recomputation.

Deterministic, contention-robust: counts distinct argument keys vs total calls for
the two pure functions the encode profile identified as dominant --
``perception_core.AC2BO`` and ``perception_core.get_UA_pairs``. If distinct << total,
memoization removes the difference and is byte-identical by construction.

Usage:
    PYTHONPATH=src .venv/bin/python tools/perf_encode_redundancy.py \
        --dataset tmCAT-tmPHOTO_xyz_dataset --molecule QIDKUL_comp_0
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


def find_xyz(dataset_dir: str, molecule: str) -> str:
    if os.path.exists(molecule):
        return molecule
    for sub in ("cat", "photo", "regression_inputs", "cohort-v0.4.5-5k", "."):
        p = os.path.join(dataset_dir, sub, f"{molecule}.xyz")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{molecule} not found under {dataset_dir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="tmCAT-tmPHOTO_xyz_dataset")
    ap.add_argument("--molecule", required=True)
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args()

    xyz = find_xyz(args.dataset, args.molecule)
    with open(xyz) as f:
        natoms = int(f.readline().strip())

    from oinsmiles.utils import perception_core as loc

    stats = {
        "ac2bo_calls": 0,
        "ac2bo_distinct": Counter(),
        "uap_calls": 0,
        "uap_calls_expensive": 0,  # non-empty bond list -> runs matching
        "uap_distinct": Counter(),
        "uap_distinct_expensive": Counter(),
        "valence_combos_iterated": [],
        "ac2bo_wall": [],
    }

    orig_ac2bo = loc.AC2BO
    orig_uap = loc.get_UA_pairs
    orig_get_bonds = loc.get_bonds

    # count how many valence combos each AC2BO call actually iterates
    combo_counter = {"n": 0}
    orig_get_UA = loc.get_UA

    def counted_get_UA(*a, **kw):
        combo_counter["n"] += 1
        return orig_get_UA(*a, **kw)

    def counted_ac2bo(AC, atoms, charge, **kw):
        stats["ac2bo_calls"] += 1
        key = hashlib.sha256(
            AC.tobytes() + repr((list(atoms), charge, sorted(kw.items()))).encode()
        ).hexdigest()[:16]
        stats["ac2bo_distinct"][key] += 1
        combo_counter["n"] = 0
        t0 = time.perf_counter()
        try:
            return orig_ac2bo(AC, atoms, charge, **kw)
        finally:
            dt = time.perf_counter() - t0
            stats["ac2bo_wall"].append(round(dt, 3))
            stats["valence_combos_iterated"].append(combo_counter["n"])

    def counted_uap(UA, AC, DU, use_graph=True):
        stats["uap_calls"] += 1
        key = repr((tuple(UA), tuple(DU), use_graph))
        stats["uap_distinct"][key] += 1
        # replicate the cheap early-out test without running it twice
        expensive = bool(UA) and len(orig_get_bonds(UA, AC)) > 0
        if expensive:
            stats["uap_calls_expensive"] += 1
            stats["uap_distinct_expensive"][key] += 1
        return orig_uap(UA, AC, DU, use_graph=use_graph)

    loc.AC2BO = counted_ac2bo
    loc.get_UA_pairs = counted_uap
    loc.get_UA = counted_get_UA
    try:
        from oinsmiles import XYZToSMILES

        t0 = time.perf_counter()
        oin = XYZToSMILES().convert(xyz)
        wall = time.perf_counter() - t0
    finally:
        loc.AC2BO = orig_ac2bo
        loc.get_UA_pairs = orig_uap
        loc.get_UA = orig_get_UA

    out = {
        "molecule": args.molecule,
        "atoms": natoms,
        "wall_s": round(wall, 2),
        "oin": oin,
        "loadavg": os.getloadavg(),
        "AC2BO": {
            "calls": stats["ac2bo_calls"],
            "distinct_args": len(stats["ac2bo_distinct"]),
            "per_key_counts": sorted(stats["ac2bo_distinct"].values(), reverse=True),
            "wall_per_call": stats["ac2bo_wall"],
            "valence_combos_iterated": stats["valence_combos_iterated"],
        },
        "get_UA_pairs": {
            "calls": stats["uap_calls"],
            "distinct_keys": len(stats["uap_distinct"]),
            "calls_with_matching": stats["uap_calls_expensive"],
            "distinct_keys_with_matching": len(stats["uap_distinct_expensive"]),
            "top_repeat_counts": sorted(stats["uap_distinct_expensive"].values(), reverse=True)[
                :10
            ],
        },
    }
    print(json.dumps(out, indent=2, default=str))
    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump(out, f, indent=2, default=str)


if __name__ == "__main__":
    main()
