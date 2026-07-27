"""How coarse a memo key can AC2BO's inner loop tolerate, and how much does it collapse?

Structural claims this probe tests (all read off the source, then measured):

1. ``get_UA_pairs(UA, AC, DU)`` reads DU **only** through the predicate ``du > 1``
   (it allocates a virtual matching node per such atom). So its output is a function of
   ``(AC, tuple(UA), tuple(du > 1 for du in DU))`` -- strictly coarser than its arguments.
2. ``get_bonds(UA, AC)`` is a function of ``(AC, tuple(UA))`` alone.
3. ``get_BO(AC, UA, DU, valences, UA_pairs)`` is a function of ``(AC, tuple(valences))``
   -- UA/DU/UA_pairs are all derived from those two.
4. ``_select_lig_mol``'s charge/carbene ladder calls ``AC2BO`` repeatedly on the *same*
   AC, so candidates recur across arms.

Prints total vs distinct for each key, i.e. the exact upper bound on what memoization
removes. Deterministic; unaffected by host load.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

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
    from oinsmiles.utils import perception_core as loc

    S = {
        "uap_total": 0,
        "uap_expensive": 0,
        "uap_key_args": set(),
        "uap_key_coarse": set(),
        "bonds_total": 0,
        "bonds_keys": set(),
        "getbo_total": 0,
        "getbo_keys": set(),
        "ac2bo": [],
    }
    ac_id = {"cur": None}

    orig_uap = loc.get_UA_pairs
    orig_bonds = loc.get_bonds
    orig_getbo = loc.get_BO
    orig_ac2bo = loc.AC2BO

    def ac_key(AC):
        return hashlib.sha256(AC.tobytes()).hexdigest()[:12]

    def w_ac2bo(AC, atoms, charge, **kw):
        ac_id["cur"] = ac_key(AC)
        t0 = time.perf_counter()
        before = (S["uap_expensive"], S["getbo_total"])
        try:
            return orig_ac2bo(AC, atoms, charge, **kw)
        finally:
            S["ac2bo"].append(
                {
                    "ac": ac_id["cur"],
                    "charge": charge,
                    "allow_carbenes": kw.get("allow_carbenes", True),
                    "wall": round(time.perf_counter() - t0, 2),
                    "uap_expensive": S["uap_expensive"] - before[0],
                    "get_BO": S["getbo_total"] - before[1],
                }
            )

    def w_uap(UA, AC, DU, use_graph=True):
        S["uap_total"] += 1
        a = ac_id["cur"] or ac_key(AC)
        S["uap_key_args"].add((a, tuple(UA), tuple(DU), use_graph))
        S["uap_key_coarse"].add((a, tuple(UA), tuple(d > 1 for d in DU), use_graph))
        if UA and len(orig_bonds(UA, AC)) > 0:
            S["uap_expensive"] += 1
        return orig_uap(UA, AC, DU, use_graph=use_graph)

    def w_bonds(UA, AC):
        S["bonds_total"] += 1
        S["bonds_keys"].add((ac_id["cur"] or ac_key(AC), tuple(UA)))
        return orig_bonds(UA, AC)

    def w_getbo(AC, UA, DU, valences, UA_pairs, use_graph=True):
        S["getbo_total"] += 1
        S["getbo_keys"].add((ac_id["cur"] or ac_key(AC), tuple(valences)))
        return orig_getbo(AC, UA, DU, valences, UA_pairs, use_graph=use_graph)

    loc.AC2BO = w_ac2bo
    loc.get_UA_pairs = w_uap
    loc.get_bonds = w_bonds
    loc.get_BO = w_getbo
    try:
        from oinsmiles import XYZToSMILES

        t0 = time.perf_counter()
        oin = XYZToSMILES().convert(xyz)
        wall = time.perf_counter() - t0
    finally:
        loc.AC2BO = orig_ac2bo
        loc.get_UA_pairs = orig_uap
        loc.get_bonds = orig_bonds
        loc.get_BO = orig_getbo

    out = {
        "molecule": args.molecule,
        "wall_s": round(wall, 2),
        "loadavg": [round(x, 2) for x in os.getloadavg()],
        "oin_sha256": hashlib.sha256(oin.encode()).hexdigest(),
        "get_UA_pairs": {
            "total_calls": S["uap_total"],
            "calls_running_matching": S["uap_expensive"],
            "distinct_by_arguments": len(S["uap_key_args"]),
            "distinct_by_coarse_key(AC,UA,du>1)": len(S["uap_key_coarse"]),
        },
        "get_bonds": {
            "total_calls": S["bonds_total"],
            "distinct_by_(AC,UA)": len(S["bonds_keys"]),
        },
        "get_BO": {
            "total_calls": S["getbo_total"],
            "distinct_by_(AC,valences)": len(S["getbo_keys"]),
        },
        "AC2BO_calls": S["ac2bo"],
    }
    print(json.dumps(out, indent=2))
    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
