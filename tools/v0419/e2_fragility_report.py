"""v0.4.19 E2 lane -- fragility before/after on a subset (or the cohort) of e_selfconsistency rows.

    python tools/v0419/e2_fragility_report.py --arm <arm.jsonl> [--shipped <shipped.jsonl>]
                                              [--sweep results-v0.4.18-release-sweep] [--names f]

Per molecule: the census's own fragility predicate (`attribution_table.attribute`'s rule, restated
here: base encoded AND any of renum0-2 / noise0-2 not byte_exact; KEY-fragile ignores
slot_renumber / rdkit_canonical), under the arm and under the shipped run, plus whether the arm's
BASE string differs from the sweep of record's `smiles_1` (the string the generator would be
handed). Prints the 2x2 and the per-axis split (renumber vs noise), lists the movers.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
FRAG = ("renum0", "renum1", "renum2", "noise0", "noise1", "noise2")
RENUM, NOISE = FRAG[:3], FRAG[3:]
STRING_ONLY = (None, "byte_exact", "key_equal/slot_renumber", "key_equal/rdkit_canonical")


def load(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.strip() and line.startswith("{"):  # skip the #DONE trailer
                r = json.loads(line)
                out[r["molecule"]] = r
    return out


def fragile(r, keys=FRAG, key_level=False):
    if not r or not r.get("base"):
        return None  # not encoded: the predicate does not apply
    rel = r.get("rel", {})
    bad = STRING_ONLY if key_level else (None, "byte_exact")
    return any(rel.get(k) not in bad for k in keys)


def sweep_smiles(sweep, mols):
    out = {}
    for m in mols:
        p = sweep / "individual_reports" / f"{m}.json"
        if p.exists():
            try:
                out[m] = json.load(open(p)).get("smiles_1")
            except Exception:  # noqa: BLE001
                out[m] = None
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--arm", type=Path, required=True)
    ap.add_argument(
        "--shipped",
        type=Path,
        default=MAIN / "results-v0.4.17-exactfold" / "e_selfconsistency_exact.jsonl",
        help="the census of record's C2 run",
    )
    ap.add_argument("--sweep", type=Path, default=MAIN / "results-v0.4.18-release-sweep")
    ap.add_argument("--names", type=Path, help="restrict to these molecules (one per line)")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    arm, shipped = load(args.arm), load(args.shipped)
    mols = sorted(arm)
    if args.names:
        keep = {ln.split()[0] for ln in args.names.read_text().splitlines() if ln.strip()}
        mols = [m for m in mols if m in keep]
    smi = sweep_smiles(args.sweep, mols)

    rows = []
    for m in mols:
        a, s = arm.get(m), shipped.get(m)
        rows.append(
            {
                "molecule": m,
                "arm_encoded": bool(a and a.get("base")),
                "shipped_encoded": bool(s and s.get("base")),
                "arm_fragile": fragile(a),
                "shipped_fragile": fragile(s),
                "arm_key_fragile": fragile(a, key_level=True),
                "shipped_key_fragile": fragile(s, key_level=True),
                "arm_renum_fragile": fragile(a, RENUM),
                "shipped_renum_fragile": fragile(s, RENUM),
                "arm_noise_fragile": fragile(a, NOISE),
                "shipped_noise_fragile": fragile(s, NOISE),
                "base_moved": bool(a and a.get("base") and smi.get(m) and a["base"] != smi[m]),
                "arm_error": (a or {}).get("base_error"),
                "arm_rel": (a or {}).get("rel"),
            }
        )

    n = len(rows)
    enc = sum(r["arm_encoded"] for r in rows)
    print(
        f"molecules {n}  arm encoded {enc}  shipped encoded {sum(r['shipped_encoded'] for r in rows)}"
    )
    two = Counter((r["shipped_fragile"], r["arm_fragile"]) for r in rows)
    print("fragile (shipped, arm):", dict(two))
    print(
        "  fixed  :",
        two[(True, False)],
        "  broken:",
        two[(False, True)],
        "  still :",
        two[(True, True)],
    )
    twok = Counter((r["shipped_key_fragile"], r["arm_key_fragile"]) for r in rows)
    print("KEY-fragile (shipped, arm):", dict(twok))
    for axis in ("renum", "noise"):
        c = Counter((r[f"shipped_{axis}_fragile"], r[f"arm_{axis}_fragile"]) for r in rows)
        print(
            f"  {axis:6s} fixed {c[(True, False)]}  broken {c[(False, True)]}  still {c[(True, True)]}"
        )
    print("base string moved vs sweep smiles_1:", sum(r["base_moved"] for r in rows), "of", enc)
    errs = Counter(str(r["arm_error"])[:60] for r in rows if not r["arm_encoded"])
    if errs:
        print("arm encode errors:", dict(errs))
    print("\nstill fragile under the arm (molecule: columns not byte_exact):")
    for r in rows:
        if r["arm_fragile"]:
            bad = {k: v for k, v in (r["arm_rel"] or {}).items() if k in FRAG and v != "byte_exact"}
            print(f"  {r['molecule']:16s} {bad}")
    print("\nnewly fragile under the arm:")
    for r in rows:
        if r["arm_fragile"] and r["shipped_fragile"] is False:
            print(f"  {r['molecule']}")
    if args.out:
        with open(args.out, "w") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")


if __name__ == "__main__":
    main()
