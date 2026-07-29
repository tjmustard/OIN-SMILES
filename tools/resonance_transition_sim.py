#!/usr/bin/env python
"""What does OIN_RESONANCE_DONOR_FOLD move, against the v0.4.13 SHIPPED state? (v0.4.14 Lane 1)

WHY THIS IS NOT ``fold_transition_sim.py`` WITH A FLAG
=====================================================
That tool's baseline arm is a pure STRING operation over the frozen corpus, because the fold it
measured was being compared against a v0.4.8 encoder that had no veto. This lever is different
in a way that changes the shape of the measurement: the thing it must be compared against is
**v0.4.13's shipped default path**, which runs the parity veto and therefore needs COORDINATES
on both sides. Both arms here re-encode from disk.

    arm BEFORE   OIN_CANONICAL_DONOR_FOLD=1  OIN_FOLD_PARITY_VETO=1  OIN_RESONANCE_DONOR_FOLD=0
    arm AFTER    OIN_CANONICAL_DONOR_FOLD=1  OIN_FOLD_PARITY_VETO=1  OIN_RESONANCE_DONOR_FOLD=1

THE POPULATION IS EXACT, NOT A SAMPLE
=====================================
Only molecules whose emitted string moves under the widening can change bucket. If neither
string moves at the string level (with the fold held on), then ``s_fold`` is identical in both
arms, so the veto sees identical inputs, reaches an identical verdict and emits an identical
string. So restricting to the string-level movers loses nothing, and every other molecule keeps
its v0.4.13 bucket by construction.

THE DIRECTION THAT MATTERS MOST IS THE BAD ONE
==============================================
A widening can LOSE points, and the mechanism is specific enough to name: a molecule the base
fold moved and the veto ALLOWED may, with the widening on, fold further, collapse a mirror pair,
and be vetoed **as a whole** -- surrendering the v0.4.13 gain as well as the new one. That is
counted here as ``byte_exact -> anything`` and reported on its own line. A release that only
counts gains cannot see it.

THE CONTROL THAT MUST NOT BE SKIPPED
====================================
Re-encode fold-OFF/veto-OFF first and require it to reproduce the frozen string. A molecule that
fails is ``drift`` and is excluded, because its transition would be describing encoder movement
since v0.4.8 rather than this lever.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/resonance_transition_sim.py \\
        --sweep <MAIN>/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest \\
        --dataset <MAIN>/tmCAT-tmPHOTO_xyz_dataset/cat \\
        --dataset <MAIN>/tmCAT-tmPHOTO_xyz_dataset/photo \\
        --baseline-byte-exact 3794 --out-json resonance_transition.json

⚠ ``--dataset`` has NO default. ``fold_transition_sim.py``'s relative defaults resolved to
nothing from a worktree, every molecule landed in ``unavailable``, and it printed a refuted
number under a correct-looking heading.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

from roundtrip_bucket_report import classify  # noqa: E402

from oinsmiles.oin import fold_parity  # noqa: E402
from oinsmiles.oin.canonical_slots import canonicalize_oin_slots  # noqa: E402

FOLD = "OIN_CANONICAL_DONOR_FOLD"
VETO = "OIN_FOLD_PARITY_VETO"
RESO = "OIN_RESONANCE_DONOR_FOLD"


def _levers(fold, veto, reso):
    os.environ[FOLD] = "1" if fold else "0"
    os.environ[VETO] = "1" if veto else "0"
    os.environ[RESO] = "1" if reso else "0"


def load_rows(sweep, limit=None):
    rows = []
    for f in sorted(glob.glob(os.path.join(sweep, "individual_reports", "*.json"))):
        try:
            with open(f) as fh:
                rows.append(json.load(fh))
        except Exception:  # noqa: BLE001 -- an unreadable report is skipped, not fatal
            continue
        if limit and len(rows) >= limit:
            break
    return rows


def bucket_of(rep, s1, s2):
    shadow = dict(rep)
    shadow["smiles_1"] = s1
    shadow["smiles_2_indep"] = s2
    b, sub = classify(shadow, score="honest")
    return f"{b}/{sub}" if sub else b


def _find_input_xyz(molecule, roots):
    for root in roots:
        p = os.path.join(root, f"{molecule}.xyz")
        if os.path.exists(p):
            return p
    return None


def _encode(path):
    """Encode under the CURRENT lever settings. ``(string, veto outcome)``."""
    from oinsmiles import XYZToSMILES

    fold_parity._state.outcome = None
    try:
        s = XYZToSMILES().convert(path)
    except Exception:  # noqa: BLE001 -- an unencodable structure is excluded, not fatal
        return None, fold_parity.last_outcome()
    return s, fold_parity.last_outcome()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--dataset", action="append", required=True, help="ABSOLUTE input-xyz root")
    ap.add_argument("--limit", type=int)
    ap.add_argument(
        "--baseline-byte-exact",
        type=int,
        required=True,
        help="byte_exact count of the state being compared against (v0.4.13: 3794)",
    )
    ap.add_argument("--out-json")
    args = ap.parse_args()

    for root in args.dataset:
        if not os.path.isdir(root):
            sys.exit(f"🔴 --dataset {root} is not a directory")

    rows = load_rows(args.sweep, args.limit)
    print(f"loaded {len(rows)} reports from {args.sweep}", flush=True)

    # Population: string-level movers of the widening, with the fold held ON in both arms.
    movers = []
    for rep in rows:
        mol = rep.get("molecule")
        moved = False
        for field in ("smiles_1", "smiles_2_indep"):
            s = rep.get(field)
            if not s:
                continue
            try:
                _levers(fold=True, veto=False, reso=False)
                a = canonicalize_oin_slots(s)
                _levers(fold=True, veto=False, reso=True)
                b = canonicalize_oin_slots(s)
            except Exception:  # noqa: BLE001 -- an unfoldable string cannot move
                continue
            if a != b:
                moved = True
        if moved:
            movers.append(mol)
    print(f"{RESO} moves a string on {len(movers)} molecules -- re-encoding both arms", flush=True)

    by_mol = {r.get("molecule"): r for r in rows}
    struct_dir = os.path.join(args.sweep, "structures")
    per_mol, drift, unavailable = [], [], []

    for i, mol in enumerate(movers, 1):
        rep = by_mol[mol]
        in_xyz = _find_input_xyz(mol, args.dataset)
        gen_xyz = os.path.join(struct_dir, f"{mol}_generated.xyz")
        if not in_xyz or not os.path.exists(gen_xyz):
            unavailable.append(mol)
            continue

        _levers(fold=False, veto=False, reso=False)
        c1, _ = _encode(in_xyz)
        c2, _ = _encode(gen_xyz)
        if c1 != rep.get("smiles_1") or c2 != rep.get("smiles_2_indep"):
            drift.append(mol)
            continue

        _levers(fold=True, veto=True, reso=False)
        b1, ob1 = _encode(in_xyz)
        b2, ob2 = _encode(gen_xyz)
        _levers(fold=True, veto=True, reso=True)
        a1, oa1 = _encode(in_xyz)
        a2, oa2 = _encode(gen_xyz)

        per_mol.append(
            {
                "molecule": mol,
                "before": bucket_of(rep, b1, b2),
                "after": bucket_of(rep, a1, a2),
                "outcome_before": [ob1, ob2],
                "outcome_after": [oa1, oa2],
            }
        )
        if i % 20 == 0:
            print(f"  [{i}/{len(movers)}] {len(per_mol)} measured", flush=True)

    # 🔴 REFUSE rather than report. With every mover excluded this prints "0 gains, 0 losses",
    # which reads as "the lever is safe and worthless" -- a clean, plausible, wrong answer.
    if movers and not per_mol:
        sys.exit(
            f"\n🔴 REFUSING TO REPORT: measured 0 of {len(movers)} movers "
            f"({len(unavailable)} unavailable, {len(drift)} drift)."
        )

    changed = [r for r in per_mol if r["before"] != r["after"]]
    trans = Counter((r["before"], r["after"]) for r in changed)
    gains = sum(c for (b, a), c in trans.items() if a == "byte_exact")
    losses = sum(c for (b, a), c in trans.items() if b == "byte_exact")
    bad = sum(
        c
        for (b, a), c in trans.items()
        if b == "byte_exact" or a in ("structural", "facmer_divergent", "encode_fail")
    )

    n = len(rows)
    be0 = args.baseline_byte_exact
    be1 = be0 + gains - losses

    print(f"\n=== {RESO}: transition against the v0.4.13 shipped state, n={n} ===")
    print(f"  movers measured : {len(per_mol)}/{len(movers)}", end="")
    print(f"  (unavailable {len(unavailable)}, drift {len(drift)})")
    print(f"  bucket changed  : {len(changed)}")
    for (b, a), c in trans.most_common():
        print(f"    {b:34s} -> {a:24s} {c:5d}")
    if not trans:
        print("    (nothing moved bucket)")

    print(
        f"\nbyte_exact {be0} -> {be1}  ({100 * be0 / n:.2f}% -> {100 * be1 / n:.2f}%, "
        f"{100 * (be1 - be0) / n:+.2f} points)"
    )
    print(
        f"moved in a BAD direction (out of byte_exact, or into structural/facmer/encode_fail): {bad}"
    )

    vc = Counter()
    for r in changed:
        vc[tuple(r["outcome_after"])] += 1
    if vc:
        print("\nveto outcome on the molecules that changed bucket:")
        for k, c in vc.most_common():
            print(f"    {str(k):64s} {c:5d}")

    out = {
        "sweep": args.sweep,
        "lever": RESO,
        "n": n,
        "baseline_byte_exact": be0,
        "byte_exact_after": be1,
        "points": round(100 * (be1 - be0) / n, 4) if n else 0.0,
        "n_movers": len(movers),
        "n_measured": len(per_mol),
        "gains": gains,
        "losses": losses,
        "bad_direction": bad,
        "transitions": {f"{b} -> {a}": c for (b, a), c in trans.items()},
        "changed": changed,
        "excluded_drift": drift,
        "excluded_unavailable": unavailable,
    }
    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nwrote {args.out_json}")


if __name__ == "__main__":
    main()
