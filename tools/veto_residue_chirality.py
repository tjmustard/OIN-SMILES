#!/usr/bin/env python
"""Are the veto-reverted molecules a GENERATOR chirality error, or an encoder gap? (v0.4.14)

THE CLAIM THIS EXISTS TO STOP ME MAKING WITHOUT EVIDENCE
========================================================
``tools/veto_outcome_audit.py`` measures that all 222 molecules the parity veto reverts are
``vetoed_collapse`` rather than ``declined_*``. It is tempting to read that as "the round trip
produced the wrong enantiomer, so those 4.44 points are a generator problem". **It does not say
that.**

``vetoed_collapse`` says only: *this structure's mirror encodes differently today, and folding
would make it encode identically.* It is a property of one structure. It says nothing about the
relationship between the INPUT and the ROUND TRIP, which is the thing the release would be
re-filing on the roadmap.

The two possibilities have completely different owners:

    (a) GENERATOR CHIRALITY ERROR -- the generated structure really is the input's mirror image.
        ``byte_exact`` failing is then CORRECT, the veto is protecting a true negative, and no
        encoder change should ever recover these. Owner: the generator.

    (b) ENCODER GAP -- the two structures are the same enantiomer, labelled two ways that this
        particular fold cannot unify without collateral damage to the mirror distinction. Then a
        *parity-aware* canonicalization could still reach them. Owner: a future encoder lane.

THE TEST
========
Both strings are re-derived from coordinates with the fold OFF (the rotation-only labeling), and
the input is additionally encoded MIRRORED:

    s1    = encode(input.xyz)
    s2    = encode(generated.xyz)          (the sweep's stored structure)
    s1_m  = encode(mirror(input.xyz))

    s1_m == s2   ->  MIRROR_MATCH   : the round trip built the enantiomer            (a)
    s1_m != s2   ->  NOT_A_MIRROR   : same enantiomer, labelled differently          (b)

``s1 == s2`` would mean the pair is not divergent at all and the molecule does not belong to
this population; it is counted as ``UNEXPECTED_IDENTICAL`` and reported rather than dropped,
because silently discarding it would flatter whichever verdict is left.

THE CONTROL
===========
The fold-OFF re-encode must reproduce the frozen ``smiles_1`` / ``smiles_2_indep``. A molecule
that fails is ``drift`` and is excluded -- otherwise the verdict describes encoder movement
since v0.4.8 rather than the round trip.

The mirror is a z-negation, matching ``fold_parity._mirror_coords`` and
``tools/mirror_audit_donor_fold.py`` line for line, so all three agree by construction rather
than by three independently-reasoned reflections.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/veto_residue_chirality.py \\
        --outcomes <results>/veto_outcomes.json \\
        --sweep <MAIN>/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest \\
        --dataset <MAIN>/tmCAT-tmPHOTO_xyz_dataset/cat \\
        --dataset <MAIN>/tmCAT-tmPHOTO_xyz_dataset/photo --out-json residue_chirality.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from collections import Counter

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

FOLD = "OIN_CANONICAL_DONOR_FOLD"
VETO = "OIN_FOLD_PARITY_VETO"
RESO = "OIN_RESONANCE_DONOR_FOLD"


def _rotation_only():
    """Rotation-only labeling: every donor-fold widening off, veto irrelevant."""
    os.environ[FOLD] = "0"
    os.environ[VETO] = "0"
    os.environ[RESO] = "0"


def read_xyz(path):
    with open(path) as fh:
        lines = fh.readlines()
    n = int(lines[0].split()[0])
    syms, coords = [], []
    for line in lines[2 : 2 + n]:
        p = line.split()
        syms.append(p[0])
        coords.append([float(x) for x in p[1:4]])
    return syms, coords


def _encode(path):
    from oinsmiles import XYZToSMILES

    try:
        return XYZToSMILES().convert(path)
    except Exception:  # noqa: BLE001 -- an unencodable structure is excluded, not fatal
        return None


def _encode_mirror(path):
    syms, coords = read_xyz(path)
    fh = tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False)
    try:
        fh.write(f"{len(syms)}\n\n")
        for sym, (x, y, z) in zip(syms, coords):
            fh.write(f"{sym:<3} {x:>16.10f} {y:>16.10f} {-z:>16.10f}\n")
        fh.close()
        return _encode(fh.name)
    finally:
        try:
            os.unlink(fh.name)
        except OSError:
            pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outcomes", required=True, help="veto_outcome_audit.py --out-json")
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--dataset", action="append", required=True, help="ABSOLUTE input-xyz root")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--out-json")
    args = ap.parse_args()

    for root in args.dataset:
        if not os.path.isdir(root):
            sys.exit(f"🔴 --dataset {root} is not a directory")

    audit = json.load(open(args.outcomes))
    reverted = [
        r["molecule"]
        for r in audit["per_molecule"]
        if r["after"] == r["before"]
        and "vetoed_collapse" in (r["outcome_input"], r["outcome_generated"])
    ]
    if args.limit:
        reverted = reverted[: args.limit]
    print(f"{len(reverted)} veto-reverted molecules to classify", flush=True)
    if not reverted:
        sys.exit("🔴 REFUSING: the outcomes file lists no vetoed_collapse molecules.")

    struct_dir = os.path.join(args.sweep, "structures")
    frozen = {}
    for mol in reverted:
        p = os.path.join(args.sweep, "individual_reports", f"{mol}.json")
        if os.path.exists(p):
            frozen[mol] = json.load(open(p))

    _rotation_only()
    rows, drift, unavailable = [], [], []
    for i, mol in enumerate(reverted, 1):
        rep = frozen.get(mol)
        in_xyz = next(
            (
                os.path.join(r, f"{mol}.xyz")
                for r in args.dataset
                if os.path.exists(os.path.join(r, f"{mol}.xyz"))
            ),
            None,
        )
        gen_xyz = os.path.join(struct_dir, f"{mol}_generated.xyz")
        if rep is None or in_xyz is None or not os.path.exists(gen_xyz):
            unavailable.append(mol)
            continue

        s1, s2 = _encode(in_xyz), _encode(gen_xyz)
        if s1 != rep.get("smiles_1") or s2 != rep.get("smiles_2_indep"):
            drift.append(mol)
            continue

        s1m = _encode_mirror(in_xyz)
        if s1m is None:
            unavailable.append(mol)
            continue

        if s1 == s2:
            verdict = "UNEXPECTED_IDENTICAL"
        elif s1m == s2:
            verdict = "MIRROR_MATCH"
        else:
            verdict = "NOT_A_MIRROR"
        rows.append({"molecule": mol, "verdict": verdict})
        if i % 25 == 0:
            print(f"  [{i}/{len(reverted)}] {len(rows)} classified", flush=True)

    # 🔴 REFUSE rather than report an all-zero table that reads like a clean verdict.
    if reverted and not rows:
        sys.exit(
            f"\n🔴 REFUSING TO REPORT: classified 0 of {len(reverted)} "
            f"({len(unavailable)} unavailable, {len(drift)} drift)."
        )

    tally = Counter(r["verdict"] for r in rows)
    n = len(rows)
    print(f"\n=== VETO-REVERTED RESIDUE, n={n} of {len(reverted)} ===")
    print(f"  excluded: unavailable {len(unavailable)}, drift {len(drift)}")
    for k, c in tally.most_common():
        print(f"  {k:22s} {c:5d}  ({100 * c / n:.1f}%)")
    print(
        "\n  MIRROR_MATCH  -> the round trip built the ENANTIOMER. byte_exact failing is CORRECT,\n"
        "                   the veto is protecting a true negative, and no encoder change should\n"
        "                   recover these. A GENERATOR lane owns them.\n"
        "  NOT_A_MIRROR  -> same enantiomer, two labelings this fold cannot unify. Still an\n"
        "                   ENCODER gap, reachable by a parity-aware canonicalization."
    )

    out = {
        "outcomes": args.outcomes,
        "sweep": args.sweep,
        "n_reverted": len(reverted),
        "n_classified": n,
        "tally": dict(tally),
        "rows": rows,
        "excluded_drift": drift,
        "excluded_unavailable": unavailable,
    }
    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nwrote {args.out_json}")


if __name__ == "__main__":
    main()
