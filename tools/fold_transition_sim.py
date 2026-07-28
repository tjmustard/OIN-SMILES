#!/usr/bin/env python
"""What does the donor fold move, and what does the parity veto give back? (v0.4.12 Lane 1)

WHY THIS EXISTS AS A TOOL
=========================
v0.4.11 measured the fold's +7.86 points with an ad-hoc script that was never committed, so
the release's headline number could not be re-derived without rewriting it. It is re-derived
here at least three times in v0.4.12 -- once as a sanity anchor before any code is written,
once after the veto lands, and once at close-out -- so it belongs in ``tools/``.

THE MEASUREMENT
===============
Offline simulation over the stored string pairs of a frozen sweep. The v0.4.8 precedent: the
SAME stored strings go through every arm, so nothing generator-side moves underneath and the
A/B confound dissolves. No generation is run and no sweep is required.

    arm base    the frozen classification, recomputed (control: must reproduce the report)
    arm fold    both strings re-canonicalized with OIN_CANONICAL_DONOR_FOLD on
    arm veto    same, plus OIN_FOLD_PARITY_VETO -- requires COORDINATES, see below

⚠ ARM ``veto`` IS NOT A STRING OPERATION, AND ASSUMING IT WAS IS A LIVE TRAP.
``canonicalize_oin_slots`` cannot apply the veto: reflection parity is not a property of the
emitted string, so the veto lives one level up in ``get_oin_string`` where the pristine
conformer is still in hand. Setting ``OIN_FOLD_PARITY_VETO=1`` and re-canonicalizing STRINGS
therefore produces output identical to the ``fold`` arm -- a clean, plausible, and completely
false "the veto costs nothing". The first draft of this tool did exactly that.

So the ``veto`` arm re-encodes FROM COORDINATES: the input XYZ for ``smiles_1`` and the
sweep's stored ``structures/*_generated.xyz`` for ``smiles_2_indep``. It is restricted to the
molecules the ``fold`` arm actually moved, which is exact rather than a sample -- the veto can
only act where the fold fires, so every other molecule is unchanged by construction.

THE CONTROL THAT MUST NOT BE SKIPPED
====================================
The veto arm re-encodes structures the frozen corpus encoded with the v0.4.8 encoder, and
v0.4.9/v0.4.10/v0.4.11 all merely CLAIM to have left the default path byte-identical. So the
arm re-encodes each molecule with the fold OFF first and requires it to reproduce the stored
string. A molecule that fails that control is reported as ``drift`` and excluded, because a
number derived from it would be measuring three releases of encoder change, not the veto.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/fold_transition_sim.py \\
        --sweep tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest \\
        --arm fold --out-json sim.json
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

from oinsmiles.oin.canonical_slots import canonicalize_oin_slots  # noqa: E402

FOLD = "OIN_CANONICAL_DONOR_FOLD"
VETO = "OIN_FOLD_PARITY_VETO"


def _set(name, on):
    os.environ[name] = "1" if on else "0"


def _recanon(s):
    """Re-canonicalize one emitted string under the CURRENT lever settings.

    Returns the input unchanged when the post-pass declines it. That is the conservative
    direction: a string the fold cannot parse simply does not move buckets, exactly as the
    shipped encoder would leave it.
    """
    if not s:
        return s
    try:
        return canonicalize_oin_slots(s)
    except Exception:  # noqa: BLE001 -- an unfoldable string keeps its original classification
        return s


def load_rows(sweep, limit=None):
    rows = []
    for f in sorted(glob.glob(os.path.join(sweep, "individual_reports", "*.json"))):
        try:
            with open(f) as fh:
                rows.append(json.load(fh))
        except Exception:  # noqa: BLE001 -- an unreadable report is counted, not fatal
            continue
        if limit and len(rows) >= limit:
            break
    return rows


def bucket_of(rep, s1, s2):
    """Classify a report with ``smiles_1``/``smiles_2_indep`` overridden.

    Reuses ``roundtrip_bucket_report.classify`` rather than re-implementing it, so the
    ``status`` gate and the ``key_equal`` sub-split cannot drift from the authoritative
    report. Note the honest score is not optional here: the frozen baseline is honest, and a
    scored arm is not comparable to it.
    """
    shadow = dict(rep)
    shadow["smiles_1"] = s1
    shadow["smiles_2_indep"] = s2
    b, sub = classify(shadow, score="honest")
    return f"{b}/{sub}" if sub else b


def _find_input_xyz(molecule, dataset_roots):
    for root in dataset_roots:
        p = os.path.join(root, f"{molecule}.xyz")
        if os.path.exists(p):
            return p
    return None


def _encode_file(path, fold, veto):
    """Encode one XYZ under an explicit lever configuration. ``None`` if it will not encode."""
    from oinsmiles import XYZToSMILES

    _set(FOLD, fold)
    _set(VETO, veto)
    try:
        return XYZToSMILES().convert(path)
    except Exception:  # noqa: BLE001 -- an unencodable structure is excluded, not fatal
        return None


def veto_arm(rows, movers, sweep, dataset_roots):
    """Re-derive the movers' buckets with the veto active, from COORDINATES.

    Returns ``(results, drift, unavailable)``. ``results`` maps molecule -> new bucket for the
    molecules that could be measured; ``drift`` lists molecules whose fold-OFF re-encode did
    not reproduce the stored string (excluded -- see the module docstring); ``unavailable``
    lists molecules whose input or generated structure is missing from disk.
    """
    by_mol = {r.get("molecule"): r for r in rows}
    results, drift, unavailable = {}, [], []
    struct_dir = os.path.join(sweep, "structures")

    for i, mol in enumerate(movers, 1):
        rep = by_mol.get(mol)
        if rep is None:
            unavailable.append(mol)
            continue
        in_xyz = _find_input_xyz(mol, dataset_roots)
        gen_xyz = os.path.join(struct_dir, f"{mol}_generated.xyz")
        if not in_xyz or not os.path.exists(gen_xyz):
            unavailable.append(mol)
            continue

        # CONTROL first: does today's encoder still reproduce the frozen strings?
        s1_off = _encode_file(in_xyz, fold=False, veto=False)
        s2_off = _encode_file(gen_xyz, fold=False, veto=False)
        if s1_off != rep.get("smiles_1") or s2_off != rep.get("smiles_2_indep"):
            drift.append(mol)
            continue

        s1_v = _encode_file(in_xyz, fold=True, veto=True)
        s2_v = _encode_file(gen_xyz, fold=True, veto=True)
        results[mol] = bucket_of(rep, s1_v, s2_v)
        if i % 25 == 0:
            print(f"  [{i}/{len(movers)}] veto arm", flush=True)
    return results, drift, unavailable


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", required=True, help="frozen sweep dir with individual_reports/")
    ap.add_argument("--arm", choices=["fold", "veto"], default="fold")
    ap.add_argument(
        "--dataset",
        action="append",
        default=[],
        help="input-xyz root, repeatable (veto arm only). Defaults to the cat/photo pair.",
    )
    ap.add_argument("--limit", type=int, help="first N reports only (smoke test)")
    ap.add_argument("--out-json")
    args = ap.parse_args()

    rows = load_rows(args.sweep, args.limit)
    print(f"loaded {len(rows)} reports from {args.sweep}", flush=True)

    _set(FOLD, False)
    _set(VETO, False)

    base = Counter()
    moved = Counter()
    per_mol = []
    base_by_mol = {}

    for rep in rows:
        s1 = rep.get("smiles_1")
        s2 = rep.get("smiles_2_indep")
        b0 = bucket_of(rep, s1, s2)
        base[b0] += 1
        base_by_mol[rep.get("molecule")] = (b0, s1, s2)

    _set(FOLD, True)
    _set(VETO, False)

    fold_bucket = {}
    for rep in rows:
        mol = rep.get("molecule")
        b0, s1, s2 = base_by_mol[mol]
        fold_bucket[mol] = bucket_of(rep, _recanon(s1), _recanon(s2))

    drift, unavailable = [], []
    if args.arm == "veto":
        movers = [m for m, b in fold_bucket.items() if b != base_by_mol[m][0]]
        print(f"veto arm: re-encoding {len(movers)} movers from coordinates", flush=True)
        roots = args.dataset or [
            "tmCAT-tmPHOTO_xyz_dataset/cat",
            "tmCAT-tmPHOTO_xyz_dataset/photo",
        ]
        veto_bucket, drift, unavailable = veto_arm(rows, movers, args.sweep, roots)
        # Excluded molecules keep their FOLD-arm bucket rather than silently reverting to
        # baseline: that is the pessimistic direction for the veto's measured cost, so the
        # surviving-gain count below can only be understated, never inflated.
        fold_bucket.update(veto_bucket)
        if drift:
            print(
                f"⚠ EXCLUDED {len(drift)} molecules: fold-OFF re-encode did not reproduce the "
                f"frozen string -- {drift[:8]}{'...' if len(drift) > 8 else ''}"
            )
        if unavailable:
            print(f"⚠ EXCLUDED {len(unavailable)} molecules: structure not on disk")

        # 🔴 REFUSE rather than report when the veto arm measured NOTHING (v0.4.13).
        #
        # Excluded movers keep their FOLD-arm bucket, which is the pessimistic direction for a
        # PARTIAL exclusion -- but it is catastrophic for a TOTAL one: with every mover excluded
        # the veto arm degenerates into the fold arm and prints the bare fold's +7.86, the number
        # v0.4.11 REFUTED, under a heading that says "veto". That is a broken instrument printing
        # a plausible result, which is the exact failure this release exists to be able to catch.
        #
        # It is not hypothetical. Run from a git worktree, `--dataset` defaults to the RELATIVE
        # pair "tmCAT-tmPHOTO_xyz_dataset/{cat,photo}", which does not exist there -- so
        # `_find_input_xyz` returned None for all 393 movers, every one landed in `unavailable`,
        # and the run printed "+7.86 points" and exited 0.
        n_measured = len(veto_bucket)
        if movers and n_measured == 0:
            sys.exit(
                f"\n🔴 REFUSING TO REPORT: the veto arm measured 0 of {len(movers)} movers "
                f"({len(unavailable)} unavailable, {len(drift)} drift).\n"
                "   With every mover excluded this arm IS the fold arm, and would print the bare\n"
                "   fold's +7.86 -- the number v0.4.11 refuted -- labelled as the veto's.\n"
                "   Most likely cause: --dataset defaults are RELATIVE and you are not in the\n"
                "   main checkout. Pass absolute roots:\n"
                "       --dataset <main>/tmCAT-tmPHOTO_xyz_dataset/cat "
                "--dataset <main>/tmCAT-tmPHOTO_xyz_dataset/photo"
            )
        if movers:
            print(f"veto arm measured {n_measured}/{len(movers)} movers")

    for rep in rows:
        mol = rep.get("molecule")
        b0 = base_by_mol[mol][0]
        b1 = fold_bucket[mol]
        moved[b1] += 1
        if b0 != b1:
            per_mol.append({"molecule": mol, "before": b0, "after": b1})

    n = len(rows)
    be0 = base.get("byte_exact", 0)
    be1 = moved.get("byte_exact", 0)

    print(f"\n=== ARM {args.arm} — transition, n={n} ===")
    trans = Counter((r["before"], r["after"]) for r in per_mol)
    for (b, a), c in trans.most_common():
        print(f"  {b:38s} -> {a:24s} {c:5d}")
    if not trans:
        print("  (nothing moved)")

    print(
        f"\nbyte_exact {be0} -> {be1}  ({100 * be0 / n:.2f}% -> {100 * be1 / n:.2f}%, "
        f"{100 * (be1 - be0) / n:+.2f} points)"
    )
    print(
        f"facmer_divergent {base.get('facmer_divergent', 0)} -> {moved.get('facmer_divergent', 0)}"
    )

    # The safety control v0.4.11 relied on, restated as an explicit number rather than a
    # sentence: anything moving OUT of byte_exact, or into structural/facmer, is a defect.
    bad = sum(
        c
        for (b, a), c in trans.items()
        if b == "byte_exact" or a in ("structural", "facmer_divergent", "encode_fail")
    )
    print(
        f"moved in a BAD direction (out of byte_exact, or into structural/facmer/encode_fail): {bad}"
    )

    out = {
        "sweep": args.sweep,
        "arm": args.arm,
        "n": n,
        "base": dict(base),
        "after": dict(moved),
        "byte_exact_before": be0,
        "byte_exact_after": be1,
        "points": round(100 * (be1 - be0) / n, 4) if n else 0.0,
        "bad_direction": bad,
        "transitions": {f"{b} -> {a}": c for (b, a), c in trans.items()},
        "moved": per_mol,
        "excluded_drift": drift,
        "excluded_unavailable": unavailable,
    }
    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nwrote {args.out_json}")


if __name__ == "__main__":
    main()
