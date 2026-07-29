#!/usr/bin/env python
"""WHY did the parity veto decline? The five outcomes, counted. (v0.4.14 Lane 1)

WHAT THIS EXISTS TO SEPARATE
============================
v0.4.13 promoted the donor fold and the parity veto together. The fold reaches **393** of the
496 ``key_equal/slot_renumber`` molecules; the veto lets **171** through and reverts **222**.
The release recorded those three numbers and stopped, because the instrument it had --
``fold_transition_sim.py`` -- reports BUCKETS, and every reverted molecule lands in the same
bucket regardless of why.

``fold_parity.resolve`` does not have one reason. It has five, and ``fold_parity.last_outcome``
already exposes them:

    vetoed_collapse            evidence says folding makes the mirror encode identically
    allowed_preexisting_fold   the shipped encoder already folds this pair (achiral, or the
                               held-off metal Delta/Lambda -- not this lever's doing)
    allowed_separation_survives the fold fired and the mirror pair stayed distinct
    declined_no_conformer      \\
    declined_no_pairs           |  NO EVIDENCE EITHER WAY. The veto did not decide the fold
    declined_no_self_encode     |  was unsafe; it never got an instrument reading at all.
    declined_reconstruction_drift/
    fold_inactive              the fold did not fire (not a decline)

The distinction is the whole point. ``vetoed_collapse`` is the veto **working** -- those
molecules encode differently from their round trip because the round trip is a different
enantiomer, and folding them would raise ``byte_exact`` by deleting that fact. A ``declined_*``
molecule is the veto **blind** -- it emits the identical conservative string, so no bucket, no
golden and no mirror audit can tell the two apart. That indistinguishability is not a
hypothetical: ``fold_parity``'s own docstring records a build of this module that declined on
18 of 18 movers while three fixture tests passed.

So: a ``vetoed_collapse`` residue is a GENERATOR finding and must not be "fixed" encoder-side.
A ``declined_*`` residue is recoverable evidence-gathering, and is worth points honestly.

THE MEASUREMENT
===============
Re-encode both strings of every fold-arm mover from COORDINATES with fold+veto on, recording
``last_outcome()`` for each encode. Offline over a frozen sweep -- the generator is not re-run,
which is exact here because the change is encoder-side (``BASELINE.md`` section 3).

    ``smiles_1``       <- the input XYZ from the dataset
    ``smiles_2_indep`` <- the sweep's stored ``structures/*_generated.xyz``

THE CONTROL THAT MUST NOT BE SKIPPED
====================================
Same as ``fold_transition_sim.py``: re-encode fold-OFF first and require it to reproduce the
frozen string. A molecule that fails is reported as ``drift`` and excluded, because its outcome
would be describing an encoder that has moved rather than the veto.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/veto_outcome_audit.py \\
        --sweep <MAIN>/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest \\
        --dataset <MAIN>/tmCAT-tmPHOTO_xyz_dataset/cat \\
        --dataset <MAIN>/tmCAT-tmPHOTO_xyz_dataset/photo \\
        --out-json veto_outcomes.json

⚠ ``--dataset`` has NO default on purpose. ``fold_transition_sim.py`` defaulted to a pair of
RELATIVE paths, which resolve to nothing from a worktree; every molecule landed in
``unavailable`` and the tool printed a refuted number under a correct-looking heading.
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

#: Outcomes that mean "the veto had no reading", as opposed to "the veto decided against".
#: Kept as an explicit set rather than a ``startswith("declined")`` test so that adding an
#: outcome to ``fold_parity`` forces a deliberate choice here instead of silently joining a
#: bucket whose meaning it may not share.
NO_EVIDENCE = frozenset(
    {
        "declined_no_conformer",
        "declined_no_pairs",
        "declined_no_self_encode",
        "declined_reconstruction_drift",
        "declined_no_mirror",
    }
)
DECIDED_AGAINST = frozenset({"vetoed_collapse"})
ALLOWED = frozenset({"allowed_preexisting_fold", "allowed_separation_survives"})


def _set(name: str, on: bool) -> None:
    os.environ[name] = "1" if on else "0"


def _recanon(s):
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


def _find_input_xyz(molecule, dataset_roots):
    for root in dataset_roots:
        p = os.path.join(root, f"{molecule}.xyz")
        if os.path.exists(p):
            return p
    return None


def _encode_with_outcome(path, fold, veto):
    """``(string, outcome)`` for one XYZ under an explicit lever configuration.

    ``outcome`` is ``fold_parity.last_outcome()`` observed immediately after the encode. It is
    reset to ``None`` first so a failed encode cannot report the PREVIOUS molecule's verdict --
    a stale thread-local reading is exactly the kind of plausible nothing this tool exists to
    refuse.
    """
    from oinsmiles import XYZToSMILES

    fold_parity._state.outcome = None
    _set(FOLD, fold)
    _set(VETO, veto)
    try:
        s = XYZToSMILES().convert(path)
    except Exception:  # noqa: BLE001 -- an unencodable structure is excluded, not fatal
        return None, fold_parity.last_outcome()
    return s, fold_parity.last_outcome()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", required=True, help="frozen sweep dir with individual_reports/")
    ap.add_argument(
        "--dataset",
        action="append",
        required=True,
        help="input-xyz root, repeatable. NO DEFAULT -- pass ABSOLUTE paths (see module docstring)",
    )
    ap.add_argument("--limit", type=int, help="first N reports only (smoke test)")
    ap.add_argument("--out-json")
    args = ap.parse_args()

    for root in args.dataset:
        if not os.path.isdir(root):
            sys.exit(f"🔴 --dataset {root} is not a directory")

    rows = load_rows(args.sweep, args.limit)
    print(f"loaded {len(rows)} reports from {args.sweep}", flush=True)

    # Baseline and fold arms are pure string operations -- cheap, and they define the movers.
    _set(FOLD, False)
    _set(VETO, False)
    base_by_mol = {}
    for rep in rows:
        mol = rep.get("molecule")
        s1, s2 = rep.get("smiles_1"), rep.get("smiles_2_indep")
        base_by_mol[mol] = (bucket_of(rep, s1, s2), s1, s2)

    _set(FOLD, True)
    _set(VETO, False)
    movers = []
    for rep in rows:
        mol = rep.get("molecule")
        b0, s1, s2 = base_by_mol[mol]
        if bucket_of(rep, _recanon(s1), _recanon(s2)) != b0:
            movers.append(mol)
    print(f"fold arm moves {len(movers)} molecules -- auditing the veto's verdict on each")

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

        s1_off, _ = _encode_with_outcome(in_xyz, fold=False, veto=False)
        s2_off, _ = _encode_with_outcome(gen_xyz, fold=False, veto=False)
        if s1_off != rep.get("smiles_1") or s2_off != rep.get("smiles_2_indep"):
            drift.append(mol)
            continue

        s1_v, o1 = _encode_with_outcome(in_xyz, fold=True, veto=True)
        s2_v, o2 = _encode_with_outcome(gen_xyz, fold=True, veto=True)
        per_mol.append(
            {
                "molecule": mol,
                "before": base_by_mol[mol][0],
                "after": bucket_of(rep, s1_v, s2_v),
                "outcome_input": o1,
                "outcome_generated": o2,
            }
        )
        if i % 25 == 0:
            print(f"  [{i}/{len(movers)}] {len(per_mol)} measured", flush=True)

    # 🔴 REFUSE rather than report when nothing was measured. An empty audit prints a clean
    # all-zero table that reads exactly like "the veto never declines" -- the most attractive
    # possible wrong answer, and the failure mode that cost v0.4.13 four instruments.
    if movers and not per_mol:
        sys.exit(
            f"\n🔴 REFUSING TO REPORT: measured 0 of {len(movers)} movers "
            f"({len(unavailable)} unavailable, {len(drift)} drift).\n"
            "   Most likely cause: --dataset roots do not contain the cohort's input XYZ files."
        )

    reverted = [r for r in per_mol if r["after"] == r["before"]]
    kept = [r for r in per_mol if r["after"] != r["before"]]

    print(f"\n=== VETO OUTCOMES, n={len(per_mol)} of {len(movers)} movers ===")
    print(f"  kept (fold survived):  {len(kept)}")
    print(f"  reverted by the veto:  {len(reverted)}")
    if unavailable:
        print(f"  ⚠ EXCLUDED unavailable: {len(unavailable)}")
    if drift:
        print(f"  ⚠ EXCLUDED drift:       {len(drift)}")

    def _tally(label, subset):
        print(f"\n--- {label} (n={len(subset)}) ---")
        c = Counter()
        for r in subset:
            c[(r["outcome_input"], r["outcome_generated"])] += 1
        for (a, b), n in c.most_common():
            print(f"  input={a!s:30s} generated={b!s:30s} {n:5d}")
        return {f"{a} | {b}": n for (a, b), n in c.items()}

    kept_t = _tally("KEPT", kept)
    rev_t = _tally("REVERTED", reverted)

    # The headline split: of the reverted molecules, how many did the veto DECIDE against and
    # how many did it merely fail to get a reading on? Only the first group is a real result.
    def _cls(r):
        outs = {r["outcome_input"], r["outcome_generated"]}
        if outs & DECIDED_AGAINST:
            return "decided_against"
        if outs & NO_EVIDENCE:
            return "no_evidence"
        return "other"

    split = Counter(_cls(r) for r in reverted)
    n_rev = len(reverted)
    print(f"\n=== WHY THE {n_rev} REVERTED MOLECULES REVERTED ===")
    for k in ("decided_against", "no_evidence", "other"):
        v = split.get(k, 0)
        pct = f"{100 * v / n_rev:.1f}%" if n_rev else "-"
        print(f"  {k:18s} {v:5d}  ({pct})")
    print(
        "\n  decided_against = the veto working: folding would collapse a mirror pair.\n"
        "                    A GENERATOR finding. Do not 'fix' this encoder-side.\n"
        "  no_evidence     = the veto never got a reading. Recoverable, and worth points."
    )

    out = {
        "sweep": args.sweep,
        "n_reports": len(rows),
        "n_movers": len(movers),
        "n_measured": len(per_mol),
        "kept": len(kept),
        "reverted": n_rev,
        "reverted_split": dict(split),
        "outcomes_kept": kept_t,
        "outcomes_reverted": rev_t,
        "per_molecule": per_mol,
        "excluded_drift": drift,
        "excluded_unavailable": unavailable,
    }
    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nwrote {args.out_json}")


if __name__ == "__main__":
    main()
