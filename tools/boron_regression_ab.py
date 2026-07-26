#!/usr/bin/env python
"""Prove the OIN_BORON_CAGE lever cannot regress a passing molecule.

Two claims, both measured rather than asserted:

1. **Lever OFF is byte-identical to pre-change `main`.** Re-encode a sample of
   currently-passing molecules on this branch with the lever unset and compare
   against the OIN stored in the frozen capstone reports. Any diff is a
   regression in code that was supposed to be inert.

2. **Lever ON does not disturb a non-cage molecule.** Re-encode the *same* sample
   with `OIN_BORON_CAGE=1`. Every gate in the change is behind a B-B-B-triangle
   motif test, so a molecule without a cage must come out byte-identical here too.
   This is the claim worth checking hardest -- it is the one that would let a
   default-ON promotion be safe, and the one an assertion would get wrong.

Also reports, as a control for the round-trip prover: the donor-H delta that
appears when an OIN fragment is re-parsed with its slot markers stripped. If
ordinary passing molecules show the same non-zero delta, that drift is a property
of the OIN format, not something the cage work introduced.

Usage:
    PYTHONPATH=src python tools/boron_regression_ab.py \
        --dataset-dir <abs> --reports-dir <abs>/results-capstone-v042/individual_reports \
        --sample 120 --seed 0
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import random
import subprocess
import sys
from collections import Counter

PER_MOL_TIMEOUT_S = 150


def worker(mol_file: str, dataset_dir: str) -> dict:
    import warnings

    warnings.filterwarnings("ignore")
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from boron_roundtrip import _boron_h, graph_fp, oin_fragments, parse_frag_tolerant

    from oinsmiles import XYZToSMILES
    from oinsmiles.utils.xyz2mol import get_tmc_mol

    path = os.path.join(dataset_dir, mol_file)
    rec = {"mol": mol_file}
    try:
        oin = XYZToSMILES().convert(path)
        rec["oin"] = oin
    except BaseException as e:  # noqa: BLE001
        rec["oin"] = None
        rec["err"] = f"{type(e).__name__}: {str(e)[:120]}"
        return rec
    # donor-H drift control: same measurement the round-trip prover makes
    try:
        tmc, _ = get_tmc_mol(path, 0, with_stereo=False)
        els_enc, _b = graph_fp(tmc)
        els_rt: Counter = Counter()
        for f in oin_fragments(oin):
            m, _mode = parse_frag_tolerant(f)
            if m is not None:
                e, _ = graph_fp(m)
                els_rt += e
        rec["total_H_delta"] = els_rt.get("H", 0) - els_enc.get("H", 0)
        rec["B_H"] = _boron_h(tmc)
    except BaseException:  # noqa: BLE001
        rec["total_H_delta"] = None
    return rec


def run_arm(mols, dataset_dir, lever_on):
    env = dict(os.environ)
    env.pop("OIN_BORON_CAGE", None)
    if lever_on:
        env["OIN_BORON_CAGE"] = "1"
    out = {}
    for i, m in enumerate(mols):
        cmd = [
            sys.executable,
            os.path.abspath(__file__),
            "--worker",
            m,
            "--dataset-dir",
            dataset_dir,
        ]
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=PER_MOL_TIMEOUT_S, env=env
            )
            line = next((ln for ln in proc.stdout.splitlines() if ln.startswith("{")), None)
            out[m] = json.loads(line) if line else {"mol": m, "oin": None, "err": "no output"}
        except subprocess.TimeoutExpired:
            out[m] = {"mol": m, "oin": None, "err": "TIMEOUT"}
        if (i + 1) % 20 == 0:
            print(f"  ...{i + 1}/{len(mols)}", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", default=None)
    ap.add_argument("--dataset-dir", required=True)
    ap.add_argument("--reports-dir", default=None)
    ap.add_argument("--sample", type=int, default=120)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="tools/boron_regression_ab.json")
    args = ap.parse_args()

    if args.worker:
        print(json.dumps(worker(args.worker, args.dataset_dir)), flush=True)
        return

    # collect molecules whose frozen report has a non-null smiles_1 (i.e. the
    # encoder already succeeded on them) -- the population a regression could harm
    stored = {}
    for rp in glob.glob(os.path.join(args.reports_dir, "*.json")):
        try:
            with open(rp) as fh:
                d = json.load(fh)
        except Exception:
            continue
        s1 = d.get("smiles_1")
        name = d.get("molecule") or os.path.basename(rp).replace(".json", "")
        if s1:
            stored[name] = s1
    print(f"{len(stored)} molecules with a stored OIN in {args.reports_dir}")

    names = sorted(stored)
    random.Random(args.seed).shuffle(names)
    picked = []
    for n in names:
        for sub in ("cat", "photo"):
            p = os.path.join(args.dataset_dir, sub, f"{n}.xyz")
            if os.path.exists(p):
                picked.append(os.path.join(sub, f"{n}.xyz"))
                break
        if len(picked) >= args.sample:
            break
    print(f"sampled {len(picked)} molecules (seed {args.seed})")

    print("arm OFF:")
    off = run_arm(picked, args.dataset_dir, lever_on=False)
    print("arm ON:")
    on = run_arm(picked, args.dataset_dir, lever_on=True)

    same_frozen = diff_frozen = off_fail = 0
    off_vs_on_same = off_vs_on_diff = 0
    diffs = []
    h_deltas = Counter()
    for m in picked:
        name = os.path.basename(m).replace(".xyz", "")
        o, n_ = off.get(m, {}), on.get(m, {})
        if o.get("oin") is None:
            off_fail += 1
        else:
            if stored.get(name) == o["oin"]:
                same_frozen += 1
            else:
                diff_frozen += 1
                diffs.append(
                    {
                        "mol": name,
                        "kind": "off_vs_frozen",
                        "frozen": stored.get(name),
                        "off": o["oin"],
                    }
                )
            if o.get("total_H_delta") is not None:
                h_deltas[o["total_H_delta"]] += 1
        if o.get("oin") == n_.get("oin"):
            off_vs_on_same += 1
        else:
            off_vs_on_diff += 1
            diffs.append(
                {"mol": name, "kind": "off_vs_on", "off": o.get("oin"), "on": n_.get("oin")}
            )

    summary = {
        "sampled": len(picked),
        "off_arm_encode_failed": off_fail,
        "off_matches_frozen": same_frozen,
        "off_differs_from_frozen": diff_frozen,
        "off_equals_on": off_vs_on_same,
        "off_differs_from_on": off_vs_on_diff,
        "donor_H_delta_histogram": dict(sorted(h_deltas.items())),
    }
    print("\n=== SUMMARY ===")
    for k, v in summary.items():
        print(f"{k:32s} {v}")
    with open(args.out, "w") as fh:
        json.dump({"summary": summary, "diffs": diffs[:40]}, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
