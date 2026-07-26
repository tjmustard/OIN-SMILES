#!/usr/bin/env python
"""Whole-encode A/B for ``OIN_VALENCE_CHARGE_FILTER``: what does the OIN string become?

The ligand-level probe (``valorder_probe.py``) shows the filter turning ``found_valid = 0``
into a real Lewis structure. This asks the question that actually matters: **does the emitted
OIN change, and is the new one right?**

"Right" is not "different". Three checks travel with every changed string:

* ``found_valid`` went from 0 to non-zero on the over-cap calls -- so the perception is a
  structure that satisfies ``BO_is_OK``, which *by construction* means every atom's bond
  orders sum to its assigned valence and the fragment's formal charges sum to the target
  charge. That is the valence/charge sanity check, not a proxy for it.
* every body fragment of the emitted OIN re-parses through ``Chem.MolFromSmiles`` and
  survives ``SanitizeMol``, and the parsed heavy-atom count and net formal charge match what
  the string claims. A string that perceives a nicer bond order but no longer re-reads is a
  regression, not a fix.
* both arms run **in one process on the same input file**, and the OFF arm is repeated last
  as a determinism self-check.

    $V tools/valorder_encode_ab.py --dataset <dir> --mols QIDKUL_comp_0 --cap 400
"""

import argparse
import contextlib
import hashlib
import json
import os
import re
import signal
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from rdkit import Chem, RDLogger  # noqa: E402

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.utils import xyz2mol_local  # noqa: E402
from oinsmiles.utils.xyz2mol_local import _CHARGE_FILTER_ENV  # noqa: E402

RDLogger.DisableLog("rdApp.*")

_BODY_SPLIT = re.compile(r"\|")
_SLOT = re.compile(r"\{\d+[<>]?\}")
_METAL_TAG = re.compile(r"\[([A-Z][a-z]?)_[A-Z]{2,3}\]")


class _Timeout(Exception):
    pass


@contextlib.contextmanager
def _time_limit(seconds):
    def _fire(signum, frame):
        raise _Timeout()

    old = signal.signal(signal.SIGALRM, _fire)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


@contextlib.contextmanager
def _silence_fds():
    with open(os.devnull, "w") as devnull:
        old_out, old_err = os.dup(1), os.dup(2)
        try:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            os.dup2(old_out, 1)
            os.dup2(old_err, 2)
            os.close(old_out)
            os.close(old_err)


def _sha(obj):
    return hashlib.sha256(str(obj).encode()).hexdigest()[:12]


def reparse_check(oin):
    """Do the OIN body fragments still re-read as molecules?"""
    if not oin:
        return None
    body = _BODY_SPLIT.split(oin)[0]
    frags = [f for f in body.split(".") if f]
    ok, bad, heavy, charge = 0, [], 0, 0
    for frag in frags:
        smi = _SLOT.sub("", frag)
        smi = _METAL_TAG.sub(r"[\1]", smi)
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            bad.append(smi)
            continue
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:  # noqa: BLE001 - a failed fragment is a result
            bad.append(f"{smi} :: {type(exc).__name__}")
            continue
        ok += 1
        heavy += mol.GetNumHeavyAtoms()
        charge += Chem.GetFormalCharge(mol)
    return {
        "fragments": len(frags),
        "reparsed": ok,
        "failed": bad,
        "heavy_atoms": heavy,
        "net_formal_charge": charge,
    }


def encode(path, filter_on, cap):
    if filter_on:
        os.environ[_CHARGE_FILTER_ENV] = "1"
    else:
        os.environ.pop(_CHARGE_FILTER_ENV, None)
    xyz2mol_local.reset_ac2bo_stats()
    t0 = time.time()
    oin, error = None, None
    try:
        with _time_limit(cap), _silence_fds():
            oin = XYZToSMILES().convert(str(path))
    except _Timeout:
        error = f"TIMEOUT after {cap}s"
    except Exception as exc:  # noqa: BLE001 - a failed arm is a result
        error = f"{type(exc).__name__}: {exc}"
    stats = dict(xyz2mol_local.AC2BO_STATS)
    return {
        "filter": "ON" if filter_on else "OFF",
        "oin": oin,
        "oin_sha": _sha(oin) if oin else None,
        "error": error,
        "wall": round(time.time() - t0, 2),
        "stats": stats,
        "reparse": reparse_check(oin),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=".")
    ap.add_argument("--mols", default="", help="refcodes resolved under --dataset")
    ap.add_argument("--files", default="", help="explicit xyz paths (the goldens live in tests/)")
    ap.add_argument("--cap", type=float, default=400.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = Path(args.dataset)
    targets = []
    for name in [x.strip() for x in args.mols.split(",") if x.strip()]:
        hits = sorted(root.glob(f"cat/{name}.xyz")) + sorted(root.glob(f"photo/{name}.xyz"))
        if not hits:
            print(f"{name}: NOT FOUND", flush=True)
            continue
        targets.append((name, hits[0]))
    for raw in [x.strip() for x in args.files.split(",") if x.strip()]:
        path = Path(raw)
        if not path.exists():
            print(f"{raw}: NOT FOUND", flush=True)
            continue
        targets.append((path.stem, path))
    if not targets:
        sys.exit("FATAL: nothing to encode -- pass --mols and/or --files")
    out = {}
    for name, path in targets:
        print(f"\n=== {name} ({path.parent.name})", flush=True)
        arms = [encode(path, False, args.cap), encode(path, True, args.cap)]
        arms.append(encode(path, False, args.cap))  # repeat OFF: determinism self-check
        for i, arm in enumerate(arms):
            st = arm["stats"]
            tag = " (repeat)" if i == 2 else ""
            print(
                f"  filter={arm['filter']:<3s}{tag} wall={arm['wall']:>7.2f}s "
                f"over_cap={st['over_cap_calls']} filtered={st['over_cap_filtered_calls']} "
                f"cands={st['candidates']} found_valid={st['found_valid']} "
                f"exhausted={st['over_cap_exhausted']} oin_sha={arm['oin_sha']}"
                + (f"  {arm['error']}" if arm["error"] else ""),
                flush=True,
            )
            if arm["reparse"]:
                rp = arm["reparse"]
                print(
                    f"        reparse {rp['reparsed']}/{rp['fragments']} frags, "
                    f"heavy={rp['heavy_atoms']} net_charge={rp['net_formal_charge']}"
                    + (f"  FAILED {rp['failed']}" if rp["failed"] else ""),
                    flush=True,
                )
        repeat_ok = arms[0]["oin_sha"] == arms[2]["oin_sha"]
        print(f"  self-check: {'REPEAT-OK' if repeat_ok else 'REPEAT-MISMATCH'}", flush=True)
        changed = arms[0]["oin_sha"] != arms[1]["oin_sha"]
        print(f"  OIN {'CHANGED' if changed else 'SAME'}", flush=True)
        if arms[0]["oin"]:
            print(f"    OFF: {arms[0]['oin']}", flush=True)
        if arms[1]["oin"]:
            print(f"    ON : {arms[1]['oin']}", flush=True)
        out[name] = {"path": str(path), "repeat_ok": repeat_ok, "changed": changed, "arms": arms}
        if args.out:
            Path(args.out).write_text(json.dumps(out, indent=1))
    if args.out:
        print(f"\nwrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
