#!/usr/bin/env python
"""Q2: does the fallback answer depend on how long you search?

For each molecule, encode it once per ``OIN_VALENCE_FALLBACK_TRIES`` budget **in one
process** (the budget is read per ``AC2BO`` call, so no subprocess is needed) and record:

* the sha of every ``AC2BO`` return value, tagged over-cap / sub-cap and
  valid-Lewis-structure / fallback-guess -- this is the ``best_BO`` half of the question;
* the **emitted OIN string** -- the answer that actually matters, since a changed
  intermediate bond order may be washed out by later normalization;
* the deterministic counters from ``xyz2mol_local.AC2BO_STATS``.

Two properties make this a fair A/B:

* Arms run **in one process against the same input file**, so nothing but the budget
  differs. Wall clock is still reported but is not the claim -- the host runs the release
  sweep and load is above 12, where wall clock is meaningless (RESTART.md §6).
* The budget list is run in the given order and the **20 000 arm is repeated last** as a
  self-check: if the first and last 20 000 arms disagree, the harness itself is
  nondeterministic and every other row is void. This is printed as ``REPEAT-OK`` /
  ``REPEAT-MISMATCH`` rather than assumed.

    $V tools/valsearch_budget_ab.py --mols QIDKUL_comp_0 --dataset <dir> \
        --budgets 20000,5000,1000,200
"""

import argparse
import contextlib
import hashlib
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

import re  # noqa: E402

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.utils import xyz2mol_local  # noqa: E402

_CHARGE_RE = re.compile(r"Charge:\s*(-?\d+)")


@contextlib.contextmanager
def _silence_fds():
    """Redirect C-level stdout/stderr to devnull (openbabel prints distance warnings)."""
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


class Recorder:
    """Wrap ``AC2BO`` to record what each call returned, and whether it was a guess."""

    def __init__(self):
        self.original = xyz2mol_local.AC2BO
        self.calls = []

    def __enter__(self):
        def wrapped(AC, atoms, charge, **kwargs):
            before = dict(xyz2mol_local.AC2BO_STATS)
            t0 = time.time()
            BO, ave = self.original(AC, atoms, charge, **kwargs)
            after = xyz2mol_local.AC2BO_STATS
            self.calls.append(
                {
                    "ac_sha": _sha(AC.tolist()),
                    "charge": charge,
                    "allow_carbenes": kwargs.get("allow_carbenes", True),
                    "over_cap": after["over_cap_calls"] > before["over_cap_calls"],
                    # A call that early-returned found a real Lewis structure; one that did
                    # not returned best_BO, i.e. a guess.
                    "found_valid": after["found_valid"] > before["found_valid"],
                    "candidates": after["candidates"] - before["candidates"],
                    "bo_sha": _sha(BO.tolist()),
                    "bo_sum": int(BO.sum()),
                    "wall": round(time.time() - t0, 3),
                }
            )
            return BO, ave

        xyz2mol_local.AC2BO = wrapped
        return self

    def __exit__(self, *exc):
        xyz2mol_local.AC2BO = self.original
        return False


def encode_at_budget(path, budget):
    """One encode with ``OIN_VALENCE_FALLBACK_TRIES=budget``. Returns a result dict."""
    os.environ["OIN_VALENCE_FALLBACK_TRIES"] = str(budget)
    xyz2mol_local.reset_ac2bo_stats()
    t0 = time.time()
    with Recorder() as rec:
        try:
            with _silence_fds():
                oin = XYZToSMILES().convert(str(path))
            error = None
        except Exception as exc:  # noqa: BLE001 - a failed arm is a result, not a crash
            oin, error = None, f"{type(exc).__name__}: {exc}"
    wall = time.time() - t0
    stats = dict(xyz2mol_local.AC2BO_STATS)
    guesses = [c for c in rec.calls if c["over_cap"] and not c["found_valid"]]
    return {
        "budget": budget,
        "oin": oin,
        "oin_sha": _sha(oin) if oin else None,
        "error": error,
        "wall": round(wall, 2),
        "stats": stats,
        "calls": rec.calls,
        # The fingerprint of the perception decisions: every over-cap fallback guess, in
        # call order. This is what a budget cut is expected to move.
        "guess_sha": _sha([c["bo_sha"] for c in guesses]),
        "n_guesses": len(guesses),
        "all_bo_sha": _sha([c["bo_sha"] for c in rec.calls]),
    }


def read_charge(path):
    with open(path) as fh:
        fh.readline()
        comment = fh.readline()
    match = _CHARGE_RE.search(comment)
    return int(match.group(1)) if match else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--mols", required=True, help="comma-separated refcodes, or @file")
    ap.add_argument("--budgets", default="20000,5000,1000,200")
    ap.add_argument("--out", default=None)
    ap.add_argument(
        "--one-arm",
        type=int,
        default=None,
        help="Encode each molecule at exactly this budget and print one line per molecule, "
        "then exit. For ligands where the 20 000 arm never returns, an in-process sweep "
        "cannot be bounded -- wrap this in `timeout N` instead, so a non-terminating arm "
        "is recorded as a timeout rather than hanging the whole comparison.",
    )
    args = ap.parse_args()

    root = Path(args.dataset)
    if args.mols.startswith("@"):
        names = [x.strip() for x in Path(args.mols[1:]).read_text().split() if x.strip()]
    else:
        names = [x.strip() for x in args.mols.split(",") if x.strip()]
    budgets = [int(x) for x in args.budgets.split(",")]
    # Repeat the first budget last: a same-process determinism self-check.
    order = budgets + [budgets[0]]

    if args.one_arm is not None:
        for name in names:
            hits = list(root.glob(f"*/{name}.xyz"))
            if not hits:
                print(f"{name}: NOT FOUND under {root}", flush=True)
                continue
            res = encode_at_budget(hits[0], args.one_arm)
            print(
                f"{name} tries={args.one_arm} wall={res['wall']:.2f}s "
                f"cands={res['stats']['candidates']} "
                f"over_cap_calls={res['stats']['over_cap_calls']} "
                f"exhausted={res['stats']['over_cap_exhausted']} "
                f"found_valid={res['stats']['found_valid']} "
                f"matching={res['stats']['matching_calls']} "
                f"guess_sha={res['guess_sha']} oin_sha={res['oin_sha']}"
                + (f" ERROR {res['error']}" if res["error"] else ""),
                flush=True,
            )
            if args.out:
                Path(args.out).write_text(json.dumps({name: res}, indent=1))
        return

    out = {}
    for name in names:
        hits = list(root.glob(f"*/{name}.xyz"))
        if not hits:
            print(f"{name}: NOT FOUND under {root}", flush=True)
            continue
        path = hits[0]
        print(f"\n=== {name}  ({path})", flush=True)
        arms = []
        for i, budget in enumerate(order):
            res = encode_at_budget(path, budget)
            arms.append(res)
            tag = "  (repeat)" if i == len(order) - 1 else ""
            print(
                f"  tries={budget:<6d} wall={res['wall']:>7.2f}s "
                f"cands={res['stats']['candidates']:<7d} "
                f"guesses={res['n_guesses']} guess_sha={res['guess_sha']} "
                f"oin_sha={res['oin_sha']}{tag}"
                + (f"  ERROR {res['error']}" if res["error"] else ""),
                flush=True,
            )
        repeat_ok = arms[0]["oin_sha"] == arms[-1]["oin_sha"] and (
            arms[0]["guess_sha"] == arms[-1]["guess_sha"]
        )
        print(f"  self-check: {'REPEAT-OK' if repeat_ok else 'REPEAT-MISMATCH'}", flush=True)
        base = arms[0]
        for res in arms[1:-1]:
            same_oin = res["oin_sha"] == base["oin_sha"]
            same_guess = res["guess_sha"] == base["guess_sha"]
            print(
                f"    {res['budget']:<6d} vs 20000: "
                f"OIN {'SAME' if same_oin else 'CHANGED'} · "
                f"best_BO {'SAME' if same_guess else 'CHANGED'}",
                flush=True,
            )
        out[name] = {"path": str(path), "repeat_ok": repeat_ok, "arms": arms}

    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
