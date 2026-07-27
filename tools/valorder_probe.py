#!/usr/bin/env python
"""Is ``found_valid = 0`` on the over-cap branch a property of the ligand, or of the ORDER?

``docs/agentic-notes/v0.4.5/VALENCE_SEARCH_v0.4.5.md`` measured that ``AC2BO``'s over-cap branch grinds 20 000
candidates and finds none valid, so it returns ``best_BO`` -- a guess. But that branch also
skips ``_ordered_valences``, the O/N/C/P/S grouping heuristic, and sub-cap ligands find a
valid Lewis structure in **3 to 8** candidates *because* of it. So the failure may be a
statement about search order rather than about the ligand.

This probe answers that **on the ligand, not through a whole encode**, which is what makes it
affordable: the largest over-cap ligands cost ~0.36 s per candidate, so a whole-encode A/B is
a timeout in both arms and measures nothing.

Two phases:

* **capture** -- run the real encoder with ``AC2BO`` wrapped, and on the first call whose
  ``valence_combo_size`` exceeds the cap, record ``(AC, atoms, charge, allow_carbenes)`` and
  abort the encode. These are the encoder's own inputs, not a reconstruction, so the probe
  cannot disagree with production about what it is searching. Cached to ``--cache`` so the
  arms below never pay for it twice.
* **probe** -- call ``AC2BO`` directly, once per arm, with a wall-clock cap per arm
  (``signal.setitimer``; the search is pure Python so the handler fires promptly). A capped
  arm still reports its counters, because ``AC2BO_STATS`` is incremented as it goes -- so
  "timed out after N candidates, none valid" is a result rather than a lost run.

The number the hypothesis lives or dies by is ``candidates`` on a ``found_valid=1`` arm: the
candidate index at which a valid Lewis structure appears.

    $V tools/valorder_probe.py --dataset <dir> --mols HICLAG_comp_0 \
        --tries 20000 --cap 120 --arms raw,ordered
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

import numpy as np  # noqa: E402

from oinsmiles.utils import perception_core  # noqa: E402
from oinsmiles.utils.perception_core import (  # noqa: E402
    _CHARGE_FILTER_ENV,
    _ORDERED_FALLBACK_ENV,
    _VALENCE_COMBO_CAP,
    possible_valences,
    valence_combo_size,
)

_CHARGE_RE = re.compile(r"Charge:\s*(-?\d+)")


class _Captured(Exception):
    """Sentinel: the over-cap AC2BO inputs are in hand, stop the encode."""


class _ArmTimeout(Exception):
    pass


@contextlib.contextmanager
def _time_limit(seconds):
    def _fire(signum, frame):
        raise _ArmTimeout()

    old = signal.signal(signal.SIGALRM, _fire)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


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


def capture_over_cap_inputs(path, cache=None):
    """The encoder's own ``(AC, atoms, charge, allow_carbenes)`` for the first over-cap call."""
    if cache and Path(cache).exists():
        data = np.load(cache, allow_pickle=False)
        return {
            "AC": data["AC"],
            "atoms": list(data["atoms"]),
            "charge": int(data["charge"]),
            "allow_carbenes": bool(data["allow_carbenes"]),
            "cached": True,
        }

    from oinsmiles import XYZToSMILES

    grabbed = {}
    original = perception_core.AC2BO

    def wrapped(AC, atoms, charge, **kwargs):
        allow_carbenes = kwargs.get("allow_carbenes", True)
        vll = possible_valences(list(AC.sum(axis=1)), atoms, allow_carbenes=allow_carbenes)
        if valence_combo_size(vll) > _VALENCE_COMBO_CAP:
            grabbed.update(
                AC=AC.copy(),
                atoms=list(atoms),
                charge=int(charge),
                allow_carbenes=bool(allow_carbenes),
                cached=False,
            )
            raise _Captured()
        return original(AC, atoms, charge, **kwargs)

    perception_core.AC2BO = wrapped
    try:
        with contextlib.suppress(Exception), _silence_fds():
            XYZToSMILES().convert(str(path))
    finally:
        perception_core.AC2BO = original

    if not grabbed:
        return None
    if cache:
        np.savez(
            cache,
            AC=grabbed["AC"],
            atoms=np.array(grabbed["atoms"]),
            charge=grabbed["charge"],
            allow_carbenes=grabbed["allow_carbenes"],
        )
    return grabbed


def run_arm(inputs, arm, tries, cap):
    """One ``AC2BO`` call under one enumeration strategy and one budget."""
    for env in (_ORDERED_FALLBACK_ENV, _CHARGE_FILTER_ENV):
        os.environ.pop(env, None)
    if arm == "ordered":
        os.environ[_ORDERED_FALLBACK_ENV] = "1"
    elif arm == "filtered":
        os.environ[_CHARGE_FILTER_ENV] = "1"
    os.environ["OIN_VALENCE_FALLBACK_TRIES"] = str(tries)

    perception_core.reset_ac2bo_stats()
    t0 = time.time()
    timed_out = False
    BO = None
    try:
        with _time_limit(cap):
            BO, _ = perception_core.AC2BO(
                inputs["AC"],
                inputs["atoms"],
                inputs["charge"],
                allow_charged_fragments=True,
                use_graph=True,
                allow_carbenes=inputs["allow_carbenes"],
            )
    except _ArmTimeout:
        timed_out = True
    wall = time.time() - t0
    stats = dict(perception_core.AC2BO_STATS)
    cands = stats["candidates"]
    return {
        "arm": arm,
        "tries": tries,
        "cap_s": cap,
        "timed_out": timed_out,
        "wall": round(wall, 2),
        "candidates": cands,
        "s_per_candidate": round(wall / cands, 4) if cands else None,
        "found_valid": stats["found_valid"],
        "over_cap_calls": stats["over_cap_calls"],
        "over_cap_ordered_calls": stats["over_cap_ordered_calls"],
        "over_cap_filtered_calls": stats["over_cap_filtered_calls"],
        "exhausted": stats["over_cap_exhausted"],
        "best_bo_improved": stats["over_cap_best_bo_improved"],
        "matching_calls": stats["matching_calls"],
        "bo_sha": _sha(BO.tolist()) if BO is not None else None,
        "bo_sum": int(BO.sum()) if BO is not None else None,
        "ac_sum": int(inputs["AC"].sum()),
        "bo_is_ac": bool(BO is not None and np.array_equal(BO, inputs["AC"])),
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
    ap.add_argument("--tries", default="20000")
    ap.add_argument("--arms", default="raw,ordered")
    ap.add_argument("--cap", type=float, default=120.0, help="wall cap per arm, seconds")
    ap.add_argument("--cache-dir", default=None, help="cache captured AC2BO inputs here")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = Path(args.dataset)
    if args.mols.startswith("@"):
        names = [x.strip() for x in Path(args.mols[1:]).read_text().split() if x.strip()]
    else:
        names = [x.strip() for x in args.mols.split(",") if x.strip()]
    budgets = [int(x) for x in args.tries.split(",")]
    arms = [a.strip() for a in args.arms.split(",") if a.strip()]

    out = {}
    for name in names:
        hits = sorted(root.glob(f"cat/{name}.xyz")) + sorted(root.glob(f"photo/{name}.xyz"))
        if not hits:
            print(f"{name}: NOT FOUND under {root}/cat or {root}/photo", flush=True)
            continue
        path = hits[0]
        cache = None
        if args.cache_dir:
            Path(args.cache_dir).mkdir(parents=True, exist_ok=True)
            cache = str(Path(args.cache_dir) / f"{name}.npz")
        t0 = time.time()
        inputs = capture_over_cap_inputs(path, cache=cache)
        if inputs is None:
            print(f"{name}: no over-cap AC2BO call reached (not an over-cap molecule?)", flush=True)
            continue
        vll = possible_valences(
            list(inputs["AC"].sum(axis=1)), inputs["atoms"], allow_carbenes=inputs["allow_carbenes"]
        )
        print(
            f"\n=== {name}  ({path.parent.name})  ligand {len(inputs['atoms'])} atoms  "
            f"charge={inputs['charge']}  carbenes={inputs['allow_carbenes']}  "
            f"combo_size>{_VALENCE_COMBO_CAP} ({valence_combo_size(vll, cap=10**18)})  "
            f"capture {time.time() - t0:.1f}s{' (cached)' if inputs['cached'] else ''}",
            flush=True,
        )
        rows = []
        for tries in budgets:
            for arm in arms:
                res = run_arm(inputs, arm, tries=tries, cap=args.cap)
                rows.append(res)
                verdict = (
                    "TIMEOUT"
                    if res["timed_out"]
                    else ("VALID" if res["found_valid"] else "no-valid")
                )
                print(
                    f"  {arm:<8s} tries={tries:<6d} {verdict:<9s} "
                    f"cands={res['candidates']:<6d} "
                    f"s/cand={res['s_per_candidate']} "
                    f"wall={res['wall']:>7.2f}s "
                    f"bo_sum={res['bo_sum']} ac_sum={res['ac_sum']} "
                    f"bo_sha={res['bo_sha']}",
                    flush=True,
                )
        out[name] = {
            "path": str(path),
            "ligand_atoms": len(inputs["atoms"]),
            "charge": inputs["charge"],
            "allow_carbenes": inputs["allow_carbenes"],
            "combo_size": valence_combo_size(vll, cap=10**18),
            "arms": rows,
        }
        if args.out:
            Path(args.out).write_text(json.dumps(out, indent=1))

    if args.out:
        print(f"\nwrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
