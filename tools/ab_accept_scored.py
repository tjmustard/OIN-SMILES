#!/usr/bin/env python
"""A/B one generator lever on runtime AND pass-rate AND structure quality.

Built for `OIN_ACCEPT_SCORED`, whose whole character is a trade: accepting the first
conformer the SCORE would credit is much faster, but it bypasses `_select_by_geometry`'s
clash-first ranking, so the returned structure can be worse. Measuring only runtime would
confirm the change and hide its cost; measuring only pass-rate would miss the point of it.

So every molecule reports three numbers, and a promotion needs all three to hold:

    elapsed_s     wall-clock for encode + generate
    passed        canonical_roundtrip_key(oin_in) == key(get_oin_string(gen.mol, coords))
                  -- the SAME predicate tools/test_dataset_roundtrip.py scores with, so a
                  pass here means a pass there
    clashes       clash.mol_clash_count(gen.mol) on the returned structure

One molecule per subprocess is deliberate: the lever is read from the environment at
predicate-construction time, and a single process cannot host both arms honestly.

Usage:
    python tools/ab_accept_scored.py --cohort cohort.json --out ab.json [--workers 2]
    python tools/ab_accept_scored.py --one <mol.xyz>        # single molecule, one arm
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))


def measure_one(xyz_path: str, timeout: float) -> dict:
    """Encode, generate, score, and count clashes for one molecule in THIS process."""
    import numpy as np

    from oinsmiles import XYZToSMILES
    from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen
    from oinsmiles.generator3d import clash
    from oinsmiles.oin.compare import canonical_roundtrip_key
    from oinsmiles.utils.xyz2mol import get_oin_string

    out: dict = {"molecule": os.path.splitext(os.path.basename(xyz_path))[0]}
    t0 = time.monotonic()
    try:
        oin_in = XYZToSMILES().convert(xyz_path)
        out["oin_in"] = oin_in
        gen = OIN3DGeneratorMetallogen(
            optimizer=None, ensemble_size=1, timeout=timeout, ff_params=None
        )
        res = gen.generate(oin_in)
        out["elapsed_s"] = round(time.monotonic() - t0, 2)

        mol = getattr(res, "mol", None)
        if mol is None:
            out["passed"] = False
            out["note"] = "generator returned no mol (eta fallback)"
            return out

        lines = res.xyz.splitlines()
        natoms = int(lines[0].strip())
        coords = np.array([[float(v) for v in lines[2 + i].split()[1:4]] for i in range(natoms)])
        oin_out = get_oin_string(mol, coords)
        out["oin_out"] = oin_out
        out["passed"] = canonical_roundtrip_key(oin_in) == canonical_roundtrip_key(oin_out)
        out["clashes"] = int(clash.mol_clash_count(mol))
    except Exception as e:
        out.setdefault("elapsed_s", round(time.monotonic() - t0, 2))
        out["passed"] = False
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def run_arm(cohort: list[dict], arm: str, env_extra: dict, timeout: float, workers: int) -> list:
    """Run every molecule of one arm, one subprocess each so the lever is read cleanly."""
    results: list[dict] = []
    pending = list(cohort)
    running: list[tuple] = []
    while pending or running:
        while pending and len(running) < workers:
            rec = pending.pop(0)
            env = dict(os.environ)
            env.update(env_extra)
            env["PYTHONPATH"] = os.path.join(ROOT, "src")
            p = subprocess.Popen(
                [
                    sys.executable,
                    os.path.abspath(__file__),
                    "--one",
                    rec["xyz"],
                    "--timeout",
                    str(timeout),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                env=env,
                text=True,
            )
            running.append((p, rec))
        time.sleep(0.5)
        for p, rec in list(running):
            if p.poll() is None:
                continue
            running.remove((p, rec))
            raw = p.stdout.read() if p.stdout else ""
            try:
                r = json.loads(raw.strip().splitlines()[-1])
            except Exception:
                r = {
                    "molecule": rec["mol"],
                    "passed": False,
                    "error": f"probe produced no parsable result (exit {p.returncode})",
                }
            r["arm"] = arm
            r["eta"] = rec.get("eta")
            r["capstone_elapsed"] = rec.get("elapsed")
            results.append(r)
            print(
                f"  [{arm}] {r['molecule']:16s} "
                f"{r.get('elapsed_s', '?'):>8}s pass={r.get('passed')} "
                f"clash={r.get('clashes', '-')}",
                flush=True,
            )
    return results


def summarize(rows: list[dict], label: str) -> dict:
    done = [r for r in rows if isinstance(r.get("elapsed_s"), (int, float))]
    el = [r["elapsed_s"] for r in done]
    cl = [r["clashes"] for r in rows if isinstance(r.get("clashes"), int)]
    s = {
        "arm": label,
        "n": len(rows),
        "passed": sum(1 for r in rows if r.get("passed")),
        "median_s": round(statistics.median(el), 2) if el else None,
        "mean_s": round(statistics.mean(el), 2) if el else None,
        "total_s": round(sum(el), 1) if el else None,
        "over_30s": sum(1 for v in el if v > 30),
        "clash_total": sum(cl) if cl else 0,
        "clash_mols_with_any": sum(1 for v in cl if v > 0),
    }
    return s


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--one", help="measure a single xyz in this process and print JSON")
    ap.add_argument("--cohort", help="JSON with {'eta': [...], 'control': [...]} records")
    ap.add_argument("--out")
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()

    if args.one:
        print(json.dumps(measure_one(args.one, args.timeout)))
        return 0

    if not args.cohort:
        ap.error("need --cohort or --one")
    raw = json.load(open(args.cohort))
    cohort = list(raw.get("eta", [])) + list(raw.get("control", []))
    cohort = [r for r in cohort if os.path.exists(r["xyz"])]
    print(f"cohort: {len(cohort)} molecules, {args.workers} workers, timeout {args.timeout}s")

    print("\n--- ARM A: default (independent confirm ON) ---")
    a = run_arm(cohort, "A-default", {"OIN_ACCEPT_SCORED": "0"}, args.timeout, args.workers)
    print("\n--- ARM B: OIN_ACCEPT_SCORED=1 ---")
    b = run_arm(cohort, "B-scored", {"OIN_ACCEPT_SCORED": "1"}, args.timeout, args.workers)

    sa, sb = summarize(a, "A-default"), summarize(b, "B-scored")
    by = {r["molecule"]: r for r in a}
    regressions = [
        r["molecule"] for r in b if not r.get("passed") and by.get(r["molecule"], {}).get("passed")
    ]
    fixes = [
        r["molecule"] for r in b if r.get("passed") and not by.get(r["molecule"], {}).get("passed")
    ]

    print("\n================ SUMMARY ================")
    for s in (sa, sb):
        print(
            f"  {s['arm']:10s} pass {s['passed']}/{s['n']}  median {s['median_s']}s  "
            f"total {s['total_s']}s  >30s: {s['over_30s']}  "
            f"clashes {s['clash_total']} over {s['clash_mols_with_any']} mols"
        )
    print(f"  PASS REGRESSIONS (A pass -> B fail): {regressions or 'none'}")
    print(f"  PASS FIXES       (A fail -> B pass): {fixes or 'none'}")

    out = {
        "arm_a": a,
        "arm_b": b,
        "summary_a": sa,
        "summary_b": sb,
        "regressions": regressions,
        "fixes": fixes,
    }
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=1)
        print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
