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
    clash_vdw     non-bonded pairs inside vdW contact, computed FROM THE RETURNED COORDINATES
    clash_severe  the subset below the severe cutoff
    worst_overlap smallest non-bonded dist/(rvdw+rvdw); continuous, so the comparison is not
                  hostage to where clash_cutoff happens to sit
    sha_in        sha256(smiles_1). CONTROL: the lever is generation-side, so this must not
                  move between arms. If it does, the run is confounded and nothing else in it
                  can be read.
    sha_out       sha256(smiles_2), where smiles_2 is `get_oin_string(gen.mol, coords)` --
                  byte-for-byte the string tools/test_dataset_roundtrip.py records.
    indep_*       the test the lever DROPS: a full XYZToSMILES().convert() of the generated
                  XYZ, compared to oin_in. `passed` above cannot see this by construction --
                  it scores with the very predicate the lever accepts on, so 18/22 -> 18/22
                  is not evidence of no loss. This is.

WHY sha_out AND NOT JUST `passed`
---------------------------------
`passed` compares `canonical_roundtrip_key`, which FOLDS slot renumbering, fragment reorder
and metal-@. The v0.4.7 wave's contract is a BYTE-identical OIN string. A conformer can
satisfy the key and still emit different bytes, so a key-only A/B cannot clear a byte gate.
Pure string equality is also load-independent, which the seconds in this script are not.

⚠ Do NOT use `clash.mol_clash_count(gen_result.mol)` here. It duck-types on `mol.atom_list`
and RETURNS 0 on AttributeError, and `gen_result.mol` is a bare `rdkit.Chem.rdchem.Mol` with
no such attribute. The first version of this script did exactly that and reported "clash 0"
for all 44 measurements in both arms -- a degenerate metric that would have certified the
quality arm without measuring it. `vdw_clash_count(positions, atomic_numbers)` takes raw
coordinates and is the honest call at this site.

One molecule per subprocess is deliberate: the lever is read from the environment at
predicate-construction time, and a single process cannot host both arms honestly.

Usage:
    python tools/ab_accept_scored.py --cohort cohort.json --out ab.json [--workers 2]
    python tools/ab_accept_scored.py --one <mol.xyz>        # single molecule, one arm
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))


def _sha(s) -> str | None:
    if not isinstance(s, str):
        return None
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:16]


def measure_one(xyz_path: str, timeout: float, dump_dir: str | None = None) -> dict:
    """Encode, generate, score, hash, re-perceive, and count clashes for one molecule here.

    ``dump_dir`` (L5-attach): additionally persist the ACCEPTED conformer's coordinates and
    the generator's own claimed metal-donor set. Nothing in this script's verdicts depends on
    it. It exists because the run JSONs record only derived scalars (sha/clash/oin), so the
    §6.5 falsification -- "does a coordinate-only donor-set predicate separate arm A's
    accepted conformer from arm B's?" -- was otherwise unanswerable without regenerating
    every molecule once per candidate predicate. With the dump, generation is paid once and
    predicate iteration is free and offline.
    """
    import tempfile

    import numpy as np
    from rdkit.Chem import GetPeriodicTable

    from oinsmiles import XYZToSMILES
    from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen
    from oinsmiles.generator3d import clash
    from oinsmiles.oin.compare import canonical_roundtrip_key
    from oinsmiles.utils.perception_tmc import get_oin_string

    out: dict = {"molecule": os.path.splitext(os.path.basename(xyz_path))[0]}
    t0 = time.monotonic()
    try:
        oin_in = XYZToSMILES().convert(xyz_path)
        out["oin_in"] = oin_in
        out["sha_in"] = _sha(oin_in)
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
        out["sha_out"] = _sha(oin_out)
        out["passed"] = canonical_roundtrip_key(oin_in) == canonical_roundtrip_key(oin_out)
        out["byte_equal_in_out"] = oin_out == oin_in

        pt = GetPeriodicTable()
        znums = [pt.GetAtomicNumber(lines[2 + i].split()[0]) for i in range(natoms)]
        cv, cs, worst = clash.vdw_clash_count(coords, znums)
        out["clash_vdw"] = int(cv)
        out["clash_severe"] = int(cs)
        out["worst_overlap"] = round(float(worst), 4)

        # --- the arm the lever's own docstring asks for -------------------------------
        # `passed` above is CIRCULAR under OIN_ACCEPT_SCORED: the lever accepts a conformer
        # on exactly the predicate `passed` scores with, so by construction it cannot detect
        # what dropping the independent re-perception costs. This runs the dropped test: a
        # full XYZToSMILES().convert() of the generated XYZ, which re-perceives connectivity
        # from coordinates alone rather than reusing the generator's own bond graph.
        # Timed separately -- it is measurement, not generation, and must not enter elapsed_s.
        ti = time.monotonic()
        try:
            with tempfile.TemporaryDirectory() as td:
                gp = os.path.join(td, "gen.xyz")
                with open(gp, "w") as fh:
                    fh.write(res.xyz)
                oin_indep = XYZToSMILES().convert(gp)
            out["oin_indep"] = oin_indep
            out["sha_indep"] = _sha(oin_indep)
            out["indep_passed"] = canonical_roundtrip_key(oin_in) == canonical_roundtrip_key(
                oin_indep
            )
            out["indep_byte_equal"] = oin_indep == oin_in
        except Exception as e:
            out["indep_passed"] = False
            out["indep_error"] = f"{type(e).__name__}: {e}"
        out["indep_s"] = round(time.monotonic() - ti, 2)

        if dump_dir:
            # The accepted conformer, plus the generator's CLAIM about which atoms are bonded
            # to the metal. The claim is a REFERENCE only -- never a measurement of attachment
            # (§6.1: a detached ligand keeps its bond object). Written so the falsification can
            # be re-scored offline against any candidate predicate.
            try:
                os.makedirs(dump_dir, exist_ok=True)
                stem = os.path.join(dump_dir, out["molecule"])
                with open(stem + ".xyz", "w") as fh:
                    fh.write(res.xyz)
                from oinsmiles.utils.perception_tmc import TRANSITION_METALS_NUM

                midx = next(
                    (
                        a.GetIdx()
                        for a in mol.GetAtoms()
                        if a.GetAtomicNum() in TRANSITION_METALS_NUM
                    ),
                    None,
                )
                claim = (
                    sorted(b.GetOtherAtomIdx(midx) for b in mol.GetAtomWithIdx(midx).GetBonds())
                    if midx is not None
                    else []
                )
                with open(stem + ".claim.json", "w") as fh:
                    json.dump(
                        {
                            "molecule": out["molecule"],
                            "input_xyz": xyz_path,
                            "metal_idx": midx,
                            "claimed_donors": claim,
                            "claimed_elements": [mol.GetAtomWithIdx(i).GetSymbol() for i in claim],
                            "natoms": mol.GetNumAtoms(),
                            "oin_in": out.get("oin_in"),
                            "oin_out": out.get("oin_out"),
                            "oin_indep": out.get("oin_indep"),
                            "indep_passed": out.get("indep_passed"),
                            "passed": out.get("passed"),
                        },
                        fh,
                        indent=1,
                    )
            except Exception as e:  # a dump failure must never change a verdict
                out["dump_error"] = f"{type(e).__name__}: {e}"
    except Exception as e:
        out.setdefault("elapsed_s", round(time.monotonic() - t0, 2))
        out["passed"] = False
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def run_arm(
    cohort: list[dict],
    arm: str,
    env_extra: dict,
    timeout: float,
    workers: int,
    hard_cap: float,
    dump_root: str | None = None,
) -> list:
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
            cmd = [
                sys.executable,
                os.path.abspath(__file__),
                "--one",
                rec["xyz"],
                "--timeout",
                str(timeout),
            ]
            if dump_root:
                cmd += ["--dump-xyz", os.path.join(dump_root, arm)]
            p = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                env=env,
                text=True,
            )
            running.append((p, rec, time.monotonic()))
        time.sleep(0.5)
        for p, rec, started in list(running):
            # HARD wall-clock kill. The generator's `timeout` is ADVISORY, not a bound:
            # embed_time_budget bounds the embed attempt loop, not the OIN-direct assembly
            # around it (measured: 60s requested, 60.7-137.9s spent, GOHWOQ 2.3x over). The
            # first run of this script relied on `timeout` and sat on one molecule for 30+
            # minutes. Only the harness's SIGKILL subprocess really enforces a budget, so this
            # reproduces that. Applied identically to both arms so the comparison stays fair.
            if p.poll() is None:
                if time.monotonic() - started > hard_cap:
                    p.kill()
                    print(
                        f"  [{arm}] {rec['mol']:16s} KILLED at hard cap {hard_cap:.0f}s", flush=True
                    )
                continue
            running.remove((p, rec, started))
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
                f"indep={r.get('indep_passed', '-')} sha={r.get('sha_out', '-')} "
                f"clash={r.get('clash_vdw', '-')} worst={r.get('worst_overlap', '-')}",
                flush=True,
            )
    return results


def summarize(rows: list[dict], label: str) -> dict:
    done = [r for r in rows if isinstance(r.get("elapsed_s"), (int, float))]
    el = [r["elapsed_s"] for r in done]
    cl = [r["clash_vdw"] for r in rows if isinstance(r.get("clash_vdw"), int)]
    sev = [r["clash_severe"] for r in rows if isinstance(r.get("clash_severe"), int)]
    wo = [r["worst_overlap"] for r in rows if isinstance(r.get("worst_overlap"), float)]
    s = {
        "arm": label,
        "n": len(rows),
        "passed": sum(1 for r in rows if r.get("passed")),
        "median_s": round(statistics.median(el), 2) if el else None,
        "mean_s": round(statistics.mean(el), 2) if el else None,
        "total_s": round(sum(el), 1) if el else None,
        "over_30s": sum(1 for v in el if v > 30),
        "clash_measured_on": len(cl),
        "clash_total": sum(cl) if cl else 0,
        "clash_mols_with_any": sum(1 for v in cl if v > 0),
        "clash_severe_total": sum(sev) if sev else 0,
        "worst_overlap_min": round(min(wo), 4) if wo else None,
        "worst_overlap_median": round(statistics.median(wo), 4) if wo else None,
        # G2: the independent re-perception the lever drops. `indep_measured_on` is the
        # DENOMINATOR and is reported for the same reason `clash_measured_on` is -- "0 fail"
        # and "0 fail over 0 measured" printed identically once already in this script.
        "indep_measured_on": sum(1 for r in rows if "indep_passed" in r),
        "indep_passed": sum(1 for r in rows if r.get("indep_passed")),
        # molecules the circular metric credits but independent re-perception does not:
        # the true price of the lever.
        "circular_only": sum(1 for r in rows if r.get("passed") and not r.get("indep_passed")),
        "byte_equal_in_out": sum(1 for r in rows if r.get("byte_equal_in_out")),
    }
    return s


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--one", help="measure a single xyz in this process and print JSON")
    ap.add_argument("--cohort", help="JSON with {'eta': [...], 'control': [...]} records")
    ap.add_argument("--out")
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument(
        "--hard-cap",
        type=float,
        default=0.0,
        help="wall-clock SIGKILL per molecule; default 1.6x --timeout, since the "
        "generator's own timeout is advisory (see run_arm)",
    )
    ap.add_argument(
        "--single-arm",
        choices=["0", "1"],
        help="run ONE arm with OIN_ACCEPT_SCORED set to this, and dump raw rows. Exists for "
        "the A-vs-A' determinism control: an A/B sha difference is only attributable to the "
        "lever if sha_out is stable across two runs of the SAME arm.",
    )
    ap.add_argument("--label", default="", help="arm label for --single-arm output")
    ap.add_argument(
        "--dump-xyz",
        help="persist each ACCEPTED conformer's coordinates + the generator's claimed metal "
        "donor set under this directory (one subdir per arm). Purely additive; no verdict in "
        "this script reads it. Exists so the §6.5 attachment falsification can be re-scored "
        "offline instead of regenerating the cohort once per candidate predicate.",
    )
    args = ap.parse_args()

    if args.one:
        print(json.dumps(measure_one(args.one, args.timeout, args.dump_xyz)))
        return 0

    if not args.cohort:
        ap.error("need --cohort or --one")
    raw = json.load(open(args.cohort))
    cohort = list(raw.get("eta", [])) + list(raw.get("control", []))
    cohort = [r for r in cohort if os.path.exists(r["xyz"])]
    print(f"cohort: {len(cohort)} molecules, {args.workers} workers, timeout {args.timeout}s")
    hard_cap = args.hard_cap if args.hard_cap else args.timeout * 1.6

    if args.single_arm is not None:
        label = args.label or f"arm-{args.single_arm}"
        rows = run_arm(
            cohort,
            label,
            {"OIN_ACCEPT_SCORED": args.single_arm},
            args.timeout,
            args.workers,
            hard_cap,
            args.dump_xyz,
        )
        s = summarize(rows, label)
        print(f"\n  {label}: {json.dumps(s)}")
        if args.out:
            with open(args.out, "w") as fh:
                json.dump({"arm": label, "rows": rows, "summary": s}, fh, indent=1)
            print(f"  wrote {args.out}")
        return 0

    print("\n--- ARM A: default (independent confirm ON) ---")
    a = run_arm(
        cohort,
        "A-default",
        {"OIN_ACCEPT_SCORED": "0"},
        args.timeout,
        args.workers,
        hard_cap,
        args.dump_xyz,
    )
    print("\n--- ARM B: OIN_ACCEPT_SCORED=1 ---")
    b = run_arm(
        cohort,
        "B-scored",
        {"OIN_ACCEPT_SCORED": "1"},
        args.timeout,
        args.workers,
        hard_cap,
        args.dump_xyz,
    )

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
            f"clash {s['clash_total']} over {s['clash_mols_with_any']}/"
            f"{s['clash_measured_on']} mols (severe {s['clash_severe_total']}, "
            f"worst_overlap min {s['worst_overlap_min']} med {s['worst_overlap_median']})"
        )
        print(
            f"  {'':10s} INDEPENDENT re-perception pass {s['indep_passed']}/"
            f"{s['indep_measured_on']} measured  |  credited by the circular metric but NOT "
            f"by independent re-perception: {s['circular_only']}  |  "
            f"byte-equal smiles_2==smiles_1: {s['byte_equal_in_out']}"
        )
    print(f"  PASS REGRESSIONS (A pass -> B fail): {regressions or 'none'}")
    print(f"  PASS FIXES       (A fail -> B pass): {fixes or 'none'}")

    indep_reg = [
        r["molecule"]
        for r in b
        if not r.get("indep_passed") and by.get(r["molecule"], {}).get("indep_passed")
    ]
    indep_fix = [
        r["molecule"]
        for r in b
        if r.get("indep_passed") and not by.get(r["molecule"], {}).get("indep_passed")
    ]
    print(f"  INDEP REGRESSIONS (A indep-pass -> B indep-fail): {indep_reg or 'none'}")
    print(f"  INDEP FIXES       (A indep-fail -> B indep-pass): {indep_fix or 'none'}")

    # --- G3: the wave's actual gate. Byte identity of smiles_2, per molecule. -----------
    # sha_in is the CONTROL: the lever is generation-side, so it must not move.
    print("\n---- G3  sha256(smiles_2) BYTE IDENTITY  (load-independent) ----")
    print(f"  {'molecule':18s} {'sha_A(out)':18s} {'sha_B(out)':18s} equal?  sha_in equal?")
    same_out = diff_out = only_one = neither = 0
    ctrl_moved = []
    for r in b:
        m = r["molecule"]
        ra = by.get(m, {})
        sa_o, sb_o = ra.get("sha_out"), r.get("sha_out")
        # Only a control VIOLATION when both arms actually produced an input hash. A
        # SIGKILLed molecule yields a synthesized row with no `oin_in` at all (exit -9), and
        # comparing None against a real hash reports "MOVED" for what is merely a missing
        # measurement. That false positive would condemn the whole run as confounded.
        ia, ib = ra.get("sha_in"), r.get("sha_in")
        if ia is None or ib is None:
            ctrl = "n/a (one arm produced no input hash)"
        elif ia == ib:
            ctrl = "OK"
        else:
            ctrl = "*** MOVED ***"
            ctrl_moved.append(m)
        if sa_o and sb_o:
            eq = "YES" if sa_o == sb_o else "NO"
            same_out += sa_o == sb_o
            diff_out += sa_o != sb_o
        elif sa_o or sb_o:
            eq = "one-arm"
            only_one += 1
        else:
            eq = "neither"
            neither += 1
        print(f"  {m:18s} {str(sa_o):18s} {str(sb_o):18s} {eq:7s} {ctrl}")
    print(
        f"  byte-identical smiles_2: {same_out} | DIFFERENT: {diff_out} | "
        f"only one arm produced a string: {only_one} | neither: {neither}"
    )
    if ctrl_moved:
        print(f"  *** CONTROL VIOLATION: sha_in moved on {ctrl_moved} -- run is confounded ***")

    out = {
        "arm_a": a,
        "arm_b": b,
        "summary_a": sa,
        "summary_b": sb,
        "regressions": regressions,
        "fixes": fixes,
        "indep_regressions": indep_reg,
        "indep_fixes": indep_fix,
        "sha_out_same": same_out,
        "sha_out_diff": diff_out,
        "sha_out_one_arm": only_one,
        "sha_out_neither": neither,
        "sha_in_control_moved": ctrl_moved,
    }
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=1)
        print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
