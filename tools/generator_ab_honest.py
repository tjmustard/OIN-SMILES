#!/usr/bin/env python
"""Does an encoder-side lever survive REAL GENERATION? (v0.4.14)

WHY THIS EXISTS, AND WHY THE THREE INSTRUMENTS ALREADY IN tools/ CANNOT ANSWER IT
================================================================================
v0.4.13 and v0.4.14 both measured an encoder lever by OFFLINE RE-SCORE: take a frozen sweep's
stored generated structures, re-encode them under each lever setting, re-classify. The licence
for that is ``tools/fold_key_invariance.py`` reading **0 comparison keys changed**, which is
argued to mean the generator returns the same conformers either way.

**That argument has a hole, and this tool exists because the hole is real.** ``accept_fn`` decides
by comparison KEY, so a key-invariant lever cannot change *which* conformer is accepted from a
given pool. But the generator's INPUT is the OIN **string**, not the key:

    gen.generate(oin_in)  ->  generation/oin_parser.OINParser  ->  ParsedOIN
                          ->  metallogen_adapter: m-SMILES, CoordMap pins donors to vertices

A slot relabeling changes that string. Different ``ParsedOIN`` -> different CoordMap -> a
different embedding. **Key-invariance bounds the ACCEPTANCE step; it says nothing about the
EMBEDDING step**, and the offline re-score holds the embedding fixed by construction.

So a lever can read GENERATOR_NEUTRAL, produce a clean one-directional offline transition matrix,
and still change what ships.

THE OTHER TWO INSTRUMENTS ARE ALSO WRONG FOR THIS
=================================================
* ``tools/gate_v047.sh arm2`` re-generates, but scores with
  ``get_oin_string(result.mol, coords)`` -- the GENERATOR'S OWN bond graph
  (``gate_arm2_roundtrip_one.py``). That is the circular predicate ``OIN_INDEP_SCORE`` replaced in
  v0.4.8 at a measured **9.6% false-positive rate (28.1% on haptic inputs)**. It is a byte-identity
  gate, not an accuracy instrument, and reading its hashes as a round-trip verdict is a mistake
  this file exists partly to stop.
* ``tools/resonance_transition_sim.py`` is the offline re-score itself.

WHAT THIS DOES
==============
For each molecule and each lever setting, the full shipped path, scored HONESTLY:

    oin_in  = XYZToSMILES().convert(input.xyz)          # encoder, under this lever setting
    res     = generate(oin_in)                          # the generator actually runs
    oin_out = XYZToSMILES().convert(written res.xyz)    # RE-PERCEIVED from coordinates
    byte_exact  <=>  oin_in == oin_out

``oin_out`` is never taken from ``res.mol``. That is the whole point.

DETERMINISM, AND WHY THIS IS NOT THE STOCHASTIC-A/B TRAP
========================================================
This project's standing rule is *never A/B by re-running a stochastic harness*. It does not bind
here, and that was checked rather than assumed: ``MetalloGenAdapter`` takes ``seed=42`` and the
same molecule re-run under the same lever setting reproduces its hashes exactly (measured on
``HEKFEL_comp_0``, both arms, two runs each). ``--repeat`` re-runs every arm N times and reports
any molecule that disagrees with itself, so the assumption is re-checked on every cohort rather
than inherited from that one measurement.

⚠ What is NOT deterministic is the pool SIZE: ``OIN3DGenerator(timeout=)`` is advisory and the
embed loop checks its deadline between attempts, so a heavily loaded box can produce a smaller
pool. Run this on a quiet machine, and treat a ``--repeat`` disagreement as a load artifact
before treating it as a code defect.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/generator_ab_honest.py \\
        --cohort-dir <MAIN>/tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k \\
        --molecules-file mols.txt --lever OIN_RESONANCE_DONOR_FOLD \\
        --out-json gen_ab.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import time
from collections import Counter

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def _run_one(xyz_path, lever, on, timeout):
    """``(oin_in, oin_out, elapsed_s, error)`` for one molecule at one lever setting."""
    os.environ[lever] = "1" if on else "0"

    from oinsmiles import XYZToSMILES
    from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen

    t0 = time.monotonic()
    try:
        oin_in = XYZToSMILES().convert(xyz_path)
    except Exception as exc:  # noqa: BLE001 -- an unencodable input is data, not a crash
        return None, None, round(time.monotonic() - t0, 2), f"encode_in:{type(exc).__name__}"

    try:
        gen = OIN3DGeneratorMetallogen(
            optimizer=None, ensemble_size=1, timeout=timeout, ff_params=None
        )
        res = gen.generate(oin_in)
    except Exception as exc:  # noqa: BLE001 -- a failed generation is data
        return oin_in, None, round(time.monotonic() - t0, 2), f"generate:{type(exc).__name__}"

    xyz = getattr(res, "xyz", None)
    if not xyz:
        return oin_in, None, round(time.monotonic() - t0, 2), "no_structure"

    fh = tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False)
    try:
        fh.write(xyz)
        fh.close()
        # 🔴 HONEST. Never `get_oin_string(res.mol, coords)` -- that is the generator's own bond
        # graph, i.e. exactly the artifact that would have to be wrong for this test to fail.
        oin_out = XYZToSMILES().convert(fh.name)
    except Exception as exc:  # noqa: BLE001
        return oin_in, None, round(time.monotonic() - t0, 2), f"encode_out:{type(exc).__name__}"
    finally:
        try:
            os.unlink(fh.name)
        except OSError:
            pass

    return oin_in, oin_out, round(time.monotonic() - t0, 2), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort-dir", required=True)
    ap.add_argument("--molecules-file", required=True, help="one molecule name per line")
    ap.add_argument("--lever", required=True)
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument(
        "--repeat", type=int, default=1, help="re-run each arm N times to check determinism"
    )
    ap.add_argument("--out-json")
    args = ap.parse_args()

    if not os.path.isdir(args.cohort_dir):
        sys.exit(f"🔴 --cohort-dir {args.cohort_dir} is not a directory")
    # ⚠ `#` comments are SKIPPED. The sample-membership files in `measurements/` carry their
    # provenance on a leading `#` line -- a rate without its sample is not reproducible -- and
    # without this the header would be read as a molecule name, land in `input_missing`, and
    # silently shift every count by one. That is the same failure mode as v0.4.14's missing
    # trailing newline, which merged two names and dropped the one molecule a run existed to
    # settle.
    names = [
        n.strip() for n in open(args.molecules_file) if n.strip() and not n.lstrip().startswith("#")
    ]
    if not names:
        sys.exit("🔴 REFUSING: molecules file is empty")

    rows, nondet = [], []
    for i, mol in enumerate(names, 1):
        p = os.path.join(args.cohort_dir, f"{mol}.xyz")
        if not os.path.exists(p):
            rows.append({"molecule": mol, "error": "input_missing"})
            continue
        arm = {}
        for on in (False, True):
            reps = [_run_one(p, args.lever, on, args.timeout) for _ in range(args.repeat)]
            if args.repeat > 1 and len({(r[0], r[1]) for r in reps}) > 1:
                nondet.append(f"{mol}:{'on' if on else 'off'}")
            arm["on" if on else "off"] = reps[0]
        (i1, o1, e1, err1) = arm["off"]
        (i2, o2, e2, err2) = arm["on"]
        rows.append(
            {
                "molecule": mol,
                "off": {"byte_exact": bool(i1 and o1 and i1 == o1), "elapsed_s": e1, "error": err1},
                "on": {"byte_exact": bool(i2 and o2 and i2 == o2), "elapsed_s": e2, "error": err2},
                "input_string_moved": bool(i1 and i2 and i1 != i2),
                "output_string_moved": bool(o1 and o2 and o1 != o2),
            }
        )
        r = rows[-1]
        print(
            f"[{i}/{len(names)}] {mol:20s} byte_exact {r['off']['byte_exact']!s:5s} -> "
            f"{r['on']['byte_exact']!s:5s}  in_moved={r['input_string_moved']} "
            f"out_moved={r['output_string_moved']}",
            flush=True,
        )

    scored = [r for r in rows if "off" in r]
    if not scored:
        sys.exit(f"\n🔴 REFUSING TO REPORT: 0 of {len(names)} molecules produced a verdict.")

    gained = [r["molecule"] for r in scored if not r["off"]["byte_exact"] and r["on"]["byte_exact"]]
    lost = [r["molecule"] for r in scored if r["off"]["byte_exact"] and not r["on"]["byte_exact"]]
    tal = Counter((r["off"]["byte_exact"], r["on"]["byte_exact"]) for r in scored)

    print(f"\n=== HONEST GENERATOR A/B — {args.lever}, n={len(scored)} of {len(names)} ===")
    for (a, b), c in sorted(tal.items()):
        print(f"  byte_exact  OFF={a!s:5s} ON={b!s:5s}  {c:4d}")
    print(f"\n  REAL gains  (OFF fail -> ON pass): {len(gained)}  {gained[:10]}")
    print(f"  REAL losses (OFF pass -> ON fail): {len(lost)}  {lost[:10]}")
    print(f"  net: {len(gained) - len(lost):+d} molecules on this sample")
    print(
        f"\n  input string moved by the lever : {sum(1 for r in scored if r['input_string_moved'])}"
    )
    print(
        f"  GENERATED output moved           : {sum(1 for r in scored if r['output_string_moved'])}"
    )
    print(
        "\n  ⚠ 'GENERATED output moved' > 0 is the finding the offline re-score cannot make:\n"
        "    the lever changed the string the GENERATOR CONSUMES, so it embedded differently.\n"
        "    Key-invariance bounds acceptance, not embedding."
    )
    if nondet:
        print(f"\n  🔴 NON-DETERMINISTIC across --repeat: {len(nondet)} arms -- {nondet[:8]}")
        print("     Check machine load before reading this as a code defect.")

    out = {
        "lever": args.lever,
        "cohort_dir": args.cohort_dir,
        "n_requested": len(names),
        "n_scored": len(scored),
        "gains": gained,
        "losses": lost,
        "net": len(gained) - len(lost),
        "input_moved": sum(1 for r in scored if r["input_string_moved"]),
        "output_moved": sum(1 for r in scored if r["output_string_moved"]),
        "nondeterministic": nondet,
        "rows": rows,
    }
    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nwrote {args.out_json}")


if __name__ == "__main__":
    main()
