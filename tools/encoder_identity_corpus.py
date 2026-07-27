"""Corpus-scale encoder byte-identity gate: re-encode every INPUT, diff against a past sweep.

WHY
---
Any re-scoring of a stored sweep rests on one assumption: **the encoder has not moved since
that sweep ran.** If it has, the "correction" being measured is a mixture of two changes and
is not attributable to either.

``tools/gate_v047.sh arm1`` checks this on 61 fixtures. That is the right shape for a
per-commit gate and the wrong shape for underwriting a re-baseline: 61 molecules cannot rule
out drift on a class none of them exercise. This tool asks the same question of all 5000, by
re-encoding each report's ``input_xyz`` and diffing the result against the ``smiles_1`` that
sweep recorded.

It is the honest counterpart to the sweep's own claim of stability. A bucket report re-run
over stored JSON proves only that the CLASSIFIER did not move; it cannot see an encoder change,
because it never encodes anything.

Cost: ~1.06 s/molecule, so ~7 minutes for 5000 on 12 workers.

Usage:
    PYTHONPATH=src python tools/encoder_identity_corpus.py \\
        --results-dir tmCAT-tmPHOTO_xyz_dataset/results-v0.4.6-sweep \\
        --out /path/to/encoder_identity.jsonl
"""

import argparse
import glob
import json
import os
import signal
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))


class _EncodeTimeout(Exception):
    pass


def _alarm(_s, _f):
    raise _EncodeTimeout()


def _worker(job):
    mol, path, expected, timeout = job
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    from oinsmiles import XYZToSMILES

    t0 = time.time()
    got, err = None, None
    signal.signal(signal.SIGALRM, _alarm)
    signal.setitimer(signal.ITIMER_REAL, timeout)
    try:
        got = XYZToSMILES().convert(path)
    except _EncodeTimeout:
        err = f"TimeoutError: encode exceeded {timeout}s"
    except Exception as e:  # noqa: BLE001 - an encode failure is a datum
        err = f"{type(e).__name__}: {e}"
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    return {
        "molecule": mol,
        "identical": (err is None and got == expected),
        "expected": expected,
        "got": got,
        "error": err,
        "elapsed_s": round(time.time() - t0, 3),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--out", required=True, help="jsonl log; appended with a #DONE sentinel")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    jobs = []
    for f in sorted(glob.glob(os.path.join(args.results_dir, "individual_reports", "*.json"))):
        rep = json.load(open(f))
        s1, xyz = rep.get("smiles_1"), rep.get("input_xyz")
        if s1 and xyz and os.path.exists(xyz):
            jobs.append((rep["molecule"], xyz, s1, args.timeout))
    if args.limit:
        jobs = jobs[: args.limit]
    print(f"re-encoding {len(jobs)} inputs on {args.workers} workers", flush=True)

    same = drift = failed = 0
    t0 = time.time()
    with open(args.out, "w") as log, ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_worker, j): j[0] for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            try:
                rec = fut.result()
            except Exception as e:  # noqa: BLE001
                rec = {"molecule": futs[fut], "identical": False, "error": f"worker died: {e}"}
            log.write(json.dumps(rec) + "\n")
            log.flush()
            if rec["identical"]:
                same += 1
            elif rec.get("error"):
                failed += 1
            else:
                drift += 1
            if i % 500 == 0:
                print(f"  {i}/{len(jobs)}  same={same} drift={drift} err={failed}", flush=True)
        log.write(f"#DONE {len(jobs)}\n")
        log.flush()

    print(f"\n{'=' * 66}")
    print(f"byte-identical : {same}/{len(jobs)}  ({100 * same / max(len(jobs), 1):.2f}%)")
    print(f"DRIFTED        : {drift}")
    print(f"encode errors  : {failed}   (these re-encode differently only because they raise)")
    print(f"elapsed        : {time.time() - t0:.0f}s")
    if drift:
        print("\n⚠ THE ENCODER MOVED. Any re-score of this sweep mixes two changes.")
        sys.exit(1)


if __name__ == "__main__":
    main()
