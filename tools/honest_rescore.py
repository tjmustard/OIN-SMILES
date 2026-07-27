"""Offline HONEST re-score of a completed sweep -- no generator run required.

WHY THIS EXISTS
---------------
``tools/test_dataset_roundtrip.py`` scores a round trip with
``get_oin_string(gen_result.mol, coords)`` -- the GENERATOR's own bond graph. That object
carries the bonds the generator *intended*, not the ones the coordinates support, so the
metric asserts coordination the geometry lacks (61 false positives / 9.6%, 28.1% of haptic
inputs) and drops stereo the geometry has (8 false negatives). It is not a bug so much as a
circularity: ``gen_result.mol`` is exactly the artifact that would have to be wrong for the
test to fail.

The fix is one call -- a full ``XYZToSMILES().convert()`` of the generated XYZ, which
re-perceives bonds AND stereo from coordinates alone. That is what ``OIN_INDEP_SCORE``
records in the live harness.

THE POINT OF THIS TOOL: that call does not need the generator.
``save_artifacts`` writes ``structures/<mol>_generated.xyz`` from the SAME ``gen_result.xyz``
string that ``_attempt_generation`` writes to the temp path the lever converts
(``test_dataset_roundtrip.py`` line ~167 vs ~221). Re-encoding the stored file is therefore
**bit-identical to what the lever would have computed**, not an approximation of it.

Consequences, both of which matter more than the speedup:

1. **Cost.** Measured 0.33 s/molecule. The full 4688-structure corpus re-scores in ~5
   minutes, against the ~55 CPU-hours a live re-sweep costs.
2. **The A/B confound disappears.** Conformers are held fixed by construction, so
   ``smiles_1`` cannot move, no second encode enters the ``--mol-timeout`` budget, and no
   marginal pass converts into a timeout for reasons unrelated to honesty. A live re-sweep
   cannot make that guarantee; this can.

WHAT IT DOES NOT DO
-------------------
It cannot detect a change in the GENERATOR, because it never runs one -- it re-scores the
conformers a past sweep produced. Pair it with a bounded live arm when the generator itself
may have moved. It also cannot help the molecules that never produced a structure at all
(timeouts, no-conformer): those are failures in both arms and are passed through untouched.

The source results directory is opened READ-ONLY and never modified. Output is a re-scored
twin containing ``individual_reports/`` plus the sidecar log, which
``tools/roundtrip_bucket_report.py --score honest`` reads directly.

Usage:
    PYTHONPATH=src python tools/honest_rescore.py \\
        --results-dir tmCAT-tmPHOTO_xyz_dataset/results-v0.4.6-sweep \\
        --output-dir  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest \\
        --workers 12
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

from oinsmiles.oin.compare import canonical_roundtrip_key  # noqa: E402

#: Transition labels. The pair (scored verdict -> honest verdict) is the unit of the
#: correction: reporting only "how many moved" hides that a false positive can land in
#: ``key_equal`` instead of falling out of ``byte_exact`` altogether, which is the single
#: most misread thing about a re-baseline.
CLASSES = (
    "byte->byte",
    "byte->key",
    "byte->FAIL",
    "key->byte",
    "key->key",
    "key->FAIL",
    "FAIL->byte",
    "FAIL->key",
    "fail->fail",
    "no_structure",
    "indep_encode_FAILED",
)


class _EncodeTimeout(Exception):
    pass


def _alarm(_signum, _frame):
    raise _EncodeTimeout()


def _verdict(s1, s2):
    """('byte' | 'key' | 'FAIL') for one pair of OIN strings."""
    if s1 is None or s2 is None:
        return "FAIL"
    if s1 == s2:
        return "byte"
    try:
        if canonical_roundtrip_key(s1) == canonical_roundtrip_key(s2):
            return "key"
    except Exception:
        return "FAIL"
    return "FAIL"


def _worker(job):
    """Re-encode one stored generated structure. Runs in a child process.

    RDKit's logger is disabled HERE rather than at import in the parent: the pool's children
    are forked/spawned per worker and a parent-side disable does not reliably survive. The
    spam is not merely noisy -- it is per-molecule and would drown the sentinel discipline
    this tool depends on. Note we silence RDKit's own logger, never the process's stderr:
    a real traceback must still reach the log file.
    """
    mol, struct_path, s1, s2, timeout = job
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    from oinsmiles import XYZToSMILES

    t0 = time.time()
    indep, err = None, None
    signal.signal(signal.SIGALRM, _alarm)
    signal.setitimer(signal.ITIMER_REAL, timeout)
    try:
        indep = XYZToSMILES().convert(struct_path)
    except _EncodeTimeout:
        err = f"TimeoutError: independent encode exceeded {timeout}s"
    except Exception as e:  # noqa: BLE001 - any encoder failure is a datum, not a crash
        err = f"{type(e).__name__}: {e}"
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)

    scored = _verdict(s1, s2)
    honest = _verdict(s1, indep)
    if err is not None:
        cls = "indep_encode_FAILED"
    elif scored == "FAIL" and honest == "FAIL":
        cls = "fail->fail"
    else:
        cls = f"{scored}->{honest}"
    return {
        "molecule": mol,
        "smiles_2_indep": indep,
        "indep_key_match": (None if indep is None else honest in ("byte", "key")),
        "indep_error": err,
        "honest_class": cls,
        "scored_verdict": scored,
        "honest_verdict": honest,
        "indep_encode_s": round(time.time() - t0, 3),
    }


def _load_done(log_path):
    """Molecules already re-scored, from a previous (possibly killed) run.

    Only lines that parse are trusted. A run killed mid-write leaves a truncated final line;
    dropping it silently is correct -- that molecule simply gets re-done.
    """
    done = {}
    if not os.path.exists(log_path):
        return done
    with open(log_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                rec = json.loads(line)
            except Exception:
                continue
            if "molecule" in rec:
                done[rec["molecule"]] = rec
    return done


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--results-dir", required=True, help="completed sweep; opened READ-ONLY")
    ap.add_argument("--output-dir", required=True, help="re-scored twin to write")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--timeout", type=float, default=120.0, help="per-molecule encode budget")
    ap.add_argument("--limit", type=int, default=0, help="debug: stop after N molecules")
    ap.add_argument("--resume", action="store_true", help="skip molecules already in the log")
    ap.add_argument(
        "--fill-coordination",
        action="store_true",
        help=(
            "also compute report['coordination'] from the two geometries. Needed for sweeps "
            "predating the diagnostic (results-v0.4.6-sweep has no such field), and it is the "
            "independent cross-check on the honest verdict -- coordinate-only, ~2.2 ms, and it "
            "consults neither bond graph."
        ),
    )
    args = ap.parse_args()

    src = os.path.abspath(args.results_dir)
    out = os.path.abspath(args.output_dir)
    if src == out:
        sys.exit("refusing to write into the source sweep: --output-dir must differ")
    indiv_src = os.path.join(src, "individual_reports")
    struct_src = os.path.join(src, "structures")
    indiv_out = os.path.join(out, "individual_reports")
    os.makedirs(indiv_out, exist_ok=True)
    with open(os.path.join(out, "SOURCE"), "w") as fh:
        fh.write(f"{src}\nre-scored by tools/honest_rescore.py -- structures NOT copied\n")

    paths = sorted(glob.glob(os.path.join(indiv_src, "*.json")))
    if args.limit:
        paths = paths[: args.limit]
    reports = {}
    for p in paths:
        with open(p) as fh:
            rep = json.load(fh)
        reports[rep["molecule"]] = rep
    print(f"source reports: {len(reports)}", flush=True)

    log_path = os.path.join(out, "honest_rescore.jsonl")
    done = _load_done(log_path) if args.resume else {}
    if done:
        print(f"resuming: {len(done)} already scored", flush=True)

    jobs, no_struct = [], []
    for mol, rep in reports.items():
        if mol in done:
            continue
        sp = os.path.join(struct_src, f"{mol}_generated.xyz")
        if not os.path.exists(sp) or not rep.get("smiles_1"):
            no_struct.append(mol)
            continue
        jobs.append((mol, sp, rep.get("smiles_1"), rep.get("smiles_2"), args.timeout))
    print(f"to re-encode: {len(jobs)}   no stored structure: {len(no_struct)}", flush=True)

    # Every result line is a REAL append+flush, not a buffer inside a long-lived process.
    # A `timeout`/kill mid-run must not be able to discard output that would then look like
    # agreement once sorted -- the discipline tools/gate_v047.sh documents.
    results = dict(done)
    t0 = time.time()
    with open(log_path, "a") as log:
        for mol in no_struct:
            rec = {
                "molecule": mol,
                "smiles_2_indep": None,
                "indep_key_match": None,
                "indep_error": None,
                "honest_class": "no_structure",
                "scored_verdict": None,
                "honest_verdict": None,
                "indep_encode_s": 0.0,
            }
            results[mol] = rec
            log.write(json.dumps(rec) + "\n")
        log.flush()

        n = 0
        if jobs:
            with ProcessPoolExecutor(max_workers=args.workers) as ex:
                futs = {ex.submit(_worker, j): j[0] for j in jobs}
                for fut in as_completed(futs):
                    mol = futs[fut]
                    try:
                        rec = fut.result()
                    except Exception as e:  # noqa: BLE001 - a dead worker is a datum
                        rec = {
                            "molecule": mol,
                            "smiles_2_indep": None,
                            "indep_key_match": None,
                            "indep_error": f"worker died: {type(e).__name__}: {e}",
                            "honest_class": "indep_encode_FAILED",
                            "scored_verdict": None,
                            "honest_verdict": None,
                            "indep_encode_s": None,
                        }
                    results[mol] = rec
                    log.write(json.dumps(rec) + "\n")
                    log.flush()
                    n += 1
                    if n % 250 == 0:
                        print(f"  {n}/{len(jobs)}  ({time.time() - t0:.0f}s)", flush=True)

        # The sentinel and its denominator. Checked BEFORE any comparison is trusted: an
        # empty or truncated results file must never be able to look like consensus.
        log.write(f"#DONE {len(results)}\n")
        log.flush()

    if args.fill_coordination:
        from oinsmiles.oin.coordination import coordination_report

        n_filled = 0
        for mol, rep in reports.items():
            if rep.get("coordination") is not None:
                continue
            sp = os.path.join(struct_src, f"{mol}_generated.xyz")
            if not os.path.exists(sp) or not os.path.exists(rep.get("input_xyz", "")):
                continue
            try:
                with open(rep["input_xyz"]) as a, open(sp) as b:
                    rep["coordination"] = coordination_report(a.read(), b.read())
                n_filled += 1
            except Exception as e:  # noqa: BLE001 - a diagnostic must not break the run
                rep["coordination"] = {
                    "intact": None,
                    "reason": f"probe failed: {type(e).__name__}: {e}",
                }
        print(f"coordination backfilled for {n_filled} molecules", flush=True)

    for mol, rep in reports.items():
        rec = results.get(mol)
        if rec:
            rep["smiles_2_indep"] = rec["smiles_2_indep"]
            rep["indep_key_match"] = rec["indep_key_match"]
            rep["honest_class"] = rec["honest_class"]
            if rec.get("indep_error"):
                rep["indep_error"] = rec["indep_error"]
            rep["indep_encode_s"] = rec["indep_encode_s"]
        with open(os.path.join(indiv_out, f"{mol}.json"), "w") as fh:
            json.dump(rep, fh, indent=2)

    counts = {}
    for rec in results.values():
        counts[rec["honest_class"]] = counts.get(rec["honest_class"], 0) + 1
    times = [r["indep_encode_s"] for r in results.values() if r.get("indep_encode_s")]
    print(f"\nre-scored {len(results)} molecules in {time.time() - t0:.0f}s", flush=True)
    if times:
        print(f"encode cost: mean {sum(times) / len(times):.2f}s  max {max(times):.1f}s")
    print("\ntransition (scored -> honest):")
    for c in CLASSES:
        if counts.get(c):
            print(f"  {counts[c]:6d}  {c}")
    for c, v in sorted(counts.items()):
        if c not in CLASSES:
            print(f"  {v:6d}  {c}  (UNEXPECTED CLASS)")
    print(f"\nwrote {indiv_out}")

    # A copy of the structures dir would double 9673 files for nothing; downstream tools
    # read individual_reports/ only. Anything that needs the geometry follows SOURCE.
    if not os.path.exists(os.path.join(out, "structures")):
        try:
            os.symlink(struct_src, os.path.join(out, "structures"))
        except OSError as e:
            print(f"note: could not link structures/ ({e}); SOURCE records the path")


if __name__ == "__main__":
    main()
