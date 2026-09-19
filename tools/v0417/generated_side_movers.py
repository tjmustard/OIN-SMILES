"""Which molecules can the exact fold change the VERDICT of without changing the generator's input?

``exact_fold_live_report.py`` gives the INPUT-side movers: ``E(x)`` changes, so the generator is
handed a different string and anything can happen. But a verdict is ``E(x) == E(G(E(x)))``, and
the second encode moves too. A molecule whose input string is untouched still flips if the
re-encode of its generated structure does -- at least 94 of the 252 ``E1_NONCANONICAL`` failures
the lever repairs are of that kind. So the A/B population is input movers UNION these.

Derived from COORDINATES, not stored strings (v0.4.14's rule): every
``<sweep>/structures/*_generated.xyz`` is encoded twice, lever OFF and lever ON, by the encoder in
this checkout. The stored ``smiles_2_indep`` is read only as a cross-check on the OFF arm -- if
today's OFF encode disagrees with it on more than the odd molecule, the sweep of record no longer
describes this code and the A/B needs a fresh OFF arm (it gets one anyway).

Runner: forked workers, one ``#START``/result pair per task, kernel ``SIGALRM`` (``SIG_DFL``) as
the hard budget, parent reconciles a killed worker and respawns the remainder -- the pattern
``tools/census/e_selfconsistency.py`` ran 110,000 encodes through.

    PYTHONPATH=$PWD/src <main>/.venv/bin/python tools/v0417/generated_side_movers.py --cpu 10
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import signal
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "census"))  # sibling module only
import e_selfconsistency as esc  # noqa: E402

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
LEVER = "OIN_EXACT_DONOR_FOLD"
ARMS = ("off", "on")


def _worker(tasks, part, structures, timeout):
    signal.signal(signal.SIGALRM, signal.SIG_DFL)  # a Python handler cannot interrupt C++
    from oinsmiles import XYZToSMILES  # pay the import outside the first alarm

    with open(part, "a") as fh:
        for mol, arm in tasks:
            fh.write(f"#START {mol} {arm}\n")
            fh.flush()
            os.environ[LEVER] = "1" if arm == "on" else "0"  # never unset: "0" is OFF
            t0 = time.time()
            signal.alarm(timeout)
            try:
                with esc._silence_fds():
                    oin = XYZToSMILES().convert(str(structures / f"{mol}_generated.xyz"))
                err = None
            except Exception as exc:  # noqa: BLE001 -- an unencodable structure is data
                oin, err = None, f"{type(exc).__name__}: {exc}"[:200]
            signal.alarm(0)
            fh.write(
                json.dumps(
                    {
                        "molecule": mol,
                        "arm": arm,
                        "oin": oin,
                        "error": err,
                        "t": round(time.time() - t0, 3),
                    }
                )
                + "\n"
            )
            fh.flush()
        os.fsync(fh.fileno())
    os._exit(0)  # skip interpreter finalizers (BLAS thread pools)


def _reconcile(part, todo, timeout):
    done, started = set(), None
    if part.exists():
        for ln in part.read_text().splitlines():
            if ln.startswith("#START "):
                started = tuple(ln.split()[1:3])
            else:
                r = json.loads(ln)
                done.add((r["molecule"], r["arm"]))
    if started is not None and started not in done:
        with open(part, "a") as fh:
            fh.write(
                json.dumps(
                    {
                        "molecule": started[0],
                        "arm": started[1],
                        "oin": None,
                        "error": "TIMEOUT",
                        "t": float(timeout),
                    }
                )
                + "\n"
            )
        done.add(started)
    return [t for t in todo if t not in done]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--sweep", type=Path, default=MAIN / "results-v0.4.14-sweep")
    ap.add_argument("--out", type=Path, default=MAIN / "results-v0.4.17-exactfold")
    ap.add_argument("--cpu", type=int, default=10)
    ap.add_argument("--timeout", type=int, default=300, help="per-ENCODE hard kill, seconds")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    import oinsmiles

    print("oinsmiles from:", oinsmiles.__file__, flush=True)
    structures = args.sweep / "structures"
    mols = sorted(p.name.removesuffix("_generated.xyz") for p in structures.glob("*_generated.xyz"))
    if args.limit:
        mols = mols[: args.limit]
    if not mols:
        sys.exit(f"ABORT: no *_generated.xyz under {structures}")
    print(f"generated structures: {len(mols)}  encodes: {2 * len(mols)}", flush=True)

    name = "generated_side_movers"
    parts = [args.out / f".{name}.part{k}" for k in range(args.cpu)]
    for p in parts:
        p.unlink(missing_ok=True)
    todo = {k: [(m, a) for m in mols[k :: args.cpu] for a in ARMS] for k in range(args.cpu)}
    procs: dict = {}
    t0 = time.time()
    while any(todo.values()) or procs:
        for k in range(args.cpu):
            p = procs.get(k)
            if p is not None and p.is_alive():
                continue
            if p is not None:
                p.join()
                del procs[k]
                todo[k] = _reconcile(parts[k], todo[k], args.timeout)
            if todo[k]:
                procs[k] = mp.Process(
                    target=_worker, args=(todo[k], parts[k], structures, args.timeout)
                )
                procs[k].start()
        time.sleep(0.2)

    by: dict = {m: {} for m in mols}
    for p in parts:
        if p.exists():
            for ln in p.read_text().splitlines():
                if ln[0] != "#":
                    r = json.loads(ln)
                    by[r["molecule"]][r["arm"]] = r
    out = args.out / f"{name}.jsonl"
    with open(out, "w") as fh:
        for m in mols:
            off, on = by[m].get("off", {}), by[m].get("on", {})
            fh.write(
                json.dumps(
                    {
                        "molecule": m,
                        "off": off.get("oin"),
                        "on": on.get("oin"),
                        "err_off": off.get("error"),
                        "err_on": on.get("error"),
                        "t_off": off.get("t"),
                        "t_on": on.get("t"),
                    }
                )
                + "\n"
            )
        fh.write(f"#DONE {len(mols)}\n")
    for p in parts:
        p.unlink(missing_ok=True)
    print(f"wrote {out}  {time.time() - t0:.0f}s wall", flush=True)

    recs = [json.loads(ln) for ln in out.read_text().splitlines() if ln[0] == "{"]
    ok = [r for r in recs if r["off"] and r["on"]]
    moved = sorted(r["molecule"] for r in ok if r["off"] != r["on"])
    print(f"\nDENOMINATOR structures={len(recs)}  encoded in both arms={len(ok)}")
    print("errors:", dict(Counter((r["err_off"] or r["err_on"] or "ok")[:30] for r in recs)))
    print(f"GENERATED-SIDE MOVERS: {len(moved)}/{len(ok)}")
    print(
        f"CPU-h off={sum(r['t_off'] or 0 for r in recs) / 3600:.2f}  "
        f"on={sum(r['t_on'] or 0 for r in recs) / 3600:.2f}"
    )
    rep = args.sweep / "bucket_report_honest.json"
    if rep.exists():
        B = {r["molecule"].removesuffix(".xyz"): r for r in json.loads(rep.read_text())}
        key = "smiles_2_indep" if any("smiles_2_indep" in v for v in B.values()) else "smiles_2"
        cmp_ = [r for r in ok if B.get(r["molecule"], {}).get(key)]
        same = sum(1 for r in cmp_ if B[r["molecule"]][key] == r["off"])
        print(f"OFF arm == sweep-of-record {key}: {same}/{len(cmp_)}  (drift = code or load)")
    (args.out / "generated_side_movers.txt").write_text("\n".join(moved) + "\n")
    print("wrote", args.out / "generated_side_movers.txt")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
