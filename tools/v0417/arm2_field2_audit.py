"""ARM 2 golden re-freeze, step 1: WHICH rows moved, and was it this lever? Encode-only, all rows.

An ARM 2 golden row is ``name, sha256(smiles_1), sha256(smiles_2), ...``. Field 2 is a pure function
of the input file and the encoder, so it can be audited for EVERY row in minutes, with no generator.
That matters because the precedent (v0.4.14) re-ran only the rows a sweep predicted would move, and
a predicted list cannot see a row that was already wrong.

Each molecule is encoded twice, each time in a FRESH PROCESS (perception memos are process-lifetime,
and several levers are read at import):

    shipped   no ``OIN_*`` variable set at all -- what a user gets
    off       ``OIN_EXACT_DONOR_FOLD=0``       -- the encoder as it was before v0.4.17

and each row lands in exactly one class:

    SAME             golden == off == shipped        nothing to do
    MOVED_BY_LEVER   golden == off != shipped        re-freeze; the reason is the exact fold
    STALE            golden != off == shipped        the golden was ALREADY wrong before this lane
    STALE_AND_MOVED  golden != off != shipped        both
    HEALED           golden != off, shipped == golden  stale before, and the lever moved it BACK:
                                                     field 2 now needs no re-freeze (NASZOY)
    SENTINEL         the golden holds NO_ENCODE@Ns   no string to compare
    ENCODE_FAILED    an arm raised or hit --timeout  reported, never silently dropped

WHAT A BROKEN VERSION OF THIS WOULD PRINT. A worker that imported the MAIN checkout's ``oinsmiles``
would report ``shipped == off`` on every row -- a flawless "nothing moved". So every worker reports
the ``oinsmiles`` it imported, the driver refuses any that is not this checkout's ``src/``, and the
run ABORTS if the lever fired on zero rows.

    <main>/.venv/bin/python tools/v0417/arm2_field2_audit.py \
        --golden tools/gate_v049_arm2_golden.tsv --cohort-dir <data>/cohort-v049-strata --out <jsonl>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve()
SRC = HERE.parents[2] / "src"
LEVER = "OIN_EXACT_DONOR_FOLD"


def _worker(xyz: str):
    sys.path.insert(0, str(SRC))  # beats PYTHONPATH, which is the point
    import oinsmiles
    from oinsmiles import XYZToSMILES

    out = {"src": str(Path(oinsmiles.__file__).resolve().parents[1])}
    try:
        s = XYZToSMILES().convert(xyz)
        out.update(oin=s, sha=hashlib.sha256(s.encode()).hexdigest(), n=len(s))
    except Exception as e:  # an encoder error is a result
        out["error"] = f"{type(e).__name__}:{e}"
    print(json.dumps(out), flush=True)


def _encode(xyz: Path, arm: str, timeout: float, lever: str = LEVER):
    env = {k: v for k, v in os.environ.items() if not k.startswith("OIN_")}
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    if arm == "off":
        env[lever] = "0"
    try:
        p = subprocess.run(
            [sys.executable, str(HERE), "--one", str(xyz)],
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"error": f"TIMEOUT@{timeout:.0f}s"}
    lines = [ln for ln in p.stdout.splitlines() if ln.startswith("{")]
    if not lines:
        return {"error": f"no output (rc={p.returncode}): {p.stderr.strip()[-200:]}"}
    return json.loads(lines[-1])


def _classify(golden, off, shipped):
    if "@" in golden or golden.startswith("NO_"):
        return "SENTINEL"
    if "sha" not in off or "sha" not in shipped:
        return "ENCODE_FAILED"
    stale, moved = golden != off["sha"], off["sha"] != shipped["sha"]
    if stale and shipped["sha"] == golden:
        return "HEALED"
    return {
        (False, False): "SAME",
        (False, True): "MOVED_BY_LEVER",
        (True, False): "STALE",
        (True, True): "STALE_AND_MOVED",
    }[(stale, moved)]


CLASSES = (
    "SAME",
    "MOVED_BY_LEVER",
    "STALE",
    "STALE_AND_MOVED",
    "HEALED",
    "SENTINEL",
    "ENCODE_FAILED",
)


def _reclassify(path: Path):
    """The shas are the measurement; the class is derived from them, so it can be re-derived."""
    lines = path.read_text().splitlines()
    if not lines or not lines[-1].startswith("#DONE"):
        sys.exit(f"ABORT: {path} has no #DONE trailer")
    rows = [json.loads(ln) for ln in lines[:-1]]
    for r in rows:
        r["class"] = _classify(r["golden_sha1"], r["off"], r["shipped"])
    path.write_text("".join(json.dumps(r) + "\n" for r in rows) + lines[-1] + "\n")
    print(f"{path.name}: DENOMINATOR {len(rows)} rows (reclassified from stored shas)")
    for k in CLASSES:
        names = sorted(r["molecule"].replace("_comp_0", "") for r in rows if r["class"] == k)
        print(f"  {len(names):4d}  {k}" + (f"   {' '.join(names)}" if k != "SAME" else ""))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--one", help=argparse.SUPPRESS)
    ap.add_argument("--golden", type=Path)
    ap.add_argument("--cohort-dir", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument(
        "--lever",
        default=LEVER,
        help="the default-ON lever being promoted; the `off` arm sets it to 0 (default: %(default)s)",
    )
    ap.add_argument(
        "--reclassify",
        action="store_true",
        help="re-derive `class` from the shas already in --out; encodes nothing",
    )
    args = ap.parse_args()
    if args.one:
        return _worker(args.one)
    if args.reclassify:
        return _reclassify(args.out)

    rows = [
        ln.split("\t")
        for ln in args.golden.read_text().splitlines()
        if ln and not ln.startswith("#")
    ]
    if not rows:
        sys.exit("ABORT: the golden has no rows -- refusing an empty denominator")

    def one(row):
        xyz = args.cohort_dir / f"{row[0]}.xyz"
        if not xyz.exists():
            return row, {"error": "input missing"}, {"error": "input missing"}
        return (
            row,
            _encode(xyz, "off", args.timeout, args.lever),
            _encode(xyz, "shipped", args.timeout, args.lever),
        )

    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(one, rows))

    foreign = {r.get("src") for _row, off, shp in results for r in (off, shp)} - {None, str(SRC)}
    if foreign:
        sys.exit(f"ABORT: a worker imported oinsmiles from {foreign}, not {SRC}")

    tally = Counter()
    with open(args.out, "w") as fh:
        for row, off, shp in sorted(results, key=lambda t: t[0][0]):
            cls = _classify(row[1], off, shp)
            tally[cls] += 1
            fh.write(
                json.dumps(
                    {
                        "molecule": row[0],
                        "class": cls,
                        "golden_sha1": row[1],
                        "off": {k: v for k, v in off.items() if k != "src"},
                        "shipped": {k: v for k, v in shp.items() if k != "src"},
                    }
                )
                + "\n"
            )
        fh.write(f"#DONE {len(results)}\n")

    fired = tally["MOVED_BY_LEVER"] + tally["STALE_AND_MOVED"] + tally["HEALED"]
    print(f"{args.golden.name}: DENOMINATOR {len(results)} rows   src={SRC}")
    for k in CLASSES:
        names = [r[0].replace("_comp_0", "") for r, o, s in results if _classify(r[1], o, s) == k]
        print(f"  {tally[k]:4d}  {k}" + (f"   {' '.join(sorted(names))}" if k != "SAME" else ""))
    if not fired:
        sys.exit(
            "ABORT: the lever fired on 0 rows -- a dead lever and a clean audit print the same"
        )
    sys.stdout.flush()


if __name__ == "__main__":
    main()
