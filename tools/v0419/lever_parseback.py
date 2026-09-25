"""v0.4.19 serializer lane, step 0: does a STRING-side change make the string describe the input?

The ruler here is the census's parse-back (``tools/census/string_sufficiency.py::parseback_one``):
the string is read back to a graph with NO 3D -- at parser level and at the generator adapter's
level -- and compared to the input xyz. It is EXACT for an encoder change (no generator, no seed,
no budget), and it is the ruler that assigns ``E1_HCOUNT`` / ``E1_GRAPH`` in the census. A round-trip
A/B cannot see a string-side repair through a generator fault, which is how OIN_H_FAITHFUL came to
be held off on "moves nothing" (levers.py) -- measured with the wrong ruler.

Arms (fresh process per molecule and arm; levers WRITTEN, never unset):

    shipped    the encoder as it ships
    hfaith     OIN_H_FAITHFUL=1                        (H1/H6 of the lane map)
    rc1off     OINDiscreteAligner._canonical_ring_signature -> None: the aligner's RC1 rank swap
               among same-mass eta fragments takes its own fail-safe            (H2, diagnostic)
    both       hfaith + rc1off
    rc1prop    OIN_RC1_PROPAGATE=1: the swap is kept and PROPAGATED to get_oin_string  (H2, the fix)
    fix2       hfaith + rc1prop -- the two levers together
    cap        OIN_CAP_IGNORES_METAL=1: a ligand atom's valence cap ignores its metal contact  (H4)
    fix3       hfaith + rc1prop + cap
    n2 / e2b / e2c   E2 lane: OIN_N_VALENCE_2 / + OIN_CANONICAL_RESONANCE / + OIN_CANONICAL_CHARGES

Per arm and population it prints: strings that CHANGED vs shipped (a canonicality lever must not
move a verified pass), and the parse-back classes -- adapter graph ISO / adapter H SAME -- with the
transitions vs shipped. A population's verified passes are the CONTROL: a repair that moves them
off ISO/SAME is a regression whatever it does elsewhere.

    <main>/.venv/bin/python tools/v0419/lever_parseback.py --names <file> --classes <json> \\
        --arms shipped,hfaith,rc1off,both --out <jsonl> [--jobs 6]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve()
SRC = HERE.parents[2] / "src"
COHORT = Path(
    "/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k"
)
ARMS = {
    "shipped": ({}, False),
    "hfaith": ({"OIN_H_FAITHFUL": "1"}, False),
    "rc1off": ({}, True),
    "both": ({"OIN_H_FAITHFUL": "1"}, True),
    "rc1prop": ({"OIN_RC1_PROPAGATE": "1"}, False),
    "fix2": ({"OIN_H_FAITHFUL": "1", "OIN_RC1_PROPAGATE": "1"}, False),
    "cap": ({"OIN_CAP_IGNORES_METAL": "1"}, False),
    "fix3": (
        {"OIN_H_FAITHFUL": "1", "OIN_RC1_PROPAGATE": "1", "OIN_CAP_IGNORES_METAL": "1"},
        False,
    ),
    # E2 lane (results-v0.4.19-e2/): perception-order levers
    "n2": ({"OIN_N_VALENCE_2": "1"}, False),
    "e2b": ({"OIN_N_VALENCE_2": "1", "OIN_CANONICAL_RESONANCE": "1"}, False),
    "e2c": (
        {
            "OIN_N_VALENCE_2": "1",
            "OIN_CANONICAL_RESONANCE": "1",
            "OIN_CANONICAL_CHARGES": "1",
        },
        False,
    ),
}
ISO = ("ISO", "ISO_MARGINAL", "ISO_CLASH")


def _worker(mol: str, arm: str):
    sys.path.insert(0, str(SRC))  # beats PYTHONPATH, which is the point
    sys.path.insert(0, str(HERE.parents[1] / "census"))
    import oinsmiles
    from oinsmiles import XYZToSMILES

    out = {"molecule": mol, "arm": arm, "src": str(Path(oinsmiles.__file__).resolve().parents[1])}
    if ARMS[arm][1]:
        from oinsmiles.utils import oin_aligner

        oin_aligner.OINDiscreteAligner._canonical_ring_signature = classmethod(
            lambda cls, smiles: None
        )
        out["rc1_forced_off"] = True
    xyz = COHORT / f"{mol}.xyz"
    t0 = time.monotonic()
    try:
        oin = XYZToSMILES().convert(str(xyz))
        out["oin"] = oin
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"[:160]
        print(json.dumps(out), flush=True)
        return
    out["encode_s"] = round(time.monotonic() - t0, 2)
    import string_sufficiency as ss

    rec = ss.parseback_one({"molecule": mol, "oin": oin, "xyz": str(xyz)})
    for k in (
        "parser_graph",
        "adapter_graph",
        "adapter_h_decoration",
        "adapter_h_diff_atoms",
        "n_unslotted",
        "adapter_error",
        "parser_error",
    ):
        if k in rec:
            out[k] = rec[k]
    print(json.dumps(out), flush=True)


def _run(mol, arm, timeout):
    env = {k: v for k, v in os.environ.items() if not k.startswith("OIN_")}
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    env.update(ARMS[arm][0])
    try:
        p = subprocess.run(
            [sys.executable, str(HERE), "--one", mol, "--arm", arm],
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"molecule": mol, "arm": arm, "error": f"TIMEOUT@{timeout:.0f}s"}
    lines = [ln for ln in p.stdout.splitlines() if ln.startswith("{")]
    if not lines:
        return {"molecule": mol, "arm": arm, "error": "NO_OUTPUT: " + p.stderr[-200:]}
    return json.loads(lines[-1])


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--one", help=argparse.SUPPRESS)
    ap.add_argument("--arm", choices=list(ARMS), help=argparse.SUPPRESS)
    ap.add_argument("--names", type=Path)
    ap.add_argument("--classes", type=Path, help="json {class: [molecules]}")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--arms", default="shipped,hfaith,rc1off,both")
    args = ap.parse_args()
    if args.one:
        return _worker(args.one, args.arm)
    arms = args.arms.split(",")
    if not set(arms) <= set(ARMS) or "shipped" not in arms:
        sys.exit(f"ABORT: arms must include 'shipped' and be from {list(ARMS)}")
    names = args.names.read_text().split()
    if not names:
        sys.exit("ABORT: no molecules -- refusing an empty denominator")
    cls = {m: k for k, v in json.loads(args.classes.read_text()).items() for m in v}
    jobs = [(m, a) for m in names for a in arms]
    with ThreadPoolExecutor(args.jobs) as pool:
        recs = list(pool.map(lambda j: _run(j[0], j[1], args.timeout), jobs))
    foreign = {r.get("src") for r in recs} - {None, str(SRC)}
    if foreign:
        sys.exit(f"ABORT: a worker imported oinsmiles from {foreign}, not {SRC}")
    with open(args.out, "w") as fh:
        for r in recs:
            fh.write(json.dumps({k: v for k, v in r.items() if k != "src"}) + "\n")
        fh.write(f"#DONE {len(recs)}\n")

    by = defaultdict(dict)
    for r in recs:
        by[r["molecule"]][r["arm"]] = r
    print(f"DENOMINATOR {len(names)} molecules x {len(arms)} arms   src={SRC}")

    def ok(r):
        return r.get("adapter_graph") in ISO and r.get("adapter_h_decoration") == "SAME"

    for k in sorted(set(cls.values())):
        ms = [m for m in names if cls.get(m) == k]
        print(f"\n== {k}  n={len(ms)}")
        for a in arms:
            rs = [by[m].get(a, {}) for m in ms]
            sh = [by[m].get("shipped", {}) for m in ms]
            changed = sum(
                1 for r, s in zip(rs, sh) if r.get("oin") and s.get("oin") and r["oin"] != s["oin"]
            )
            errs = sum(1 for r in rs if r.get("error"))
            graph = Counter(r.get("adapter_graph", "?") for r in rs)
            hdec = Counter(r.get("adapter_h_decoration", "?") for r in rs)
            n_ok = sum(ok(r) for r in rs)
            gain = sum(ok(r) and not ok(s) for r, s in zip(rs, sh))
            loss = sum(ok(s) and not ok(r) for r, s in zip(rs, sh))
            print(
                f"   {a:8s} string changed vs shipped {changed:3d}  errors {errs:2d}  "
                f"ISO&SAME {n_ok:3d} (+{gain}/-{loss} vs shipped)   graph {dict(graph.most_common(4))}"
                f"   H {dict(hdec.most_common(3))}"
            )
    sys.stdout.flush()


if __name__ == "__main__":
    main()
