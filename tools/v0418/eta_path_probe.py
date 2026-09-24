"""L2 DETACHED, step 1: does the OTHER construction path attach the ring? Two arms, fresh processes.

``eta_distance_audit.py`` measured the defect: on the default path the generator puts eta carbons a
median +0.73 A too far in the DETACHED class (+0.17 A even in verified passes), and the generated
M-C(eta) distance barely depends on the metal. The generator has two ways to build an eta complex:

    dg      the DEFAULT -- OIN-direct assembly feeding the dummy-atom distance-geometry embed: a
            dummy atom bonded to every ring atom is pinned at (r_M + r_dummy) x 0.4..0.7
    rigid   ``ff_params={"oin_direct": True}`` -- the ligand is built alone and its eta face is
            PLACED so that M-atom = (r_M + r_atom) x scale  (``embed._place_haptic``). Default OFF:
            it "regressed eta" in v0.4.4, judged by a metric that could not see the structure.

This runs both on a sample and reports, per arm: did it build, how far are the eta atoms, is the
coordination intact (``oin.coordination``), and does the HONEST round trip hold -- an independent
``XYZToSMILES().convert`` of the generated coordinates against the input string, never the
generator's own bond graph.

⚠ A SAMPLE, not a verdict: it says which path to take seriously. The verdict is a real A/B through
the harness with the neutral ruler (VERIFIED), over the whole class and its passing controls.

    <main>/.venv/bin/python tools/v0418/eta_path_probe.py --names <file> --out <jsonl> [--jobs 6]
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
#: arm -> (ff_params, environment). Levers are WRITTEN, never deleted: unset means the default.
ARMS = {
    "dg": (None, {}),
    "rigid": ({"oin_direct": True}, {}),
    "cov": (None, {"OIN_ETA_COVALENT_TARGET": "1"}),
    "cov_unscaled": (None, {"OIN_ETA_COVALENT_TARGET": "1", "OIN_ETA_TARGET_UNSCALED": "1"}),
    "exempt": (None, {"OIN_VDW_EXEMPT_BINDING": "1"}),
    "cov_exempt": (None, {"OIN_ETA_COVALENT_TARGET": "1", "OIN_VDW_EXEMPT_BINDING": "1"}),
    "all3": (
        None,
        {
            "OIN_ETA_COVALENT_TARGET": "1",
            "OIN_ETA_TARGET_UNSCALED": "1",
            "OIN_VDW_EXEMPT_BINDING": "1",
        },
    ),
}


def _worker(mol: str, arm: str, timeout: float):
    sys.path.insert(0, str(SRC))  # beats PYTHONPATH, which is the point
    sys.path.insert(0, str(HERE.parent))
    import tempfile

    import eta_distance_audit as audit

    import oinsmiles
    from oinsmiles import XYZToSMILES
    from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen
    from oinsmiles.oin.coordination import coordination_report

    out = {"molecule": mol, "arm": arm, "src": str(Path(oinsmiles.__file__).resolve().parents[1])}
    in_text = (COHORT / f"{mol}.xyz").read_text()
    t0 = time.monotonic()
    try:
        oin1 = XYZToSMILES().convert(str(COHORT / f"{mol}.xyz"))
        gen = OIN3DGeneratorMetallogen(
            optimizer=None, ensemble_size=1, timeout=timeout, ff_params=ARMS[arm][0]
        )
        res = gen.generate(oin1)
        out["elapsed_s"] = round(time.monotonic() - t0, 1)
        if res is None or not getattr(res, "xyz", None):
            out["status"] = "NO_STRUCTURE"
            print(json.dumps(out), flush=True)
            return
        out["status"] = "BUILT"
        out["xyz"] = res.xyz
        cmp_ = audit.compare(in_text, res.xyz)
        out["groups"] = [g for g in cmp_.get("groups", []) if g["kind"] == "eta"]
        out["sigma_inflation"] = [
            round(g["inflation"], 3) for g in cmp_.get("groups", []) if g["kind"] == "sigma"
        ]
        out["intact"] = bool(coordination_report(in_text, res.xyz).get("intact"))
        with tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False) as fh:
            fh.write(res.xyz)
        try:
            oin2 = XYZToSMILES().convert(fh.name)
        finally:
            os.unlink(fh.name)
        out["honest_byte_exact"] = oin2 == oin1
    except Exception as e:  # a generator error is a result
        out["status"] = f"ERROR:{type(e).__name__}"
        out["error"] = str(e)[:160]
        out["elapsed_s"] = round(time.monotonic() - t0, 1)
    print(json.dumps(out), flush=True)


def _run(mol, arm, timeout):
    env = {k: v for k, v in os.environ.items() if not k.startswith("OIN_")}
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    env.update(ARMS[arm][1])
    try:
        p = subprocess.run(
            [sys.executable, str(HERE), "--one", mol, "--arm", arm, "--timeout", str(timeout)],
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout * 1.5,
        )
    except subprocess.TimeoutExpired:
        return {"molecule": mol, "arm": arm, "status": "KILLED", "elapsed_s": timeout * 1.5}
    lines = [ln for ln in p.stdout.splitlines() if ln.startswith("{")]
    if not lines:
        return {"molecule": mol, "arm": arm, "status": "NO_OUTPUT", "error": p.stderr[-200:]}
    return json.loads(lines[-1])


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--one", help=argparse.SUPPRESS)
    ap.add_argument("--arm", choices=list(ARMS), help=argparse.SUPPRESS)
    ap.add_argument("--names", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--arms", default="dg,rigid", help="comma-separated, from: " + ",".join(ARMS))
    args = ap.parse_args()
    if args.one:
        return _worker(args.one, args.arm, args.timeout)
    arms = args.arms.split(",")
    if not set(arms) <= set(ARMS):
        sys.exit(f"ABORT: unknown arm in {arms}")

    names = args.names.read_text().split()
    if not names:
        sys.exit("ABORT: no molecules -- refusing an empty denominator")
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
    print(f"{'molecule':16s} {'metal':5s} | " + " | ".join(f"{a:>34s}" for a in arms))
    for m in names:
        cells, metal = [], ""
        for a in arms:
            r = by[m].get(a, {})
            gs = r.get("groups") or []
            if gs:
                metal = gs[0]["metal"]
                k = sum(g["k"] for g in gs)
                infl = sum(g["inflation"] * g["k"] for g in gs) / k
                inb = sum(g["in_band"] for g in gs)
                cells.append(
                    f"{infl:+.2f}A {inb:2d}/{k:<2d} "
                    f"{'INTACT' if r.get('intact') else 'detach'} "
                    f"{'RT' if r.get('honest_byte_exact') else '--'} {r.get('elapsed_s', 0):5.0f}s"
                )
            else:
                cells.append(f"{r.get('status', '?')[:22]:22s} {r.get('elapsed_s', 0):5.0f}s")
        print(f"{m:16s} {metal:5s} | " + " | ".join(f"{c:>34s}" for c in cells))
    for a in arms:
        rs = [by[m].get(a, {}) for m in names]
        c = Counter("built" if r.get("status") == "BUILT" else "no structure" for r in rs)
        print(
            f"{a:6s} built {c['built']}/{len(rs)}   intact {sum(bool(r.get('intact')) for r in rs)}"
            f"   honest byte-exact {sum(bool(r.get('honest_byte_exact')) for r in rs)}"
            f"   sum elapsed {sum(r.get('elapsed_s', 0) for r in rs):.0f}s"
        )
    sys.stdout.flush()


if __name__ == "__main__":
    main()
