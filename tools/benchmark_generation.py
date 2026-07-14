#!/usr/bin/env python3
"""Serial, seed-aware, per-stage profiler for the MetalloGen OIN -> 3D generator.

This is P0 of the v0.4.0-perf wave. It produces the measurement ground truth every
other phase's acceptance gate ("N.Nx faster, non-overlapping IQR") cites; see the
v0.4.0 entry in ``CHANGELOG.md`` for the headline numbers it measured.

Why this exists / what it does differently from a naive timer:

* **cProfile cannot see into RDKit.** ``AllChem.EmbedMolecule`` is a Boost.Python
  callable; cProfile folds its wall-time into the *calling* Python function's
  ``tottime`` (on fac-Ir(ppy)3, ``get_embedding`` shows ~48 s of "self time" that is
  actually ETKDG). We attribute cost by *wrapping the callable and timing it*, split
  by caller source line. ``EmbedMolecule`` returns ``-1`` on failure rather than
  raising -- we count that.
* **Generation is deterministic (S6, seed=42).** The same OIN yields byte-identical
  XYZ across runs, so remaining wall-clock variance is machine jitter. We still report
  median + IQR over N runs and mark a stage ``UNSTABLE`` when its run-to-run
  IQR > 20% of its median.
* **A hang is a finding, not a crash.** Each molecule's runs execute in a spawned
  subprocess under a hard wall-clock deadline; ``signal.alarm`` cannot interrupt a
  native RDKit/PuLP hang, so we ``SIGKILL`` on expiry (pattern lifted from
  ``tools/test_dataset_roundtrip.py::_supervise`` -- copied, not imported, because P6
  owns that file).
* **Always FF-only.** We construct ``OIN3DGenerator(engine="metallogen",
  optimizer=None)``. The default ``optimizer="xtb"`` costs a ``shutil.which`` +
  ``deepcopy`` per conformer and refines nothing when no xtb binary is present.

Heavy imports (rdkit / oinsmiles) are deferred into the subprocess and the encode
helpers so this module imports cheaply for unit testing.

Usage::

    uv run python tools/benchmark_generation.py --goldens --runs 5 --json goldens.json
    uv run python tools/benchmark_generation.py --sample 60 --stratify-cn \\
        --dataset-dir /path/tmCAT-tmPHOTO_xyz_dataset --runs 5 --json sample.json
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
import sys
import time

# --------------------------------------------------------------------------------------
# The four golden OINs -- the perf-wave A/B set (cisplatin / ferrocene / fac-Ir(ppy)3 /
# PdCl2-BINAP)
# --------------------------------------------------------------------------------------
GOLDENS: dict[str, str] = {
    "cisplatin": "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
    "ferrocene": (
        "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
    ),
    "fac-Ir(ppy)3": (
        "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{1}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{4}1"
    ),
    "PdCl2-BINAP": (
        "[Pd_SPL].c1ccc(P{0}(c2ccccc2)c2ccc3ccccc3c2-c2c(P{1}(c3ccccc3)c3ccccc3)"
        "ccc3ccccc23)cc1.[Cl]{2}.[Cl]{3}"
    ),
}

# OIN geo-code -> coordination number. Mirrors
# oinsmiles.generation.metallogen_adapter.OIN_TO_METALLOGEN_GEO /
# _expected_coordination_number (the CN is the integer prefix of the MetalloGen geo
# name, e.g. SPL -> "4_square_planar" -> 4). Kept local so unit tests need no heavy
# imports; keep in sync if the source table changes.
GEO_CN: dict[str, int] = {
    "LIN": 2,
    "TPL": 3,
    "SQP": 4,
    "SPL": 4,
    "TET": 4,
    "TPY": 4,
    "SPY": 5,
    "TBP": 5,
    "OCT": 6,
    "PBP": 7,
    "SQA": 8,
}

# Stage under whose wall an internal EmbedMolecule is *also* counted -- documented
# overlap, never summed into a synthetic total.
UNSTABLE_FRACTION = 0.20
SUBMS_SECONDS = 1e-3  # stages with median below this are reported count-first


# ======================================================================================
# Per-stage recorder + monkeypatch wrappers  (installed inside the subprocess only)
# ======================================================================================
class StageRecorder:
    """Accumulates ``{stage: {wall, count, fail}}`` for one ``generate()`` call."""

    def __init__(self) -> None:
        self.stages: dict[str, dict[str, float]] = {}
        self.depth = 0  # reentrancy counter for _timed_outermost

    def reset(self) -> None:
        self.stages = {}
        self.depth = 0

    def add(self, key: str, dt: float, fail: int = 0) -> None:
        s = self.stages.get(key)
        if s is None:
            s = self.stages[key] = {"wall": 0.0, "count": 0, "fail": 0}
        s["wall"] += dt
        s["count"] += 1
        s["fail"] += fail

    def snapshot(self) -> dict[str, dict[str, float]]:
        return {k: dict(v) for k, v in self.stages.items()}


REC = StageRecorder()

# Saved originals (populated by install_wrappers, restored by uninstall_wrappers).
_ORIG: dict[str, object] = {}


def embed_failed(r) -> int:
    """AllChem.EmbedMolecule signals failure by RETURNING -1 (it does not raise), so a
    counter that only catches exceptions reports 100% success on a half-failing embed."""
    return int(r == -1)


def _timed_outermost(rec: StageRecorder, key: str, fn, args, kwargs, fail_pred):
    """Time ``fn(*args)`` under ``key``, reentrancy-safe: a nested (recursive) call
    while this stage is already on the stack runs through untimed and uncounted, so a
    self-recursing target (e.g. get_valid_molecule -> get_valid_molecule) is attributed
    exactly once. Returns fn's result."""
    if rec.depth > 0:
        return fn(*args, **kwargs)
    rec.depth += 1
    t = time.perf_counter()
    r = None
    try:
        r = fn(*args, **kwargs)
        return r
    finally:
        rec.depth -= 1
        rec.add(key, time.perf_counter() - t, fail=int(bool(fail_pred(r))))


def _install_wrappers() -> None:
    """Patch the eight stage chokepoints. All targets verified call-by-reference so a
    module/class attribute patch intercepts the pipeline's own calls."""
    import pulp.apis.coin_api as coin_api  # noqa: PLC0415
    from rdkit.Chem import AllChem  # noqa: PLC0415

    from oinsmiles.generation import metallogen_adapter as MA  # noqa: PLC0415
    from oinsmiles.generator3d import chem as GC  # noqa: PLC0415
    from oinsmiles.generator3d import clean_geometry as CG  # noqa: PLC0415
    from oinsmiles.utils import oin_aligner as OA  # noqa: PLC0415

    perf = time.perf_counter

    # 1) embed.ETKDG -- wrap the shared AllChem.EmbedMolecule; attribute by caller
    #    source line (embed.py:547 primary / :553 metal-swap / :583 rebuild retry, and
    #    clean_geometry.py:239). EmbedMolecule returns -1 on failure (does not raise).
    _ORIG["embed"] = AllChem.EmbedMolecule

    def embed_spy(*a, **k):
        fr = sys._getframe(1)
        loc = f"{os.path.basename(fr.f_code.co_filename)}:{fr.f_lineno}"
        t = perf()
        r = _ORIG["embed"](*a, **k)
        dt = perf() - t
        fail = embed_failed(r)
        REC.add(f"embed.ETKDG[{loc}]", dt, fail)
        REC.add("embed.ETKDG", dt, fail)
        return r

    AllChem.EmbedMolecule = embed_spy

    # 2) perception.PuLP -- Molecule.get_valid_molecule; recurses (chem.py:2298/2339),
    #    so time only the outermost (depth-0) call to avoid double counting.
    _ORIG["gvm"] = GC.Molecule.get_valid_molecule

    def gvm_spy(self, *a, **k):
        return _timed_outermost(
            REC, "perception.PuLP", _ORIG["gvm"], (self, *a), k, lambda r: r is None
        )

    GC.Molecule.get_valid_molecule = gvm_spy

    # 3) perception.CBC -- count CBC subprocess spawns (PULP_CBC_CMD inherits solve_CBC
    #    unoverridden from the COIN_CMD base).
    _ORIG["cbc"] = coin_api.COIN_CMD.solve_CBC

    def cbc_spy(self, *a, **k):
        t = perf()
        r = _ORIG["cbc"](self, *a, **k)
        REC.add("perception.CBC", perf() - t)
        return r

    coin_api.COIN_CMD.solve_CBC = cbc_spy

    # 4) ff.clean_geometry -- TMCOptimizer.clean_geometry (bool return). Its wall
    #    INCLUDES the clean_geometry.py:239 embed and its internal PuLP calls;
    #    overlap is expected and documented, never summed into total.
    _ORIG["clean"] = CG.TMCOptimizer.clean_geometry

    def clean_spy(self, *a, **k):
        t = perf()
        r = _ORIG["clean"](self, *a, **k)
        REC.add("ff.clean_geometry", perf() - t, fail=int(not r))
        return r

    CG.TMCOptimizer.clean_geometry = clean_spy

    # 5) select.geometry -- the stage function, plus the _map_to_template chokepoint
    #    (its call count exposes P5's double-match).
    _ORIG["select"] = MA._select_by_geometry

    def select_spy(*a, **k):
        t = perf()
        r = _ORIG["select"](*a, **k)
        REC.add("select.geometry", perf() - t)
        return r

    MA._select_by_geometry = select_spy

    _ORIG["map"] = OA.OINDiscreteAligner._map_to_template

    def map_spy(self, *a, **k):
        t = perf()
        r = _ORIG["map"](self, *a, **k)
        REC.add("select.map_to_template", perf() - t)
        return r

    OA.OINDiscreteAligner._map_to_template = map_spy

    # 6) contract_mol -- module-level build_contract_mol (Mol | None).
    _ORIG["contract"] = MA.build_contract_mol

    def contract_spy(*a, **k):
        t = perf()
        r = _ORIG["contract"](*a, **k)
        REC.add("contract_mol", perf() - t, fail=int(r is None))
        return r

    MA.build_contract_mol = contract_spy

    # 7) reencode -- module-level _reencode_oin_fast (str | None).
    _ORIG["reencode"] = MA._reencode_oin_fast

    def reencode_spy(*a, **k):
        t = perf()
        r = _ORIG["reencode"](*a, **k)
        REC.add("reencode", perf() - t, fail=int(r is None))
        return r

    MA._reencode_oin_fast = reencode_spy


def _uninstall_wrappers() -> None:
    """Restore originals (used by tests; the subprocess just exits)."""
    if not _ORIG:
        return
    import pulp.apis.coin_api as coin_api  # noqa: PLC0415
    from rdkit.Chem import AllChem  # noqa: PLC0415

    from oinsmiles.generation import metallogen_adapter as MA  # noqa: PLC0415
    from oinsmiles.generator3d import chem as GC  # noqa: PLC0415
    from oinsmiles.generator3d import clean_geometry as CG  # noqa: PLC0415
    from oinsmiles.utils import oin_aligner as OA  # noqa: PLC0415

    AllChem.EmbedMolecule = _ORIG["embed"]
    GC.Molecule.get_valid_molecule = _ORIG["gvm"]
    coin_api.COIN_CMD.solve_CBC = _ORIG["cbc"]
    CG.TMCOptimizer.clean_geometry = _ORIG["clean"]
    MA._select_by_geometry = _ORIG["select"]
    OA.OINDiscreteAligner._map_to_template = _ORIG["map"]
    MA.build_contract_mol = _ORIG["contract"]
    MA._reencode_oin_fast = _ORIG["reencode"]
    _ORIG.clear()


# ======================================================================================
# Subprocess worker: encode-once (sample) then N timed generations, streamed back
# ======================================================================================
# Message tags on the result queue.
MSG_ENCODED = "encoded"
MSG_META = "meta"
MSG_RUN = "run"
MSG_ERROR = "error"
MSG_DONE = "done"

_EMPTY = object()  # sentinel returned by the get_msg callable when the queue is empty


def _collect_runs(get_msg, is_alive, runs, mol_timeout, monotonic, initial_oin):
    """Drain streamed child messages into up to ``runs`` records under a per-run
    wall-clock deadline (reset on every message). Pure -- callables are injected so the
    watchdog is testable without a live subprocess.

    ``get_msg()`` returns a message tuple or ``_EMPTY``; ``is_alive()`` -> bool;
    ``monotonic()`` -> float clock.
    """
    records: list[dict] = []
    oin = initial_oin
    timed_out = False
    error = None
    deadline = monotonic() + mol_timeout
    while len(records) < runs:
        msg = get_msg()
        if msg is _EMPTY:
            if monotonic() > deadline:
                timed_out = True
                break
            if not is_alive():
                error = error or "child died with no result"
                break
            continue
        tag = msg[0]
        if tag == MSG_ENCODED:
            oin = msg[1]
            deadline = monotonic() + mol_timeout
        elif tag == MSG_RUN:
            records.append(msg[1])
            deadline = monotonic() + mol_timeout
        elif tag == MSG_ERROR:
            error = f"{msg[1]['type']}: {msg[1]['msg']}"
            break
        elif tag == MSG_DONE:
            break
    return {"records": records, "oin": oin, "timed_out": timed_out, "error": error}


def _child_entry(queue, payload: str, is_oin: bool, runs: int, seed: int) -> None:
    """Run inside a spawned subprocess. Installs wrappers, encodes (if given an XYZ
    path), then times ``runs`` generations, streaming one message per run so a later
    hang does not lose earlier results."""
    import contextlib  # noqa: PLC0415
    import hashlib  # noqa: PLC0415
    import random  # noqa: PLC0415

    from rdkit import RDLogger  # noqa: PLC0415

    RDLogger.DisableLog("rdApp.*")
    random.seed(seed)  # partial guarantee until P1 threads --seed into the generator
    devnull = open(os.devnull, "w")  # noqa: SIM115

    try:
        _install_wrappers()
        from oinsmiles.generation.engine import OIN3DGenerator  # noqa: PLC0415

        if is_oin:
            oin = payload
        else:
            from oinsmiles import XYZToSMILES  # noqa: PLC0415

            with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
                oin = XYZToSMILES().convert(payload)
            queue.put((MSG_ENCODED, oin))

        for _ in range(runs):
            REC.reset()
            gen = OIN3DGenerator(engine="metallogen", optimizer=None)
            t0 = time.perf_counter()
            with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
                result = gen.generate(oin)
            total = time.perf_counter() - t0
            xyz = result.xyz or ""
            atoms = max(0, len(xyz.strip().splitlines()) - 2)
            queue.put(
                (
                    MSG_RUN,
                    {
                        "total": total,
                        "stages": REC.snapshot(),
                        "atoms": atoms,
                        "xyz_sha": hashlib.sha256(xyz.encode()).hexdigest()[:16],
                    },
                )
            )
        queue.put((MSG_DONE,))
    except Exception as exc:  # a real generation failure for this molecule
        import traceback as _tb  # noqa: PLC0415

        queue.put(
            (MSG_ERROR, {"type": type(exc).__name__, "msg": str(exc), "tb": _tb.format_exc()})
        )


def _supervise(payload: str, is_oin: bool, runs: int, seed: int, mol_timeout: float) -> dict:
    """Spawn ``_child_entry`` and collect ``runs`` streamed run-records under a per-run
    wall-clock deadline. SIGKILL on expiry -- a timeout is recorded, not raised.

    Returns a dict: {records, oin, timed_out, error}. ``records`` may be shorter than
    ``runs`` if the child hung or died partway.
    """
    import multiprocessing  # noqa: PLC0415
    import queue as queue_mod  # noqa: PLC0415

    ctx = multiprocessing.get_context("spawn")
    q = ctx.Queue()
    proc = ctx.Process(target=_child_entry, args=(q, payload, is_oin, runs, seed))
    proc.start()

    def get_msg():
        try:
            return q.get(timeout=0.5)
        except queue_mod.Empty:
            return _EMPTY

    try:
        out = _collect_runs(
            get_msg, proc.is_alive, runs, mol_timeout, time.monotonic, payload if is_oin else None
        )
    finally:
        if proc.is_alive():
            proc.kill()
        proc.join(timeout=30)
        if proc.is_alive():  # pragma: no cover -- kill already sent
            proc.kill()
            proc.join()
    return out


# ======================================================================================
# Pure helpers (unit-tested)
# ======================================================================================
def geo_to_cn(metal_geo: str | None) -> int | None:
    """Coordination number from a registry ``metal_geo`` token (e.g. ``"Pd_SPL"``)."""
    if not metal_geo or "_" not in metal_geo:
        return GEO_CN.get(metal_geo) if metal_geo else None
    return GEO_CN.get(metal_geo.rsplit("_", 1)[-1])


def summarize(values: list[float]) -> dict:
    """median, p25, p75, IQR, min, max, iqr_pct, unstable for a list of wall-times."""
    vals = [float(v) for v in values]
    n = len(vals)
    med = statistics.median(vals) if n else 0.0
    if n >= 2:
        q = statistics.quantiles(vals, n=4, method="inclusive")  # [p25, p50, p75]
        p25, p75 = q[0], q[2]
    else:
        p25 = p75 = med
    iqr = p75 - p25
    iqr_pct = (iqr / med * 100.0) if med > 0 else 0.0
    # A stage with a truly negligible median is reported count-first, not flagged.
    unstable = med >= SUBMS_SECONDS and iqr > UNSTABLE_FRACTION * med
    return {
        "median": med,
        "p25": p25,
        "p75": p75,
        "iqr": iqr,
        "iqr_pct": iqr_pct,
        "min": min(vals) if n else 0.0,
        "max": max(vals) if n else 0.0,
        "n": n,
        "unstable": unstable,
        "subms": 0 < med < SUBMS_SECONDS,
    }


def aggregate_runs(records: list[dict]) -> dict:
    """Fold N per-run records into total + per-stage summaries. Counts/fails are taken
    from run 0 (deterministic pipeline) and flagged if they drift across runs."""
    if not records:
        return {"runs": 0}
    totals = [r["total"] for r in records]
    stage_keys: set[str] = set()
    for r in records:
        stage_keys.update(r["stages"].keys())

    stages = {}
    for key in sorted(stage_keys):
        walls = [r["stages"].get(key, {}).get("wall", 0.0) for r in records]
        counts = [r["stages"].get(key, {}).get("count", 0) for r in records]
        fails = [r["stages"].get(key, {}).get("fail", 0) for r in records]
        s = summarize(walls)
        s["count"] = counts[0]
        s["fail"] = fails[0]
        s["count_stable"] = len(set(counts)) == 1
        s["fail_stable"] = len(set(fails)) == 1
        stages[key] = s

    shas = {r.get("xyz_sha") for r in records}
    return {
        "runs": len(records),
        "atoms": records[0].get("atoms"),
        "total": summarize(totals),
        "stages": stages,
        "byte_identical": len(shas) == 1,
    }


def find_conflicts(pgrep_output: dict[str, str], my_pids: set[int]) -> list[str]:
    """Given {pattern: pgrep_stdout}, return conflict descriptions (excluding my_pids).

    Factored out so the serial guard is testable without a live process table.
    """
    conflicts = []
    for pattern, out in pgrep_output.items():
        pids = [int(p) for p in out.split() if p.strip().isdigit() and int(p) not in my_pids]
        if pids:
            conflicts.append(f"{pattern} (pids: {', '.join(map(str, pids))})")
    return conflicts


def stratified_sample(
    registry: list[dict],
    path_index: dict[str, str],
    n: int,
    cn_lo: int,
    cn_hi: int,
    seed: int,
) -> tuple[list[dict], dict[int, int]]:
    """CN-stratified draw of ``n`` molecules from the registry.

    ``registry`` = list of {molecule, metal_geo, ...}; ``path_index`` maps molecule
    name -> xyz path (only molecules present here are eligible). Allocation across CN
    buckets follows the observed distribution (largest-remainder), with >=1 per
    non-empty in-range bucket so small CNs stay represented. Deterministic under seed.

    Returns (selected, observed_distribution).
    """
    import random  # noqa: PLC0415

    buckets: dict[int, list[dict]] = {}
    for rec in registry:
        cn = geo_to_cn(rec.get("metal_geo"))
        if cn is None or not (cn_lo <= cn <= cn_hi):
            continue
        name = rec.get("molecule")
        if name not in path_index:
            continue
        buckets.setdefault(cn, []).append({"molecule": name, "cn": cn, "path": path_index[name]})

    observed = {cn: len(v) for cn, v in sorted(buckets.items())}
    total_avail = sum(observed.values())
    if total_avail == 0:
        return [], observed
    n = min(n, total_avail)

    # Largest-remainder allocation proportional to bucket size, floor 1 per bucket.
    raw = {cn: (cnt / total_avail) * n for cn, cnt in observed.items()}
    alloc = {cn: min(observed[cn], max(1, int(raw[cn]))) for cn in observed}
    # Reconcile to exactly n, capped by availability.
    while sum(alloc.values()) > n:
        cn = max(alloc, key=lambda c: (alloc[c] - raw[c], c))
        if alloc[cn] > 1:
            alloc[cn] -= 1
        else:
            break
    while sum(alloc.values()) < n:
        cand = [c for c in observed if alloc[c] < observed[c]]
        if not cand:
            break
        cn = max(cand, key=lambda c: (raw[c] - alloc[c], c))
        alloc[cn] += 1

    rng = random.Random(seed)
    selected: list[dict] = []
    for cn in sorted(buckets):
        pool = sorted(buckets[cn], key=lambda d: d["molecule"])
        rng.shuffle(pool)
        selected.extend(pool[: alloc.get(cn, 0)])
    selected.sort(key=lambda d: (d["cn"], d["molecule"]))
    return selected, observed


def build_path_index(dataset_dir: str) -> dict[str, str]:
    """molecule name (basename w/o .xyz, excluding *_generated) -> absolute xyz path."""
    index: dict[str, str] = {}
    for sub in ("cat", "photo"):
        base = os.path.join(dataset_dir, sub)
        if not os.path.isdir(base):
            continue
        for fn in os.listdir(base):
            if fn.endswith(".xyz") and not fn.endswith("_generated.xyz"):
                index.setdefault(fn[:-4], os.path.join(base, fn))
    return index


def load_registry(dataset_dir: str) -> list[dict]:
    reg_path = os.path.join(dataset_dir, "20260707-results", "case_registry.json")
    with open(reg_path) as f:
        return json.load(f)


# ======================================================================================
# Environment / metadata
# ======================================================================================
def _git_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.abspath(__file__)),
            check=False,
        )
        return out.stdout.strip() or "unknown"
    except Exception:  # pragma: no cover
        return "unknown"


def _rdkit_version() -> str:
    try:
        import rdkit  # noqa: PLC0415

        return rdkit.__version__
    except Exception:  # pragma: no cover
        return "unknown"


def _xtb_status() -> str:
    import shutil  # noqa: PLC0415

    return shutil.which("xtb") or "absent"


def _ancestor_pids() -> set[int]:
    """This process and its full parent chain (python, its uv/bash launchers, ...).

    ``pgrep -f benchmark_generation`` matches every process whose command line contains
    the script path -- which includes the ``uv run`` and shell wrappers that launched us.
    Excluding only pid+ppid would misread those launchers as a competing run, so we
    exclude the whole ancestry.
    """
    pids: set[int] = set()
    pid = os.getpid()
    for _ in range(64):
        pids.add(pid)
        try:
            with open(f"/proc/{pid}/stat") as f:
                data = f.read()
            ppid = int(data[data.rindex(")") + 2 :].split()[1])
        except Exception:
            break
        if ppid <= 0 or ppid in pids:
            break
        pid = ppid
    return pids


def _assert_serial() -> None:
    """Refuse to run alongside a dataset sweep or another benchmark (contention
    fabricates no_conformers failures and corrupts every timing number)."""
    my_pids = _ancestor_pids()
    patterns = ("test_dataset_roundtrip", "benchmark_generation")
    outputs = {}
    for pat in patterns:
        r = subprocess.run(["pgrep", "-f", pat], capture_output=True, text=True, check=False)
        outputs[pat] = r.stdout
    conflicts = find_conflicts(outputs, my_pids)
    if conflicts:
        sys.exit(
            "REFUSING TO RUN -- a competing process is active; contention corrupts "
            "timings:\n  " + "\n  ".join(conflicts)
        )


# ======================================================================================
# Reporting
# ======================================================================================
def _fmt_stage_line(key: str, s: dict) -> str:
    flag = " UNSTABLE" if s.get("unstable") else ""
    cnt = f" x{s['count']}" if "count" in s else ""
    fail = f" fail={s['fail']}" if s.get("fail") else ""
    if s.get("subms"):
        return f"    {key:34s} {s['median'] * 1000:7.2f}ms (sub-ms){cnt}{fail}"
    return (
        f"    {key:34s} {s['median']:7.3f}s  IQR {s['iqr']:6.3f}s "
        f"({s['iqr_pct']:4.1f}%){cnt}{fail}{flag}"
    )


def _print_molecule(name: str, agg: dict, timed_out: bool, error: str | None) -> None:
    if timed_out:
        print(
            f"\n{name}: TIMED OUT (recorded as finding; {len(agg.get('stages', {}))} partial stages)"
        )
    if error and not agg.get("runs"):
        print(f"\n{name}: FAILED -- {error}")
        return
    if not agg.get("runs"):
        print(f"\n{name}: no runs completed")
        return
    t = agg["total"]
    bi = "byte-identical" if agg["byte_identical"] else "XYZ VARIES"
    tflag = " UNSTABLE" if t["unstable"] else ""
    print(
        f"\n{name}  atoms={agg['atoms']}  runs={agg['runs']}  {bi}\n"
        f"    {'total':34s} {t['median']:7.3f}s  IQR {t['iqr']:6.3f}s "
        f"({t['iqr_pct']:4.1f}%)  min {t['min']:.2f} max {t['max']:.2f}{tflag}"
    )
    for key, s in agg["stages"].items():
        print(_fmt_stage_line(key, s))


def _benchmark_one(name, payload, is_oin, runs, seed, mol_timeout) -> dict:
    print(f"  benchmarking {name} ...", flush=True)
    sup = _supervise(payload, is_oin, runs, seed, mol_timeout)
    agg = aggregate_runs(sup["records"])
    entry = {
        "name": name,
        "oin": sup["oin"],
        "timed_out": sup["timed_out"],
        "error": sup["error"],
        **agg,
    }
    _print_molecule(name, agg, sup["timed_out"], sup["error"])
    return entry


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="benchmark_generation",
        description="Serial per-stage profiler for the MetalloGen OIN->3D generator.",
    )
    p.add_argument("--goldens", action="store_true", help="Benchmark the four golden OINs.")
    p.add_argument("--sample", type=int, default=0, help="Benchmark N dataset molecules.")
    p.add_argument(
        "--stratify-cn", action="store_true", help="CN-stratify the --sample draw (CN 2-7)."
    )
    p.add_argument("--runs", type=int, default=5, help="Runs per molecule (default 5).")
    p.add_argument(
        "--seed", type=int, default=42, help="Seed Python random (partial until P1 plumbs it)."
    )
    p.add_argument("--sample-seed", type=int, default=1234, help="Deterministic sampling seed.")
    p.add_argument("--mol-timeout", type=float, default=240.0, help="Per-run wall-clock cap (s).")
    p.add_argument(
        "--dataset-dir",
        default="/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset",
        help="Dataset root (main-checkout only; gitignored).",
    )
    p.add_argument("--cn-lo", type=int, default=2)
    p.add_argument("--cn-hi", type=int, default=7)
    p.add_argument("--json", dest="json_out", default=None, help="Write machine-readable output.")
    p.add_argument(
        "--no-serial-guard", action="store_true", help="Skip the pgrep serial guard (tests only)."
    )
    args = p.parse_args(argv)

    if not args.goldens and args.sample <= 0:
        p.error("nothing to do: pass --goldens and/or --sample N")
    if not args.no_serial_guard:
        _assert_serial()

    meta = {
        "rdkit": _rdkit_version(),
        "commit": _git_commit(),
        "runs": args.runs,
        "seed": args.seed,
        "optimizer": None,
        "engine": "metallogen",
        "host_cpu_count": os.cpu_count(),
        "xtb": _xtb_status(),
        "mol_timeout": args.mol_timeout,
    }
    print("== benchmark_generation ==")
    print(
        f"   rdkit={meta['rdkit']} commit={meta['commit']} cores={meta['host_cpu_count']} "
        f"runs={args.runs} xtb={meta['xtb']}"
    )

    results = {"meta": meta, "goldens": [], "sample": [], "sample_stratification": {}}

    if args.goldens:
        print("\n-- Four golden OINs --")
        for name, oin in GOLDENS.items():
            results["goldens"].append(
                _benchmark_one(name, oin, True, args.runs, args.seed, args.mol_timeout)
            )

    if args.sample > 0:
        print(f"\n-- CN-stratified sample (target {args.sample}) --")
        registry = load_registry(args.dataset_dir)
        path_index = build_path_index(args.dataset_dir)
        selected, observed = stratified_sample(
            registry, path_index, args.sample, args.cn_lo, args.cn_hi, args.sample_seed
        )
        chosen_dist: dict[int, int] = {}
        for d in selected:
            chosen_dist[d["cn"]] = chosen_dist.get(d["cn"], 0) + 1
        results["sample_stratification"] = {
            "available_by_cn": observed,
            "chosen_by_cn": chosen_dist,
            "selected": [d["molecule"] for d in selected],
        }
        print(
            f"   available_by_cn={observed}  chosen_by_cn={chosen_dist}  ({len(selected)} molecules)"
        )
        for i, d in enumerate(selected):
            tag = f"{d['molecule']}[CN{d['cn']}] ({i + 1}/{len(selected)})"
            results["sample"].append(
                _benchmark_one(tag, d["path"], False, args.runs, args.seed, args.mol_timeout)
            )

    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump(results, f, indent=2, default=str)
        print(f"\nWrote {args.json_out}")

    # Summary of instability / timeouts for quick scanning.
    unstable_total = [
        e["name"]
        for grp in ("goldens", "sample")
        for e in results[grp]
        if e.get("total", {}).get("unstable")
    ]
    timed_out = [
        e["name"] for grp in ("goldens", "sample") for e in results[grp] if e.get("timed_out")
    ]
    print(f"\n== done ==  unstable-total={len(unstable_total)}  timed-out={len(timed_out)}")
    if timed_out:
        print("   timed out:", ", ".join(timed_out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
