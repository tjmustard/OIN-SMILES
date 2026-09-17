"""I1 of the census: judge every GENERATED structure against its INPUT, with no encoder.

    python tools/census/g_vs_input.py --control identity   # must be 100% ISO + SAME/NO_STEREO
    python tools/census/g_vs_input.py --control mirror     # every ruler-chiral input -> MIRROR
    python tools/census/g_vs_input.py --control squeeze    # ligands pushed 0.35 A inward -> still ISO
    python tools/census/g_vs_input.py --control detach     # must be 0% ISO
    python tools/census/g_vs_input.py --control delete     # must be 0% ISO
    python tools/census/g_vs_input.py                      # the real run
    python tools/census/g_vs_input.py --crosstab           # verdict x honest round-trip bucket

WHAT A BROKEN VERSION WOULD PRINT
---------------------------------
A ruler that ignored coordinates would pass ``identity`` and also pass ``mirror`` as all-SAME;
one that ignored the graph would pass ``identity`` and fail ``detach``/``delete``. The four
controls are chosen so no single failure mode passes all of them. Every control structure is
also randomly ROTATED and its atoms randomly PERMUTED, so the controls exercise the claim that
the verdict does not depend on atom order or on which isomorphism the matcher returns.
The real run is void unless all four controls pass -- run them first, read them first.

Paths default to the MAIN checkout (absolute): the dataset and the sweep results are gitignored
and do not exist in a worktree. A relative ``--dataset`` from a worktree once made an instrument
silently exclude its whole population; hence absolute defaults and a printed denominator.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import signal
import sys
import time
import zlib
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))  # sibling module only; no oinsmiles import
import neutral_graph as ng  # noqa: E402

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
DEF_COHORT = MAIN / "cohort-v0.4.5-5k"
DEF_SWEEP = MAIN / "results-v0.4.14-sweep"
DEF_OUT = MAIN / "results-census"


def _rot(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array(
        [
            [a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
            [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],
            [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d],
        ]
    )


def make_control(mode, syms, xyz, rng):
    """A synthetic 'generated' structure whose correct verdict is known by construction."""
    syms, xyz = list(syms), np.array(xyz, dtype=float)
    if mode == "mirror":
        xyz[:, 2] *= -1.0
    elif mode == "squeeze":
        # mimic the generator's short M-L bonds: every ligand slides 0.35 A toward the metal.
        # Topology is untouched, so the verdict MUST stay ISO -- this is the tolerance control
        # the other four lack (they test exact copies and outright breakage, nothing between).
        g = ng.NeutralGraph(syms, xyz)
        if len(g.metals) != 1:
            return None
        m = g.heavy[g.metals[0]]
        hv = np.array(g.heavy)
        owner = {
            int(i): int(hv[int(np.argmin(np.linalg.norm(xyz[hv] - xyz[i], axis=1)))])
            for i in np.where(g.z == 1)[0]
        }
        for comp in g.ligand_components():
            donors = [g.heavy[a] for a in comp if g.metals[0] in g.nbrs[a]]
            if not donors:
                continue
            atoms = {g.heavy[a] for a in comp}
            idx = sorted(atoms | {h for h, o in owner.items() if o in atoms})
            v = xyz[m] - xyz[donors].mean(axis=0)
            n = np.linalg.norm(v)
            if n > 1e-6:
                xyz[idx] += 0.35 * v / n
    elif mode in ("detach", "delete"):
        g = ng.NeutralGraph(syms, xyz)
        bound = [
            c for c in g.ligand_components() if any(m in g.nbrs[a] for a in c for m in g.metals)
        ]
        if not bound or not g.metals:
            return None
        comp = min(bound, key=len) if mode == "delete" else max(bound, key=len)
        atoms = {g.heavy[a] for a in comp}
        if mode == "delete":
            keep = [i for i in range(len(syms)) if i not in atoms]
            syms, xyz = [syms[i] for i in keep], xyz[keep]
        else:
            # carry the ligand's own hydrogens along, so only the METAL attachment changes
            hv = np.array(g.heavy)
            for i in np.where(g.z == 1)[0]:
                k = int(np.argmin(np.linalg.norm(xyz[hv] - xyz[i], axis=1)))
                if hv[k] in atoms:
                    atoms.add(int(i))
            idx = sorted(atoms)
            m = g.heavy[g.metals[0]]
            v = xyz[idx].mean(axis=0) - xyz[m]
            n = np.linalg.norm(v)
            v = v / n if n > 1e-6 else np.array([1.0, 0, 0])
            xyz[idx] += 3.0 * v
    perm = rng.permutation(len(syms))
    xyz = xyz[perm] @ _rot(rng).T + rng.normal(size=3) * 5.0
    return [syms[i] for i in perm], xyz


def _one(mol, cohort, sweep, control, seed):
    inp = cohort / f"{mol}.xyz"
    if control:
        syms, xyz = ng.parse_xyz(inp)
        if not syms:
            return {"molecule": mol, "graph": "UNREADABLE", "stereo": "NA"}
        gen = make_control(
            control, syms, xyz, np.random.default_rng([seed, zlib.crc32(mol.encode())])
        )
        if gen is None:
            return {"molecule": mol, "graph": "CONTROL_NA", "stereo": "NA"}
        rec = ng.compare(inp, gen=gen)
    else:
        rec = ng.compare(inp, sweep / "structures" / f"{mol}_generated.xyz")
    rec["molecule"] = mol
    return rec


def _worker(mols, part, cohort, sweep, control, seed, timeout):
    # SIGALRM stays at SIG_DFL on purpose: a Python handler cannot interrupt a C++ graph match,
    # but the kernel's default action kills the process whatever it is doing.
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    with open(part, "a") as fh:
        for mol in mols:
            fh.write(f"#START {mol}\n")
            fh.flush()
            signal.alarm(timeout)
            try:
                rec = _one(mol, cohort, sweep, control, seed)
            except Exception as exc:  # a crash is a verdict about the instrument; record it
                rec = {
                    "molecule": mol,
                    "graph": "ERROR",
                    "stereo": "NA",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            signal.alarm(0)
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
        os.fsync(fh.fileno())
    os._exit(0)  # skip interpreter finalizers (BLAS thread pools) on the way out


def _reconcile(part, todo):
    """After a worker exits: drop finished molecules; a #START with no verdict was killed."""
    done, started = set(), None
    if part.exists():
        for ln in part.read_text().splitlines():
            if ln.startswith("#START "):
                started = ln[7:]
            else:
                done.add(json.loads(ln)["molecule"])
    if started is not None and started not in done:
        with open(part, "a") as fh:
            fh.write(json.dumps({"molecule": started, "graph": "TIMEOUT", "stereo": "NA"}) + "\n")
        done.add(started)
    return [m for m in todo if m not in done]


def run(args):
    control = args.control
    gen_have = {
        p.name[: -len("_generated.xyz")]
        for p in (args.sweep / "structures").glob("*_generated.xyz")
    }
    cohort = sorted(p.stem for p in args.cohort.glob("*.xyz"))
    dangling = [m for m in cohort if not (args.cohort / f"{m}.xyz").exists()]
    if dangling:
        sys.exit(
            f"ABORT: {len(dangling)} dangling cohort symlinks -- the dataset has vanished; restore it first"
        )
    mols = cohort if control else [m for m in cohort if m in gen_have]
    if args.n:
        mols = list(
            np.random.default_rng(args.seed).choice(
                mols, size=min(args.n, len(mols)), replace=False
            )
        )
    print(
        f"cohort={len(cohort)}  with_generated_xyz={len(gen_have & set(cohort))}  scoring={len(mols)}  control={control}"
    )
    args.out.mkdir(parents=True, exist_ok=True)
    name = f"g_verdict_control-{control}" if control else "g_verdict"
    parts = [args.out / f".{name}.part{k}" for k in range(args.cpu)]
    t0 = time.time()
    # Workers are forked from the MAIN thread and polled: forking from helper threads deadlocked
    # the children at exit (inherited locks), with every verdict already written.
    todo = {k: mols[k :: args.cpu] for k in range(args.cpu)}
    procs: dict = {}
    for p in parts:
        p.unlink(missing_ok=True)
    while any(todo.values()) or procs:
        for k in range(args.cpu):
            p = procs.get(k)
            if p is not None and p.is_alive():
                continue
            if p is not None:
                p.join()
                del procs[k]
                todo[k] = _reconcile(parts[k], todo[k])
            if todo[k]:
                procs[k] = mp.Process(
                    target=_worker,
                    args=(
                        todo[k],
                        parts[k],
                        args.cohort,
                        args.sweep,
                        control,
                        args.seed,
                        args.timeout,
                    ),
                )
                procs[k].start()
        time.sleep(0.2)
    recs = []
    for p in parts:
        if p.exists():
            recs += [json.loads(ln) for ln in p.read_text().splitlines() if not ln.startswith("#")]
            p.unlink()
    recs.sort(key=lambda r: r["molecule"])
    out = args.out / f"{name}.jsonl"
    with open(out, "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
        fh.write(f"#DONE {len(recs)}\n")
    print(f"wrote {out}  n={len(recs)}  {time.time() - t0:.0f}s")
    summarize(recs, control)


def summarize(recs, control):
    n = len(recs)
    g = Counter(r["graph"] for r in recs)
    s = Counter(r["stereo"] for r in recs)
    print(f"\nDENOMINATOR n={n}")
    print("graph :", dict(g.most_common()))
    print("stereo:", dict(s.most_common()))
    if not control:
        return
    iso = ("ISO", "ISO_MARGINAL", "ISO_CLASH")
    usable = [r for r in recs if r["graph"] not in ("CONTROL_NA", "UNREADABLE")]
    if control == "squeeze":
        bad = [
            r
            for r in usable
            if r["graph"] not in iso or r["stereo"] not in ("SAME", "NO_STEREO", "SAME_PARTIAL")
        ]
    elif control == "identity":
        bad = [
            r
            for r in usable
            if r["graph"] != "ISO"
            or r["stereo"] not in ("SAME", "NO_STEREO")
            or r.get("sphere") != "SAME"
            or r.get("h_decoration") != "SAME"
        ]
    elif control == "mirror":
        bad = []
        for r in usable:
            want = ("MIRROR",) if r.get("input_chiral_by_ruler") else ("SAME", "NO_STEREO")
            if r["graph"] != "ISO" or r["stereo"] not in want or r.get("sphere") != "SAME":
                bad.append(r)
        ch = sum(1 for r in usable if r.get("input_chiral_by_ruler"))
        print(f"ruler-chiral inputs: {ch}/{len(usable)} = {100 * ch / max(1, len(usable)):.1f}%")
    else:
        bad = [r for r in usable if r["graph"] in iso]
    print(
        f"CONTROL {control}: usable={len(usable)}  FAIL={len(bad)}  ->  {'PASS' if not bad and usable else 'FAIL'}"
    )
    for r in bad[:12]:
        print("   ", r["molecule"], r["graph"], r["stereo"], r.get("sphere"), r.get("h_decoration"))


def crosstab(args):
    recs = [
        json.loads(ln)
        for ln in (args.out / "g_verdict.jsonl").read_text().splitlines()
        if not ln.startswith("#")
    ]
    verdict = {r["molecule"]: r for r in recs}
    rows = json.loads((args.sweep / "bucket_report_honest.json").read_text())
    rows = rows["records"] if isinstance(rows, dict) else rows
    print(f"DENOMINATOR: bucket rows={len(rows)}  ruler verdicts={len(verdict)}")

    def cell(r):
        v = verdict.get(r["molecule"].removesuffix(".xyz"))
        if v is None:
            return "NO_STRUCTURE"
        if v["graph"] not in ("ISO", "ISO_MARGINAL", "ISO_CLASH"):
            return v["graph"]
        st = v["stereo"]
        tag = {"ISO": "ISO", "ISO_MARGINAL": "ISO~", "ISO_CLASH": "ISO!"}[v["graph"]]
        if v.get("sphere") == "ARRANGEMENT_DIFF":
            return f"{tag}+ARRANGEMENT_DIFF"
        return f"{tag}+{st}"

    tab: dict = {}
    for r in rows:
        b = r["bucket"] if r["bucket"] != "key_equal" else f"key_equal/{r.get('subclass')}"
        tab.setdefault(b, Counter())[cell(r)] += 1
    cols = sorted({c for t in tab.values() for c in t})
    print("| bucket | n | " + " | ".join(cols) + " |")
    print("|---|---|" + "---|" * len(cols))
    for b in sorted(tab, key=lambda b: -sum(tab[b].values())):
        print(
            f"| {b} | {sum(tab[b].values())} | "
            + " | ".join(str(tab[b].get(c, 0)) for c in cols)
            + " |"
        )
    (args.out / "g_crosstab.json").write_text(
        json.dumps({b: dict(t) for b, t in tab.items()}, indent=1, sort_keys=True)
    )


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--cohort", type=Path, default=DEF_COHORT)
    ap.add_argument("--sweep", type=Path, default=DEF_SWEEP)
    ap.add_argument("--out", type=Path, default=DEF_OUT)
    ap.add_argument("--control", choices=["identity", "mirror", "squeeze", "detach", "delete"])
    ap.add_argument("--crosstab", action="store_true")
    ap.add_argument("--n", type=int, default=0, help="random subsample (0 = all)")
    ap.add_argument("--cpu", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--timeout", type=int, default=90, help="per-molecule hard kill, seconds")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    for p in (args.cohort, args.sweep):
        if not p.is_absolute():
            sys.exit(f"ABORT: {p} is relative -- pass an absolute path (worktrees have no dataset)")
    crosstab(args) if args.crosstab else run(args)


if __name__ == "__main__":
    main()
