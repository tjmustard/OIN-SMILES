"""I2 of the census: is the ENCODER self-consistent on every input? No generator anywhere.

One molecule, presented eleven ways; a canonical lossless hash has a KNOWN answer for each:

    base     E(x)                                  the reference string
    again    E(x) a second time                    MUST equal base      (determinism)
    rewrite  E(x re-written by this script)        MUST equal base      (controls the file writer)
    rot0     E(R x), R a proper rotation           MUST equal base
    renum0-2 E(x, atoms permuted)                  MUST equal base
    noise0-2 E(x + N(0, sigma)), sigma = 0.02 A    SHOULD equal base    (perception robustness)
    mirror   E(x, z -> -z)                         achiral-by-ruler: MUST equal base
                                                   chiral-by-ruler:  MUST differ

    python tools/census/e_selfconsistency.py --only FEQFIS_comp_0,ABAZIO_comp_0   # smoke
    python tools/census/e_selfconsistency.py                                       # all 5,000
    python tools/census/e_selfconsistency.py --summarize                           # tables + the 2x2

The mirror column is joined to the I1 ruler's ``g_verdict_control-mirror.jsonl::
input_chiral_by_ruler`` (a verdict that never imported ``oinsmiles``) to give the 2x2:
achiral & differs = NON-CANONICAL; chiral & same = NON-INJECTIVE.

WHAT A BROKEN VERSION WOULD PRINT
---------------------------------
A probe whose variants were never really transformed prints "100% stable" -- the same line a
perfect encoder prints. So: every variant is asserted to BE its transform before it is encoded
(non-identity permutation, det = +1 / -1, 0 < max displacement); ``rewrite`` separates "the
encoder moved" from "my file writer moved it"; the chiral half of the 2x2 shows the mirror column
FIRES; ``--summarize`` reports agreement with C1's independent ``mirror_probe.json``; and the
renumber column has a must-drift control -- a v0.4.5 drifter with the levers that fixed it OFF:

    OIN_STABLE_STEREO=0 OIN_CANONICAL_BODY=0 OIN_CANONICAL_SLOTS=0 OIN_CANONICAL_PERCEPTION=0 \\
      python tools/census/e_selfconsistency.py --only FEQFIS_comp_0 --out <scratch>

This DOES import ``oinsmiles`` (it measures the encoder) and prints ``oinsmiles.__file__`` so a
run against the wrong source tree is visible: ``sys.path.insert`` beats ``PYTHONPATH``, and six
A/B arms once ran main's code and printed a flawless null. Run with the MAIN checkout's venv:

    PYTHONPATH=$PWD/src <main>/.venv/bin/python tools/census/e_selfconsistency.py ...

Encodes can hang inside C++ (max observed 185 s, v0.4.8), so each encode runs under SIGALRM at
SIG_DFL in a worker forked from the MAIN thread; a killed encode is recorded as TIMEOUT and the
worker is restarted on the next task (pattern copied from ``g_vs_input.py``).
"""

from __future__ import annotations

import argparse
import contextlib
import json
import multiprocessing as mp
import os
import re
import shutil
import signal
import sys
import time
import zlib
from collections import Counter
from pathlib import Path

import numpy as np

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
DEF_COHORT = MAIN / "cohort-v0.4.5-5k"
DEF_SWEEP = MAIN / "results-v0.4.14-sweep"
DEF_OUT = MAIN / "results-census"

TAGS = ["base", "again", "rewrite", "rot0"]
TAGS += [f"renum{i}" for i in range(3)] + [f"noise{i}" for i in range(3)] + ["mirror"]
SIGMA = 0.02
_SLOT_RE = re.compile(r"\{(\d+)([<>^]?)\}")


@contextlib.contextmanager
def _silence_fds():
    """Redirect C-level stdout/stderr to devnull (openbabel prints distance warnings)."""
    with open(os.devnull, "w") as devnull:
        old_out, old_err = os.dup(1), os.dup(2)
        try:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            os.dup2(old_out, 1)
            os.dup2(old_err, 2)
            os.close(old_out)
            os.close(old_err)


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    syms, xyz = [], []
    for ln in lines[2 : 2 + n]:
        p = ln.split()
        syms.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, np.asarray(xyz, dtype=float), (lines[1] if len(lines) > 1 else "")


def write_xyz(path, syms, xyz, comment):
    with open(path, "w") as fh:
        fh.write(f"{len(syms)}\n{comment}\n")
        for s, c in zip(syms, xyz):
            fh.write(f"{s:<3} {c[0]:>14.8f} {c[1]:>14.8f} {c[2]:>14.8f}\n")


def _rot(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


class InstrumentError(Exception):
    """A variant was not the transform its tag names (an encoder ``assert`` is NOT this)."""


def make_variant(tag, syms, xyz, rng):
    """The transformed structure, ASSERTED to be the transform its tag names."""
    syms, out = list(syms), np.array(xyz, dtype=float)
    if tag == "rewrite":
        pass
    elif tag.startswith("rot"):
        R = _rot(rng)
        assert abs(np.linalg.det(R) - 1.0) < 1e-9 and not np.allclose(R, np.eye(3), atol=1e-3)
        out = out @ R.T
    elif tag.startswith("renum"):
        perm = rng.permutation(len(syms))
        if len(syms) > 1:
            while np.array_equal(perm, np.arange(len(syms))):
                perm = rng.permutation(len(syms))
        orig = list(syms)
        syms, out = [orig[i] for i in perm], out[perm]
        assert sorted(syms) == sorted(orig) and len(out) == len(xyz), "renumber lost an atom"
        assert len(syms) < 2 or not np.allclose(out, xyz), "renumber produced the identity"
    elif tag.startswith("noise"):
        out = out + rng.normal(scale=SIGMA, size=out.shape)
        d = np.linalg.norm(out - xyz, axis=1).max()
        assert 0.0 < d < 0.25, f"noise displacement {d}"
    elif tag == "mirror":
        out[:, 2] *= -1.0
        assert np.allclose(out[:, :2], xyz[:, :2]) and np.allclose(out[:, 2], -xyz[:, 2])
    else:
        raise ValueError(tag)
    return syms, out


def _encode_task(mol, tag, cohort, tmpdir, seed):
    from oinsmiles import XYZToSMILES

    src = cohort / f"{mol}.xyz"
    if tag in ("base", "again"):
        path = src
    else:
        syms, xyz, comment = read_xyz(src)
        rng = np.random.default_rng([seed, zlib.crc32(mol.encode()), zlib.crc32(tag.encode())])
        try:
            vs, vx = make_variant(tag, syms, xyz, rng)
        except AssertionError as exc:  # the INSTRUMENT is broken, not the encoder
            raise InstrumentError(str(exc)) from exc
        # same basename as the input, so the file name is not a variable
        path = tmpdir / f"{mol}.xyz"
        write_xyz(path, vs, vx, comment)
    with _silence_fds():
        return XYZToSMILES().convert(str(path))


def _worker(tasks, part, cohort, tmpdir, seed, timeout):
    # SIGALRM stays at SIG_DFL on purpose: a Python handler cannot interrupt C++, but the
    # kernel's default action kills the process whatever it is doing.
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    from oinsmiles import XYZToSMILES  # noqa: F401  pay the import OUTSIDE the first alarm

    tmpdir.mkdir(parents=True, exist_ok=True)
    no_base = set()
    with open(part, "a") as fh:
        for mol, tag in tasks:
            if mol in no_base:
                continue
            fh.write(f"#START {mol} {tag}\n")
            fh.flush()
            t0 = time.time()
            signal.alarm(timeout)
            try:
                oin, err = _encode_task(mol, tag, cohort, tmpdir, seed), None
            except InstrumentError as exc:
                oin, err = None, f"INSTRUMENT: {exc}"
            except Exception as exc:
                oin, err = None, f"{type(exc).__name__}: {exc}"[:300]
            signal.alarm(0)
            if tag == "base" and not oin:
                no_base.add(mol)
            fh.write(
                json.dumps(
                    {
                        "molecule": mol,
                        "tag": tag,
                        "oin": oin,
                        "error": err,
                        "t": round(time.time() - t0, 3),
                    }
                )
                + "\n"
            )
            fh.flush()
        os.fsync(fh.fileno())
    os._exit(0)  # skip interpreter finalizers (BLAS thread pools) on the way out


def _reconcile(part, todo, timeout):
    """After a worker exits: drop finished tasks; a #START with no result line was killed."""
    done, started, no_base = set(), None, set()
    if part.exists():
        for ln in part.read_text().splitlines():
            if ln.startswith("#START "):
                started = tuple(ln.split()[1:3])
            else:
                r = json.loads(ln)
                done.add((r["molecule"], r["tag"]))
                if r["tag"] == "base" and not r["oin"]:
                    no_base.add(r["molecule"])
    if started is not None and started not in done:
        with open(part, "a") as fh:
            fh.write(
                json.dumps(
                    {
                        "molecule": started[0],
                        "tag": started[1],
                        "oin": None,
                        "error": "TIMEOUT",
                        "t": float(timeout),
                    }
                )
                + "\n"
            )
        done.add(started)
        if started[1] == "base":
            no_base.add(started[0])
    return [t for t in todo if t not in done and t[0] not in no_base]


def _relation(a, b):
    """Sweep taxonomy for a pair of strings (``roundtrip_bucket_report.py::classify``)."""
    from oinsmiles.oin.compare import canonical_roundtrip_key, normalize_oin_for_comparison

    if a == b:
        return "byte_exact"
    try:
        k1, k2 = canonical_roundtrip_key(a), canonical_roundtrip_key(b)
    except Exception:
        return "key_error"
    if k1 != k2:
        return "facmer_divergent" if (k1[0] == k2[0] and k1[2] == k2[2]) else "structural"
    n1, n2 = normalize_oin_for_comparison(a.strip()), normalize_oin_for_comparison(b.strip())
    f1, f2 = [f for f in n1.split(".") if f], [f for f in n2.split(".") if f]
    if sorted(f1) == sorted(f2) and f1 != f2:
        return "key_equal/fragment_reorder"
    if _SLOT_RE.sub(r"{\2}", n1) == _SLOT_RE.sub(r"{\2}", n2):
        return "key_equal/slot_renumber"
    sw = re.compile(r"\{(\d+)[<>^]\}")
    if sw.sub(r"{\1}", n1) == sw.sub(r"{\1}", n2):
        return "key_equal/winding_star_drift"
    return "key_equal/rdkit_canonical"


def assemble(lines, mols):
    """Per-(molecule, tag) lines -> one record per molecule."""
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")  # the key parses fragments; its warnings are not findings
    by: dict = {m: {} for m in mols}
    for r in lines:
        by[r["molecule"]][r["tag"]] = r
    recs = []
    for m in mols:
        got = by[m]
        base = got.get("base", {})
        rec = {
            "molecule": m,
            "base": base.get("oin"),
            "base_error": base.get("error"),
            "t_base": base.get("t"),
            "t_total": round(sum(r["t"] for r in got.values()), 3),
            "rel": {},
            "diff": {},
        }
        if rec["base"]:
            for tag in TAGS[1:]:
                r = got.get(tag)
                if r is None:
                    rec["rel"][tag] = "MISSING"
                elif not r["oin"]:
                    rec["rel"][tag] = (
                        "TIMEOUT"
                        if r["error"] == "TIMEOUT"
                        else "INSTRUMENT"
                        if (r["error"] or "").startswith("INSTRUMENT")
                        else "encode_fail"
                    )
                    rec["diff"][tag] = r["error"]
                else:
                    rec["rel"][tag] = _relation(rec["base"], r["oin"])
                    if r["oin"] != rec["base"]:
                        rec["diff"][tag] = r["oin"]
        recs.append(rec)
    return recs


def run(args):
    import oinsmiles

    print("oinsmiles from:", oinsmiles.__file__)
    print("levers in env :", {k: v for k, v in sorted(os.environ.items()) if k.startswith("OIN_")})
    cohort = sorted(p.stem for p in args.cohort.glob("*.xyz"))
    dangling = [m for m in cohort if not (args.cohort / f"{m}.xyz").exists()]
    if dangling:
        sys.exit(f"ABORT: {len(dangling)} dangling cohort symlinks -- restore the dataset first")
    mols = cohort
    if args.only:
        want = {m.strip().removesuffix(".xyz") for m in args.only.split(",") if m.strip()}
        mols = [m for m in cohort if m in want]
        if want - set(mols):
            sys.exit(f"ABORT: --only named molecules not in the cohort: {sorted(want - set(mols))}")
    elif args.n:
        rng = np.random.default_rng(args.seed)
        mols = sorted(rng.choice(cohort, size=min(args.n, len(cohort)), replace=False))
    print(f"cohort={len(cohort)}  probing={len(mols)}  encodes={len(mols) * len(TAGS)}")
    args.out.mkdir(parents=True, exist_ok=True)
    name = args.name
    parts = [args.out / f".{name}.part{k}" for k in range(args.cpu)]
    tmproot = args.out / f".{name}.tmp"
    t0 = time.time()
    # molecule-major order: a molecule's eleven encodes stay on one worker, base first
    todo = {k: [(m, t) for m in mols[k :: args.cpu] for t in TAGS] for k in range(args.cpu)}
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
                todo[k] = _reconcile(parts[k], todo[k], args.timeout)
            if todo[k]:
                procs[k] = mp.Process(
                    target=_worker,
                    args=(
                        todo[k],
                        parts[k],
                        args.cohort,
                        tmproot / str(k),
                        args.seed,
                        args.timeout,
                    ),
                )
                procs[k].start()
        time.sleep(0.2)
    lines = []
    for p in parts:
        if p.exists():
            lines += [json.loads(ln) for ln in p.read_text().splitlines() if ln[0] != "#"]
    recs = assemble(lines, mols)
    out = args.out / f"{name}.jsonl"
    with open(out, "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
        fh.write(f"#DONE {len(recs)}\n")
    for p in parts:
        p.unlink(missing_ok=True)
    shutil.rmtree(tmproot, ignore_errors=True)
    print(f"wrote {out}  n={len(recs)}  {time.time() - t0:.0f}s wall")
    summarize(args)


def _load(path):
    rows = [ln for ln in Path(path).read_text().splitlines()]
    if not rows or not rows[-1].startswith("#DONE"):
        sys.exit(f"ABORT: {path} has no #DONE sentinel -- the run is unfinished or was killed")
    return [json.loads(ln) for ln in rows if ln[0] != "#"]


def _group(rec, prefix):
    """Worst relation over a transform's k trials: (byte_stable, key_stable, relations)."""
    rels = [v for t, v in rec["rel"].items() if t.startswith(prefix)]
    byte = all(v == "byte_exact" for v in rels)
    key = all(v == "byte_exact" or v.startswith("key_equal") for v in rels)
    return byte, key, rels


def summarize(args):
    recs = _load(args.out / f"{args.name}.jsonl")
    n = len(recs)
    enc = [r for r in recs if r["base"]]
    print(f"\nDENOMINATOR n={n}  base encoded={len(enc)}  base failed={n - len(enc)}")
    print(
        "base failures:",
        dict(Counter((r["base_error"] or "?")[:40] for r in recs if not r["base"]).most_common(8)),
    )
    bad_inst = [r["molecule"] for r in enc if "INSTRUMENT" in r["rel"].values()]
    miss = [r["molecule"] for r in enc if "MISSING" in r["rel"].values()]
    print(
        f"INSTRUMENT assertion failures={len(bad_inst)}  MISSING tags={len(miss)}  (both must be 0)"
    )
    print(f"CPU-h={sum(r['t_total'] for r in recs) / 3600:.2f}")

    print("\n| transform | k | byte-stable | key-stable | not byte-stable: worst relations |")
    print("|---|---|---|---|---|")
    flags = {}
    for label, prefix, k in (
        ("again (determinism)", "again", 1),
        ("rewrite (writer control)", "rewrite", 1),
        ("rotate", "rot", 1),
        ("renumber", "renum", 3),
        (f"noise sigma={SIGMA}", "noise", 3),
        ("mirror (ALL molecules)", "mirror", 1),
    ):
        b = kk = 0
        c = Counter()
        for r in enc:
            byte, key, rels = _group(r, prefix)
            flags.setdefault(r["molecule"], {})[prefix] = (byte, key)
            b += byte
            kk += key
            if not byte:
                c.update({v for v in rels if v != "byte_exact"})
        print(
            f"| {label} | {k} | {b}/{len(enc)} ({100 * b / max(1, len(enc)):.2f}%) | "
            f"{kk}/{len(enc)} ({100 * kk / max(1, len(enc)):.2f}%) | {dict(c.most_common())} |"
        )

    # cross-process determinism, free: the v0.4.14 sweep encoded the same inputs in another
    # process (another hash seed) and an older tree. Read the JSON-derived report, never *.oin.
    rep = args.sweep / "bucket_report_honest.json"
    B = {}
    if rep.exists():
        B = {r["molecule"].removesuffix(".xyz"): r for r in json.loads(rep.read_text())}
        both = [r for r in enc if B.get(r["molecule"], {}).get("smiles_1")]
        same = sum(1 for r in both if B[r["molecule"]]["smiles_1"] == r["base"])
        print(
            f"\nbase == v0.4.14 sweep smiles_1: {same}/{len(both)}  (differences = code drift OR cross-process nondeterminism)"
        )

    ruler = args.out / "g_verdict_control-mirror.jsonl"
    if not ruler.exists():
        print(f"\n(no {ruler.name}; run g_vs_input.py --control mirror for the 2x2)")
        return
    chiral = {r["molecule"]: r.get("input_chiral_by_ruler") for r in _load(ruler)}
    tab = Counter()
    cells: dict = {}
    for r in enc:
        rel = r["rel"].get("mirror")
        col = (
            "same"
            if rel == "byte_exact"
            else "differs"
            if rel not in ("TIMEOUT", "encode_fail", "MISSING", "INSTRUMENT")
            else "no_mirror_encode"
        )
        row = {True: "chiral", False: "achiral", None: "no_ruler_verdict"}[
            chiral.get(r["molecule"])
        ]
        tab[(row, col)] += 1
        cells.setdefault((row, col), []).append(r)
    print(f"\nTHE 2x2 (ruler chirality x E(mirror) vs E(x)); joined={sum(tab.values())}/{len(enc)}")
    print("| ruler | E(mirror)==E(x) | E(mirror)!=E(x) | no mirror encode |")
    print("|---|---|---|---|")
    for row in ("achiral", "chiral", "no_ruler_verdict"):
        print(
            f"| {row} | {tab[(row, 'same')]} | {tab[(row, 'differs')]} | {tab[(row, 'no_mirror_encode')]} |"
        )
    nc, ni = cells.get(("achiral", "differs"), []), cells.get(("chiral", "same"), [])
    ach = tab[("achiral", "same")] + tab[("achiral", "differs")]
    chi = tab[("chiral", "same")] + tab[("chiral", "differs")]
    print(
        f"NON-CANONICAL (achiral & differs): {len(nc)}/{ach} = {100 * len(nc) / max(1, ach):.1f}% of achiral"
    )
    print(
        f"NON-INJECTIVE (chiral & same)    : {len(ni)}/{chi} = {100 * len(ni) / max(1, chi):.1f}% of chiral"
    )
    print(f"mirror column FIRES on chiral    : {tab[('chiral', 'differs')]}/{chi}")
    print(
        "  achiral & differs, by relation :",
        dict(Counter(r["rel"]["mirror"] for r in nc).most_common()),
    )

    def other(r):
        return all(
            flags[r["molecule"]][p][0] for p in ("again", "rewrite", "rot", "renum", "noise")
        )

    print(
        f"  ...of which stable under EVERY other transform (a pure reflection defect): {sum(map(other, nc))}/{len(nc)}"
    )

    if B:
        print(
            "\n| honest bucket (v0.4.14) | n enc | any proper-transform byte drift | achiral & mirror differs | chiral & mirror same |"
        )
        print("|---|---|---|---|---|")
        per: dict = {}
        for r in enc:
            b = B.get(r["molecule"])
            if not b:
                continue
            name = b["bucket"] if b["bucket"] != "key_equal" else f"key_equal/{b.get('subclass')}"
            per.setdefault(name, []).append(r)
        ncs, nis = {r["molecule"] for r in nc}, {r["molecule"] for r in ni}
        for name, rs in sorted(per.items(), key=lambda kv: -len(kv[1])):
            drift = sum(1 for r in rs if not other(r))
            print(
                f"| {name} | {len(rs)} | {drift} ({100 * drift / len(rs):.1f}%) | "
                f"{sum(1 for r in rs if r['molecule'] in ncs)} | {sum(1 for r in rs if r['molecule'] in nis)} |"
            )

    c1 = args.out / "mirror_probe.json"
    if c1.exists():
        rows = json.loads(c1.read_text())["rows"]
        mine = {r["molecule"]: r["rel"].get("mirror") == "byte_exact" for r in enc}
        cmp = [
            (m, v["same"], mine[m]) for m, v in rows.items() if v["same"] is not None and m in mine
        ]
        agree = sum(1 for _, a, b in cmp if a == b)
        print(f"\nagreement with C1 mirror_probe.json (independent script): {agree}/{len(cmp)}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--cohort", type=Path, default=DEF_COHORT)
    ap.add_argument("--sweep", type=Path, default=DEF_SWEEP)
    ap.add_argument("--out", type=Path, default=DEF_OUT)
    ap.add_argument(
        "--name", default="e_selfconsistency", help="output stem (use another for smokes)"
    )
    ap.add_argument("--only", help="comma-separated molecule ids")
    ap.add_argument("--n", type=int, default=0, help="random subsample (0 = all)")
    ap.add_argument("--summarize", action="store_true")
    ap.add_argument("--cpu", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--timeout", type=int, default=300, help="per-ENCODE hard kill, seconds")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    for p in (args.cohort, args.sweep, args.out):
        if not p.is_absolute():
            sys.exit(f"ABORT: {p} is relative -- pass an absolute path (worktrees have no dataset)")
    summarize(args) if args.summarize else run(args)


if __name__ == "__main__":
    main()
