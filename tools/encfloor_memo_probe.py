#!/usr/bin/env python
"""Does the ``AC2BO`` memo hit across two *differently numbered* copies of one molecule?

``docs/agentic-notes/v0.4.5/V045_STATUS_2026-07-25.md`` measured that a generated conformer's atom order never
matches the input (0/36) and concluded the memo "cannot hit on the generated side, because
its key is a permutation away". That measurement predates the canonicalising wrapper: today
``AC2BO`` computes ``perm = _canonical_atom_permutation(AC, atoms)`` and calls
``_AC2BO_core`` on the **permuted** matrix, and ``_ac2bo_memo_anchor`` tags on the bytes of
*that permuted* matrix. If the canonical permutation is graph-determined, two encodes of the
same graph should share a tag and the second should reuse the first's entries.

Renumbering the input XYZ is the **generator-free upper bound** on that question: it holds
the graph exactly fixed (unlike a generated conformer, whose perceived graph matched the
input only 16/36 = 44 % of the time) while destroying the atom order completely. If the memo
does not hit here, it cannot hit on the generated side either, and LEAD 3 is dead without
spending generator time. If it does hit, 44 % bounds the population that could benefit.

Counters, not seconds: ``AC2BO_STATS['matching_calls']`` and ``['candidates']`` are exact and
load-independent. Two arms per molecule, in the SAME process (a memo hit is the thing being
measured, so a fresh process would destroy it):

  arm REUSE  encode original, then encode renumbered copy WITHOUT clearing the memo
  arm COLD   encode original, then clear the memo, then encode the renumbered copy

``matching_calls`` on encode 2 is the payload. REUSE << COLD means the memo hits.
The OIN strings are compared too: canonical perception claims they should be equal.

    $V tools/encfloor_memo_probe.py --dataset $DS --molecules QIDKIZ_comp_0 --seed 7
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


def find_xyz(dataset_dir: str, molecule: str) -> str:
    if os.path.exists(molecule):
        return molecule
    for sub in ("cat", "photo", "regression_inputs", "."):
        p = os.path.join(dataset_dir, sub, f"{molecule}.xyz")
        if os.path.exists(p):
            return p
    p = os.path.join("tests/fixtures", f"{molecule}.xyz")
    if os.path.exists(p):
        return p
    raise FileNotFoundError(f"{molecule} not found under {dataset_dir}")


def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0].split()[0])
    syms, coords = [], []
    for line in lines[2 : 2 + n]:
        p = line.split()
        syms.append(p[0])
        coords.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, np.asarray(coords, dtype=float), (lines[1].rstrip("\n") if len(lines) > 1 else "")


def write_xyz(path, syms, coords, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n{comment}\n")
        for s, c in zip(syms, coords):
            f.write(f"{s:<3} {c[0]:>14.8f} {c[1]:>14.8f} {c[2]:>14.8f}\n")


def encode(path, loc):
    """Encode one file; return (sha, stats delta, wall). Wall is ADVISORY."""
    from oinsmiles import XYZToSMILES

    before = dict(loc.AC2BO_STATS)
    t0 = time.perf_counter()
    err = None
    oin = None
    try:
        oin = XYZToSMILES().convert(path)
    except Exception as exc:  # noqa: BLE001
        err = f"{type(exc).__name__}: {exc}"
    wall = time.perf_counter() - t0
    delta = {k: loc.AC2BO_STATS[k] - before.get(k, 0) for k in loc.AC2BO_STATS}
    sha = None if oin is None else hashlib.sha256(oin.encode()).hexdigest()[:16]
    return sha, delta, wall, err


def probe(path, name, seed, tmpdir):
    from oinsmiles.utils import perception_core as loc

    syms, coords, comment = read_xyz(path)
    rng = random.Random(seed)
    order = list(range(len(syms)))
    rng.shuffle(order)
    ren = os.path.join(tmpdir, f"{name}_renum.xyz")
    write_xyz(ren, [syms[i] for i in order], coords[order], comment)

    out = {"molecule": name, "n_atoms": len(syms), "seed": seed}
    for arm in ("REUSE", "COLD"):
        loc._ac2bo_memo_clear()
        loc.reset_ac2bo_stats()
        sha1, d1, w1, e1 = encode(path, loc)
        if arm == "COLD":
            loc._ac2bo_memo_clear()
        sha2, d2, w2, e2 = encode(ren, loc)
        out[arm] = {
            "sha_orig": sha1,
            "sha_renum": sha2,
            "oin_equal": sha1 == sha2,
            "enc1": {
                "matching": d1["matching_calls"],
                "candidates": d1["candidates"],
                "ac2bo": d1["ac2bo_calls"],
                "wall": w1,
                "err": e1,
            },
            "enc2": {
                "matching": d2["matching_calls"],
                "candidates": d2["candidates"],
                "ac2bo": d2["ac2bo_calls"],
                "wall": w2,
                "err": e2,
            },
            "memo_entries_after": loc._ac2bo_memo_entries(),
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="tmCAT-tmPHOTO_xyz_dataset")
    ap.add_argument("--molecules", required=True)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args()

    names = [x.strip() for x in args.molecules.split(",") if x.strip()]
    print(f"loadavg_at_start={os.getloadavg()}", flush=True)
    results = []
    with tempfile.TemporaryDirectory() as td:
        for name in names:
            try:
                path = find_xyz(args.dataset, name)
            except FileNotFoundError as exc:
                print(f"SKIP {name}: {exc}", flush=True)
                continue
            r = probe(path, name, args.seed, td)
            results.append(r)
            print(f"\n=== {name} ({r['n_atoms']} atoms) ===", flush=True)
            for arm in ("REUSE", "COLD"):
                a = r[arm]
                print(
                    f"  {arm:6s} enc1 matching={a['enc1']['matching']:>7d} "
                    f"cand={a['enc1']['candidates']:>7d} | "
                    f"enc2 matching={a['enc2']['matching']:>7d} "
                    f"cand={a['enc2']['candidates']:>7d} | "
                    f"oin_equal={a['oin_equal']} entries={a['memo_entries_after']} "
                    f"(wall1={a['enc1']['wall']:.2f}s wall2={a['enc2']['wall']:.2f}s ADVISORY)",
                    flush=True,
                )
            reuse2 = r["REUSE"]["enc2"]["matching"]
            cold2 = r["COLD"]["enc2"]["matching"]
            verdict = "HIT" if cold2 and reuse2 < cold2 else ("no-op" if cold2 == 0 else "MISS")
            print(f"  -> memo on encode2: {verdict}  ({reuse2} vs {cold2} matching calls)")
            print(f"#DONE {len(results)}", flush=True)

    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump(results, f, indent=1)
        print(f"wrote {args.json_out}")
    print(f"#DONE_ALL {len(results)}/{len(names)}", flush=True)


if __name__ == "__main__":
    main()
