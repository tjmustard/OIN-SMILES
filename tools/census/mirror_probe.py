"""Encoder-only probe: is ``E(mirror(x)) == E(x)`` for molecules the ruler calls ACHIRAL?

A canonical hash must give an achiral molecule ONE string: its mirror image is the same
molecule. This encodes each input and its z-reflection and compares bytes, for three groups
picked from the census verdicts (run ``g_vs_input.py`` and its ``--control mirror`` first):

* ruler-achiral molecules that FAILED the round trip as ``key_equal/slot_renumber``
* ruler-achiral molecules that PASSED                       (the at-risk base rate)
* ruler-CHIRAL molecules that passed                        (control: the string MUST differ)

Unlike the ruler this DOES import ``oinsmiles`` -- it measures the encoder. It prints
``oinsmiles.__file__`` so a run against the wrong source tree is visible, not silent:

    PYTHONPATH=$PWD/src <main>/.venv/bin/python tools/census/mirror_probe.py
"""

import json
import os
import random
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import TimeoutError as FT

D = "/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset"


def enc_pair(m):
    import oinsmiles
    from oinsmiles import XYZToSMILES

    src = f"{D}/cohort-v0.4.5-5k/{m}.xyz"
    lines = open(src).read().splitlines()
    n = int(lines[0])
    out = lines[:2]
    for ln in lines[2 : 2 + n]:
        p = ln.split()
        out.append(f"{p[0]} {p[1]} {p[2]} {-float(p[3]):.8f}")
    with tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False) as fh:
        fh.write("\n".join(out) + "\n")
        tmp = fh.name
    try:
        a = XYZToSMILES().convert(src)
        b = XYZToSMILES().convert(tmp)
    except Exception as e:
        return m, None, None, f"{type(e).__name__}"
    finally:
        os.unlink(tmp)
    return m, a, b, oinsmiles.__file__


if __name__ == "__main__":
    V = {}
    for ln in open(f"{D}/results-census/g_verdict.jsonl"):
        if ln[0] != "#":
            r = json.loads(ln)
            V[r["molecule"]] = r
    Mi = {}
    for ln in open(f"{D}/results-census/g_verdict_control-mirror.jsonl"):
        if ln[0] != "#":
            r = json.loads(ln)
            Mi[r["molecule"]] = r
    B = {
        r["molecule"].removesuffix(".xyz"): r
        for r in json.load(open(f"{D}/results-v0.4.14-sweep/bucket_report_honest.json"))
    }
    achiral = [m for m, r in Mi.items() if not r.get("input_chiral_by_ruler")]
    rng = random.Random(42)
    fails = [
        m
        for m in achiral
        if B[m]["bucket"] == "key_equal" and B[m].get("subclass") == "slot_renumber"
    ]
    passes = rng.sample([m for m in achiral if B[m]["bucket"] == "byte_exact"], 600)
    chiral_pass = rng.sample(
        [
            m
            for m, r in Mi.items()
            if r.get("input_chiral_by_ruler") and B[m]["bucket"] == "byte_exact"
        ],
        200,
    )
    groups = {
        "achiral & slot_renumber FAIL": fails,
        "achiral & byte_exact PASS": passes,
        "CHIRAL & byte_exact PASS (control: must differ)": chiral_pass,
    }
    allm = sorted(set(fails + passes + chiral_pass))
    res = {}
    ex = ProcessPoolExecutor(10)
    futs = {m: ex.submit(enc_pair, m) for m in allm}
    for m, f in futs.items():
        try:
            res[m] = f.result(timeout=240)
        except FT:
            res[m] = (m, None, None, "TIMEOUT")
    print("oinsmiles from:", next((r[3] for r in res.values() if r[1]), "?"))
    out = {}
    for name, ms in groups.items():
        c = Counter()
        for m in ms:
            _, a, b, _ = res[m]
            c[
                "encode_error"
                if a is None or b is None
                else ("E(mirror)==E(x)" if a == b else "E(mirror)!=E(x)")
            ] += 1
        print(f"{name:50s} n={len(ms):4d}  {dict(c)}")
        out[name] = dict(c)
    json.dump(
        {
            "groups": out,
            "rows": {
                m: {"same": (r[1] == r[2]) if r[1] and r[2] else None} for m, r in res.items()
            },
        },
        open(f"{D}/results-census/mirror_probe.json", "w"),
    )
    os._exit(0)
