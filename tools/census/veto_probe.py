"""WHY does ``E`` give an ACHIRAL molecule two strings? Test one named mechanism, with controls.

HYPOTHESIS. ``oin/fold_parity.py::resolve`` vetoes the donor fold when

    S_rot(x) != S_rot(mirror x)   and   S_fold(x) == S_fold(mirror x)

and reads the left conjunct as "x and its mirror are an enantiomer pair". For an ACHIRAL molecule
that premise is false: the rotation-only labeling differing between the two hands IS the
canonicality defect, the fold is the repair, and the veto undoes it -- on ONE hand only (the
other hand's raw labeling already equals the folded form, so the fold is inactive there and
nothing is vetoed). One molecule, two strings.

Three arms, each encoding x and its z-mirror, recording ``fold_parity.last_outcome()``:

    default          fold ON,  veto ON    (shipped)
    veto_off         fold ON,  veto OFF   the fold alone
    fold_off         fold OFF, veto OFF   the rotation-only labeling

Populations (ruler chirality from ``g_verdict_control-mirror.jsonl``; mirror verdicts from
``e_selfconsistency.jsonl`` when it has its #DONE, else C1's ``mirror_probe.json``):

    achiral & E(mirror) != E(x)     the defect. Hypothesis: veto_off makes them EQUAL.
    achiral & E(mirror) == E(x)     control: no arm should matter.
    chiral  & E(mirror) != E(x)     control: WHAT A BROKEN VERSION WOULD PRINT -- if veto_off
                                    never collapsed a chiral pair here, the veto lever would not
                                    be reaching the encoder and the first row would mean nothing.

    PYTHONPATH=$PWD/src <main>/.venv/bin/python tools/census/veto_probe.py [--cpu 2]

``--mode renumber`` asks the same question of the RENUMBER column: for every molecule whose
three renumberings drift at slot level (``key_equal/*``), re-encode base + the same three
permutations (identical seeds to ``e_selfconsistency.py``) under the three arms. If the drift
vanishes with the veto off, it is the same asymmetry: one presentation's raw labeling happens to
equal the folded form (fold inactive), another's does not and is vetoed.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import TimeoutError as FT
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # sibling module only
import e_selfconsistency as esc  # noqa: E402

ARMS = {
    "default": {},
    "veto_off": {"OIN_FOLD_PARITY_VETO": "0"},
    "fold_off": {"OIN_FOLD_PARITY_VETO": "0", "OIN_CANONICAL_DONOR_FOLD": "0"},
}


def _enc(path, levers):
    from oinsmiles import XYZToSMILES
    from oinsmiles.oin import fold_parity as fp

    old = {k: os.environ.get(k) for k in levers}
    os.environ.update(levers)
    fp._state.outcome = None  # last_outcome() is sticky; a veto-OFF encode never writes it
    try:
        with esc._silence_fds():
            return XYZToSMILES().convert(str(path)), fp.last_outcome()
    except Exception as exc:
        return None, f"{type(exc).__name__}"
    finally:
        for k, v in old.items():
            os.environ.pop(k, None) if v is None else os.environ.__setitem__(k, v)


def probe(mol, mode="mirror"):
    import zlib

    import numpy as np

    import oinsmiles

    src = esc.DEF_COHORT / f"{mol}.xyz"
    syms, xyz, comment = esc.read_xyz(src)
    tags = ["mirror"] if mode == "mirror" else ["renum0", "renum1", "renum2"]
    with tempfile.TemporaryDirectory() as td:
        paths = []
        for i, tag in enumerate(tags):
            rng = np.random.default_rng([42, zlib.crc32(mol.encode()), zlib.crc32(tag.encode())])
            vs, vx = esc.make_variant(tag, syms, xyz, rng)
            d = Path(td) / str(i)
            d.mkdir()
            esc.write_xyz(d / f"{mol}.xyz", vs, vx, comment)
            paths.append(d / f"{mol}.xyz")
        rec = {"molecule": mol, "src": oinsmiles.__file__}
        for arm, lev in ARMS.items():
            a, oa = _enc(src, lev)
            vs = [_enc(p, lev) for p in paths]
            rec[arm] = {
                "same": all(a == b for b, _ in vs) if a and all(b for b, _ in vs) else None,
                "outcome_x": oa,
                "outcome_m": vs[0][1] if mode == "mirror" else [o for _, o in vs],
                "x": a,
                "m": vs[0][0] if mode == "mirror" else [b for b, _ in vs],
                "rel": [esc._relation(a, b) if a and b else "encode_fail" for b, _ in vs],
            }
    return rec


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--n-control", type=int, default=150)
    ap.add_argument("--mode", choices=["mirror", "renumber"], default="mirror")
    ap.add_argument(
        "--chiral-all",
        action="store_true",
        help="mirror mode: EVERY ruler-chiral & differs molecule, no other group -- how many real "
        "pairs does the veto protect? Writes veto_probe_chiral.jsonl",
    )
    args = ap.parse_args()
    out = esc.DEF_OUT
    if args.mode == "renumber":
        return run_renumber(args, out)
    chiral = {
        r["molecule"]: r["input_chiral_by_ruler"]
        for r in esc._load(out / "g_verdict_control-mirror.jsonl")
    }
    full = out / "e_selfconsistency.jsonl"
    if full.exists() and full.read_text().rstrip().splitlines()[-1].startswith("#DONE"):
        same = {
            r["molecule"]: r["rel"]["mirror"] == "byte_exact"
            for r in esc._load(full)
            if r["base"]
            and r["rel"].get("mirror", "").split("/")[0]
            in ("byte_exact", "key_equal", "structural", "facmer_divergent")
        }
        source = full.name
    else:
        rows = json.loads((out / "mirror_probe.json").read_text())["rows"]
        same = {m: v["same"] for m, v in rows.items() if v["same"] is not None}
        source = "mirror_probe.json (C1 SAMPLE -- re-run once e_selfconsistency.jsonl has #DONE)"
    rng = random.Random(42)
    pick = lambda ms: sorted(rng.sample(ms, min(args.n_control, len(ms))))  # noqa: E731
    groups = {
        "achiral & differs (THE DEFECT)": sorted(
            m for m, s in same.items() if not s and not chiral[m]
        ),
        "achiral & same (control)": pick(sorted(m for m, s in same.items() if s and not chiral[m])),
        "chiral & differs (control: veto_off MUST collapse some)": pick(
            sorted(m for m, s in same.items() if not s and chiral[m])
        ),
    }
    if args.chiral_all:
        groups = {
            "chiral & differs (ALL): pairs the veto protects = vetoed_collapse on a hand": sorted(
                m for m, s in same.items() if not s and chiral[m]
            )
        }
    print("mirror verdicts from:", source)
    print("populations:", {k: len(v) for k, v in groups.items()})
    allm = sorted({m for ms in groups.values() for m in ms})
    ex = ProcessPoolExecutor(args.cpu)
    res = _collect(ex, allm, "mirror")
    with open(
        out / ("veto_probe_chiral.jsonl" if args.chiral_all else "veto_probe.jsonl"), "w"
    ) as fh:
        for name, ms in groups.items():
            for m in ms:
                if res[m]:
                    fh.write(json.dumps({"group": name, **res[m]}) + "\n")
        fh.write(f"#DONE {len(allm)}\n")
    for name, ms in groups.items():
        rs = [res[m] for m in ms if res[m]]
        print(f"\n== {name}   n={len(ms)}  probed={len(rs)}")
        for arm in ARMS:
            c = Counter(
                {True: "same", False: "differs", None: "encode_fail"}[r[arm]["same"]] for r in rs
            )
            print(f"   {arm:9s} E(mirror) vs E(x): {dict(c)}")
        c = Counter(
            tuple(sorted((str(r["default"]["outcome_x"]), str(r["default"]["outcome_m"]))))
            for r in rs
        )
        print("   default-arm veto outcomes {x, mirror}:")
        for k, v in c.most_common(8):
            print(f"      {v:4d}  {k}")
        if "DEFECT" in name:
            vet = [
                r
                for r in rs
                if "vetoed_collapse" in (r["default"]["outcome_x"], r["default"]["outcome_m"])
            ]
            fixed = [r for r in vet if r["veto_off"]["same"]]
            rest = [r for r in rs if r not in vet]
            print(
                f"   MECHANISM: a hand was vetoed in {len(vet)}/{len(rs)}; of those veto_off makes E(mirror)==E(x) in {len(fixed)}/{len(vet)}"
            )
            print(
                f"   NOT this mechanism: {len(rest)}/{len(rs)}  (veto_off same: {sum(1 for r in rest if r['veto_off']['same'])}, "
                f"fold_off same: {sum(1 for r in rest if r['fold_off']['same'])})"
            )
    sys.stdout.flush()  # os._exit skips the interpreter's buffers; a piped run printed NOTHING
    os._exit(0)


def _collect(ex, mols, mode):
    futs = {m: ex.submit(probe, m, mode) for m in mols}
    res = {}
    for m, f in futs.items():
        try:
            res[m] = f.result(timeout=1800)
        except FT:
            res[m] = None
    print("oinsmiles from:", next((r["src"] for r in res.values() if r), "?"), flush=True)
    return res


def run_renumber(args, out):
    recs = [r for r in esc._load(out / "e_selfconsistency.jsonl") if r["base"]]
    order = ["TIMEOUT", "encode_fail", "structural", "facmer_divergent"]

    def worst(r):
        rels = [r["rel"][t] for t in ("renum0", "renum1", "renum2")]
        if any(v in order for v in rels):
            return "key_changes"
        return "slot_level" if any(v != "byte_exact" for v in rels) else "stable"

    pop = {"slot_level": [], "key_changes": [], "stable": []}
    for r in recs:
        pop[worst(r)].append(r["molecule"])
    rng = random.Random(42)
    groups = {
        "renumber drifts at SLOT level (THE QUESTION)": sorted(pop["slot_level"]),
        "renumber changes the KEY (control: not a labeling problem)": sorted(
            rng.sample(pop["key_changes"], min(args.n_control, len(pop["key_changes"])))
        ),
        "renumber stable (control)": sorted(rng.sample(pop["stable"], args.n_control)),
    }
    print("populations:", {k: len(v) for k, v in groups.items()}, flush=True)
    allm = sorted({m for ms in groups.values() for m in ms})
    res = _collect(ProcessPoolExecutor(args.cpu), allm, "renumber")
    with open(out / "veto_probe_renumber.jsonl", "w") as fh:
        for name, ms in groups.items():
            for m in ms:
                if res[m]:
                    fh.write(json.dumps({"group": name, **res[m]}) + "\n")
        fh.write(f"#DONE {len(allm)}\n")
    for name, ms in groups.items():
        rs = [res[m] for m in ms if res[m]]
        print(f"\n== {name}   n={len(ms)}  probed={len(rs)}")
        for arm in ARMS:
            c = Counter(
                {True: "all 3 same", False: "drifts", None: "encode_fail"}[r[arm]["same"]]
                for r in rs
            )
            print(f"   {arm:9s} 3 renumberings vs base: {dict(c)}")
        c = Counter(
            "vetoed on some presentation"
            if "vetoed_collapse" in [r["default"]["outcome_x"], *r["default"]["outcome_m"]]
            else "never vetoed"
            for r in rs
        )
        print("   default-arm:", dict(c))
        if "QUESTION" in name:
            vet = [
                r
                for r in rs
                if "vetoed_collapse" in [r["default"]["outcome_x"], *r["default"]["outcome_m"]]
            ]
            print(
                f"   MECHANISM: vetoed on some presentation {len(vet)}/{len(rs)}; of those veto_off makes all 3 equal in {sum(1 for r in vet if r['veto_off']['same'])}/{len(vet)}"
            )
            rest = [r for r in rs if r not in vet]
            print(
                f"   NOT this mechanism: {len(rest)}/{len(rs)} (veto_off all-same: {sum(1 for r in rest if r['veto_off']['same'])}, fold_off all-same: {sum(1 for r in rest if r['fold_off']['same'])})"
            )
    sys.stdout.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
