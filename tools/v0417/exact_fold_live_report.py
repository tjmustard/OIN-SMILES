"""The exact fold, LIVE: C2's self-consistency run under ``OIN_EXACT_DONOR_FOLD=1`` vs shipped.

Both files come from ``tools/census/e_selfconsistency.py`` with the same seed, so the eleven
presentations of every molecule are the same files in both arms and every column below is
comparable molecule by molecule. This script only reads; it encodes nothing.

Printed, in the order they should be believed:

1. DENOMINATORS and instrument health for both arms (a missing ``#DONE`` aborts).
2. OFFLINE == LIVE: the audit's predicted ``E_exact(x)`` / ``E_exact(mirror x)`` against the
   strings the live encoder emitted. This is what licenses reading the offline audit at all --
   and what a lever that never reached the encoder would fail (it would match the SHIPPED column).
3. The 2x2 (ruler chirality x ``E(mirror) == E(x)``), shipped vs exact.
4. Renumber drift at slot level vs key level, shipped vs exact. The key level is perception and
   must NOT move; if it does, the lever is touching something it has no business touching.
5. The MOVERS -- molecules whose ``E(x)`` changed -- by census bucket and fault. This is the
   population a generator A/B owes, derived from coordinates, not from stored strings.

    <main>/.venv/bin/python tools/v0417/exact_fold_live_report.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
RENUM = ("renum0", "renum1", "renum2")
STABLE = ("again", "rewrite", "rot0")
KEY_LEVEL = ("structural", "facmer_divergent", "encode_fail", "TIMEOUT")


def _load(path):
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not rows or not rows[-1].startswith("#DONE"):
        sys.exit(f"ABORT: {path} has no #DONE trailer -- the run is not finished")
    return {r["molecule"]: r for r in map(json.loads, rows[:-1])}


def _variant(rec, tag):
    """The string a presentation encoded to: ``diff`` holds it only when it differs from base."""
    rel = rec["rel"].get(tag)
    if rel == "byte_exact":
        return rec["base"]
    return rec["diff"].get(tag) if rel not in ("TIMEOUT", "encode_fail", "MISSING") else None


def _renum_class(rec):
    rels = [rec["rel"].get(t) for t in RENUM]
    if any(r in KEY_LEVEL for r in rels):
        return "key_level"
    return "slot_level" if any(r != "byte_exact" for r in rels) else "stable"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--shipped", type=Path, default=MAIN / "results-census/e_selfconsistency.jsonl")
    ap.add_argument(
        "--exact",
        type=Path,
        default=MAIN / "results-v0.4.17-exactfold/e_selfconsistency_exact.jsonl",
    )
    ap.add_argument(
        "--audit", type=Path, default=MAIN / "results-v0.4.17-exactfold/autofold_audit.json"
    )
    ap.add_argument(
        "--ruler", type=Path, default=MAIN / "results-census/g_verdict_control-mirror.jsonl"
    )
    ap.add_argument("--table", type=Path, default=MAIN / "results-census/attribution_table.tsv")
    ap.add_argument("--movers-out", type=Path, default=None)
    args = ap.parse_args()

    S, X = _load(args.shipped), _load(args.exact)
    print(f"1. DENOMINATORS  shipped n={len(S)}  exact n={len(X)}")
    for name, d in (("shipped", S), ("exact", X)):
        enc = [r for r in d.values() if r["base"]]
        rels = Counter(v for r in enc for v in r["rel"].values())
        print(
            f"   {name:8s} base encoded={len(enc)}  TIMEOUT={rels['TIMEOUT']}  "
            f"encode_fail={rels['encode_fail']}  INSTRUMENT={rels['INSTRUMENT']}  "
            f"MISSING={rels['MISSING']}  CPU-h={sum(r['t_total'] for r in d.values()) / 3600:.2f}"
        )
        bad = sum(1 for r in enc if any(r["rel"].get(t) != "byte_exact" for t in STABLE))
        print(f"   {name:8s} again/rewrite/rot0 not byte-stable: {bad}  (determinism floor)")
    both = sorted(m for m in S if S[m]["base"] and X.get(m, {}).get("base"))
    print(f"   encoded in BOTH arms: {len(both)}")

    print("\n2. OFFLINE AUDIT == LIVE ENCODER")
    audit = json.loads(args.audit.read_text())
    seen, agree = {}, Counter()
    for group, rows in audit.items():
        if "PROTECTED" in group:
            continue
        for w in rows:
            seen[w["molecule"]] = w
    for m, w in seen.items():
        if m not in X or not X[m]["base"]:
            agree["not_encoded_live"] += 1
            continue
        agree[("x", X[m]["base"] == w["auto_x"])] += 1
        lm = _variant(X[m], "mirror")
        if lm is not None:
            agree[("m", lm == w["auto_m"])] += 1
    print(f"   audited molecules: {len(seen)}")
    print(
        f"   E_exact(x)      offline == live: {agree[('x', True)]}/{agree[('x', True)] + agree[('x', False)]}"
    )
    print(
        f"   E_exact(mirror) offline == live: {agree[('m', True)]}/{agree[('m', True)] + agree[('m', False)]}"
    )
    shipped_like = sum(
        1 for m, w in seen.items() if m in both and X[m]["base"] == S[m]["base"] != w["auto_x"]
    )
    print(
        f"   live == SHIPPED where the audit predicted a change (a dead lever prints >0): {shipped_like}"
    )

    print("\n3. THE 2x2  (ruler chirality x E(mirror)==E(x))")
    chiral = {
        r["molecule"]: r.get("input_chiral_by_ruler")
        for r in map(json.loads, (ln for ln in args.ruler.read_text().splitlines() if ln[0] == "{"))
    }
    for name, d in (("shipped", S), ("exact", X)):
        tab = Counter()
        for m in both:
            rel = d[m]["rel"].get("mirror")
            if rel in ("TIMEOUT", "encode_fail", "MISSING") or chiral.get(m) is None:
                tab["unusable"] += 1
                continue
            tab[("chiral" if chiral[m] else "achiral", rel == "byte_exact")] += 1
        print(
            f"   {name:8s} achiral&differs (NON-CANONICAL)={tab[('achiral', False)]:4d}  "
            f"achiral&same={tab[('achiral', True)]:4d}  chiral&same (NON-INJECTIVE)={tab[('chiral', True)]:4d}  "
            f"chiral&differs={tab[('chiral', False)]:4d}  unusable={tab['unusable']}"
        )

    print("\n4. RENUMBER DRIFT (three renumberings vs base)")
    for name, d in (("shipped", S), ("exact", X)):
        c = Counter(_renum_class(d[m]) for m in both)
        print(
            f"   {name:8s} stable={c['stable']}  slot_level={c['slot_level']}  key_level={c['key_level']}"
        )
    moved_key = [
        m
        for m in both
        if (_renum_class(S[m]) == "key_level") != (_renum_class(X[m]) == "key_level")
    ]
    print(f"   key-level membership changed between arms (perception; expect ~0): {len(moved_key)}")

    print("\n5. MOVERS  E_exact(x) != E_shipped(x)")
    T = {r["molecule"]: r for r in csv.DictReader(open(args.table), delimiter="\t")}
    movers = [m for m in both if S[m]["base"] != X[m]["base"]]
    print(f"   movers: {len(movers)}/{len(both)}")
    print("   by bucket:", dict(Counter(T[m]["bucket"] for m in movers).most_common()))
    print("   by fault :", dict(Counter(T[m]["fault"] for m in movers).most_common()))
    e1 = [m for m in both if T[m]["fault"] == "E1_NONCANONICAL"]
    fixed = [m for m in e1 if X[m]["rel"].get("mirror") == "byte_exact"]
    print(f"   E1_NONCANONICAL failures now mirror-stable: {len(fixed)}/{len(e1)}")
    dest = args.movers_out or (args.exact.parent / "exact_fold_movers.txt")
    dest.write_text("\n".join(movers) + "\n")
    print("   wrote", dest)
    sys.stdout.flush()


if __name__ == "__main__":
    main()
