"""v0.4.18 -- the promotion sweep against the A/B that predicted it, and against the old baseline.

Reads only. The A/B (``eta_ab_report.py``) projected 86.04% self-consistent / 77.60% VERIFIED onto
the 5,000 on two claims this script CHECKS rather than repeats:

  1. a NON-eta molecule cannot reach either lever, so it is unchanged BY CONSTRUCTION
     -> every non-eta structure must be byte-identical to the sweep this one replaces, and
  2. generation is deterministic at this load (the A/B's OFF arm == that sweep, 1,063/1,063)
     -> every eta structure must be byte-identical to the A/B's ON arm's.

Parameterised for the release sweep (--sweep, --record, --ab-on, --predicted); the defaults are
the L2 promotion's (results-v0.4.18-sweep against results-v0.4.17-sweep and ab_on).

A broken promotion prints differently in each half: levers that never reached the shipped default
leave the eta half identical to the OLD sweep instead of to ``ab_on`` (counted and printed); an
encoder that moved shows up as ``smiles_1`` differences (must be 0).

    <main>/.venv/bin/python tools/v0418/sweep_vs_ab.py [--out report.json]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
ISO = ("ISO", "ISO_MARGINAL", "ISO_CLASH")
STRING_FAULTS = ("P_E1_COVERAGE", "DATA_MULTI", "P_DETACHED", "E1_GRAPH", "E1_HCOUNT")
BAD_STEREO = ("MIRROR", "MIRROR_PARTIAL", "DIFFERENT")
PREDICTED = "4136+166,3730+150"  # the A/B's projection: base passes + net, self then verified


def _sum(expr: str) -> int:
    """'4305+45' -> 4350; no eval."""
    return sum(int(t) for t in expr.replace("-", "+-").split("+") if t)


def _jsonl(path):
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not rows or not rows[-1].startswith("#DONE"):
        sys.exit(f"ABORT: {path} has no #DONE trailer")
    return {r["molecule"]: r for r in map(json.loads, rows[:-1])}


def _buckets(d):
    return {
        r["molecule"].removesuffix(".xyz"): r
        for r in json.loads((d / "bucket_report_honest.json").read_text())
    }


def _sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--sweep", type=Path, default=MAIN / "results-v0.4.18-sweep")
    ap.add_argument("--record", type=Path, default=MAIN / "results-v0.4.17-sweep")
    ap.add_argument("--ab-on", type=Path, default=MAIN / "results-v0.4.18-eta-detached/ab_on")
    ap.add_argument("--eta-cohort", type=Path, default=MAIN / "cohort-v0.4.18-eta")
    ap.add_argument(
        "--table", type=Path, default=MAIN / "results-v0.4.17-reattribution/attribution_table.tsv"
    )
    ap.add_argument(
        "--predicted",
        default=PREDICTED,
        help="'<self base>+<net>,<verified base>+<net>' -- what the A/B projected (default: L2's)",
    )
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    predicted = dict(zip(("self", "verified"), (_sum(x) for x in args.predicted.split(","))))

    T = {r["molecule"]: r for r in csv.DictReader(open(args.table), delimiter="\t")}
    B, R, A = _buckets(args.sweep), _buckets(args.record), _buckets(args.ab_on)
    G, GR = _jsonl(args.sweep / "g_verdict.jsonl"), _jsonl(args.record / "g_verdict.jsonl")
    eta = {p.stem for p in args.eta_cohort.glob("*.xyz")}
    if len(B) != 5000 or set(B) != set(R) or len(eta) != 1146 or not eta <= set(B):
        sys.exit(f"ABORT: denominators -- sweep {len(B)}, record {len(R)}, eta {len(eta)}")
    print(f"DENOMINATOR 5000 = {len(eta)} eta-bound + {5000 - len(eta)} non-eta")

    def s1(X, m):
        v = X[m].get("smiles_1")
        return v if v and v != "None" else None

    moved = [m for m in B if s1(B, m) and s1(R, m) and s1(B, m) != s1(R, m)]
    one_sided = [m for m in B if bool(s1(B, m)) != bool(s1(R, m))]
    print(
        f"\n1. INPUT SIDE   smiles_1 differs from the v0.4.17 sweep: {len(moved)} (MUST be 0)   "
        f"encoded in one sweep only: {len(one_sided)}"
    )

    def struct(d, m):
        return _sha(d / "structures" / f"{m}_generated.xyz")

    rep = {"smiles_1_moved": len(moved)}
    print("\n2. THE TWO CLAIMS THE PROJECTION STOOD ON")
    for name, ms, ref_dir, ref_B, label in (
        ("non-eta", sorted(set(B) - eta), args.record, R, "v0.4.17 sweep"),
        ("eta", sorted(eta), args.ab_on.parent / "ab_on", A, "A/B ON arm"),
    ):
        cmp_ = [(m, struct(args.sweep, m), struct(ref_dir, m)) for m in ms]
        both = [(m, a, b) for m, a, b in cmp_ if a and b]
        same = sum(a == b for _, a, b in both)
        only = sum(1 for _, a, b in cmp_ if bool(a) != bool(b))
        bsame = sum(B[m]["bucket"] == ref_B[m]["bucket"] for m in ms)
        print(
            f"   {name:8s} n={len(ms):4d}  structure == {label}: {same}/{len(both)}  (built in one "
            f"only: {only})   bucket ==: {bsame}/{len(ms)}"
        )
        rep[name] = {
            "n": len(ms),
            "struct_same": same,
            "struct_cmp": len(both),
            "bucket_same": bsame,
        }
    old = [(m, struct(args.sweep, m), struct(args.record, m)) for m in sorted(eta)]
    n_old = sum(1 for _, a, b in old if a and b and a == b)
    print(
        f"   eta structures still identical to the OLD sweep: {n_old}  (all of them = dead levers)"
    )

    def verified(m, X, GX):
        if X[m]["bucket"] != "byte_exact" or T[m]["fault"] in STRING_FAULTS:
            return False
        g = GX.get(m)
        return bool(g) and g.get("graph") in ISO and g.get("stereo") not in BAD_STEREO

    print("\n3. THE TWO NUMBERS   (VERIFIED with the string faults CARRIED from the census table;")
    print("   two_numbers.txt re-derives them from this sweep's own parse-back -- they must agree)")
    for key, fn in (
        ("self", lambda m, X, GX: X[m]["bucket"] == "byte_exact"),
        ("verified", verified),
    ):
        new = {m: fn(m, B, G) for m in B}
        was = {m: fn(m, R, GR) for m in B}
        g = [m for m in B if new[m] and not was[m]]
        l_ = [m for m in B if was[m] and not new[m]]
        n = sum(new.values())
        print(
            f"   {key:9s} {sum(was.values())} -> {n} = {n / 50:.2f}%   predicted {predicted[key]} "
            f"({predicted[key] / 50:.2f}%)   off by {n - predicted[key]:+d}"
        )
        print(
            f"             gains {len(g)} (eta {sum(m in eta for m in g)}, non-eta "
            f"{sum(m not in eta for m in g)})   losses {len(l_)} (eta {sum(m in eta for m in l_)}, "
            f"non-eta {sum(m not in eta for m in l_)})"
        )
        rep[key] = {"was": sum(was.values()), "now": n, "gains": len(g), "losses": len(l_)}
    print("\n4. BUCKETS  ", dict(Counter(r["bucket"] for r in B.values()).most_common()))
    print("   was      ", dict(Counter(r["bucket"] for r in R.values()).most_common()))
    el = [float(r.get("elapsed_s") or 0) for r in B.values()]
    elr = [float(r.get("elapsed_s") or 0) for r in R.values()]
    print(
        f"\n5. RUNTIME  sum elapsed_s {sum(elr) / 3600:.2f} h -> {sum(el) / 3600:.2f} h   >30 s "
        f"{sum(x > 30 for x in elr)} -> {sum(x > 30 for x in el)}   max {max(elr):.0f} -> {max(el):.0f} s"
    )
    rep["over_30s"] = sum(x > 30 for x in el)
    if args.out:
        args.out.write_text(json.dumps(rep, indent=1, sort_keys=True))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
