"""v0.4.17 L1b -- read the two harness arms and say what OIN_EXACT_DONOR_FOLD does to the verdict.

Reads only. Each arm directory (``ab_off``, ``ab_on``) must already hold

    bucket_report_honest.json   tools/roundtrip_bucket_report.py --results-dir <arm> --score honest
    g_verdict.jsonl             tools/census/g_vs_input.py --cohort <movers> --sweep <arm> --out <arm>

and the script prints those two commands and stops if either is missing.

TWO NUMBERS, because the census showed one is not enough. ``byte_exact`` is a SELF-AGREEMENT rate:
``E(x) == E(G(E(x)))``. A canonicality lever can raise it by teaching the second encode to agree
with the first, whatever was built. So every transition is also scored VERIFIED, with the census's
own predicate (``attribution_table.py`` fault ``NONE``): the pass must also have a string that
describes the input (parse-back rules clean -- a slot relabeling cannot change these, so the
census verdict is carried) and a GENERATED STRUCTURE the neutral ruler finds graph-isomorphic to
the input with no stereo element mirrored or changed -- recomputed per arm from that arm's own
structures.

Printed in the order it should be believed:

1. COMPLETENESS and provenance per arm (report count, tree the harness loaded, xtb).
2. DEAD-LEVER CHECK: molecules whose ``smiles_1`` differs between the arms. Zero means the lever
   never reached the ON arm and nothing below means anything.
3. NOISE FLOOR: the OFF arm against the sweep of record on the same molecules. The generator is
   seeded but its budget is advisory, so load moves verdicts; gains and losses smaller than this
   count are not distinguishable from it.
4. The transition table, gains, losses -- self-consistent and verified.
5. The projection onto the 5,000-molecule headline. Non-movers are unchanged BY CONSTRUCTION
   (same input string, same re-encode), so the projection is exact up to (3).

    <main>/.venv/bin/python tools/v0417/generator_ab_report.py
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
ISO = ("ISO", "ISO_MARGINAL", "ISO_CLASH")
STRING_FAULTS = ("P_E1_COVERAGE", "DATA_MULTI", "P_DETACHED", "E1_GRAPH", "E1_HCOUNT")
BAD_STEREO = ("MIRROR", "MIRROR_PARTIAL", "DIFFERENT")
N_COHORT, PASS_SELF, PASS_VERIFIED = 5000, 3858, 3462  # sweep of record; census C4


def _jsonl(path):
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not rows or not rows[-1].startswith("#DONE"):
        sys.exit(f"ABORT: {path} has no #DONE trailer")
    return {r["molecule"]: r for r in map(json.loads, rows[:-1])}


def _arm(d, cohort):
    need = {
        "bucket_report_honest.json": f"tools/roundtrip_bucket_report.py --results-dir {d} --score honest",
        "g_verdict.jsonl": f"tools/census/g_vs_input.py --cohort {cohort} --sweep {d} --out {d}",
    }
    missing = [cmd for f, cmd in need.items() if not (d / f).exists()]
    if missing:
        sys.exit(
            "ABORT: run first (PYTHONPATH=$PWD/src <main>/.venv/bin/python ...):\n  "
            + "\n  ".join(missing)
        )
    B = {
        r["molecule"].removesuffix(".xyz"): r
        for r in json.loads((d / "bucket_report_honest.json").read_text())
    }
    return B, _jsonl(d / "g_verdict.jsonl")


def _verified(m, B, G, T):
    if B[m]["bucket"] != "byte_exact" or T[m]["fault"] in STRING_FAULTS:
        return False
    g = G.get(m)
    return bool(g) and g.get("graph") in ISO and g.get("stereo") not in BAD_STEREO


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--ab", type=Path, default=MAIN / "results-v0.4.17-exactfold")
    ap.add_argument("--cohort", type=Path, default=MAIN / "cohort-v0.4.17-exactfold-movers")
    ap.add_argument("--record", type=Path, default=MAIN / "results-v0.4.14-sweep")
    ap.add_argument("--table", type=Path, default=MAIN / "results-census/attribution_table.tsv")
    args = ap.parse_args()

    cohort = sorted(p.stem for p in args.cohort.glob("*.xyz"))
    T = {r["molecule"]: r for r in csv.DictReader(open(args.table), delimiter="\t")}
    arms = {}
    print(f"1. COMPLETENESS   cohort={len(cohort)}")
    for a in ("off", "on"):
        d = args.ab / f"ab_{a}"
        B, G = _arm(d, args.cohort)
        arms[a] = (B, G)
        reps = list((d / "individual_reports").glob("*.json"))
        one = json.loads(reps[0].read_text()) if reps else {}
        tree = {
            m.group(1)
            for p in d.glob("shard*.log")
            for m in [re.search(r"oinsmiles loaded from: (\S+)", p.read_text()[:4000])]
            if m
        }
        print(
            f"   {a:3s} reports={len(reps)}  bucketed={len(B)}  ruler verdicts={len(G)}  "
            f"xtb_available={one.get('xtb_available')}  commit={one.get('commit_id')}  tree={sorted(tree)}"
        )
    (Bf, Gf), (Bn, Gn) = arms["off"], arms["on"]
    both = [m for m in cohort if m in Bf and m in Bn]
    print(f"   scored in BOTH arms: {len(both)}/{len(cohort)}")

    moved_in = [m for m in both if Bf[m].get("smiles_1") != Bn[m].get("smiles_1")]
    print(
        f"\n2. DEAD-LEVER CHECK   smiles_1 differs between arms: {len(moved_in)}  (0 = lever never fired)"
    )

    R = {
        r["molecule"].removesuffix(".xyz"): r
        for r in json.loads((args.record / "bucket_report_honest.json").read_text())
    }
    flip = [
        m for m in both if (R[m]["bucket"] == "byte_exact") != (Bf[m]["bucket"] == "byte_exact")
    ]
    up = sum(1 for m in flip if Bf[m]["bucket"] == "byte_exact")
    print(
        f"\n3. NOISE FLOOR   OFF arm vs sweep of record, pass/fail flips: {len(flip)}/{len(both)}  "
        f"(now pass {up}, now fail {len(flip) - up})"
    )

    print("\n4. TRANSITIONS  OFF -> ON")
    for label, fn in (
        ("self-consistent (byte_exact)", lambda m, B, G: B[m]["bucket"] == "byte_exact"),
        ("VERIFIED (ruler on the structure)", lambda m, B, G: _verified(m, B, G, T)),
    ):
        off = {m: fn(m, Bf, Gf) for m in both}
        on = {m: fn(m, Bn, Gn) for m in both}
        gains = sorted(m for m in both if not off[m] and on[m])
        losses = sorted(m for m in both if off[m] and not on[m])
        print(f"   {label}")
        print(f"      pass OFF={sum(off.values())}  pass ON={sum(on.values())}  of {len(both)}")
        print(
            f"      GAINS {len(gains)}   LOSSES {len(losses)}   NET {len(gains) - len(losses):+d}"
        )
        mv = set(moved_in)
        print(
            f"      gains: input moved {sum(1 for m in gains if m in mv)} / generated side only "
            f"{sum(1 for m in gains if m not in mv)};  losses: input moved "
            f"{sum(1 for m in losses if m in mv)} / generated side only {sum(1 for m in losses if m not in mv)}"
        )
        print(
            "      gains by census fault :",
            dict(Counter(T[m]["fault"] for m in gains).most_common()),
        )
        print(
            "      losses by census fault:",
            dict(Counter(T[m]["fault"] for m in losses).most_common()),
        )
        print(
            "      loss transitions      :",
            dict(Counter(Bn[m]["bucket"] for m in losses).most_common()),
        )
        base = PASS_SELF if label.startswith("self") else PASS_VERIFIED
        net = len(gains) - len(losses)
        print(
            f"      PROJECTION on n={N_COHORT}: {100 * base / N_COHORT:.2f}% -> "
            f"{100 * (base + net) / N_COHORT:.2f}%  ({net / (N_COHORT / 100):+.2f} pts; noise floor "
            f"+/-{len(flip) / (N_COHORT / 100):.2f})"
        )
        if label.startswith("self"):
            (args.ab / "ab_losses.txt").write_text("\n".join(losses) + "\n")
            (args.ab / "ab_gains.txt").write_text("\n".join(gains) + "\n")
    tab = Counter(
        (Bf[m]["bucket"], Bn[m]["bucket"]) for m in both if Bf[m]["bucket"] != Bn[m]["bucket"]
    )
    print("\n   bucket transitions (off -> on):")
    for (x, y), c in tab.most_common():
        print(f"      {c:4d}  {x} -> {y}")
    el = {a: sum(float(arms[a][0][m].get("elapsed_s") or 0) for m in both) / 3600 for a in arms}
    print(f"\n5. RUNTIME  sum elapsed_s  OFF={el['off']:.2f} h  ON={el['on']:.2f} h")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
