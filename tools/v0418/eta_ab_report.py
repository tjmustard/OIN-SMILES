"""v0.4.18 L2 -- read the two harness arms and say what the eta construction levers do. Reads only.

    OFF   OIN_ETA_COVALENT_TARGET=0  OIN_VDW_EXEMPT_BINDING=0
    ON    OIN_ETA_COVALENT_TARGET=1  OIN_VDW_EXEMPT_BINDING=1

Each arm directory (``ab_off``, ``ab_on``; ``tools/v0418/post_eta_ab.sh`` makes both files) holds

    bucket_report_honest.json   tools/roundtrip_bucket_report.py --results-dir <arm> --score honest
    g_verdict.jsonl             tools/census/g_vs_input.py --cohort <eta cohort> --sweep <arm> --out <arm>

TWO NUMBERS. ``byte_exact`` is SELF-AGREEMENT, ``E(x) == E(G(E(x)))``. VERIFIED is the census's
fault ``NONE``: the pass must also carry a string that describes the input (a property of the INPUT
string, identical in both arms and in the sweep of record -- carried from the re-attribution table)
and a GENERATED STRUCTURE the neutral ruler finds graph-isomorphic to the input with no stereo
element mirrored or changed -- recomputed per arm from that arm's own structures.

THIS IS A GENERATOR-SIDE LEVER, so v0.4.17's dead-lever check (``smiles_1`` differs between arms)
is INVERTED here: the input strings must be IDENTICAL -- a difference is contamination, and the
report aborts on it. The lever is shown to fire by counting GENERATED STRUCTURES that differ.

Printed in the order it should be believed:

1. COMPLETENESS and provenance per arm.
2. INPUT CONTROL (must be 0) and DEAD-LEVER CHECK (must not be 0).
3. NOISE FLOOR: the OFF arm against the sweep of record on the same molecules -- pass/fail flips
   and byte-identical structures. The generator is seeded but its budget is advisory, so load
   moves verdicts; gains and losses below this count are not distinguishable from it.
4. Transitions, gains, losses -- self-consistent and verified -- by census fault, by metal.
5. The structures themselves: coordination intact (the harness's own ``coordination`` block) and
   the eta M-C inflation (``eta_distance_audit.compare``), per arm, by census class.
6. Runtime, and the projection onto the 5,000-molecule headline. A non-eta molecule is unchanged
   BY CONSTRUCTION (both levers are read only inside a complex that has a multi-atom binding
   group), so the projection is exact up to (3).

    <main>/.venv/bin/python tools/v0418/eta_ab_report.py [--out report.json]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import statistics as st
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "src"))  # beats PYTHONPATH, which is the point
sys.path.insert(0, str(HERE.parent))

import eta_distance_audit as audit  # noqa: E402

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
ISO = ("ISO", "ISO_MARGINAL", "ISO_CLASH")
STRING_FAULTS = ("P_E1_COVERAGE", "DATA_MULTI", "P_DETACHED", "E1_GRAPH", "E1_HCOUNT")
BAD_STEREO = ("MIRROR", "MIRROR_PARTIAL", "DIFFERENT")
#: results-v0.4.17-sweep, the baseline of record: self-consistent / VERIFIED of 5,000
N_COHORT, PASS_SELF, PASS_VERIFIED = 5000, 4136, 3730


def _jsonl(path):
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not rows or not rows[-1].startswith("#DONE"):
        sys.exit(f"ABORT: {path} has no #DONE trailer")
    return {r["molecule"]: r for r in map(json.loads, rows[:-1])}


def _buckets(d: Path):
    return {
        r["molecule"].removesuffix(".xyz"): r
        for r in json.loads((d / "bucket_report_honest.json").read_text())
    }


def _arm(d: Path, cohort: Path):
    need = {
        "bucket_report_honest.json": f"tools/roundtrip_bucket_report.py --results-dir {d} --score honest",
        "g_verdict.jsonl": f"tools/census/g_vs_input.py --cohort {cohort} --sweep {d} --out {d}",
    }
    missing = [cmd for f, cmd in need.items() if not (d / f).exists()]
    if missing:
        sys.exit("ABORT: run tools/v0418/post_eta_ab.sh first:\n  " + "\n  ".join(missing))
    return _buckets(d), _jsonl(d / "g_verdict.jsonl")


def _sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def _verified(m, B, G, T):
    if B[m]["bucket"] != "byte_exact" or T[m]["fault"] in STRING_FAULTS:
        return False
    g = G.get(m)
    return bool(g) and g.get("graph") in ISO and g.get("stereo") not in BAD_STEREO


def _cls(row):
    if row["fault"] == "G_CONSTRUCTION" and row["sub"] == "DETACHED":
        return "DETACHED"
    if row["fault"] == "NONE":
        return "VERIFIED_PASS"
    return "other_FAIL" if row["outcome"] == "FAIL" else "other_false_PASS"


def _median(xs):
    return f"{st.median(xs):+.3f}" if xs else "  n/a"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--ab", type=Path, default=MAIN / "results-v0.4.18-eta-detached")
    ap.add_argument("--cohort", type=Path, default=MAIN / "cohort-v0.4.18-eta")
    ap.add_argument("--record", type=Path, default=MAIN / "results-v0.4.17-sweep")
    ap.add_argument(
        "--table",
        type=Path,
        default=MAIN / "results-v0.4.17-reattribution/attribution_table.tsv",
    )
    ap.add_argument("--out", type=Path)
    ap.add_argument(
        "--smoke",
        action="store_true",
        help="tool debugging only: print past the two control aborts, every line after them tagged",
    )
    args = ap.parse_args()

    def _abort(msg):
        if not args.smoke:
            sys.exit("ABORT: " + msg)
        print(f"   !!! SMOKE MODE, CONTROL FAILED ({msg}) -- NOTHING BELOW IS A RESULT !!!")

    cohort = sorted(p.stem for p in args.cohort.glob("*.xyz"))
    if not cohort:
        sys.exit(f"ABORT: {args.cohort} is empty -- refusing an empty denominator")
    T = {r["molecule"]: r for r in csv.DictReader(open(args.table), delimiter="\t")}
    report = {"n_cohort": len(cohort)}

    # ---- 1. completeness -------------------------------------------------------------------
    arms, reps = {}, {}
    print(f"1. COMPLETENESS   cohort={len(cohort)}   {(args.ab / 'AB_COMMIT').read_text().strip()}")
    for a in ("off", "on"):
        d = args.ab / f"ab_{a}"
        arms[a] = _arm(d, args.cohort)
        reps[a] = {
            p.stem: json.loads(p.read_text()) for p in (d / "individual_reports").glob("*.json")
        }
        one = next(iter(reps[a].values()), {})
        tree = {
            m.group(1)
            for p in d.glob("shard*.log")
            for m in [re.search(r"oinsmiles loaded from: (\S+)", p.read_text()[:4000])]
            if m
        }
        print(
            f"   {a:3s} reports={len(reps[a])}  bucketed={len(arms[a][0])}  ruler verdicts="
            f"{len(arms[a][1])}  xtb_available={one.get('xtb_available')}  "
            f"commit={one.get('commit_id')}  tree={sorted(tree)}"
        )
    (Bf, Gf), (Bn, Gn) = arms["off"], arms["on"]
    both = [m for m in cohort if m in Bf and m in Bn]
    print(f"   scored in BOTH arms: {len(both)}/{len(cohort)}")
    if len(both) != len(cohort):
        sys.exit("ABORT: an arm is incomplete -- nothing below is printed")

    # ---- 2. input control + dead-lever check -----------------------------------------------
    R = _buckets(args.record)

    def _s1(B, m):
        v = B[m].get("smiles_1")
        return v if v and v != "None" else None

    # a string on one side and NONE on the other is an encode that ran out of budget under load
    # (timing, counted); two DIFFERENT strings is an encoder that is not the same encoder (fatal)
    in_moved = [m for m in both if _s1(Bf, m) and _s1(Bn, m) and _s1(Bf, m) != _s1(Bn, m)]
    in_vs_rec = [m for m in both if _s1(Bf, m) and _s1(R, m) and _s1(Bf, m) != _s1(R, m)]
    one_sided = [m for m in both if bool(_s1(Bf, m)) != bool(_s1(Bn, m))]
    print(
        f"\n2. INPUT CONTROL   smiles_1 differs OFF vs ON: {len(in_moved)}   OFF vs sweep of record: "
        f"{len(in_vs_rec)}   (both MUST be 0: the levers are generator-side)   encoded in one arm "
        f"only: {len(one_sided)} (encode budget under load)"
    )
    if in_moved or in_vs_rec:
        _abort(f"the input side moved ({(in_moved + in_vs_rec)[:5]}) -- contamination")
    sha = {
        a: {m: _sha(args.ab / f"ab_{a}" / "structures" / f"{m}_generated.xyz") for m in both}
        for a in ("off", "on")
    }
    built_both = [m for m in both if sha["off"][m] and sha["on"][m]]
    differ = [m for m in built_both if sha["off"][m] != sha["on"][m]]
    only_on = [m for m in both if sha["on"][m] and not sha["off"][m]]
    only_off = [m for m in both if sha["off"][m] and not sha["on"][m]]
    print(
        f"   DEAD-LEVER CHECK   structure in both arms {len(built_both)}, of which DIFFERENT "
        f"{len(differ)}  (0 = the levers never fired);  built only ON {len(only_on)}, only OFF "
        f"{len(only_off)}"
    )
    if not differ:
        _abort("no generated structure differs between the arms -- the levers never fired")
    report.update(
        n_structures_differ=len(differ), built_only_on=len(only_on), built_only_off=len(only_off)
    )

    # ---- 3. noise floor --------------------------------------------------------------------
    flip = [
        m for m in both if (R[m]["bucket"] == "byte_exact") != (Bf[m]["bucket"] == "byte_exact")
    ]
    up = sum(1 for m in flip if Bf[m]["bucket"] == "byte_exact")
    rec_sha = {m: _sha(args.record / "structures" / f"{m}_generated.xyz") for m in both}
    n_cmp = [m for m in both if rec_sha[m] and sha["off"][m]]
    n_same = sum(1 for m in n_cmp if rec_sha[m] == sha["off"][m])
    print(
        f"\n3. NOISE FLOOR   OFF arm vs sweep of record: pass/fail flips {len(flip)}/{len(both)} "
        f"(now pass {up}, now fail {len(flip) - up});  structures byte-identical {n_same}/{len(n_cmp)}"
    )
    report.update(noise_flips=len(flip), off_vs_record_identical=[n_same, len(n_cmp)])

    # ---- 4. transitions --------------------------------------------------------------------
    print("\n4. TRANSITIONS  OFF -> ON")
    for label, key, fn, base in (
        (
            "self-consistent (byte_exact)",
            "self",
            lambda m, B, G: B[m]["bucket"] == "byte_exact",
            PASS_SELF,
        ),
        (
            "VERIFIED (ruler on the structure)",
            "verified",
            lambda m, B, G: _verified(m, B, G, T),
            PASS_VERIFIED,
        ),
    ):
        off = {m: fn(m, Bf, Gf) for m in both}
        on = {m: fn(m, Bn, Gn) for m in both}
        gains = sorted(m for m in both if not off[m] and on[m])
        losses = sorted(m for m in both if off[m] and not on[m])
        net = len(gains) - len(losses)
        print(f"   {label}")
        print(f"      pass OFF={sum(off.values())}  pass ON={sum(on.values())}  of {len(both)}")
        print(f"      GAINS {len(gains)}   LOSSES {len(losses)}   NET {net:+d}")
        for tag, ms in (("gains ", gains), ("losses", losses)):
            print(
                f"      {tag} by census fault:",
                dict(
                    Counter(f"{T[m]['fault']}/{T[m]['sub']}".rstrip("/") for m in ms).most_common(8)
                ),
            )
            print(
                f"      {tag} by metal       :",
                dict(Counter(T[m]["metal"] for m in ms).most_common(12)),
            )
        print(
            "      loss transitions      :",
            dict(Counter(f"{Bf[m]['bucket']}->{Bn[m]['bucket']}" for m in losses).most_common()),
        )
        print(
            f"      PROJECTION on n={N_COHORT}: {100 * base / N_COHORT:.2f}% -> "
            f"{100 * (base + net) / N_COHORT:.2f}%  ({net / (N_COHORT / 100):+.2f} pts; noise floor "
            f"+/-{len(flip) / (N_COHORT / 100):.2f})"
        )
        (args.ab / f"ab_gains_{key}.txt").write_text("\n".join(gains) + "\n")
        (args.ab / f"ab_losses_{key}.txt").write_text("\n".join(losses) + "\n")
        report[key] = {
            "pass_off": sum(off.values()),
            "pass_on": sum(on.values()),
            "gains": len(gains),
            "losses": len(losses),
            "projection_pct": round(100 * (base + net) / N_COHORT, 2),
        }
    tab = Counter(
        (Bf[m]["bucket"], Bn[m]["bucket"]) for m in both if Bf[m]["bucket"] != Bn[m]["bucket"]
    )
    print("\n   bucket transitions (off -> on):")
    for (x, y), c in tab.most_common():
        print(f"      {c:4d}  {x} -> {y}")

    # ---- 5. the structures -----------------------------------------------------------------
    print("\n5. THE STRUCTURES   by census class (the v0.4.17 sweep's attribution)")
    print(
        f"   {'class':17s} {'n':>5s} | {'intact OFF':>10s} {'intact ON':>10s} | "
        f"{'eta infl OFF':>12s} {'eta infl ON':>12s} | {'in band OFF':>11s} {'in band ON':>11s}"
    )
    intact = {
        a: {m: bool((reps[a][m].get("coordination") or {}).get("intact")) for m in both}
        for a in arms
    }
    infl = {a: defaultdict(list) for a in arms}
    band = {a: defaultdict(lambda: [0, 0]) for a in arms}
    unmatched = Counter()
    for a in arms:
        for m in both:
            g = args.ab / f"ab_{a}" / "structures" / f"{m}_generated.xyz"
            if not g.exists():
                continue
            cmp_ = audit.compare((args.cohort / f"{m}.xyz").read_text(), g.read_text())
            unmatched[a] += cmp_.get("unmatched_ligands", 0)
            for grp in cmp_.get("groups", []):
                if grp["kind"] == "eta":
                    infl[a][_cls(T[m])].append(grp["inflation"])
                    band[a][_cls(T[m])][0] += grp["in_band"]
                    band[a][_cls(T[m])][1] += grp["k"]
    by_cls = Counter(_cls(T[m]) for m in both)
    for c, n in by_cls.most_common():
        ms = [m for m in both if _cls(T[m]) == c]
        pct = lambda a: (  # noqa: E731
            f"{100 * band[a][c][0] / band[a][c][1]:5.1f}%" if band[a][c][1] else "  n/a"
        )
        print(
            f"   {c:17s} {n:5d} | {sum(intact['off'][m] for m in ms):10d} "
            f"{sum(intact['on'][m] for m in ms):10d} | {_median(infl['off'][c]):>12s} "
            f"{_median(infl['on'][c]):>12s} | {pct('off'):>11s} {pct('on'):>11s}"
        )
    print(
        f"   all               {len(both):5d} | {sum(intact['off'].values()):10d} "
        f"{sum(intact['on'].values()):10d}   (unmatched ligands, counted not dropped: "
        f"OFF {unmatched['off']}, ON {unmatched['on']})"
    )
    report["intact"] = {a: sum(intact[a].values()) for a in arms}

    # ---- 6. runtime ------------------------------------------------------------------------
    el = {a: [float(arms[a][0][m].get("elapsed_s") or 0) for m in both] for a in arms}
    print(
        f"\n6. RUNTIME  sum elapsed_s  OFF={sum(el['off']) / 3600:.2f} h  ON={sum(el['on']) / 3600:.2f} h"
        f"   >30 s: OFF {sum(x > 30 for x in el['off'])}  ON {sum(x > 30 for x in el['on'])}"
        f"   max: OFF {max(el['off']):.0f} s  ON {max(el['on']):.0f} s"
    )
    report["elapsed_h"] = {a: round(sum(el[a]) / 3600, 2) for a in arms}
    report["over_30s"] = {a: sum(x > 30 for x in el[a]) for a in arms}
    if args.out:
        args.out.write_text(json.dumps(report, indent=1, sort_keys=True))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
