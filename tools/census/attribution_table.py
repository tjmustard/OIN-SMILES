"""I4 of the census: ONE table, 5,000 rows, one fault per molecule -- passes included.

    python tools/census/attribution_table.py            # writes attribution_table.tsv + summary
    python tools/census/attribution_table.py --gz       # also gzips the per-molecule census files
                                                        # for /freeze-measurements (census_*.jsonl.gz)

Joins, per molecule: the honest bucket (v0.4.14 sweep), the I1 ruler verdict on the generated
structure (C1, ``g_verdict.jsonl``) and on the mirrored input (chirality, ``g_verdict_control-
mirror.jsonl``), the encoder's self-consistency (C2, ``e_selfconsistency.jsonl``), the parse-back
and perception flags (C3, ``parseback.jsonl``, ``parseback_gen.jsonl``, ``pflags.jsonl``,
``collisions.json``), and the attach class (``attach_class_audit.json``).

THE FAULT IS THE FIRST RULE THAT FIRES -- in this order, which is the order of the pipeline
--------------------------------------------------------------------------------------------
 1  no string at all                                   P_E1_COVERAGE
 2  the string does not describe the input (C3):
      a fragment with no slot, unbound in the input    DATA_MULTI      (two molecules in one file)
      a fragment with no slot, bound in the input      P_DETACHED      (perception dropped a ligand)
      adapter-level graph != input                     E1_GRAPH
      adapter-level H count != input                   E1_HCOUNT
 3  no generated structure                             G_NOTHING       (sub: TIMEOUT / NOCONF / ...)
 4  generated graph != input (C1)                      G_CONSTRUCTION  (sub: DETACHED / ...)
 5  graph =, chirality MIRROR, E separates the hands   G_HANDEDNESS
    graph =, chirality MIRROR, E(mirror x) == E(x)     E1_NONINJECTIVE (G was never told)
    graph =, chirality DIFFERENT (a diastereomer)      G_DIASTEREOMER  (sub: ARRANGEMENT / CENTRES)
    either of the two G_ verdicts on a PASS             E1_NOT_ENCODED  (the string never carried
                                                                        the element G changed)
 6  graph = and chirality =, and the round trip FAILED:
      H decoration differs on the generated structure  G_HCOUNT
      achiral by the ruler and E(mirror x) != E(x)     E1_NONCANONICAL (one molecule, two strings)
      string changes under renumbering / 0.02 A noise  E2_P_FRAGILE
      else                                             E2_P_OTHER      (bond orders / labels on a
                                                                        correct structure)
 7  PASSED and no rule above fired                     NONE
 8  FAILED and no rule above fired                     UNATTRIBUTED    (reported, never folded)

A PASS that carries a fault from rules 2, 4 or 5 is a METRIC FALSE PASS: the round trip agreed
with itself about a wrong structure or a wrong string. Every flag is also kept as its own column,
so the overlaps the first-match rule hides can be read back.

WHAT A BROKEN VERSION WOULD PRINT
---------------------------------
A join keyed on the wrong name field silently yields 5,000 rows of ``UNATTRIBUTED``; a rule that
never fires yields a plausible table with one class missing. So: the bucket column must reproduce
``bucket_report_honest`` exactly, the classes must sum to 1,142 FAIL + 3,858 PASS, each rule's
count is printed, the DETACHED cross-check against the attach-class audit is printed, and
``UNATTRIBUTED`` is printed as a number.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
DEF_SWEEP = MAIN / "results-v0.4.14-sweep"
DEF_OUT = MAIN / "results-census"

#: absolute path of a checkout on this box -> ``<checkout:name>`` in anything that leaves it
_LOCAL_PATH = re.compile(r"/home/[^/\"]+/Documents/GitHub/([^/\"]+)")

ISO = ("ISO", "ISO_MARGINAL", "ISO_CLASH")
STEREO_SAME = ("SAME", "NO_STEREO", "SAME_PARTIAL")
#: What the bucket column must reproduce, PER SWEEP, and how many rows must come out ``NONE``.
#: Neither number is computed here: the buckets are each sweep's frozen ``bucket_report_honest``
#: and ``none`` is the VERIFIED count ``tools/v0417/sweep_two_numbers.py`` printed for that sweep
#: from a six-line re-implementation of rules 1-5. Two tools, one answer, or the run aborts.
EXPECTED = {
    "results-v0.4.14-sweep": {
        "buckets": {
            "byte_exact": 3858,
            "structural": 484,
            "key_equal": 365,
            "hard_fail": 266,
            "facmer_divergent": 15,
            "encode_fail": 12,
        },
        "none": 3462,
    },
    "results-v0.4.17-sweep": {
        "buckets": {
            "byte_exact": 4136,
            "structural": 479,
            "key_equal": 122,
            "hard_fail": 236,
            "facmer_divergent": 16,
            "encode_fail": 11,
        },
        "none": 3730,
    },
}
EXPECTED_BUCKETS = EXPECTED["results-v0.4.14-sweep"]["buckets"]
PT_PER_MOL = 100.0 / 5000

COLUMNS = [
    "molecule", "bucket", "subclass", "outcome", "fault", "sub", "rule",
    "metal", "geo", "cn_in", "eta", "n_atoms", "elapsed_s", "charge_stated", "os_in_allowed",
    "n_radical", "n_charged_c", "attach_class", "cluster",
    "g_graph", "g_stereo", "g_sphere", "g_hdeco", "ruler_chiral",
    "e2_base_ok", "e2_mirror_differs", "e2_fragile", "e2_key_fragile",
    "pb_parser_graph", "pb_adapter_graph", "pb_adapter_h", "pb_unslotted", "pb_str_valence_fail",
    "gen_parser_graph", "gen_unslotted", "probe_equal", "error",
]  # fmt: skip


def _jl(path):
    # A missing instrument used to read as "no rows", and a join over no rows is a plausible table
    # with one class quietly absent. Every input is required; say which one is not there.
    if not path.exists():
        sys.exit(f"ABORT: instrument file missing: {path}")
    return [json.loads(ln) for ln in path.read_text().splitlines() if not ln.startswith("#")]


def _by_mol(recs):
    return {r["molecule"]: r for r in recs}


def _g_nothing_sub(err: str) -> str:
    e = err or ""
    if "TimeoutException" in e and "generating" in e:
        return "TIMEOUT"
    if "TimeoutException" in e:
        return "TIMEOUT_OTHER"
    if "Atom count mismatch" in e:
        return "ATOMCOUNT"
    if "Generation/Verification failed" in e or "child process died" in e:
        return "NOCONF"
    return "BLANK" if not e else "OTHER"


def attribute(row, g1, gm, e2, pb, pg, pf, attach, cluster) -> dict:
    """One molecule -> one TSV row; the fault is the first rule that fires."""
    m = row["molecule"]
    passed = row["bucket"] == "byte_exact"
    rel = (e2 or {}).get("rel", {})
    e2_base_ok = bool((e2 or {}).get("base"))
    e2_mirror_differs = e2_base_ok and rel.get("mirror") not in (None, "byte_exact")
    frag_keys = ("renum0", "renum1", "renum2", "noise0", "noise1", "noise2")
    e2_fragile = e2_base_ok and any(rel.get(k) not in (None, "byte_exact") for k in frag_keys)
    e2_key_fragile = e2_base_ok and any(
        rel.get(k)
        not in (None, "byte_exact", "key_equal/slot_renumber", "key_equal/rdkit_canonical")
        for k in frag_keys
    )
    ruler_chiral = bool((gm or {}).get("input_chiral_by_ruler"))
    g_graph = (g1 or {}).get("graph")
    g_stereo = (g1 or {}).get("stereo")
    g_sphere = (g1 or {}).get("sphere")
    g_hdeco = (g1 or {}).get("h_decoration")
    pb_parser = (pb or {}).get("parser_graph")
    pb_adapter = (pb or {}).get("adapter_graph")
    pb_h = (pb or {}).get("adapter_h_decoration")
    pb_unslotted = int((pb or {}).get("n_unslotted") or 0)
    gen_unslotted = int((pg or {}).get("n_unslotted") or 0)

    fault, sub, rule = None, "", 0
    if not row.get("smiles_1"):
        fault, sub, rule = "P_E1_COVERAGE", _g_nothing_sub(row.get("error")), 1
    elif pb_unslotted:
        rule = 2
        fault = "DATA_MULTI" if pb_parser in ISO else "P_DETACHED"
        sub = f"components in={pb.get('n_components_in')} str={pb.get('n_components_str')}"
    elif pb_adapter is not None and pb_adapter not in ISO:
        fault, sub, rule = "E1_GRAPH", pb_adapter, 2
    elif pb_h == "DIFF":
        fault, sub, rule = "E1_HCOUNT", f"h_diff={pb.get('adapter_h_diff_atoms')}", 2
    elif g1 is None:
        fault, sub, rule = "G_NOTHING", _g_nothing_sub(row.get("error")), 3
    elif g_graph not in ISO:
        rule = 4
        fault = "G_CONSTRUCTION"
        sub = "DETACHED" if (gen_unslotted or attach == "DETACHED") else g_graph
    elif g_stereo in ("MIRROR", "MIRROR_PARTIAL"):
        rule = 5
        fault = "G_HANDEDNESS" if e2_mirror_differs else "E1_NONINJECTIVE"
        sub = g_stereo
    elif g_stereo == "DIFFERENT":
        rule = 5
        fault = "G_DIASTEREOMER"
        sub = "ARRANGEMENT" if g_sphere == "ARRANGEMENT_DIFF" else "CENTRES"
    if rule == 5 and passed and fault in ("G_HANDEDNESS", "G_DIASTEREOMER"):
        # the generator changed a stereo element and the round trip still agreed byte for byte:
        # the string never carried that element (ligand sp3 centres are not enforced by design)
        fault, sub = "E1_NOT_ENCODED", f"{g_stereo}/{sub}"
    if rule:
        pass  # rules 1-5 fired above
    elif not passed:
        rule = 6
        if g_hdeco == "DIFF":
            fault, sub = "G_HCOUNT", ""
        elif e2_mirror_differs and not ruler_chiral:
            fault, sub = "E1_NONCANONICAL", "veto one-hand (C2)"
        elif e2_fragile:
            fault, sub = "E2_P_FRAGILE", "key" if e2_key_fragile else "string"
        else:
            fault, sub = "E2_P_OTHER", g_stereo or ""
        if fault == "UNATTRIBUTED":
            rule = 8
    else:
        fault, sub, rule = "NONE", "", 7
    return {
        "molecule": m,
        "bucket": row["bucket"],
        "subclass": row.get("subclass") or "",
        "outcome": "PASS" if passed else "FAIL",
        "fault": fault,
        "sub": sub,
        "rule": rule,
        "metal": ((pf or {}).get("metal") or [""])[0],
        "geo": (pb or {}).get("geo", ""),
        "cn_in": ((pb or {}).get("cn_in") or [""])[0],
        "eta": int(bool(row.get("eta"))),
        "n_atoms": (pf or {}).get("n_atoms", (pb or {}).get("n_heavy_in", "")),
        "elapsed_s": row.get("elapsed_s", ""),
        "charge_stated": (pf or {}).get("charge_stated", ""),
        "os_in_allowed": (pf or {}).get("os_in_allowed", ""),
        "n_radical": (pf or {}).get("n_radical_atoms", ""),
        "n_charged_c": (pf or {}).get("n_charged_carbon", ""),
        "attach_class": attach or "",
        "cluster": cluster or "",
        "g_graph": g_graph or "",
        "g_stereo": g_stereo or "",
        "g_sphere": g_sphere or "",
        "g_hdeco": g_hdeco or "",
        "ruler_chiral": int(ruler_chiral),
        "e2_base_ok": int(e2_base_ok),
        "e2_mirror_differs": int(e2_mirror_differs),
        "e2_fragile": int(e2_fragile),
        "e2_key_fragile": int(e2_key_fragile),
        "pb_parser_graph": pb_parser or "",
        "pb_adapter_graph": pb_adapter or "",
        "pb_adapter_h": pb_h or "",
        "pb_unslotted": pb_unslotted,
        "pb_str_valence_fail": (pb or {}).get("str_valence_fail", ""),
        "gen_parser_graph": (pg or {}).get("parser_graph", ""),
        "gen_unslotted": gen_unslotted,
        "probe_equal": (pf or {}).get("probe_equal", ""),
        "error": (row.get("error") or "")[:60],
    }


def build(args):
    """Join every instrument, attribute, verify, write."""
    if args.sweep.name not in EXPECTED:
        sys.exit(
            f"ABORT: no EXPECTED entry for {args.sweep.name}. Add its frozen bucket counts and its "
            "VERIFIED count -- from that sweep's own reports, never from this tool's output."
        )
    expect = EXPECTED[args.sweep.name]
    src = {
        "g_verdict": args.g_verdict or args.out / "g_verdict.jsonl",
        "mirror": args.mirror or args.out / "g_verdict_control-mirror.jsonl",
        "e2": args.e2 or args.out / "e_selfconsistency.jsonl",
        "parseback": args.parseback or args.out / "parseback.jsonl",
        "parseback_gen": args.parseback_gen or args.out / "parseback_gen.jsonl",
        "pflags": args.pflags or args.out / "pflags.jsonl",
        "attach": args.attach or args.sweep / "attach_class_audit.json",
        "collisions": args.collisions or args.out / "collisions.json",
    }
    print(f"sweep: {args.sweep.name}")
    for k, v in src.items():
        shown = _LOCAL_PATH.sub(r"<checkout:\1>", str(v))  # this print is frozen: no home dir
        print(f"  {k:14s} <- {shown}")
    rows = json.loads((args.sweep / "bucket_report_honest.json").read_text())
    g1 = _by_mol(_jl(src["g_verdict"]))
    gm = _by_mol(_jl(src["mirror"]))
    e2 = _by_mol(_jl(src["e2"]))
    pb = _by_mol(_jl(src["parseback"]))
    pg = _by_mol(_jl(src["parseback_gen"]))
    pf = _by_mol(_jl(src["pflags"]))
    if not src["attach"].exists() or not src["collisions"].exists():
        sys.exit(f"ABORT: instrument file missing: {src['attach']} or {src['collisions']}")
    aca = json.loads(src["attach"].read_text())
    attach = {m: cls for t in aca["table"].values() for cls, ms in t.items() for m in ms}
    col = json.loads(src["collisions"].read_text())
    cluster = {m: c["cluster"] for c in col["clusters"] for m in c["members"]}
    print(
        f"inputs: rows={len(rows)} g1={len(g1)} mirror={len(gm)} e2={len(e2)} parseback={len(pb)} "
        f"parseback_gen={len(pg)} pflags={len(pf)} attach={len(attach)} clustered={len(cluster)}"
    )
    table = [
        attribute(
            r, g1.get(r["molecule"]), gm.get(r["molecule"]), e2.get(r["molecule"]),
            pb.get(r["molecule"]), pg.get(r["molecule"]), pf.get(r["molecule"]),
            attach.get(r["molecule"]), cluster.get(r["molecule"]),
        )
        for r in rows
    ]  # fmt: skip
    table.sort(key=lambda t: t["molecule"])

    # ---- verification, before anything is written ------------------------------------------
    n = len(table)
    bc = Counter(t["bucket"] for t in table)
    oc = Counter(t["outcome"] for t in table)
    print(f"\nDENOMINATOR rows={n}")
    print(
        "bucket column:",
        dict(bc),
        " reproduces bucket_report_honest:",
        bc == Counter(expect["buckets"]),
    )
    n_pass = expect["buckets"]["byte_exact"]
    print("outcome:", dict(oc), f" (expected PASS {n_pass} / FAIL {5000 - n_pass})")
    assert n == 5000 and bc == Counter(expect["buckets"]), (
        "the join does not reproduce the bucket table"
    )
    fc = Counter(t["fault"] for t in table)
    print(
        f"NONE = {fc.get('NONE', 0)}   VERIFIED by sweep_two_numbers.py = {expect['none']}   "
        f"{'AGREE' if fc.get('NONE', 0) == expect['none'] else 'DISAGREE'}"
    )
    if fc.get("NONE", 0) != expect["none"]:
        sys.exit("ABORT: this table and sweep_two_numbers.py disagree about how many passes verify")
    print("fault classes:", dict(fc.most_common()))
    print("rule fired:", dict(sorted(Counter(t["rule"] for t in table).items())))
    print(
        "UNATTRIBUTED:",
        fc.get("UNATTRIBUTED", 0),
        " UNDECIDED-stereo rows:",
        sum(1 for t in table if t["g_stereo"].startswith("UNDECIDED")),
    )

    # ---- cross-checks the plan asked for -----------------------------------------------------
    det = {m for m, c in attach.items() if c == "DETACHED"}
    r4 = {t["molecule"] for t in table if t["fault"] == "G_CONSTRUCTION"}
    print(
        f"\nrule 4 (G_CONSTRUCTION) n={len(r4)}  ∩ audit DETACHED (all buckets {len(det)}) = {len(r4 & det)}"
    )
    fails_det = {m for m in det if m in {t["molecule"] for t in table if t["outcome"] == "FAIL"}}
    print(
        f"  DETACHED among FAIL = {len(fails_det)}; of those G_CONSTRUCTION={len(fails_det & r4)}, "
        f"other faults={dict(Counter(t['fault'] for t in table if t['molecule'] in fails_det - r4).most_common(6))}"
    )
    sr = [t for t in table if t["subclass"] == "slot_renumber"]
    print(
        f"slot_renumber n={len(sr)} by fault:", dict(Counter(t["fault"] for t in sr).most_common())
    )

    # ---- crosstabs ---------------------------------------------------------------------------
    xb = defaultdict(Counter)
    for t in table:
        xb[t["fault"]][t["bucket"]] += 1
    order = sorted(xb, key=lambda f: -sum(xb[f].values()))
    print("\nfault x bucket (n, pts)")
    buckets = [
        "byte_exact",
        "structural",
        "key_equal",
        "hard_fail",
        "facmer_divergent",
        "encode_fail",
    ]
    print(f"  {'fault':16s} {'n':>5s} {'pts':>6s}  " + "  ".join(f"{b[:10]:>10s}" for b in buckets))
    for f in order:
        tot = sum(xb[f].values())
        print(
            f"  {f:16s} {tot:5d} {tot * PT_PER_MOL:6.2f}  "
            + "  ".join(f"{xb[f].get(b, 0):10d}" for b in buckets)
        )
    fails = [t for t in table if t["outcome"] == "FAIL"]
    passes = [t for t in table if t["outcome"] == "PASS"]
    print(f"\nFAIL partition (n={len(fails)}, {len(fails) * PT_PER_MOL:.2f} pts):")
    for f, c in Counter(t["fault"] for t in fails).most_common():
        subs = Counter(t["sub"] for t in fails if t["fault"] == f).most_common(4)
        print(f"  {f:16s} {c:5d} {c * PT_PER_MOL:6.2f}   {subs}")
    fp = [t for t in passes if t["fault"] != "NONE"]
    print(
        f"\nMETRIC FALSE PASS (PASS with a fault): {len(fp)} = {100 * len(fp) / len(passes):.1f}% of passes"
    )
    for f, c in Counter(t["fault"] for t in fp).most_common():
        print(f"  {f:16s} {c:5d}")
    by_owner = Counter()
    for t in fails:
        by_owner[t["fault"].split("_")[0]] += 1
    print(
        "\nFAIL by owner prefix:",
        {k: (v, round(v * PT_PER_MOL, 2)) for k, v in by_owner.most_common()},
    )

    # ---- write -------------------------------------------------------------------------------
    dest = args.write_to or args.out
    dest.mkdir(parents=True, exist_ok=True)
    tsv = dest / "attribution_table.tsv"
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        for t in table:
            w.writerow(t)
    summary = {
        "n_rows": n,
        "bucket": dict(bc),
        "outcome": dict(oc),
        "fault": dict(fc),
        "fault_x_bucket": {f: dict(c) for f, c in xb.items()},
        "fail_partition": {f: c for f, c in Counter(t["fault"] for t in fails).items()},
        "fail_partition_pts": {
            f: round(c * PT_PER_MOL, 2) for f, c in Counter(t["fault"] for t in fails).items()
        },
        "fail_sub": {f: dict(Counter(t["sub"] for t in fails if t["fault"] == f)) for f in fc},
        "metric_false_pass": dict(Counter(t["fault"] for t in fp)),
        "metric_false_pass_n": len(fp),
        "fail_by_owner": dict(by_owner),
        "rule_fired": dict(Counter(t["rule"] for t in table)),
        "unattributed": fc.get("UNATTRIBUTED", 0),
        "crosscheck": {
            "G_CONSTRUCTION_n": len(r4),
            "audit_DETACHED_all": len(det),
            "audit_DETACHED_among_fail": len(fails_det),
            "G_CONSTRUCTION_and_DETACHED_fail": len(fails_det & r4),
            "slot_renumber_by_fault": dict(Counter(t["fault"] for t in sr)),
        },
        "pts_per_molecule": PT_PER_MOL,
    }
    (dest / "attribution_summary.json").write_text(json.dumps(summary, indent=1))
    print(f"\nwrote {tsv} ({n} rows) and attribution_summary.json")
    if args.gz:
        with open(tsv, "rb") as fi, gzip.open(args.out / "attribution_table.tsv.gz", "wb") as fo:
            shutil.copyfileobj(fi, fo)
        for name in (
            "g_verdict", "g_verdict_control-mirror", "e_selfconsistency", "parseback",
            "parseback_gen", "pflags", "veto_probe", "veto_probe_chiral", "veto_probe_renumber",
        ):  # fmt: skip
            src = args.out / f"{name}.jsonl"
            if not src.exists():
                continue
            # measurements/ is PUBLIC. The veto probes record which oinsmiles they imported (a
            # `src` field with an absolute path: the proof against the sys.path trap). Keep the
            # proof, drop the home directory.
            with open(src) as fi, gzip.open(args.out / f"census_{name}.jsonl.gz", "wt") as fo:
                for ln in fi:
                    fo.write(_LOCAL_PATH.sub(r"<checkout:\1>", ln))
        print("gzipped per-molecule census files as census_*.jsonl.gz + attribution_table.tsv.gz")


def main():
    """CLI."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--sweep", type=Path, default=DEF_SWEEP)
    ap.add_argument("--out", type=Path, default=DEF_OUT)
    ap.add_argument("--gz", action="store_true", help="also write the gzipped freeze copies")
    ap.add_argument("--write-to", type=Path, help="write the table here instead of into --out")
    for flag, what in (
        ("--g-verdict", "C1 ruler verdicts on THIS sweep's generated structures"),
        ("--mirror", "C1 ruler on the mirrored INPUT (input-only: reusable across sweeps)"),
        ("--e2", "C2 encoder self-consistency, under the encoder THIS sweep shipped"),
        ("--parseback", "C3 parse-back of THIS sweep's smiles_1"),
        ("--parseback-gen", "C3 parse-back of THIS sweep's smiles_2 / generated structures"),
        ("--pflags", "C3 perception flags (reads THIS sweep's smiles_1)"),
        ("--attach", "attach_class_audit.json of THIS sweep"),
        ("--collisions", "natural-twin clusters (input-only: reusable across sweeps)"),
    ):
        ap.add_argument(flag, type=Path, help=f"{what} (default: the census layout)")
    args = ap.parse_args()
    for p in (args.sweep, args.out):
        if not p.exists():
            sys.exit(f"missing: {p}")
    build(args)


if __name__ == "__main__":
    main()
