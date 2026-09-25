"""v0.4.19 serializer lane -- what the two string-side levers do to the cohort, in the order the
evidence should be believed. Reads only.

    OFF   shipped defaults (OIN_H_FAITHFUL, OIN_RC1_PROPAGATE both off) = results-v0.4.18-release-sweep
    ON    OIN_H_FAITHFUL=1  OIN_RC1_PROPAGATE=1

An ENCODER lever moves ``smiles_1``. That has three consequences the v0.4.18 generator-side report
cannot express, and each is a subcommand here:

  changed   The complete CHANGED-STRING SET, from the full-cohort ``e_selfconsistency`` run under the
            levers (base column) against the sweep of record's ``smiles_1``. This is the population a
            live arm must run; every other molecule's generated structure is unchanged by
            construction (generation is byte-deterministic at this load, v0.4.18) and is re-scored
            offline. Also attributes each changed string to a lever when the single-lever encodes
            exist (tools/v0419/lever_parseback.py --arms shipped,hfaith,rc1prop over the set).
  audit     CANONICALITY, whole cohort: determinism / rewrite / rotation / 3 renumberings / 3 noise /
            mirror, lever run vs shipped run, per transform -- and the molecules that became stable
            or unstable. A canonicality lever must not make a string depend on atom order.
  offline   The UNCHANGED-string rows: the sweep's generated xyz re-encoded under the levers
            (tools/honest_rescore.py) against the lever's own base string. A pass that stops agreeing
            is a LOSS the lever causes with no generator involved; the reverse is a gain.
  ab        The live arm over the changed set vs the sweep's own rows, scored self-consistent
            (``byte_exact``, honest) AND VERIFIED -- the census's fault ``NONE``, computed per arm by
            the census's own ``attribute()`` with that arm's bucket, ruler verdict, parse-back and
            self-consistency records. CONTROL: on the OFF arm the recomputed fault must equal the
            release census table's fault, row for row.
  all       everything above, then the projection onto the 5,000 headline.

    <main>/.venv/bin/python tools/v0419/serializer_ab_report.py all [--out report.json]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1] / "census"))

from attribution_table import attribute  # noqa: E402

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
SWEEP = MAIN / "results-v0.4.18-release-sweep"
CENSUS = MAIN / "results-v0.4.18-release-census"
LANE = MAIN / "results-v0.4.19-serializer"
SHIPPED_ESC = MAIN / "results-v0.4.17-exactfold" / "e_selfconsistency_exact.jsonl"
RULER_MIRROR = MAIN / "results-v0.4.17-exactfold" / "g_verdict_control-mirror.jsonl"
N_COHORT, PASS_SELF, PASS_VERIFIED = 5000, 4347, 3920
TRANSFORMS = (
    "again",
    "rewrite",
    "rot0",
    "renum0",
    "renum1",
    "renum2",
    "noise0",
    "noise1",
    "noise2",
    "mirror",
)


def _jsonl(path: Path, need_done=True):
    rows = [ln for ln in path.read_text().splitlines() if ln.strip()]
    if need_done and (not rows or not rows[-1].startswith("#DONE")):
        sys.exit(f"ABORT: {path} has no #DONE trailer -- still running or died")
    recs = [json.loads(r) for r in rows if r.startswith("{")]
    return {r["molecule"]: r for r in recs}


def _buckets(d: Path):
    p = d / "bucket_report_honest.json"
    if not p.exists():
        sys.exit(
            f"ABORT: {p} missing -- run tools/roundtrip_bucket_report.py --score honest on {d}"
        )
    return {r["molecule"].removesuffix(".xyz"): r for r in json.loads(p.read_text())}


def _table():
    with open(CENSUS / "attribution_table.tsv") as fh:
        return {r["molecule"]: r for r in csv.DictReader(fh, delimiter="\t")}


def _reports(d: Path, mols):
    out = {}
    for m in mols:
        p = d / "individual_reports" / f"{m}.json"
        out[m] = json.loads(p.read_text()) if p.exists() else None
    return out


# ----------------------------------------------------------------------------- changed
def cmd_changed(args, res):
    esc = _jsonl(args.esc)
    T = _table()
    if len(esc) != N_COHORT:
        sys.exit(f"ABORT: lever self-consistency run has {len(esc)} molecules, not {N_COHORT}")
    S = _buckets(SWEEP)
    changed, enc_fail, same = [], [], 0
    for m, r in esc.items():
        s1 = (S.get(m) or {}).get("smiles_1")
        if r.get("base_error") or not r.get("base"):
            enc_fail.append((m, r.get("base_error"), bool(s1)))
        elif r["base"] != (s1 or ""):
            changed.append(m)
        else:
            same += 1
    print(f"== CHANGED-STRING SET  ({args.esc.name} vs {SWEEP.name} smiles_1)")
    print(
        f"   unchanged {same}   changed {len(changed)}   encode failed under the levers {len(enc_fail)}"
    )
    lost = [x for x in enc_fail if x[2]]
    if lost:
        print(
            f"   🔴 {len(lost)} molecules the shipped encoder encodes and the lever run does NOT:"
        )
        for m, e, _ in lost[:20]:
            print(f"      {m}: {e}")
    no_string_both = [x for x in enc_fail if not x[2]]
    print(f"   no string in either ({len(no_string_both)}): P_E1_COVERAGE, unchanged")
    print("   changed, by census fault:")
    for (f, o), n in sorted(
        Counter((T[m]["fault"], T[m]["outcome"]) for m in changed).items(), key=lambda x: -x[1]
    ):
        print(f"      {f:18s} {o:4s} {n:4d}")
    print(
        "   changed, eta:",
        dict(Counter("eta" if T[m]["eta"] == "1" else "non-eta" for m in changed)),
    )
    # per-lever attribution when the single-lever encodes exist
    if args.single and args.single.exists():
        by = defaultdict(dict)
        for ln in args.single.read_text().splitlines():
            if ln.startswith("{"):
                r = json.loads(ln)
                by[r["molecule"]][r["arm"]] = r.get("oin")
        attr = Counter()
        for m in changed:
            a = by.get(m, {})
            h = a.get("hfaith") is not None and a.get("hfaith") != a.get("shipped")
            p = a.get("rc1prop") is not None and a.get("rc1prop") != a.get("shipped")
            attr[
                "H_FAITHFUL only"
                if h and not p
                else "RC1_PROPAGATE only"
                if p and not h
                else "both"
                if h and p
                else "neither(!)"
            ] += 1
        print("   attribution (single-lever encodes):", dict(attr))
    res["changed"] = {
        "n_changed": len(changed),
        "n_unchanged": same,
        "encode_failed": enc_fail,
        "molecules": changed,
    }
    if args.write_set:
        args.write_set.write_text(
            "# v0.4.19 changed-string set: base under OIN_H_FAITHFUL=1 OIN_RC1_PROPAGATE=1 != release-sweep smiles_1\n"
            + "\n".join(changed)
            + "\n"
        )
        print(f"   wrote {args.write_set} ({len(changed)} molecules)")
    return set(changed)


# ----------------------------------------------------------------------------- audit
def cmd_audit(args, res):
    lever, ship = _jsonl(args.esc), _jsonl(SHIPPED_ESC)
    ruler = _jsonl(RULER_MIRROR, need_done=False)
    T = _table()
    print(
        f"\n== CANONICALITY AUDIT  {len(lever)} molecules, lever run vs shipped run ({SHIPPED_ESC.name})"
    )
    print(
        f"   {'transform':10s} {'shipped byte-stable':>20s} {'lever byte-stable':>18s} {'became unstable':>16s} {'became stable':>14s}"
    )
    detail = {}
    for t in TRANSFORMS:
        s_ok = {m for m, r in ship.items() if r.get("base") and r["rel"].get(t) == "byte_exact"}
        l_ok = {m for m, r in lever.items() if r.get("base") and r["rel"].get(t) == "byte_exact"}
        both = set(ship) & set(lever)
        worse = sorted((both & s_ok) - l_ok)
        better = sorted((both & l_ok) - s_ok)
        detail[t] = {"shipped": len(s_ok), "lever": len(l_ok), "worse": worse, "better": better}
        print(f"   {t:10s} {len(s_ok):20d} {len(l_ok):18d} {len(worse):16d} {len(better):14d}")

    # the 2x2 against the ruler: achiral & mirror differs = NON-CANONICAL; chiral & same = NON-INJECTIVE
    def two_by_two(E):
        c = Counter()
        for m, r in E.items():
            if not r.get("base") or m not in ruler:
                continue
            chiral = bool(ruler[m].get("input_chiral_by_ruler"))
            differs = r["rel"].get("mirror") not in (None, "byte_exact")
            c[("chiral" if chiral else "achiral", "differs" if differs else "same")] += 1
        return c

    s2, l2 = two_by_two(ship), two_by_two(lever)
    print("   mirror 2x2 (ruler chirality x string):   shipped  ->  lever")
    for k in (
        ("achiral", "same"),
        ("achiral", "differs"),
        ("chiral", "differs"),
        ("chiral", "same"),
    ):
        tag = {
            ("achiral", "differs"): "  NON-CANONICAL",
            ("chiral", "same"): "  NON-INJECTIVE",
        }.get(k, "")
        print(f"      {k[0]:8s} {k[1]:8s} {s2.get(k, 0):6d}  ->  {l2.get(k, 0):6d}{tag}")

    def fragile(r):
        return bool(r.get("base")) and any(
            r["rel"].get(k) not in (None, "byte_exact")
            for k in ("renum0", "renum1", "renum2", "noise0", "noise1", "noise2")
        )

    sf = {m for m, r in ship.items() if fragile(r)}
    lf = {m for m, r in lever.items() if fragile(r)}
    print(
        f"   E2-fragile (renumber/noise moves the string): shipped {len(sf)} -> lever {len(lf)}   became fragile {len(lf - sf)}   became stable {len(sf - lf)}"
    )
    newly = sorted(lf - sf)
    if newly:
        print(
            "   🔴 became fragile under the levers:",
            ", ".join(f"{m}({T[m]['fault']})" for m in newly[:25]),
        )
    res["audit"] = {
        "per_transform": detail,
        "mirror_2x2": {
            "shipped": {f"{a}/{b}": n for (a, b), n in s2.items()},
            "lever": {f"{a}/{b}": n for (a, b), n in l2.items()},
        },
        "fragile": {
            "shipped": len(sf),
            "lever": len(lf),
            "became_fragile": newly,
            "became_stable": sorted(sf - lf),
        },
    }


# ----------------------------------------------------------------------------- offline
def cmd_offline(args, res, changed):
    """Unchanged-string rows: the string is the sweep's, the structure is the sweep's (deterministic),
    so the only things the levers can move are the honest re-encode of that structure and the
    self-consistency columns. Recompute the census fault with exactly those two swapped in."""
    esc_l, esc_s = _jsonl(args.esc), _jsonl(SHIPPED_ESC)
    S = _buckets(SWEEP)
    rs = args.rescore / "honest_rescore.jsonl"
    if not rs.exists():
        sys.exit(f"ABORT: {rs} missing")
    R = {}
    for ln in rs.read_text().splitlines():
        if ln.startswith("{"):
            r = json.loads(ln)
            R[r["molecule"]] = r
    T = _table()
    G = _jsonl(SWEEP / "g_verdict.jsonl")
    PB = _jsonl(SWEEP / "parseback.jsonl")
    PG = _jsonl(CENSUS / "parseback_gen.jsonl") if (CENSUS / "parseback_gen.jsonl").exists() else {}
    GM = _jsonl(RULER_MIRROR, need_done=False)
    PF = (
        _jsonl(CENSUS / "pflags.jsonl", need_done=False)
        if (CENSUS / "pflags.jsonl").exists()
        else {}
    )
    print(
        f"\n== OFFLINE RE-SCORE of the UNCHANGED-string rows ({args.rescore.name}: generated xyz re-encoded under the levers)"
    )
    trans, tv = Counter(), Counter()
    losses, gains, errs, ctrl_mism = [], [], [], []
    v_losses, v_gains = [], []
    n = 0
    for m, row in S.items():
        if m in changed or not esc_l.get(m, {}).get("base"):
            continue
        row = dict(row)
        row.setdefault("eta", T[m]["eta"] == "1")
        # control: the recomputation with the SHIPPED esc and bucket must give the table's fault
        f0 = attribute(
            row, G.get(m), GM.get(m), esc_s.get(m), PB.get(m), PG.get(m), PF.get(m), None, None
        )["fault"]
        if f0 != T[m]["fault"]:
            ctrl_mism.append((m, T[m]["fault"], f0))
        r = R.get(m)
        was = row["bucket"] == "byte_exact"
        if r is None or r.get("honest_class") == "no_structure":
            trans["no_structure"] += 1
            now = was  # nothing to re-encode; the verdict cannot move
        elif r.get("indep_error") or not r.get("smiles_2_indep"):
            errs.append((m, r.get("indep_error")))
            trans["encode_error_under_lever"] += 1
            now = False
        else:
            n += 1
            now = r["smiles_2_indep"] == esc_l[m]["base"]
            trans[("pass" if was else "fail") + "->" + ("pass" if now else "fail")] += 1
        if was and not now:
            losses.append(m)
        if now and not was:
            gains.append(m)
        row_l = dict(row)
        if now != was:
            row_l["bucket"] = "byte_exact" if now else "structural"
        f1 = attribute(
            row_l, G.get(m), GM.get(m), esc_l.get(m), PB.get(m), PG.get(m), PF.get(m), None, None
        )["fault"]
        v0, v1 = f0 == "NONE", f1 == "NONE"
        tv[("pass" if v0 else "fail") + "->" + ("pass" if v1 else "fail")] += 1
        if v0 and not v1:
            v_losses.append((m, f1))
        if v1 and not v0:
            v_gains.append((m, f0))
    print(
        f"   CONTROL recomputed fault (shipped inputs) == census table: {len(S) - len(changed) - len(ctrl_mism)} / {len(S) - len(changed)}"
    )
    for m, a, b in ctrl_mism[:8]:
        print(f"      MISMATCH {m}: table {a}  recomputed {b}")
    print(f"   rows re-encoded {n};  self-consistent transitions {dict(trans)}")
    print(
        f"   SELF-CONSISTENT losses (pass -> fail, structure held fixed): {len(losses)}   gains: {len(gains)}"
    )
    for m in losses[:20]:
        print(f"      LOSS {m}  fault={T[m]['fault']}  eta={T[m]['eta']}")
    for m in gains[:20]:
        print(f"      GAIN {m}  fault={T[m]['fault']}/{T[m]['sub']}  eta={T[m]['eta']}")
    print(
        f"   VERIFIED transitions (fault recomputed with the lever's self-consistency columns) {dict(tv)}"
    )
    print(f"   VERIFIED losses {len(v_losses)}   gains {len(v_gains)}")
    for m, f in v_losses[:20]:
        print(f"      V-LOSS {m}  table={T[m]['fault']} -> {f}")
    for m, f in v_gains[:20]:
        print(f"      V-GAIN {m}  table={f} -> NONE")
    if errs:
        print(
            f"   🔴 re-encode errors under the levers: {len(errs)} (a shipped structure the lever run cannot encode)"
        )
        for m, e in errs[:10]:
            print(f"      {m}: {e}")
    res["offline"] = {
        "n": n,
        "transitions": dict(trans),
        "losses": losses,
        "gains": gains,
        "errors": errs,
        "verified_transitions": dict(tv),
        "verified_losses": v_losses,
        "verified_gains": v_gains,
        "control_mismatch": ctrl_mism,
    }


# ----------------------------------------------------------------------------- ab
def _arm_fault(d: Path, esc_path: Path, cohort: set, T):
    """Recompute the census fault per molecule from THIS arm's instruments."""
    B = _buckets(d)
    G = _jsonl(d / "g_verdict.jsonl")
    PB = _jsonl(d / "parseback.jsonl")
    PG = (
        _jsonl(d / "parseback_gen.jsonl", need_done=True)
        if (d / "parseback_gen.jsonl").exists()
        else {}
    )
    E = _jsonl(esc_path)
    GM = _jsonl(RULER_MIRROR, need_done=False)
    PF = (
        _jsonl(CENSUS / "pflags.jsonl", need_done=False)
        if (CENSUS / "pflags.jsonl").exists()
        else {}
    )
    out = {}
    for m in sorted(cohort):
        row = B.get(m)
        if row is None:
            sys.exit(f"ABORT: {d} has no bucket row for {m}")
        row = dict(row)
        row.setdefault("eta", T[m]["eta"] == "1")
        out[m] = attribute(
            row, G.get(m), GM.get(m), E.get(m), PB.get(m), PG.get(m), PF.get(m), None, None
        )
    return out, B


def cmd_ab(args, res, changed):
    T = _table()
    off, on = args.ab / "ab_off", args.ab / "ab_on"
    for d in (off, on):
        if not (d / "DONE").exists():
            sys.exit(f"ABORT: {d}/DONE missing")
        n = len(list((d / "individual_reports").glob("*.json")))
        if n != len(changed):
            sys.exit(f"ABORT: {d.name} has {n} reports for a changed set of {len(changed)}")
    print(f"\n== LIVE A/B over the changed-string set ({len(changed)} molecules)")
    for d in (off, on):
        print(f"   {d.name}: {(d / 'AB_COMMIT').read_text().strip()[:100]}")
    F_off, B_off = _arm_fault(off, SHIPPED_ESC, changed, T)
    F_on, B_on = _arm_fault(on, args.esc, changed, T)
    # CONTROL: the recomputed OFF fault must be the census table's fault
    mism = [
        (m, T[m]["fault"], F_off[m]["fault"]) for m in changed if T[m]["fault"] != F_off[m]["fault"]
    ]
    print(
        f"   CONTROL recomputed OFF fault == census table: {len(changed) - len(mism)} / {len(changed)}"
    )
    for m, a, b in mism[:10]:
        print(f"      MISMATCH {m}: table {a}  recomputed {b}")
    # the same ON strings read by the SHIPPED reader (OIN_H_FAITHFUL gates the adapter too)
    sr = on / "parseback_shippedreader.jsonl"
    if sr.exists():
        PB_on, PB_sr = _jsonl(on / "parseback.jsonl"), _jsonl(sr)
        dep = [
            m
            for m in changed
            if (PB_on.get(m) or {}).get("adapter_h_decoration")
            != (PB_sr.get(m) or {}).get("adapter_h_decoration")
        ]

        def h_ok(P):
            return sum(
                1
                for m in changed
                if (P.get(m) or {}).get("adapter_graph") in ("ISO", "ISO_MARGINAL", "ISO_CLASH")
                and (P.get(m) or {}).get("adapter_h_decoration") == "SAME"
            )

        print(
            f"   READER COUPLING: ON strings parse back ISO&SAME on {h_ok(PB_on)} rows under the ON reader, {h_ok(PB_sr)} under the SHIPPED reader; the H verdict depends on the reader's lever on {len(dep)} rows"
        )
        res.setdefault("ab", {})["reader_coupling"] = {
            "on_reader_iso_same": h_ok(PB_on),
            "shipped_reader_iso_same": h_ok(PB_sr),
            "reader_dependent": dep,
        }
    # dead-lever check, the ENCODER way: smiles_1 must differ on every row of this cohort
    same_s1 = [m for m in changed if B_off[m].get("smiles_1") == B_on[m].get("smiles_1")]
    print(
        f"   DEAD-LEVER CHECK (encoder): smiles_1 identical across arms on {len(same_s1)} rows (must be 0)"
    )
    R_off, R_on = _reports(off, changed), _reports(on, changed)
    # runtime and structure identity
    t_off = sum((R_off[m] or {}).get("metrics", {}).get("elapsed_s", 0) or 0 for m in changed)
    t_on = sum((R_on[m] or {}).get("metrics", {}).get("elapsed_s", 0) or 0 for m in changed)
    print(f"   runtime sum elapsed_s: off {t_off:.0f} s   on {t_on:.0f} s")

    def score(F, B):
        s = {m: B[m]["bucket"] == "byte_exact" for m in changed}
        v = {m: F[m]["fault"] == "NONE" for m in changed}
        return s, v

    s_off, v_off = score(F_off, B_off)
    s_on, v_on = score(F_on, B_on)

    def trans(a, b):
        c = Counter()
        for m in changed:
            c[("pass" if a[m] else "fail") + "->" + ("pass" if b[m] else "fail")] += 1
        return c

    ts, tv = trans(s_off, s_on), trans(v_off, v_on)
    print(
        f"   SELF-CONSISTENT: off {sum(s_off.values())} -> on {sum(s_on.values())}   {dict(ts)}   net {sum(s_on.values()) - sum(s_off.values()):+d}"
    )
    print(
        f"   VERIFIED:        off {sum(v_off.values())} -> on {sum(v_on.values())}   {dict(tv)}   net {sum(v_on.values()) - sum(v_off.values()):+d}"
    )
    print("   ON-arm fault of the rows still failing VERIFIED, by fault:")
    for f, n in Counter(
        F_on[m]["fault"] + ("/" + F_on[m]["sub"] if F_on[m]["fault"].startswith("G_C") else "")
        for m in changed
        if not v_on[m]
    ).most_common(12):
        print(f"      {f:36s} {n:4d}")
    print("   transitions by census fault of the row (VERIFIED):")
    by = defaultdict(Counter)
    for m in changed:
        by[T[m]["fault"]][
            ("pass" if v_off[m] else "fail") + "->" + ("pass" if v_on[m] else "fail")
        ] += 1
    for f, c in sorted(by.items(), key=lambda x: -sum(x[1].values())):
        print(f"      {f:18s} {dict(c)}")
    losses_v = [m for m in changed if v_off[m] and not v_on[m]]
    print(f"   VERIFIED losses ({len(losses_v)}):")
    for m in losses_v[:25]:
        print(
            f"      {m}  {T[m]['metal']} eta={T[m]['eta']}  on-fault {F_on[m]['fault']}/{F_on[m]['sub']}  on-bucket {B_on[m]['bucket']}"
        )
    res["ab"] = {
        **res.get("ab", {}),
        "n": len(changed),
        "control_mismatch": mism,
        "dead_lever_same_s1": same_s1,
        "self": {"off": sum(s_off.values()), "on": sum(s_on.values()), "transitions": dict(ts)},
        "verified": {"off": sum(v_off.values()), "on": sum(v_on.values()), "transitions": dict(tv)},
        "verified_losses": losses_v,
        "verified_gains": [m for m in changed if v_on[m] and not v_off[m]],
        "self_losses": [m for m in changed if s_off[m] and not s_on[m]],
        "self_gains": [m for m in changed if s_on[m] and not s_off[m]],
        "on_fault": {m: F_on[m]["fault"] for m in changed},
        "runtime": {"off": t_off, "on": t_on},
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("cmd", choices=["changed", "audit", "offline", "ab", "all"])
    ap.add_argument("--esc", type=Path, default=LANE / "e_selfconsistency_fix2.jsonl")
    ap.add_argument(
        "--single",
        type=Path,
        default=LANE / "changed_single.jsonl",
        help="lever_parseback.py output over the changed set, arms shipped,hfaith,rc1prop",
    )
    ap.add_argument("--rescore", type=Path, default=LANE / "rescore_fix2")
    ap.add_argument("--ab", type=Path, default=LANE / "ab")
    ap.add_argument("--write-set", type=Path, default=None)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    res = {}
    changed = cmd_changed(args, res) if args.cmd in ("changed", "offline", "ab", "all") else None
    if args.cmd in ("audit", "all"):
        cmd_audit(args, res)
    if args.cmd in ("offline", "all"):
        cmd_offline(args, res, changed)
    if args.cmd in ("ab", "all"):
        cmd_ab(args, res, changed)
    if args.cmd == "all" and "ab" in res and "offline" in res:
        d_self = (
            res["ab"]["self"]["on"]
            - res["ab"]["self"]["off"]
            + len(res["offline"]["gains"])
            - len(res["offline"]["losses"])
        )
        d_ver = (
            res["ab"]["verified"]["on"]
            - res["ab"]["verified"]["off"]
            + len(res["offline"]["verified_gains"])
            - len(res["offline"]["verified_losses"])
        )
        print(
            "\n== PROJECTION onto the 5,000 (live arm + offline re-score; unchanged rows exact by determinism)"
        )
        print(
            f"   self-consistent {PASS_SELF} -> {PASS_SELF + d_self}  ({100 * PASS_SELF / N_COHORT:.2f}% -> {100 * (PASS_SELF + d_self) / N_COHORT:.2f}%)"
        )
        print(
            f"   VERIFIED        {PASS_VERIFIED} -> {PASS_VERIFIED + d_ver}  ({100 * PASS_VERIFIED / N_COHORT:.2f}% -> {100 * (PASS_VERIFIED + d_ver) / N_COHORT:.2f}%)"
        )
        res["projection"] = {"self": PASS_SELF + d_self, "verified": PASS_VERIFIED + d_ver}
    if args.out:
        args.out.write_text(json.dumps(res, indent=1, default=list))
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
