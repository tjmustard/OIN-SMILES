"""Stage the v0.4.19 SERIALIZER lane for ``/freeze-measurements``. Reads only.

Reuses ``tools/v0418/freeze_stage.py``'s scrub -> gzip(mtime=0) -> verify-what-was-written machinery
with this lane's directory and prefix (``v0419_ser_``), its own release ``v0.4.19-serializer``.

What is frozen, and what each file answers:

    pop1_report.txt / pop1b_report.txt / pop1*.jsonl.gz     the parse-back probe, 296 x arms (S1)
    changed_set.txt / changed_single.jsonl.gz               the 80 changed strings and per-lever attribution
    esc_fix2_summary.txt                                    the full-cohort self-consistency run's table
    e_selfconsistency_fix2.jsonl.gz                         the run itself (5,000 x 11)
    rescore_fix2.jsonl.gz                                   the offline re-score under the levers
    serializer_ab_report.txt / .json                        everything, in the order to believe it
    ab_on_{bucket_report_honest.json,g_verdict.jsonl,parseback.jsonl,parseback_shippedreader.jsonl}.gz
    ab_commits.tsv

    <main>/.venv/bin/python tools/v0419/freeze_stage.py
    <main>/.venv/bin/python tools/v0419/freeze_stage.py --verify <main>/measurements/v0.4.19-serializer
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1] / "v0418"))

import freeze_stage as fs  # noqa: E402

fs.R = fs.DATA / "results-v0.4.19-serializer"
fs.STAGE = fs.R / "freeze"
fs.P = "v0419_ser_"
EXPECTED_ARMS = {  # rule v3 of OIN_CAP_IGNORES_METAL; cap alone, and all three levers
    # offline_v_losses: OHUTIV (NONE -> G_HCOUNT, a generated structure's H count); fragile: 2 failing rows
    "cap": {
        "changed": 31,
        "self": (5, 2),
        "verified": (6, 0),
        "offline": (13, 5),
        "offline_v_losses": 1,
        "fragile": 2,
    },
    "fix3": {
        "changed": 109,
        "self": (33, 3),
        "verified": (30, 0),
        "offline": (6, 5),
        "offline_v_losses": 1,
        "fragile": 2,
    },
}
EXPECTED = {
    "self": (29, 1),
    "verified": (24, 0),
    "changed": 80,
    "offline_gains": 5,
    "offline_losses": 0,
}

TEXT = [
    "pop1_report.txt",
    "pop1b_report.txt",
    "changed_set.txt",
    "serializer_ab_report.txt",
    "serializer_ab_report.json",
]
GZ = {
    "pop1.jsonl.gz": "pop1.jsonl",
    "pop1b.jsonl.gz": "pop1b.jsonl",
    "changed_single.jsonl.gz": "changed_single.jsonl",
    "e_selfconsistency_fix2.jsonl.gz": "e_selfconsistency_fix2.jsonl",
    "rescore_fix2.jsonl.gz": "rescore_fix2/honest_rescore.jsonl",
    "ab_on_bucket_report_honest.json.gz": "ab/ab_on/bucket_report_honest.json",
    "ab_on_g_verdict.jsonl.gz": "ab/ab_on/g_verdict.jsonl",
    "ab_on_parseback.jsonl.gz": "ab/ab_on/parseback.jsonl",
    "ab_on_parseback_shippedreader.jsonl.gz": "ab/ab_on/parseback_shippedreader.jsonl",
}


ARMS = ("cap", "fix3")


def _esc_summary(name: str = "esc_fix2.log") -> str:
    log = (fs.R / name).read_text()
    keep = [
        ln
        for ln in log.splitlines()
        if not ln.startswith(("!!!", "[", "  ")) and "Warning" not in ln
    ]
    return "\n".join(keep) + "\n"


def stage():
    if fs.STAGE.exists():
        shutil.rmtree(fs.STAGE)  # a stale staged file would be harvested as if it were current
    fs.STAGE.mkdir(parents=True)
    total = n = 0
    for f in TEXT:
        if not (fs.R / f).exists():
            sys.exit(f"ABORT: {fs.R / f} is missing -- refusing to stage a partial release")
        total += fs._write(fs.P + f, (fs.R / f).read_text(), False)
        n += 1
    total += fs._write(fs.P + "esc_fix2_summary.txt", _esc_summary(), False)
    n += 1
    for name, src in GZ.items():
        p = fs.R / src
        if not p.exists():
            sys.exit(f"ABORT: {p} is missing -- refusing to stage a partial release")
        txt = p.read_text()
        if src.endswith(".jsonl") and not txt.rstrip().splitlines()[-1].startswith("#DONE"):
            sys.exit(f"ABORT: {p} has no #DONE trailer")
        total += fs._write(fs.P + name, txt, True)
        n += 1
    commits = (
        "off\t"
        + (fs.R / "ab/ab_off/AB_COMMIT").read_text().strip()
        + "\non\t"
        + (fs.R / "ab/ab_on/AB_COMMIT").read_text().strip()
        + "\n"
    )
    total += fs._write(fs.P + "ab_commits.tsv", "arm\tlaunched_from\n" + commits, False)
    n += 1
    # the later arms (cap = OIN_CAP_IGNORES_METAL alone; fix3 = all three), each with its own
    # audit, re-score, changed set and live arm -- staged only once its report exists, and then
    # completely or not at all
    if (fs.R / "pop1c_report.txt").exists():
        total += fs._write(
            fs.P + "pop1c_report.txt", (fs.R / "pop1c_report.txt").read_text(), False
        )
        total += fs._write(fs.P + "pop1c.jsonl.gz", (fs.R / "pop1c.jsonl").read_text(), True)
        n += 2
    for arm in ARMS:
        rep = fs.R / f"serializer_ab_report_{arm}.txt"
        if not rep.exists():
            print(f"  (arm {arm}: no report yet -- not staged)")
            continue
        text = [
            f"serializer_ab_report_{arm}.txt",
            f"serializer_ab_report_{arm}.json",
            f"changed_set_{arm}.txt",
        ]
        gz = {
            f"changed_single_{arm}.jsonl.gz": f"changed_single_{arm}.jsonl",
            f"e_selfconsistency_{arm}.jsonl.gz": f"e_selfconsistency_{arm}.jsonl",
            f"rescore_{arm}.jsonl.gz": f"rescore_{arm}/honest_rescore.jsonl",
            f"ab_{arm}_on_bucket_report_honest.json.gz": f"ab_{arm}/ab_on/bucket_report_honest.json",
            f"ab_{arm}_on_g_verdict.jsonl.gz": f"ab_{arm}/ab_on/g_verdict.jsonl",
            f"ab_{arm}_on_parseback.jsonl.gz": f"ab_{arm}/ab_on/parseback.jsonl",
            f"ab_{arm}_on_parseback_shippedreader.jsonl.gz": f"ab_{arm}/ab_on/parseback_shippedreader.jsonl",
        }
        for f in text:
            if not (fs.R / f).exists():
                sys.exit(f"ABORT: {fs.R / f} is missing -- refusing to stage a partial arm {arm}")
            total += fs._write(fs.P + f, (fs.R / f).read_text(), False)
            n += 1
        total += fs._write(fs.P + f"esc_{arm}_summary.txt", _esc_summary(f"esc_{arm}.log"), False)
        n += 1
        for name, src in gz.items():
            q = fs.R / src
            if not q.exists():
                sys.exit(f"ABORT: {q} is missing -- refusing to stage a partial arm {arm}")
            txt = q.read_text()
            if src.endswith(".jsonl") and not txt.rstrip().splitlines()[-1].startswith("#DONE"):
                sys.exit(f"ABORT: {q} has no #DONE trailer")
            total += fs._write(fs.P + name, txt, True)
            n += 1
        commits = (
            "off\t"
            + (fs.R / f"ab_{arm}/ab_off/AB_COMMIT").read_text().strip()
            + "\non\t"
            + (fs.R / f"ab_{arm}/ab_on/AB_COMMIT").read_text().strip()
            + "\n"
        )
        total += fs._write(fs.P + f"ab_commits_{arm}.tsv", "arm\tlaunched_from\n" + commits, False)
        n += 1
    print(f"staged {n} files, {total / 1024:.0f} KB -> {fs.STAGE}")


def verify(frozen: Path):
    """Re-derive the headline from the frozen tree alone -- no results directory, no worktree."""
    P = fs.P
    rep = json.loads((frozen / f"{P}serializer_ab_report.json").read_text())
    changed = [
        ln
        for ln in (frozen / f"{P}changed_set.txt").read_text().splitlines()
        if ln and not ln.startswith("#")
    ]
    ok = True

    def check(label, got, want):
        nonlocal ok
        flag = "OK " if got == want else "🔴 "
        ok &= got == want
        print(f"  {flag}{label}: {got}  (expected {want})")

    check("changed-string set", len(changed), EXPECTED["changed"])
    # recompute the A/B from the frozen arm files, not from the report's own numbers
    B = {
        r["molecule"].removesuffix(".xyz"): r
        for r in json.loads(
            gzip.open(frozen / f"{P}ab_on_bucket_report_honest.json.gz", "rt").read()
        )
    }
    on_pass = sum(1 for m in changed if B[m]["bucket"] == "byte_exact")
    check(
        "ON arm byte_exact passes (from the frozen bucket report)", on_pass, rep["ab"]["self"]["on"]
    )
    check(
        "self-consistent gains/losses",
        (len(rep["ab"]["self_gains"]), len(rep["ab"]["self_losses"])),
        EXPECTED["self"],
    )
    check(
        "VERIFIED gains/losses",
        (len(rep["ab"]["verified_gains"]), len(rep["ab"]["verified_losses"])),
        EXPECTED["verified"],
    )
    check(
        "offline gains/losses",
        (len(rep["offline"]["gains"]), len(rep["offline"]["losses"])),
        (EXPECTED["offline_gains"], EXPECTED["offline_losses"]),
    )
    check("offline VERIFIED losses", len(rep["offline"]["verified_losses"]), 0)
    check(
        "canonicality: molecules that became fragile",
        len(rep["audit"]["fragile"]["became_fragile"]),
        0,
    )
    check(
        "mirror 2x2 unchanged",
        rep["audit"]["mirror_2x2"]["shipped"] == rep["audit"]["mirror_2x2"]["lever"],
        True,
    )
    check(
        "controls: OFF fault == census table (A/B, offline)",
        (len(rep["ab"]["control_mismatch"]), len(rep["offline"]["control_mismatch"])),
        (0, 0),
    )
    check("dead-lever check: identical smiles_1 rows", len(rep["ab"]["dead_lever_same_s1"]), 0)
    rc = rep["ab"].get("reader_coupling", {})
    print(
        f"  reader coupling: ON strings ISO&SAME under ON reader {rc.get('on_reader_iso_same')} / under SHIPPED reader {rc.get('shipped_reader_iso_same')}; reader-dependent rows {len(rc.get('reader_dependent', []))}"
    )
    proj = rep.get("projection", {})
    print(
        f"  projection: self-consistent {proj.get('self')} / 5000 = {100 * proj.get('self', 0) / 5000:.2f}%   VERIFIED {proj.get('verified')} / 5000 = {100 * proj.get('verified', 0) / 5000:.2f}%"
    )
    for arm in ARMS:
        rp = frozen / f"{P}serializer_ab_report_{arm}.json"
        if not rp.exists():
            continue
        r = json.loads(rp.read_text())
        exp = EXPECTED_ARMS.get(arm)
        print(f"  -- arm {arm}")
        ch = [
            ln
            for ln in (frozen / f"{P}changed_set_{arm}.txt").read_text().splitlines()
            if ln and not ln.startswith("#")
        ]
        Bm = {
            x["molecule"].removesuffix(".xyz"): x
            for x in json.loads(
                gzip.open(frozen / f"{P}ab_{arm}_on_bucket_report_honest.json.gz", "rt").read()
            )
        }
        check(
            f"  {arm} ON arm byte_exact passes (frozen bucket report)",
            sum(1 for m in ch if Bm[m]["bucket"] == "byte_exact"),
            r["ab"]["self"]["on"],
        )
        if exp:
            check(f"  {arm} changed-string set", len(ch), exp["changed"])
            check(
                f"  {arm} self-consistent gains/losses",
                (len(r["ab"]["self_gains"]), len(r["ab"]["self_losses"])),
                exp["self"],
            )
            check(
                f"  {arm} VERIFIED gains/losses",
                (len(r["ab"]["verified_gains"]), len(r["ab"]["verified_losses"])),
                exp["verified"],
            )
            check(
                f"  {arm} offline gains/losses",
                (len(r["offline"]["gains"]), len(r["offline"]["losses"])),
                exp["offline"],
            )
        check(
            f"  {arm} offline VERIFIED losses",
            len(r["offline"]["verified_losses"]),
            (exp or {}).get("offline_v_losses", 0),
        )
        check(
            f"  {arm} became fragile",
            len(r["audit"]["fragile"]["became_fragile"]),
            (exp or {}).get("fragile", 0),
        )
        check(
            f"  {arm} controls (A/B, offline)",
            (len(r["ab"]["control_mismatch"]), len(r["offline"]["control_mismatch"])),
            (0, 0),
        )
        pj = r.get("projection", {})
        print(
            f"     projection: self-consistent {pj.get('self')} ({100 * pj.get('self', 0) / 5000:.2f}%)   VERIFIED {pj.get('verified')} ({100 * pj.get('verified', 0) / 5000:.2f}%)"
        )
    print("VERIFY", "OK" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--verify", type=Path, help="a frozen release dir; re-derive the headline")
    args = ap.parse_args()
    if args.verify:
        sys.exit(0 if verify(args.verify) else 1)
    stage()


if __name__ == "__main__":
    main()
