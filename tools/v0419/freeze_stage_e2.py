"""Stage the v0.4.19 E2 lane (perception order: OIN_N_VALENCE_2, OIN_CANONICAL_RESONANCE) for
``/freeze-measurements``. Reads only.

Reuses ``tools/v0418/freeze_stage.py``'s scrub -> gzip(mtime=0) -> verify-what-was-written
machinery with this lane's directory (``results-v0.4.19-e2``), prefix (``v0419_e2_``) and release
(``v0.4.19-e2``). The lane stands ON v0.4.19: its baseline is the v0.4.19 release sweep and the
"shipped" self-consistency run is the serializer lane's fix3 run (the three promoted levers ON).

What is frozen, in the order the lane happened:

    fragile105.txt, esc105_e2b_final.jsonl.gz, e2b_final_fragility.jsonl.gz
        the 105 census E2_P_FRAGILE rows under the first cut (88 stable, renumbering 94/100)
    n2_valence_scan.json.gz            how many ligands carry a two-coordinate N, product sizes
    e_selfconsistency_e2b.jsonl.gz     CUT 1 (valence 2 appended for every two-coordinate N) on the
    canonicality_audit_e2b.{txt,json}  whole cohort: 343 strings moved, 167 of them VERIFIED passes;
    e2b_cohort_fragility.txt           fragile fixed 215 / broken 53 (noise +5/-44) -- REFUTED
    changed_set_e2b.txt, changed_single_e2b.jsonl.gz, became_fragile_e2b.txt, movers_e2b.txt
    e_selfconsistency_e2d_movers.jsonl.gz, e2d_movers_vs_{shipped,e2b}.txt
                                       CUT 2 (a second pass inside AC2BO) on the 443 movers:
                                       251 still moved -- REFUTED
    movers_single_e2e.jsonl.gz         CUT 3 (a ladder-level fallback on a bond-order guess):
                                       177 of 443 moved, 59 verified passes -- REFUTED
    movers_single_e2f.jsonl.gz, n2_scoped_movers.json
                                       CUT 4 (the option scoped to pyrrolide N): 84 moved, 66 already
                                       fragile
    e_selfconsistency_{e2f,res}.jsonl.gz, {e2f,res}_cohort_fragility.txt, changed_{e2f,res}.txt,
    changed_set_{e2f,res}.txt          the whole cohort under cut 4 + RES, and RES alone
    rescore_{e2f,res}.jsonl.gz         the offline re-score of the unchanged-string rows
    serializer_ab_report_{e2f,res}.{txt,json}, ab_{e2f,res}_on_*.gz, commits.tsv
                                       the live A/B over each changed set, and the projection

    <main>/.venv/bin/python tools/v0419/freeze_stage_e2.py
    <main>/.venv/bin/python tools/harvest_measurements.py --release v0.4.19-e2 \\
        --from <results-v0.4.19-e2>/freeze --dry-run
    <main>/.venv/bin/python tools/v0419/freeze_stage_e2.py --verify <main>/measurements/v0.4.19-e2

NOT staged: the derived per-molecule fragility jsonl of the whole cohort (2.9 MB; --verify
recomputes the 2x2 from the frozen self-consistency runs), the 105-row intermediate passes.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
# v0418 FIRST: this directory has its own freeze_stage.py (the serializer lane's), which must not
# shadow the shared machinery. e2_fragility_report is found through the script's own directory.
sys.path.insert(0, str(HERE.parents[1] / "v0418"))

import freeze_stage as fs  # noqa: E402
from e2_fragility_report import fragile  # noqa: E402

fs.R = fs.DATA / "results-v0.4.19-e2"
fs.STAGE = fs.R / "freeze"
fs.P = "v0419_e2_"

TEXT = [
    "fragile105.txt",
    "movers_e2b.txt",
    "became_fragile_e2b.txt",
    "changed_set_e2b.txt",
    "canonicality_audit_e2b.txt",
    "canonicality_audit_e2b.json",
    "e2b_cohort_fragility.txt",
    "e2d_movers_vs_shipped.txt",
    "e2d_movers_vs_e2b.txt",
    "n2_scoped_movers.json",
    "e2f_cohort_fragility.txt",
    "res_cohort_fragility.txt",
    "changed_e2f.txt",
    "changed_res.txt",
    "changed_set_e2f.txt",
    "changed_set_res.txt",
    "serializer_ab_report_e2f.txt",
    "serializer_ab_report_e2f.json",
    "serializer_ab_report_res.txt",
    "serializer_ab_report_res.json",
]
# run logs reduced to their summary lines (the rest is RDKit/aligner warnings per molecule)
LOGS = {
    "esc_e2b_summary.txt": "esc_e2b.log",
    "esc_e2d_movers_summary.txt": "esc_e2d_movers.log",
    "esc_e2f_summary.txt": "esc_e2f.log",
    "esc_res_summary.txt": "esc_res.log",
    "changed_single_e2b_summary.txt": "changed_single_e2b.log",
    "movers_single_e2e_summary.txt": "movers_single_e2e.log",
    "movers_single_e2f_summary.txt": "movers_single_e2f.log",
}
GZ = {
    "esc105_e2b_final.jsonl.gz": "esc105_e2b_final.jsonl",
    "e2b_final_fragility.jsonl.gz": "e2b_final_fragility.jsonl",
    "n2_valence_scan.json.gz": "n2_valence_scan.json",
    "e_selfconsistency_e2b.jsonl.gz": "e_selfconsistency_e2b.jsonl",
    "e_selfconsistency_e2d_movers.jsonl.gz": "e_selfconsistency_e2d_movers.jsonl",
    "e_selfconsistency_e2f.jsonl.gz": "e_selfconsistency_e2f.jsonl",
    "e_selfconsistency_res.jsonl.gz": "e_selfconsistency_res.jsonl",
    "changed_single_e2b.jsonl.gz": "changed_single_e2b.jsonl",
    "movers_single_e2e.jsonl.gz": "movers_single_e2e.jsonl",
    "movers_single_e2f.jsonl.gz": "movers_single_e2f.jsonl",
    "rescore_e2f.jsonl.gz": "rescore_e2f/honest_rescore.jsonl",
    "rescore_res.jsonl.gz": "rescore_res/honest_rescore.jsonl",
}
AB_ARMS = ("e2f", "res")
AB_FILES = (
    "bucket_report_honest.json",
    "g_verdict.jsonl",
    "parseback.jsonl",
    "parseback_shippedreader.jsonl",
)
# files whose producer writes a #DONE trailer -- an unfinished run must not be frozen
NEEDS_DONE = (
    "e_selfconsistency_",
    "changed_single_",
    "movers_single_",
    "honest_rescore",
    "esc105_",
)

EXPECTED = {
    # arm -> (changed strings, fragile fixed, fragile broken) vs the shipped (fix3) run
    "e2b": (343, 215, 53),
    "e2f": (136, 120, 5),
    "res": (56, 81, 1),
}
EXPECTED_AB = {
    # arm -> self (gains, losses), VERIFIED (gains, losses), offline self, offline VERIFIED,
    #        projection (self, VERIFIED), dead-lever rows (RAXJEH: its encode sits on the resonance
    #        isolation's 120 CPU-s budget under RES; the harness and the audit fell on opposite sides)
    "e2f": ((69, 7), (66, 7), (17, 2), (11, 1), (4466, 4025), 1),
    "res": ((23, 3), (20, 3), (18, 1), (12, 0), (4426, 3985), 1),
}


def _summary(name: str) -> str:
    keep = [
        ln
        for ln in (fs.R / name).read_text().splitlines()
        if ln.strip()
        and not ln.startswith(("!!!", "[", "{", "/"))
        and "Warning" not in ln
        and "warn(" not in ln
        and not ln.lstrip().startswith(("dummy_mol", "R, rmsd", "rotation"))
    ]
    return "\n".join(keep) + "\n"


def _put(name: str, src: Path, gz: bool) -> int:
    if not src.exists():
        sys.exit(f"ABORT: {src} is missing -- refusing to stage a partial lane")
    txt = src.read_text()
    if src.suffix == ".jsonl" and src.name.startswith(NEEDS_DONE):
        if not txt.rstrip().splitlines()[-1].startswith("#DONE"):
            sys.exit(f"ABORT: {src} has no #DONE trailer -- an unfinished run must not be frozen")
    return fs._write(fs.P + name, txt, gz)


def stage():
    if fs.STAGE.exists():
        shutil.rmtree(fs.STAGE)  # a stale staged file would be harvested as if it were current
    fs.STAGE.mkdir(parents=True)
    total = n = 0
    for f in TEXT:
        total += _put(f, fs.R / f, False)
        n += 1
    for name, src in LOGS.items():
        if not (fs.R / src).exists():
            sys.exit(f"ABORT: {fs.R / src} is missing")
        total += fs._write(fs.P + name, _summary(src), False)
        n += 1
    for name, src in GZ.items():
        total += _put(name, fs.R / src, True)
        n += 1
    rows = ["run\tcommit"]
    for p in sorted(fs.R.glob("E2_COMMIT_*")):
        rows.append(f"{p.name.removeprefix('E2_COMMIT_')}\t{p.read_text().strip()}")
    for arm in AB_ARMS:
        for f in AB_FILES:
            total += _put(f"ab_{arm}_on_{f}.gz", fs.R / f"ab_{arm}" / "ab_on" / f, True)
            n += 1
        for side in ("off", "on"):
            rows.append(
                f"ab_{arm}_{side}\t{(fs.R / f'ab_{arm}' / f'ab_{side}' / 'AB_COMMIT').read_text().strip()}"
            )
    total += fs._write(fs.P + "commits.tsv", "\n".join(rows) + "\n", False)
    n += 1
    print(f"staged {n} files, {total / 1024:.0f} KB -> {fs.STAGE}")


def _jsonl_gz(p: Path) -> dict:
    out = {}
    for ln in gzip.open(p, "rt"):
        if ln.startswith("{"):
            r = json.loads(ln)
            out[r["molecule"]] = r
    return out


def verify(frozen: Path) -> bool:
    """Re-derive the lane's numbers from the frozen trees alone -- this release, the v0.4.19
    release sweep's frozen bucket report (smiles_1) and the serializer lane's frozen fix3 run (the
    shipped self-consistency columns). No results directory, no worktree."""
    P = fs.P
    root = frozen.parent
    ship = _jsonl_gz(root / "v0.4.19-serializer" / "v0419_ser_e_selfconsistency_fix3.jsonl.gz")
    sweep = {
        r["molecule"].removesuffix(".xyz"): r.get("smiles_1")
        for r in json.loads(
            gzip.open(
                root / "v0.4.19-release-sweep" / "v0419_rsweep_bucket_report_honest.json.gz", "rt"
            ).read()
        )
    }
    ok = True

    def check(label, got, want):
        nonlocal ok
        good = got == want
        ok &= good
        print(f"  {'OK ' if good else '🔴 '}{label}: {got}  (expected {want})")

    for arm, (n_changed, n_fixed, n_broken) in EXPECTED.items():
        esc = _jsonl_gz(frozen / f"{P}e_selfconsistency_{arm}.jsonl.gz")
        check(f"{arm}: molecules in the self-consistency run", len(esc), 5000)
        changed = sorted(m for m, r in esc.items() if r.get("base") and r["base"] != sweep.get(m))
        check(f"{arm}: base string != v0.4.19 sweep smiles_1", len(changed), n_changed)
        listed = [
            ln
            for ln in (frozen / f"{P}changed_set_{arm}.txt").read_text().splitlines()
            if ln and not ln.startswith("#")
        ]
        check(f"{arm}: frozen changed set == recomputed", sorted(listed) == changed, True)
        fixed = broken = 0
        for m, r in esc.items():
            a, s = fragile(r), fragile(ship.get(m))
            fixed += bool(s is True and a is False)
            broken += bool(s is False and a is True)
        check(f"{arm}: fragile fixed / broken vs shipped", (fixed, broken), (n_fixed, n_broken))

    for arm, (self_gl, ver_gl, off_self, off_ver, proj, dead) in EXPECTED_AB.items():
        rep = json.loads((frozen / f"{P}serializer_ab_report_{arm}.json").read_text())
        ch = [
            ln
            for ln in (frozen / f"{P}changed_set_{arm}.txt").read_text().splitlines()
            if ln and not ln.startswith("#")
        ]
        B = {
            x["molecule"].removesuffix(".xyz"): x
            for x in json.loads(
                gzip.open(frozen / f"{P}ab_{arm}_on_bucket_report_honest.json.gz", "rt").read()
            )
        }
        check(
            f"{arm}: ON arm byte_exact passes (frozen bucket report)",
            sum(1 for m in ch if B[m]["bucket"] == "byte_exact"),
            rep["ab"]["self"]["on"],
        )
        ab, off = rep["ab"], rep["offline"]
        check(
            f"{arm}: live self gains/losses",
            (len(ab["self_gains"]), len(ab["self_losses"])),
            self_gl,
        )
        check(
            f"{arm}: live VERIFIED gains/losses",
            (len(ab["verified_gains"]), len(ab["verified_losses"])),
            ver_gl,
        )
        check(
            f"{arm}: offline self gains/losses", (len(off["gains"]), len(off["losses"])), off_self
        )
        check(
            f"{arm}: offline VERIFIED gains/losses",
            (len(off["verified_gains"]), len(off["verified_losses"])),
            off_ver,
        )
        check(
            f"{arm}: controls (A/B, offline)",
            (len(ab["control_mismatch"]), len(off["control_mismatch"])),
            (0, 0),
        )
        check(f"{arm}: dead-lever rows (RAXJEH, explained)", len(ab["dead_lever_same_s1"]), dead)
        pj = rep["projection"]
        check(f"{arm}: projection (self, VERIFIED)", (pj["self"], pj["verified"]), proj)
        print(
            f"     = {100 * pj['self'] / 5000:.2f}% self-consistent / {100 * pj['verified'] / 5000:.2f}% VERIFIED"
            "  (v0.4.19: 87.78% / 79.12%)"
        )
    print("VERIFY", "OK" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--verify", type=Path, help="a frozen release dir; re-derive the lane")
    args = ap.parse_args()
    if args.verify:
        sys.exit(0 if verify(args.verify) else 1)
    stage()


if __name__ == "__main__":
    main()
