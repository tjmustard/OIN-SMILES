"""ARM 2 golden re-freeze, steps 3-4: read the gate's rows, say WHY each moved, splice. No chemistry.

    diff     merge the `on` shards, apply the gate's own comparison to every golden row, cross-check
             field 2 against tools/v0417/arm2_field2_audit.py, and build the control cohorts
             (symlink dirs) holding exactly the rows `on` did not reproduce
    splice   read the `off` control for those rows, give each ONE reason, and write the goldens

THE REASON IS THE POINT. A row that differs from the golden is not thereby a row the v0.4.17 lever
moved. Generation is seeded, so with ``OIN_EXACT_DONOR_FOLD=0`` a healthy golden row comes back
byte-for-byte. Each unreproduced row is therefore exactly one of:

    LEVER          `off` reproduces the golden row; `on` does not       -> the exact fold moved it
    STALE          `off` does not reproduce it and `off` == `on`        -> wrong BEFORE this lane
    STALE+LEVER    `off` does not reproduce it and `off` != `on`        -> both

and a row nothing gates (``NO_ENCODE@``) or gates on field 2 only (``NO_STRUCTURE@``) is never
re-frozen for a field the gate does not read.

WHAT IS SPLICED. Fields 1-6 (name, sha_in, sha_out, len_in, len_out, eta) come from the fresh `on`
row. Fields 7+ are PRESERVED: in a v0.4.9 golden field 7 is the runtime band ``--band`` filters on,
while a fresh row carries ``xyz_sha`` there (v0.4.14 found this by breaking it). A fresh row with no
``sha_out`` is NOT spliced by this tool -- whether that is a budget or a defect is a judgement, and
it aborts so that a person makes it.

``# MANIFEST_SHA256`` is recomputed: ARM 2 never verifies it, so a stale one goes unseen.

    <main>/.venv/bin/python tools/v0417/arm2_refreeze.py diff
    <main>/.venv/bin/python tools/v0417/arm2_refreeze.py splice [--write]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
OUT = DATA / "results-v0.4.17-exactfold" / "arm2_refreeze"
SHARDS = 6
GOLDENS = {
    "v047": (REPO / "tools/gate_v047_arm2_golden.tsv", DATA / "cohort-v047-slow100"),
    "v049": (REPO / "tools/gate_v049_arm2_golden.tsv", DATA / "cohort-v049-strata"),
}


def control_cohort(tag: str) -> Path:
    return DATA / f"cohort-v0.4.17-arm2-control-{tag}"


def golden_rows(path: Path):
    return [ln.split("\t") for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]


def fresh_rows(arm: str, tag: str):
    """Every shard must carry its ``#DONE`` -- a killed shard and a clean one look the same."""
    rows, done = {}, 0
    for i in range(1, SHARDS + 1):
        p = OUT / arm / f"{tag}_{i}.tsv"
        if not p.exists():
            sys.exit(f"ABORT: {p} is missing -- arm '{arm}' has not finished")
        lines = p.read_text().splitlines()
        if not lines or not lines[-1].startswith("#DONE"):
            sys.exit(f"ABORT: {p} has no #DONE trailer -- refusing a partial shard")
        done += int(lines[-1].split()[1])
        for ln in lines[:-1]:
            f = ln.split("\t")
            if f[0] in rows:
                sys.exit(f"ABORT: {f[0]} appears in two shards of {arm}/{tag}")
            rows[f[0]] = f
    if done != len(rows):
        sys.exit(f"ABORT: {arm}/{tag} shards say #DONE {done} but hold {len(rows)} rows")
    return rows


def gate_verdict(g, f):
    """``gate_v047.sh::run_arm2``'s comparison, in the same order. Returns (verdict, field)."""
    g1, g2 = g[1], g[2]
    f1, f2 = (f[1], f[2]) if f else ("", "")
    if f is None:
        return "MISSING", None
    if g1.startswith("NO_ENCODE@"):
        return "OBSERVE", None
    if f1 != g1:
        return "MISMATCH", 2
    if g2.startswith("NO_STRUCTURE@"):
        return "SHA1_ONLY", None
    if g2 == "NO_STRUCTURE_DET":
        return ("MISMATCH", 3) if f2 else ("DET", None)
    return ("MISMATCH", 3) if f2 != g2 else ("COMPARED", None)


def _audit(tag):
    p = OUT / f"field2_audit_{tag}.jsonl"
    lines = p.read_text().splitlines()
    if not lines[-1].startswith("#DONE"):
        sys.exit(f"ABORT: {p} has no #DONE trailer")
    return {r["molecule"]: r for r in map(json.loads, lines[:-1])}


def cmd_diff(_args):
    for tag, (gpath, cohort) in GOLDENS.items():
        G, F, A = golden_rows(gpath), fresh_rows("on", tag), _audit(tag)
        tally, owed = Counter(), []
        for g in G:
            v, field = gate_verdict(g, F.get(g[0]))
            tally[v if field is None else f"MISMATCH_f{field}"] += 1
            if v in ("MISMATCH", "MISSING"):
                owed.append((g[0], field, F.get(g[0], [""] * 8)[-1]))
        extra = sorted(set(F) - {g[0] for g in G})
        print(
            f"\n=== {tag}: {len(G)} golden rows, {len(F)} fresh rows ({len(extra)} ungated: {extra})"
        )
        print("   ", dict(sorted(tally.items())))

        # Two instruments, one question. The gate's field-2 mismatches and the encode-only audit's
        # "shipped != golden" were produced by different code paths and must name the same rows.
        by_gate = {m for m, field, _s in owed if field == 2}
        by_audit = {
            m
            for m, r in A.items()
            if r["class"] != "SENTINEL" and r["shipped"].get("sha") != r["golden_sha1"]
        }
        if by_gate != by_audit:
            sys.exit(
                f"ABORT: field 2 disagrees. gate only: {sorted(by_gate - by_audit)}  "
                f"audit only: {sorted(by_audit - by_gate)}"
            )
        print(f"    field 2: gate and encode-only audit agree on all {len(by_gate)} rows")

        cdir = control_cohort(tag)
        cdir.mkdir(exist_ok=True)
        for old in cdir.iterdir():
            old.unlink()  # a stale member would be run, and spliced, as if it were owed
        for m, field, status in owed:
            os.symlink(os.path.realpath(cohort / f"{m}.xyz"), cdir / f"{m}.xyz")
            print(f"    owed  {m:22s} field {field}   fresh status: {status[:70]}")
        (OUT / f"owed_{tag}.txt").write_text("".join(f"{m}\n" for m, _f, _s in owed))
        print(f"    control cohort: {len(owed)} rows -> {cdir.name}")
    sys.stdout.flush()


def cmd_splice(args):
    report = []
    for tag, (gpath, _cohort) in GOLDENS.items():
        G, ON, OFF, A = (
            golden_rows(gpath),
            fresh_rows("on", tag),
            fresh_rows("off", tag),
            _audit(tag),
        )
        owed = (OUT / f"owed_{tag}.txt").read_text().split()
        if set(OFF) != set(owed):
            sys.exit(
                f"ABORT: {tag} control ran {sorted(set(OFF) ^ set(owed))} differently from owed"
            )
        new, reasons = [], {}
        for g in G:
            m = g[0]
            if m not in owed:
                new.append(g)
                continue
            on, off = ON[m], OFF[m]
            if not on[1] or not on[2]:
                sys.exit(f"ABORT: {m} has no fresh sha_out ({on[-1][:80]}). Decide it by hand.")
            off_ok = gate_verdict(g, off)[0] != "MISMATCH"
            moved = (on[1], on[2]) != (off[1], off[2])
            reason = "LEVER" if off_ok else ("STALE+LEVER" if moved else "STALE")
            if off_ok and not moved:
                sys.exit(
                    f"ABORT: {m}: `off` reproduces the golden AND equals `on`, yet `on` "
                    "mismatched. The run is not deterministic; nothing here can be trusted."
                )
            healed = A[m]["class"] != "SENTINEL" and A[m]["shipped"].get("sha") == g[1]
            reasons[m] = (reason, 2 if on[1] != g[1] else 3, healed)
            new.append(on[:6] + g[6:])
        report.append((tag, gpath, G, new, reasons))
        print(
            f"\n=== {tag}: {len(reasons)} of {len(G)} rows re-frozen   "
            f"{dict(Counter(r for r, _f, _h in reasons.values()))}"
        )
        for m, (reason, field, _h) in sorted(reasons.items()):
            print(f"    {m:22s} {reason:12s} first gated field that moved: {field}")

    if not args.write:
        print("\n(dry run -- pass --write to rewrite the goldens)")
        return
    for tag, gpath, _G, new, reasons in report:
        old_comments = [
            ln
            for ln in gpath.read_text().splitlines()
            if ln.startswith("#") and not ln.startswith(("# MANIFEST_SHA256=", "#DONE"))
        ]
        block = Path(args.comment_file).read_text().rstrip("\n").splitlines()
        block += [
            f"#   {m}\t{reason}\tfield {field}"
            for m, (reason, field, _h) in sorted(reasons.items())
        ]
        data = ["\t".join(r) for r in new]
        manifest = "\n".join(data)
        digest = hashlib.sha256(manifest.encode()).hexdigest()
        text = "\n".join(block + old_comments + data + [f"# MANIFEST_SHA256={digest}"])
        gpath.write_text(text + f"\n#DONE {len(data)}\n")
        print(f"wrote {gpath.name}: {len(data)} rows, MANIFEST_SHA256={digest}")
    rows = [
        f"{tag}\t{m}\t{reason}\t{field}\t" + "\t".join(next(r for r in new if r[0] == m)[:6])
        for tag, _p, _G, new, reasons in report
        for m, (reason, field, _h) in sorted(reasons.items())
    ]
    (OUT / "arm2_refreeze_rows.tsv").write_text(
        "# golden\tmolecule\treason\tfirst_moved_field\tfields 1-6 of the re-frozen row\n"
        + "\n".join(rows)
        + f"\n#DONE {len(rows)}\n"
    )


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("diff")
    sp = sub.add_parser("splice")
    sp.add_argument("--write", action="store_true")
    sp.add_argument("--comment-file", help="the '# v0.4.17: ...' block to put at the top")
    args = ap.parse_args()
    {"diff": cmd_diff, "splice": cmd_splice}[args.cmd](args)


if __name__ == "__main__":
    main()
