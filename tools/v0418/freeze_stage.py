"""Stage the L2 eta lane for ``/freeze-measurements``: scrub, THEN gzip, one prefix. Reads only.

Same discipline as ``tools/census/stage_reattribution.py``, for the same reasons: the harvester
copies a ``.gz`` VERBATIM into a PUBLIC tree, so text is scrubbed before it is compressed and the
decompressed bytes are checked afterwards; gzip is written with ``mtime=0`` so a re-stage is
byte-identical; every file carries one prefix (``v0418_l2_``).

THE ACCEPTANCE TEST of this freeze is ``--verify <measurements/v0.4.18-l2>``: both headline pairs
(+195/-29 and +189/-39 for M1+M2; +203/-51 and +179/-59 for all three) re-derived from the frozen
files and the frozen census table ALONE, no results directory on the path.

Two files are MADE here rather than copied, because the claim they carry lives in 3,438 structure
files and 3,438 harness reports that are far too large to freeze:

    ab_rows.tsv.gz              molecule x {off, on, on3}: sha256 of the generated structure (the
                                dead-lever check), the harness's ``coordination.intact``, elapsed_s;
                                plus the sweep of record's structure sha256 (the noise floor)
    eta_distance_audit.jsonl.gz the audit's rows with floats rounded to 3 dp -- 582 KB -> under the
                                harvester's 512 KB cap. The .txt beside it is the unrounded run.

NOT staged: the probes' ``.jsonl`` (each record embeds a whole generated xyz); their ``.txt`` is.

    <main>/.venv/bin/python tools/v0418/freeze_stage.py            # writes <results>/freeze/
    <main>/.venv/bin/python tools/v0418/freeze_stage.py --verify <main>/measurements/v0.4.18-l2
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import sys
from pathlib import Path

DATA = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
R = DATA / "results-v0.4.18-eta-detached"
RECORD = DATA / "results-v0.4.17-sweep"
STAGE = R / "freeze"
PER_FILE_CAP = 512 * 1024
P = "v0418_l2_"
ARMS = ("off", "on", "on3")
EXPECTED = {  # arm -> (self gains, self losses, verified gains, verified losses)
    "on": (195, 29, 189, 39),
    "on3": (203, 51, 179, 59),
}
ISO = ("ISO", "ISO_MARGINAL", "ISO_CLASH")
STRING_FAULTS = ("P_E1_COVERAGE", "DATA_MULTI", "P_DETACHED", "E1_GRAPH", "E1_HCOUNT")
BAD_STEREO = ("MIRROR", "MIRROR_PARTIAL", "DIFFERENT")

_CHECKOUT = re.compile(r"/home/[^/\s\"']+/Documents/GitHub/([A-Za-z0-9._-]+)")
_HOME = re.compile(r"/home/[^/\s\"']+/")
_SCRATCH = re.compile(r"/tmp/claude-\d+/[^\s\"')]*")
_TMPFILE = re.compile(r"/tmp/tmp[A-Za-z0-9_]+")

TEXT = [
    "eta_distance_identity.txt",
    "eta_distance_audit.txt",
    "eta_path_probe.txt",
    "eta_target_probe.txt",
    "eta_exempt_probe.txt",
    "eta_scoped_probe.txt",
    "eta_rescoped_probe.txt",
    "noneta_identity_probe.txt",
    "eta_all3_probe.txt",
    "probe3_classes.json",
    "probe4_classes.json",
    "eta_ab_report.txt",
    "eta_ab_report.json",
    "eta_ab_report_on3.txt",
    "eta_ab_report_on3.json",
]


def scrub(text: str) -> str:
    text = _CHECKOUT.sub(r"<checkout:\1>", text)
    text = _SCRATCH.sub("<SCRATCH>", text)
    text = _TMPFILE.sub("<TMPFILE>", text)
    return _HOME.sub("<HOME>/", text)


def _sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""


def _round(o):
    if isinstance(o, float):
        return round(o, 3)
    if isinstance(o, list):
        return [_round(x) for x in o]
    if isinstance(o, dict):
        return {k: _round(v) for k, v in o.items()}
    return o


def _audit_rounded() -> str:
    lines = (R / "eta_distance_audit.jsonl").read_text().splitlines()
    if not lines[-1].startswith("#DONE"):
        sys.exit("ABORT: eta_distance_audit.jsonl has no #DONE trailer")
    out = [json.dumps(_round(json.loads(ln)), separators=(",", ":")) for ln in lines[:-1]]
    return "\n".join(out + [lines[-1]]) + "\n"


def _ab_rows() -> str:
    cohort = sorted(p.stem for p in (DATA / "cohort-v0.4.18-eta").glob("*.xyz"))
    if len(cohort) != 1146:
        sys.exit(f"ABORT: cohort has {len(cohort)} molecules, not 1146")
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t", lineterminator="\n")
    head = ["molecule", "sha256_record"]
    for a in ARMS:
        head += [f"sha256_{a}", f"intact_{a}", f"elapsed_s_{a}"]
    w.writerow(head)
    for m in cohort:
        row = [m, _sha(RECORD / "structures" / f"{m}_generated.xyz")]
        for a in ARMS:
            rep = json.loads((R / f"ab_{a}" / "individual_reports" / f"{m}.json").read_text())
            row += [
                _sha(R / f"ab_{a}" / "structures" / f"{m}_generated.xyz"),
                int(bool((rep.get("coordination") or {}).get("intact"))),
                (rep.get("metrics") or {}).get("elapsed_s", ""),
            ]
        w.writerow(row)
    return buf.getvalue() + f"#DONE {len(cohort)}\n"


def _write(name: str, raw: str, gz: bool):
    clean = scrub(raw)
    if "/home/" in clean or "/tmp/claude" in clean:
        sys.exit(f"ABORT: a local path survives scrubbing in {name}")
    dest = STAGE / name
    if gz:
        with open(dest, "wb") as raw_fh:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, mtime=0) as fh:
                fh.write(clean.encode())
        back = gzip.open(dest, "rt").read()  # check what was WRITTEN, not what was intended
        if back != clean or "/home/" in back or "/tmp/claude" in back:
            sys.exit(f"ABORT: {dest} does not decompress to the scrubbed text")
    else:
        dest.write_text(clean)
    size = dest.stat().st_size
    print(f"  {size:8d}  {name}{'   (scrubbed)' if clean != raw else ''}")
    if size > PER_FILE_CAP:
        sys.exit(f"ABORT: {name} is {size} bytes -- over the harvester's cap, it would be SKIPPED")
    return size


def stage():
    if STAGE.exists():
        shutil.rmtree(STAGE)  # a stale staged file would be harvested as if it were current
    STAGE.mkdir(parents=True)
    n = total = 0
    for f in TEXT:
        if not (R / f).exists():
            sys.exit(f"ABORT: {R / f} is missing -- refusing to stage a partial release")
        total += _write(P + f, (R / f).read_text(), False)
        n += 1
    for a in ARMS:
        d = R / f"ab_{a}"
        gv = (d / "g_verdict.jsonl").read_text()
        if not gv.rstrip().splitlines()[-1].startswith("#DONE"):
            sys.exit(f"ABORT: {d}/g_verdict.jsonl has no #DONE trailer")
        total += _write(f"{P}ab_{a}_g_verdict.jsonl.gz", gv, True)
        total += _write(
            f"{P}ab_{a}_bucket_report_honest.json.gz",
            (d / "bucket_report_honest.json").read_text(),
            True,
        )
        n += 2
    total += _write(P + "ab_rows.tsv.gz", _ab_rows(), True)
    total += _write(P + "eta_distance_audit.jsonl.gz", _audit_rounded(), True)
    commits = "".join(
        f"{a}\t{(c.read_text().strip() if c.exists() else 'see off')}\n"
        for a in ARMS
        for c in [R / f"ab_{a}" / "AB_COMMIT" if a == "on3" else R / "AB_COMMIT"]
    )
    total += _write(P + "ab_commits.tsv", "arm\tlaunched_from\n" + commits, False)
    print(f"staged {n + 3} files, {total / 1024:.0f} KB -> {STAGE}")


def verify(frozen: Path, expected=None, noise=(1063, 1063)):
    """Both headline pairs from the FROZEN files alone. No results directory is read."""
    expected = EXPECTED if expected is None else expected
    census = frozen.parent / "v0.4.18-census" / "v0418_census_attribution_table.tsv.gz"
    T = {r["molecule"]: r for r in csv.DictReader(gzip.open(census, "rt"), delimiter="\t")}

    def arm(a):
        B = {
            r["molecule"].removesuffix(".xyz"): r
            for r in json.load(gzip.open(frozen / f"{P}ab_{a}_bucket_report_honest.json.gz", "rt"))
        }
        G = {}
        for ln in gzip.open(frozen / f"{P}ab_{a}_g_verdict.jsonl.gz", "rt"):
            if not ln.startswith("#"):
                r = json.loads(ln)
                G[r["molecule"]] = r
        return B, G

    def self_ok(m, B, G):
        return B[m]["bucket"] == "byte_exact"

    def ver_ok(m, B, G):
        if not self_ok(m, B, G) or T[m]["fault"] in STRING_FAULTS:
            return False
        g = G.get(m)
        return bool(g) and g.get("graph") in ISO and g.get("stereo") not in BAD_STEREO

    Bf, Gf = arm("off")
    ok = True
    for a, want in expected.items():
        Bn, Gn = arm(a)
        got = []
        for fn in (self_ok, ver_ok):
            off = {m: fn(m, Bf, Gf) for m in Bf}
            on = {m: fn(m, Bn, Gn) for m in Bf}
            got += [sum(on[m] and not off[m] for m in Bf), sum(off[m] and not on[m] for m in Bf)]
        flag = "OK" if tuple(got) == want else "MISMATCH"
        ok &= tuple(got) == want
        print(
            f"  arm {a:3s} n={len(Bf)}  self +{got[0]}/-{got[1]}  verified +{got[2]}/-{got[3]}"
            f"   expected {want}   {flag}"
        )
    rows = list(csv.DictReader(gzip.open(frozen / f"{P}ab_rows.tsv.gz", "rt"), delimiter="\t"))
    rows = [r for r in rows if not r["molecule"].startswith("#")]
    cmp_ = [r for r in rows if r["sha256_record"] and r["sha256_off"]]
    same = sum(r["sha256_record"] == r["sha256_off"] for r in cmp_)
    print(
        f"  noise floor: OFF == sweep of record on {same}/{len(cmp_)} structures "
        f"(expected {noise[0]}/{noise[1]})"
    )
    ok &= (same, len(cmp_)) == tuple(noise)
    if not ok:
        sys.exit("ABORT: the frozen tree does NOT reproduce the headline")
    print("  the frozen tree reproduces the headline")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--verify", type=Path, help="a frozen release dir; re-derive the headline")
    args = ap.parse_args()
    verify(args.verify) if args.verify else stage()


if __name__ == "__main__":
    main()
