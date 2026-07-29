#!/usr/bin/env python
"""Freeze a release's comparison numbers into the tracked ``measurements/`` tree.

WHY THIS EXISTS
===============
``tmCAT-tmPHOTO_xyz_dataset/`` is **gitignored in its entirety** (``.gitignore:85``) -- the
2.1 GB dataset and every ``results-*/`` sweep. ``spec/handoffs/`` is gitignored too. So for six
consecutive releases the project wrote excellent prose about numbers that lived on exactly one
disk, and:

    🔴 v0.4.11's mirror-audit JSON -- the 19/250 measurement that chartered the whole of
       v0.4.12 -- is GONE. It was recoverable only because the TOOL was committed and
       deterministic, at ~50 minutes of CPU to re-run.

``CLOSEOUT.md`` §4 mandated git-durable *prose* and said nothing about *data*, which is exactly
how the lapse went unnoticed: v0.4.3 committed 2 data files, v0.4.6 committed 5, and
v0.4.7-v0.4.12 committed **zero**.

This tool is the other half of that ritual step. It is meant to be run automatically from
CLOSEOUT §4b and manually via ``/freeze-measurements``.

⚠ THE TREE IS PUBLIC. ``origin`` is a public GitHub repo and a sibling session pushes ``main``
regularly, so anything harvested here becomes public and permanent. Two consequences: local
absolute paths must never leak (checked below), and size discipline is not cosmetic -- hence
the caps.

WHERE IT WRITES, AND WHY NOT ``$PWD``
=====================================
Always ``<main checkout>/measurements/``, resolved from ``git rev-parse --git-common-dir``,
**never** the current worktree. A worktree's files vanish on ``git worktree remove`` -- which is
how this session nearly lost its own data, having created and then removed
``../oin-v0412-release``. Untracked files also do not survive ``git clean -fd``. The main
checkout plus a commit is the only durable destination.

Same resolution trick v0.4.9 used to stop ``gate_v047.sh`` silently selecting a sibling
project's venv.

SELECTION IS AN ALLOWLIST, NOT A SIZE SWEEP
===========================================
Measured: a naive "every small .json/.md under a cap" rule selects **12,756 files / 19 MB** of
per-molecule probe dumps and case registries. The allowlist below yields **29 files / 838 KB**
across the whole historical tree -- the artifacts a *later release actually diffs*.

Usage
=====
    python tools/harvest_measurements.py --backfill --dry-run
    python tools/harvest_measurements.py --release v0.4.12 --from <dir> [--from <dir2>]
    python tools/harvest_measurements.py --link          # read-access symlink in a worktree
"""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

#: Filename patterns worth preserving forever. Each is a *comparison* artifact: a later release
#: diffs it to learn whether it regressed. Report-backing evidence for a single write-up belongs
#: in ``docs/agentic-notes/<release>/`` instead -- see ``measurements/README.md`` for the test.
ALLOW = [
    "bucket_report*.md",
    "FROZEN.md",
    "SOURCE",
    "MANIFEST*",
    "VALIDATION.md",
    "README.md",
    "triage_hard_fails.md",
    "transitions.json",
    "*_DONE",
    "*TREND.tsv",
    "CASE_REGISTRY.md",
    # v0.4.12 onward: the per-release instruments' own output, harvested from a scratchpad.
    "mirror_audit*.json",
    "transition_*.json",
    "ab_*.json",
    "cohort_*.json",
    "audit_*.json",
    # v0.4.13. Added because this release's instruments are EXACTLY the artifact class whose
    # loss motivated this tool -- the two mirror arms are the same kind of file as v0.4.11's
    # vanished 19/250 audit, and `fold_transition_veto.json` names all 171 molecules behind the
    # +3.42 headline, which no prose in `docs/` reproduces.
    # ⚠ `mirror_*` deliberately, NOT `mirror_arm*`. The narrower pattern was written first and
    # silently dropped `mirror_cat_{promoted,noveto}.json` -- the arms that reproduce v0.4.12's
    # PUBLISHED gate (19 -> 0), i.e. half the evidence the v0.4.13 promotion rests on. Caught only
    # because the dry-run's file list was read line by line against the scratchpad. A too-narrow
    # allowlist fails exactly like a broken instrument: it prints a plausible total.
    "mirror_*.json",
    "fold_*.json",
    "attach_class_audit.json",
    "prefilter_*.json",
]

#: Directories that are raw inputs or bulk per-molecule output. Never harvested.
PRUNE_DIRS = {
    "individual_reports",
    "structures",
    "cat",
    "photo",
    "rebaseline_inputs",
    "regression_inputs",
    "reports",
    "__pycache__",
    ".git",
}

#: Extensions/names that are never comparison artifacts regardless of the allowlist.
DENY = ["*.log", "*.xyz", "*.oin", "*.bak", "*.pre-rebuild-bak", "worker_pids.txt", "*.pyc"]

PER_FILE_CAP = 512 * 1024
TOTAL_CAP = 5 * 1024 * 1024

#: Provenance for names this project's tools produce. Anything unmatched is reported
#: ``UNKNOWN`` rather than given a plausible-looking guess -- a figure without its command is an
#: order of magnitude, not a measurement.
PROVENANCE = [
    (
        r"^mirror_audit_seed(\d+)_veto_(on|off)\.json$",
        "tools/mirror_audit_donor_fold.py --dataset <cat> --n 250 --seed {0}"
        "   (OIN_FOLD_PARITY_VETO={1})",
    ),
    (
        r"^transition_(fold|veto)\.json$",
        "tools/fold_transition_sim.py --sweep <frozen sweep> --arm {0}",
    ),
    (
        r"^ab_.*\.json$",
        "tools/ab_accept_scored.py --cohort <cohort> --lever <lever> --timeout 150 --hard-cap 240",
    ),
    (r"^cohort_.*\.json$", "cohort manifest, built from the frozen sweep"),
    # v0.4.13's instruments.
    (
        r"^mirror_arm([AB])_(\w+)\.json$",
        "tools/mirror_audit_donor_fold.py --dataset <cohort-v0.4.5-5k> --n 250 --seed 7"
        "   (arm {0}: {1}; the noveto arm sets OIN_FOLD_PARITY_VETO=0) -- mixed cat+photo draw,"
        " reads 33 collapses -> 0",
    ),
    (
        r"^mirror_cat_(\w+)\.json$",
        "tools/mirror_audit_donor_fold.py --dataset <cat> --n 250 --seed 7"
        "   ({0}; the noveto arm sets OIN_FOLD_PARITY_VETO=0) -- CAT-ONLY draw, reproduces"
        " v0.4.12's published 19 -> 0 with achiral unmoved at 157",
    ),
    (
        r"^fold_transition_(fold|veto)\.json$",
        "tools/fold_transition_sim.py --sweep <frozen sweep> --arm {0}"
        "   --dataset <cat> --dataset <photo>   (ABSOLUTE roots: the relative default"
        " silently excludes every mover from a worktree)",
    ),
    (
        r"^fold_key_invariance\.json$",
        "tools/fold_key_invariance.py --sweep <frozen sweep>"
        "   (does the fold ever change the round-trip KEY? 0 => generator-neutral"
        " => an offline re-score is exact and no sweep is owed)",
    ),
    (
        r"^attach_class_audit\.json$",
        "tools/attach_class_audit.py --results-dir <frozen sweep>"
        "   (MEDZUR/GAVSED split, with the byte_exact control arm)",
    ),
    (
        r"^prefilter_.*\.json$",
        "tools/prefilter_prevalence.py --xyz <input> | --cohort <dir>"
        "   (OIN_PREFILTER_ADVISORY two-arm; needs a QUIET machine — it reports a latency cost)",
    ),
    (r"^bucket_report.*\.md$", "tools/roundtrip_bucket_report.py --results-dir <dir>"),
    (r"^FROZEN\.md$|^SOURCE$", "hand-written provenance record for a frozen sweep"),
    (r"^triage_hard_fails\.md$", "tools/triage_hard_fails.py"),
    (r"^CASE_REGISTRY\.md$|^.*TREND\.tsv$", "tools/rebuild_summary.py / milestone_report.py"),
]


def repo_root() -> Path:
    """The MAIN checkout, resolved from anywhere -- including inside a linked worktree.

    ``--git-common-dir`` is the shared ``.git`` directory; from a worktree it is an absolute
    path into the main checkout, and from the main checkout it is the relative ``.git``. Its
    parent is the main working tree in both cases.
    """
    common = subprocess.check_output(["git", "rev-parse", "--git-common-dir"], text=True).strip()
    return Path(common).resolve().parent


def in_worktree() -> bool:
    git_dir = subprocess.check_output(["git", "rev-parse", "--git-dir"], text=True).strip()
    return "worktrees" in Path(git_dir).resolve().parts


def _wanted(name: str) -> bool:
    if any(fnmatch.fnmatch(name, d) for d in DENY):
        return False
    return any(fnmatch.fnmatch(name, a) for a in ALLOW)


def collect(src: Path) -> list[tuple[int, Path]]:
    """Allowlisted files under ``src``, each within the per-file cap."""
    out: list[tuple[int, Path]] = []
    for dirpath, dirnames, files in os.walk(src):
        dirnames[:] = [d for d in dirnames if d not in PRUNE_DIRS]
        for f in files:
            if not _wanted(f):
                continue
            p = Path(dirpath) / f
            try:
                size = p.stat().st_size
            except OSError:
                continue
            if size <= PER_FILE_CAP:
                out.append((size, p))
    return sorted(out, key=lambda t: str(t[1]))


def provenance_for(name: str) -> str:
    for pattern, template in PROVENANCE:
        m = re.match(pattern, name)
        if m:
            try:
                return template.format(*m.groups())
            except (IndexError, KeyError):
                return template
    return "UNKNOWN"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def scrub(text: str, root: Path) -> tuple[str, bool]:
    """Rewrite local absolute paths to repo-relative form. Returns ``(text, still_leaks)``.

    ⚠ This REWRITES rather than refuses, and the distinction matters. The tools that produced
    these artifacts wrote absolute paths into their own headers, so a refuse-on-`/home/` rule
    rejects exactly the files worth keeping -- measured: it dropped every `bucket_report_*.md`
    and the v0.4.8 `SOURCE`, i.e. the frozen baseline this whole tree exists to preserve, while
    happily keeping the `ARM_A_DONE` sentinels. Losing the data to protect the path is the wrong
    trade when the path is trivially removable.

    ``<repo>/`` is stripped (the reference stays meaningful, since the tree lives in the repo)
    and any surviving home directory is masked. Anything still matching afterwards is a genuine
    leak and the caller drops the file -- so the guarantee is kept, just not by throwing the
    baby out.
    """
    out = text.replace(str(root) + "/", "").replace(str(root), ".")
    out = re.sub(r"/home/[^/\s]+/", "<HOME>/", out)
    out = re.sub(r"/tmp/claude-\d+/[^\s\"')]*", "<SCRATCH>", out)
    still = "/home/" in out or "/tmp/claude" in out
    return out, still


#: Renames applied on the way in. The v0.4.12 scratchpad called the same kind of run
#: ``audit_baseline_seed7`` and ``audit_base_seed11``, and "baseline" is ambiguous anyway --
#: mirror_audit_donor_fold.py runs BOTH fold arms internally, so what differed between those
#: runs was the ambient veto, not the fold.
RENAME = {
    "audit_baseline_seed7.json": "mirror_audit_seed07_veto_off.json",
    "audit_veto_seed7.json": "mirror_audit_seed07_veto_on.json",
    "audit_base_seed11.json": "mirror_audit_seed11_veto_off.json",
    "audit_veto_seed11.json": "mirror_audit_seed11_veto_on.json",
    "transition_fold_baseline.json": "transition_fold.json",
    "ab_pilot_eta.json": "ab_eta_accept_stale_cohort.json",
    "ab_pilot2_eta.json": "ab_eta_accept_realpop.json",
    "cohort_pilot.json": "cohort_pilot_stale.json",
}

#: Scratch intermediates: real files, but not artifacts anyone will diff.
SKIP_NAMES = {"movers6.json", "movers40.json"}


def write_release(
    dest: Path, release: str, picks: list[tuple[int, Path, str]], dry: bool, root: Path
):
    rel_dir = dest / release
    rows = []
    for size, src, newname in picks:
        rows.append((newname, size, src))
        if dry:
            continue
        rel_dir.mkdir(parents=True, exist_ok=True)
        try:
            text = src.read_text(errors="strict")
        except (OSError, UnicodeDecodeError):
            shutil.copy2(src, rel_dir / newname)  # binary/unreadable: copy verbatim
            continue
        cleaned, _still = scrub(text, root)
        (rel_dir / newname).write_text(cleaned)

    if dry:
        return rows

    lines = [
        f"# `measurements/{release}` — frozen comparison artifacts",
        "",
        "Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.",
        "",
        "| file | bytes | sha256 | produced by |",
        "|---|---:|---|---|",
    ]
    for name, size, src in sorted(rows):
        lines.append(f"| `{name}` | {size} | `{sha256(rel_dir / name)}` | {provenance_for(name)} |")
    # ⚠ The generated index is itself published, so it gets the same scrubbing as the harvested
    # files. Missing this leaked the full scratchpad path (`/tmp/claude-<uid>/-home-<user>-...`)
    # into a public repo on the first run -- the guard has to cover what the guard writes.
    src_lines, _ = scrub("\n".join(f"{n}  <-  {s}" for n, _sz, s in sorted(rows)), root)
    lines += ["", "Source paths at harvest time:", "", "```", src_lines, "```"]
    (rel_dir / "README.md").write_text("\n".join(lines) + "\n")
    return rows


def infer_release(dirname: str) -> str:
    """``results-v0.4.8-honest`` -> ``v0.4.8-honest``; ``results-capstone-v042`` -> ``capstone-v042``.

    The ``results-`` prefix is always stripped -- it is noise inside a tree whose every entry is
    a result. Names that carry no version at all keep their remaining text rather than being
    forced into a version-shaped slot they do not fit.
    """
    return re.sub(r"^results-", "", dirname)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", help="target release folder, e.g. v0.4.12")
    ap.add_argument(
        "--from", dest="sources", action="append", default=[], help="source directory (repeatable)"
    )
    ap.add_argument(
        "--backfill",
        action="store_true",
        help="harvest every tmCAT-tmPHOTO_xyz_dataset/results-* dir",
    )
    ap.add_argument("--dry-run", action="store_true", help="print the selection, write nothing")
    ap.add_argument(
        "--link", action="store_true", help="create .measurements-main in this worktree and exit"
    )
    ap.add_argument("--force", action="store_true", help="proceed past the total-size cap")
    args = ap.parse_args()

    root = repo_root()
    dest = root / "measurements"
    print(f"main checkout : {root}")
    print(f"destination   : {dest}")
    if in_worktree():
        print("⚠ running inside a WORKTREE — writing to the main checkout above, not here.")

    if args.link:
        link = Path.cwd() / ".measurements-main"
        if link.is_symlink() or link.exists():
            print(f"already present: {link}")
            return 0
        link.symlink_to(dest)
        print(f"linked {link} -> {dest}")
        return 0

    jobs: list[tuple[str, Path]] = []
    if args.backfill:
        ds = root / "tmCAT-tmPHOTO_xyz_dataset"
        for d in sorted(ds.glob("results-*")):
            if d.is_dir():
                jobs.append((infer_release(d.name), d))
    if args.sources:
        if not args.release:
            ap.error("--from requires --release")
        for s in args.sources:
            jobs.append((args.release, Path(s).expanduser()))
    if not jobs:
        ap.error("nothing to do: pass --backfill and/or --release with --from")

    grand_total, grand_files, unresolved = 0, 0, []
    for release, src in jobs:
        if not src.is_dir():
            print(f"  ⚠ missing source, skipped: {src}")
            continue
        picks = []
        for size, p in collect(src):
            if p.name in SKIP_NAMES:
                continue
            try:
                text = p.read_text(errors="strict")
            except (OSError, UnicodeDecodeError):
                picks.append((size, p, RENAME.get(p.name, p.name)))
                continue
            _cleaned, still = scrub(text, root)
            if still:
                # Scrubbing could not make it safe -- that IS a refusal, and a rare one.
                print(f"  🔴 REFUSED (local path survives scrubbing, tree is public): {p}")
                continue
            picks.append((size, p, RENAME.get(p.name, p.name)))
        if not picks:
            continue
        total = sum(s for s, _p, _n in picks)
        grand_total += total
        grand_files += len(picks)
        print(f"\n{release}: {len(picks)} files, {total / 1024:.0f} KB   <- {src}")
        for size, p, newname in picks:
            tag = "" if provenance_for(newname) != "UNKNOWN" else "   [provenance UNKNOWN]"
            rename = f"  (was {p.name})" if newname != p.name else ""
            print(f"    {size:>8}  {newname}{rename}{tag}")
            if provenance_for(newname) == "UNKNOWN":
                unresolved.append(f"{release}/{newname}")
        write_release(dest, release, picks, args.dry_run, root)

    print(f"\nTOTAL: {grand_files} files, {grand_total / 1024:.0f} KB")
    if unresolved:
        print(
            f"⚠ {len(unresolved)} file(s) with UNKNOWN provenance — recorded as such, not guessed."
        )
    if grand_total > TOTAL_CAP and not args.force:
        print(
            f"🔴 REFUSED: {grand_total / 1024 / 1024:.1f} MB exceeds the {TOTAL_CAP // 1024 // 1024} MB cap."
        )
        print("   measurements/ is a comparison tree in a PUBLIC repo, not a data dump.")
        print("   Narrow the selection, or pass --force if this is genuinely warranted.")
        return 1
    if args.dry_run:
        print("(dry run — nothing written)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
