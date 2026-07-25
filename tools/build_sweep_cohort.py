"""Freeze a reproducible round-trip sweep cohort as a deduped symlink dir + manifest.

WHY A SYMLINK DIR AND NOT `--limit`
===================================
``test_dataset_roundtrip.py`` walks a dataset dir and keys every report by the file's
**basename**. The tmCAT-tmPHOTO tree has 26,230 ``.xyz`` files but only **25,197 unique
basenames** -- 1,033 names exist in BOTH ``cat/`` and ``photo/``. Pointing the harness at
the raw tree makes those 1,033 double-match and race each other's report writes. So a
cohort is materialized as a directory of symlinks, one per selected molecule, pointing at
the exact copy that was selected. This is the same pattern as the v0.4.4
``regression_inputs/`` cohort.

(For the record: ``run_regression_sweep.sh`` says "239 names"; that is stale. Two
independent counts -- ``uniq -d`` over basenames, and 26,230 - 25,197 -- both give 1,033.)

WHY NOT `--random`
==================
``test_dataset_roundtrip.py --random`` calls ``random.seed()`` with **system time**, so
two arms would not see the same molecules. Any A/B needs a frozen list. This tool seeds
explicitly (default 42, the project-wide convention) and writes the selected names to a
manifest so the cohort is auditable and re-creatable.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/build_sweep_cohort.py \
        --n 5000 --seed 42 --out tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k

Then point the sweep at ``--dataset-dir <out>``.
"""

import argparse
import json
import os
import random
import subprocess
import sys
from datetime import datetime

DEFAULT_DATASET = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "tmCAT-tmPHOTO_xyz_dataset")
)


def collect_unique(dataset_dir, subdirs):
    """Map basename -> absolute path, deduped deterministically.

    Walks ``subdirs`` in the given order and keeps the FIRST occurrence of each basename,
    so the dedup choice is a stated rule (``cat/`` wins over ``photo/`` by default) rather
    than filesystem order. Skips harness output (``*_generated.xyz``).
    """
    chosen: dict[str, str] = {}
    dupes: dict[str, list[str]] = {}
    for sub in subdirs:
        root = os.path.join(dataset_dir, sub)
        if not os.path.isdir(root):
            sys.exit(f"error: {root} not found")
        for dirpath, _dirnames, filenames in sorted(os.walk(root)):
            for fn in sorted(filenames):
                if not fn.endswith(".xyz") or fn.endswith("_generated.xyz"):
                    continue
                full = os.path.join(dirpath, fn)
                if fn in chosen:
                    dupes.setdefault(fn, [chosen[fn]]).append(full)
                else:
                    chosen[fn] = full
    return chosen, dupes


def main():
    ap = argparse.ArgumentParser(description="Freeze a sweep cohort as a symlink dir + manifest.")
    ap.add_argument("--n", type=int, required=True, help="Cohort size")
    ap.add_argument("--seed", type=int, default=42, help="Sampling seed (default 42)")
    ap.add_argument("--dataset-dir", default=DEFAULT_DATASET)
    ap.add_argument("--subdirs", default="cat,photo", help="Comma-separated, dedup priority order")
    ap.add_argument("--out", required=True, help="Cohort symlink dir to create")
    ap.add_argument(
        "--overlap-with",
        action="append",
        default=[],
        help="Results dir whose molecule set to report overlap against (repeatable)",
    )
    args = ap.parse_args()

    dataset_dir = os.path.abspath(args.dataset_dir)
    out = os.path.abspath(args.out)
    subdirs = [s.strip() for s in args.subdirs.split(",") if s.strip()]

    chosen, dupes = collect_unique(dataset_dir, subdirs)
    names = sorted(chosen)
    print(f"unique basenames: {len(names)}  (duplicated across subdirs: {len(dupes)})")

    if args.n > len(names):
        sys.exit(f"error: --n {args.n} exceeds {len(names)} available")

    # Sample from the SORTED name list so the draw depends only on (names, seed).
    picked = sorted(random.Random(args.seed).sample(names, args.n))

    if os.path.exists(out):
        sys.exit(f"error: {out} already exists -- refusing to overwrite a frozen cohort")
    os.makedirs(out)
    for name in picked:
        os.symlink(chosen[name], os.path.join(out, name))

    # Overlap against prior sweeps: the only slice that is identical-molecule diffable.
    overlaps = {}
    for res in args.overlap_with:
        indiv = os.path.join(os.path.abspath(res), "individual_reports")
        if not os.path.isdir(indiv):
            print(f"warning: {indiv} not found, skipping overlap", file=sys.stderr)
            continue
        prior = {fn[: -len(".json")] + ".xyz" for fn in os.listdir(indiv) if fn.endswith(".json")}
        shared = sorted(set(picked) & prior)
        overlaps[os.path.basename(os.path.abspath(res))] = {
            "prior_molecules": len(prior),
            "overlap": len(shared),
            "molecules": shared,
        }
        print(f"overlap with {os.path.basename(res)}: {len(shared)} / {len(picked)}")

    try:
        commit = (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                cwd=os.path.dirname(os.path.abspath(__file__)),
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
    except Exception:
        commit = "unknown"

    manifest = {
        "built_at": datetime.now().isoformat(timespec="seconds"),
        "commit_id": commit,
        "dataset_dir": dataset_dir,
        "subdirs_dedup_priority": subdirs,
        "unique_basenames_available": len(names),
        "duplicated_basenames": len(dupes),
        "seed": args.seed,
        "n": args.n,
        "cohort_dir": out,
        "overlaps": overlaps,
        "molecules": picked,
    }
    manifest_path = os.path.join(out, "..", f"{os.path.basename(out)}_manifest.json")
    manifest_path = os.path.abspath(manifest_path)
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"\ncohort: {len(picked)} symlinks in {out}")
    print(f"manifest: {manifest_path}")


if __name__ == "__main__":
    main()
