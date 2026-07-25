"""Generator-free canonicality A/B: re-encode a completed sweep's stored geometry pairs.

WHAT THIS MEASURES -- and what it cannot
========================================
Every round-trip report from ``test_dataset_roundtrip.py`` records an absolute
``input_xyz`` and (on success) a stored ``structures/<mol>_generated.xyz``. That pair is
two *different geometries of the same isomer*: the crystal input and the generated
structure. Re-encoding BOTH with the current code and comparing the two OIN strings is
therefore a direct measurement of **encoder canonicality** -- exactly the quantity
v0.4.5 Lanes 1-3 move.

Because it re-encodes stored geometry it **never invokes the 3D generator**, so it is
immune to the confound that dominates the batch harness (``BASELINE.md`` §6): 67.4% of
round-trip failures are generation timeouts, so any code change that shifts runtime moves
the pass-rate for reasons unrelated to the notation. Here the geometry is frozen, both
arms see byte-identical inputs, and the only variable is the encoder.

The flip side: it is **blind to generator changes**. Lanes 4/5/6 add descriptors the
generator must reproduce; measure those with ``tools/injectivity/*`` instead, never here.

It is also a *lower bound* on the drift a real sweep sees, because it can only compare
molecules whose generated structure was stored -- rows that timed out or produced nothing
have no second geometry and are excluded (and counted).

CONTRAST WITH ``recalculate_oin_smiles.py``
==========================================
That tool does the same re-encode but **mutates the source reports in place** and applies
conservative geometry gates so a row is only promoted to success when RMSD/atom-count also
pass. This one is **non-mutating** (writes a fresh results dir) and deliberately drops the
geometry gates: every row with both strings present is stamped ``success`` so that
``roundtrip_bucket_report.classify()`` routes it purely on the STRINGS, into
``byte_exact`` / ``key_equal`` / ``facmer_divergent`` / ``structural``. Geometry quality is
not what this instrument is for.

Usage
=====
    CAP=tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042
    for i in 1 2 3 4 5 6; do
      PYTHONPATH=src .venv/bin/python tools/reencode_ab.py \
        --results-dir $CAP --out /path/armA --shard $i:6 &
    done; wait
    PYTHONPATH=src .venv/bin/python tools/roundtrip_bucket_report.py --results-dir /path/armA

Shards write disjoint files into one ``--out`` dir, so no merge step is needed --
``roundtrip_bucket_report.py`` just globs them. ``--shard`` is **1-based** (``1:6`` .. ``6:6``),
matching ``test_dataset_roundtrip.py``.
"""

import argparse
import contextlib
import glob
import json
import os
import subprocess
import sys
import time
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from oinsmiles import XYZToSMILES  # noqa: E402


@contextlib.contextmanager
def _silence_fds():
    """Redirect C-level stdout/stderr to devnull (openbabel prints distance warnings)."""
    with open(os.devnull, "w") as devnull:
        old_out, old_err = os.dup(1), os.dup(2)
        try:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            os.dup2(old_out, 1)
            os.dup2(old_err, 2)
            os.close(old_out)
            os.close(old_err)


def _find_generated_xyz(results_dir, basename):
    """Stored generated structure for one molecule, or None.

    Same two locations ``recalculate_oin_smiles.py:37`` looks in.
    """
    for candidate in (
        os.path.join(results_dir, "structures", f"{basename}_generated.xyz"),
        os.path.join(results_dir, "test_failures", basename, "last_generated.xyz"),
    ):
        if os.path.exists(candidate):
            return candidate
    return None


def _provenance(results_dir):
    from rdkit import rdBase

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    try:
        commit = (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                cwd=repo_root,
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
        if subprocess.call(
            ["git", "diff", "--quiet", "HEAD"], cwd=repo_root, stderr=subprocess.DEVNULL
        ):
            commit += "-dirty"
    except Exception:
        commit = "unknown"
    # Capture every OIN_* lever so an arm's configuration is recorded in its own output.
    levers = {k: v for k, v in sorted(os.environ.items()) if k.startswith("OIN_")}
    return {
        "reencode_at": datetime.now().isoformat(timespec="seconds"),
        "reencode_commit_id": commit,
        "reencode_rdkit_version": rdBase.rdkitVersion,
        "reencode_source_dir": os.path.abspath(results_dir),
        "reencode_levers": levers,
    }


def _parse_shard(spec):
    """``I:N`` 1-based -> (I, N). Raises ValueError on a malformed spec."""
    i_str, n_str = spec.split(":")
    i, n = int(i_str), int(n_str)
    if n < 1 or not (1 <= i <= n):
        raise ValueError(f"--shard {spec}: need 1 <= I <= N and N >= 1")
    return i, n


def main():
    ap = argparse.ArgumentParser(
        description="Non-mutating generator-free re-encode of a sweep's stored geometry pairs."
    )
    ap.add_argument(
        "--results-dir",
        required=True,
        help="Source results dir with individual_reports/ and structures/",
    )
    ap.add_argument("--out", required=True, help="Destination results dir (created if absent)")
    ap.add_argument("--shard", help="1-based shard spec I:N, e.g. 3:6")
    ap.add_argument("--limit", type=int, help="Process at most N reports (after --shard)")
    ap.add_argument("--only", help="Comma-separated molecule ids")
    ap.add_argument(
        "--skip-existing",
        action="store_true",
        help="Resume: skip molecules already written to --out (safe only within one arm)",
    )
    args = ap.parse_args()

    src = os.path.abspath(args.results_dir)
    out = os.path.abspath(args.out)
    indiv_src = os.path.join(src, "individual_reports")
    indiv_out = os.path.join(out, "individual_reports")

    if not os.path.isdir(indiv_src):
        sys.exit(f"error: {indiv_src} not found")
    os.makedirs(indiv_out, exist_ok=True)

    reports = sorted(glob.glob(os.path.join(indiv_src, "*.json")))
    if args.only:
        wanted = {m.strip().removesuffix(".xyz") for m in args.only.split(",") if m.strip()}
        reports = [p for p in reports if os.path.basename(p)[: -len(".json")] in wanted]
    if args.shard:
        i, n = _parse_shard(args.shard)
        reports = reports[i - 1 :: n]
    if args.limit:
        reports = reports[: args.limit]

    tag = args.shard or "all"
    prov = _provenance(src)
    conv = XYZToSMILES()

    n_written = 0
    skipped_no_input = 0
    skipped_no_generated = 0
    skipped_existing = 0
    t0 = time.time()

    for pos, path in enumerate(reports, 1):
        if args.skip_existing:
            stem = os.path.basename(path)[: -len(".json")]
            if os.path.exists(os.path.join(indiv_out, f"{stem}.json")):
                skipped_existing += 1
                continue
        try:
            with open(path) as f:
                rep = json.load(f)
        except Exception as e:  # noqa: BLE001
            print(f"[{tag}] warning: unreadable {path}: {e}", file=sys.stderr)
            continue

        basename = rep.get("molecule") or os.path.basename(path)[: -len(".json")]
        input_xyz = rep.get("input_xyz")
        gen_xyz = _find_generated_xyz(src, basename)

        # Exclusions are counted, never silently dropped: a molecule with no stored
        # generated structure has no second geometry, so it carries no canonicality signal.
        if not input_xyz or not os.path.exists(input_xyz):
            skipped_no_input += 1
            continue
        if gen_xyz is None:
            skipped_no_generated += 1
            continue

        s1 = s2 = None
        err = None
        try:
            with _silence_fds():
                s1 = conv.convert(input_xyz)
        except Exception as e:  # noqa: BLE001
            err = f"Input re-encode failed: {type(e).__name__}: {e}"
        if s1 is not None:
            try:
                with _silence_fds():
                    s2 = conv.convert(gen_xyz)
            except Exception as e:  # noqa: BLE001
                err = f"Generated re-encode failed: {type(e).__name__}: {e}"

        # Stamp `success` whenever both strings exist so that
        # roundtrip_bucket_report.classify() buckets on the STRINGS alone. A missing
        # smiles_1 lands in `encode_fail`; a missing smiles_2 lands in `hard_fail` via
        # classify's defensive branch -- both honest signals about the encoder.
        new = {
            "molecule": basename,
            "input_xyz": input_xyz,
            "generated_xyz": gen_xyz,
            "smiles_1": s1,
            "smiles_2": s2,
            "status": "success" if (s1 is not None and s2 is not None) else "failed",
            "error": err,
            # Inherited from the SOURCE sweep -- generation time, not re-encode time.
            # Kept so the bucket report's eta/elapsed percentile sections still populate.
            "metrics": rep.get("metrics"),
            "source_status": rep.get("status"),
            "source_smiles_1": rep.get("smiles_1"),
            "source_smiles_2": rep.get("smiles_2"),
            **prov,
        }
        with open(os.path.join(indiv_out, f"{basename}.json"), "w") as f:
            json.dump(new, f, indent=2)
        n_written += 1

        if pos % 200 == 0:
            rate = pos / max(time.time() - t0, 1e-9)
            print(f"[{tag}] {pos}/{len(reports)}  {rate:.1f} mol/s", flush=True)

    elapsed = time.time() - t0
    stats = {
        "shard": tag,
        "reports_seen": len(reports),
        "written": n_written,
        "skipped_no_input_xyz": skipped_no_input,
        "skipped_no_generated_xyz": skipped_no_generated,
        "skipped_existing": skipped_existing,
        "elapsed_s": round(elapsed, 1),
        **prov,
    }
    shard_name = tag.replace(":", "of")
    with open(os.path.join(out, f"reencode_stats_{shard_name}.json"), "w") as f:
        json.dump(stats, f, indent=2)

    print(
        f"[{tag}] done: {n_written} written, "
        f"{skipped_no_input} no-input, {skipped_no_generated} no-generated, "
        f"{skipped_existing} pre-existing, {elapsed:.0f}s"
    )


if __name__ == "__main__":
    main()
