import argparse
import json
import logging
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import traceback
from datetime import datetime


class TimeoutException(Exception):
    pass


def timeout_handler(signum, frame):
    raise TimeoutException("Molecule processing timed out")


logger = logging.getLogger(__name__)


def _get_git_commit_id() -> str:
    """Return the short git commit hash, with '-dirty' suffix if the working tree has changes.

    Returns:
        A string like 'a1b2c3d' or 'a1b2c3d-dirty'. Returns 'unknown' if git
        is unavailable or the directory is not a git repository.
    """
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
        # Check for uncommitted changes
        dirty = subprocess.call(
            ["git", "diff", "--quiet", "HEAD"],
            cwd=repo_root,
            stderr=subprocess.DEVNULL,
        )
        if dirty != 0:
            commit += "-dirty"
        return commit
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        logger.warning("Could not determine git commit ID.")
        return "unknown"


# Add src and tests/integration to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../tests/integration")))

from rdkit import Chem
from rmsd_utils import calculate_tmc_rmsd
from verify_roundtrip import (
    canonical_roundtrip_key,
    normalize_oin_for_comparison,
    read_atom_count,
)

from oinsmiles import XYZToSMILES
from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen as OIN3DGenerator

# Environment fields stamped into every report alongside commit_id, so each row
# in summary_roundtrip.json can be attributed to the code + env that produced it.
# Populated once in main(); merged into reports by save_artifacts().
RUN_ENV = {}


def _build_run_env(args) -> dict:
    from rdkit import rdBase

    return {
        "rdkit_version": rdBase.rdkitVersion,
        "quick": bool(args.quick),
        "xtb_available": shutil.which("xtb") is not None,
    }


def _attempt_generation(tier_name, generator, oin1_string, xyz_path, report):
    """Attempt 1D -> 3D -> 1D for a single generator config. Returns (success_bool, last_gen_xyz_content)"""
    tmp_dir = tempfile.mkdtemp()
    gen_xyz_path = os.path.join(tmp_dir, "gen.xyz")
    last_gen_xyz_content = None
    xyz_to_smiles = XYZToSMILES()

    try:
        # OIN(1) -> XYZ(Gen)
        gen_result = generator.generate(oin1_string)

        with open(gen_xyz_path, "w") as f:
            f.write(gen_result.xyz)

        last_gen_xyz_content = gen_result.xyz
        mol_gen_bonded = gen_result.mol

        # XYZ(Gen) -> OIN(2)
        if mol_gen_bonded is not None:
            try:
                from oinsmiles.utils.xyz2mol import get_oin_string

                with open(gen_xyz_path, "r") as f:
                    xyz_lines = f.readlines()
                natoms = int(xyz_lines[0].strip())
                xyz_coords = []
                for i in range(2, 2 + natoms):
                    parts = xyz_lines[i].split()
                    xyz_coords.append([float(x) for x in parts[1:4]])
                import numpy as np

                xyz_coords = np.array(xyz_coords)
                oin2_string = get_oin_string(mol_gen_bonded, xyz_coords)
            except Exception:
                oin2_string = xyz_to_smiles.convert(gen_xyz_path)
        else:
            oin2_string = xyz_to_smiles.convert(gen_xyz_path)

        report["smiles_2"] = oin2_string

        # Verification: compare by structure-level canonical key (collapses
        # chemically-meaningless notation drift -- implicit-H, carbene, symmetric
        # donor, fragment order -- while still catching genuinely different
        # connectivity, metal/geometry, or eta winding). The normalized strings
        # are kept only for the human-readable diagnostic message.
        s1 = normalize_oin_for_comparison(oin1_string.strip())
        s2 = normalize_oin_for_comparison(oin2_string.strip())

        if canonical_roundtrip_key(oin1_string) != canonical_roundtrip_key(oin2_string):
            report["error"] = f"String mismatch at {tier_name}. Exp: {s1}, Got: {s2}"
            return False, last_gen_xyz_content

        # Geometric fidelity
        mol_orig = Chem.MolFromXYZFile(xyz_path)
        mol_gen_xyz = Chem.MolFromXYZFile(gen_xyz_path)
        if mol_gen_bonded is None:
            mol_gen_bonded = mol_gen_xyz

        if mol_orig and mol_gen_xyz:
            rmsd = calculate_tmc_rmsd(mol_orig, mol_gen_xyz, mol2_bonded=mol_gen_bonded)
            report["metrics"]["rmsd"] = round(rmsd, 4)
            if rmsd >= 1.0:
                report["error"] = f"High RMSD at {tier_name}: {rmsd:.4f}"
                return False, last_gen_xyz_content

        # Atom count
        atom_count_input = read_atom_count(xyz_path)
        atom_count_generated = read_atom_count(gen_xyz_path)
        if atom_count_input != atom_count_generated:
            report["error"] = (
                f"Atom count mismatch at {tier_name}. Input {atom_count_input} != Gen {atom_count_generated}"
            )
            return False, last_gen_xyz_content

        # If we got here, it's a success
        report["status"] = "success"
        report["tier_passed"] = tier_name
        report["error"] = None
        return True, last_gen_xyz_content

    except Exception as e:
        report["error"] = (
            f"Generation/Verification failed at {tier_name}: {type(e).__name__}: {e}\n{traceback.format_exc()}"
        )
        return False, last_gen_xyz_content
    finally:
        shutil.rmtree(tmp_dir)


def initialize_report(xyz_path: str, commit_id: str = "unknown") -> dict:
    """Create a fresh report dict for a molecule.

    Args:
        xyz_path: Absolute path to the input XYZ file.
        commit_id: Git commit hash to stamp into the report.

    Returns:
        A dict with all standard report fields.
    """
    basename = os.path.splitext(os.path.basename(xyz_path))[0]
    return {
        "molecule": basename,
        "input_xyz": xyz_path,
        "commit_id": commit_id,
        "status": "pending",
        "tier_passed": None,
        "metrics": {},
        "smiles_1": None,
        "smiles_2": None,
        "error": None,
    }


def save_artifacts(report, last_xyz, output_dir, is_final=False):
    basename = report["molecule"]

    report.update(RUN_ENV)
    report["saved_at"] = datetime.now().isoformat(timespec="seconds")

    # Save individual JSON
    indiv_path = os.path.join(output_dir, "individual_reports", f"{basename}.json")
    with open(indiv_path, "w") as f:
        json.dump(report, f, indent=2)

    # A molecule that now passes must not keep stale failure forensics around
    # (they read as current failures long after the underlying bug is fixed).
    if report["status"] == "success":
        stale_fail_dir = os.path.join(output_dir, "test_failures", basename)
        if os.path.isdir(stale_fail_dir):
            shutil.rmtree(stale_fail_dir, ignore_errors=True)

    # Save structures and OINs for inspection if successful or if it's the final pass
    if report["status"] == "success" or is_final:
        struct_dir = os.path.join(output_dir, "structures")
        os.makedirs(struct_dir, exist_ok=True)

        if report["smiles_2"] is not None:
            with open(os.path.join(struct_dir, f"{basename}.oin"), "w") as f:
                f.write(report["smiles_2"])
        elif report["smiles_1"] is not None:
            with open(os.path.join(struct_dir, f"{basename}.oin"), "w") as f:
                f.write(report["smiles_1"])

        if last_xyz is not None:
            with open(os.path.join(struct_dir, f"{basename}_generated.xyz"), "w") as f:
                f.write(last_xyz)

    # Save failure forensics if it's a total failure
    if report["status"] == "failed" and is_final:
        fail_dir = os.path.join(output_dir, "test_failures", basename)
        os.makedirs(fail_dir, exist_ok=True)
        with open(os.path.join(fail_dir, "report.json"), "w") as f:
            json.dump(report, f, indent=2)
        if last_xyz is not None:
            with open(os.path.join(fail_dir, "last_generated.xyz"), "w") as f:
                f.write(last_xyz)


def main():
    parser = argparse.ArgumentParser(
        description="Dataset Roundtrip Testing (Multi-Pass Architecture)"
    )
    parser.add_argument(
        "--dataset-dir",
        type=str,
        default="../tmCAT-tmPHOTO_xyz_dataset",
        help="Path to dataset directory",
    )
    parser.add_argument(
        "--output-dir", type=str, default="../results", help="Directory to save results"
    )
    parser.add_argument("--limit", type=int, default=None, help="Limit number of molecules to test")
    parser.add_argument("--cpu", action="store_true", help="Force CPU execution")
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run with a 60-second timeout for g-xTB and limited UFF pool size",
    )
    parser.add_argument(
        "--continue",
        dest="continue_run",
        action="store_true",
        help="Continue from previous run (skip already processed molecules in summary_roundtrip.json and append new results)",
    )
    parser.add_argument(
        "--rerun-failed",
        action="store_true",
        help="Only run on molecules that previously failed (requires existing summary_roundtrip.json)",
    )
    parser.add_argument(
        "--only",
        type=str,
        default=None,
        help="Comma-separated molecule names (e.g. ABAFOZ_comp_0) to process, bypassing "
        "--continue/--rerun-failed filtering. For targeted repro of specific cases.",
    )
    parser.add_argument(
        "--shard",
        type=str,
        default=None,
        help="I:N — process only the I-th of N deterministic slices of the (filtered, "
        "sorted) molecule list, e.g. --shard 2:5. For parallel workers.",
    )
    parser.add_argument(
        "--no-summary",
        action="store_true",
        help="Write individual reports only; skip rewriting summary_roundtrip.json "
        "(parallel workers use this, then tools/rebuild_summary.py merges).",
    )
    parser.add_argument(
        "--random",
        action="store_true",
        help="Randomly shuffle the dataset before selecting molecules (useful with --limit)",
    )
    parser.add_argument(
        "--mol-timeout",
        type=int,
        default=0,
        help="Global timeout in seconds for processing a single molecule (0 to disable).",
    )
    args = parser.parse_args()

    if args.cpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    dataset_dir = os.path.abspath(args.dataset_dir)
    output_dir = os.path.abspath(args.output_dir)

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "test_failures"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "individual_reports"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "structures"), exist_ok=True)

    xyz_files = []
    output_dir_abs = os.path.abspath(output_dir)
    for root, dirs, files in os.walk(dataset_dir):
        # Prevent recursing into the output directory if it's nested inside the dataset directory
        dirs[:] = [
            d for d in dirs if not os.path.abspath(os.path.join(root, d)).startswith(output_dir_abs)
        ]

        for f in files:
            if f.endswith(".xyz") and not f.endswith("_generated.xyz"):
                xyz_files.append(os.path.join(root, f))

    # Sort for deterministic order
    xyz_files = sorted(xyz_files)

    old_report = []
    summary_path = os.path.join(output_dir, "summary_roundtrip.json")
    # Always load an existing summary: the final write merges old + new by
    # molecule (newest wins), so a targeted run never clobbers other rows.
    if os.path.exists(summary_path):
        with open(summary_path, "r") as f:
            old_report = json.load(f)

    # Filter the molecule list
    if args.only:
        wanted = {n.strip().removesuffix(".xyz") for n in args.only.split(",") if n.strip()}
        xyz_files = [f for f in xyz_files if os.path.splitext(os.path.basename(f))[0] in wanted]
        print(f"Only mode: matched {len(xyz_files)} of {len(wanted)} requested molecules.")
    elif args.rerun_failed:
        if old_report:
            failed_mols = {r["molecule"] for r in old_report if r["status"] == "failed"}
            xyz_files = [
                f for f in xyz_files if os.path.splitext(os.path.basename(f))[0] in failed_mols
            ]
            print(f"Rerun-failed: filtered to {len(xyz_files)} previously failed molecules.")
        else:
            print(f"Warning: --rerun-failed specified but {summary_path} not found. Running all.")
    elif args.continue_run:
        if old_report:
            processed_mols = {r["molecule"] for r in old_report}
            xyz_files = [
                f
                for f in xyz_files
                if os.path.splitext(os.path.basename(f))[0] not in processed_mols
            ]
            print(f"Continue mode: skipping {len(processed_mols)} already processed molecules.")
        else:
            print(f"Note: --continue specified but {summary_path} not found. Starting fresh.")

    if args.shard:
        try:
            shard_i, shard_n = (int(x) for x in args.shard.split(":"))
        except ValueError:
            parser.error(f"--shard must be I:N (got {args.shard!r})")
        if not (1 <= shard_i <= shard_n):
            parser.error(f"--shard index out of range: {args.shard!r}")
        xyz_files = xyz_files[shard_i - 1 :: shard_n]
        print(f"Shard {shard_i}/{shard_n}: {len(xyz_files)} molecules in this slice.")

    if args.random:
        import random

        random.seed()  # Use system time
        random.shuffle(xyz_files)
        print("Randomly shuffling the dataset order.")

    if args.limit:
        xyz_files = xyz_files[: args.limit]

    print(f"Found {len(xyz_files)} XYZ files to process.")

    if args.mol_timeout > 0:
        signal.signal(signal.SIGALRM, timeout_handler)

    global_report = []
    requires_g_xtb = []

    commit_id = _get_git_commit_id()
    print(f"Git commit: {commit_id}")
    RUN_ENV.update(_build_run_env(args))
    print(f"Env: rdkit {RUN_ENV['rdkit_version']}, xtb_available={RUN_ENV['xtb_available']}")

    xyz_to_smiles = XYZToSMILES()

    # Determine quick settings
    timeout_val = 60 if args.quick else 300
    ff_params_fast = {"uff_pool_size": 2, "max_attempts": 10} if args.quick else None

    print("\n--- PASS 1: UFF FAST-PASS ---")
    gen_uff = OIN3DGenerator(
        optimizer=None, ensemble_size=1, timeout=timeout_val, ff_params=ff_params_fast
    )

    for i, xyz_path in enumerate(xyz_files, 1):
        basename = os.path.splitext(os.path.basename(xyz_path))[0]
        print(f"[{i}/{len(xyz_files)}] UFF Pass: {basename}...", end=" ", flush=True)

        report = initialize_report(xyz_path, commit_id=commit_id)

        try:
            oin1_string = xyz_to_smiles.convert(xyz_path)
            report["smiles_1"] = oin1_string
        except Exception as e:
            report["status"] = "failed"
            report["error"] = (
                f"XYZToSMILES failed: {type(e).__name__}: {e}\n{traceback.format_exc()}"
            )
            save_artifacts(report, None, output_dir, is_final=True)
            global_report.append(report)
            print("FAILED (1D conversion)")
            continue

        try:
            if args.mol_timeout > 0:
                signal.alarm(args.mol_timeout)
            success, last_xyz = _attempt_generation("UFF_1", gen_uff, oin1_string, xyz_path, report)
        except TimeoutException as e:
            report["status"] = "failed"
            report["error"] = f"TimeoutException: {e}"
            success, last_xyz = False, None
        finally:
            if args.mol_timeout > 0:
                signal.alarm(0)

        if success:
            save_artifacts(report, last_xyz, output_dir, is_final=True)
            global_report.append(report)
            print("SUCCESS")
        else:
            err = report.get("error", "")
            if err and (
                err.startswith("Generation/Verification failed at")
                or err.startswith("TimeoutException")
            ):
                report["status"] = "failed"
                save_artifacts(report, last_xyz, output_dir, is_final=True)
                global_report.append(report)
                print("FAILED (Hard failure, skipping g-xTB)")
            else:
                report["status"] = "pending_g-xtb"
                save_artifacts(report, last_xyz, output_dir, is_final=False)
                requires_g_xtb.append((xyz_path, oin1_string, report))
                print("FAILED (queued for g-xTB)")

    if requires_g_xtb:
        print(f"\n--- PASS 2: g-xTB PASS ({len(requires_g_xtb)} files) ---")
        gen_g_xtb_1 = OIN3DGenerator(
            optimizer="g-xtb", ensemble_size=1, timeout=timeout_val, ff_params=ff_params_fast
        )
        gen_g_xtb_5 = OIN3DGenerator(
            optimizer="g-xtb", ensemble_size=5, timeout=timeout_val, ff_params=ff_params_fast
        )

        for i, (xyz_path, oin1_string, report) in enumerate(requires_g_xtb, 1):
            basename = report["molecule"]
            print(f"[{i}/{len(requires_g_xtb)}] g-xTB Pass: {basename}...", flush=True)

            # Attempt g-xTB_1
            print("  -> Trying g-xTB_1...", end=" ", flush=True)
            try:
                if args.mol_timeout > 0:
                    signal.alarm(args.mol_timeout)
                success, last_xyz = _attempt_generation(
                    "g-xTB_1", gen_g_xtb_1, oin1_string, xyz_path, report
                )
            except TimeoutException as e:
                report["status"] = "failed"
                report["error"] = f"TimeoutException at g-xTB_1: {e}"
                success, last_xyz = False, None
            finally:
                if args.mol_timeout > 0:
                    signal.alarm(0)

            if success:
                save_artifacts(report, last_xyz, output_dir, is_final=True)
                global_report.append(report)
                print("SUCCESS")
                continue

            err = report.get("error", "")
            if err and (
                err.startswith("Generation/Verification failed at")
                or err.startswith("TimeoutException")
            ):
                report["status"] = "failed"
                save_artifacts(report, last_xyz, output_dir, is_final=True)
                global_report.append(report)
                print("FAILED (Hard failure, skipping g-xTB_5)")
                continue

            print("FAILED (Soft failure)")

            # Attempt g-xTB_5
            print("  -> Trying g-xTB_5...", end=" ", flush=True)
            try:
                if args.mol_timeout > 0:
                    signal.alarm(args.mol_timeout)
                success, last_xyz = _attempt_generation(
                    "g-xTB_5", gen_g_xtb_5, oin1_string, xyz_path, report
                )
            except TimeoutException as e:
                report["status"] = "failed"
                report["error"] = f"TimeoutException at g-xTB_5: {e}"
                success, last_xyz = False, None
            finally:
                if args.mol_timeout > 0:
                    signal.alarm(0)

            if success:
                print("SUCCESS")
            else:
                report["status"] = "failed"
                print("FAILED (Total Failure)")

            save_artifacts(report, last_xyz, output_dir, is_final=True)
            global_report.append(report)

    # Save global report: merge old + new by molecule, newest wins. This both
    # fixes the old --rerun-failed duplicate-row bug (old failed row + new row
    # for the same molecule) and stops a targeted run from clobbering the rows
    # of molecules it did not process.
    run_successes = sum(1 for r in global_report if r["status"] == "success")
    print(f"\nFinished processing {len(xyz_files)} files.")
    print(f"This run: {run_successes} successes / {len(global_report)} processed.")

    if args.no_summary:
        print("--no-summary: individual reports written; summary left untouched.")
        print("Merge later with: uv run python tools/rebuild_summary.py --output-dir <dir>")
        return

    merged = {r["molecule"]: r for r in old_report}
    merged.update({r["molecule"]: r for r in global_report})
    final_report = sorted(merged.values(), key=lambda r: r.get("molecule", ""))

    global_path = os.path.join(output_dir, "summary_roundtrip.json")
    with open(global_path, "w") as f:
        json.dump(final_report, f, indent=2)

    successes = sum(1 for r in final_report if r["status"] == "success")
    print(f"Summary totals: {successes} successes, {len(final_report) - successes} failures.")
    print(f"Global report saved to {global_path}")


if __name__ == "__main__":
    main()
