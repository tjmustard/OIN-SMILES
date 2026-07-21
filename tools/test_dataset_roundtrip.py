import argparse
import json
import logging
import multiprocessing
import os
import queue as queue_mod
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from datetime import datetime

logger = logging.getLogger(__name__)

# Progress marker a pass-1 child sends the moment the encode succeeds, so a later
# SIGKILL during generation does not lose smiles_1.
_ENCODED = "__encoded__"


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
from rmsd_utils import calculate_tmc_rmsd_detailed
from verify_roundtrip import (
    canonical_roundtrip_key,
    normalize_oin_for_comparison,
    read_atom_count,
)

from oinsmiles import XYZToSMILES
from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen as OIN3DGenerator
from oinsmiles.generator3d.ml_optimizer import resolve_xtb_binary

# Environment fields stamped into every report alongside commit_id, so each row
# in summary_roundtrip.json can be attributed to the code + env that produced it.
# Populated once in main(); merged into reports by save_artifacts().
RUN_ENV = {}

# Coordination-sphere mean-RMSD pass threshold (Angstrom). A single named constant
# so the gate value is legible in both the gate check and the per-report run-env
# provenance -- the value is unchanged (this is not a threshold change), it is just
# no longer a bare literal buried at the gate. See the FF-floor note at the gate.
RMSD_GATE = 1.0


def _build_run_env(args) -> dict:
    from rdkit import rdBase

    xtb_available = resolve_xtb_binary() is not None
    return {
        "rdkit_version": rdBase.rdkitVersion,
        "quick": bool(args.quick),
        # The per-molecule wall-clock cap and the RMSD gate the row was judged
        # against. Stamped so a failed row is self-describing: a high_rmsd under
        # optimizer_effective="ff" is an FF-floor row (string already matched,
        # geometry only just over the gate), and a timeout under quick=True with a
        # small mol_timeout is a budget artifact -- neither is a chemistry defect.
        "mol_timeout": args.mol_timeout,
        "rmsd_gate": RMSD_GATE,
        "xtb_available": xtb_available,
        # FF-only is the honest default when the g-xTB binary is absent; PASS 2
        # then runs an FF re-roll rather than any semi-empirical optimization.
        "optimizer_effective": "g-xtb" if xtb_available else "ff",
    }


def _pass2_config(xtb_available: bool):
    """Return (optimizer, tier1_name, tier5_name) for the PASS-2 recovery ladder.

    PASS 2 recovers PASS-1 soft failures two ways that are independent of g-xTB:
    a fresh stochastic re-embed and a wider ensemble (1 -> 5). When the 'xtb'
    binary is absent the g-xTB optimizer is a no-op, so we run FF-only and name
    the tiers for what actually happened rather than mislabeling them g-xTB.
    """
    if xtb_available:
        return "g-xtb", "g-xTB_1", "g-xTB_5"
    return None, "FF_reroll_1", "FF_reroll_5"


def _honesty_breakdown(reports):
    """Split failed reports into honest buckets so artifacts stop reading as defects.

    Returns a dict of counts. The two artifact buckets are the ones S7 makes legible:

    * ``ff_floor_high_rmsd`` -- the string round-trip matched and only the geometric
      tightness gate failed, under FF-only (no ``xtb``). A chemically-correct
      round-trip, not an accuracy defect.
    * ``quick_timeout`` -- a ``TimeoutException`` under ``--quick`` (the 30 s
      ``--mol-timeout`` budget), i.e. a budget artifact, not a chemistry hang.

    Everything else failed is a ``real_failure`` (string mismatch, atom-count,
    encode crash, generation exception) plus non-quick ``timeout`` kept separate.
    Purely descriptive: it changes no status and no gate.
    """
    counts = {
        "ff_floor_high_rmsd": 0,
        "quick_timeout": 0,
        "timeout_full_budget": 0,
        "real_failure": 0,
    }
    for r in reports:
        if r.get("status") != "failed":
            continue
        err = r.get("error") or ""
        if err.startswith("High RMSD") and (
            r.get("ff_floor") or r.get("optimizer_effective") == "ff"
        ):
            counts["ff_floor_high_rmsd"] += 1
        elif "TimeoutException" in err:
            counts["quick_timeout" if r.get("quick") else "timeout_full_budget"] += 1
        else:
            counts["real_failure"] += 1
    return counts


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
            rmsd, reason = calculate_tmc_rmsd_detailed(
                mol_orig, mol_gen_xyz, mol2_bonded=mol_gen_bonded
            )
            # A sphere the metric cannot map is not bad geometry. Reporting it as
            # "High RMSD: 996.0" made 62 chemically-correct round-trips read as the
            # worst structures in the dataset.
            if rmsd is None:
                report["metrics"]["rmsd"] = None
                report["metrics"]["rmsd_mapping_reason"] = reason
                report["error"] = f"RMSD mapping failed at {tier_name}: {reason}"
                return False, last_gen_xyz_content
            report["metrics"]["rmsd"] = round(rmsd, 4)
            if rmsd >= RMSD_GATE:
                report["error"] = f"High RMSD at {tier_name}: {rmsd:.4f}"
                # The string round-trip already matched above, so this is a
                # chemically-correct round-trip that only missed the geometric
                # tightness gate. Under FF-only that is the FF floor, not a defect;
                # flag it so triage/backlog stop conflating it with real failures.
                if RUN_ENV.get("optimizer_effective") == "ff":
                    report["ff_floor"] = True
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


def _encode_and_attempt_in_child(result_queue, tier_name, gen_kwargs, xyz_path, report):
    """Child entry point for pass 1: encode XYZ -> OIN, then run one tier.

    The encode lives in here, not in the parent, because it is where UGUHAH_comp_0
    actually wedges -- ``XYZToSMILES.convert`` on that 97-atom porphyrinoid runs for
    minutes. A watchdog wrapped only around generation would never have caught it.
    """
    oin1_string = None
    try:
        oin1_string = XYZToSMILES().convert(xyz_path)
        report["smiles_1"] = oin1_string
    except BaseException as e:  # noqa: BLE001 - the parent must learn about anything fatal
        report["status"] = "failed"
        report["error"] = f"XYZToSMILES failed: {type(e).__name__}: {e}\n{traceback.format_exc()}"
        result_queue.put((False, report, None, None))
        return

    # Hand the encode back immediately: if generation is later SIGKILLed, the parent
    # still knows smiles_1, and can tell "hung while encoding" from "hung while
    # generating".
    result_queue.put((_ENCODED, oin1_string))

    try:
        generator = OIN3DGenerator(**gen_kwargs)
        success, last_xyz = _attempt_generation(tier_name, generator, oin1_string, xyz_path, report)
        result_queue.put((success, report, last_xyz, oin1_string))
    except BaseException as e:  # noqa: BLE001
        report["error"] = (
            f"Generation/Verification failed at {tier_name}: {type(e).__name__}: {e}\n"
            f"{traceback.format_exc()}"
        )
        result_queue.put((False, report, None, oin1_string))


def _attempt_in_child(result_queue, tier_name, gen_kwargs, oin1_string, xyz_path, report):
    """Child entry point for pass 2: build a generator, run one tier, ship the report back.

    The generator is constructed here rather than passed in, because a spawned child
    cannot inherit it.
    """
    try:
        generator = OIN3DGenerator(**gen_kwargs)
        success, last_xyz = _attempt_generation(tier_name, generator, oin1_string, xyz_path, report)
        result_queue.put((success, report, last_xyz))
    except BaseException as e:  # noqa: BLE001
        report["error"] = (
            f"Generation/Verification failed at {tier_name}: {type(e).__name__}: {e}\n"
            f"{traceback.format_exc()}"
        )
        result_queue.put((False, report, None))


def _supervise(proc, result_queue, timeout, tier_name, report, progress=None, stage="generating"):
    """Run ``proc`` under a hard wall-clock cap, SIGKILLing it on expiry.

    ``signal.alarm`` cannot interrupt a hang inside native C++ or a tight loop that
    never yields to the interpreter: the signal is queued until control returns to
    Python, which never happens. UGUHAH_comp_0 wedged a Phase-0 shard for 35+ minutes
    despite a 420 s cap. Only an OS kill clears that.

    Args:
        progress: Optional dict; ``_ENCODED`` markers from the child are recorded here
            rather than mistaken for the final result.
        stage: What the child starts out doing, named in the timeout message.

    Returns:
        The child's payload tuple, or None if it timed out or died. On None, ``report``
        carries the reason and ``progress`` whatever the child managed to send first.
    """
    proc.start()

    # Drain the queue before joining: a child blocks on exit until its pipe is read.
    payload = None
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            message = result_queue.get(timeout=0.5)
        except queue_mod.Empty:
            if not proc.is_alive():
                break  # died without reporting -- segfault, OOM kill, ...
            continue
        if message and message[0] == _ENCODED:
            if progress is not None:
                progress["oin1"] = message[1]
            stage = "generating"
            continue
        payload = message
        break

    if payload is not None:
        proc.join(timeout=30)
        if proc.is_alive():
            proc.kill()
            proc.join()
        return payload

    report["status"] = "failed"
    if proc.is_alive():
        proc.kill()
        proc.join()
        report["error"] = (
            f"TimeoutException at {tier_name}: exceeded {timeout}s while {stage} (hard kill)"
        )
    else:
        report["error"] = (
            f"Generation/Verification failed at {tier_name}: "
            f"child process died with exit code {proc.exitcode} and no result"
        )
    return None


def _encode_and_attempt(tier_name, generator, gen_kwargs, xyz_path, report, timeout, converter):
    """Pass 1 for one molecule: encode then generate, under the watchdog if asked.

    Spawning costs a few seconds (rdkit + MetalloGen re-import), so an untimed run stays
    in-process and pays nothing.

    Returns:
        tuple: (success_bool, last_gen_xyz_content, oin1_string_or_None)
    """
    if timeout <= 0:
        try:
            oin1_string = converter.convert(xyz_path)
            report["smiles_1"] = oin1_string
        except Exception as e:
            report["status"] = "failed"
            report["error"] = (
                f"XYZToSMILES failed: {type(e).__name__}: {e}\n{traceback.format_exc()}"
            )
            return False, None, None
        success, last_xyz = _attempt_generation(tier_name, generator, oin1_string, xyz_path, report)
        return success, last_xyz, oin1_string

    ctx = multiprocessing.get_context("spawn")
    result_queue = ctx.Queue()
    proc = ctx.Process(
        target=_encode_and_attempt_in_child,
        args=(result_queue, tier_name, gen_kwargs, xyz_path, report),
    )
    progress = {}
    payload = _supervise(proc, result_queue, timeout, tier_name, report, progress, stage="encoding")
    if payload is None:
        # A kill during generation still leaves us the encode the child sent back.
        oin1_string = progress.get("oin1")
        report["smiles_1"] = oin1_string
        return False, None, oin1_string

    success, child_report, last_xyz, oin1_string = payload
    report.clear()
    report.update(child_report)
    return success, last_xyz, oin1_string


def _run_attempt(tier_name, generator, gen_kwargs, oin1_string, xyz_path, report, timeout):
    """Pass 2 for one molecule: run one generator tier, under the watchdog if asked.

    Returns:
        tuple: (success_bool, last_gen_xyz_content)
    """
    if timeout <= 0:
        return _attempt_generation(tier_name, generator, oin1_string, xyz_path, report)

    ctx = multiprocessing.get_context("spawn")
    result_queue = ctx.Queue()
    proc = ctx.Process(
        target=_attempt_in_child,
        args=(result_queue, tier_name, gen_kwargs, oin1_string, xyz_path, report),
    )
    payload = _supervise(proc, result_queue, timeout, tier_name, report)
    if payload is None:
        return False, None

    success, child_report, last_xyz = payload
    report.clear()
    report.update(child_report)
    return success, last_xyz


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
        help="Run with a 30-second timeout for g-xTB and limited UFF pool size",
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
        help="Hard wall-clock timeout in seconds per molecule per tier (0 to disable). "
        "Non-zero runs each attempt in a subprocess that is SIGKILLed on expiry, so a "
        "hang inside native code cannot wedge the run. Costs a few seconds per tier.",
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
        print(
            f"Watchdog: each molecule runs in a killable subprocess, {args.mol_timeout}s cap "
            "(costs a few seconds per tier to re-import rdkit + MetalloGen)."
        )

    global_report = []
    requires_g_xtb = []

    commit_id = _get_git_commit_id()
    print(f"Git commit: {commit_id}")
    RUN_ENV.update(_build_run_env(args))
    print(
        f"Env: rdkit {RUN_ENV['rdkit_version']}, xtb_available={RUN_ENV['xtb_available']}, "
        f"optimizer_effective={RUN_ENV['optimizer_effective']}"
    )

    xyz_to_smiles = XYZToSMILES()

    # Determine quick settings
    timeout_val = 30 if args.quick else 300
    ff_params_fast = {"uff_pool_size": 2, "max_attempts": 10} if args.quick else None

    print("\n--- PASS 1: UFF FAST-PASS ---")
    uff_kwargs = {
        "optimizer": None,
        "ensemble_size": 1,
        "timeout": timeout_val,
        "ff_params": ff_params_fast,
    }
    gen_uff = OIN3DGenerator(**uff_kwargs)

    for i, xyz_path in enumerate(xyz_files, 1):
        basename = os.path.splitext(os.path.basename(xyz_path))[0]
        print(f"[{i}/{len(xyz_files)}] UFF Pass: {basename}...", end=" ", flush=True)

        report = initialize_report(xyz_path, commit_id=commit_id)

        # The encode runs inside the watchdog too: UGUHAH_comp_0 hangs in
        # XYZToSMILES.convert(), not in the generator.
        t0 = time.monotonic()
        success, last_xyz, oin1_string = _encode_and_attempt(
            "UFF_1", gen_uff, uff_kwargs, xyz_path, report, args.mol_timeout, xyz_to_smiles
        )
        # Stamp wall-clock spent here. Set *after* the call so the subprocess path's
        # report.clear()/update(child_report) cannot wipe it. PASS 2 (if reached) adds
        # its own tier time to this figure.
        report.setdefault("metrics", {})["elapsed_s"] = round(time.monotonic() - t0, 3)

        if oin1_string is None and not success:
            save_artifacts(report, None, output_dir, is_final=True)
            global_report.append(report)
            err = report.get("error", "")
            print("FAILED (1D conversion)" if err.startswith("XYZToSMILES") else "FAILED (hard)")
            continue

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
        # With no 'xtb' binary the g-xTB optimizer is a no-op, so PASS 2 runs an
        # FF re-roll (fresh re-embed + ensemble 1 -> 5) and is named honestly.
        xtb_available = RUN_ENV.get("xtb_available", resolve_xtb_binary() is not None)
        opt2, tier1, tier5 = _pass2_config(xtb_available)
        print(f"\n--- PASS 2: {tier1}/{tier5} PASS ({len(requires_g_xtb)} files) ---")
        pass2_1_kwargs = {
            "optimizer": opt2,
            "ensemble_size": 1,
            "timeout": timeout_val,
            "ff_params": ff_params_fast,
        }
        pass2_5_kwargs = {**pass2_1_kwargs, "ensemble_size": 5}
        gen_pass2_1 = OIN3DGenerator(**pass2_1_kwargs)
        gen_pass2_5 = OIN3DGenerator(**pass2_5_kwargs)

        for i, (xyz_path, oin1_string, report) in enumerate(requires_g_xtb, 1):
            basename = report["molecule"]
            print(f"[{i}/{len(requires_g_xtb)}] PASS 2: {basename}...", flush=True)

            # Carry PASS-1 wall-clock forward, then accumulate each tier's time onto it,
            # so elapsed_s ends up as the molecule's total across both passes. Read now,
            # before any _run_attempt subprocess-path report.clear() can drop it.
            elapsed_s = report.get("metrics", {}).get("elapsed_s", 0.0)

            # Attempt tier 1 (ensemble_size=1)
            print(f"  -> Trying {tier1}...", end=" ", flush=True)
            t0 = time.monotonic()
            success, last_xyz = _run_attempt(
                tier1,
                gen_pass2_1,
                pass2_1_kwargs,
                oin1_string,
                xyz_path,
                report,
                args.mol_timeout,
            )
            elapsed_s += time.monotonic() - t0
            report.setdefault("metrics", {})["elapsed_s"] = round(elapsed_s, 3)

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
                print(f"FAILED (Hard failure, skipping {tier5})")
                continue

            print("FAILED (Soft failure)")

            # Attempt tier 5 (ensemble_size=5)
            print(f"  -> Trying {tier5}...", end=" ", flush=True)
            t0 = time.monotonic()
            success, last_xyz = _run_attempt(
                tier5,
                gen_pass2_5,
                pass2_5_kwargs,
                oin1_string,
                xyz_path,
                report,
                args.mol_timeout,
            )
            elapsed_s += time.monotonic() - t0
            report.setdefault("metrics", {})["elapsed_s"] = round(elapsed_s, 3)

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
    print(f"Optimizer effective: {RUN_ENV.get('optimizer_effective', 'unknown')}")

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
    # Honest failure split: FF-floor high_rmsd and quick-timeouts are budget/FF
    # artifacts, not accuracy defects -- surfaced here so the backlog stops
    # conflating them with real failures.
    bd = _honesty_breakdown(final_report)
    print(
        "Failure breakdown: "
        f"{bd['real_failure']} real, "
        f"{bd['ff_floor_high_rmsd']} FF-floor high_rmsd (string matched, geometry-only), "
        f"{bd['quick_timeout']} quick-timeout, "
        f"{bd['timeout_full_budget']} full-budget timeout."
    )
    print(f"Global report saved to {global_path}")


if __name__ == "__main__":
    main()
