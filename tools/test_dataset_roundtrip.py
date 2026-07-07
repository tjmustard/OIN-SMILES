import argparse
import json
import os
import shutil
import sys
import tempfile
import traceback

# Add src and tests/integration to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../tests/integration")))

from rdkit import Chem
from rmsd_utils import calculate_tmc_rmsd
from verify_roundtrip import normalize_oin_for_comparison, read_atom_count

from oinsmiles import XYZToSMILES
from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen as OIN3DGenerator


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

        # Verification
        s1 = normalize_oin_for_comparison(oin1_string.strip())
        s2 = normalize_oin_for_comparison(oin2_string.strip())

        if s1 != s2:
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


def initialize_report(xyz_path):
    basename = os.path.splitext(os.path.basename(xyz_path))[0]
    return {
        "molecule": basename,
        "input_xyz": xyz_path,
        "status": "pending",
        "tier_passed": None,
        "metrics": {},
        "smiles_1": None,
        "smiles_2": None,
        "error": None,
    }


def save_artifacts(report, last_xyz, output_dir, is_final=False):
    basename = report["molecule"]

    # Save individual JSON
    indiv_path = os.path.join(output_dir, "individual_reports", f"{basename}.json")
    with open(indiv_path, "w") as f:
        json.dump(report, f, indent=2)

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
    for root, _, files in os.walk(dataset_dir):
        for f in files:
            if f.endswith(".xyz"):
                xyz_files.append(os.path.join(root, f))

    # Sort for deterministic order
    xyz_files = sorted(xyz_files)

    if args.limit:
        xyz_files = xyz_files[: args.limit]

    print(f"Found {len(xyz_files)} XYZ files to process.")

    global_report = []
    requires_mace = []

    xyz_to_smiles = XYZToSMILES()

    print("\n--- PASS 1: UFF FAST-PASS ---")
    gen_uff = OIN3DGenerator(optimizer=None, ensemble_size=1)

    for i, xyz_path in enumerate(xyz_files, 1):
        basename = os.path.splitext(os.path.basename(xyz_path))[0]
        print(f"[{i}/{len(xyz_files)}] UFF Pass: {basename}...", end=" ", flush=True)

        report = initialize_report(xyz_path)

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

        success, last_xyz = _attempt_generation("UFF_1", gen_uff, oin1_string, xyz_path, report)

        if success:
            save_artifacts(report, last_xyz, output_dir, is_final=True)
            global_report.append(report)
            print("SUCCESS")
        else:
            report["status"] = "pending_mace"
            save_artifacts(report, last_xyz, output_dir, is_final=False)
            requires_mace.append((xyz_path, oin1_string, report))
            print("FAILED (queued for MACE)")

    if requires_mace:
        print(f"\n--- PASS 2: MACE HEAVY-PASS ({len(requires_mace)} files) ---")
        gen_mace_1 = OIN3DGenerator(optimizer="mace-omol-0-extra-large-1024", ensemble_size=1)
        gen_mace_5 = OIN3DGenerator(optimizer="mace-omol-0-extra-large-1024", ensemble_size=5)

        for i, (xyz_path, oin1_string, report) in enumerate(requires_mace, 1):
            basename = report["molecule"]
            print(f"[{i}/{len(requires_mace)}] MACE Pass: {basename}...", flush=True)

            # Attempt MACE_1
            print("  -> Trying MACE_1...", end=" ", flush=True)
            success, last_xyz = _attempt_generation(
                "MACE_1", gen_mace_1, oin1_string, xyz_path, report
            )
            if success:
                save_artifacts(report, last_xyz, output_dir, is_final=True)
                global_report.append(report)
                print("SUCCESS")
                continue
            print("FAILED")

            # Attempt MACE_5
            print("  -> Trying MACE_5...", end=" ", flush=True)
            success, last_xyz = _attempt_generation(
                "MACE_5", gen_mace_5, oin1_string, xyz_path, report
            )
            if success:
                print("SUCCESS")
            else:
                report["status"] = "failed"
                print("FAILED (Total Failure)")

            save_artifacts(report, last_xyz, output_dir, is_final=True)
            global_report.append(report)

    # Save global report
    global_path = os.path.join(output_dir, "summary_roundtrip.json")
    with open(global_path, "w") as f:
        json.dump(global_report, f, indent=2)

    # Print simple summary
    successes = sum(1 for r in global_report if r["status"] == "success")
    print(f"\nFinished processing {len(xyz_files)} files.")
    print(f"Successes: {successes}")
    print(f"Failures: {len(xyz_files) - successes}")
    print(f"Global report saved to {global_path}")


if __name__ == "__main__":
    main()
