import argparse
import json
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../tests/integration")))

from verify_roundtrip import canonical_roundtrip_key, normalize_oin_for_comparison

from oinsmiles import XYZToSMILES


def recalculate_smiles(xyz_to_smiles, xyz_path):
    if not xyz_path or not os.path.exists(xyz_path):
        return None
    try:
        return xyz_to_smiles.convert(xyz_path)
    except Exception as e:
        print(f"Failed to convert {xyz_path}: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Recalculate OIN SMILES for existing results")
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Path to the results directory containing summary_roundtrip.json",
    )
    args = parser.parse_args()

    output_dir = os.path.abspath(args.output_dir)
    summary_path = os.path.join(output_dir, "summary_roundtrip.json")
    indiv_dir = os.path.join(output_dir, "individual_reports")

    if not os.path.exists(summary_path):
        print(f"Error: {summary_path} not found.")
        return

    with open(summary_path, "r") as f:
        global_report = json.load(f)

    xyz_to_smiles = XYZToSMILES()

    print(f"Loaded {len(global_report)} reports from {summary_path}")
    print("Recalculating OIN SMILES using the latest code...")

    updated_count = 0
    changed_status = 0

    for report in global_report:
        basename = report["molecule"]
        input_xyz = report.get("input_xyz")

        # Locate generated xyz
        gen_xyz = os.path.join(output_dir, "structures", f"{basename}_generated.xyz")
        if not os.path.exists(gen_xyz):
            gen_xyz = os.path.join(output_dir, "test_failures", basename, "last_generated.xyz")
            if not os.path.exists(gen_xyz):
                gen_xyz = None

        new_smiles_1 = recalculate_smiles(xyz_to_smiles, input_xyz)
        new_smiles_2 = recalculate_smiles(xyz_to_smiles, gen_xyz) if gen_xyz else None

        report["smiles_1"] = new_smiles_1
        report["smiles_2"] = new_smiles_2

        old_status = report.get("status")

        if not new_smiles_1:
            report["status"] = "failed"
            report["error"] = "1D conversion failed on recalculation"
        elif gen_xyz and new_smiles_2:
            key1 = canonical_roundtrip_key(new_smiles_1)
            key2 = canonical_roundtrip_key(new_smiles_2)

            if key1 != key2:
                s1 = normalize_oin_for_comparison(new_smiles_1.strip())
                s2 = normalize_oin_for_comparison(new_smiles_2.strip())
                report["status"] = "failed"
                report["error"] = f"String mismatch on recalculation. Exp: {s1}, Got: {s2}"
            else:
                # If it failed previously only due to string mismatch, it might now succeed.
                # However, if it failed due to RMSD or Atom Count, we should preserve that.
                if old_status == "failed" and "String mismatch" in (report.get("error") or ""):
                    report["status"] = "success"
                    report["error"] = None
                # If it was a success, it remains a success
                elif old_status == "success":
                    report["status"] = "success"
                    report["error"] = None
        else:
            # No generated structure to compare
            if old_status != "pending_xtb":
                report["status"] = "failed"
                report["error"] = report.get(
                    "error", "No generated structure found during recalculation"
                )

        if old_status != report["status"]:
            changed_status += 1
            print(f"[{basename}] Status changed: {old_status} -> {report['status']}")

        # Write individual report back
        indiv_path = os.path.join(indiv_dir, f"{basename}.json")
        if os.path.exists(indiv_dir):
            with open(indiv_path, "w") as f:
                json.dump(report, f, indent=2)

        updated_count += 1

    # Write summary back
    with open(summary_path, "w") as f:
        json.dump(global_report, f, indent=2)

    successes = sum(1 for r in global_report if r.get("status") == "success")
    failures = len(global_report) - successes

    print(f"\nFinished updating {updated_count} reports.")
    print(f"Statuses changed: {changed_status}")
    print(f"New Totals -> Successes: {successes} | Failures: {failures}")
    print(f"Updated global report saved to {summary_path}")


if __name__ == "__main__":
    main()
