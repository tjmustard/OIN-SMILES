"""Recalculate OIN strings for existing round-trip results with the current code.

Re-encodes each report's input and generated XYZ files and re-judges the string
comparison with the current canonical comparator. CONSERVATIVE by design:

* A row is only flipped failed -> success when the canonical keys match AND the
  stored atom counts match AND a previously recorded RMSD metric already passed
  (< 1.0). String keys alone say nothing about geometry -- the original version
  of this tool flipped string-mismatch failures to "success" without re-running
  RMSD/atom-count, which created unverified-success rows.
* When keys match but geometry was never verified, the row is marked failed with
  a "geometry unverified" error so ``--rerun-failed`` picks it up.
* A failed recalculation keeps the old SMILES strings for forensics instead of
  overwriting them with null.
* Every touched row is stamped with recalculation provenance.
"""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../tests/integration")))

from verify_roundtrip import (
    canonical_roundtrip_key,
    normalize_oin_for_comparison,
    read_atom_count,
)

from oinsmiles import XYZToSMILES


def _find_generated_xyz(output_dir, basename):
    for candidate in (
        os.path.join(output_dir, "structures", f"{basename}_generated.xyz"),
        os.path.join(output_dir, "test_failures", basename, "last_generated.xyz"),
    ):
        if os.path.exists(candidate):
            return candidate
    return None


def _recalc_provenance():
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
    return {
        "recalculated_at": datetime.now().isoformat(timespec="seconds"),
        "recalc_commit_id": commit,
        "recalc_rdkit_version": rdBase.rdkitVersion,
    }


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

    provenance = _recalc_provenance()
    xyz_to_smiles = XYZToSMILES()

    print(f"Loaded {len(global_report)} reports from {summary_path}")
    print("Recalculating OIN SMILES using the latest code...")

    changed_status = 0

    for report in global_report:
        basename = report["molecule"]
        input_xyz = report.get("input_xyz")
        old_status = report.get("status")
        gen_xyz = _find_generated_xyz(output_dir, basename)

        report.update(provenance)

        # Re-encode the ORIGINAL input. On failure keep the stored strings so
        # the row remains inspectable; record the crash as the error.
        try:
            new_smiles_1 = xyz_to_smiles.convert(input_xyz) if input_xyz else None
            encode_error = None
        except Exception as e:
            new_smiles_1 = None
            encode_error = f"XYZToSMILES failed on recalculation: {type(e).__name__}: {e}"

        if new_smiles_1 is None:
            report["status"] = "failed"
            report["error"] = encode_error or "Input XYZ missing on recalculation"
            if old_status != report["status"]:
                changed_status += 1
                print(f"[{basename}] Status changed: {old_status} -> failed (encode)")
            _write_individual(indiv_dir, basename, report)
            continue

        report["smiles_1"] = new_smiles_1

        if gen_xyz is None:
            # Nothing generated to compare against: a success without its
            # structure is unverifiable; a failure stays a failure.
            if old_status == "success":
                report["status"] = "failed"
                report["error"] = (
                    "No generated structure found on recalculation; "
                    "geometry unverified -- rerun required"
                )
                changed_status += 1
                print(f"[{basename}] Status changed: success -> failed (no structure)")
            _write_individual(indiv_dir, basename, report)
            continue

        try:
            new_smiles_2 = xyz_to_smiles.convert(gen_xyz)
        except Exception as e:
            new_smiles_2 = None
            report["error"] = (
                f"Generated-XYZ re-encode failed on recalculation: {type(e).__name__}: {e}"
            )

        if new_smiles_2 is None:
            report["status"] = "failed"
            if old_status != "failed":
                changed_status += 1
                print(f"[{basename}] Status changed: {old_status} -> failed (re-encode)")
            _write_individual(indiv_dir, basename, report)
            continue

        report["smiles_2"] = new_smiles_2

        if canonical_roundtrip_key(new_smiles_1) != canonical_roundtrip_key(new_smiles_2):
            s1 = normalize_oin_for_comparison(new_smiles_1.strip())
            s2 = normalize_oin_for_comparison(new_smiles_2.strip())
            report["status"] = "failed"
            report["error"] = f"String mismatch on recalculation. Exp: {s1}, Got: {s2}"
            if old_status != "failed":
                changed_status += 1
                print(f"[{basename}] Status changed: {old_status} -> failed (mismatch)")
            _write_individual(indiv_dir, basename, report)
            continue

        # Keys match. Only call it a success if geometry is ALSO verified:
        # matching atom counts plus an already-recorded passing RMSD.
        atoms_in = read_atom_count(input_xyz)
        atoms_gen = read_atom_count(gen_xyz)
        rmsd = (report.get("metrics") or {}).get("rmsd")

        if atoms_in != atoms_gen:
            report["status"] = "failed"
            report["error"] = (
                f"Atom count mismatch on recalculation. Input {atoms_in} != Gen {atoms_gen}"
            )
        elif rmsd is not None and rmsd < 1.0:
            report["status"] = "success"
            report["error"] = None
        else:
            report["status"] = "failed"
            report["error"] = (
                "OIN keys match on recalculation but geometry is unverified "
                "(no passing RMSD on record) -- rerun required"
            )

        if old_status != report["status"]:
            changed_status += 1
            print(f"[{basename}] Status changed: {old_status} -> {report['status']}")
        _write_individual(indiv_dir, basename, report)

    with open(summary_path, "w") as f:
        json.dump(global_report, f, indent=2)

    successes = sum(1 for r in global_report if r.get("status") == "success")
    print(f"\nFinished updating {len(global_report)} reports.")
    print(f"Statuses changed: {changed_status}")
    print(f"New Totals -> Successes: {successes} | Failures: {len(global_report) - successes}")
    print(f"Updated global report saved to {summary_path}")


def _write_individual(indiv_dir, basename, report):
    if os.path.isdir(indiv_dir):
        with open(os.path.join(indiv_dir, f"{basename}.json"), "w") as f:
            json.dump(report, f, indent=2)


if __name__ == "__main__":
    main()
