import argparse
import glob
import json
import os


def main():
    parser = argparse.ArgumentParser(
        description="Rebuild summary_roundtrip.json from individual reports"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Path to the results directory containing individual_reports/",
    )
    args = parser.parse_args()

    output_dir = os.path.abspath(args.output_dir)
    indiv_dir = os.path.join(output_dir, "individual_reports")
    summary_path = os.path.join(output_dir, "summary_roundtrip.json")

    if not os.path.exists(indiv_dir):
        print(f"Error: Directory {indiv_dir} does not exist.")
        return

    print(f"Scanning {indiv_dir} for individual reports...")

    global_report = []

    # Iterate through all JSON files in the individual_reports directory
    for report_file in glob.glob(os.path.join(indiv_dir, "*.json")):
        try:
            with open(report_file, "r") as f:
                report = json.load(f)
                global_report.append(report)
        except Exception as e:
            print(f"Warning: Failed to parse {report_file}: {e}")

    if not global_report:
        print("No valid reports found. Exiting.")
        return

    # Sort the global report by molecule name for consistency
    global_report.sort(key=lambda x: x.get("molecule", ""))

    print(f"Found {len(global_report)} reports. Writing to {summary_path}...")

    with open(summary_path, "w") as f:
        json.dump(global_report, f, indent=2)

    successes = sum(1 for r in global_report if r.get("status") == "success")
    failures = len(global_report) - successes

    print("\nRebuild complete.")
    print(f"Successes: {successes}")
    print(f"Failures: {failures}")


if __name__ == "__main__":
    main()
