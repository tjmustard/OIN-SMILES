#!/usr/bin/env python3
import argparse
import csv
import json
import os
import subprocess
import time
from typing import Dict, List, Any


def run_verify(
    output_dir: str,
    limit: int = None,
    only: str = None,
    ff_preset: str = None,
    optimizer: str = None,
    ensemble_size: int = None,
    uff_pool_size: int = None,
    cpu: bool = False,
) -> Dict[str, Any]:
    """Runs verify_roundtrip.py as a subprocess and parses its JSON output."""
    cmd = ["uv", "run", "python", "tests/integration/verify_roundtrip.py"]
    
    cmd.extend(["--output-dir", output_dir])
    if limit:
        cmd.extend(["--limit", str(limit)])
    if only:
        cmd.extend(["--only", only])
    if ff_preset:
        cmd.extend(["--ff-preset", ff_preset])
    if optimizer:
        cmd.extend(["--optimizer", optimizer])
    if ensemble_size is not None:
        cmd.extend(["--ensemble-size", str(ensemble_size)])
    if uff_pool_size is not None:
        cmd.extend(["--uff-pool-size", str(uff_pool_size)])
    if cpu:
        cmd.append("--cpu")

    print(f"Running: {' '.join(cmd)}")
    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end_time = time.time()
    
    # Parse the summary
    summary_path = os.path.join(output_dir, "summary_roundtrip.json")
    if not os.path.exists(summary_path):
        print(f"ERROR: {summary_path} not found. Subprocess failed.")
        print(result.stdout)
        print(result.stderr)
        return {
            "total": 0, "passed": 0, "failed": 0, "pass_rate": 0.0,
            "avg_rmsd": 0.0, "avg_time_sec": 0.0, "total_time": end_time - start_time,
            "error": "Summary not found"
        }

    with open(summary_path, "r") as f:
        summary = json.load(f)

    total = summary.get("total", 0)
    passed = summary.get("passed", 0)
    pass_rate = (passed / total * 100.0) if total > 0 else 0.0

    rmsds = []
    times = []
    for res in summary.get("results", []):
        metrics = res.get("metrics", {})
        if metrics.get("rmsd") is not None:
            rmsds.append(metrics["rmsd"])
        if metrics.get("time_seconds") is not None:
            times.append(metrics["time_seconds"])

    avg_rmsd = (sum(rmsds) / len(rmsds)) if rmsds else 0.0
    avg_time = (sum(times) / len(times)) if times else 0.0

    return {
        "total": total,
        "passed": passed,
        "failed": summary.get("failed", 0),
        "pass_rate": pass_rate,
        "avg_rmsd": avg_rmsd,
        "avg_time_sec": avg_time,
        "total_time": end_time - start_time,
    }


def main():
    parser = argparse.ArgumentParser(description="Run Optimization Grid Search")
    parser.add_argument("--phase", type=str, default="all", choices=["1", "2", "3", "all"],
                        help="Which phase of the grid search to run.")
    args = parser.parse_args()

    base_out_dir = "verification_artifacts_opt"
    os.makedirs(base_out_dir, exist_ok=True)

    experiments = []

    # Phase 1: UFF Grid Search (Cheapest dimensions, small subset)
    if args.phase in ["1", "all"]:
        print("\n--- [Phase 1: UFF Optimization] ---")
        for uff_pool in [5, 10, 20]:
            for ff_preset in ["loose", "default", "tight"]:
                experiments.append({
                    "phase": "1",
                    "name": f"P1_uff{uff_pool}_{ff_preset}",
                    "kwargs": {
                        "limit": 5, # Small subset
                        "optimizer": "ff",
                        "ff_preset": ff_preset,
                        "uff_pool_size": uff_pool,
                        "ensemble_size": 1,
                        "cpu": True, # FF is purely CPU
                    }
                })

    # Phase 2: MACE Grid Search (More difficult dimensions, small subset)
    if args.phase in ["2", "all"]:
        print("\n--- [Phase 2: MACE Optimization] ---")
        for ens_size in [1, 5, 10]:
            for cpu in [False, True]:
                hw = "cpu" if cpu else "gpu"
                experiments.append({
                    "phase": "2",
                    "name": f"P2_mace_ens{ens_size}_{hw}",
                    "kwargs": {
                        "limit": 5, # Small subset
                        "optimizer": "xtb",
                        "ff_preset": "default", # Assuming default is best, can adjust later
                        "uff_pool_size": 10,
                        "ensemble_size": ens_size,
                        "cpu": cpu,
                    }
                })

    # Phase 3: Full suite verification at top combinations
    if args.phase in ["3", "all"]:
        print("\n--- [Phase 3: Full Suite] ---")
        # Example top combinations to test on full suite
        combos = [
            {"optimizer": "ff", "ff_preset": "default", "uff_pool_size": 10, "ensemble_size": 1, "cpu": True},
            {"optimizer": "xtb", "ff_preset": "default", "uff_pool_size": 10, "ensemble_size": 1, "cpu": False},
        ]
        for i, combo in enumerate(combos, 1):
            hw = "cpu" if combo["cpu"] else "gpu"
            opt_short = "ff" if combo["optimizer"] == "ff" else "xtb"
            experiments.append({
                "phase": "3",
                "name": f"P3_full_{opt_short}_ens{combo['ensemble_size']}_{hw}",
                "kwargs": {
                    "limit": None, # Full suite
                    **combo
                }
            })

    results = []

    for exp in experiments:
        out_dir = os.path.join(base_out_dir, exp["name"])
        os.makedirs(out_dir, exist_ok=True)
        
        print(f"\n>> Starting Experiment: {exp['name']}")
        metrics = run_verify(output_dir=out_dir, **exp["kwargs"])
        
        row = {
            "Phase": exp["phase"],
            "Experiment": exp["name"],
            "Optimizer": exp["kwargs"].get("optimizer"),
            "FF_Preset": exp["kwargs"].get("ff_preset"),
            "UFF_Pool": exp["kwargs"].get("uff_pool_size"),
            "Ensemble": exp["kwargs"].get("ensemble_size"),
            "Hardware": "CPU" if exp["kwargs"].get("cpu") else "GPU",
            "Limit": exp["kwargs"].get("limit") or "All",
            "Total_Tests": metrics["total"],
            "Pass_Rate(%)": f"{metrics['pass_rate']:.1f}",
            "Avg_RMSD": f"{metrics['avg_rmsd']:.4f}",
            "Avg_Time_Sec": f"{metrics['avg_time_sec']:.2f}",
            "Total_Time_Sec": f"{metrics['total_time']:.1f}",
        }
        results.append(row)
        print(f"Finished {exp['name']}: {row['Pass_Rate(%)']}% passed, Avg Time: {row['Avg_Time_Sec']}s")

    # Write CSV
    csv_path = os.path.join(base_out_dir, "optimization_results.csv")
    if results:
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nWritten results to {csv_path}")

    # Write Markdown Table
    md_path = os.path.join(base_out_dir, "optimization_results.md")
    if results:
        with open(md_path, "w") as f:
            f.write("# Optimization Results\n\n")
            # Create table header
            headers = list(results[0].keys())
            f.write("| " + " | ".join(headers) + " |\n")
            f.write("|" + "|".join(["---" for _ in headers]) + "|\n")
            
            # Create table rows (sort by Phase, then Pass Rate descending)
            sorted_results = sorted(results, key=lambda x: (x["Phase"], -float(x["Pass_Rate(%)"])))
            for row in sorted_results:
                f.write("| " + " | ".join(str(row[h]) for h in headers) + " |\n")
        print(f"Written results to {md_path}")


if __name__ == "__main__":
    main()
