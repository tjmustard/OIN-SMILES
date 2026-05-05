"""Compare DG strategy performance across curated examples.

Evaluates "single", "ensemble", and "directed" strategies by running the full
round-trip (XYZ→OIN→3D→OIN) on each of the curated examples. Records pass/fail
+ geometric quality metrics (min_dist, RMSD) for side-by-side comparison.
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generation.molassembler_adapter import _min_inter_atomic_dist
from oinsmiles import XYZToSMILES
import numpy as np
from rdkit import Chem
import tempfile
import argparse
import time
import json

from verify_xyz_to_oin import get_examples
from verify_roundtrip import normalize_oin_for_comparison, read_atom_count
from rmsd_utils import calculate_tmc_rmsd


def _xyz_positions(xyz_path: str) -> np.ndarray:
    """Return (N, 3) float array of atom positions parsed from an XYZ file."""
    with open(xyz_path, "r") as f:
        lines = f.readlines()
    if len(lines) < 3:
        return np.zeros((0, 3))
    atom_count = int(lines[0].strip())
    positions = []
    for line in lines[2 : 2 + atom_count]:
        parts = line.split()
        if len(parts) >= 4:
            try:
                positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
            except ValueError:
                pass
    return np.array(positions, dtype=float) if positions else np.zeros((0, 3))


def run_strategy_test(
    example,
    generator: OIN3DGenerator,
    xyz_to_smiles,
    strategy: str = "single",
    output_dir: str = None,
    example_index: int = None,
) -> dict:
    """Run a single strategy test for one example.

    Parameters
    ----------
    example : Example
        The test example with xyz_content.
    generator : OIN3DGenerator
        Generator configured for a specific strategy.
    xyz_to_smiles : XYZToSMILES
        Converter for XYZ→OIN.
    strategy : str
        Strategy name ("single", "ensemble", or "directed") for output naming.
    output_dir : str, optional
        If provided, save generated XYZ files here with suffix.
    example_index : int, optional
        Index of example (1-based) for filename prefix (e.g., Ex1_, Ex2_).

    Returns
    -------
    dict
        Keys: ok, exception, min_dist, rmsd, oin_stable, elapsed_s, xyz_saved
    """
    result = {
        "ok": False,
        "exception": None,
        "min_dist": None,
        "rmsd": None,
        "oin_stable": False,
        "elapsed_s": None,
        "xyz_saved": None,
    }

    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            input_xyz_path = f.name
            f.write(example.xyz_content)

        try:
            oin1_string = xyz_to_smiles.convert(input_xyz_path)

            gen_xyz_path = None
            try:
                t0 = time.perf_counter()
                gen_result = generator.generate(oin1_string)
                result["elapsed_s"] = time.perf_counter() - t0

                gen_xyz_content = gen_result.xyz

                # Save to output_dir if provided
                if output_dir and example_index is not None:
                    import re
                    # Sanitize example name: remove special chars, keep alphanumeric and dash
                    safe_name = re.sub(r'[^\w\-]', '', example.name)
                    # Use Ex1_ExampleName prefix
                    prefix = f"Ex{example_index}_{safe_name}_"

                    # Save original XYZ only once (on single strategy)
                    if strategy == "single":
                        orig_path = os.path.join(output_dir, f"{prefix}original.xyz")
                        with open(orig_path, 'w') as f:
                            f.write(example.xyz_content)

                    xyz_out_path = os.path.join(output_dir, f"{prefix}{strategy}.xyz")
                    with open(xyz_out_path, 'w') as f:
                        f.write(gen_xyz_content)
                    result["xyz_saved"] = xyz_out_path

                with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as gen_f:
                    gen_xyz_path = gen_f.name
                    gen_f.write(gen_xyz_content)

                oin2_string = xyz_to_smiles.convert(gen_xyz_path)

                result["min_dist"] = _min_inter_atomic_dist(_xyz_positions(gen_xyz_path))

                mol_orig = Chem.MolFromXYZFile(input_xyz_path)
                mol_gen_xyz = Chem.MolFromXYZFile(gen_xyz_path)  # topology-free, for RMSD
                # Bonded mol for MOL/SDF output (from generator when available)
                mol_gen_bonded = gen_result.mol if gen_result.mol is not None else mol_gen_xyz

                # Write MOL and SDF files — mol_gen_bonded has bonds if generator produced them
                if output_dir and example_index is not None:
                    import re as _re_mol
                    safe_name = _re_mol.sub(r'[^\w\-]', '', example.name)
                    prefix = f"Ex{example_index}_{safe_name}_"

                    if mol_orig and strategy == "single":
                        try:
                            orig_mol_path = os.path.join(output_dir, f"{prefix}original.mol")
                            Chem.MolToMolFile(mol_orig, orig_mol_path)
                            orig_sdf_path = os.path.join(output_dir, f"{prefix}original.sdf")
                            writer = Chem.SDWriter(orig_sdf_path)
                            writer.write(mol_orig)
                            writer.close()
                        except Exception:
                            pass
                    if mol_gen_bonded:
                        try:
                            gen_mol_path = os.path.join(output_dir, f"{prefix}{strategy}.mol")
                            Chem.MolToMolFile(mol_gen_bonded, gen_mol_path)
                            gen_sdf_path = os.path.join(output_dir, f"{prefix}{strategy}.sdf")
                            writer = Chem.SDWriter(gen_sdf_path)
                            writer.write(mol_gen_bonded)
                            writer.close()
                        except Exception:
                            pass

                result["rmsd"] = calculate_tmc_rmsd(mol_orig, mol_gen_xyz, mol2_bonded=mol_gen_bonded) if (mol_orig and mol_gen_xyz) else None

                s1 = normalize_oin_for_comparison(oin1_string.strip())
                s2 = normalize_oin_for_comparison(oin2_string.strip())
                result["oin_stable"] = s1 == s2

                atom_count_match = read_atom_count(input_xyz_path) == read_atom_count(gen_xyz_path)
                rmsd_ok = result["rmsd"] is not None and result["rmsd"] < 1.0
                result["ok"] = result["oin_stable"] and rmsd_ok and atom_count_match

            finally:
                if gen_xyz_path:
                    try:
                        os.unlink(gen_xyz_path)
                    except OSError:
                        pass

        finally:
            try:
                os.unlink(input_xyz_path)
            except OSError:
                pass

    except Exception as e:
        result["exception"] = str(e)

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Compare DG strategies on curated examples"
    )
    parser.add_argument("--output-dir", type=str, help="Directory to save detailed results")
    parser.add_argument("--limit", type=int, help="Limit number of examples to run")
    parser.add_argument("--include-tmqm", action="store_true", help="Include tmQM examples (slow)")
    args = parser.parse_args()

    output_dir = args.output_dir
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    examples = get_examples(include_tmqm=args.include_tmqm)
    if args.limit:
        examples = examples[: args.limit]

    xyz_to_smiles = XYZToSMILES()
    strategies = ["single", "ensemble", "directed"]
    generators = {s: OIN3DGenerator(dg_strategy=s) for s in strategies}

    print(f"\n{'='*100}")
    print(f"DG Strategy Comparison on {len(examples)} Examples")
    print(f"{'='*100}\n")

    all_results = {}

    for i, example in enumerate(examples, 1):
        print(f"\n[{i}/{len(examples)}] {example.name}")
        print(f"  {'─' * 80}")

        example_results = {}
        for strategy in strategies:
            result = run_strategy_test(
                example,
                generators[strategy],
                xyz_to_smiles,
                strategy=strategy,
                output_dir=output_dir,
                example_index=i,
            )
            example_results[strategy] = result

            status = "✓ PASS" if result["ok"] else "✗ FAIL"
            details = []
            if result["exception"]:
                details.append(f"Exception: {result['exception']}")
            else:
                if not result["oin_stable"]:
                    details.append("OIN unstable")
                if result["rmsd"] is None:
                    details.append("RMSD: N/A")
                elif result["rmsd"] >= 1.0:
                    details.append(f"RMSD: {result['rmsd']:.2f} (high)")
                else:
                    details.append(f"RMSD: {result['rmsd']:.2f}")
                if result["min_dist"] is not None:
                    details.append(f"Min-dist: {result['min_dist']:.2f}Å")
                if result["elapsed_s"] is not None:
                    details.append(f"Time: {result['elapsed_s']:.2f}s")

            detail_str = " | ".join(details) if details else "OK"
            print(f"  {strategy:10s}: {status:10s}  {detail_str}")

        all_results[example.name] = example_results

    print(f"\n{'='*100}")
    print("Summary Table (Pass/Fail Counts)")
    print(f"{'='*100}\n")

    summary = {strategy: {"pass": 0, "fail": 0} for strategy in strategies}
    for strat_results in all_results.values():
        for strategy, result in strat_results.items():
            summary[strategy]["pass" if result["ok"] else "fail"] += 1

    print(f"{'Strategy':<15} {'Pass':<8} {'Fail':<8} {'Pass Rate':<12}")
    print("─" * 45)
    for strategy in strategies:
        p = summary[strategy]["pass"]
        fa = summary[strategy]["fail"]
        total = p + fa
        rate = f"{100.0 * p / total:.1f}%" if total > 0 else "N/A"
        print(f"{strategy:<15} {p:<8} {fa:<8} {rate:<12}")

    if output_dir:
        md_path = os.path.join(output_dir, "comparison.md")
        with open(md_path, "w") as f:
            f.write("# DG Strategy Comparison\n\n")
            f.write(f"Tested {len(examples)} examples with 3 strategies.\n\n")
            f.write("## Generated Structures\n\n")
            f.write("For each example, generated XYZ files are saved with indexed example names:\n")
            f.write("- `Ex1_ExampleName_original.xyz`, `Ex2_ExampleName_original.xyz`, ... — original input structures\n")
            f.write("- `Ex1_ExampleName_single.xyz`, `Ex2_ExampleName_single.xyz`, ... — generated by single strategy\n")
            f.write("- `Ex1_ExampleName_ensemble.xyz`, `Ex2_ExampleName_ensemble.xyz`, ... — generated by ensemble strategy\n")
            f.write("- `Ex1_ExampleName_directed.xyz`, `Ex2_ExampleName_directed.xyz`, ... — generated by directed strategy\n\n")
            f.write("These can be opened in any molecular viewer (PyMOL, Jmol, etc.) for visual QA.\n\n")

            f.write("## Per-Strategy Pass Rates\n\n")
            f.write("| Strategy | Pass | Fail | Rate |\n")
            f.write("|----------|------|------|------|\n")
            for strategy in strategies:
                p = summary[strategy]["pass"]
                fa = summary[strategy]["fail"]
                total = p + fa
                rate = f"{100.0 * p / total:.1f}%" if total > 0 else "N/A"
                f.write(f"| {strategy} | {p} | {fa} | {rate} |\n")

            f.write("\n## Per-Example Results\n\n")
            f.write("| Example | Single | Ensemble | Directed |\n")
            f.write("|---------|--------|----------|----------|\n")
            for example_name, strat_results in all_results.items():
                cols = [example_name]
                for strategy in strategies:
                    result = strat_results[strategy]
                    if result["ok"]:
                        cols.append(
                            f"✓ ({result['rmsd']:.2f}Å, {result['min_dist']:.2f}Å)"
                            if result["rmsd"] is not None
                            else "✓"
                        )
                    else:
                        reason = (result.get("exception") or "validation failed").split("\n")[0][:30]
                        cols.append(f"✗ ({reason})")
                f.write("| " + " | ".join(cols) + " |\n")

        print(f"\nMarkdown summary written to: {md_path}")

        json_path = os.path.join(output_dir, "comparison_strategies.json")
        with open(json_path, "w") as f:
            json.dump({
                "examples_tested": len(examples),
                "summary": summary,
                "results": all_results
            }, f, indent=2)
        print(f"JSON artifact written to: {json_path}")


if __name__ == "__main__":
    main()
