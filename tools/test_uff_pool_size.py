import argparse
import contextlib
import itertools
import os
import sys
import time

# Ensure tests/integration is accessible to load built-in examples
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "tests", "integration"))
)
from verify_xyz_to_oin import get_examples

from oinsmiles.generation.metallogen_adapter import convert_oin_to_msmiles
from oinsmiles.generator3d import calculate_heavy_atom_rmsd, embed
from oinsmiles.generator3d.clean_geometry import TMCOptimizer
from oinsmiles.generator3d.om import get_om_from_modified_smiles


def load_test_oins(additional_oins=None):
    """Loads OIN strings from built-in examples and user input."""
    test_oins = []

    # 1. Load built-in examples from integration tests (excluding tmQM by default)
    for ex in get_examples(include_tmqm=False):
        oin = getattr(ex, "expected_oin_string", None) or getattr(ex, "oin_string", None)
        if oin:
            test_oins.append(oin)

    # 2. Add user-specified OINs
    if additional_oins:
        test_oins.extend(additional_oins)

    return test_oins


def run_benchmark(oin: str, max_attempts: int, cleaner: TMCOptimizer):
    """Runs generation and deduplication benchmark for a single OIN."""
    start_time = time.time()

    msmiles = convert_oin_to_msmiles(oin)
    metal_complex = get_om_from_modified_smiles(msmiles)

    successful_mols = []
    options = [0, 1, 2]
    scales = [0.8, 0.9, 1.0, 1.1, 1.2]
    combinations = list(itertools.product(scales, options))

    # Generate conformers in a loop, suppressing console spam
    with contextlib.redirect_stdout(sys.stderr):
        for i in range(max_attempts):
            scale, option = combinations[i % len(combinations)]

            try:
                positions = embed.get_embedding(
                    metal_complex, scale, option, align=True, use_random=True
                )
            except Exception:
                positions = None

            if positions is not None:
                tmp_complex = metal_complex.copy()
                tmp_complex.set_position(positions)

                if cleaner.clean_geometry(tmp_complex, scale):
                    mol = tmp_complex.get_molecule()
                    mol.generation_index = i + 1  # Tag with iteration number
                    successful_mols.append(mol)

    if not successful_mols:
        return "FAILED TO GENERATE", "N/A", time.time() - start_time

    # Sort by energy (handling None safely)
    successful_mols.sort(
        key=lambda m: (
            getattr(m, "energy", None) if getattr(m, "energy", None) is not None else float("inf")
        )
    )

    # Deduplicate based on heavy atom RMSD and Energy
    dedup_mols = []
    rmsd_threshold = 0.5
    energy_threshold = 2.0

    for mol in successful_mols:
        is_unique = True
        for acc_mol in dedup_mols:
            rmsd = calculate_heavy_atom_rmsd(mol, acc_mol)
            e1 = getattr(mol, "energy", 0.0) or 0.0
            e2 = getattr(acc_mol, "energy", 0.0) or 0.0

            if rmsd < rmsd_threshold and abs(e1 - e2) <= energy_threshold:
                is_unique = False
                break

        if is_unique:
            dedup_mols.append(mol)

    # Calculate benchmark metrics
    unique_confs = len(dedup_mols)
    max_index_survived = max((mol.generation_index for mol in dedup_mols), default=0)
    lowest_energy = getattr(dedup_mols[0], "energy", 0.0) or 0.0

    found_str = f"{unique_confs} / {max_index_survived}"
    return found_str, f"{lowest_energy:.4f}", time.time() - start_time


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark UFF pool sizes against unique conformers found."
    )
    parser.add_argument("--oins", nargs="*", default=[], help="Additional OIN strings to test")
    parser.add_argument(
        "--max-attempts", type=int, default=50, help="Number of iterations to generate structures"
    )
    args = parser.parse_args()

    test_oins = load_test_oins(args.oins)
    print(f"\nLoaded {len(test_oins)} OINs for testing.\n")

    cleaner = TMCOptimizer()
    results = []

    for i, oin in enumerate(test_oins, 1):
        print(f"[{i}/{len(test_oins)}] Benchmarking OIN: {oin[:30]}...", file=sys.stderr)
        try:
            # Truncate OIN for clean display
            display_oin = oin if len(oin) <= 47 else oin[:47] + "..."

            found_str, energy_str, duration = run_benchmark(oin, args.max_attempts, cleaner)
            results.append(
                f"{display_oin:<50} | {found_str:<25} | {energy_str:<20} | {duration:<10.2f}"
            )

        except Exception as e:
            # Print a cleaned up error snippet
            error_snip = str(e).split("\n")[0]
            error_snip = error_snip if len(error_snip) <= 17 else error_snip[:17] + "..."
            results.append(f"{display_oin:<50} | {'ERROR':<25} | {error_snip:<20} | {0.00:<10.2f}")

    # Now print the clean table at the very end, safe from C++ logs!
    print("\n\n" + "=" * 115)
    print("FINAL BENCHMARK RESULTS")
    print("=" * 115)
    print(f"{'OIN':<50} | {'Found / Max Attempt':<25} | {'Lowest Energy':<20} | {'Time (s)':<10}")
    print("-" * 115)
    for res in results:
        print(res)
    print("=" * 115 + "\n")


if __name__ == "__main__":
    main()
