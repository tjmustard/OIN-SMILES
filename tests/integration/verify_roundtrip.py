import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

import argparse
import shutil
import tempfile
import time

from rdkit import Chem
from rmsd_utils import calculate_tmc_rmsd

from oinsmiles import XYZToSMILES
from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen as OIN3DGenerator
from oinsmiles.generation.molassembler_adapter import (
    _build_connected_smiles,
    _compute_expected_trans_sym_pairs,
    _get_binding_sym,
    _pick_masm_permutation,
)
from oinsmiles.generation.oin_parser import OINParser as _OINParser

# Comparison helpers live in the package (lightweight: rdkit + OINInlineHandler
# only) so they can be reused without importing the 3D-generation stack. They are
# re-exported here for backward compatibility with existing callers that do
# `from verify_roundtrip import normalize_oin_for_comparison`.
from oinsmiles.oin.compare import (  # noqa: F401
    canonical_roundtrip_key,
    normalize_oin_for_comparison,
    winding_canonical_key,
)


def read_atom_count(xyz_path: str) -> int:
    """Read the atom count from the first line of an XYZ file.

    Args:
        xyz_path: Path to the XYZ file.

    Returns:
        Integer atom count from line 1.

    Raises:
        ValueError: If the first line cannot be parsed as an integer.
    """
    with open(xyz_path, "r") as f:
        first_line = f.readline().strip()
    return int(first_line)


def _log_step2_inputs(oin_string: str) -> None:
    """Log diagnostic information about what will be sent to Molassembler.

    Parses the OIN string and displays the geometry, fragments, vectors,
    and computed Molassembler inputs (connected SMILES, permutation, etc).
    """
    try:
        parser = _OINParser()
        parsed = parser.parse(oin_string)

        print("  [Parsed OIN]")
        print(f"  geo_code:  {parsed.geo_code or 'NON'}")
        print(f"  fragments: {parsed.fragments}")
        print(f"             metal → fragment[{parsed.metal_fragment_idx}]")

        # Print vectors (slot assignments)
        if parsed.vectors:
            print("  vectors:")
            for vec in parsed.vectors:
                sym = (
                    _get_binding_sym(parsed.fragments[vec.fragment_idx], vec.atom_in_fragment_idx)
                    or "?"
                )
                print(
                    f"    frag[{vec.fragment_idx}] {sym:<2}  "
                    f"slot({vec.vector[0]:7.3f}, {vec.vector[1]:7.3f}, {vec.vector[2]:7.3f})"
                )

        print("\n  [Molassembler Inputs]")
        connected_smiles = _build_connected_smiles(parsed)
        print(f"  connected SMILES: {connected_smiles}")

        perm_idx = _pick_masm_permutation(parsed)
        perm_label = "TRANS" if perm_idx == 1 else "CIS/default"
        print(f"  perm_idx:         {perm_idx}  ({perm_label})")

        trans_pairs = _compute_expected_trans_sym_pairs(parsed)
        print(f"  trans sym pairs:  {trans_pairs or 'None'}")

        # Expected bindings
        expected_bindings = []
        for vec in parsed.vectors:
            if vec.fragment_idx == parsed.metal_fragment_idx:
                continue
            sym = _get_binding_sym(parsed.fragments[vec.fragment_idx], vec.atom_in_fragment_idx)
            if sym:
                expected_bindings.append((sym, tuple(vec.vector)))

        if expected_bindings:
            print(f"  expected bindings: {len(expected_bindings)} atom(s)")
            for sym, vec in expected_bindings:
                print(f"    {sym:<2}  slot({vec[0]:7.3f}, {vec[1]:7.3f}, {vec[2]:7.3f})")

    except Exception as e:
        print(f"  [Parse diagnostic error (non-fatal): {e}]")


from reporting import VerificationReporter
from verify_xyz_to_oin import get_examples


def main():
    parser = argparse.ArgumentParser(description="Verify OIN Round-Trip")
    parser.add_argument("--output-dir", type=str, help="Directory to save verification artifacts")
    parser.add_argument(
        "--limit", type=int, help="Limit number of examples to run (for fast testing)"
    )
    parser.add_argument(
        "--only",
        type=str,
        help="Run only examples whose name contains this substring (case-insensitive).",
    )
    parser.add_argument(
        "--ff-preset",
        type=str,
        default=None,
        help="FF convergence preset for the MetalloGen engine (loose/default/tight/very_tight).",
    )
    parser.add_argument(
        "--optimizer",
        type=str,
        default="xtb",
        help=(
            "Post-FF optimizer for the MetalloGen engine. Default: 'xtb' (g-xTB, "
            "subprocess wrapper with graceful FF fallback). Use "
            "'mace-omol-0-extra-large-1024' for the accurate MACE sign-off, or 'ff' for "
            "FF-only."
        ),
    )
    parser.add_argument(
        "--ensemble-size",
        type=int,
        default=1,
        help="Number of conformers to generate and optimize. Default: 1.",
    )
    parser.add_argument(
        "--cpu",
        action="store_true",
        help="Force CPU execution.",
    )
    parser.add_argument(
        "--uff-pool-size",
        type=int,
        default=None,
        help="Override the default UFF conformer pool size.",
    )
    args = parser.parse_args()

    if args.cpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    output_dir = args.output_dir
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Saving artifacts to: {output_dir}")

    print("Verifying Round-Trip Tests (V2.3 String Identity)")

    reporter = VerificationReporter("Round-Trip Verification Report")

    xyz_to_smiles = XYZToSMILES()
    ff_params = {}
    if args.uff_pool_size is not None:
        ff_params["uff_pool_size"] = args.uff_pool_size

    generator = OIN3DGenerator(
        ff_preset=args.ff_preset,
        optimizer=args.optimizer,
        ensemble_size=args.ensemble_size,
        ff_params=ff_params if ff_params else None,
    )
    if args.ff_preset or args.optimizer or args.ensemble_size or args.uff_pool_size:
        print(
            f"MetalloGen engine: ff_preset={args.ff_preset!r} optimizer={args.optimizer!r} ensemble_size={args.ensemble_size} uff_pool_size={args.uff_pool_size}"
        )

    examples = get_examples()
    if args.only:
        needles = [n.strip().lower() for n in args.only.split(",")]
        examples = [e for e in examples if any(n in e.name.lower() for n in needles)]
        print(f"Filtering to {len(examples)} example(s) matching '{args.only}'.")
    if args.limit:
        print(f"Limiting to first {args.limit} examples.")
        examples = examples[: args.limit]

    print(f"Loaded {len(examples)} examples.")

    for i, example in enumerate(examples, 1):
        print("\n==================================================")
        print(f"Running Round-Trip for Example {i}: {example.name}")
        print("==================================================")

        safe_name = "".join(x for x in example.name if x.isalnum() or x in (" ", "-", "_")).replace(
            " ", "_"
        )

        # --- Unified Round-Trip Test ---
        # Flow: XYZ -> OIN(1) -> XYZ(Gen) -> OIN(2)
        # Checks: RMSD(XYZ, XYZ_Gen) < 1.0  AND  OIN(1) == OIN(2)

        test_name = f"{example.name} (Unified Round-Trip)"
        print("\n--- Unified Round-Trip (XYZ -> OIN -> XYZ -> OIN) ---")

        try:
            # 0. Setup Paths
            base_name = f"Ex{i}_{safe_name}"

            if output_dir:
                input_xyz_path = os.path.join(output_dir, f"{base_name}_input.xyz")
                gen_xyz_path = os.path.join(output_dir, f"{base_name}_generated.xyz")
                oin1_path = os.path.join(output_dir, f"{base_name}_step1.oin")
                oin2_path = os.path.join(output_dir, f"{base_name}_step2.oin")
            else:
                tmp_dir = tempfile.mkdtemp()
                input_xyz_path = os.path.join(tmp_dir, "input.xyz")
                gen_xyz_path = os.path.join(tmp_dir, "gen.xyz")
                oin1_path = os.path.join(tmp_dir, "step1.oin")
                oin2_path = os.path.join(tmp_dir, "step2.oin")

            # -------------------------------------------------------------
            # Step 1: START (Determine Input)
            # -------------------------------------------------------------
            start_oin = None

            if example.xyz_content:
                # Flow A: Start from XYZ
                with open(input_xyz_path, "w") as f:
                    f.write(example.xyz_content)

                # XYZ -> OIN(1)
                print("Step 1: Convert Input XYZ -> OIN(1)")
                oin1_string = xyz_to_smiles.convert(input_xyz_path)
                print(f"  OIN(1): {oin1_string}")

                if output_dir:
                    with open(oin1_path, "w") as f:
                        f.write(oin1_string)

                start_oin = oin1_string
            else:
                # Flow B: Start from OIN (if no XYZ provided)
                # Use example.oin_string or expected
                start_oin = example.oin_string or example.expected_oin_string
                if not start_oin:
                    print("Skipping: No XYZ content and no OIN string provided.")
                    reporter.log_failure(test_name, "No data provided")
                    continue

                print(f"Step 1: Start from provided OIN: {start_oin}")
                oin1_string = start_oin
                # No input XYZ path to compare RMSD against later (unless we generate one?)
                # We'll skip RMSD check if started from OIN.
                input_xyz_path = None

            # -------------------------------------------------------------
            # Step 2: OIN(1) -> XYZ(Gen) (Molassembler)
            # -------------------------------------------------------------
            print("Step 2: Generate Structure OIN(1) -> XYZ(Gen)")
            _log_step2_inputs(oin1_string)
            start_time = time.time()
            gen_result = generator.generate(oin1_string)
            end_time = time.time()
            generation_time = end_time - start_time
            print(f"  Generation Time: {generation_time:.4f}s")
            with open(gen_xyz_path, "w") as f:
                f.write(gen_result.xyz)

            # -------------------------------------------------------------
            # Step 3: XYZ(Gen) -> OIN(2)
            # -------------------------------------------------------------
            print("Step 3: Convert XYZ(Gen) -> OIN(2)")
            # If we have a bonded mol from the generator, use it instead of inferring
            if gen_result.mol is not None:
                try:
                    from oinsmiles.utils.xyz2mol import get_oin_string

                    # Extract coordinates from XYZ
                    with open(gen_xyz_path, "r") as f:
                        xyz_lines = f.readlines()
                    natoms = int(xyz_lines[0].strip())
                    xyz_coords = []
                    for i in range(2, 2 + natoms):
                        parts = xyz_lines[i].split()
                        xyz_coords.append([float(x) for x in parts[1:4]])
                    import numpy as np

                    xyz_coords = np.array(xyz_coords)
                    oin2_string = get_oin_string(gen_result.mol, xyz_coords)
                except Exception as e:
                    # Fall back to xyz2mol if get_oin_string fails
                    print(f"    Note: get_oin_string failed ({e}), falling back to xyz2mol")
                    oin2_string = xyz_to_smiles.convert(gen_xyz_path)
            else:
                oin2_string = xyz_to_smiles.convert(gen_xyz_path)
            print(f"  OIN(2): {oin2_string}")

            if output_dir:
                with open(oin2_path, "w") as f:
                    f.write(oin2_string)

            # -------------------------------------------------------------
            # Step 4: Verification
            # -------------------------------------------------------------
            print("\n--- Verification Results ---")
            passed = True
            details = []

            # Check 1: String Identity (OIN 1 vs OIN 2)
            # Normalize: strip atom-ordering-dependent @SP/@OH/@TB descriptors
            # from the metal fragment before comparing — the slot assignments
            # already encode the isomer; the @XY## label is xyz-order-dependent.
            metrics: dict = {"time_seconds": round(generation_time, 4)}
            s1 = normalize_oin_for_comparison(oin1_string.strip())
            s2 = normalize_oin_for_comparison(oin2_string.strip())

            if winding_canonical_key(s1) == winding_canonical_key(s2):
                msg = "[PASS] OIN Stability: Strings Identical (normalized)"
                print(msg)
                details.append(msg)
            else:
                passed = False
                msg = "[FAIL] OIN Stability: Mismatch"
                print(msg)
                print(f"  Expected: {s1}")
                print(f"  Got:      {s2}")
                details.append(f"{msg}<br>Exp: `{s1}`<br>Got: `{s2}`")

            # Check 2: Geometric Fidelity (RMSD) - Only if we started from XYZ
            if input_xyz_path and os.path.exists(input_xyz_path):
                mol_orig = Chem.MolFromXYZFile(input_xyz_path)
                mol_gen_xyz = Chem.MolFromXYZFile(gen_xyz_path)  # topology-free, for RMSD
                # Bonded mol for MOL/SDF output (from generator when available)
                mol_gen_bonded = gen_result.mol if gen_result.mol is not None else mol_gen_xyz

                # If we have a bonded mol from the generator, use it to recompute OIN(2)
                # instead of inferring bonding from coordinates (which can fail for ansa-metallocenes)
                if gen_result.mol is not None:
                    try:
                        from oinsmiles.utils.xyz2mol import get_oin_string

                        # Extract coordinates from XYZ
                        with open(gen_xyz_path, "r") as f:
                            xyz_lines = f.readlines()
                        natoms = int(xyz_lines[0].strip())
                        xyz_coords = []
                        for i in range(2, 2 + natoms):
                            parts = xyz_lines[i].split()
                            xyz_coords.append([float(x) for x in parts[1:4]])
                        import numpy as np

                        xyz_coords = np.array(xyz_coords)
                        oin2_string = get_oin_string(gen_result.mol, xyz_coords)
                    except Exception as e:
                        # Fall back to xyz2mol if get_oin_string fails
                        print(f"    Note: get_oin_string failed ({e}), falling back to xyz2mol")
                        oin2_string = xyz_to_smiles.convert(gen_xyz_path)

                if mol_orig and mol_gen_xyz:
                    # Write MOL and SDF files — mol_gen_bonded has bonds if generator produced them
                    if output_dir:
                        try:
                            orig_mol_path = os.path.splitext(input_xyz_path)[0] + ".mol"
                            Chem.MolToMolFile(mol_orig, orig_mol_path)
                            orig_sdf_path = os.path.splitext(input_xyz_path)[0] + ".sdf"
                            writer = Chem.SDWriter(orig_sdf_path)
                            writer.write(mol_orig)
                            writer.close()
                        except Exception:
                            pass
                        if mol_gen_bonded:
                            try:
                                gen_mol_path = os.path.splitext(gen_xyz_path)[0] + ".mol"
                                Chem.MolToMolFile(mol_gen_bonded, gen_mol_path)
                                gen_sdf_path = os.path.splitext(gen_xyz_path)[0] + ".sdf"
                                writer = Chem.SDWriter(gen_sdf_path)
                                writer.write(mol_gen_bonded)
                                writer.close()
                            except Exception:
                                pass

                    rmsd = calculate_tmc_rmsd(mol_orig, mol_gen_xyz, mol2_bonded=mol_gen_bonded)
                    print(f"  RMSD Input vs Generated (coord sphere): {rmsd:.4f}")
                    metrics["rmsd"] = round(rmsd, 4)

                    if rmsd < 1.0:
                        msg = f"[PASS] Geometry: RMSD {rmsd:.4f} < 1.0"
                        print(msg)
                        details.append(f"RMSD: {rmsd:.4f}")
                    else:
                        passed = False
                        msg = f"[FAIL] Geometry: High RMSD {rmsd:.4f}"
                        print(msg)
                        details.append(f"<b>High RMSD: {rmsd:.4f}</b>")
                else:
                    msg = (
                        "[WARN] RDKit failed to support RMSD calc (atom mismatch or parsing error)"
                    )
                    print(msg)
                    details.append(msg)
                    metrics["rmsd"] = None
            else:
                details.append("(Skipped RMSD - No Input XYZ)")
                metrics["rmsd"] = None

            # Check 3: Atom Count Fidelity (detect missing H atoms)
            # Only applicable when we started from an input XYZ (Flow A).
            if input_xyz_path and os.path.exists(input_xyz_path):
                atom_count_input = read_atom_count(input_xyz_path)
                atom_count_generated = read_atom_count(gen_xyz_path)
                metrics["atom_count_input"] = atom_count_input
                metrics["atom_count_generated"] = atom_count_generated
                print(
                    f"  Atom count — Input: {atom_count_input}, Generated: {atom_count_generated}"
                )
                if atom_count_input != atom_count_generated:
                    passed = False
                    msg = (
                        f"[FAIL] Atom Count: Input {atom_count_input} "
                        f"!= Generated {atom_count_generated} (missing atoms)"
                    )
                    print(msg)
                    details.append(msg)
                else:
                    msg = f"[PASS] Atom Count: {atom_count_input}"
                    print(msg)
                    details.append(msg)
            else:
                details.append("(Skipped Atom Count - No Input XYZ)")

            if passed:
                reporter.log_success(test_name, " | ".join(details), metrics=metrics)
            else:
                reporter.log_failure(
                    test_name, "Validations Failed", got="<br>".join(details), metrics=metrics
                )

            # Cleanup
            if not output_dir:
                if "tmp_dir" in locals() and os.path.exists(tmp_dir):
                    shutil.rmtree(tmp_dir)

        except Exception as e:
            print(f"Unified Test FAILED: {e}")
            import traceback

            traceback.print_exc()
            if output_dir:
                # Persist a forensic trail for crashed examples (e.g. TiCat3/4
                # leave no step2.oin). Prefer base_name (computed at the top of
                # the try with the correct outer index) over rebuilding from i,
                # which the step-3 coordinate loops shadow.
                try:
                    stem = base_name if "base_name" in locals() else f"Ex{i}_{safe_name}"
                    err_path = os.path.join(output_dir, f"{stem}_error.txt")
                    with open(err_path, "w") as fh:
                        fh.write(f"{type(e).__name__}: {e}\n\n{traceback.format_exc()}")
                except Exception:
                    pass
            reporter.log_failure(test_name, f"Exception: {str(e)}")

    # Final Summary
    reporter.print_summary()
    if output_dir:
        json_path = os.path.join(output_dir, "summary_roundtrip.json")
        reporter.write_summary_json(json_path)
        print(f"JSON summary written to: {json_path}")


if __name__ == "__main__":
    main()
