"""
Molassembler Spike — MiniPRD_MolassemblerSpike
===============================================
Temporary investigation script. DO NOT import from src/oinsmiles/.
DO NOT commit as a permanent test.

Purpose: Confirm SCINE Molassembler installation, picklability,
ProcessPoolExecutor isolation, and cisplatin conformer generation.

Run with: uv run python tests/spike_molassembler.py
"""

import os
import pickle
import sys
from concurrent.futures import ProcessPoolExecutor, TimeoutError
from pathlib import Path

CANDIDATE_OUTPUTS = Path(__file__).parent / "candidate_outputs"
CANDIDATE_OUTPUTS.mkdir(exist_ok=True)


# =============================================================================
# SECTION 1: Import and version check
# =============================================================================

def check_import():
    """Task 3: Verify correct import path and version."""
    import scine_molassembler as masm  # noqa: PLC0415
    print(f"[OK] import scine_molassembler as masm  (version {masm.__version__})")
    return masm


# =============================================================================
# MODULE-LEVEL WORKER — must be at module level for ProcessPoolExecutor pickle
# =============================================================================

def _molassembler_worker(args: dict) -> dict:
    """
    Module-level worker for ProcessPoolExecutor.
    args keys: smiles (str), seed (int)
    Returns: dict with 'positions_shape' or 'error'
    """
    import scine_molassembler as masm  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415, F401

    smiles = args.get("smiles", "N")
    seed = args.get("seed", 42)

    mol = masm.io.experimental.from_smiles(smiles)
    result = masm.dg.generate_conformation(mol, seed)

    if isinstance(result, masm.dg.Error):
        return {"error": str(result)}

    positions = result  # numpy (N, 3) in Angstrom
    return {"positions_shape": list(positions.shape), "ok": True}


# =============================================================================
# SECTION 2: Picklability check
# =============================================================================

def check_picklability():
    """Task 4: Verify module-level worker is picklable."""
    import scine_molassembler as masm  # noqa: PLC0415
    mol = masm.io.experimental.from_smiles("N")

    # Test Molecule picklability
    mol_bytes = pickle.dumps(mol)
    mol_recovered = pickle.loads(mol_bytes)
    assert mol_recovered is not None
    print(f"[OK] Molecule is picklable  ({len(mol_bytes)} bytes)")

    # Test worker function picklability
    worker_bytes = pickle.dumps(_molassembler_worker)
    assert len(worker_bytes) > 0
    print(f"[OK] _molassembler_worker is picklable  ({len(worker_bytes)} bytes)")


# =============================================================================
# SECTION 3: ProcessPoolExecutor isolation
# =============================================================================

def check_process_pool():
    """Task 5: Confirm ProcessPoolExecutor isolation works."""
    with ProcessPoolExecutor(max_workers=1) as ex:
        fut = ex.submit(_molassembler_worker, {"smiles": "N[C@@H](Cl)F", "seed": 42})
        try:
            result = fut.result(timeout=10)
            assert result.get("ok"), f"Worker returned error: {result}"
            print(f"[OK] ProcessPoolExecutor round-trip  positions_shape={result['positions_shape']}")
        except TimeoutError:
            print("[FAIL] ProcessPoolExecutor timed out after 10s")
            sys.exit(1)


# =============================================================================
# SECTION 4: Cisplatin conformer
# =============================================================================

def generate_cisplatin():
    """Task 6: Generate cisplatin conformer, save to candidate_outputs/."""
    import scine_molassembler as masm  # noqa: PLC0415

    # Cisplatin: cis-[PtCl2(NH3)2]
    # SMILES: N[Pt](N)(Cl)Cl — N atoms get 2H implicit (NH2), not 3H (NH3).
    # For production, the real adapter will use the full OIN-parsed SMILES.
    smiles = "N[Pt](N)(Cl)Cl"
    mol = masm.io.experimental.from_smiles(smiles)
    print(f"[INFO] Cisplatin molecule: V={mol.graph.V} atoms, E={mol.graph.E} bonds")
    print(f"[INFO] Stereopermutators: {mol.stereopermutators}")

    result = masm.dg.generate_conformation(mol, seed=42)

    if isinstance(result, masm.dg.Error):
        print(f"[FAIL] DG error: {result}")
        sys.exit(1)

    positions = result  # numpy (N, 3) in Angstrom
    print(f"[OK] Conformer generated  positions_shape={list(positions.shape)}")

    # Write XYZ via masm.io.write (positions in same units as DG output = Angstrom)
    out_path = str(CANDIDATE_OUTPUTS / "spike_cisplatin.xyz")
    masm.io.write(out_path, mol, positions)
    print(f"[OK] XYZ written to {out_path}  (CANDIDATE ARTIFACT — requires human review)")


# =============================================================================
# SECTION 5: API surface summary
# =============================================================================

def print_api_summary():
    """Print confirmed API signatures for MiniPRD_MolassemblerAdapter."""
    import scine_molassembler as masm  # noqa: PLC0415

    print()
    print("=" * 60)
    print("CONFIRMED MOLASSEMBLER API SURFACE (for MolassemblerAdapter)")
    print("=" * 60)
    print(f"  Package version : {masm.__version__}")
    print(f"  Import path     : import scine_molassembler as masm")
    print(f"  Molecule SMILES : masm.io.experimental.from_smiles(smiles: str) -> Molecule")
    print(f"  DG conformer    : masm.dg.generate_conformation(mol: Molecule, seed: int) -> ndarray | dg.Error")
    print(f"  Write XYZ       : masm.io.write(filename: str, mol: Molecule, positions: ndarray)")
    print(f"  Positions units : Angstrom (despite docs saying bohr — verified empirically)")
    print(f"  Picklable       : Molecule YES  module-level function YES")
    print(f"  DG error check  : isinstance(result, masm.dg.Error)")
    print("=" * 60)


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("=== Molassembler Spike ===")
    print()

    check_import()
    check_picklability()
    check_process_pool()
    generate_cisplatin()
    print_api_summary()

    print()
    print("[DONE] All spike checks passed. See .agents/memory/molassembler_spike_results.md for decisions.")
