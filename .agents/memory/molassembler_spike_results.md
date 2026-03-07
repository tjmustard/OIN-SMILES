# Molassembler Spike Results
Date: 2026-03-04
Status: ALL CHECKS PASSED — proceed to MiniPRD_MolassemblerAdapter

---

## 1. Installation

```
uv add "scine-molassembler>=2.0.0"
→ Installed: scine-molassembler==3.0.1, scine-utilities==10.1.0
```

**Note:** Version 3.0.1 installed (not 2.x). The `>=2.0.0` pin resolves to v3.

---

## 2. Correct Import Path

```python
import scine_molassembler as masm  # underscore, NOT dot
# NOT: import scine.molassembler  (ModuleNotFoundError)
```

---

## 3. Picklability

| Target | Picklable | Bytes |
|--------|-----------|-------|
| `masm.Molecule` | ✅ YES (has `__getstate__`) | 209 |
| Module-level function | ✅ YES | 48 |
| Class method / lambda | ❌ UNTESTED — do not use |

**Conclusion:** ProcessPoolExecutor pattern is fully viable. Module-level `_molassembler_worker` confirmed picklable.

---

## 4. Confirmed API Surface

```python
import scine_molassembler as masm

# Create Molecule from SMILES
mol = masm.io.experimental.from_smiles(smiles: str) -> masm.Molecule

# Generate conformer (distance geometry)
result = masm.dg.generate_conformation(mol: masm.Molecule, seed: int) -> np.ndarray | masm.dg.Error

# Check for error
if isinstance(result, masm.dg.Error):
    error_message = str(result)  # human-readable

# Positions (success case)
positions = result  # numpy.ndarray shape (N_atoms, 3), units: ANGSTROM

# Write XYZ file
masm.io.write(filename: str, mol: masm.Molecule, positions: np.ndarray)
# Note: io.write expects positions in same units as DG output (Angstrom).
# Documentation says "bohr" but empirical check (Pt-N distance 2.09 Å ✓) confirms Angstrom.

# Molecule graph properties
mol.graph.V  # int — number of atoms
mol.graph.E  # int — number of bonds
mol.stereopermutators  # StereopermutatorList

# Molecule element type (WARNING: requires scine_utilities; fails if setuptools missing)
# mol.graph.element_type(i)  # returns Scine::Utils::ElementType (not a Python str)
# Use masm.io.write() to get element symbols — it handles this internally.
```

---

## 5. ProcessPoolExecutor Pattern (confirmed)

```python
from concurrent.futures import ProcessPoolExecutor, TimeoutError as FuturesTimeout

def _molassembler_worker(args: dict) -> dict:
    """Module-level — picklable. args: smiles, seed."""
    import scine_molassembler as masm
    mol = masm.io.experimental.from_smiles(args["smiles"])
    result = masm.dg.generate_conformation(mol, args.get("seed", 42))
    if isinstance(result, masm.dg.Error):
        return {"error": str(result)}
    return {"positions": result.tolist(), "ok": True}

with ProcessPoolExecutor(max_workers=1) as ex:
    fut = ex.submit(_molassembler_worker, {"smiles": "N[Pt](N)(Cl)Cl", "seed": 42})
    try:
        result = fut.result(timeout=60)
    except FuturesTimeout:
        raise MolassemblerTimeoutError("timed out after 60s")
```

---

## 6. Cisplatin Test

- SMILES: `N[Pt](N)(Cl)Cl` → V=9, E=8, 3 stereopermutators
- DG conformer: generated successfully (seed=42)
- XYZ saved to: `tests/candidate_outputs/spike_cisplatin.xyz` (CANDIDATE ARTIFACT)

---

## 7. Known Limitations

- `scine_utilities` requires `setuptools` at import time — `import scine_utilities` fails if setuptools not present. This is not needed for our use case (we use `masm.io.write` which handles element types internally).
- `mol.graph.element_type(i)` is not directly usable from Python (C++ type not registered). Use `masm.io.write()` for XYZ output instead of manual element-type lookups.
- The `io.write` positions "in bohr" documentation appears inaccurate. Pass Angstrom positions (raw DG output) directly.

---

## 8. Decisions for MiniPRD_MolassemblerAdapter

1. **Import**: `import scine_molassembler as masm`
2. **Worker**: Module-level `_molassembler_worker(args: dict) -> dict` — pass args as dict, return dict with `positions` list or `error` string
3. **XYZ output**: Use `masm.io.write(tmpfile, mol, positions)` with `tempfile.NamedTemporaryFile`; read back the XYZ string
4. **Timeout**: `fut.result(timeout=self.timeout)` with `concurrent.futures.TimeoutError` → `MolassemblerTimeoutError`
5. **Version note**: Actual installed version is 3.0.1 (not 2.x). Update SuperPRD NFR-1 note if needed.
