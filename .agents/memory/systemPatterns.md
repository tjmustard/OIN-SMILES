---
name: System Architecture & Design Patterns
description: OIN-SMILES technical architecture, pipeline design, and recurring patterns used across the codebase
type: reference
---

# OIN-SMILES System Architecture

## Overview
OIN-SMILES is a dual-pipeline system for lossless conversion between 3D molecular structures (XYZ) and 1D SMILES representations for Transition Metal Complexes (TMCs), using Open Isomer Notation (OIN v3.6) as the intermediate canonical format.

## Pipeline 1: XYZ → OIN (3D Structure to SMILES)

**Flow:**
```
XYZToSMILES.convert(xyz_file)
  ↓
xyz2mol.get_tmc_mol() — Graph generation using Jensen Group algorithm
  ↓
CIPAssigner.assign_all() — Stereocenters (P/N atoms) before fragmentation
  ↓
OINDiscreteAligner — Geometry detection, slot assignment, fragmentation
  ↓
OINSanitizer.generate_robust_smiles() — RDKit sanitization, properties preserved
  ↓
ChiralityRecoveryUtility.recover() — Re-apply @/@@ after sanitization
  ↓
OINInlineHandler.generate_inline_string() — Final V3.6 format
  ↓
[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}
```

**Key Architectural Decisions:**
1. **Early chirality assignment** — CIPAssigner runs on the full (pre-fragmentation) mol to capture 3D stereochemistry context
2. **Property preservation** — Atom properties (CIP codes, chiral tags) survive RDKit sanitization through custom sanitizer
3. **Regex-only parsing** — `parse_inline_string()` uses regex only (no MolFromSmiles round-trip) to preserve @/@@ markers
4. **Graceful fallback** — PseudoAtomStrategy provides wildcards (*) for uncomputable stereocenters instead of crashing

## Pipeline 2: OIN → XYZ (SMILES to 3D Structure)

**Flow:**
```
OIN3DGenerator.generate(oin_string, timeout=60)
  ↓
generation/OINParser — Inline string parsing → ParsedOIN object
  ↓
molassembler_adapter._template_generate() — Template-based placement for known geometries
  ↓
molassembler_adapter._stitch_fragment() — Attach ligands (aromatic η-ligands use ETKDG)
  ↓
molassembler_adapter.generate_conformation() — Distance geometry (DG) fallback
  ↓
RDKit mol assembly — CombineMols + dative bonds + conformer
  ↓
GeneratedStructure(xyz: str, mol: Optional[Chem.Mol])
```

**Key Architectural Decisions:**
1. **Template-first strategy** — Known geometries (SP4, SP3, OC-6, etc.) use hand-coded templates for speed/accuracy
2. **DG fallback** — Unplaced ligands or missing templates use SCINE Molassembler distance geometry (GIL-safe via ProcessPoolExecutor)
3. **ETKDG for aromatic η-ligands** — 5-membered aromatic rings (Cp, indenyl) use de-aromatization + ETKDG embedding (avoids RDKit kekulization failures)
4. **Dual output** — Returns both XYZ string (for geometry) and RDKit mol (for bond topology); mol is None for fallback cases

## Molassembler Integration

**Module:** `src/oinsmiles/generation/molassembler_adapter.py`

**API Surface:**
```python
import scine_molassembler as masm

# Create mol from SMILES
mol = masm.io.experimental.from_smiles("N[Pt](N)(Cl)Cl")

# Generate conformer (distance geometry)
result = masm.dg.generate_conformation(mol, seed=42)
# → numpy.ndarray (N_atoms × 3), units: Angstrom
# or masm.dg.Error if failed

# Write XYZ file
masm.io.write("complex.xyz", mol, positions)  # positions in Angstrom
```

**Timeout Pattern:**
```python
from concurrent.futures import ProcessPoolExecutor, TimeoutError as FuturesTimeout

def _molassembler_worker(args: dict) -> dict:
    import scine_molassembler as masm
    # ... work with masm ...
    return {"positions": result.tolist()}

with ProcessPoolExecutor(max_workers=1) as executor:
    future = executor.submit(_molassembler_worker, args)
    try:
        result = future.result(timeout=60)  # GIL-safe timeout
    except FuturesTimeout:
        raise MolassemblerTimeoutError(f"timed out after 60s")
```

**Key Property:**
- `Molecule` objects are picklable ✅
- Module-level worker functions are picklable ✅
- Ideal for `ProcessPoolExecutor` + timeout enforcement

## OIN Format (v3.6)

**Canonical inline format:**
```
[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}
```

**Components:**
- `[Pt@SP1_SPL]` — Metal atom with chirality tag (@) and geometry code (SP1 = square planar)
- `.[Cl]{0}` — Fragment (Cl ligand) with slot marker {0}
- `N{2}` — Atom-level slot marker (nitrogen at slot 2)
- Geometry codes: `SP1`, `SP3`, `OC-6`, `TBPY-5`, `PBP-1`, etc.
- Chirality tags: `@` (R), `@@` (S), or absent (achiral)
- Winding direction: `>` (clockwise), `<` (counterclockwise) — encoded in slot markers for v3.6

## Design Patterns

### Pattern 1: Property Preservation Through Pipelines
**Problem:** RDKit sanitization destroys custom atom properties (CIP codes, chiral tags)
**Solution:** Custom `OINSanitizer` that:
1. Copies atom properties before sanitization
2. Re-applies properties post-sanitization
3. Calls `CIPAssigner` AFTER `Chem.SanitizeMol()` (hard precondition)

### Pattern 2: Graceful Degradation for Non-Standard Stereocenters
**Problem:** Some P/N atoms have non-standard valence (can't assign CIP codes)
**Solution:** `PseudoAtomStrategy`:
1. Attempt CIP assignment on all P/N atoms
2. If no `_CIPCode` found: replace with wildcard atom (`*`, atomic_num=0)
3. Call `strip_pseudo_atoms()` before OIN serialization

### Pattern 3: Template-First with Fallback
**Problem:** 3D generation can fail for complex or unknown geometries
**Solution:** `molassembler_adapter`:
1. Try template-based placement (hand-coded coordinates)
2. If template missing: invoke DG (SCINE Molassembler)
3. If DG times out: return best-effort XYZ block with partial 3D

### Pattern 4: Aromatic Ligand De-aromatization
**Problem:** RDKit can't kekulize 5-membered aromatic rings (5π violates Hückel)
**Solution:** `_stitch_eta_fragment()`:
1. Extract aromatic ring SMILES
2. De-aromatize: change aromatic bonds→SINGLE, clear aromatic flags
3. Run ETKDG embedding on full fragment (not just ring)
4. Extract ring coordinates from conformer

## Code Organization

```
src/oinsmiles/
├── __init__.py              — Public API (XYZToSMILES, SMILESToXYZ)
├── cli.py                   — CLI entry point (oin-smiles command)
├── core/
│   ├── chirality.py         — CIPAssigner, ChiralityRecoveryUtility, PseudoAtomStrategy
│   ├── graph.py             — Graph utility functions
│   └── translator.py        — XYZToSMILES, SMILESToXYZ wrappers
├── generation/
│   ├── engine.py            — OIN3DGenerator (main API, timeout handling)
│   ├── molassembler_adapter.py  — Template placement, DG fallback, ETKDG logic
│   └── oin_parser.py        — Inline string parsing
├── oin/
│   ├── inline.py            — OIN v3.6 inline format handler
│   ├── parser.py            — Alternative parser (returns tuple, not used in v0.2.0+)
│   └── writer.py            — OIN generation from geometry data
└── utils/
    ├── xyz2mol.py           — Graph generation (Jensen Group algorithm)
    ├── xyz2mol_local.py     — Local fork with modifications
    └── oin_aligner.py       — OINDiscreteAligner, OINSanitizer
```

## Testing Strategy

**Unit tests:** `tests/unit/` — Individual component functionality
**Integration tests:** `tests/integration/verify_roundtrip.py` — End-to-end XYZ→OIN→XYZ with RMSD validation
**Regression tests:** Known complexes (Pt, Fe, Ir, Ti) validated against baseline RMSD thresholds

**Test Suite:** `uv run python -m unittest discover tests`
