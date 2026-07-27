
# System Patterns

## Purpose
Documents the "How" — architectural decisions, design patterns, tech stack, and conventions for OIN-SMILES.

## Architecture: Two Independent Pipelines

### XYZ → OIN (Forward)
```
XYZToSMILES.convert()
  → perception_tmc.get_tmc_mol()  # Jensen Group graph algorithm
  → CIPAssigner.assign_all()      # 3D-derived CIP codes for P/N atoms
  → OINDiscreteAligner            # slot assignment & geometry template
  → OINSanitizer                  # SMILES canonicalization
  → OIN v3.6 inline string
```

### OIN → XYZ (Reverse)
```
OIN3DGenerator.generate()
  → generation/OINParser          # returns ParsedOIN dataclass
  → MetalloGenAdapter
      → dummy-metal + RDKit CoordMap embed
      → constrained MMFF/UFF cleanup + optimizer refinement (g-xTB / MACE / FF)
  → GeneratedStructure(xyz, mol)
```

## Metal Center Invariant
The metal center is **always `fragments[0]`** in both pipelines. This is a load-bearing canonical-form property. Never reorder fragment lists.

## Two OINParser Classes (Not Interchangeable)
| File | Class | Returns | Used by |
|---|---|---|---|
| `src/oinsmiles/oin/parser.py` | `OINParser` | `(str, Dict)` | `SMILESToXYZ` (incomplete) |
| `src/oinsmiles/generation/oin_parser.py` | `OINParser` | `ParsedOIN` | `OIN3DGenerator` |

## Tech Stack
- **Language**: Python ≥ 3.10
- **Package manager**: `uv` (`uv sync`, `uv run`)
- **3D generation**: vendored MetalloGen engine (`oinsmiles.generator3d`); dummy-metal + RDKit `CoordMap` embed
- **Graph operations**: RDKit (`Chem`, `AllChem`)
- **Graph construction**: perception_core (Jensen Group xyz2mol algorithm, vendored)
- **Build**: `uv build`, entry point `oin-smiles` registered in `pyproject.toml`

## Key API (v0.2.0+)
- `XYZToSMILES().convert(path)` → `str` (OIN v3.6 inline)
- `OIN3DGenerator(timeout=60).generate(oin_str)` → `GeneratedStructure(xyz: str, mol: Optional[Chem.Mol])`
- `GeneratedStructure.xyz` always available; `.mol` is None for eta fallback cases

## Design Patterns
- **Dummy-metal CoordMap embed**: MetalloGen replaces the metal with a dummy centre, embeds ligands against coordination-geometry-matched slot vectors, then relaxes with a constrained force field + optimizer
- **CIPAssigner before fragmentation**: 3D-derived CIP codes must be computed on the full TMC mol, then propagated via atom properties
- **`AssignAtomChiralTagsFromStructure` precedes `AssignStereochemistry`** to get @/@@ tags in SMILES

## SMILES Handling Conventions
- `Chem.SanitizeMol()` is a hard precondition before `CIPAssigner.assign_all()`
- Zone A P/N atoms (direct metal binders, `total_degree < 4` in fragment) → chiral tag cleared
- `parse_inline_string()` uses regex only — no RDKit round-trip (preserves @/@@)
- `OINSanitizer.generate_robust_smiles()` returns `(smiles, kmol)` with properties preserved

## OIN Format Versions
| Version | Format |
|---|---|
| V2.4 (sidecar, obsolete) | `[Pt].[Cl] \|g:SPL\|w:1.0:0;2.0:1\|` |
| V3.0+ (inline) | `[Pt_SPL].N{0}.[Cl]{1}` |
| V3.4 | Heading atom `{0>}` |
| V3.6 (canonical) | Winding direction: `>` CW, `<` CCW |

## Technical Debt (Known)
- **TD-001**: `XYZToSMILES.convert()` defined twice — second shadows first
- **TD-002**: `OINInlineHandler.generate_inline_string()` has `pass` stub
- **TD-003**: `SMILESToXYZ` in translator.py is incomplete (dummy atoms)
- **TD-005**: `TEMPLATES`/`TEMPLATE_SPECS` duplicated across two files

## Anti-Patterns
- Do not reorder fragment lists — metal-first invariant is load-bearing
- Do not call `Chem.SanitizeMol()` after `AssignAtomChiralTagsFromStructure` — it clears stereo tags

## Conventions
- Testing: `uv run python -m unittest discover tests`
- Integration verification: `tests/integration/verify_roundtrip.py`
- RMSD metric: use **mean** RMSD across atoms, not max-per-atom (conformational flexibility makes per-atom thresholds noise-prone)
