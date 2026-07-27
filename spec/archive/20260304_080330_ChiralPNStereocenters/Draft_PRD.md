# Draft PRD: Chiral P/N Stereocenter Support

## Metadata
- **Project Name**: OIN-SMILES — Chiral Phosphorus & Nitrogen Support
- **Version**: 1.0.0
- **Status**: Draft
- **Owner**: Architect Agent (2026-03-03)
- **Parent SuperPRD**: `spec/compiled/SuperPRD.md`

---

## 1. Introduction & Goals

### 1.1 Problem Statement
OIN-SMILES currently encodes transition metal coordination geometry (the metal center's stereochemistry) but silently drops the chirality of **non-metal stereocenters** — specifically chiral phosphorus (P) and chiral nitrogen (N) atoms in TMC ligands.

**Failure modes today:**
1. **XYZ → OIN (encoding):** The SMILES fragment for a chiral phosphine ligand (e.g., PPh₂Me with R configuration) is generated without `@`/`@@` markers. CIP information is lost.
2. **OIN → XYZ (decoding):** Architector ignores `@`/`@@` markers in the SMILES fragments. Even if markers were present, the generated 3D structure would not respect them.
3. **Round-trip failure:** `OIN(XYZ₁) ≠ OIN(XYZ₂)` for chiral complexes because the second conversion cannot reproduce the original chirality without knowing it was chiral.

**Examples of affected systems:**
- Monodentate chiral phosphines (e.g., Pd-PR₃ with 3 different R groups)
- Chelating chiral diphosphines (e.g., BINAP-Rh, chiral P bonded to metal)
- Pendant chiral centers (e.g., chiral N in an aminoacid chelate backbone, not bonded to metal)

### 1.2 Solution Overview
Two independent sub-problems, each requiring distinct engineering:

**A. Encoding (XYZ → OIN):**
- Assign CIP codes (R/S) to chiral P/N atoms on the **full complex** before fragmentation.
- For **pendant (non-Zone A) chiral centers:** `isomericSmiles=True` already passes them through; verify and add regression tests.
- For **Zone A chiral centers (P/N bonded to metal):** Use a metal pseudo-atom strategy — substitute the metal with a dummy isotope before generating the fragment SMILES, so RDKit can determine CIP with full neighbor context. Best-effort: log a warning if assignment fails.

**B. Decoding (OIN → XYZ):**
- **Molassembler** (SCINE, BSD-3-Clause) replaces Architector as the **sole** 3D structure generator for all OIN→XYZ conversions — chiral and achiral alike.
- Molassembler uses rigorous graph-theoretical stereopermutators and 4D distance geometry to generate structures that respect `@`/`@@` when present and produce correct coordination geometry for achiral complexes.
- `ArchitectorAdapter` and `ArchitectorWrapper` are deprecated and removed. There is no dual-backend or fallback path.
- **Rationale:** A single backend eliminates routing complexity, legacy string ambiguity (achiral OIN strings that were actually chiral), and the maintenance burden of two separate 3D generation pathways.

### 1.3 Target Audience
Computational chemists and software engineers using OIN-SMILES for TMC databases, catalyst screening, or stereochemical analysis.

---

## 2. Confidence Mandate
**Confidence Score**: 7/10

**Open Questions (for Red Team):**
- [ ] Can Molassembler's Python API accept a full TMC SMILES (metal + dative bonds + @/@@) directly, or does it require manual stereopermutator assignment?
- [ ] Does the metal pseudo-atom strategy reliably produce correct @/@@ for Zone A P/N, or does RDKit's CIP algorithm fail for non-standard atom environments?
- [ ] Is `scine-molassembler` pip-installable on Linux without extra system deps (SCINE has a complex C++ build chain)?
- [ ] What seed API does Molassembler expose for reproducibility?
- [ ] Which OIN geometry codes does Molassembler natively support? (LIN, TPL, SPL, TET, TPY, TBP, SPY, OCT, PBP — all must be covered since Architector is fully removed.)
- [ ] What is the expected performance envelope for Molassembler on large achiral complexes (e.g., Ferrocene, Ir(ppy)₃) compared to the retired Architector baseline?

---

## 3. Scope

### 3.1 In-Scope
- CIP assignment and @/@@ marker propagation for chiral P and N atoms (all three center types)
- Pendant chiral center encoding verification and regression tests
- Zone A chiral center encoding via metal pseudo-atom strategy (best-effort)
- **Full replacement of Architector with Molassembler** for all OIN→XYZ conversions (chiral and achiral)
- Molassembler must correctly reproduce all OIN geometry codes: LIN, TPL, SPL, TET, TPY, TBP, SPY, OCT, PBP
- Deprecation and removal of `ArchitectorAdapter` and `ArchitectorWrapper`
- Verification that all existing achiral round-trip tests (Cisplatin, Ferrocene, etc.) pass through Molassembler
- New `oin-smiles` CLI entry point (single-file XYZ → OIN conversion)
- Unit tests: known chiral P/N complexes with ground-truth OIN strings
- Integration tests: round-trip OIN string stability for both chiral and achiral complexes

### 3.2 Out-of-Scope
- Axial chirality (atropisomers, e.g., BINAP's biaryl axial chirality — defer to separate PRD)
- Other stereocenters (C, Si, S, As) — this PRD targets P and N only
- Metal-centered chirality (Λ/Δ octahedral, this is handled by the OIN geometry tags)
- Charged complex support (TD-008 in baseline — separate debt item)
- Batch/directory CLI processing (single-file only in this PRD)

---

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
| :--- | :--- | :--- | :--- |
| US-001 | As a chemist, I want `XYZToSMILES.convert()` on a complex with a pendant chiral N/P to produce an OIN string with the correct `@`/`@@` marker, so that CIP information is preserved. | 1. OIN string contains `@` or `@@` on chiral atom.<br>2. CIP code (R/S) matches RDKit assignment on input XYZ.<br>3. No regression on existing non-chiral test cases. | High |
| US-002 | As a chemist, I want `XYZToSMILES.convert()` on a complex with a Zone A chiral P (bonded to metal) to attempt CIP assignment using a pseudo-atom strategy. | 1. System tries pseudo-atom substitution.<br>2. If CIP can be assigned, @/@@ appears in OIN.<br>3. If assignment fails, a warning is logged and no @/@@ is emitted (no crash). | High |
| US-003 | As a developer, I want `OIN3DGenerator.generate()` to use Molassembler for all OIN strings, so that chiral and achiral complexes share a single generation backend. | 1. All OIN strings — with or without @/@@ — are processed by Molassembler.<br>2. Achiral complexes (Cisplatin, Ferrocene, Ir(ppy)₃) produce geometrically correct structures.<br>3. Chiral complexes produce structures with correct CIP on P/N atoms. | High |
| US-004 | As a developer, I want the round-trip `OIN(XYZ₁) == OIN(XYZ₂)` to hold for chiral complexes. | 1. `OIN(XYZ₁) == OIN(Molassembler(OIN(XYZ₁)))` — string-identical including @/@@ markers. | High |
| US-005 | As a user, I want a CLI entry point so I can run `oin-smiles convert complex.xyz` from the shell. | 1. `uv run oin-smiles convert <file.xyz>` prints OIN string to stdout.<br>2. Exit code 0 on success, 1 on failure. | Medium |
| US-006 | As a developer, I want unit tests covering ground-truth chiral OIN strings for Rh-BINAP and a pendant chiral N complex. | 1. Test with known expected OIN string passes.<br>2. Test for pendant chiral center (non-Zone A) passes.<br>3. Test for Zone A chiral P best-effort passes. | High |

---

## 5. Technical Specifications (The Blueprint)

### 5.1 Architecture & Resolved Trade-offs

#### Data Flow: XYZ → OIN (Encoding)

```
XYZ File
  → [xyz2mol.get_tmc_mol()] → RDKit Mol (full complex, metal included)
  → [NEW: CIPAssigner.assign_all(mol)] → Mol with _CIPCode properties set on all atoms
  → [Existing fragment separation] → Metal fragment + Ligand fragments
  → [OINSanitizer.generate_robust_smiles() — MODIFIED]:
      For each Zone A atom:
        if atom is chiral P/N:
          → replace metal neighbor with pseudo-atom [Zz] (heavy dummy)
          → call RDKit AssignStereochemistry on modified ligand
          → extract @/@@ marker
          → remove [Zz] from SMILES, transplant @/@@ to correct position
          → on failure: log WARNING, proceed without @/@@
      For each non-Zone A chiral P/N:
        → isomericSmiles=True handles this; verify CIPCode matches
  → [OINDiscreteAligner] → canonical slot assignment
  → [OINInlineHandler] → final OIN string with @/@@ markers intact
```

#### Data Flow: OIN → XYZ (Decoding)

```
OIN String (chiral OR achiral)
  → [OINParser.parse()] → ParsedOIN (fragment SMILES may contain @/@@)
  → [NEW: MolassemblerAdapter.convert(parsed_oin)] → Molassembler molecule
      Inputs: metal symbol, ligand SMILES (with @/@@ if chiral), slot vectors from OIN
      Molassembler: assigns stereopermutators from graph + @/@@ markers
      Molassembler DG (4D): → 3D structure with correct geometry + chirality
  → [xTB optimization] → refined structure
  → [Post-optimization CIP check]: if chiral atoms present, verify CIP preserved;
      if inverted → skip xTB, return Molassembler pre-optimization structure
```

#### Resolved Trade-offs

- **Issue:** Architector ignores @/@@ — cannot be patched to handle chirality.
  - **Options:** (A) Patch Architector | (B) Post-process Architector output | (C) Molassembler for chiral cases | (D) Molassembler for all cases
  - **Resolution:** **Option D** — Molassembler replaces Architector entirely. A single backend eliminates routing logic, removes the "legacy achiral OIN string that is actually chiral" failure mode, and eliminates the maintenance burden of two parallel adapters. `ArchitectorAdapter` and `ArchitectorWrapper` are removed.

- **Issue:** Zone A chiral atoms lose their metal neighbor during fragmentation, breaking CIP assignment.
  - **Options:** (A) Skip Zone A chirality | (B) Pseudo-atom substitution | (C) CIP from full complex, transplant to fragment
  - **Resolution:** Option C as primary: assign CIP on full complex before fragmentation, store `_CIPCode`, then call `recover_chirality_tag()` AFTER `OINSanitizer` runs to apply the correct @/@@ to the sanitized fragment SMILES. Option B (pseudo-atom) is NOT used as the primary mechanism — see Red Team C-1 finding. Best-effort — log and skip if recovery fails.

- **Issue:** Molassembler has a complex C++ build chain (SCINE project).
  - **Options:** (A) Required dep | (B) Optional soft dep | (C) Docker wrapper
  - **Resolution:** Option A — add `scine-molassembler` to `pyproject.toml` as a required dependency. Since Architector is removed, there is no fallback. A clear `ImportError` with installation instructions is raised if the package is missing.

### 5.2 System Graph Blast Radius

Files modified:

| File | Change Type | Rationale |
|------|-------------|-----------|
| `src/oinsmiles/utils/oin_aligner.py` | Modified | `OINSanitizer.generate_robust_smiles()` gets Zone A CIP preservation logic; `recover_chirality_tag()` called post-sanitization |
| `src/oinsmiles/utils/perception_tmc.py` | Modified | Pre-fragmentation `CIPAssigner.assign_all()` call added; `_CIPCode` stored per atom |
| `src/oinsmiles/core/translator.py` | Modified | `XYZToSMILES.convert()` wired to `CIPAssigner`; `SMILESToXYZ.convert()` wired to `MolassemblerAdapter` |
| `src/oinsmiles/generation/engine.py` | Modified | `OIN3DGenerator` wired exclusively to `MolassemblerAdapter`; all Architector references removed |
| `src/oinsmiles/generation/oin_parser.py` | Modified | Verify @/@@ markers survive `parse_inline_string()` and `MolToSmiles()` calls |
| `src/oinsmiles/oin/inline.py` | Modified | Verify `generate_inline_string()` preserves @/@@ during {slot} injection |
| `pyproject.toml` | Modified | Add `scine-molassembler`; remove `architector` if it is a direct dep |

Files removed (deprecated):

| File | Rationale |
|------|-----------|
| `src/oinsmiles/generation/architector_adapter.py` | Replaced by `MolassemblerAdapter` |
| `src/oinsmiles/generation/wrapper.py` | Architector wrapper no longer used |

New files:

| File | Purpose |
|------|---------|
| `src/oinsmiles/core/chirality.py` | `CIPAssigner`, `has_chiral_pn()`, `recover_chirality_tag()` utilities |
| `src/oinsmiles/generation/molassembler_adapter.py` | `MolassemblerAdapter` — converts ParsedOIN to Molassembler API for all complexes |
| `src/oinsmiles/cli.py` | CLI entry point using `argparse` |
| `tests/unit/test_chirality.py` | Unit tests: CIP assignment, @/@@ preservation through inline injection |
| `tests/fixtures/chiral_complexes/` | XYZ fixtures: pendant chiral N complex, Zone A chiral phosphine |

### 5.3 Execution Checklist (MiniPRDs)

- [ ] `spec/compiled/MiniPRD_MolassemblerSpike.md` — Pre-work: install Molassembler, verify API shape and geometry coverage, confirm @/@@ survival through inline handler
- [ ] `spec/compiled/MiniPRD_ChiralEncoding.md` — CIPAssigner + OINSanitizer modifications (XYZ→OIN)
- [ ] `spec/compiled/MiniPRD_MolassemblerAdapter.md` — MolassemblerAdapter as full Architector replacement (OIN→XYZ for all cases); deprecate ArchitectorAdapter and ArchitectorWrapper
- [ ] `spec/compiled/MiniPRD_CLI.md` — CLI entry point
- [ ] `spec/compiled/MiniPRD_ChiralTests.md` — Test fixtures, unit tests, integration round-trip tests (chiral and achiral via Molassembler)

### 5.4 API Contracts / Schema

#### `CIPAssigner` (new module: `src/oinsmiles/core/chirality.py`)

```python
class CIPAssigner:
    @staticmethod
    def assign_all(mol: Chem.Mol) -> Chem.Mol:
        """
        Assigns CIP codes to all stereocenters on the full complex mol.
        Returns mol with '_CIPCode' property set on each chiral atom.
        Must be called BEFORE fragmentation.
        """

    @staticmethod
    def has_chiral_pn(mol: Chem.Mol) -> bool:
        """Returns True if any P or N atom has _CIPCode property set."""

class PseudoAtomStrategy:
    PSEUDO_ATOMIC_NUM: int = 118  # Oganesson — highest priority by CIP rules

    @staticmethod
    def substitute_metal(mol: Chem.Mol, metal_idx: int) -> Tuple[Chem.Mol, int]:
        """Replace metal with Og pseudo-atom. Returns (modified_mol, pseudo_idx)."""

    @staticmethod
    def recover_chirality_tag(
        ligand_smiles: str,
        chiral_atom_local_idx: int,
        cip_code: str  # 'R' or 'S'
    ) -> str:
        """
        Given a ligand SMILES (without metal), the local idx of the chiral atom,
        and the known CIP code, return the SMILES with correct @/@@ applied.
        Returns original smiles with warning logged on failure.
        """
```

#### `MolassemblerAdapter` (new module: `src/oinsmiles/generation/molassembler_adapter.py`)

```python
class MolassemblerAdapter:
    def __init__(self, seed: Optional[int] = 42, timeout_seconds: int = 120):
        """
        seed: fixed seed for DG if Molassembler API supports it; None = non-deterministic.
        timeout_seconds: max wall-clock time for DG generation before TimeoutError.
        """

    def convert(self, parsed_oin: ParsedOIN) -> GeneratedStructure:
        """
        Converts ParsedOIN to a 3D structure via Molassembler for ALL input OIN strings
        (chiral and achiral). Replaces ArchitectorAdapter + ArchitectorWrapper entirely.

        Steps:
          1. Build Molassembler molecule from metal + ligand SMILES + slot vectors.
          2. Assign stereopermutators: metal center from OIN geometry tag; P/N from @/@@.
          3. Run 4D distance geometry within timeout_seconds.
          4. Run xTB optimization.
          5. If any chiral P/N: verify CIP preserved post-xTB; revert to pre-xTB if inverted.
          6. Return GeneratedStructure.

        Raises:
          ImportError: if scine-molassembler is not installed.
          TimeoutError: if DG exceeds timeout_seconds.
          MolassemblerGeometryError: if the OIN geometry code is unsupported.
        """

class GeneratedStructure:
    """Thin wrapper around Molassembler output to expose a consistent interface."""
    def write_file(self, path: str) -> None: ...
    def get_atoms(self) -> List[Tuple[str, float, float, float]]: ...
    def get_xyz_string(self) -> str: ...
```

#### CLI (`src/oinsmiles/cli.py`)

```
Usage: oin-smiles <command> [options]

Commands:
  convert <xyz_file>        Convert XYZ to OIN-SMILES string (stdout)
  generate <oin_string>     Generate XYZ from OIN string (stdout or --output file)

Options:
  --charge INT              Molecular charge (default: 0)
  --output FILE             Write XYZ to file instead of stdout
  --verbose                 Enable debug logging
```

Entry point in `pyproject.toml`:
```toml
[project.scripts]
oin-smiles = "oinsmiles.cli:main"
```

### 5.5 Dependencies

**New dependency:**
- `scine-molassembler` (BSD-3-Clause) — Reiher Group, ETH Zurich

**Existing (already in pyproject.toml):**
- `rdkit` ≥ 2025.9.3 — CIP assignment, SMILES manipulation
- `xtb` ≥ 22.1 — post-Molassembler geometry refinement (optional)
- `numpy`, `scipy` — no change

---

## 6. Negative Constraints (The "Do NOTs")

- **DO NOT** modify the OIN string syntax. No new tags, delimiters, or keywords. The `@`/`@@` markers are standard SMILES — no OIN-specific chirality tag is introduced.
- **DO NOT** break existing achiral round-trips. All currently passing test cases (Cisplatin, Transplatin, Ferrocene, Ir(ppy)₃) must pass through Molassembler and produce correct structures.
- **DO NOT** retain Architector as a fallback. There is one backend: Molassembler. If Molassembler cannot handle a geometry, the error must be explicit (`MolassemblerGeometryError`), not a silent re-route.
- **DO NOT** crash on Zone A chiral center CIP failure. Log a warning and omit `@`/`@@` — graceful degradation is mandatory.
- **DO NOT** add GPL-licensed code. Only BSD/MIT/Apache dependencies are permitted (established by the baseline and confirmed by OpenSourceTMCBuilderReport.md analysis).
- **DO NOT** hardcode chiral P/N ligand lists (analogous to the `SYMMETRIC_LIGANDS` antipattern in `oin_aligner.py`). Chirality detection must be algorithmic.
- **DO NOT** use Oganesson (Z=118) or any other high-Z pseudo-atom to substitute the metal for CIP assignment. Use CIP codes stored from the full complex pre-fragmentation and recover @/@@ post-sanitization via `recover_chirality_tag()` only.

---

## 7. Risks & Mitigation

| Risk | Severity | Mitigation |
|------|----------|------------|
| `scine-molassembler` does not support all OIN geometry codes (LIN/TPL/SPL/TET/TPY/TBP/SPY/OCT/PBP) | **Critical** | Resolve in Spike MiniPRD before all other MiniPRDs. If any geometry is unsupported, scope must be adjusted or a geometry-mapping layer added. |
| `scine-molassembler` has complex C++ build deps and may not install via pip on all platforms | High | Verify install in project `.venv` during Spike. Since there is no fallback, install failure is a hard blocker. |
| Molassembler performance regression for achiral cases (Ferrocene, Ir(ppy)₃) vs. retired Architector | High | Benchmark Molassembler on all existing test fixtures during Spike. SLA: ≤ 120 seconds per `generate()` call. Expose `timeout_seconds` parameter. |
| CIP assignment from full complex not recoverable as @/@@ in fragment SMILES after `OINSanitizer` | High | `recover_chirality_tag()` must be tested on ≥ 3 known chiral P/N fragments. If systematic failure, Zone A chirality must be descoped. |
| @/@@ markers destroyed during `OINInlineHandler.generate_inline_string()` injection step | High | Verify empirically during Spike (required pre-work). Add explicit regression test in MiniPRD_ChiralTests. |
| Molassembler's DG produces incorrect stereogeometry even with @/@@ provided | Medium | Post-generation CIP verification in `MolassemblerAdapter.convert()`. If CIP wrong, raise `MolassemblerGeometryError`. |
| xTB optimization post-Molassembler inverts a chiral center or distorts achiral geometry | Medium | Post-xTB CIP check for all P/N atoms. If inverted → use pre-xTB structure. For achiral: compare RMSD with pre-xTB; if > 2Å, flag as suspect. |
| `scine-molassembler` version API breakage across releases | Medium | Pin to `scine-molassembler>=X.Y.Z,<X+1` after Spike confirms working version. |
| Atom ordering in Molassembler output differs from input XYZ, breaking OIN round-trip stability | Medium | Existing issue for Ir(ppy)₃ (~5Å RMSD). Verify canonical OIN string is identical regardless of Molassembler atom ordering. |
| CLI error handling: unformatted Python traceback on Molassembler failure | Low | All exceptions caught in `cli.py`; user-facing messages must be human-readable. Exit code 1 with stderr message only. |

---

## 8. Success Metrics

1. **Achiral regression**: All currently passing integration tests (`Cisplatin`, `Transplatin`, `Ferrocene`, `fac/mer-Ir(ppy)₃`) pass through Molassembler with OIN string stability — strings identical to pre-migration baselines.
2. **Chiral round-trip stability**: `OIN(XYZ₁) == OIN(Molassembler(OIN(XYZ₁)))` for ≥ 2 chiral test complexes (one pendant chiral N + one Zone A chiral phosphine), including `@`/`@@` markers.
3. **CIP correctness**: For pendant chiral centers, CIP code in OIN string matches RDKit's CIP assignment on the input XYZ in ≥ 95% of well-defined test cases (where RDKit confirms CIP can be assigned).
4. **Zone A best-effort coverage**: At least 1 Zone A chiral phosphine test case where CIP is correctly assigned and the OIN carries the correct `@`/`@@` — verified against a known crystal structure or published CIP assignment.
5. **Geometry coverage**: Molassembler successfully generates structures for all 9 OIN geometry codes (verified in Spike).
6. **Performance**: All `generate()` calls complete within 120 seconds for the test fixture set.
7. **CLI functional**: `oin-smiles convert tests/fixtures/cisplatin.xyz` returns the known Cisplatin OIN string. `oin-smiles generate "<oin_string>"` produces a valid XYZ.
