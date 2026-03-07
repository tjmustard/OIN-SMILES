# SuperPRD — OIN-SMILES: Chiral P/N Stereocenter Support

## Metadata
- **Project Name**: OIN-SMILES
- **Version**: 1.1.0 (Chiral P/N Stereocenter Feature & Baseline Update)
- **Status**: Approved — System of Record
- **Owner**: Resolution Agent (2026-03-04) / Baseline Agent

---

## 1. Current State Baseline

Before detailing the new feature additions, it is important to establish the ground truth of the system's current functionality (Baseline).

### 1.1 Current Core Value Loop
OIN-SMILES is a toolkit designed for the representation, translation, and 3D structure generation of Transition Metal Complexes (TMCs). The core value loop consists of bidirectional translation and 3D coordinate generation:
1.  **XYZ to OIN-SMILES:** Parses an XYZ coordinate file (using `xyz2mol` heuristics to detect connectivity and haptics), canonically aligns the complex based on Principal Axes of Inertia (PAI), fragments the ligands, canonicalizes the ligands into robust SMILES (with strict explicit hydrogen tracking), determines optimal geometry mapping against predefined standard templates, and emits a canonical Open Isomer Notation (OIN) string describing the 3D topology.
2.  **OIN-SMILES to 3D Generation:** Parses an OIN string, extracts the metal and ligand fragments along with their topological binding slots (haptics and directionality), scales the geometry template vectors based on covalent radii, and constructs an `inputDict` payload. This payload is handed off to the backend engine to generate a realistic 3D structure.

### 1.2 Existing Data Models
*   **TMCGraph / RDKit Mol:** Transition Metal Complexes are primarily handled as RDKit `Mol` and `RWMol` instances. Custom properties (e.g., `__origIdx`, explicit `AtomMapNum` for tagging slots) are heavily utilized to track atoms through canonicalization.
*   **ParsedOIN (Dataclass):** The parser distills an OIN string into a structure containing `smiles`, `fragments` (with metal at index 0), `geometry` (e.g., `"SPL"`, `"OCT"`), and `vectors` mapping specific ligand atoms to template slots (`slot_idx`), capturing haptic heading and winding direction.
*   **Geometry Templates:** A standard dictionary (`TEMPLATES`) defining ideal polyhedra (e.g., Octahedral `OCT`, Square Planar `SPL`). Each slot contains a `pos` vector and a `ref` vector for relative winding determination.

### 1.3 Established Patterns
*   **Zone A Sanitization:** To prevent RDKit from arbitrarily interpreting implicit hydrogens upon re-parsing, explicit hydrogen counts are heavily enforced for all atoms coordinating to the metal.
*   **Mass-First Canonical Sorting:** Ligands are canonically sorted first by mass, then by binding atom mass, then by SMILES string.
*   **Alignment via PAI:** Molecules are universally translated to the origin and rotated such that their Principal Axes of Inertia align with Cartesian axes prior to topological mapping.
*   **Haptic Expansion & Resolution:** Haptic bonds are abstracted to a centroid "Virtual Atom" for geometry detection and expanded into a cone of vectors during generation.

### 1.4 Identified Technical Debt (Baseline)
*   **RDKit Map Number Overloading:** `AtomMapNum` is heavily overloaded to pass data between steps (e.g., forcing brackets, or sneaking slot/winding data through `MolToSmiles`).
*   **SMILES Re-parsing Redundancy:** Modules frequently generate SMILES only to immediately re-parse them with `SANITIZE=False` to apply heuristics.

---

## 2. Introduction & Goals

### 2.1 Problem Statement
The OIN-SMILES v3.6 baseline correctly encodes coordination geometry (cis/trans, facial/meridional, haptic ring orientation) but ignores stereocenters on ligand atoms. Specifically, P and N atoms in TMC ligands can bear `@`/`@@` chirality markers in SMILES, but:
1. The current XYZ→OIN pipeline destroys these markers during `OINInlineHandler.parse_inline_string()` via an RDKit canonicalization round-trip.
2. The pipeline performs no CIP assignment on Zone A atoms (those directly bonded to metal), because fragmentation removes the metal context required for RDKit CIP resolution.
3. The 3D generation backend is currently not stereo-aware for P/N centers.

Chiral phosphines (e.g., BINAP derivatives, chiral monodentate phosphines) and chiral amines are common in asymmetric catalysis TMC datasets. Without @/@@ encoding, enantiomers are indistinguishable in OIN strings.

### 2.2 Solution Overview
OIN-SMILES v1.1.0 adds:
- **Pre-fragmentation CIP assignment** (`CIPAssigner`) on the full sanitized mol, before the metal fragment is separated.
- **Chirality recovery** (`ChiralityRecoveryUtility`) that re-applies `@`/`@@` to ligand SMILES after `OINSanitizer` runs, using stored `_CIPCode` atom properties.
- **PseudoAtomStrategy** (fallback only): wildcard `*` sentinel (`PSEUDO_ATOMIC_NUM=0`) for non-standard-valence P/N atoms where CIP pre-assignment fails. Stripped from final OIN output.
- **`OINInlineHandler.parse_inline_string()` fix**: regex-only slot marker stripping — no RDKit round-trip, preserving `@`/`@@`.
- **Molassembler backend**: replaces Architector entirely. Uses SCINE graph-theoretical stereopermutators + 4D distance geometry. `ProcessPoolExecutor` isolation for GIL-safe timeout.
- **CLI entry point**: `oin-smiles xyz2oin <path>` and `oin-smiles oin2xyz <oin>`.

### 2.3 Target Audience
Computational chemists working with asymmetric TMC datasets who need lossless P/N stereocenter encoding in OIN strings for database storage and similarity search.

---

## 3. Confidence Mandate

**Confidence Score**: 9/10
**Resolved Red Team Flags**: 3 Critical, 4 High, 5 Medium, 4 Low
**Remaining uncertainty**: Molassembler Python binding picklability (confirmed by MolassemblerSpike before MolassemblerAdapter implementation begins).

---

## 4. Scope

### 4.1 In-Scope
- Pre-fragmentation CIP assignment for P and N stereocenters (`CIPAssigner`)
- Chirality recovery after `OINSanitizer` (`ChiralityRecoveryUtility`)
- `PseudoAtomStrategy` as fallback (wildcard `*` sentinel, `PSEUDO_ATOMIC_NUM=0`)
- `OINInlineHandler.parse_inline_string()` regex-only fix (no re-canonicalization)
- Molassembler integration as sole 3D generation backend (`MolassemblerAdapter`)
- `ProcessPoolExecutor`-based timeout isolation (default 60s)
- CLI entry point (`cli.py`)
- P-chiral and N-chiral test fixtures with RDKit CIP oracle validation
- BINAP as structural stability fixture
- Axial-chiral ligand SMILES must encode chirality for those specific atoms in the output SMILES

### 4.2 Out-of-Scope
- Automatic axial chirality detection
- Charged complex support (TD-008 remains open)
- Multi-metal complex support
- Architector backend (fully removed)
- Chiral metal centers (only P/N ligand stereocenters)

---

## 5. Technical Specifications

### 5.1 Architecture & Resolved Trade-offs

**Data flow (XYZ → OIN, with chirality):**
```
XYZ File
  → [get_tmc_mol()] → Full unsanitized mol
  → [Chem.SanitizeMol(mol)]                 ← HARD PRECONDITION for CIP
  → [CIPAssigner.assign_all(full_mol)]       ← BEFORE fragmentation
      → Chem.AssignStereochemistry(cleanIt=True, force=True)
      → stores '_CIPCode' ('R'/'S') on P/N atoms
  → [Fragment Separation] → Metal fragment + Ligand fragments
  → [OINSanitizer.generate_robust_smiles()]  ← Zone A explicit H brackets
  → [ChiralityRecoveryUtility.recover()]     ← reads _CIPCode, re-applies @/@@
      → if no _CIPCode: PseudoAtomStrategy(PSEUDO_ATOMIC_NUM=0)
  → [OINDiscreteAligner]                     ← geometry detection, slot mapping
  → [OINInlineHandler.generate_inline_string] → inline OIN string
```

**`OINInlineHandler.parse_inline_string()` fix (C-2):**
Strip slot markers and metal prefix via `SLOT_REGEX` and `METAL_REGEX` substitution only. No `MolFromSmiles → MolToSmiles` round-trip. Output is non-canonical SMILES but `@`/`@@`-safe.

**Resolved Trade-offs Log:**
- **C-1 (Zone A encoding contradiction):** Option C (pre-fragmentation CIP) is primary. `PseudoAtomStrategy` kept as fallback. `PSEUDO_ATOMIC_NUM = 0` (RDKit wildcard `*`). All `[Zz]`/Oganesson references purged.
- **C-2 (@/@@ corruption in parse_inline_string):** Option B — regex-only stripping. No RDKit round-trip in `parse_inline_string()`. Output SMILES is non-canonical but chirality-safe.
- **C-3 (Mol sanitization prerequisite):** `Chem.SanitizeMol(mol)` is a hard precondition called before `CIPAssigner.assign_all()`.
- **H-1 (BINAP fixture):** BINAP kept as structural stability fixture. Dedicated P-chiral fixture added with explicit `@`/`@@` literal assertion verified by RDKit CIP oracle.
- **H-2 (Timeout GIL):** `ProcessPoolExecutor(max_workers=1)` for all Molassembler calls. Default `timeout = 60` seconds. Raises `MolassemblerTimeoutError`.

### 5.2 System Graph Blast Radius

15 of 36 nodes affected (42% of graph). Reference: `spec/compiled/architecture.yml`.

**Modified (8):**
- `atom_XYZToSMILES` — calls `CIPAssigner.assign_all()`; calls `ChiralityRecoveryUtility.recover()`
- `atom_OINSanitizer` — must not destroy `_CIPCode` atom properties
- `atom_OINDiscreteAligner` — receives chiral mol with `@`/`@@` intact
- `atom_OIN3DGenerator` — rewired to call `MolassemblerAdapter`
- `atom_OINInlineHandler` — `parse_inline_string()` rewritten to regex-only stripping
- `atom_OINParser_oin` — must not destroy `@`/`@@`
- `atom_SMILESToXYZ` — rewired to `OIN3DGenerator`
- `mod_generation` — Architector dependency removed; SCINE/Molassembler added

**Removed (2):**
- `atom_ArchitectorAdapter`
- `atom_ArchitectorWrapper`

**Added (5):**
- `atom_CIPAssigner` → `src/oinsmiles/core/chirality.py`
- `atom_ChiralityRecoveryUtility` → `src/oinsmiles/core/chirality.py`
- `atom_PseudoAtomStrategy` → `src/oinsmiles/core/chirality.py`
- `atom_MolassemblerAdapter` → `src/oinsmiles/generation/molassembler_adapter.py`
- `atom_CLI` → `src/oinsmiles/cli.py`

### 5.3 Execution Checklist (MiniPRDs)
- [x] `spec/compiled/MiniPRD_MolassemblerSpike.md` — validate SCINE import + `ProcessPoolExecutor` picklability
- [x] `spec/compiled/MiniPRD_ChiralEncoding.md` — `CIPAssigner`, `ChiralityRecoveryUtility`, `PseudoAtomStrategy`, `parse_inline_string()` fix
- [x] `spec/compiled/MiniPRD_MolassemblerAdapter.md` — `MolassemblerAdapter` + `OIN3DGenerator` rewire + delete Architector
- [ ] `spec/compiled/MiniPRD_CLI.md` — `cli.py` entry point
- [ ] `spec/compiled/MiniPRD_ChiralTests.md` — test fixtures + RDKit CIP oracle

### 5.4 API Contracts

```python
# src/oinsmiles/core/chirality.py
class CIPAssigner:
    def assign_all(self, mol: Chem.Mol) -> Chem.Mol:
        # Precondition: mol must be sanitized (Chem.SanitizeMol called by caller).
        # Stores '_CIPCode' ('R' or 'S') as an atom integer property.
        pass

class ChiralityRecoveryUtility:
    def recover(self, mol: Chem.Mol) -> Chem.Mol:
        # Restores @/@@ chirality tags on P/N atoms using '_CIPCode'.
        pass

class PseudoAtomStrategy:
    PSEUDO_ATOMIC_NUM: int = 0

# src/oinsmiles/generation/molassembler_adapter.py
class MolassemblerTimeoutError(RuntimeError): pass

class MolassemblerAdapter:
    def generate(self, parsed_oin: ParsedOIN) -> str:
        # Runs Molassembler via ProcessPoolExecutor
        pass

def _molassembler_worker(args: dict) -> str:
    # Picklable worker submitted to ProcessPoolExecutor
    pass

# src/oinsmiles/generation/engine.py
class OIN3DGenerator:
    def generate(self, oin_string: str) -> str:
        # Uses MolassemblerAdapter
        pass

# src/oinsmiles/cli.py
def main() -> None:
    # oin-smiles xyz2oin <path> | oin-smiles oin2xyz <oin>
    pass
```

### 5.5 Dependencies
**Added:** `scine-molassembler >= 2.0.0`, `concurrent.futures`
**Removed:** `architector`, `xTB`
**Retained:** `RDKit >= 2025.9.3`, `NumPy`, `SciPy`, `mendeleev`, `networkx`, `openbabel-wheel`

---

## 6. Negative Constraints
- **DO NOT** use Oganesson (Z=118) or `[Zz]` as a pseudo-atom sentinel.
- **DO NOT** use `threading.Thread` for Molassembler timeout — use `ProcessPoolExecutor`.
- **DO NOT** pass SMILES through `MolFromSmiles → MolToSmiles` in `parse_inline_string()`.
- **DO NOT** call `CIPAssigner.assign_all()` on an unsanitized mol.
- **DO NOT** allow `PseudoAtomStrategy` wildcard atoms (`*`) to appear in final OIN string.
- **DO NOT** assert `@`/`@@` on BINAP P atoms in tests.
- **DO NOT** implement `_molassembler_worker` as a class method or lambda.

---

## 7. Risks & Mitigation
- **Risk: Molassembler Python bindings not picklable** → **Mitigation:** `MolassemblerSpike` confirms picklability via `pickle.dumps(worker_fn)`.
- **Risk: SCINE build complexity on CI** → **Mitigation:** Pin `>=2.0.0`; add to CI via conda-forge. Document build instructions in `README`.
- **Risk: `PseudoAtomStrategy` wildcard (`*`) conflicts** → **Mitigation:** Strip wildcard atoms before final OIN serialization.

---

## 8. Success Metrics
- All 6 existing verified complexes pass OIN stability check after refactor.
- P-chiral and N-chiral fixtures pass the RDKit CIP oracle test.
- Molassembler conformer generated within 60s for all existing fixtures.
- `oin-smiles xyz2oin tests/fixtures/cisplatin.xyz` returns correct OIN string from CLI.
- `MolassemblerTimeoutError` correctly raised on artificial infinite worker.
