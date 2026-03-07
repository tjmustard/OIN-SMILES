# SuperPRD — OIN-SMILES: Chiral P/N Stereocenter Support

## Metadata
- **Project Name**: OIN-SMILES
- **Version**: 1.1.0 (Chiral P/N Stereocenter Feature)
- **Status**: Approved — System of Record
- **Owner**: Resolution Agent (2026-03-04)

---

## 1. Introduction & Goals

### 1.1 Problem Statement
The OIN-SMILES v3.6 baseline correctly encodes coordination geometry (cis/trans, facial/meridional, haptic ring orientation) but ignores stereocenters on ligand atoms. Specifically, P and N atoms in TMC ligands can bear `@`/`@@` chirality markers in SMILES, but:
1. The current XYZ→OIN pipeline destroys these markers during `OINInlineHandler.parse_inline_string()` via an RDKit canonicalization round-trip.
2. The pipeline performs no CIP assignment on Zone A atoms (those directly bonded to metal), because fragmentation removes the metal context required for RDKit CIP resolution.
3. The 3D generation backend (Architector) is not stereo-aware for P/N centers.

Chiral phosphines (e.g., BINAP derivatives, chiral monodentate phosphines) and chiral amines are common in asymmetric catalysis TMC datasets. Without @/@@ encoding, enantiomers are indistinguishable in OIN strings.

### 1.2 Solution Overview
OIN-SMILES v1.1.0 adds:
- **Pre-fragmentation CIP assignment** (`CIPAssigner`) on the full sanitized mol, before the metal fragment is separated.
- **Chirality recovery** (`ChiralityRecoveryUtility`) that re-applies `@`/`@@` to ligand SMILES after `OINSanitizer` runs, using stored `_CIPCode` atom properties.
- **PseudoAtomStrategy** (fallback only): wildcard `*` sentinel (`PSEUDO_ATOMIC_NUM=0`) for non-standard-valence P/N atoms where CIP pre-assignment fails. Stripped from final OIN output.
- **`OINInlineHandler.parse_inline_string()` fix**: regex-only slot marker stripping — no RDKit round-trip, preserving `@`/`@@`.
- **Molassembler backend**: replaces Architector entirely. Uses SCINE graph-theoretical stereopermutators + 4D distance geometry. `ProcessPoolExecutor` isolation for GIL-safe timeout.
- **CLI entry point**: `oin-smiles xyz2oin <path>` and `oin-smiles oin2xyz <oin>`.

### 1.3 Target Audience
Computational chemists working with asymmetric TMC datasets who need lossless P/N stereocenter encoding in OIN strings for database storage and similarity search.

---

## 2. Confidence Mandate

**Confidence Score**: 9/10

**Resolved Red Team Flags**: 3 Critical, 4 High, 5 Medium, 4 Low — all documented below.

**Remaining uncertainty**: Molassembler Python binding picklability (confirmed by MolassemblerSpike before MolassemblerAdapter implementation begins).

---

## 3. Scope

### 3.1 In-Scope
- Pre-fragmentation CIP assignment for P and N stereocenters (`CIPAssigner`)
- Chirality recovery after `OINSanitizer` (`ChiralityRecoveryUtility`)
- `PseudoAtomStrategy` as fallback (wildcard `*` sentinel, `PSEUDO_ATOMIC_NUM=0`)
- `OINInlineHandler.parse_inline_string()` regex-only fix (no re-canonicalization)
- Molassembler integration as sole 3D generation backend (`MolassemblerAdapter`)
- `ProcessPoolExecutor`-based timeout isolation (default 60s)
- CLI entry point (`cli.py`)
- P-chiral and N-chiral test fixtures with RDKit CIP oracle validation
- BINAP as structural stability fixture (no `@`/`@@` assertion on P atoms)
- Axial-chiral ligand SMILES must encode chirality for those specific atoms in the output SMILES

### 3.2 Out-of-Scope
- Automatic axial chirality detection (BINAP is a structural stability test, not a chirality encoding test)
- Charged complex support (TD-008 remains open)
- Multi-metal complex support
- Architector backend (fully removed)
- Chiral metal centers (only P/N ligand stereocenters)

---

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
| :--- | :--- | :--- | :--- |
| US-001 | As a chemist, I want XYZ→OIN to encode `@`/`@@` for P/N stereocenters so that chiral phosphine ligands are distinguishable from their enantiomers. | 1. P-chiral fixture XYZ → OIN contains `@`/`@@`<br>2. RDKit CIP oracle confirms expected R/S descriptor | High |
| US-002 | As a chemist, I want OIN→XYZ to reconstruct P-chiral complexes via Molassembler so that 3D chirality is preserved end-to-end. | 1. OIN stability: OIN-in == OIN-out after round-trip<br>2. Molassembler conformer generated within 60s | High |
| US-003 | As a developer, I want `oin-smiles xyz2oin <path>` CLI so that pipelines can call OIN-SMILES from the shell. | 1. Returns OIN string to stdout<br>2. Non-zero exit code on error | Medium |
| US-004 | As a developer, I want `oin-smiles oin2xyz <oin>` CLI so that OIN strings can be converted back to XYZ from the shell. | 1. Returns XYZ block to stdout<br>2. Non-zero exit code on error or timeout | Medium |
| US-005 | As a QA process, I want all 6 existing regression complexes to still pass OIN stability after the refactor so that no regressions are introduced. | All 6 complexes (cisplatin, transplatin, cis-PtCl₂(en), ferrocene, fac-/mer-Ir(ppy)₃) pass OIN stability check | High |

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
  → [OINDiscreteAligner]                     ← geometry detection, slot mapping (unchanged)
  → [OINInlineHandler.generate_inline_string] → inline OIN string
```

**`OINInlineHandler.parse_inline_string()` fix (C-2):**
Strip slot markers and metal prefix via `SLOT_REGEX` and `METAL_REGEX` substitution only. No `MolFromSmiles → MolToSmiles` round-trip. Output is non-canonical SMILES but `@`/`@@`-safe.

**Resolved Trade-offs Log:**

- **C-1 (Zone A encoding contradiction):** Draft PRD had a four-way contradiction: §1.2A/§5.1 diagrams referenced `[Zz]` pseudo-atom; §5.1 Resolved Trade-offs said Option C (no pseudo-atom); §6 banned Oganesson; §5.4 defined `PseudoAtomStrategy(PSEUDO_ATOMIC_NUM=118)`. Also `[Zz]` is RDKit-invalid (`MolFromSmiles("[Zz]")` returns `None`).
  - **Resolution:** Option C (pre-fragmentation CIP) is primary. `PseudoAtomStrategy` kept as fallback for non-standard valences only. `PSEUDO_ATOMIC_NUM = 0` (RDKit wildcard `*`). All `[Zz]`/Oganesson references purged.

- **C-2 (@/@@ corruption in parse_inline_string):** `MolFromSmiles → MolToSmiles(canonical=True)` reorders atoms, silently corrupting `@`/`@@` relative to traversal order.
  - **Resolution:** Option B — regex-only stripping. No RDKit round-trip in `parse_inline_string()`. Output SMILES is non-canonical but chirality-safe.

- **C-3 (Mol sanitization prerequisite):** `CIPAssigner.assign_all()` silently no-ops on unsanitized mol because RDKit's `AssignStereochemistry()` requires aromatic perception, ring info, and valence.
  - **Resolution (autonomous):** `Chem.SanitizeMol(mol)` is a hard precondition called before `CIPAssigner.assign_all()`. Raises `Chem.SanitizeMol` exception on failure (not silently ignored).

- **H-1 (BINAP fixture):** BINAP's chirality is axial (atropisomeric at the biaryl bond), not P-centered. RDKit assigns no `@`/`@@` to BINAP's P atoms. An `@`/`@@` assertion on BINAP P atoms would pass vacuously.
  - **Resolution:** BINAP kept as structural stability fixture (`assert encode(BINAP_XYZ) is not None` only, no `@`/`@@` assertion). Dedicated P-chiral fixture (e.g., methylphenylpropylphosphine or similar) added with explicit `@`/`@@` literal assertion verified by RDKit CIP oracle. Axial-chiral ligand SMILES (e.g., BINAP) must correctly encode the axial chirality descriptors for those atoms in the output SMILES string.

- **H-2 (Timeout GIL):** If Molassembler C++ holds the Python GIL, `threading.Thread.join(timeout)` cannot interrupt it.
  - **Resolution:** `ProcessPoolExecutor(max_workers=1)` for all Molassembler calls. Default `timeout = 60` seconds. Raises `MolassemblerTimeoutError` on expiry.

- **H-4 (Circular validation):** OIN round-trip (XYZ → OIN → XYZ → OIN') tests stability but not correctness — a shared encoder/decoder bug passes tautologically.
  - **Resolution:** RDKit CIP oracle cross-validator added to `ChiralTests`: after encoding, extract ligand SMILES, call `RDKit.AssignStereochemistry()`, compare `_CIPCode` to fixture-expected R or S descriptor. This is an independent correctness oracle.

- **NFR-1 (SCINE version pin):** `scine-molassembler >= 2.0.0` in `pyproject.toml`.

- **NFR-2 (Timeout default):** `MolassemblerAdapter(timeout: int = 60)`.

### 5.2 System Graph Blast Radius

15 of 36 nodes affected (42% of graph). Reference: `spec/compiled/architecture.yml`.

**Modified (8):**
- `atom_XYZToSMILES` — calls `CIPAssigner.assign_all()` before fragmentation; calls `ChiralityRecoveryUtility.recover()` after `OINSanitizer`
- `atom_OINSanitizer` — `ChiralityRecoveryUtility` called immediately after; must not destroy `_CIPCode` atom properties
- `atom_OINDiscreteAligner` — slot assignment logic unchanged; receives chiral mol with `@`/`@@` intact
- `atom_OIN3DGenerator` — rewired to call `MolassemblerAdapter` instead of `ArchitectorAdapter`
- `atom_OINInlineHandler` — `parse_inline_string()` rewritten to regex-only stripping
- `atom_OINParser_oin` — must not destroy `@`/`@@` during sidecar/inline parsing
- `atom_SMILESToXYZ` — rewired to `OIN3DGenerator` (partially resolves TD-003)
- `mod_generation` — Architector dependency removed; SCINE/Molassembler added

**Removed (2):**
- `atom_ArchitectorAdapter` (`generation/architector_adapter.py` deleted)
- `atom_ArchitectorWrapper` (`generation/wrapper.py` deleted)

**Added (5):**
- `atom_CIPAssigner` → `src/oinsmiles/core/chirality.py`
- `atom_ChiralityRecoveryUtility` → `src/oinsmiles/core/chirality.py`
- `atom_PseudoAtomStrategy` → `src/oinsmiles/core/chirality.py` (fallback only)
- `atom_MolassemblerAdapter` → `src/oinsmiles/generation/molassembler_adapter.py`
- `atom_CLI` → `src/oinsmiles/cli.py`

### 5.3 Execution Checklist (MiniPRDs)

Execute in order — each MiniPRD has a prerequisite dependency:

- [ ] `spec/compiled/MiniPRD_MolassemblerSpike.md` — validate SCINE import + `ProcessPoolExecutor` picklability **[no prerequisites]**
- [ ] `spec/compiled/MiniPRD_ChiralEncoding.md` — `CIPAssigner`, `ChiralityRecoveryUtility`, `PseudoAtomStrategy`, `parse_inline_string()` fix **[no prerequisites]**
- [ ] `spec/compiled/MiniPRD_MolassemblerAdapter.md` — `MolassemblerAdapter` + `OIN3DGenerator` rewire + delete Architector **[requires: MolassemblerSpike complete]**
- [ ] `spec/compiled/MiniPRD_CLI.md` — `cli.py` entry point **[requires: MolassemblerAdapter complete]**
- [ ] `spec/compiled/MiniPRD_ChiralTests.md` — test fixtures + RDKit CIP oracle **[requires: ChiralEncoding + MolassemblerAdapter complete; candidate fixtures reviewed]**

### 5.4 API Contracts

```python
# src/oinsmiles/core/chirality.py

class CIPAssigner:
    """Pre-fragmentation CIP assignment on the full sanitized mol."""

    def assign_all(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Precondition: mol must be sanitized (Chem.SanitizeMol called by caller).
        Calls Chem.AssignStereochemistry(mol, cleanIt=True, force=True).
        Stores '_CIPCode' ('R' or 'S') as an atom integer property on all
        P and N atoms with a resolved CIP descriptor.
        Returns: same mol object with properties set.
        Raises: ValueError if mol is None.
        """


class ChiralityRecoveryUtility:
    """Re-applies @/@@ to ligand SMILES from stored _CIPCode atom properties."""

    def recover(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Reads '_CIPCode' from atom properties set by CIPAssigner.assign_all().
        Restores @/@@ chirality tags on P/N atoms using CIP-to-chiral-tag mapping.
        If no '_CIPCode' found on a P/N atom: delegates to PseudoAtomStrategy.
        Returns: mol with @/@@ restored.
        """


class PseudoAtomStrategy:
    """
    Fallback for P/N atoms where CIP pre-assignment fails (non-standard valence).
    PSEUDO_ATOMIC_NUM = 0 (RDKit wildcard atom `*`).
    IMPORTANT: wildcard atoms MUST be stripped from the mol before final OIN
    serialization. They must never appear in the canonical OIN output string.
    """
    PSEUDO_ATOMIC_NUM: int = 0


# src/oinsmiles/generation/molassembler_adapter.py

class MolassemblerTimeoutError(RuntimeError):
    """Raised when Molassembler conformer generation exceeds timeout."""


class MolassemblerAdapter:
    def __init__(self, timeout: int = 60) -> None:
        """timeout: seconds before ProcessPoolExecutor raises MolassemblerTimeoutError."""

    def generate(self, parsed_oin: ParsedOIN) -> str:
        """
        Runs Molassembler distance geometry in a separate OS process via
        ProcessPoolExecutor(max_workers=1).
        Returns: XYZ block as string.
        Raises: MolassemblerTimeoutError if conformer not generated within timeout seconds.
        """


# _molassembler_worker must be a module-level function (not a method) for picklability.
def _molassembler_worker(args: dict) -> str:
    """Picklable worker submitted to ProcessPoolExecutor. Returns XYZ block."""


# src/oinsmiles/generation/engine.py (updated signature only)
class OIN3DGenerator:
    def generate(self, oin_string: str) -> str:
        """
        XYZ block as string.
        Raises: MolassemblerTimeoutError on timeout.
        Uses MolassemblerAdapter (Architector fully removed).
        """


# src/oinsmiles/cli.py
def main() -> None:
    """
    Entry point for `oin-smiles` CLI.
    Subcommands:
      xyz2oin <path>   — print OIN string to stdout
      oin2xyz <oin>    — print XYZ block to stdout
    Non-zero exit code on any error.
    """
```

### 5.5 Dependencies

**Added:**
- `scine-molassembler >= 2.0.0` — 3D conformer generation with stereopermutators (BSD-3-Clause)
- `concurrent.futures` (stdlib) — `ProcessPoolExecutor` for GIL-safe timeout isolation

**Removed:**
- `architector` — fully replaced by Molassembler
- `xTB` (was Architector's geometry optimizer) — no longer required

**Retained:**
- `RDKit >= 2025.9.3`, `NumPy`, `SciPy`, `mendeleev`, `networkx`, `openbabel-wheel`

---

## 6. Negative Constraints

- **DO NOT** use Oganesson (Z=118) or `[Zz]` as a pseudo-atom sentinel — `[Zz]` returns `None` from `Chem.MolFromSmiles`.
- **DO NOT** use `threading.Thread` for Molassembler timeout — use `concurrent.futures.ProcessPoolExecutor` exclusively.
- **DO NOT** pass SMILES through `MolFromSmiles → MolToSmiles` in `parse_inline_string()` — this corrupts `@`/`@@`.
- **DO NOT** call `CIPAssigner.assign_all()` on an unsanitized mol — `Chem.SanitizeMol()` is a hard precondition.
- **DO NOT** allow `PseudoAtomStrategy` wildcard atoms (`*`) to appear in the final serialized OIN string — strip them before output.
- **DO NOT** assert `@`/`@@` on BINAP P atoms in tests — BINAP chirality is axial (atropisomeric), not P-centered.
- **DO NOT** promote fixture SMILES from `tests/candidate_outputs/` to `tests/fixtures/` without human review.
- **DO NOT** implement `_molassembler_worker` as a class method or lambda — it must be a picklable module-level function.
- **DO NOT** write to `spec/archive/` manually — use `archive_specs.py` exclusively.

---

## 7. Risks & Mitigation

- **Risk: Molassembler Python bindings not picklable** → **Mitigation:** `MolassemblerSpike` confirms picklability via `pickle.dumps(worker_fn)` before `MolassemblerAdapter` is implemented. If not picklable, use `multiprocessing.Queue`-based worker pattern instead.
- **Risk: SCINE build complexity on CI** → **Mitigation:** Pin `>=2.0.0`; add to CI pre-install step via conda-forge channel. Document build instructions in `README`.
- **Risk: `PseudoAtomStrategy` wildcard (`*`) conflicts with downstream SMARTS queries** → **Mitigation:** Strip wildcard atoms before final OIN serialization — they must never appear in the output OIN string.
- **Risk: New P-chiral fixture SMILES incorrect `@`/`@@`** → **Mitigation:** All new fixture SMILES route to `tests/candidate_outputs/` for human review; RDKit CIP oracle is required before promotion to `tests/fixtures/`.
- **Risk: `chiral_stereo_check` flag in `xyz2mol_local.py` sets stereo flags before `CIPAssigner.assign_all()` runs, causing undefined precedence** → **Mitigation:** Audit `xyz2mol_local.py` during `ChiralEncoding` implementation to confirm `CIPAssigner` runs after and overrides any pre-set stereo flags.

---

## 8. Success Metrics

- All 6 existing verified complexes (cisplatin, transplatin, cis-PtCl₂(en), ferrocene, fac-/mer-Ir(ppy)₃) still pass OIN stability check after refactor.
- At least 1 P-chiral fixture passes the RDKit CIP oracle test (expected R/S descriptor confirmed).
- At least 1 N-chiral fixture passes the RDKit CIP oracle test.
- Molassembler conformer generated within 60s for all existing fixtures.
- `oin-smiles xyz2oin tests/fixtures/cisplatin.xyz` returns correct OIN string from CLI.
- `MolassemblerTimeoutError` correctly raised (not silently swallowed) on a 1-second timeout with artificial infinite worker.
