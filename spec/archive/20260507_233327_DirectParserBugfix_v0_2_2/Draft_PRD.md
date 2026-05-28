# Draft SuperPRD: Direct Parser Bugfix & Coexistence Backend (v0.2.2)

**Version:** 0.0.1-draft
**Target Release:** v0.2.2
**Audit Reference:** `spec/audit/DirectParser_IntegrationAudit_20260506.md`
**Supersedes:** `spec/archive/MiniPRD_DirectParser_Integration_v0.2.1_DEFERRED.md`
**Date:** 2026-05-07

---

## 1. Introduction & Goals

### Problem Statement
The `parse_oin_direct()` function in `src/oinsmiles/generation/engine.py` was specified by `MiniPRD_DirectParser_Integration.md` as the new OIN→XYZ pipeline replacing the legacy `MolassemblerAdapter`. However, the May 6, 2026 audit revealed it is **non-functional** — fragment-rank-to-atom-index mapping is missing, polydentate ligands fail, isomerism selection is absent, and eta-bond translation is broken. Five distinct integration blockers prevent its use in production. The legacy adapter remains the working pipeline and ships with v0.2.1.

### Solution Overview
Refactor `parse_oin_direct()` into a fully working pipeline by sequentially fixing the 5 blockers identified in the audit. Introduce a `backend` parameter on `OIN3DGenerator` so both pipelines coexist: `"direct"` (new default) and `"legacy"` (preserved escape hatch). Both backends must pass `verify_roundtrip.py` test parity. Legacy backend retains all current functionality but gets fail-soft CI status (warnings, not blockers) to allow graceful drift in v0.3.0+.

### Target Audience
- **Primary:** Internal OIN-SMILES library development team
- **Secondary:** Library users who depend on `OIN3DGenerator` (current API contract preserved)
- **Tertiary:** Future contributors extending OIN→XYZ to new ligand types (cleaner pipeline = easier extension)

---

## 2. Confidence Mandate

**Confidence Score:** 8/10

### Rationale
- ✅ Audit document provides exhaustive technical detail for all 5 blockers
- ✅ Legacy adapter implementation provides reference for correct mapping logic
- ✅ Test parity bar is well-defined (existing `verify_roundtrip.py`)
- ✅ Sequential MiniPRD execution prevents big-bang failures
- ⚠️ Cross-backend XYZ comparison may surface unexpected geometry differences (mitigated by 1.0Å threshold)
- ⚠️ Polydentate ligand mapping (Blocker 2) requires understanding of vector data structure — may need design iteration
- ⚠️ Permutation selection (Blocker 3) is the most complex blocker; legacy `_pick_masm_permutation()` logic may not transplant cleanly

### Open Clarifying Questions
1. **Q-OPEN-1:** What's the canonical "metal" identification rule when `fragments[0]` is not the metal? (audit assumes always position 0 — does this hold?)
2. **Q-OPEN-2:** Should `parse_oin_direct()` accept v2.4 sidecar format, or only v3.6 inline? (audit examples show v3.6 only)
3. **Q-OPEN-3:** When `backend="direct"` fails on a complex the legacy adapter handles, is that a P0 bug or acceptable scope limitation for v0.2.2?

---

## 3. Scope

### In-Scope
- Rewrite `extract_oin_constraints()` with new return signature including fragment→atom mapping
- Rewrite `parse_oin_direct()` to use atom-indexed constraints and proper SMILES connectivity
- Implement polydentate ligand handling (multi-atom binding)
- Port `_pick_masm_permutation()` logic from `MolassemblerAdapter` into direct parser
- Implement eta-bond atom-index translation
- Add comprehensive test suite for direct parser (unit + integration)
- Add `backend` parameter to `OIN3DGenerator.__init__()` and `OIN3DGenerator.generate()`
- Re-export `parse_oin_direct` from `oinsmiles.__init__.py` as public utility
- Update CI to run both backends with legacy as fail-soft
- Update CHANGELOG with breaking-change notice for `extract_oin_constraints()` signature
- Bump version to v0.2.2

### Out-of-Scope
- Removing legacy `MolassemblerAdapter` (deferred to v0.3.0+ per Q1)
- Novel geometry types not currently supported by legacy adapter (per Q9)
- v2.4 sidecar format support in direct parser (Q-OPEN-2; assume v3.6 only)
- Performance optimization beyond what audit identifies as blocker
- Migration of `OINInlineHandler` or other unrelated parsers
- Documentation overhaul beyond CHANGELOG and inline docstrings
- Bit-identical XYZ output between backends (Q8: 1.0Å tolerance accepted)

---

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
|---|---|---|---|
| US-001 | As `OIN3DGenerator`, I want a `backend` parameter accepting `"direct"` or `"legacy"` so users can choose the pipeline | 1. Constructor accepts `backend` kwarg; `"direct"` is default<br>2. `generate()` accepts `backend` kwarg overriding constructor<br>3. Invalid value raises `ValueError("backend must be 'direct' or 'legacy'")`<br>4. Existing callers without backend kwarg get direct pipeline | High |
| US-002 | As `extract_oin_constraints()`, I want to return a fragment→atom-index mapping so downstream stages can translate vertex_indices correctly | 1. Return signature: `(stripped_smiles, constraints, fragment_to_atom_mapping)`<br>2. `fragment_to_atom_mapping[frag_rank] = [atom_indices]` after connected SMILES is built<br>3. All 17 test call sites updated to handle new tuple<br>4. Documented as breaking change in CHANGELOG | High |
| US-003 | As `parse_oin_direct()`, I want to handle polydentate ligands (chelates) so complexes like Pt(en)Cl₂ work correctly | 1. Multi-atom ligands (en, phen, BINAP, BDPP) generate valid Molassembler molecules<br>2. All `tests/fixtures/oin/PdCl2-RR-BD*.json` round-trip successfully<br>3. RMSD < 1.0Å vs. input XYZ | High |
| US-004 | As `parse_oin_direct()`, I want correct isomerism (cis/trans, fac/mer) so symmetric complexes produce the right geometry | 1. Cisplatin produces cis configuration; Transplatin produces trans<br>2. fac-Ir(ppy)₃ ≠ mer-Ir(ppy)₃<br>3. Permutation index selection logic ported from legacy `_pick_masm_permutation()`<br>4. RMSD < 1.0Å vs. legacy adapter on all symmetric fixtures | High |
| US-005 | As `parse_oin_direct()`, I want eta-ligand vertex indices translated to atom indices so ferrocene and TiCp₂Me₂ work | 1. η⁵-cyclopentadienyl ligands generate correct sandwich geometry<br>2. Ferrocene round-trips with `is_eta=True` bonds preserved<br>3. TiCat1/3/4 fixtures pass<br>4. No "bond to self" errors (tracking issue from audit) | High |
| US-006 | As CI/test infrastructure, I want both backends tested on `verify_roundtrip.py` so parity is enforced | 1. CI runs `verify_roundtrip.py --backend direct` (blocking)<br>2. CI runs `verify_roundtrip.py --backend legacy` (fail-soft / warning)<br>3. New test file `tests/unit/test_parse_oin_direct.py` with ≥8 unit tests covering the 9 geometry codes<br>4. Cross-backend RMSD comparison < 1.0Å for all current fixtures | High |
| US-007 | As `oinsmiles` package user, I want `parse_oin_direct` importable from the top-level so I can use it as a low-level utility | 1. `from oinsmiles import parse_oin_direct` works<br>2. Listed in `__all__` in `__init__.py`<br>3. Docstring documents it as low-level direct parser | Medium |

---

## 5. Technical Specifications

### Architecture & Resolved Trade-offs

#### Resolved Decisions
| ID | Decision | Rationale |
|---|---|---|
| Q1 | Coexistence (B): `backend="direct"\|"legacy"`, direct default | Permanent reversibility; legacy as escape hatch |
| Q2 | Five sequential MiniPRDs (B) | Foundational blocker (#1) must land first; each blocker independently auditable |
| Q3 | Functional parity (A): pass `verify_roundtrip.py` | Existing test suite is canonical; no new bar invented |
| Q4 | Backend on both constructor + method (C) | Stateful default + per-call override |
| Q5 | Break `extract_oin_constraints()` signature (A) | Only 1 production caller (will be rewritten); 17 test call sites updatable |
| Q6 | Public utility (A): re-exported from `oinsmiles` | Power users gain low-level access; advertised as stable |
| Q7 | Both backends CI, legacy fail-soft (C) | Pragmatic; reflects "legacy on its way out" reality |
| Q8 | Cross-backend RMSD < 1.0Å (A) | Matches existing test threshold; avoids overspec |
| Q9 | No novel geometries in v0.2.2 (A) | Tight scope; defer extension to v0.3.0+ |

#### Pipeline Architecture (Target State)

```
User Input (OIN-SMILES)
         │
         ▼
┌────────────────────────┐
│  OIN3DGenerator        │
│  (backend="direct" |   │
│   backend="legacy")    │
└──────────┬─────────────┘
           │
    ┌──────┴──────┐
    │             │
    ▼             ▼
┌─────────┐   ┌──────────────────────┐
│ DIRECT  │   │ LEGACY               │
│ (new)   │   │ (preserved)          │
└────┬────┘   └──────────┬───────────┘
     │                   │
     ▼                   ▼
parse_oin_direct()  OINParser.parse() →
                    MolassemblerAdapter.generate()
     │                   │
     ▼                   ▼
1. extract_oin_constraints()
   → (smiles, constraints, frag_to_atom)
2. _build_connected_smiles_with_mapping()
3. tokenize_unsanitized_smiles()
4. _translate_vertex_indices()
   (frag_rank → atom_idx via mapping)
5. _select_permutation()
   (ported from legacy)
6. construct_molassembler_mol()
   (now receives atom-indexed constraints)
7. _dg_worker() (timeout 10s)
     │                   │
     └──────────┬────────┘
                ▼
        GeneratedStructure(xyz, mol)
```

### System Graph Blast Radius

#### Modified Nodes (architecture.yml)
- `atom_oin3d_generator` (status: dirty) — adds `backend` param to `__init__` and `generate()`
- `atom_direct_parser_regex` (status: dirty) — return signature changes from 2-tuple to 3-tuple
- `atom_direct_parser_masm` (status: dirty) — input format changes (atom indices, not fragment ranks)
- `atom_direct_parser_integration` (status: dirty) — full rewrite of `parse_oin_direct()`
- `atom_api_exports` (status: dirty) — new `parse_oin_direct` export

#### New Nodes (to be added)
- `atom_direct_parser_polydentate` — polydentate ligand connectivity logic (Blocker 2)
- `atom_direct_parser_permutation` — permutation selection logic (Blocker 3)

#### Unchanged Nodes (review only)
- `atom_molassembler_adapter` — preserved as `legacy_backend`; no code changes
- `atom_oin_parser_gen` — `OINParser.parse()` still used by legacy backend
- `atom_generated_structure` — `GeneratedStructure` dataclass unchanged
- `atom_smiles_to_xyz` — `SMILESToXYZ` API unchanged

### Execution Checklist (Sequential MiniPRDs)

The audit's "Optimal Implementation Order" is **#1 → #3 → #4 → #2 → #5**. We follow that order:

1. **MiniPRD_DirectParser_FragmentMapping_v0.2.2.md** (Blocker #1, P0)
   - Add fragment→atom-index mapping in `extract_oin_constraints()` and connected-SMILES builder
   - Update return signatures
   - Foundation for all other work
   - **Estimated effort:** 2-3 hours

2. **MiniPRD_DirectParser_Permutation_v0.2.2.md** (Blocker #3, P0)
   - Port `_pick_masm_permutation()` logic to direct parser
   - Map OIN geometry → Molassembler permutation index
   - Handle symmetric cases (octahedral, TBP)
   - **Estimated effort:** 6-8 hours

3. **MiniPRD_DirectParser_EtaBonds_v0.2.2.md** (Blocker #4, P1)
   - Translate eta-bond vertex indices to atom indices using mapping from #1
   - Handle multi-atom eta bonds (η⁵ cyclopentadienyl)
   - Use `TEMPLATES` data for slot resolution
   - **Estimated effort:** 4-6 hours

4. **MiniPRD_DirectParser_Polydentate_v0.2.2.md** (Blocker #2, P1)
   - Multi-atom ligand binding (en, phen, diphosphines)
   - Multiple bonds from metal to each binding atom
   - Use vector data to determine binding atoms
   - **Estimated effort:** 4-6 hours

5. **MiniPRD_DirectParser_Tests_v0.2.2.md** (Blocker #5, P2)
   - Create `tests/unit/test_parse_oin_direct.py` with ≥8 unit tests
   - Add backend parameter to `verify_roundtrip.py`
   - CI integration: direct blocking, legacy fail-soft
   - Backend selection wiring in `OIN3DGenerator` (US-001)
   - Public API export (US-007)
   - **Estimated effort:** 2-3 hours

**Total estimated effort:** ~20-26 hours (~3-4 working days focused)

### API Contracts / Schema

#### Public API (oinsmiles.__init__.py)
```python
from oinsmiles import (
    XYZToSMILES,        # unchanged
    SMILESToXYZ,        # unchanged
    parse_oin_direct,   # NEW (Q6=A)
)
```

#### OIN3DGenerator (engine.py)
```python
class OIN3DGenerator:
    def __init__(
        self,
        timeout: int = 60,
        dg_strategy: str = "single",
        ensemble_size: int = 10,
        backend: Literal["direct", "legacy"] = "direct",  # NEW
    ) -> None: ...

    def generate(
        self,
        oin_string: str,
        backend: Optional[Literal["direct", "legacy"]] = None,  # NEW
    ) -> GeneratedStructure: ...
```

#### Free Function (engine.py)
```python
def parse_oin_direct(oin_smiles: str) -> GeneratedStructure:
    """Low-level direct parser for OIN-SMILES → 3D structure.
    Uses regex → AST → Molassembler → DG pipeline.
    For most use cases, prefer OIN3DGenerator."""
```

#### Breaking Change: extract_oin_constraints()
```python
# Before (v0.2.1):
def extract_oin_constraints(oin_smiles: str) -> Tuple[str, Dict[int, Dict]]:
    return (stripped_smiles, constraints)

# After (v0.2.2 — BREAKING):
def extract_oin_constraints(oin_smiles: str) -> Tuple[str, Dict[int, Dict], Dict[int, List[int]]]:
    return (stripped_smiles, constraints, fragment_to_atom_mapping)
```

CHANGELOG entry required. Migration: any external caller must accept the third return value (or use `_, _, _ = extract_oin_constraints(...)` pattern).

### Dependencies
- **No new dependencies.** All work uses existing libraries:
  - `scine_molassembler >= 2.0.0` (already pinned)
  - `rdkit` (already pinned)
  - `numpy` (already pinned)
- **Python:** ≥3.10 (already required)

---

## 6. Negative Constraints

- **DO NOT** modify `MolassemblerAdapter` or any function in `molassembler_adapter.py` (legacy backend must remain byte-for-byte stable)
- **DO NOT** delete `_stitch_multi_eta_fragment()`, `_template_generate()`, `_stitch_fragment()`, or `_embed_fragment()` (deferred to v0.3.0+)
- **DO NOT** change `GeneratedStructure` dataclass signature (v0.2.0 API contract)
- **DO NOT** change `OIN3DGenerator.generate()` return type or semantics (must remain `GeneratedStructure`)
- **DO NOT** introduce a third backend or pluggable backend registry (out of scope; YAGNI)
- **DO NOT** add v2.4 sidecar format support to direct parser (out of scope; v3.6 inline only)
- **DO NOT** silently fall back from direct to legacy on error (Q1=B; explicit backend selection only)
- **DO NOT** modify `OINParser.parse()` signature or behavior (used by legacy backend)
- **DO NOT** make `parse_oin_direct()` accept anything other than a string (no Mol objects, no parsed dataclasses)
- **DO NOT** introduce dynamic backend dispatch via env vars or config files (constructor + method kwargs only)
- **DO NOT** ship any MiniPRD with failing `verify_roundtrip.py --backend direct` tests
- **DO NOT** mark legacy backend as deprecated in v0.2.2 (deprecation messaging deferred to v0.3.0+)

---

## 7. Risks & Mitigation

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| Polydentate ligand mapping requires architectural change to vector data | Medium | High | Spike during MiniPRD #4; if blocked, escalate; legacy backend unaffected |
| Permutation selection logic doesn't transplant cleanly from legacy adapter | Medium | High | Read `_pick_masm_permutation()` end-to-end before MiniPRD #2; document required adaptations |
| Cross-backend RMSD exceeds 1.0Å on edge cases (e.g., ansa-metallocenes) | Medium | Medium | Per-fixture exclusion list with documented justification; revisit threshold per Q-OPEN-3 |
| Breaking `extract_oin_constraints()` signature breaks downstream consumers | Low | Medium | Confirmed only 1 production + 17 test call sites; CHANGELOG entry; semver-major-style notice |
| MiniPRD #2 (permutation) blocks #3 and #4 due to dependency | High | Medium | Explicit dependency declared in MiniPRD blocked-by; sequential execution prevents parallel-merge issues |
| Direct backend produces "valid but different" XYZ that fails fixture comparison | Medium | Medium | Use loose threshold (1.0Å); add fixture-specific overrides if needed |
| User sets `backend="direct"` but feature fails silently on unsupported geometry | Low | High | Explicit `NotImplementedError` with message naming the geometry and suggesting `backend="legacy"` |
| CI pipeline doesn't support fail-soft for legacy backend | Low | Low | Use `continue-on-error: true` in GitHub Actions; document in CI config |
| Audit assumed `fragments[0]` is always the metal — may not hold | Low | High | Add explicit metal-detection logic; raise clear error if no metal in fragment 0 (Q-OPEN-1) |

---

## 8. Success Metrics

### Acceptance (Hard Gates)
- ✅ All 5 sequential MiniPRDs land with green audits
- ✅ `verify_roundtrip.py --backend direct` passes for all current fixtures (Cisplatin, Transplatin, Cis-PtCl2(en), Ferrocene, fac/mer-Ir(ppy)₃, PdCl2-R-BINAP, PdCl2-RR-BDNN, PdCl2-RR-BDPP, TiCat1/3/4)
- ✅ `verify_roundtrip.py --backend legacy` passes (continues current state)
- ✅ Cross-backend RMSD < 1.0Å on all current fixtures
- ✅ Unit test coverage ≥80% for `parse_oin_direct()` and helpers
- ✅ `OIN3DGenerator(backend="direct")` and `OIN3DGenerator(backend="legacy")` both work
- ✅ Public API: `from oinsmiles import parse_oin_direct` succeeds
- ✅ CHANGELOG documents breaking change to `extract_oin_constraints()`

### Quality (Soft Goals)
- 📊 Code in `parse_oin_direct()` is < 200 lines (legacy adapter is ~600; we should be cleaner)
- 📊 No new lint errors or type-check failures introduced
- 📊 Doctests in module pass
- 📊 Each MiniPRD's audit completes within ≤2 cycles (no major rework)

### v0.3.0 Roadmap Inputs
- 📈 If direct backend handles ≥95% of fixtures with RMSD < 0.5Å, plan legacy removal for v0.3.0
- 📈 If polydentate handling has fewer than 3 fixture-specific overrides, declare architecture sound
- 📈 Document any geometry types where direct < legacy parity for v0.3.0 follow-up

---

## Appendix A: Audit Document Cross-Reference

This Draft PRD operationalizes the recommendations in:
**`spec/audit/DirectParser_IntegrationAudit_20260506.md`** (855 lines)

Direct mappings:
- Audit "Blockers 1–5" → MiniPRDs in Execution Checklist (Section 5)
- Audit "Option C" → This entire SuperPRD
- Audit "Recommended Path" → Sequential execution policy (Q2=B)
- Audit "Glossary" → Reused for shared terminology
- Audit "Example 1: Mapping Problem Illustrated" → Primary technical reference for MiniPRD #1
- Audit "Example 3: Legacy Adapter Comparison" → Primary technical reference for MiniPRD #2

---

## Appendix B: Open Questions for Red Team

The following Q-OPEN items should be specifically interrogated during `/hyper-redteam`:

1. **Q-OPEN-1:** Metal detection assumption (`fragments[0]` is metal) — robustness?
2. **Q-OPEN-2:** v2.4 sidecar format support in direct parser — explicit error or silent passthrough to legacy?
3. **Q-OPEN-3:** When direct backend fails on a complex legacy handles, is that P0 or scope-limit? (Today's answer: P0; needs adversarial review)
4. **Coexistence trap:** Does Q1=B (permanent coexistence) introduce a maintenance burden the team will regret? Should we pre-commit to v0.3.0 deletion?
5. **Test fairness:** Q7=C (legacy fail-soft) — is this an admission that we expect legacy to break, or a pragmatic concession?
6. **Permutation portability:** Q3 functional parity assumes `_pick_masm_permutation()` ports cleanly. What if it doesn't?
7. **Public API stability:** Q6=A re-exports `parse_oin_direct`. If the implementation changes for v0.3.0, is that a breaking change?

---

**End of Draft PRD.**

**Next Step:** Start a new conversation and run `/hyper-redteam` to perform adversarial analysis on this Draft PRD.
