# SuperPRD: Direct Parser Bugfix & Coexistence Backend (v0.2.2)

## Metadata
- **Project Name:** OIN-SMILES
- **Feature:** Repair `parse_oin_direct()` (5 audit-identified blockers) and ship a `backend` parameter on `OIN3DGenerator` so direct and legacy pipelines coexist
- **Version:** 0.2.2
- **Status:** Resolved — Ready for Implementation
- **Resolution Agent:** hyper-resolve (2026-05-07)
- **Supersedes:** `spec/compiled/SuperPRD_DirectParser.md` (v0.2.1, design-only) and `spec/archive/MiniPRD_DirectParser_Integration_v0.2.1_DEFERRED.md`
- **Audit Reference:** `spec/audit/DirectParser_IntegrationAudit_20260506.md` (855 lines)
- **Cross-Reference:** `spec/compiled/architecture.yml` (last-audited 2026-05-06)

---

## 1. Introduction & Goals

### Problem Statement
The `parse_oin_direct()` function in `src/oinsmiles/generation/engine.py` was specified by `MiniPRD_DirectParser_Integration.md` (v0.2.1) as the new OIN→XYZ pipeline replacing the legacy `MolassemblerAdapter`. The 2026-05-06 audit found it **non-functional** across five distinct integration blockers:

1. **Blocker #1 (P0):** Fragment-rank ↔ atom-index mapping missing → "bond to self" errors.
2. **Blocker #2 (P1):** Polydentate ligand connectivity not handled (en, BINAP, BDPP fail).
3. **Blocker #3 (P0):** Permutation/isomerism selection absent (cis vs. trans, fac vs. mer indistinguishable).
4. **Blocker #4 (P1):** Eta-bond translation to atom indices broken.
5. **Blocker #5 (P2):** Test coverage missing for direct parser.

Production v0.2.1 ships with the legacy `MolassemblerAdapter` exclusively; the direct parser is dead code.

### Solution Overview
Refactor `parse_oin_direct()` into a working pipeline by sequentially fixing each blocker. Introduce a `backend: Literal["direct", "legacy"]` parameter on `OIN3DGenerator` so both pipelines coexist. The `direct` backend is the new default in v0.2.2 (per Resolved Decision D1). The `legacy` backend is preserved as an explicit escape hatch with a defined deletion path in v0.3.0 (per Resolved Decision D7).

Both backends MUST pass `tests/integration/verify_roundtrip.py` for **all current fixtures** as a hard gate. Cross-backend mean RMSD MUST be < 1.0Å on every current fixture (per D2). Novel fixtures introduced post-v0.2.2 may invoke the per-fixture exclusion list; existing fixtures may not.

### Target Audience
- **Primary:** Internal OIN-SMILES library development team.
- **Secondary:** Library users of `OIN3DGenerator` (current API contract preserved; new `backend` kwarg is additive).
- **Tertiary:** Future contributors extending OIN→XYZ to new ligand types — the direct pipeline is intended to be the cleaner extension surface in v0.3.0+.

---

## 2. Confidence Mandate

**Confidence Score:** 9/10 (post-resolution; up from 8/10 in Draft)

**Rationale:**
- ✅ Audit document provides exhaustive technical detail for all 5 blockers.
- ✅ Legacy adapter implementation provides reference for correct mapping logic.
- ✅ All Q-OPEN questions resolved (D3–D7); previously load-bearing assumptions now either gated by spike (D3) or codified as project invariants (metal-first canonical form).
- ✅ Test parity bar is well-defined: `verify_roundtrip.py` mean RMSD < 1.0Å.
- ✅ Sequential MiniPRD execution prevents big-bang failures; each blocker is independently auditable.
- ⚠️ Cross-backend XYZ comparison may surface conformational differences on flexible ligands (mitigated by mean-RMSD threshold per D2).
- ⚠️ Permutation port (Blocker #3) is the most complex; gated by mandatory pre-MiniPRD-#2 write-up of `_pick_masm_permutation()` (per D4).

**Aggregated confidence rule (per A11):** Ship confidence = `min(MiniPRD confidences)`. Any MiniPRD landing below 7/10 triggers Architect re-audit before subsequent MiniPRDs merge.

### Resolved Clarifying Questions
| ID | Question | Resolution |
|---|---|---|
| Q-OPEN-1 | `fragments[0]` always metal? | **Invariant by canonical-form design** (per project memory `project_oin_invariants.md`). Pre-MiniPRD-#1 verification spike runs over all fixtures as belt-and-suspenders. |
| Q-OPEN-2 | v2.4 sidecar support in direct parser? | **No.** Raise `OINFormatError("v3.6 inline only — use OIN3DGenerator(backend='legacy') for v2.4")` at the top of `parse_oin_direct`. |
| Q-OPEN-3 | Direct fails where legacy succeeds — P0 or scope? | **Triage matrix:** known gap (in `_exclusions.yml`) = P2 with clear fallback message; unknown failure = P0 ship-blocker. |
| Q-OPEN-4 | Permanent coexistence? | **No.** `FutureWarning` on `backend='legacy'` starting v0.2.3; legacy removed in v0.3.0 with 1-version deprecation. |
| Q-OPEN-5 | Legacy fail-soft scope? | **Existing fixtures: blocking on both backends.** Fail-soft applies only to novel fixtures introduced after v0.2.2. |
| Q-OPEN-6 | Permutation port cleanly? | **Gated.** 1-page write-up of `_pick_masm_permutation()` (call graph, state deps) is a mandatory deliverable attached to MiniPRD #2. |
| Q-OPEN-7 | `parse_oin_direct` API stability? | **`@experimental` for v0.2.2** via docstring tag + `__experimental__` module attribute. Promoted to stable in v0.3.0. |

---

## 3. Scope

### In-Scope
- Rewrite `_extract_oin_constraints()` (private; renamed from `extract_oin_constraints`) to return fragment→atom mapping.
- Rewrite `parse_oin_direct()` to consume atom-indexed constraints and produce valid Molassembler molecules.
- Implement polydentate ligand handling (multi-atom binding) for en, phen, BINAP, BDPP, BDNN.
- Port `_pick_masm_permutation()` semantics from `MolassemblerAdapter` into the direct parser.
- Implement eta-bond vertex→atom-index translation with multi-eta support (η⁵ Cp).
- Add `backend` parameter to `OIN3DGenerator.__init__()` and `OIN3DGenerator.generate()`.
- Add `random_seed: Optional[int] = None` to `OIN3DGenerator.__init__()` for DG determinism (R7).
- Add `timeout: int = 60` parameter to `parse_oin_direct()` (R2).
- Add `MAX_OIN_STRING_LEN = 100_000` input-size ceiling (R3).
- Defer `scine_molassembler` import to first `OIN3DGenerator` instantiation (R8).
- Re-export `parse_oin_direct` from `oinsmiles.__init__.py` with `@experimental` markings (D6).
- Add `backend` kwarg to `tests/integration/verify_roundtrip.py`; CI runs both with the policy from D7.
- Add `tools/audit_extract_calls.py` (CI-runnable) to enumerate `_extract_oin_constraints` call sites (A8).
- Add `tests/fixtures/oin/_exclusions.yml` schema and governance doc (A12).
- Bump version to v0.2.2 in `pyproject.toml`; pin `scine_molassembler >= 3.0.0, < 4.0.0` (D5).
- Update CHANGELOG.md with the silent-default-flip notice (US-001), the legacy-deprecation roadmap, and the operational backout commitment (R10).

### Out-of-Scope
- Removing legacy `MolassemblerAdapter` (deferred to v0.3.0 per D7).
- Novel geometry types not currently supported by legacy adapter.
- v2.4 sidecar format support in direct parser (Q-OPEN-2: hard-fail by design).
- Performance optimization beyond audit-identified blockers.
- Migration of `OINInlineHandler` or other unrelated parsers.
- Documentation overhaul beyond CHANGELOG and inline docstrings.
- Bit-identical XYZ output between backends (D2: mean RMSD < 1.0Å is acceptance bar).
- Thread-safety of `OIN3DGenerator` (R4: documented as not-thread-safe; no locking work).

---

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
|---|---|---|---|
| **US-001** | As `OIN3DGenerator`, I want a `backend` parameter accepting `"direct"` or `"legacy"` so users can choose the pipeline | 1. Constructor accepts `backend` kwarg; `"direct"` is default<br>2. `generate(backend=…)` overrides constructor when not `None`<br>3. Invalid value → `ValueError("backend must be 'direct' or 'legacy'")`<br>4. Existing callers without backend kwarg get the direct pipeline (silent default-flip; documented in CHANGELOG)<br>5. Typing: `Literal["direct","legacy"] \| None = None` on `generate`; precedence-comment in docstring (A1) | High |
| **US-002** | As `_extract_oin_constraints()`, I want to return a fragment→atom-index mapping so downstream stages can translate `vertex_indices` correctly | 1. New return signature: `(stripped_smiles, constraints, fragment_to_atom_mapping)`<br>2. `fragment_to_atom_mapping[frag_rank] -> List[int]` (atom indices in the connected SMILES)<br>3. Underscore-prefixed (private; was never in `__all__`)<br>4. All call sites updated; live count produced by `tools/audit_extract_calls.py` (A8) | High |
| **US-003** | As `parse_oin_direct()`, I want to handle polydentate ligands so `Pt(en)Cl₂` and similar chelates work | 1. Multi-atom ligands (en, phen, BINAP, BDPP, BDNN) generate valid Molassembler molecules<br>2. All `tests/fixtures/oin/PdCl2-RR-BD*.json` round-trip successfully<br>3. Mean RMSD < 1.0Å vs. input XYZ | High |
| **US-004** | As `parse_oin_direct()`, I want correct isomerism so symmetric complexes produce the right geometry | 1. Cisplatin → cis; Transplatin → trans<br>2. fac-Ir(ppy)₃ ≠ mer-Ir(ppy)₃<br>3. Permutation index logic semantically equivalent to legacy `_pick_masm_permutation()`<br>4. Cross-backend mean RMSD < 1.0Å on all symmetric fixtures<br>5. **Gating deliverable:** `MiniPRD_Permutation` carries a 1-page write-up of legacy permutation logic (D4) | High |
| **US-005** | As `parse_oin_direct()`, I want eta-ligand vertex indices translated to atom indices so ferrocene and TiCp₂Me₂ work | 1. η⁵-cyclopentadienyl ligands generate correct sandwich geometry<br>2. Ferrocene round-trips with `is_eta=True` bonds preserved<br>3. TiCat1/3/4 fixtures pass<br>4. Dedicated regression test `test_no_bond_to_self_after_eta_translation` (A9) pins the audit failure mode | High |
| **US-006** | As CI/test infrastructure, I want both backends tested on `verify_roundtrip.py` so parity is enforced | 1. CI runs `verify_roundtrip.py --backend direct` (blocking)<br>2. CI runs `verify_roundtrip.py --backend legacy` (blocking on **existing** fixtures; fail-soft only on novel post-v0.2.2 fixtures per D7)<br>3. New `tests/unit/test_parse_oin_direct.py` with ≥1 test per geometry code (9 codes) per A5<br>4. Cross-backend mean RMSD < 1.0Å for every current fixture; divergences route to `tests/fixtures/oin/_exclusions.yml` with written justification | High |
| **US-007** | As `oinsmiles` package user, I want `parse_oin_direct` importable as a low-level utility | 1. `from oinsmiles import parse_oin_direct` works<br>2. Listed in `__all__` in `__init__.py` (additions only; existing entries preserved per P6.3)<br>3. Marked `@experimental` via docstring + `__experimental__ = True` module attribute (A3)<br>4. Docstring documents v0.3.0 promotion-to-stable contract | Medium |
| **US-008** | As a user with a v2.4-sidecar OIN string, I want a clear error pointing me to `backend='legacy'` | 1. `parse_oin_direct` detects non-v3.6 input via regex pre-check<br>2. Raises `OINFormatError("v3.6 inline only — use OIN3DGenerator(backend='legacy') for v2.4")` (D6)<br>3. Error message includes the input string truncated to 200 chars | Medium |
| **US-009** | As a user, I want `OIN3DGenerator` to be importable even when `scine_molassembler` import fails (e.g., on an unsupported platform) so XYZ→OIN remains usable | 1. `import oinsmiles` succeeds without Molassembler installed<br>2. `XYZToSMILES` works without Molassembler<br>3. First `OIN3DGenerator(backend='direct')` instantiation raises `ImportError("scine_molassembler is required for OIN→XYZ; install via …")` (R8)<br>4. Test in `tests/unit/test_import_resilience.py` simulates missing Molassembler via `sys.modules` patching | Medium |
| **US-010** | As an operator, I want a deterministic option for DG so reproducibility issues can be debugged | 1. `OIN3DGenerator(random_seed=42)` produces identical XYZ on repeated calls (within FP noise)<br>2. Default `random_seed=None` preserves current stochastic ensemble behavior (R7)<br>3. Documented in `OIN3DGenerator` docstring | Low |

---

## 5. Technical Specifications

### Architecture & Resolved Trade-offs

#### Resolved Decisions Log
| ID | Decision | Rationale |
|---|---|---|
| **D1** | Default backend = `direct` in v0.2.2 | User accepts silent default-flip risk in service of moving forward; mitigated by R10 backout commitment |
| **D2** | Cross-backend acceptance = mean RMSD < 1.0Å (no per-atom cap) | Conformational flexibility (methyl rotations, free torsions) makes per-atom thresholds noise-prone (saved as feedback memory) |
| **D3** | Pre-MiniPRD-#1 verification spike on metal-position assumption | Belt-and-suspenders; the metal-first invariant is by-design but the spike confirms it across all current fixtures |
| **D4** | Pre-MiniPRD-#2 1-page write-up of `_pick_masm_permutation()` | Confidence-evidence; redesign happens in writing, not in code |
| **D5** | Pin `scine_molassembler >= 3.0.0, < 4.0.0` | Matches CI-tested 3.0.1; fail-loud on hypothetical 4.x major-version drift |
| **D6** | `OINFormatError` on v2.4 input + `@experimental` tag on `parse_oin_direct` | Fail-loud at the boundary; signal that the surface is not yet locked |
| **D7** | Time-boxed coexistence: legacy blocking on existing fixtures, fail-soft only on novel; v0.3.0 removal with `FutureWarning` from v0.2.3 | Lifecycle is explicit; engineers and users know the runway |
| Q1 | Coexistence (B): `backend="direct"\|"legacy"`, direct default | Reversibility; legacy as escape hatch (refines D1) |
| Q2 | Five sequential MiniPRDs | Foundational blocker (#1) must land first |
| Q3 | Functional parity: pass `verify_roundtrip.py` | Existing test suite is canonical |
| Q4 | Backend on both constructor + method | Stateful default + per-call override |
| Q5 | Underscore-prefix `_extract_oin_constraints` | Was never `__all__`-exported; not a breaking change |
| Q6 | Public utility: `parse_oin_direct` re-exported with `@experimental` | Power users gain low-level access; surface is signal-locked |
| Q9 | No novel geometries in v0.2.2 | Tight scope; defer extension to v0.3.0+ |

#### Pipeline Architecture (Target State)
```
User Input (OIN-SMILES, v3.6 inline)
         │
         ▼
┌────────────────────────┐
│  OIN3DGenerator        │
│  backend = "direct"|   │
│  "legacy"              │
│  timeout=60            │
│  random_seed=None      │
└──────────┬─────────────┘
           │
   ┌───────┴────────┐
   ▼                ▼
┌─────────┐    ┌──────────────────────┐
│ DIRECT  │    │ LEGACY               │
│ (new)   │    │ (preserved, byte-    │
│         │    │  stable per P6.4)    │
└────┬────┘    └──────────┬───────────┘
     │                    │
     ▼                    ▼
parse_oin_direct(s, t)   OINParser.parse() →
     │                   MolassemblerAdapter.generate()
     │                   │
     ├─ length-check (R3)│
     ├─ v3.6-format-check (D6, US-008)
     ├─ _extract_oin_constraints
     │     → (smiles, constraints, frag_to_atom_map)
     ├─ _build_connected_smiles_with_mapping
     ├─ tokenize_unsanitized_smiles (existing)
     ├─ _translate_vertex_indices (frag_rank → atom_idx)
     ├─ _select_permutation (ported per D4)
     ├─ construct_molassembler_mol (existing, atom-indexed inputs)
     └─ _dg_worker (timeout = min(60-t_elapsed, 30); seed-aware)
                        │
                        ▼
              GeneratedStructure(xyz, mol)
```

### System Graph Blast Radius

#### Modified Nodes (architecture.yml — status: dirty)
- `atom_oin3d_generator` — adds `backend`, `random_seed` params; defers Molassembler import.
- `atom_direct_parser_regex` — return signature 2-tuple → 3-tuple; renames public→private.
- `atom_direct_parser_masm` — input format (atom indices, not fragment ranks).
- `atom_direct_parser_integration` — full rewrite of `parse_oin_direct()`.
- `atom_api_exports` — new `parse_oin_direct` export with `@experimental`.
- `atom_smiles_to_xyz` — set `needs_review` (transitive dep on `atom_oin3d_generator`).

#### New Nodes (to be added)
- `atom_direct_parser_polydentate` — `src/oinsmiles/generation/polydentate.py` (or section in `engine.py` if <100 LOC; per A4).
- `atom_direct_parser_permutation` — `src/oinsmiles/generation/permutation.py` (or section, same rule).
- `atom_oin_format_validator` — input length + v3.6 detection (R3, D6).
- `atom_exclusion_governance` — `tests/fixtures/oin/_exclusions.yml` schema + loader.

#### Unchanged Nodes
- `atom_molassembler_adapter` — preserved verbatim per P6.4.
- `atom_oin_parser_gen` — `OINParser.parse()` still used by legacy backend.
- `atom_generated_structure` — dataclass unchanged; optional fields may be added with defaults (P6.6 clarification).

### Execution Checklist (Sequential MiniPRDs)

Audit order **#1 → #3 → #4 → #2 → #5** (from optimal-implementation-order analysis). MiniPRD filenames append `_v0.2.2` to disambiguate from earlier work.

1. **`MiniPRD_DirectParser_FragmentMapping_v0.2.2.md`** (Blocker #1, P0) — 2–3h
2. **`MiniPRD_DirectParser_Permutation_v0.2.2.md`** (Blocker #3, P0; gated by D4 write-up) — 6–8h
3. **`MiniPRD_DirectParser_EtaBonds_v0.2.2.md`** (Blocker #4, P1) — 4–6h
4. **`MiniPRD_DirectParser_Polydentate_v0.2.2.md`** (Blocker #2, P1) — 4–6h
5. **`MiniPRD_DirectParser_Tests_v0.2.2.md`** (Blocker #5 + US-001/006/007/008/009/010 wiring) — 4–6h

**Total:** ~20–29h with 25% contingency (per P5.4); equates to ~3–4 working days.

### API Contracts

#### `oinsmiles.__init__.py` (additions only, per P6.3)
```python
from oinsmiles import (
    XYZToSMILES,        # unchanged
    SMILESToXYZ,        # unchanged
    parse_oin_direct,   # NEW (Q6=A; @experimental)
)
```

#### `OIN3DGenerator` (engine.py)
```python
class OIN3DGenerator:
    def __init__(
        self,
        timeout: int = 60,
        dg_strategy: str = "single",
        ensemble_size: int = 10,
        backend: Literal["direct", "legacy"] = "direct",       # NEW (D1)
        random_seed: Optional[int] = None,                     # NEW (R7)
    ) -> None: ...

    def generate(
        self,
        oin_string: str,
        backend: Literal["direct", "legacy"] | None = None,    # None → use self._backend (A1)
    ) -> GeneratedStructure: ...
```

#### Free Function (engine.py)
```python
def parse_oin_direct(
    oin_smiles: str,
    timeout: int = 60,                                          # NEW (R2)
) -> GeneratedStructure:
    """Low-level direct parser for OIN-SMILES → 3D structure.

    > **Status:** Experimental — signature stable for v0.2.2;
    > semantics may change in v0.3.0. Marked via `__experimental__`.

    Uses regex → AST → Molassembler → DG pipeline. For most use
    cases, prefer `OIN3DGenerator(backend='direct')`.
    """

__experimental__ = True  # module attribute (A3)
```

#### Private function (was public)
```python
# v0.2.1:
def extract_oin_constraints(oin_smiles: str) -> Tuple[str, Dict[int, Dict]]: ...

# v0.2.2 (PRIVATE, expanded — not a breaking-change announcement):
def _extract_oin_constraints(
    oin_smiles: str,
) -> Tuple[str, Dict[int, Dict], Dict[int, List[int]]]: ...
```

#### Constants
```python
MAX_OIN_STRING_LEN: int = 100_000  # R3 — DoS protection
SUPPORTED_GEOMETRY_CODES: Set[str] = {"SPL", "TBP", "OCT", ...}  # R9
```

### Dependencies
- `scine_molassembler >= 3.0.0, < 4.0.0` (tightened from `>= 2.0.0` per D5)
- `rdkit` (already pinned)
- `numpy` (already pinned)
- **Python:** ≥3.10 (already required)
- **No new runtime dependencies.** `radon` may be added as a dev-only dependency for cyclomatic-complexity check (A6); optional.

---

## 6. Negative Constraints

- **DO NOT** modify `MolassemblerAdapter` or any function in `molassembler_adapter.py` (legacy backend must remain byte-stable per P6.4).
- **DO NOT** delete `_stitch_multi_eta_fragment()`, `_template_generate()`, `_stitch_fragment()`, or `_embed_fragment()` (deferred to v0.3.0+).
- **DO NOT** change `GeneratedStructure` dataclass signature; *adding* optional fields with safe defaults is allowed (P6.6).
- **DO NOT** change `OIN3DGenerator.generate()` return type; must remain `GeneratedStructure`.
- **DO NOT** introduce a third backend or pluggable backend registry (out of scope; YAGNI).
- **DO NOT** add v2.4 sidecar format support to direct parser (Q-OPEN-2; v3.6 inline only).
- **DO NOT** silently fall back from direct to legacy on error (D6/D7; explicit selection only).
- **DO NOT** modify `OINParser.parse()` signature or behavior (legacy depends on it).
- **DO NOT** make `parse_oin_direct()` accept anything other than a string (no `Mol` objects, no dataclasses).
- **DO NOT** introduce dynamic backend dispatch via env vars or config files.
- **DO NOT** ship any MiniPRD with failing `verify_roundtrip.py --backend direct` on existing fixtures.
- **DO NOT** mark legacy as deprecated **in v0.2.2** (deprecation messaging starts v0.2.3 per D7).
- **DO NOT** modify any file outside the `engine.py` / `oin_parser.py` / `__init__.py` / `verify_roundtrip.py` / `tests/` / `polydentate.py` / `permutation.py` / `pyproject.toml` / `CHANGELOG.md` / `_exclusions.yml` perimeter without explicit Architect approval (P6.1).
- **DO NOT** change the default value of `OIN3DGenerator.timeout` (must remain `60`) (P6.2).
- **DO NOT** remove existing entries from `oinsmiles.__all__` (additions only) (P6.3).
- **DO NOT** swallow exceptions from direct backend; all failures must propagate with the input OIN string truncated to 200 chars and a recommendation to try `backend='legacy'` (P6.5).
- **DO NOT** raise the cross-backend RMSD threshold above 1.0Å mean unilaterally; divergences route to `_exclusions.yml` with justification (P7.2).
- **DO NOT** modify `xyz2mol.py` or `core/chirality.py` while in this MiniPRD chain (W6.1; XYZ→OIN side-effects out of scope).

---

## 7. Risks & Mitigation

Probability tiers calibrated per P7.1: **Low ≈ 10%, Medium ≈ 30%, High ≈ 60%**.

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| Polydentate ligand mapping requires architectural change to vector data | Medium | High | Spike during MiniPRD #4; if blocked, escalate to Architect; legacy backend unaffected. |
| Permutation logic doesn't transplant cleanly from legacy adapter | Medium | High | D4 mandate: 1-page write-up read-through *before* MiniPRD #2 starts; redesign happens in prose if needed. |
| Cross-backend mean RMSD exceeds 1.0Å on edge cases (e.g., ansa-metallocenes) | Medium | Medium | Per-fixture exclusion list with written justification per entry (A12); auto-files v0.3.0 issue (P7.6). |
| Privatizing `extract_oin_constraints` (rename to `_extract_oin_constraints`) breaks an unknown external consumer | Low | Low | Underscore prefix communicates "private"; CHANGELOG mentions; not a `__all__` member, so no public-API contract was made. |
| MiniPRD #2 (permutation) blocks #3, #4 in audit order — but execution order #1→#3→#4→#2→#5 places #2 fourth | Low | Medium | The "blocks" claim was inverted in the Draft; execution order reflects the audit's optimal sequencing. Reconciled (P7.4). |
| Direct backend produces "valid but different" XYZ that fails fixture comparison | Medium | Medium | Investigate divergence per fixture; do *not* raise threshold unilaterally (P7.2). |
| User sets `backend='direct'` but feature fails on unsupported geometry | Low | High | Explicit `NotImplementedError(f"Geometry {code} not supported by direct backend; try backend='legacy'")` (R9). |
| CI fail-soft for novel-fixture legacy doesn't behave correctly | Low | Low | Use `continue-on-error: true` keyed by fixture-tag in GitHub Actions matrix; document. |
| Audit assumed `fragments[0]` is always metal — may not hold | Very Low | High | Project invariant per `project_oin_invariants.md`. Pre-MiniPRD-#1 spike confirms across all current fixtures (D3). |
| Hidden risk: `scine_molassembler` C++ extension segfaults during subprocess timeout | Low | Medium | Catch `subprocess.TimeoutExpired` *and* `BrokenProcessPool`; surface both with diagnostic info (Q7.3 (a)). |
| RDKit aromatic-perception drift between minor versions | Low | Medium | Add CI matrix entry pinning current RDKit; failure → investigate before bumping (Q7.3 (b)). |
| Python 3.10 → 3.12 `Literal` typing semantics shift | Low | Low | Test against both 3.10 and 3.12 in CI; type-check with `mypy --python-version 3.10 3.12` (Q7.3 (c)). |
| `ProcessPoolExecutor` deprecation path or CI fork-resource exhaustion | Low | Medium | Per-call executor (R5); document batch-from-outside; no global pool (Q7.3 (d)). |
| Operational backout needed within 14 days of v0.2.2 | Low | High | R10: ship v0.2.3 with default reverted to `'legacy'`. Pre-committed in CHANGELOG. |

---

## 8. Success Metrics

### Acceptance (Hard Gates — all must pass)
- ✅ All 5 sequential MiniPRDs land with green audits (≤3 audit cycles per MiniPRD per P8.5; ≥3 cycles triggers Architect escalation, not silent retry).
- ✅ `verify_roundtrip.py --backend direct` passes on all current fixtures: Cisplatin, Transplatin, Cis-PtCl₂(en), Ferrocene, fac/mer-Ir(ppy)₃, PdCl₂-R-BINAP, PdCl₂-RR-BDNN, PdCl₂-RR-BDPP, TiCat1/3/4.
- ✅ `verify_roundtrip.py --backend legacy` passes on all current fixtures (blocking; novel-fixture fail-soft per D7 only applies post-v0.2.2).
- ✅ Cross-backend mean RMSD < 1.0Å on every current fixture (no per-fixture exclusions allowed for v0.2.2 ship).
- ✅ Unit-test coverage ≥80% line / ≥70% branch on `parse_oin_direct` and helpers; ≥1 test fixture per geometry code (9 codes), with at least one polydentate fixture per ligand class (en, BINAP, BDPP) per A5.
- ✅ `OIN3DGenerator(backend="direct")` and `OIN3DGenerator(backend="legacy")` both work on the full fixture suite.
- ✅ `from oinsmiles import parse_oin_direct` succeeds; module attribute `__experimental__ == True`.
- ✅ `import oinsmiles` succeeds when `scine_molassembler` import fails (R8 import-resilience test passes).
- ✅ `OIN3DGenerator(random_seed=42)` produces reproducible XYZ across two consecutive calls (FP-tolerance equality).
- ✅ Submitting an OIN string > `MAX_OIN_STRING_LEN` raises `OINFormatError` *before* regex runs (R3 DoS test).
- ✅ Submitting a v2.4 sidecar string raises `OINFormatError` with `backend='legacy'` recommendation (D6, US-008).
- ✅ Regression test `test_no_bond_to_self_after_eta_translation` passes (A9).
- ✅ CHANGELOG documents: silent default-flip in US-001, legacy-deprecation roadmap, R10 backout commitment, `_extract_oin_constraints` privatization.
- ✅ `pyproject.toml` pins `scine_molassembler >= 3.0.0, < 4.0.0`.

### Quality (Soft Goals)
- 📊 Cyclomatic complexity ≤ 10 per function (via `radon cc`); no function > 50 lines (replaces "<200 lines" per A6).
- 📊 No new lint or type-check failures introduced (`mypy --python-version 3.10 3.12`).
- 📊 No new `# type: ignore` comments without an inline justification.

### v0.3.0 Roadmap Inputs
- 📈 If direct backend handles ≥95% of fixtures with mean RMSD < 0.5Å, plan legacy removal for v0.3.0.
- 📈 If polydentate handling has ≤2 entries in `_exclusions.yml`, declare architecture sound for v0.3.0 default.
- 📈 Document any geometry code where direct < legacy parity for v0.3.0 follow-up.
- 📈 v0.2.3 ships with `FutureWarning` on `backend='legacy'` selection (D7).

### Operational Tracking
- 🛡️ R10: monitor v0.2.2 release issues for 14 days; default-revert in v0.2.3 if regressions surface.
- 🛡️ A11: ship confidence = `min(MiniPRD confidences)`; track per-MiniPRD audit cycle count.
- 🛡️ A10: any MiniPRD started >14 days post-audit triggers re-validation of relevant audit section.

---

## Appendix A: Audit Cross-Reference

This SuperPRD operationalizes recommendations in `spec/audit/DirectParser_IntegrationAudit_20260506.md`:
- Audit "Blockers 1–5" → MiniPRDs in Execution Checklist (§5).
- Audit "Option C" → This entire SuperPRD.
- Audit "Recommended Path" → Sequential execution policy (Q2=B).
- Audit "Example 1: Mapping Problem Illustrated" → Primary technical reference for MiniPRD #1.
- Audit "Example 3: Legacy Adapter Comparison" → Primary technical reference for MiniPRD #2.

## Appendix B: Red Team Coverage Map

Every Red Team finding (`spec/active/RedTeam_Report.md`, archived alongside this SuperPRD) is addressed:
- **P0 items (5):** D1, D2, D3, D4, D5 — all resolved.
- **P1 items (7):** D6, R3, A1, R10, A10, P7.4 (execution order), missing risk entries Q7.3 (a–d) — all addressed.
- **P2 items:** D6 (`@experimental`), R4/R6/R8 (NFRs), A6 (cyclomatic complexity), structured-logging deferred to v0.2.3 (out of scope).
- **P3 items:** R7 (random_seed) accepted; LRU cache deferred to v0.2.3+ (out of scope).

---

**End of SuperPRD.**

**Implementation Order:** MiniPRD #1 (FragmentMapping) → #2 (Permutation) → #3 (EtaBonds) → #4 (Polydentate) → #5 (Tests).
