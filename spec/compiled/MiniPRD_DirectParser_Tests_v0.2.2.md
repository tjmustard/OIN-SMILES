# MiniPRD: Direct Parser — Test Suite, Backend Wiring, and Public-API (v0.2.2 Blocker #5)

**Hypergraph Node ID:** `module_direct_parser_tests_v022`
**Parent Node:** `system_direct_parser_v022`
**Parent SuperPRD:** `spec/compiled/SuperPRD_DirectParser_v0.2.2.md`
**Audit Reference:** `spec/audit/DirectParser_IntegrationAudit_20260506.md` — "Blocker #5: Test Coverage Missing"
**Execution Order:** 5 of 5 (final integration + ship preparation)
**Blocked By:** MiniPRD #1, #2, #3, #4 (all)
**Estimated Effort:** 4–6 hours
**Priority:** High (closes ship-readiness)

---

## 1. The Confidence Mandate

**Confidence Score:** 9/10

**Rationale:**
- ✅ Each prior MiniPRD lands its own focused tests; this MiniPRD is composition + cross-cutting concerns + ship-prep.
- ✅ `verify_roundtrip.py` is established infrastructure; `--backend` flag is a small extension.
- ✅ Coverage tooling (`coverage` package, `pytest`/`unittest`) is standard.
- ⚠️ Cross-backend RMSD comparison across all 11 fixtures will surface any latent divergences not caught by per-MiniPRD tests — surprises are possible here.

---

## 2. Atomic User Stories

- **US-001:** As `OIN3DGenerator`, I want a `backend` constructor parameter and a per-call `backend=` override on `generate()`, so callers can choose the pipeline (matches SuperPRD US-001).
- **US-002:** As `OIN3DGenerator`, I want a `random_seed` constructor parameter that produces reproducible XYZ when set, so determinism is testable (R7).
- **US-003:** As `parse_oin_direct()`, I want `MAX_OIN_STRING_LEN` enforcement and v3.6-only format detection, so DoS attempts and v2.4-input misuse fail-loud at the boundary (R3, D6, US-008).
- **US-004:** As `oinsmiles` package user, I want `parse_oin_direct` re-exported with `@experimental` markings (US-007).
- **US-005:** As `oinsmiles` package, I want to remain importable when `scine_molassembler` import fails, so XYZ→OIN-only users on unsupported platforms aren't blocked (R8, US-009).
- **US-006:** As CI infrastructure, I want both backends tested per the D7 policy (legacy blocking on existing fixtures; fail-soft on novel post-v0.2.2 fixtures only).
- **US-007:** As a release manager, I want CHANGELOG, version bump, and dependency pin all landed in this MiniPRD so v0.2.2 is shippable on close.

---

## 3. Implementation Plan (Task List)

### Backend wiring (US-001, US-002)
- [ ] **Task 1:** Add `backend: Literal["direct","legacy"] = "direct"` param to `OIN3DGenerator.__init__`. Validate; store as `self._backend`.
- [ ] **Task 2:** Add `random_seed: Optional[int] = None` to `OIN3DGenerator.__init__`. Store as `self._random_seed`. Pass through to DG worker.
- [ ] **Task 3:** Add `backend: Literal["direct","legacy"] | None = None` to `OIN3DGenerator.generate()`. If `None`, use `self._backend`. Document precedence in docstring (A1).
- [ ] **Task 4:** Wire `self._backend` to dispatch: `"direct"` → `parse_oin_direct(...)`; `"legacy"` → `OINParser.parse(...) + MolassemblerAdapter.generate(...)`. Both return `GeneratedStructure`.
- [ ] **Task 5:** Add `ValueError("backend must be 'direct' or 'legacy'")` for invalid values (constructor and per-call).

### Boundary safety (US-003, R3, D6)
- [ ] **Task 6:** Add `MAX_OIN_STRING_LEN = 100_000` constant in `engine.py`. At the top of `parse_oin_direct`, raise `OINFormatError(f"OIN string too long: {len(s)} > {MAX_OIN_STRING_LEN}")` before any regex.
- [ ] **Task 7:** Add v3.6 format detection at the top of `parse_oin_direct`. If the regex doesn't match v3.6 inline pattern, raise `OINFormatError("v3.6 inline only — use OIN3DGenerator(backend='legacy') for v2.4")` with input truncated to 200 chars.

### Public API (US-004, R8, US-009)
- [ ] **Task 8:** Re-export `parse_oin_direct` from `src/oinsmiles/__init__.py`. Add to `__all__` (additions only per P6.3). Set module attribute `__experimental__ = True` in `engine.py`.
- [ ] **Task 9:** Update `parse_oin_direct` docstring with the `@experimental` block (per A3 spec in SuperPRD).
- [ ] **Task 10:** Defer `import scine_molassembler` from module-top to inside `OIN3DGenerator.__init__` (or first use). Wrap the import in try/except that re-raises `ImportError("scine_molassembler is required for OIN→XYZ; install via 'uv add scine_molassembler' or 'pip install scine_molassembler')` with a remediation hint.

### Test infrastructure (US-006, broad coverage)
- [ ] **Task 11:** Add `--backend {direct,legacy}` CLI flag to `tests/integration/verify_roundtrip.py`. Default = `direct` to mirror runtime default.
- [ ] **Task 12:** Update CI config: run `verify_roundtrip.py --backend direct` (blocking on all fixtures) and `--backend legacy` (blocking on existing fixtures; fail-soft on novel fixtures introduced post-v0.2.2). Implement via fixture tags (`existing` vs `novel`) in fixture JSON metadata.
- [ ] **Task 13:** Create `tests/fixtures/oin/_exclusions.yml` (initially empty; schema documented in inline comments). Add a `tests/unit/test_exclusions_loader.py` smoke test that asserts the file is parseable and that any entry has the required fields (`fixture_id`, `backend`, `justification`, `tracking_issue`).

### Unit-test coverage (US-006, A5)
- [ ] **Task 14:** Create `tests/unit/test_parse_oin_direct.py`. Add ≥1 test per `SUPPORTED_GEOMETRY_CODES` entry (9 codes total). Each test parses an OIN string for a representative complex of that geometry; asserts `GeneratedStructure` returned with non-empty `xyz`.
- [ ] **Task 15:** Add cross-backend RMSD tests: for each of the 11 current fixtures, run both backends; assert mean RMSD < 1.0Å. (This is the core D2 / Q8 acceptance gate.)
- [ ] **Task 16:** Add timeout tests: `parse_oin_direct(s, timeout=1)` on a fixture known to take >1s raises `TimeoutError`.
- [ ] **Task 17:** Add R10 backout-readiness smoke: a script `tools/revert_default_to_legacy.py` (one-line patch generator for v0.2.3) lives in tools/. Test asserts the script's output applies cleanly to `engine.py`.

### Determinism + import resilience (US-002, US-005)
- [ ] **Task 18:** `tests/unit/test_determinism.py::test_random_seed_reproducible`. Two `OIN3DGenerator(random_seed=42).generate(s)` calls produce identical XYZ (FP-tolerance via `numpy.allclose`).
- [ ] **Task 19:** `tests/unit/test_determinism.py::test_random_seed_none_stochastic`. Two `OIN3DGenerator(random_seed=None).generate(s)` calls produce different XYZ (stochastic ensemble preserved).
- [ ] **Task 20:** `tests/unit/test_import_resilience.py::test_import_without_molassembler`. Patch `sys.modules["scine_molassembler"] = None` (or use `importlib` magic); assert `import oinsmiles` succeeds and `XYZToSMILES` is callable; assert `OIN3DGenerator()` raises `ImportError` with the remediation message.

### Coverage gate (A5)
- [ ] **Task 21:** Run `coverage run -m unittest discover tests && coverage report --include='src/oinsmiles/generation/*'`. Assert ≥80% line, ≥70% branch on `engine.py`, `permutation.py`, `polydentate.py` (whichever exist). Add a CI step that fails if coverage drops.

### Quality gate (A6)
- [ ] **Task 22:** Run `radon cc src/oinsmiles/generation/engine.py src/oinsmiles/generation/permutation.py src/oinsmiles/generation/polydentate.py --min=B`. Assert no function exceeds cyclomatic complexity 10 or 50 lines. (Add `radon` to dev-dependencies in `pyproject.toml`.)

### Ship prep (US-007)
- [ ] **Task 23:** Bump version in `pyproject.toml` from `0.2.1` → `0.2.2`.
- [ ] **Task 24:** Update `pyproject.toml` dependency: `scine_molassembler` from `>=2.0.0` to `>=3.0.0,<4.0.0` (D5).
- [ ] **Task 25:** Update `CHANGELOG.md` with v0.2.2 section. Required entries:
  - **Default backend flip:** `OIN3DGenerator()` now defaults to `backend='direct'` (was implicitly legacy in v0.2.1). Existing callers may opt back into legacy via `OIN3DGenerator(backend='legacy')`.
  - **Operational backout commitment (R10):** if `direct` regresses production within 14 days, v0.2.3 will revert default to `legacy`.
  - **Legacy deprecation roadmap (D7):** `FutureWarning` on `backend='legacy'` selection planned for v0.2.3; legacy backend removed in v0.3.0.
  - **Dependency pin:** `scine_molassembler >=3.0.0,<4.0.0`.
  - **Privatization:** `extract_oin_constraints` → `_extract_oin_constraints` (was undocumented helper; not a public-API breaking change).
  - **New public API:** `from oinsmiles import parse_oin_direct` (`@experimental`).
  - **New constructor params:** `OIN3DGenerator(random_seed=...)`, `OIN3DGenerator(backend=...)`.
- [ ] **Task 26:** Run final acceptance suite: `uv sync && uv run python -m unittest discover tests && uv run python tests/integration/verify_roundtrip.py --backend direct && uv run python tests/integration/verify_roundtrip.py --backend legacy`. All must pass.

---

## 4. The Negative Space (Constraints)

- **DO NOT** modify `MolassemblerAdapter`, `OINParser.parse()`, `xyz2mol.py`, or `core/chirality.py` — all out of perimeter (P6.1, P6.4).
- **DO NOT** change `OIN3DGenerator.timeout` default of 60 seconds (P6.2).
- **DO NOT** remove existing entries from `oinsmiles.__all__` (P6.3).
- **DO NOT** ship if `verify_roundtrip.py --backend direct` fails on any current fixture.
- **DO NOT** ship if cross-backend mean RMSD ≥ 1.0Å on any current fixture (no exclusions for v0.2.2 ship per D7).
- **DO NOT** add structured logging or LRU cache in this MiniPRD (deferred to v0.2.3).
- **DO NOT** mark the legacy backend as `@deprecated` in v0.2.2 — `FutureWarning` starts in v0.2.3 only (D7).
- **DO NOT** suppress test failures with `@unittest.skip` to make CI green; instead, route to `_exclusions.yml` with justification.
- **DO NOT** change `GeneratedStructure` dataclass; *adding* optional fields with defaults is permitted (P6.6) but not required by this MiniPRD.

---

## 5. Integration Tests & Verification

### Test 1 (Deterministic — full fixture suite, both backends):
- **Input:** All 11 current fixtures (Cisplatin, Transplatin, Cis-PtCl₂(en), Ferrocene, fac/mer-Ir(ppy)₃, PdCl₂-R-BINAP, PdCl₂-RR-BDNN, PdCl₂-RR-BDPP, TiCat1/3/4).
- **Expected Output:** `verify_roundtrip.py --backend direct` passes all (mean RMSD < 1.0Å each). `verify_roundtrip.py --backend legacy` passes all (existing-fixture blocking per D7).

### Test 2 (Deterministic — cross-backend parity):
- **Input:** Each of the 11 fixtures, parsed by both backends.
- **Expected Output:** Pairwise mean RMSD between direct-XYZ and legacy-XYZ < 1.0Å.

### Test 3 (Deterministic — public API):
- **Input:** `import oinsmiles; oinsmiles.parse_oin_direct`.
- **Expected Output:** Resolves to the function. `oinsmiles.parse_oin_direct.__doc__` contains "Experimental". `oinsmiles.engine.__experimental__ is True`.

### Test 4 (Deterministic — DoS protection):
- **Input:** OIN string of length `MAX_OIN_STRING_LEN + 1`.
- **Expected Output:** `OINFormatError` raised before regex runs (timing assertion: < 10ms).

### Test 5 (Deterministic — v2.4 rejection):
- **Input:** A v2.4 sidecar OIN string.
- **Expected Output:** `OINFormatError` containing the substring `"backend='legacy'"`.

### Test 6 (Deterministic — random seed reproducibility):
- **Input:** Two `OIN3DGenerator(random_seed=42).generate(cisplatin_oin)` calls.
- **Expected Output:** Identical XYZ within FP tolerance.

### Test 7 (Deterministic — import resilience):
- **Input:** `import oinsmiles` with `scine_molassembler` removed from `sys.modules`.
- **Expected Output:** Import succeeds. `XYZToSMILES` callable. `OIN3DGenerator()` raises `ImportError` with remediation message.

### Test 8 (Deterministic — coverage gate):
- **Input:** `coverage run -m unittest discover tests && coverage report`.
- **Expected Output:** ≥80% line, ≥70% branch on direct-parser modules.

### Test 9 (Deterministic — cyclomatic complexity gate):
- **Input:** `radon cc src/oinsmiles/generation/{engine,permutation,polydentate}.py`.
- **Expected Output:** No function rated worse than B (≤10).

### Test 10 (Deterministic — exclusions schema):
- **Input:** Empty `tests/fixtures/oin/_exclusions.yml`.
- **Expected Output:** `test_exclusions_loader` passes; any future entry must have `fixture_id`, `backend`, `justification`, `tracking_issue`.

### Candidate Artifact Routing:
- Cross-backend RMSD divergence between 1.0Å and 1.5Å on a current fixture: **hard block for v0.2.2 ship**. (Existing-fixture blocking per D7.)
- Cyclomatic complexity violation on a single function: refactor before merge; do not exclude.
- Coverage shortfall (e.g., 78% line): add tests until threshold met; do not lower threshold.
- A new fixture introduced *after* v0.2.2 release that fails on legacy-only is the *first* legitimate `_exclusions.yml` candidate.
