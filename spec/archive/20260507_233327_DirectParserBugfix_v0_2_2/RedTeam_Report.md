# RedTeam Report: Direct Parser Bugfix & Coexistence Backend (v0.2.2)

**Subject:** `spec/active/Draft_PRD.md` (v0.0.1-draft, dated 2026-05-07)
**Reviewer:** Red Team Agent
**Date of Analysis:** 2026-05-07
**Posture:** Adversarial / hostile-but-constructive
**Cross-Reference:** `spec/compiled/architecture.yml` (last-audited 2026-05-06)

---

## Executive Summary

The Draft PRD is well-structured and operationalizes the audit findings clearly. However, it carries several latent risks that weaken its claim of being a controlled v0.2.2 release:

1. **Coexistence is asymmetric**: direct = blocking, legacy = fail-soft, but legacy is the *current production*. This inverts the usual "trust new only after burn-in" pattern.
2. **Silent default flip** in US-001 (no `backend` kwarg → direct) is a backwards-incompatible behavioral change disguised as additive.
3. **Sub-process semantics, timeouts, memory, and concurrency** are unspecified — the diagram mentions a 10s timeout that contradicts the 60s default in `__init__`.
4. **Public API stability claim for `parse_oin_direct`** (Q6=A) conflicts with the fact that it is being rewritten in v0.2.2 itself.
5. **Test parity is undefined**: cross-backend RMSD < 1.0Å is loose enough that geometrically wrong outputs may pass.
6. **Three open questions are deferred to Red Team** without proposing a specific resolution path or fallback if the answer is unfavorable.

This document attacks each section and proposes hardening actions.

---

## Section 1: Introduction & Goals — Analysis

### Clarifying Questions
- **Q1.1:** The Problem Statement enumerates *five* blockers and frames them as exhaustive. What evidence justifies "five and only five"? Was a fault-injection or adversarial walk of `parse_oin_direct()` performed, or only a static read of audit findings? If only static, what is the discovery rate of new blockers we should expect?
- **Q1.2:** Solution Overview states "Both backends must pass `verify_roundtrip.py` test parity" — yet Q7=C makes legacy fail-soft. Which is authoritative? If "must pass" is rhetorical, restate it; if it is binding, drop the fail-soft policy.
- **Q1.3:** Tertiary audience names "future contributors extending OIN→XYZ to new ligand types" but Section 6 forbids any modification to the legacy adapter and Q9 explicitly excludes new geometries from v0.2.2. How does v0.2.2 produce extension benefits the Tertiary audience can use?
- **Q1.4:** The phrase "graceful drift in v0.3.0+" appears nowhere else. What does "graceful drift" mean operationally — accepting parity regressions? Tolerating diverging XYZ output? Define.

### What-If Scenarios
- **W1.1 — Sixth blocker discovered post-merge:** Direct backend ships, a user submits a fixture (e.g., μ-bridged dimer or zero-valent carbonyl) that exposes a sixth blocker. Because direct is the default, the user gets a hard error. What is the rollback procedure? Is there a kill-switch env var? Currently there is none specified.
- **W1.2 — Audit document corrupted/contradicts code:** The audit at `spec/audit/DirectParser_IntegrationAudit_20260506.md` is treated as authoritative. If the audit's interpretation of legacy `_pick_masm_permutation()` is wrong, all five MiniPRDs inherit the error. There is no specified counter-check (e.g., diff against running code, independent re-derivation).
- **W1.3 — Legacy backend rots silently:** Q7=C makes legacy fail-soft. Contributors stop maintaining it. By v0.3.0 deletion target, legacy parity is unknown, and no migration data exists for the v0.2.0 → v0.3.0 transition.
- **W1.4 — User reads "v3.6 inline only" too late:** A user with v2.4 sidecar OIN strings hits direct backend, gets a cryptic regex-failure error, switches to legacy, and works around it. The bug report sits unanswered; the team thinks adoption is fine.

### Points for Improvement
- **P1.1:** Rephrase "five blockers" as "five known blockers per audit; additional may emerge during integration" and add an NFR: any newly discovered blocker post-MiniPRD-#1 triggers an explicit Architect-level review before ship.
- **P1.2:** Remove the contradiction between Solution Overview ("must pass") and Q7 ("legacy fail-soft"). Recommend: legacy fail-soft applies only to *novel* fixtures added in v0.3.0; existing fixtures must remain blocking on both backends.
- **P1.3:** Add an explicit "Operational Backout" subsection: "If `backend='direct'` regresses production within 14 days of v0.2.2 release, change default back to `'legacy'` in a v0.2.3 patch." Make this a written commitment, not implicit.
- **P1.4:** Define "graceful drift" or remove the phrase — undefined hedges are technical debt magnets.

---

## Section 2: Confidence Mandate — Analysis

### Clarifying Questions
- **Q2.1:** The Confidence Score is 8/10, but three Q-OPEN questions remain. The Confidence Mandate normally requires score ≥9 for execution. What gates the move from 8→9? Q-OPEN-1, Q-OPEN-2, Q-OPEN-3 all need answers — does each unblock 0.33 confidence points? Make the gating explicit.
- **Q2.2:** The bullet "Polydentate ligand mapping (Blocker 2) requires understanding of vector data structure — may need design iteration" — what is the trigger for declaring "design iteration was insufficient" and escalating? Time-box, fixture-failure-count, or reviewer judgment?
- **Q2.3:** "Permutation selection (Blocker 3) is the most complex blocker" — has anyone read `_pick_masm_permutation()` end-to-end before signing this PRD? If not, the 8/10 is overstated.

### What-If Scenarios
- **W2.1 — Q-OPEN answers contradict architecture:** Suppose Q-OPEN-1's answer is "no, fragments[0] is not always the metal — depends on canonical SMILES order." That invalidates audit examples relying on `fragments[0]` and re-opens MiniPRD #1's design. The Draft PRD has no contingency plan for "Q-OPEN answer reshapes scope."
- **W2.2 — Confidence inflation cascade:** MiniPRD #1 lands at confidence 7. MiniPRD #2 inherits #1's foundation and is graded 8 because "#1 already proved it." MiniPRD #3 is 8.5 by similar logic. By the end, the *aggregated* confidence of the v0.2.2 ship is below 7, but no MiniPRD records it that way.
- **W2.3 — Audit document age:** The audit is dated 2026-05-06 and the PRD 2026-05-07. The codebase may shift during the 3-4 day execution window. There is no policy that the audit be re-validated after each MiniPRD lands.

### Points for Improvement
- **P2.1:** Tie confidence to a checklist: each Q-OPEN must be resolved by a named owner with a written answer in the PRD before any MiniPRD is unblocked. Currently they are "for Red Team" — not for *resolution*.
- **P2.2:** Add an aggregated-confidence rule: ship confidence = min(MiniPRD confidences). If any drops below 7, all subsequent MiniPRDs require re-audit before merge.
- **P2.3:** Mandate a pre-MiniPRD-#2 read-through and write-up of `_pick_masm_permutation()` as a gating deliverable, captured in the MiniPRD itself. Without it, MiniPRD #2 cannot start.
- **P2.4:** Establish an audit-staleness policy: if the audit is more than N days old at the start of any MiniPRD, re-validate the relevant section before proceeding (where N is small, e.g., 14 days).

---

## Section 3: Scope — Analysis

### Clarifying Questions
- **Q3.1:** "Update CHANGELOG with breaking-change notice for `extract_oin_constraints()` signature" — but the function is not in the public `__all__`. Is it considered public-API by virtue of being importable, or private-with-underscore-equivalent? If private, why a breaking-change notice at all?
- **Q3.2:** Out-of-scope says "Bit-identical XYZ output between backends (Q8: 1.0Å tolerance accepted)." 1.0Å is *huge* — typical bond lengths are 1.0–2.5Å, so atoms could be off by ~50% of a bond length and still pass. What's the justification for 1.0Å rather than, say, 0.1Å?
- **Q3.3:** "Migration of `OINInlineHandler`" is out of scope. But `OINInlineHandler` is the only path for v3.4+ headings/winding metadata into the OIN string. If direct parser doesn't honor those tokens, is that a parser bug or a scope exclusion?
- **Q3.4:** Where in the scope is the policy for environments where `scine_molassembler` is unavailable (e.g., Windows, Apple Silicon edge cases)? Does the direct backend gracefully degrade or hard-fail at import time?

### What-If Scenarios
- **W3.1 — Tolerance gaming:** A direct-backend output is 0.95Å off legacy on Cisplatin (under threshold but visibly wrong). It passes CI. No one notices until a downstream user reports DFT energies diverging by 5 kcal/mol because Pt-Cl bonds are systematically stretched.
- **W3.2 — v2.4 input passed to direct:** A user pipes a v2.4-sidecar OIN string into `OIN3DGenerator(backend="direct")`. The regex preprocessor fails silently or produces a malformed mapping. Out-of-scope says "no support" but Negative Constraints don't enforce a clear error.
- **W3.3 — `extract_oin_constraints()` drift:** External consumer (e.g., a chemistry plugin not in the repo) imports `extract_oin_constraints` from a `_internal` module. Breaking the signature breaks them silently. We don't know who they are.
- **W3.4 — Performance regression:** Adding fragment→atom mapping doubles parse time. Out of scope says "no perf optimization beyond audit-identified blocker." A 10× slowdown ships and gets blamed on Molassembler.

### Points for Improvement
- **P3.1:** Either declare `extract_oin_constraints` private (rename to `_extract_oin_constraints`) and remove the "breaking change" framing, or commit to maintaining the new 3-tuple signature for at least one major version. Pick one.
- **P3.2:** Tighten the cross-backend tolerance: 1.0Å for full-structure RMSD is too loose. Recommend per-atom max-displacement < 0.5Å OR full-structure RMSD < 0.3Å, whichever is stricter. If 1.0Å is genuinely the right number, justify with a citation to existing fixture variance.
- **P3.3:** Add an explicit fail-loud policy for v2.4 sidecar input: "`parse_oin_direct()` raises `OINFormatError` if input is not v3.6 inline; error message points users to `OIN3DGenerator(backend='legacy')`."
- **P3.4:** Add a benchmarking checkpoint: pre-v0.2.2-tag, run `verify_roundtrip.py` on both backends and capture timing. If direct is >2× slower than legacy, flag for v0.2.3 perf work. Currently there's no perf telemetry at all.
- **P3.5:** Add a scope item: emit a `DeprecationWarning` if `extract_oin_constraints` is imported from outside the package, OR underscore-prefix it. Don't have a "breaking change" for an undeclared-private function with no migration story.

---

## Section 4: User Stories — Analysis

### Clarifying Questions
- **Q4.1 (US-001):** Acceptance criterion 4 says "Existing callers without backend kwarg get direct pipeline." This is a *silent* default change. v0.2.1 callers will get different behavior on a patch-version bump from 0.2.1 → 0.2.2. Does our semver allow this? Justify or change to: "default = legacy in v0.2.2; default = direct in v0.3.0."
- **Q4.2 (US-002):** "All 17 test call sites updated" — was this number derived from `grep extract_oin_constraints`? What about indirect callers via `partial(...)` or `getattr(module, 'extract_oin_constraints')`? Is the count audited?
- **Q4.3 (US-003):** "RMSD < 1.0Å vs. input XYZ" — input XYZ is the ground truth for round-trip testing, but for *novel* OIN strings (no input XYZ), what is the validation? The story is silent.
- **Q4.4 (US-004):** "Permutation index selection logic ported from legacy `_pick_masm_permutation()`." But the legacy code may use heuristics tied to internal state of `MolassemblerAdapter`. How do we verify the *port* is correct, not just structurally similar?
- **Q4.5 (US-005):** "No 'bond to self' errors (tracking issue from audit)" — is there a dedicated regression test that asserts this specific error never re-occurs? If not, add one.
- **Q4.6 (US-006):** "verify_roundtrip.py --backend legacy (fail-soft / warning)" — what is the *escalation procedure*? If legacy fails 3 fixtures in a row across 3 PRs, do we block? Define a numerical threshold.
- **Q4.7 (US-007):** "Docstring documents it as low-level direct parser" — but the implementation is being rewritten *in this same release*. How is "low-level utility" stable if the implementation is volatile?

### What-If Scenarios
- **W4.1 — Default-change ambush:** A downstream user pins `oinsmiles==0.2.*`, expects patch-level stability, and gets a backend swap on `pip install -U`. A production pipeline that worked yesterday produces different XYZ today. They blame *us*, not the version pin.
- **W4.2 — Test-site count drift:** During refactor, MiniPRD #1 lands updates to "all 17 call sites." Between MiniPRD #1 and #5, a new test fixture adds an 18th call site. The 18th breaks but only on `--backend direct` — silent failure if test isn't tagged.
- **W4.3 — Polydentate edge case (US-003):** PdCl2-RR-BDPP has a P-chiral center with 3D constraints. Direct backend produces a valid structure with *opposite* chirality (R becomes S). RMSD passes the 1.0Å threshold by accident due to symmetry. Test is green, science is wrong.
- **W4.4 — Legacy-only fixture (US-006):** A complex like `TiCat4` works on legacy but fails on direct. CI marks direct as failing → release blocked. The fix is to add a per-fixture exclusion. Now the exclusion list grows over time and becomes a hidden catalog of "things direct can't do" with no SLA.
- **W4.5 — Public API surface lock-in (US-007):** v0.2.2 ships `from oinsmiles import parse_oin_direct`. v0.3.0 needs to refactor it. We have to deprecate-and-replace, doubling the surface area.

### Points for Improvement
- **P4.1:** Change US-001 default to `backend="legacy"` for v0.2.2; flip default to `"direct"` in v0.3.0 after one full release cycle of soak time. Document the migration in `CHANGELOG.md` and add a single-emission `FutureWarning` when default is implicitly used.
- **P4.2:** Add an automated audit script: `tools/audit_extract_calls.py` that greps the codebase and tests, lists all call sites, and produces a count. Re-run as a CI step. Static "we know there are 17" is brittle.
- **P4.3:** For US-004, mandate a *unit* test that constructs a known-symmetric input (e.g., regular octahedron of identical Cl atoms) and verifies the permutation index *value* matches the legacy adapter's choice — not just the resulting RMSD.
- **P4.4:** For US-005, add a regression test named `test_no_bond_to_self_after_eta_translation` that pins down the exact failure mode from the audit.
- **P4.5:** For US-006, add an explicit escalation rule: "If legacy backend fails any blocking fixture for 2 consecutive merges, legacy gets demoted to 'archived' status and a tracking issue opened to remove it."
- **P4.6:** For US-007, mark `parse_oin_direct` as `@experimental` (Python decorator + docstring tag) for v0.2.2; promote to "stable" only in v0.3.0 once shape is settled.

---

## Section 5: Technical Specifications — Analysis

### Clarifying Questions
- **Q5.1 (Pipeline diagram):** The diagram annotates `_dg_worker() (timeout 10s)` but the API contract specifies `timeout: int = 60` on `OIN3DGenerator.__init__`. Which is authoritative? If 10s is a hardcoded sub-step within a 60s budget, document the budget partition.
- **Q5.2 (Blast Radius):** `atom_smiles_to_xyz` is currently `needs_review` per `architecture.yml`. It depends on `atom_oin3d_generator`. Why isn't it in the Modified Nodes list? A `backend` kwarg change to `OIN3DGenerator.generate()` ripples to anything that calls it.
- **Q5.3 (New Nodes):** `atom_direct_parser_polydentate` and `atom_direct_parser_permutation` are listed as new but no `associated_file` is given. Will they be new files, or new functions inside `oin_parser.py` / `engine.py`? If functions, why are they distinct nodes vs. expansions of existing atomic nodes?
- **Q5.4 (Sequential MiniPRDs):** "Total estimated effort: ~20-26 hours." This presumes zero rework. Audit-discovery rate of new blockers, integration friction, and test-suite churn are not budgeted. What contingency is built in?
- **Q5.5 (API Contract):** `def generate(self, oin_string: str, backend: Optional[Literal["direct", "legacy"]] = None)`. If `backend=None`, which is used — `self._backend` from `__init__` or the constructor default? If the latter, the Literal type lies (None is a valid value).
- **Q5.6 (Dependencies):** "scine_molassembler >= 2.0.0 (already pinned)." But CLAUDE.md memory says current version is 3.0.1. Is `>= 2.0.0` the right floor? Major-version drift between 2.x and 3.x could break the API.
- **Q5.7 (Free Function):** `parse_oin_direct(oin_smiles: str) -> GeneratedStructure` has no timeout parameter. How does it behave on infinite-loop input? Does it inherit from a global, raise, or hang?

### What-If Scenarios
- **W5.1 — ProcessPoolExecutor leak:** `parse_oin_direct` is called in a tight loop for batch processing. ProcessPoolExecutor spins up subprocesses each call. Memory grows unboundedly. No pooling strategy is specified.
- **W5.2 — Concurrent backend mutation:** A long-running service constructs `OIN3DGenerator()` once and calls `.generate(oin, backend="direct")` and `.generate(oin, backend="legacy")` from different threads. State pollution? Race conditions? The PRD doesn't mention thread safety.
- **W5.3 — Subprocess pickling failure on `parse_oin_direct`:** ProcessPoolExecutor requires picklable args. If `parse_oin_direct` closes over RDKit `Chem.Mol` objects (un-picklable), timeouts will fail at submit, not at execution. CLAUDE.md notes "confirmed picklability" for Molassembler — was the same done for the new direct path?
- **W5.4 — Geometry types not in TEMPLATES:** The Pipeline diagram step 5 says "_select_permutation()." If the OIN geometry code is not in `TEMPLATES` (e.g., PBP-7, hypothetical 8-coord), what happens? `KeyError`? `NotImplementedError`? Currently undefined.
- **W5.5 — RDKit aromatic perception drift:** The AST tokenization stage relies on RDKit's aromatic flag handling. If a future RDKit upgrade changes Kekulé perception, both backends drift but only direct breaks.
- **W5.6 — Polydentate-fragment-rank ambiguity:** A bidentate ligand (en) has two binding atoms. `vertex_indices` may list both. The fragment→atom mapping needs a 1-to-many relationship in some cases. Does `Dict[int, List[int]]` cover the inverse direction (atom → fragment) too? Probably needed for blocker #2.

### Points for Improvement
- **P5.1:** Make timeout semantics explicit. Recommend a contract: `OIN3DGenerator(timeout=N)` budgets total time; `_dg_worker` gets `min(N - parsing_time, 30)`. Hardcoded 10s is suspect.
- **P5.2:** Update the Blast Radius to include `atom_smiles_to_xyz` — even if no code change, mark `needs_review` because its dependency footprint changes.
- **P5.3:** Either materialize `atom_direct_parser_polydentate` and `atom_direct_parser_permutation` as new files (e.g., `polydentate.py`, `permutation.py`) or merge them into existing nodes. New nodes without files are tracking-debt.
- **P5.4:** Add a 25% time contingency: rebudget to ~25-32 hours. Cap any single MiniPRD at 1.5× its estimate; trigger Architect review on overrun.
- **P5.5:** Fix the `Optional[Literal[...]]` typing: either drop `None` and require explicit choice, or add `None` to the Literal and document the precedence rule clearly.
- **P5.6:** Pin `scine_molassembler == 3.0.*` (or `>= 3.0.0, < 4.0.0`) since CLAUDE.md confirms 3.0.1 is in use. The `>= 2.0.0` floor is too loose.
- **P5.7:** Specify behavior for `parse_oin_direct` with no timeout: either inherit a module-level constant, or raise `TypeError` if no timeout is in scope. Hanging is unacceptable.
- **P5.8:** Add a thread-safety NFR: declare whether `OIN3DGenerator` is thread-safe. If not, document. If yes, add a thread-stress test.
- **P5.9:** Add a memory NFR: `parse_oin_direct` must not retain references to `Chem.Mol` after returning. Add a regression test using `weakref.finalize`.
- **P5.10:** Define behavior for unsupported geometries: `raise NotImplementedError(f"Geometry {code} not supported by direct backend; try backend='legacy'")`.

---

## Section 6: Negative Constraints — Analysis

### Clarifying Questions
- **Q6.1:** "DO NOT silently fall back from direct to legacy on error" — combined with "default = direct" (US-001), how does a user discover that legacy is the right escape hatch? Is the error message specified? Required content?
- **Q6.2:** "DO NOT mark legacy backend as deprecated in v0.2.2." But Section 8 ("v0.3.0 Roadmap Inputs") explicitly plans legacy removal. What is the user-facing signal in v0.2.2 that legacy is on its way out? Without one, v0.3.0's removal is a surprise.
- **Q6.3:** "DO NOT change `GeneratedStructure` dataclass signature." Defensive — but what about adding *optional* fields with defaults (e.g., `backend_used: str = "direct"`)? Is that allowed?
- **Q6.4:** "DO NOT introduce a third backend or pluggable backend registry." But the Sequential MiniPRDs include separate atom nodes (`polydentate`, `permutation`). Are those "internal pluggable" or just structured code? Define.

### What-If Scenarios
- **W6.1 — Constraint omission:** No constraint forbids modifying `xyz2mol.py` or `core/chirality.py`. A MiniPRD author might "improve" them while in the engine.py blast radius. Boom: hidden side-effect on XYZ→OIN pipeline.
- **W6.2 — DG worker timeout regression:** No constraint preserves the existing 60s default timeout. A MiniPRD might lower it to 10s "for snappier UX" and break long-running fixtures.
- **W6.3 — Legacy import path removed:** A MiniPRD author renames `MolassemblerAdapter.generate()` → `MolassemblerAdapter.legacy_generate()` for "clarity," breaking the legacy backend even though the *file* wasn't modified.
- **W6.4 — `__all__` accidental contraction:** Adding `parse_oin_direct` to `__all__` could push the list onto a single line, accidentally dropping `XYZToSMILES` or `SMILESToXYZ`. No constraint protects existing exports.
- **W6.5 — DeprecationWarning suppressed:** The constraint forbids deprecation messaging on legacy. But what if the MiniPRD adds a `FutureWarning` on the *default* selection mechanism? Allowed?

### Points for Improvement
- **P6.1:** Add: "DO NOT modify any file outside the engine.py / oin_parser.py / __init__.py / verify_roundtrip.py / test files perimeter without explicit Architect approval."
- **P6.2:** Add: "DO NOT change the default value of `OIN3DGenerator.timeout` (must remain 60s)."
- **P6.3:** Add: "DO NOT remove existing entries from `oinsmiles.__all__`; only additions allowed."
- **P6.4:** Add: "DO NOT modify `MolassemblerAdapter` public method signatures or behavior; legacy backend wiring depends on byte-stable API."
- **P6.5:** Add: "DO NOT swallow exceptions from direct backend; all failures must propagate with a message that includes the input OIN string and a recommendation to try `backend='legacy'`."
- **P6.6:** Clarify: "Adding optional fields to `GeneratedStructure` with safe defaults is allowed; removing or repurposing existing fields is not."
- **P6.7:** Replace "DO NOT mark legacy as deprecated" with "Legacy is non-deprecated in v0.2.2; the v0.3.0 RFC must include a 1-version deprecation period before removal." This converts the constraint into a forward-looking commitment.

---

## Section 7: Risks & Mitigation — Analysis

### Clarifying Questions
- **Q7.1:** Probabilities ("Medium," "High") and Impact tiers are not calibrated. What is "Medium" — 30%? 50%? Without numerics, the risk register is decorative.
- **Q7.2:** Risk: "Direct backend produces 'valid but different' XYZ that fails fixture comparison" — mitigation is "loose threshold (1.0Å)." That's not mitigation, that's *causing* the W3.1 scenario above.
- **Q7.3:** No risk identified for: (a) Molassembler license/version changes, (b) RDKit aromatic perception drift, (c) Python 3.10 → 3.12 typing changes (`Literal` semantics), (d) ProcessPoolExecutor deprecation paths.
- **Q7.4:** "MiniPRD #2 (permutation) blocks #3 and #4" — but the Execution Checklist orders work as #1→#3→#4→#2→#5. So #2 is *fourth*, blocking nothing downstream of itself except #5. Inconsistency.
- **Q7.5:** "Audit assumed `fragments[0]` is always the metal — may not hold." Listed as Probability=Low, Impact=High. Mitigation: "raise clear error if no metal in fragment 0." But Q-OPEN-1 says the team has not yet confirmed the assumption. Probability=Low is unjustified — it's actually Probability=Unknown.

### What-If Scenarios
- **W7.1 — Combinatorial risk:** Multiple "Medium probability, Medium impact" risks compound. Probability of *at least one* ship-blocking issue is much higher than any single risk.
- **W7.2 — Mitigation by-pass:** "Per-fixture exclusion list with documented justification" as mitigation for cross-backend RMSD differences becomes a graveyard list that no one prunes. By v0.3.0, half the fixtures are excluded.
- **W7.3 — Hidden risk: subprocess crash on import:** `scine_molassembler` is a C++-extension. A bad install in CI causes the subprocess to segfault before the timeout fires. The error is "subprocess died" not "Molassembler import failed." Misdiagnosed.
- **W7.4 — Hidden risk: Windows/MSVC build:** If we ever ship to Windows, scine_molassembler may need conda-forge. Direct backend's import-time dependency could make `oinsmiles` unimportable on Windows even when only XYZToSMILES is needed.

### Points for Improvement
- **P7.1:** Quantify probabilities: Low=10%, Medium=30%, High=60%. Recompute the joint risk of "any one of the High/Medium-probability items materializes" — likely >80%, which should change posture.
- **P7.2:** Replace "loose threshold" mitigation with "investigate and document each cross-backend divergence; do not raise threshold unilaterally."
- **P7.3:** Add the four missing risk entries (Q7.3 above).
- **P7.4:** Reconcile the Execution Checklist order with the dependency claim: either "permutation must land before polydentate" (and reorder Section 5) or "polydentate doesn't depend on permutation" (and remove the risk entry).
- **P7.5:** Reclassify Q-OPEN-1 risk as Probability=Unknown; gate MiniPRD #1 on resolving the question first. "Unknown" risks are not safely ship-able.
- **P7.6:** Add a master "exclusion-list watch" risk: any time a fixture is added to a per-backend exclusion list, an issue is filed for v0.3.0.
- **P7.7:** Add an "import-time degradation" risk: if Molassembler fails to import, `oinsmiles` should still be importable for XYZ→OIN use cases. Verify with a dedicated test.

---

## Section 8: Success Metrics — Analysis

### Clarifying Questions
- **Q8.1:** "Unit test coverage ≥80%." 80% of what — lines, branches, functions? On `parse_oin_direct()` only or the whole `engine.py`? What's the baseline today?
- **Q8.2:** "Code in `parse_oin_direct()` is < 200 lines." Why 200? Is shorter better, or is there a concrete maintainability reason? Could be a vanity metric.
- **Q8.3:** "If direct backend handles ≥95% of fixtures with RMSD < 0.5Å, plan legacy removal for v0.3.0." But the acceptance bar is RMSD < 1.0Å. Two thresholds — one looser (acceptance) and one stricter (roadmap). Why?
- **Q8.4:** No SLO for: (a) build time, (b) import time, (c) memory peak, (d) per-call latency. Are these explicitly out of scope for v0.2.2, or just forgotten?
- **Q8.5:** "Doctests in module pass" — does `parse_oin_direct` even have doctests? If yes, where? If not, why is this a metric?

### What-If Scenarios
- **W8.1 — Coverage-gaming:** Author writes 8 unit tests that hit easy paths (e.g., `XYZToSMILES.convert("Pt")`-like trivial cases). Coverage hits 82%. Real bugs in polydentate handling untested. Metric green; product red.
- **W8.2 — Line count gaming:** "<200 lines" achieved by extracting helpers into a hidden `_internal.py` module that adds 400 lines elsewhere. Public function looks clean; total system is more complex.
- **W8.3 — Hidden audit cycle:** "Each MiniPRD's audit completes within ≤2 cycles (no major rework)." If MiniPRD #1 takes 4 audit cycles, do we ship anyway? Soft goal vs. hard gate is unclear.
- **W8.4 — Roadmap creep:** If direct hits 94% (just below 95%), do we delay legacy removal or push the threshold? Soft thresholds invite negotiation.

### Points for Improvement
- **P8.1:** Replace "≥80% coverage" with: "≥80% line coverage on `parse_oin_direct` and helpers, ≥70% branch coverage, with explicit coverage of all 9 geometry codes including at least one polydentate fixture per ligand class (en, BINAP, BDPP)." Make the metric specific and gameable-resistant.
- **P8.2:** Drop the "<200 lines" metric or replace with a maintainability proxy (cyclomatic complexity ≤ 10 per function, no function > 50 lines).
- **P8.3:** Reconcile the two RMSD thresholds: pick one acceptance bar and stick to it. Recommend: acceptance = < 0.5Å mean, < 1.0Å max-per-atom. Roadmap promotion to v0.3.0 default uses the same metric.
- **P8.4:** Add explicit performance SLOs: import time < 2s on a clean venv; per-call latency on Cisplatin < 5s p95; memory delta < 100MB per call. Capture before/after numbers in the release notes.
- **P8.5:** Make audit-cycle limit a hard gate: ≥3 audit cycles on any MiniPRD triggers Architect-level escalation, not silent retry.
- **P8.6:** Add a regression metric: "Number of fixtures with direct-backend XYZ that diverges from legacy by > 0.3Å (hard cap: 0; tracked per release)."
- **P8.7:** Drop "doctests in module pass" unless doctests exist and are non-trivial.

---

## Section 9: Open Questions (Appendix B) — Analysis

The Draft PRD defers seven open questions to Red Team. Each gets a hostile reading.

### Q-OPEN-1: Metal detection assumption (`fragments[0]` always metal)
- **Adversarial position:** This is a *load-bearing* assumption used throughout the audit. Without empirical confirmation across all fixtures, MiniPRD #1 cannot start safely.
- **Mitigation required:** Write a small standalone script that runs the regex preprocessor on every fixture and asserts `fragments[0]` is the metal. Run before MiniPRD #1 begins. If it fails on any fixture, redesign MiniPRD #1.
- **Severity:** P0 (blocking).

### Q-OPEN-2: v2.4 sidecar format support
- **Adversarial position:** "Assume v3.6 only" is a silent-failure setup. Users with v2.4 strings will get cryptic errors.
- **Mitigation required:** Explicit format detection at the top of `parse_oin_direct`: if input doesn't match v3.6 pattern, raise `OINFormatError(f"Direct parser supports v3.6 inline only; received {detected_version}. Use OIN3DGenerator(backend='legacy') for v2.4.")`. Document.
- **Severity:** P1.

### Q-OPEN-3: Direct fails where legacy succeeds — P0 or scope-limit?
- **Adversarial position:** The PRD says "today's answer: P0; needs adversarial review." Adversarial answer: it depends on whether the failure is a *known* gap (per per-fixture exclusion list) or *unknown* (regression). Knowns can be deferred; unknowns must block.
- **Mitigation required:** Define a triage matrix:
  - Known gap (in exclusion list): user gets clear error, falls back to legacy manually. P2.
  - Unknown failure: blocks v0.2.2 release until either fixed or added to exclusion list with justification. P0.
- **Severity:** Requires both decisions; P0 + P2.

### Q-OPEN-4: Coexistence trap
- **Adversarial position:** Permanent coexistence creates a maintenance pact that history shows is rarely honored. Without a deletion commitment, legacy will haunt v0.3.0 and v0.4.0.
- **Mitigation required:** Pre-commit to legacy deletion in v0.3.0 with a documented date (e.g., 2026-09-01). Add an `__init__.py` `FutureWarning` triggered when `backend="legacy"` is selected, starting v0.2.3.
- **Severity:** Architectural — affects Q1 resolution.

### Q-OPEN-5: Legacy fail-soft as admission of expected breakage
- **Adversarial position:** Q7=C is rationalized as pragmatic, but it's also a *signal* that the team expects legacy to break. If we expect it to break, it shouldn't be the production fallback.
- **Mitigation required:** Either keep legacy blocking on existing fixtures (and fail-soft only on novel future fixtures), OR commit to a v0.3.0 removal date. Don't have it both ways.
- **Severity:** P1.

### Q-OPEN-6: Permutation portability
- **Adversarial position:** The Confidence Mandate rates this 8/10 without anyone having read `_pick_masm_permutation()` end-to-end. Confidence is overstated.
- **Mitigation required:** Pre-MiniPRD-#2 deliverable: a 1-page write-up of the legacy permutation logic and its data dependencies. If it relies on internal `MolassemblerAdapter` state, MiniPRD #2 redesigns; if pure function, ports cleanly.
- **Severity:** P0 (gates MiniPRD #2).

### Q-OPEN-7: Public API stability vs. implementation volatility
- **Adversarial position:** Re-exporting a function that's being rewritten is contradictory. The signature may be stable but error semantics, performance, and edge-case behavior are not.
- **Mitigation required:** Mark `parse_oin_direct` as `@experimental` for v0.2.2; lock signature only. Promote to "stable" in v0.3.0 with explicit semver guarantee.
- **Severity:** P2.

---

## Section 10: Cross-Cutting / "Unknown Unknowns"

Issues not tied to any single Draft PRD section.

### Missing NFRs
- **Rate limiting:** None specified. If a service wraps `OIN3DGenerator`, what's the per-second call ceiling? `ProcessPoolExecutor` will exhaust fork resources without throttling.
- **TTL / Cache semantics:** None specified. `parse_oin_direct` repeatedly called with the same input recomputes from scratch. Should there be an LRU cache?
- **Telemetry / Logging:** No structured logging. When a user reports "direct backend produced wrong geometry," what's the diagnostic info we need? Currently nothing is logged at the boundary.
- **Security:** What's the input-size ceiling for `oin_string`? A 10MB OIN string DoSes the regex preprocessor. No max-length enforced.
- **Determinism:** Distance Geometry is stochastic. Is there a `random_seed` parameter? If two users with the same input get different XYZ, is that a bug or feature?
- **Reproducibility:** No mention of locking RDKit/Molassembler/numpy versions in CI. A version drift mid-release reproduces W7.3.

### Missing Documentation
- No mention of updating `README.md` or user-facing docs for the new `backend` parameter.
- No mention of updating CLAUDE.md project memory for the v0.2.2 changes.
- No mention of updating the v3.6 inline format spec if direct parsing exposes ambiguities.

### Missing Operational Hooks
- No `migrate_v021_to_v022.py` script (even informational).
- No "smoke test" command for users to verify their install works post-upgrade.
- No release-checklist item to run a benchmark and capture numbers.

### Missing Architectural Hooks
- The `architecture.yml` doesn't yet have nodes for Polydentate or Permutation atoms. They're listed in the Draft PRD as "to be added" but no MiniPRD is responsible for adding them.
- No MiniPRD owns the `architecture.yml` reconciliation step. After 5 MiniPRDs land dirty, who runs `/hyper-audit` and reconciles?

---

## Section 11: Recommendations Priority Matrix

| Priority | Recommendation | Section |
|---|---|---|
| **P0 (Block ship)** | Resolve Q-OPEN-1 (metal detection) before MiniPRD #1 | §9 |
| **P0 (Block ship)** | Resolve Q-OPEN-6 (permutation port) before MiniPRD #2 | §9 |
| **P0 (Block ship)** | Default backend = `"legacy"` in v0.2.2; flip in v0.3.0 | P4.1 |
| **P0 (Block ship)** | Tighten cross-backend RMSD threshold: < 0.5Å mean OR < 1.0Å max-per-atom | P3.2, P8.3 |
| **P0 (Block ship)** | Pin `scine_molassembler == 3.0.*` (not `>= 2.0.0`) | P5.6 |
| **P1 (Pre-merge)** | Add explicit error for v2.4 input on direct backend | P3.3, Q-OPEN-2 |
| **P1 (Pre-merge)** | Add input-size ceiling for `oin_string` (DoS protection) | §10 |
| **P1 (Pre-merge)** | Reconcile Optional[Literal] typing in `generate(backend=...)` | P5.5 |
| **P1 (Pre-merge)** | Add operational backout policy (revert default in patch) | P1.3 |
| **P1 (Pre-merge)** | Resolve audit-staleness policy + audit recheck after each MiniPRD | P2.4 |
| **P1 (Pre-merge)** | Reorder Execution Checklist to match risk-section dependency claim | P7.4 |
| **P1 (Pre-merge)** | Add the four missing risk-register entries (RDKit, Python, Molassembler license, ProcessPoolExecutor) | §10 |
| **P2 (Pre-release)** | Mark `parse_oin_direct` as `@experimental` for v0.2.2 | P4.6, Q-OPEN-7 |
| **P2 (Pre-release)** | Specify thread-safety, memory, and import-time NFRs | §10 |
| **P2 (Pre-release)** | Replace "<200 lines" with cyclomatic complexity metric | P8.2 |
| **P2 (Pre-release)** | Add structured logging at backend-selection boundary | §10 |
| **P3 (Nice-to-have)** | Add `random_seed` parameter for DG determinism | §10 |
| **P3 (Nice-to-have)** | Add LRU cache option for `parse_oin_direct` | §10 |

---

## Conclusion

The Draft PRD is a competent operationalization of a known-bad subsystem fix, but it ships with three structural weaknesses:

1. **Asymmetric coexistence**: direct = blocking, legacy = fail-soft, default = direct. Inverts the safety pattern of "trust new only after burn-in."
2. **Unresolved foundational questions**: Q-OPEN-1 and Q-OPEN-6 are gating but not gated.
3. **Loose acceptance thresholds**: 1.0Å RMSD, ≥80% coverage, "<200 lines" — gameable and not science-grade.

The recommended path forward is to apply the P0 items (5 changes) before any MiniPRD starts; treat the P1 items as merge-gates; defer P2/P3 to v0.2.3 or v0.3.0.

If P0 items are addressed, this PRD is shippable. If they are not, ship will produce a v0.2.2 with a hidden default-flip, two unresolved load-bearing assumptions, and a maintenance pact that v0.3.0 will struggle to honor.

---

**End of Red Team Report.**

**Next Step:** Run `/hyper-resolve` to triage these vulnerabilities, generate the final SuperPRD, and compile the executable MiniPRDs.
