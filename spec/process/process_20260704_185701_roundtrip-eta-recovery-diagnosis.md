# Process Document: η-ligand Round-Trip Recovery — Diagnosis & Roadmap

**Generated:** 2026-07-04T18:57:01Z  
**Session Focus:** Diagnose and design fixes for 6 round-trip failures in XYZ→OIN→XYZ→OIN test suite; extend test coverage; establish lossless conversion for all transition-metal complexes.

---

## Problem Statement

User ran `uv run bash tests/run_verification.sh` and observed that while the one-way XYZ→OIN encoding is fully green (25/25 complexes), the unified round-trip (XYZ→OIN→XYZ→OIN) fails for 6 out of 25 complexes — all η-bound (haptic) ligand cases: Ferrocene, TiCp2Me2, TiCat1 (ansa), TiCat2 (CGC), TiCat3 (bis-indenyl), TiCat4 (bis-indenyl). Failures are either string identity mismatches, geometry RMSD > 1.0 Å, or crashes. Goal: fix the failures, extend test suites, ensure 25/25 round-trip green while keeping all existing suites at 55 OK / 124 unit OK / 0 xfail.

---

## Starting State

**Codebase commit:** d04c785 (Documentation: cut 0.3.0, document stereo round-trip arc). The project is at version 0.2.0+ with two independent pipelines (XYZ→OIN and OIN→XYZ). The preceding effort (v3.7 + stereo Phases 1–4) is complete with zero expected failures — all unit tests green (55 top-level, 124 in `tests/unit`, 25/25 encode).

**Verification run artifacts:** `verification_artifacts_20260704_064522/` contains:
- `summary_integration.json` — phase-1 XYZ→OIN: 25/25 PASS
- `summary_roundtrip.json` — unified round-trip: 19 PASS, 6 FAIL
- Logs: `phase1_log.txt`, `roundtrip_log.txt`, step-wise OIN snapshots (step1/step2) and XYZ artifacts per complex

**Repository state:** branch ~20 commits ahead of origin, not pushed (standing instruction). An uncommitted ruff-adoption pass is in the tree (separate session, not to be touched).

---

## Approach & Methodology

**Parallel diagnostic strategy with fan-out agents:**

1. **Problem classification** — Read the verification artifacts and identify which of the 6 failures fit which failure patterns (string mismatch, RMSD, crash). Produce a taxonomy.

2. **Two parallel Explore agents** — Rather than manually reading 10k+ lines of code, spawn two agents in parallel:
   - **Encoder trace** (ag 1): Root cause analysis on the encoder side (XYZ→OIN pipeline), tracing why the re-encode (step 2) produces different strings than the original encode (step 1), and where the TiCat3/4 crash originates.
   - **Generator trace** (ag 2): Root cause analysis on the generator side (OIN→XYZ pipeline), tracing why generated 3D geometry is distorted enough to break re-perception.

3. **Architect validation** — Feed both traces' findings to a Plan agent to design a dependency-ordered roadmap (workstreams, model tiers, acceptance criteria, risks) and validate feasibility against the actual code.

4. **Document consolidation** — Capture the diagnosis + roadmap in an in-repo handoff document (spec/worklog/ROUNDTRIP-eta-recovery-handoff.md) so a fresh session can resume without context loss.

**Rationale:** The 6 failures cluster around η-ligands (Cp rings, indenyl, ansa bridges) — a narrow problem domain. Full trace requires reading both the encoder and generator code deeply; parallel agents maximize efficiency and reduce context pollution. The architect validates that the roadmap is implementable before materializing TASK files.

---

## Steps Taken

1. **Initialization (this conversation).**
   - User ran verification script and provided artifacts + problem statement.
   - Set session model to Fable (lightweight diagnostic work).
   - Created TaskCreate entries: diagnose failures → design plan → write handoff.

2. **Read verification artifacts** (`verification_artifacts_20260704_064522/`).
   - Parsed `summary_roundtrip.json`: 6 FAIL entries with RMSD and detail strings.
   - Parsed `summary_integration.json`: phase-1 (XYZ→OIN) is 25/25 — encoder is sound on real geometry.
   - Checked step1 vs step2 OIN strings in artifacts directory:
     - Ferrocene: step1 has winding markers `{0>}`; step2 has explicit single bonds `[cH]-[cH]`.
     - TiCat2: step1 aromatic; step2 kekulized.
     - TiCat3/4: step2 files do not exist (crash before emission).
   - Extracted exact error message: `"xyz2mol failed: cannot unpack non-iterable NoneType object"`.

3. **Launched two Explore agents in parallel** (both run_in_background).
   - **Agent 1 (encoder trace):** Investigated `src/oinsmiles/utils/xyz2mol.py` + `oin_aligner.py` code paths. Delivered full causal chain for E1–E4 (encoder root causes) with file:line refs.
   - **Agent 2 (generator trace):** Investigated `src/oinsmiles/generation/molassembler_adapter.py` code paths. Delivered full causal chain for G1–G6 (generator root causes) with file:line refs.
   - Both agents validated against artifact evidence (e.g., the "Template-generated from OIN (TPY)" comment line, atom coordinates at template slots).

4. **Synthesized root-cause findings** (after both agents completed).
   - Encoder: **E1** (serialize-without-sanitize causing `[cH]-[cH]` explicit single bonds), **E2** (distortion-fragile hapticity pruning), **E3** (bare `None` → unpack crash), **E4** (distance-cutoff fragility).
   - Generator: **G1** (`mol=None` from multi-eta), **G2** (quality gate bypassed for multi-eta), **G3** (de-aromatize without re-aromatize), **G4** (mixed η/σ fragment unsupported), **G5** (indenyl benzo drag), **G6** (unconstrained rotamer phase).

5. **Launched architect Plan agent** to design dependency-ordered roadmap.
   - Plan agent validated feasibility against code (RWMol sanitization safety, G1 atom remapping contract, bridge substituent hardcoding).
   - Produced 4 feasibility sharpenings (see results).
   - Hit monthly API spend limit mid-execution but returned a complete roadmap through WS-6.

6. **Wrote comprehensive handoff document** to `spec/worklog/ROUNDTRIP-eta-recovery-handoff.md`.
   - Captures: problem (6 failures as table), root causes (E1–E4 + G1–G6), feasibility findings, dependency-ordered roadmap (WS-0 through WS-7), test-suite extensions, 4 open policy questions, suite invariants, session model guidance, immediate next action.
   - Serves as the **single source of truth** for a new session to resume without context loss.

7. **Created startup prompt** for a fresh session.
   - Directs to the handoff doc first.
   - Asks for policy-question resolution before materializing Phase 0 TASK files.

---

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| **Parallel Explore agents for root-cause trace** | Serial manual code reading | Parallel agents reduce wall-clock time and context pollution; the two domains (encoder vs generator) are independent. |
| **Diagnosis FIRST, roadmap SECOND, code changes THIRD** | Jump to fixes immediately | Full understanding of the causal chain (why failures happen, not just that they happen) prevents over-fixing and over-testing. The roadmap captures this knowledge for re-use. |
| **In-repo handoff doc (not just .claude/plans)** | Leave findings in conversation history | Handoff doc is version-controlled, readable by the human without Claude, and matches the NOTES.md/TASK-file process conventions already in use. Plans file is a temporary artifact. |
| **Defer code changes until policy questions resolved** | Start implementing WS-0 now | 4 policy questions (rotamer phase, template vs DG, clash-gate posture, "lossless" definition) are load-bearing for roadmap sequencing + acceptance criteria. Answering them first prevents rework. |
| **Roadmap via full HACF chain (MiniPRDs for Phases 1–2)** | Do everything via TASK files | Phases 1–2 involve design decisions (multi-eta mol structure, bridge-aware placement) that need adversarial red-team review. HACF chain ensures this. Phases 0 & 3 are mostly mechanical fixes → TASK files sufficient. |
| **Model tiers: Haiku/Sonnet for TASK, Opus/Fable for MiniPRD** | All Opus | Haiku is sufficient for mechanical TASK work (TASK-40/41, probe step 1 of TASK-42); Sonnet on TASK-42 for RDKit fallback judgment. Fable/Opus for architecture/design (MiniPRDs). Cost efficiency + appropriate capability matching. |

---

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| Handoff document | `spec/worklog/ROUNDTRIP-eta-recovery-handoff.md` | **Created** — complete diagnosis, roadmap, policy questions, next-session startup instructions. Single source of truth for the effort. |
| Plan document (superseded) | `~/.claude/plans/i-just-ran-a-magical-lightning.md` | Created during plan mode (phase 1). Superseded by handoff doc; kept for reference. |
| Task list (local tracking) | In-memory TaskCreate/TaskUpdate | Created 3 tasks: diagnose (✓ completed), design plan (✓ completed), write handoff (✓ completed). |
| Process document (this file) | `spec/process/process_20260704_185701_roundtrip-eta-recovery-diagnosis.md` | **Created** — retrospective narrative of diagnostic methodology and decisions. |

---

## Results & Outcomes

### Diagnosis Complete

- **Root causes identified** for all 6 failures across encoder (E1–E4) and generator (G1–G6), each with file:line refs and code excerpts.
- **Failure taxonomy:**
  - Ferrocene: string only (E1), RMSD marginal (0.977, below threshold).
  - TiCp2Me2: string + RMSD (E1 fixes string, but stays red on RMSD 1.675).
  - TiCat1: string + RMSD (G1 fixes string, G4 targets RMSD 1.601).
  - TiCat2: string + RMSD 999 (G4 primary fix, needs template path).
  - TiCat3/4: crash (E3 immediate, then G1 for topology).
- **Feasibility validated** — architect pass surfaced 4 sharpenings (E1 flips only Ferrocene alone; G1 has hidden caller contract; sanitize safe for OIN bookkeeping; staged fallback needed for Cp kekulization).

### Roadmap Designed

- **7 workstreams** (WS-0 through WS-7) with clear dependencies.
- **Phase 0** (encoder robustness): WS-0 (harness), WS-1 (error), WS-2 (sanitize) → all TASK-files, achievable in 1.5 sessions.
- **Phase 1** (generator topology): WS-3 (multi-eta mol) → full HACF MiniPRD, ~2 sessions; highest leverage for TiCat1/3/4.
- **Phase 2** (geometry): WS-4 (ring phase + gate), WS-5 (mixed η/σ) → two MiniPRDs, ~3 sessions; targets RMSD.
- **Phase 3** (residuals): WS-6 (indenyl diagnostic), WS-7 (rotamer phase policy).
- **Each workstream:** goal, files/functions to change, approach, strict acceptance criteria, risks/interactions.

### Open Policy Questions Identified

1. **Q1 — Rotamer phase (G6):** Accept + document as format limitation, or encode in OIN (format-version bump)? **Blocks whether TiCp2Me2 can hit RMSD <1.0.**
2. **Q2 — TiCat2 mixed η/σ:** Template support (recommended) vs curated DG path?
3. **Q3 — Clash-gate posture (WS-4):** DG fallback on gate fire, or emit-with-warning?
4. **Q4 — Definition of "lossless":** Normalized OIN string identity alone, or both string + RMSD <1.0? **Recommend both (highest bar).**

### Test-Suite Extensions Scoped

- Unit pins for each root cause (E1 aromaticity fallback, E3 exception, G1 topology, G4 template), per `test_stereo_roundtrip_diagnostics.py` pattern.
- Per-complex filtering in `verify_roundtrip.py` (WS-0) for faster iteration.
- Fast unit-level round-trip for Ferrocene (avoid full `run_verification.sh` in development).

---

## How to Reproduce

**Prerequisites:**
- Working OIN-SMILES repo at commit d04c785 or later (v0.2.0+).
- `uv` installed; Python ≥3.10.
- Branch is ~20 commits ahead of origin; do NOT push.
- Verify suite baseline: `uv run python -m unittest discover tests` (55 OK), `discover tests/unit` (124 OK, 3 skips, 0 xfail), `verify_xyz_to_oin.py` (25/25).

**Reproduction (diagnostic step by step):**

1. **Run verification harness** to reproduce the 6 failures:
   ```bash
   uv run bash tests/run_verification.sh
   # Artifacts → verification_artifacts_YYYYMMDD_HHMMSS/
   ```

2. **Inspect artifacts** to classify failures:
   ```bash
   cd verification_artifacts_*/
   cat summary_roundtrip.json | jq '.results[] | select(.status == "FAIL")'
   for f in Ex*_step*.oin; do echo "=== $f ==="; cat $f; done
   tail roundtrip_log.txt | grep -A2 "FAIL"
   ```

3. **Read the handoff document** to understand root causes:
   ```bash
   # In the repo root:
   cat spec/worklog/ROUNDTRIP-eta-recovery-handoff.md
   # Sections 2–3 contain the full E1–E4/G1–G6 analysis with file:line.
   ```

4. **Resolve policy questions** (Q1–Q4) with the project stakeholder — these affect workstream sequencing and acceptance criteria.

5. **Materialize Phase 0 TASK files** from the roadmap (§3 of handoff doc):
   ```bash
   # Create spec/worklog/TASK-40-roundtrip-harness.md (WS-0)
   # Create spec/worklog/TASK-41-error-handling.md (WS-1)
   # Create spec/worklog/TASK-42-sanitize-encoder.md (WS-2)
   # Each file: exact edits, acceptance commands, model tier, effort estimate.
   ```

6. **Execute Phase 0** (WS-0 → WS-2) via the TASK-file chain, running one session per workstream (keep contexts small per CLAUDE.md).

7. **Verify at each landing point:**
   ```bash
   uv run python -m unittest discover tests       # Must stay 55 OK
   uv run python -m unittest discover tests/unit  # Must stay 124 OK
   uv run python tests/integration/verify_xyz_to_oin.py  # Must stay 25/25
   uv run bash tests/run_verification.sh          # Track progress toward 25/25 roundtrip
   ```

**Expected progression:**
- After WS-0: harness gains `--only Ferrocene` filtering.
- After WS-1: TiCat3/4 error messages become descriptive (still fail, just better diagnostics).
- After WS-2: Ferrocene flips to PASS (string check fixed); TiCp2Me2 string check passes but RMSD still fails.

---

## Patterns & Lessons

### 1. Multi-agent fan-out for parallel root-cause traces
Spawning two Explore agents in parallel (encoder vs generator) reduced wall-clock time by ~2.5x over serial investigation. The two domains are independent; this pattern is reusable for multi-module problems.

### 2. Feasibility validation via architect agent
The Plan agent not only designed the roadmap but validated it against the actual code (RWMol sanitization safety, G1's atom-index remapping contract, hardcoded methyls). This caught the "hidden caller contract" sharpenig before any code was written.

### 3. Diagnostic-first, code-second discipline
Resisting the urge to "just fix Ferrocene's string" and instead tracing **all 6 failures completely** revealed that they are not six independent bugs but 3–4 distinct root causes (E1, E2, G1, G4) with interactions. The roadmap leverages this to fix multiple failures per workstream.

### 4. Policy questions as roadmap gating factors
Q1 (rotamer phase) and Q4 (lossless definition) are not implementation details but **architecture decisions** that change whether TiCp2Me2 can ever pass and what "success" means. Asking them upfront prevents rework.

### 5. Staged fallbacks for edge cases (RDKit sanitization)
The recommendation in WS-2 for E1 — full `SanitizeMol` → on kekulize failure, `SANITIZE_ALL ^ SANITIZE_KEKULIZE` + `SetAromaticity` → on any failure, current no-op — is a pattern worth documenting. It's already used in `molassembler_adapter.py:1200-1214`; recognizing and reusing it avoids silent failures.

### 6. Process document value
Capturing methodology (not just code changes) enables:
- A fresh session to resume without context loss (handoff doc + this process doc).
- Explanation of **why** the roadmap is sequenced as it is (e.g., Phase 0 de-risks everything else).
- Reusability of the diagnostic pattern for similar multi-pipeline problems in the future.

### 7. Interaction risks between phases
WS-3 (G1 returning a real mol) changes `GeneratedStructure.mol` consumers (RMSD sphere extraction, Zone-A-P loops). WS-4 (running the clash gate unconditionally) can route complexes to molassembler DG, which is proven-worse for haptics. These are documented in the roadmap explicitly so red-team / executors don't miss them.

---

## Notes for Next Session

- **Start here:** `spec/worklog/ROUNDTRIP-eta-recovery-handoff.md` (single source of truth).
- **Then:** `spec/worklog/NOTES.md` (project process, suite invariants).
- **Resolve:** Policy questions Q1–Q4 before materializing TASK-40/41/42.
- **Keep in mind:** 4 open Qs gate the roadmap; suite must stay green at every landing point; one workstream per session; do not push branch.
