# MiniPRD: Zone-A P Stereocenter Encoding — Generation Enforcement (MiniPRD-B)
**Hypergraph Node ID:** atom_molassembler_adapter
**Parent Node:** mod_generation (`src/oinsmiles/generation/molassembler_adapter.py`)
**Parent SuperPRD:** `spec/compiled/SuperPRD_StereoPhase4_ZoneA_P.md` (v1.0.0)
**Execution order:** SECOND — consumes the `_OIN_CIPCode_LP` / `[P@]` contract
established by `MiniPRD_ZoneA_P_Encode.md`. Blocked until MiniPRD-A lands.

## 1. The Confidence Mandate
**Agent Instruction:** Before generating any plans or writing code, analyze this
document and output a Confidence Score (1-10). If the score is below 9, list
strictly the clarifying questions needed to reach 10.

Context an executor must load: `generation/molassembler_adapter.py` —
`_template_generate` (~line 977, assembled-complex stage) and
`_molassembler_worker` (**line 1489** — the Draft PRD's `:1418` was stale),
SuperPRD §5.1 (data flow + Resolved Trade-offs B2/B3, B8) and §6. The dummy-metal
copy + rdCIPLabeler helpers from MiniPRD-A (`core/chirality.py`) are reused, not
reimplemented.

## 2. Atomic User Stories
* **US-B1:** As the generation pipeline, the input `[P@]`/`[P@@]` tag is
  preserved through parse (already true — verify only) and enforced after
  assembly: the assembled complex's lone-pair-convention label (dummy-metal
  copy + `rdCIPLabeler`) must match the input tag; on mismatch, re-embed the
  fragment with a new ETKDG seed (max 3 attempts, **never** a mirror); on
  persistent mismatch, emit the structure + `OINStereoWarning`. Paths with no
  assembled RDKit mol (eta fallback, DG fallback) skip enforcement + warn.
* **US-B2:** As the round-trip contract, DIPAMP is lossless: XYZ→OIN→XYZ→OIN is
  byte-stable on the OIN string (pinned RDKit version), and `rdCIPLabeler`
  CIP-from-3D on the regenerated metal-present complex matches the original
  (R,R); flipping `@↔@@` in the OIN inverts both regenerated labels.
* **US-B3:** As the fallback path, molassembler either honors trivalent `[P@]`
  via `from_smiles`, or `_molassembler_worker` sets an explicit atom
  stereopermutator — decided by an in-MiniPRD investigation (Q3), with the
  finding recorded in the worklog.

## 3. Implementation Plan (Task List) — 100% complete
- [x] Task 1: Verify (audit, not edit) that `parse_inline_string` and
      `OINParser.parse` pass `[P@]{0}` / `[P@@]{1>}` through intact — MiniPRD-A
      Test 8 covers the unit level; here assert at `ParsedOIN` level.
      Confirmed via `TestZoneAPParsedOINPassthrough` (new test file) — no
      production code change required, `parse_inline_string`'s `SLOT_REGEX`
      only strips `{slot}` markers, `@`/`@@` inside bracket atoms untouched.
- [x] Task 2: In `_template_generate`, at the assembled-complex stage
      (re-located: the `combined_mol` construction block, now factored into
      `_assemble_combined_mol`), added
      `_verify_zone_a_p(assembled_mol, fragment_inputs) -> list[int]`: for
      each input fragment P that carried a chiral tag, builds the dummy-metal
      copy (reused MiniPRD-A helper from `core/chirality.py`, imported not
      copied), gets the `rdCIPLabeler` label via `_lp_cip_label`, compares to
      the input tag's expected label (`_zone_a_p_expected_labels`, a
      graph-based recompute off the OIN fragment SMILES's own tag); returns
      mismatched P indices.
- [x] Task 3: Added the bounded enforcement loop in `_template_generate`:
      while mismatches exist and attempts < 3 → re-run ETKDG on the
      offending fragment(s) with a NEW seed (`_stitch_fragment(..., seed=…)`),
      re-place, re-verify via `_assemble_combined_mol` + `_verify_zone_a_p`.
      Hard cap 3 attempts total.
- [x] Task 4: On persistent mismatch after 3 attempts: structure is emitted
      and `warnings.warn(OINStereoWarning, f"... atom {idx} ...")` fires once
      per offending atom. Verified worst-case runtime empirically (persistent
      forced-mismatch test): well under 1s, far below the 60s
      ProcessPoolExecutor budget.
- [x] Task 5: On paths where no assembled RDKit mol exists (eta fallback →
      `GeneratedStructure.mol is None`, DG fallback): enforcement is skipped
      (gated on `combined_mol is not None`) and `_warn_zone_a_p_fallback()`
      (called from `MolassemblerAdapter.generate()`, both the
      no-assembled-mol template branch and the DG-fallback branch) emits
      `OINStereoWarning` ("stereo unenforced on fallback path …") — RISK-9 is
      observable, never silent.
- [x] Task 6: Q3 investigation done — see `spec/worklog/NOTES.md` Log entry
      "MiniPRD_ZoneA_P_GenEnforce.md (MiniPRD-B) executed". Finding:
      `masm.io.experimental.from_smiles` does not silently drop the
      stereopermutator on trivalent `[P@]` — it raises `RuntimeError:
      Mismatched shape for set chiral data`, already caught by
      `_molassembler_worker`'s pre-existing bare `except Exception:` which
      routes to `_rdkit_etkdg_fallback` (correctly respects the chiral tag).
      No `_molassembler_worker` code change made (mol object never exists to
      set a stereopermutator on); `_warn_zone_a_p_fallback` added instead so
      the "never actually verified" gap on this path is visible.
- [x] Task 7: Added the test-only injection point:
      `_stitch_fragment(..., _test_flip_chiral_idx=<local_p_idx>)` flips ONE
      atom's chiral tag before ETKDG embeds — a genuine, localized mis-embed,
      never a whole-fragment mirror. Exercised by
      `TestZoneAPForcedMisEmbedCorrection` and `TestZoneAPBoundedFailure`.
- [x] Task 8: Tests added — `tests/unit/test_zone_a_p_genenforce.py` (11
      tests, §5 Tests 1-6). `_stitch_fragment` gained an explicit `seed`
      parameter (default 42, unchanged behaviour) used by all generation
      assertions. Full suite green:
      `uv run python -m unittest discover tests/unit` → 112 run, OK
      (skipped=3, expected failures=2, unchanged from before this MiniPRD);
      `uv run python -m unittest discover tests` → 55 run, OK.

## 4. The Negative Space (Constraints)
* **DO NOT** apply a mirror or any improper transform to a fragment — the
  correction mechanism is re-embed-with-new-seed ONLY (Resolved B2/B3: a global
  mirror inverts every co-resident stereocenter → diastereomer manufacture).
* **DO NOT** exceed 3 re-embed attempts or loop unboundedly — persistent
  mismatch emits + warns; a stereo warning must never become a timeout death.
* **DO NOT** place enforcement in `OIN3DGenerator` (engine level) —
  `GeneratedStructure.mol` is `Optional` and post-hoc; the adapter's
  assembled-complex stage is the pinned layer (Resolved Q2/RISK-2).
* **DO NOT** reimplement the dummy-metal copy or CIP labeling — import the
  MiniPRD-A helpers; `rdCIPLabeler` end-to-end, never legacy `_CIPCode`.
* **DO NOT** compare labels cross-convention — enforcement compares
  lone-pair-sense vs lone-pair-sense (input tag vs dummy-copy label).
* **DO NOT** assert on raw 3D coordinates in tests — assert on derived
  metal-present CIP (deterministic oracle over non-deterministic ETKDG).
* **DO NOT** silently skip enforcement anywhere — every skipped path warns.
* **DO NOT** regress carbon `@/@@` or achiral generation paths — the adapter is
  shared; TASK-10 + baseline generation tests must stay green.
* **DO NOT** touch direct-parser nodes, `oin/parser.py`, or the OIN grammar.

## 5. Integration Tests & Verification
* **Test 1 (Deterministic — enantiomer discrimination):** generate from the
  DIPAMP OIN (fixed seed) → `rdCIPLabeler` metal-present CIP on both P = (R,R).
  Flip `@↔@@` on both P in the OIN → regenerate → both labels inverted (S,S).
* **Test 2 (Deterministic — lossless round-trip):** XYZ→OIN→XYZ→OIN on
  `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz` → second OIN string byte-identical to
  the first (pinned RDKit version); regenerated complex CIP matches original.
* **Test 3 (Deterministic — forced mis-embed):** via the Task-7 injection point,
  supply a mirrored initial P pyramid → enforcement corrects it within ≤3
  re-embed attempts; assert NO mirror transform was applied (co-resident
  stereocenter in the test ligand retains its configuration).
* **Test 4 (Deterministic — bounded failure):** injection that can never
  satisfy the tag → generation completes, structure emitted, exactly one
  `OINStereoWarning` per offending atom (`catch_warnings(record=True)`),
  wall-clock far below the 60 s timeout.
* **Test 5 (Deterministic — fallback observability):** an eta-fallback complex
  (mol=None path) with a tagged P in the OIN → generation completes +
  `OINStereoWarning` stating enforcement was skipped.
* **Test 6 (Deterministic — no regression):** baseline generation suite
  (cisplatin, ferrocene, Ir(ppy)3, TASK-10 carbon set) unchanged; clean
  fixtures pass under `-W error::OINStereoWarning`.
* **Test 7 (Novel — molassembler Q3):** the Q3 investigation output (does
  `from_smiles` keep trivalent `[P@]`?) → **Candidate Artifact routing
  protocol triggered** — finding recorded in `spec/worklog/NOTES.md` and the
  chosen fallback mechanism noted for `/hyper-audit`.
