# Process Document: MiniPRD-D (WS-3) Draft PRD via /hyper-architect

**Generated:** 2026-07-04T16:52:35-07:00
**Session Focus:** Author the Draft PRD for MiniPRD-D — Phase 1 · WS-3 of the η-ligand round-trip recovery effort (G1 + G3: make `_stitch_multi_eta_fragment` return a real bonded RDKit mol) using the HACF `/hyper-architect` phase.

## Problem Statement

OIN-SMILES round-trips 3D metal-complex structures XYZ→OIN→XYZ→OIN. Six η-bound (Cp/indenyl) complexes fail that round trip; Phase 0 fixed one (Ferrocene), leaving five. The three ansa/bis-indenyl metallocenes (TiCat1/3/4) fail because the OIN→XYZ generator emits **no bond topology** for multi-eta fragments (`GeneratedStructure.mol is None`), forcing the re-encoder to re-perceive bonds from distorted coordinates — which garbles the string (TiCat1) or crashes (TiCat3/4). This session produced the requirements spec (Draft PRD) to fix that, without writing production code.

## Starting State

- **HEAD:** `d04c785` (“Docs: cut 0.3.0, document stereo round-trip arc”). Branch `main`, ~20 commits ahead of origin, not pushed (standing instruction).
- **Uncommitted working tree (Phase 0 · WS-0/1/2, landed earlier 2026-07-04, staged-for-review):** modified `src/oinsmiles/utils/oin_aligner.py` (WS-2 re-aromatizer), `src/oinsmiles/utils/xyz2mol.py` (WS-1 real error), `tests/integration/verify_roundtrip.py` (WS-0 `--only` + artifacts), `spec/worklog/NOTES.md`; untracked `TASK-40/41/42-*.md`, `ROUNDTRIP-eta-recovery-handoff.md`, `tests/fixtures/ticat3_generated_broken.xyz`, and two new unit test files.
- **Round-trip status:** 20/25 (`verify_roundtrip.py`). Remaining failures: TiCp2Me2, TiCat1–4.
- **Suites green:** `discover tests` 55 OK; `discover tests/unit` 127 OK (skip=3, xfail=0); `verify_xyz_to_oin.py` 25/25.
- **`spec/active/`:** empty (only `.gitkeep`) — no clobber risk.
- **Key defect (pre-session):** `OIN3DGenerator.generate(<TiCat1 OIN>).mol is None`; multi-eta fragments produce no combined RDKit mol.
- **Guiding docs:** `spec/worklog/ROUNDTRIP-eta-recovery-handoff.md` (diagnosis + roadmap; WS-3 = MiniPRD-D, scope G1+G3), `spec/worklog/NOTES.md` (2026-07-04 Phase-0 log entry).

## Approach & Methodology

Spec-driven, single HACF phase: `/hyper-architect` (Phase 1 of `architect → redteam → resolve → execute → audit`). The Architect role is codebase-first — resolve everything derivable from source before interviewing, and empirically ground claims rather than infer. Sequencing: (1) exhaustive read of the target function, its caller, and every downstream consumer; (2) a runtime probe to capture ground truth (current `mol=None`, real fragment SMILES, aromatic state, atom-count arithmetic); (3) a tight three-question interview limited to genuine judgment calls the brief flagged; (4) generate the Draft PRD in the repo’s house style. No production code was written — the deliverable is a requirements document that feeds `/hyper-redteam`.

## Steps Taken

1. **Read the skill + guiding docs.** Read `.agents/skills/hyper-architect/SKILL.md` (interview rules), `ROUNDTRIP-eta-recovery-handoff.md` §2 (root causes G1/G3), the NOTES.md 2026-07-04 Phase-0 entry, and `AGENTS.md`/`CLAUDE.md` (HACF is toolchain, metal is always `fragments[0]`). Purpose: establish scope (G1+G3 only; WS-4/5/6/7 out) before touching code.

2. **Mapped the generator module.** `grep` for function anchors in `src/oinsmiles/generation/molassembler_adapter.py`, then read `_stitch_multi_eta_fragment` (`:591-1133`) in full, `_assemble_combined_mol` (`:556-588`), the caller loop in `_template_generate` (`:2041-2370`), the single-eta sibling `_stitch_eta_fragment` (`:1136+`), and `generate()` (`:3356+`). Outcome: confirmed the `mol=None` return at `:1133`, the `heavy_atom_map` emission mapping (`:1050-1067`), and the de-aromatization in `_embed_fragment` (`:719-728`, = G3).

3. **Found the hidden caller contract.** `_assemble_combined_mol` (`:577-580`) wires each metal→binding DATIVE bond by indexing `frag_mol`’s own atoms with the fragment-SMILES `bidx` the caller passes (`:2088`). Since the reconstructed mol is in emission order (a permutation of `mol_h`), the function must **also return remapped emission-space binding indices** — the central design constraint.

4. **Traced the DoD consumers.** Read `get_oin_string` (`xyz2mol.py:751`), the round-trip step-2/3 branch (`verify_roundtrip.py:264-294`), and the WS-2 re-aromatizer (`oin_aligner.py:50-58`). Discovered WS-2 fires **only when both endpoint atoms are already `IsAromatic`** — which `_embed_fragment` clears — resolving the brief’s “is G3 still needed?” question: atom-flag restoration is mandatory regardless of WS-2.

5. **Ran a runtime probe** (`scratchpad/probe_multieta.py`, `uv run python`). Captured ground truth: `generate(<TiCat1 OIN>).mol is None`; parsed fragments `[Ti]` / `C[Si](C)(c1[cH][cH][cH][cH]1)c1[cH][cH][cH][cH]1` / `[CH3]` / `[CH3]`; the ansa fragment’s `sanitize=False` parse retains 10 aromatic atoms + 10 aromatic bonds (so `mol` at `:650` is a sound G3 source); 27 emission atoms == 27 `mol_h` atoms (bijection feasible).

6. **Blast-radius grep** across `verify_roundtrip.py`/`rmsd_utils.py`. Key finding: the RMSD path **already consumes `gen_result.mol`** (`:333` `mol_gen_bonded = gen_result.mol if not None else mol_gen_xyz` → `:381` `calculate_tmc_rmsd(..., mol2_bonded=...)`). Making the mol non-None silently reroutes coordination-sphere extraction through the reconstructed bonded mol — so a binding-index error would corrupt RMSD, not just the string.

7. **Interview — Q1 (G3 scope).** Recommended full aromatic restoration (atom flags + bond types) over the minimal atom-flags-only variant. **User accepted.**

8. **Interview — Q2 (degrade posture).** Recommended `logger.warning` on coverage-guard miss (never `warnings.warn`/`OINStereoWarning`, to avoid the `-W error` coupling that bit TASK-31/32) plus soft-degrade (`mol=None`, XYZ kept, never bare `None`). **User accepted.**

9. **Interview — Q3 (acceptance model).** Recommended enforced `tests/unit` invariant pins (mol-not-None, 12 DATIVE bonds, aromatic ring bonds, conformer/atom-count alignment, `normalize(OIN1)==normalize(OIN2)`) + integration confirmation via `verify_roundtrip.py`, RMSD measured-not-gated, and no candidate artifacts (all outputs deterministic or self-checking). **User accepted.**

10. **Matched house style.** Read the deprecated schema templates (`.agents/schemas/{SuperPRD,MiniPRD}_template.md`), the most recent Draft PRD (`spec/archive/20260704_000832_EtaRingCanonicalization/Draft_PRD.md`), and `grep`ed `spec/compiled/architecture.yml` for the correct node IDs (`atom_molassembler_adapter`, `atom_generated_structure`, `atom_oin3d_generator`, `atom_xyz2mol`).

11. **Wrote the Draft PRD** to `spec/active/Draft_PRD_MiniPRD_D_MultiEtaMol.md` (WS-3-specific filename, no clobber): 9 sections including a Resolved Trade-offs Log (D-1/D-2/D-3), a return-signature contract (4-tuple → 5-tuple), blast radius, 6 user stories, risks R1–R6, and a §9 “Notes for /hyper-redteam” seeding the next phase. No `.py` files edited, so `ruff format` did not apply.

12. **Appended a NOTES.md Log entry** recording the codebase-first findings, the empirical probe, the RMSD-consumer finding, and the three decisions, ending with the `/hyper-redteam` handoff.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| G3 = **full** aromatic restoration (atom `IsAromatic` + `AROMATIC` bond types) from the `mol` parse | Atom-flags-only, lean on WS-2 re-aromatizer for bond types | WS-2 only fires when atoms are already aromatic, and `_embed_fragment` clears them — atom flags are mandatory regardless. Full restore is one extra `SetBondType` per bond, makes `gen.mol` valid for every consumer (RMSD sphere, SDF writers), and turns WS-2 into a provable no-op instead of a hidden dependency. |
| Graceful degrade = `logger.warning` + soft return (`mol=None`, XYZ kept) | (a) silent `mol=None`; (b) `warnings.warn`; (c) bare `None` | `logger.*` can’t trip the `-W error::OINStereoWarning` gate that bit TASK-31/32; soft return preserves the XYZ (bare `None` would abort placement → DG fallback, worse for haptics). Never emit a wrong mol. |
| Return signature 4-tuple → 5-tuple (add emission-space binding indices) | Return the full `heavy_atom_map`; keep 4-tuple and have caller recompute | Plain tuple matches sibling conventions; returning remapped indices directly is the minimal correct fix for `_assemble_combined_mol`’s hidden `frag_start + bidx` contract. |
| Acceptance = enforced unit pins + integration confirmation; no candidate artifacts; RMSD measured-not-gated | Integration-only (`verify_roundtrip.py`); add RMSD gate | Topology is deterministic and string acceptance is a self-checking identity, so it belongs in the fast `discover tests/unit` loop; RMSD quality is WS-4’s concern. |
| Coverage guard: exactly-2-rings + SiMe2 bridge + total `output→mol_h` bijection + `etkdg_ok` | General/looser bridge acceptance | The command scoped MiniPRD-D to SiMe2 bis-Cp/indenyl; fail-closed on anything else prevents a mis-mapped mol. |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| Draft PRD — MiniPRD-D (WS-3) | `spec/active/Draft_PRD_MiniPRD_D_MultiEtaMol.md` | created |
| Worklog log entry | `spec/worklog/NOTES.md` | updated (appended 2026-07-04 architect entry) |
| Probe script (throwaway) | `scratchpad/probe_multieta.py` | created (session scratchpad, not committed) |

No production/source code (`src/`) or tests changed this session. Nothing staged, committed, or pushed.

## Results & Outcomes

- A complete, code-grounded Draft PRD for MiniPRD-D exists at `spec/active/Draft_PRD_MiniPRD_D_MultiEtaMol.md`, ready for `/hyper-redteam`.
- All three genuinely-open design decisions are resolved with the user and logged (§5.1 D-1/D-2/D-3).
- The specification is empirically anchored (probe confirmed `mol=None`, the fragment SMILES, aromatic state, and the 27=27 atom bijection) rather than inferred.
- A previously-unstated hazard is captured for the Red Team: the RMSD path already reads `gen_result.mol`, so binding-index correctness now also governs RMSD, not just the re-encoded string.
- Definition of Done is restated precisely: TiCat1/3/4 `gen.mol is not None` with correct topology, step-2 re-encode via `get_oin_string` (which also removes the TiCat3/4 crash), normalized string identity — **not** RMSD (WS-4). Ferrocene and all suites stay green.

## How to Reproduce

Prerequisites: repo at `main` with the uncommitted Phase-0 (WS-0/1/2) tree present; `uv` toolchain; Python ≥3.10; RDKit + scipy installed (`uv sync`).

1. Read, in order: `spec/worklog/ROUNDTRIP-eta-recovery-handoff.md` (§2 G1/G3, §3 WS-3, §5 policy), the NOTES.md 2026-07-04 Phase-0 entry, `CLAUDE.md`/`AGENTS.md`.
2. Read the target function and its collaborators in `src/oinsmiles/generation/molassembler_adapter.py`: `_stitch_multi_eta_fragment` (`:591-1133`), `_assemble_combined_mol` (`:556-588`), the caller in `_template_generate` (`:2041-2370`), plus `get_oin_string` (`src/oinsmiles/utils/xyz2mol.py:751`) and the WS-2 re-aromatizer (`src/oinsmiles/utils/oin_aligner.py:50-58`).
3. Probe ground truth: run a script that calls `OIN3DGenerator().generate("[Ti_TET].C[Si](C)(c{0}1[cH]{0}[cH]{0<}[cH]{0}[cH]{0}1)c{1}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1.[CH3]{2}.[CH3]{3}")` and prints `.mol is None` and the parsed fragments (via `generation.oin_parser.OINParser().parse(...)`). Expect `mol is None: True`, 36-atom XYZ, and 10 aromatic atoms/bonds on the `sanitize=False` ansa parse.
4. Confirm the RMSD consumer: `grep -n "gen_result.mol\|mol2_bonded" tests/integration/verify_roundtrip.py` → lines 333/381.
5. Invoke `/hyper-architect` with the WS-3 goal. Answer three interview questions: G3 scope → full restoration; degrade → `logger.warning` + soft return; acceptance → unit pins + integration, no candidate artifacts.
6. Confirm the generated Draft PRD lands at `spec/active/Draft_PRD_MiniPRD_D_MultiEtaMol.md` and a NOTES.md log entry is appended.

Gotchas / order-dependencies:
- Line numbers in the handoff are approximate — re-verify anchors against current source (they had drifted from the WS-0/1/2 edits).
- `mkdir` is unnecessary for `spec/process`/`spec/active` (they exist); `archive_specs.py` needs `uv run python` (bare `python` isn’t on PATH).
- Keep `git add` scoped; do not push or commit (standing instruction). The draft is left unstaged for the chain.

## Patterns & Lessons

- **Codebase-first + empirical probe beats inference.** Reading the function is necessary but not sufficient; a 20-line probe converted “the mol is probably `None`” and “the fragment is probably aromatic” into measured facts (27=27 atoms, 10 aromatic bonds) that the PRD cites as ground truth and forbids re-deriving.
- **Hunt the hidden caller contract.** The load-bearing design constraint (remapped binding indices) was invisible in the target function — it lived in how `_assemble_combined_mol` indexes `frag_start + bidx`. Always trace how a return value is *consumed*, not just produced.
- **Re-examine “subsumed” claims mechanically.** The brief suggested WS-2 might subsume G3 for the string; reading WS-2’s exact predicate (`GetBeginAtom().GetIsAromatic()`) proved atom-flag restoration is still mandatory. Verify the trigger condition, don’t trust the summary.
- **Follow the mol to every reader.** Flipping one field (`mol=None` → real mol) changed a *silent* downstream path (RMSD coordination-sphere extraction), not just the obvious one (string re-encode). Enumerate all consumers of a value before changing its nullability.
- **Interview only genuine judgment calls.** The command pre-decided scope and most mechanics; the interview was three questions, each an either/or the codebase couldn’t answer, each with a recommended default — keeping the human turn count minimal while still surfacing the decisions the Red Team must see.
- **Match house style from the newest example, not the deprecated template.** The schema templates are marked deprecated; the most recent archived Draft PRD is the authoritative structure to mirror.
