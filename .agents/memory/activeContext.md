
# Active Context

## Purpose
Captures the current state of OIN-SMILES development. Updated after significant task completions. Read first to understand where to pick up.

## Current State (as of 2026-07-08)

### Release Status
- **v0.2.0** — Released 2026-03-07: Molassembler backend, P/N stereocenter encoding, CLI, OIN v3.6
- **v0.3.0** — In `CHANGELOG.md [0.3.0] - 2026-07-04` (pyproject synced): OIN v3.7 descriptor-free metal token; stereo round-trip arc (winding preservation, haptic-face control, Zone-A P `[P@]`/`[P@@]` encoding + square-planar enforcement); eta-ring canonicalization; bidentate incompatible-bite → DG routing; plus the prior v0.2.1 work (CLI fix, P/N fixtures, TiCat ETKDG fix, Direct Parser audits). Not yet pushed/tagged.
- **v0.3.1** — In `CHANGELOG.md [0.3.1] - 2026-07-05`: MACE MLIP optimizer support, FF convergence knobs, Distinct RMSD error codes, and stability fixes for Ir drift, BDPP/BDNN, FeCO5/FeH2CO4. (Note: `pyproject.toml` was NOT bumped for 0.3.1 — reconciled to 0.3.2 below.)
- **v0.3.2** — In `CHANGELOG.md [0.3.2] - 2026-07-06` (pyproject bumped 0.3.0 → 0.3.2): geometry-code-aware conformer selection (fixes stochastic BDNN square-plane; BDNN 5/5, full round-trip 25/25); eta RMSD recovery for TiCat1–4 (coord-sphere element-key drop + robust RMSD; TiCat2 Cp aromaticity); wider deduplicated UFF pool; generator energy-None sort crash fix; **eta-ligand winding rac/meso diastereomer fix + generalization to N rings** (TiCat3/4 swap; TiCat5/6 fixtures). Committed on `feature/metallogen-3d-generator` (`04ad753`, `5d6260f`, `118b82c`, `7f0e880`), not pushed.
- **v0.3.3** — In `CHANGELOG.md [0.3.3] - 2026-07-06` (pyproject 0.3.2 → 0.3.3): **metallogen is now the default `OIN3DGenerator`/CLI engine** (was legacy Molassembler) with **g-xTB (`optimizer="xtb"`) as the default optimizer** — a subprocess `xtb --gxtb --opt` wrapper that falls back to FF if the binary is absent. MACE (`mace-omol-*`, fails loudly w/o torch+weights) and legacy Molassembler (`engine="legacy"`, reference for Zone-A P enforcement) are opt-in. Added `oin2xyz --engine/--optimizer`, `tools/install_gxtb.sh`, `run_optimization_grid.py`. `verify_roundtrip.py --optimizer` default → `xtb`. Legacy-specific real-generation unit tests pinned to `engine="legacy"` to keep the fast suite deterministic + heavy-optimizer-free. History squashed to `10e18a1` (consolidated) + `0e264d7` (removed `spec/`) and **pushed** to `origin/feature/metallogen-3d-generator`; the g-xTB layer landed as a follow-on local commit.
- **v0.3.4** — In `CHANGELOG.md [0.3.4] - 2026-07-07` (pyproject 0.3.3 → 0.3.4): default g-xTB optimizer crash fixed (`ASEOptimizer.optimize()` read energy from a non-existent calculator on the subprocess path); dead frag-based charge/BO paths rerouted through the PuLP solver; **MACE/PyTorch moved to an opt-in `--extra mace`** (default install is lightweight FF + g-xTB, no torch); configurable g-xTB `timeout` (default 300s) with FF fallback; `ruff` + CI added (repo lint-clean); packaging metadata + `docs/OPTIMIZERS.md`; repo clutter removed.
- **v0.3.5** — In `CHANGELOG.md [0.3.5] - 2026-07-08` (pyproject 0.3.4 → 0.3.5): **encoder Track A1 — canonical symmetric-donor binding slot** (carboxylate two-oxygen `{slot}` ambiguity, `3b9b6ae`); **η³-allyl double-bond-loss fix** (`_flatten_template` now clears radicals so `build_contract_mol`'s connectivity match transfers bond orders; 5 dataset allyl cases ABAZEK/ABETIK/ABETOQ/ACALOI/AGOVOK key-match, `86aec45`/`b2270f3`); **Zone-A chiral P donor lone-pair stereo recovery** on OIN→3D→OIN (primes `_OIN_CIPCode_LP` from `rdCIPLabeler`, ACUWUT, `65255a1`); **`--quick` roundtrip crash fix** (`max_attempts` filtered out of `ff_params` before `TMCOptimizer`) + **optimizer rename `xtb` → `g-xtb`** across CLI/harness (both spellings still accepted, `f3cab8f`); dataset roundtrip harness enhancements + `recalculate_oin_smiles.py`/`rebuild_summary.py` (`1123423`). On `feature/metallogen-3d-generator`.
- **v0.2.2 (Direct Parser bugfixes)** — still planned/deferred (5 P0/P1 blockers, see below); superseded numbering-wise by 0.3.0 but the work itself is untouched.

### No Active Sprint
The MetalloGen backend is now the **default** generation engine (v0.3.3 cutover, g-xTB default optimizer) and the stereo diagnostic backlog is fully closed — test suite green with **zero expected failures** (`discover tests/unit` 153 OK skip=3; `discover tests` 55 OK; `verify_xyz_to_oin.py` 27/27; full MACE round-trip 25/25). NOTE: the `spec/` tree (worklog/process/compiled) was **removed from the repo** in `0e264d7`, so prior session-persistent notes no longer live in-repo — rely on git history + this file + the auto-memory. Deferred, non-blocking follow-ups: Zone-A N encoding (needs Option-C out-of-band marker); a real compatible-bite bidentate 3D fixture; DG-path set-based enforcement limitation; extending Zone-A P stereo enforcement to the metallogen engine (currently legacy-Molassembler-only); and the Direct Parser bugfixes below. (The TiCat1/3 `[Ti_TET]`↔`[Ti_TPY]` eta string drift was **closed** in `118b82c`; the eta winding rac/meso diastereomer swap in `7f0e880`.)

## Direct Parser — Deferred to v0.2.2
**Audit doc**: `spec/audit/DirectParser_IntegrationAudit_20260506.md`

Integration is blocked by 5 issues in `src/oinsmiles/generation/oin_parser.py`:
1. **(P0 Blocker)** Fragment rank ↔ atom index mapping missing → "bond to self" errors
2. **(P1)** Polydentate ligand connectivity not handled
3. **(P0 Blocker)** No permutation/isomerism selection (cis vs. trans, fac vs. mer)
4. **(P1)** Eta bond translation to atom indices broken
5. **(P2)** Missing test coverage for direct parser

**Production pipeline (current)**: Uses legacy `OINParser.parse()` + `MolassemblerAdapter.generate()` — all integration tests pass.

## Recent Completions

### v0.3.2 Geometry-Aware Selection + Eta RMSD Recovery (2026-07-06)
- **Geometry-fit-ranked conformer selection** — `classify_coordination_geometry()` / `coordination_geometry_fit()` (`utils/oin_aligner.py`) + `_select_by_geometry()` (`generation/metallogen_adapter.py`). Picks the pool conformer with the tightest fit to the requested geometry template (energy ties), fixing the stochastic BDNN `[Pd_SPL]`→`[Pd_TPY]` distortion. Haptic/η gated out (non-regressive). BDNN 5/5, full MACE round-trip 25/25.
- **Eta RMSD-99x recovery (TiCat1–4)** — coord-sphere element-key drop + `_compute_robust_rmsd` (`tests/integration/rmsd_utils.py`); TiCat2 Cp aromaticity via `RemoveHs(sanitize=False)`.
- **Concurrent generator work (bundled)** — wider deduplicated UFF pool (`uff_pool_size`/`rmsd_threshold`/`energy_threshold`), generator energy-None sort crash fix, `TMCOptimizer` returns `(success, energy)`.
- **Closed follow-up (`118b82c`)** — TiCat1/3 eta `[Ti_TET]`↔`[Ti_TPY]` string drift. Extended geometry-fit selection to haptic ligands: `_coordination_vectors` now reduces hapticity to centroid donors (`_reduce_haptic_positions`, <1.6 Å clustering, only when the group count equals the expected coordination number), so TiCat (14→4) / ferrocene (10→2) / Zeise η²-alkene (5→4) become eligible instead of falling back to lowest-energy. Non-regressive; TiCat1/3 now hold `[Ti_TET]` deterministically (8/8 FF).

### v0.3.0 Stereo Round-Trip Arc (2026-07-02 → 07-04)
- **OIN v3.7 descriptor-free metal token** — fixed a stale `is_metal` bug that leaked RDKit's `@SP1` into the metal token (`xyz2mol.py`).
- **Winding preservation (Phase 1)** — `{n>}`/`{n<}` markers now parse through to `ParsedOIN.winding_by_slot` (`oin/inline.py`, `generation/oin_parser.py`).
- **Haptic-face control (Phase 3)** — winding marker steers the generated ring face; per-ring signed-circulation mirror correction (`generation/molassembler_adapter.py`, `oin/winding.py`).
- **Zone-A P encoding + SPL enforcement (Phase 4 / MiniPRD-C)** — `[P@]`/`[P@@]` on metal-bound P; dummy-metal embed makes both enantiomers reachable on square-planar (`core/chirality.py`, `generation/molassembler_adapter.py`).
- **Eta-ring canonicalization** — multi-substituted η-rings round-trip byte-stably: canonical ring-SMILES fragment order + lowest-`CanonicalRankAtoms` heading atom (`utils/oin_aligner.py`).
- **Bidentate incompatible-bite → DG routing** — DIPAMP-class chelates fall back to distance geometry instead of colliding on the template path.

### v0.2.1 Work (folded into [0.3.0])
- **Direct Parser MiniPRD audits (4 of 5)**: RegexPreprocessor, ASTTokenization, MolassemblerInstantiation, FragmentMapping — all audited, archived. Integration MiniPRD deferred.
- **CLI fix**: `oin2xyz` command updated to access `.xyz` from `GeneratedStructure` return type
- **P/N stereocenters fixtures**: PdCl2-R-BINAP, PdCl2-RR-BDNN, PdCl2-RR-BDPP — all pass round-trip RMSD < 1.0 Å
- **TiCat1/3/4 3D generation fix**: ETKDG + de-aromatization strategy for aromatic η5 ligands. See `docs/ETKDG_AROMATIC_FIX.md`.
- **`_extract_oin_constraints()` rename**: 31 call sites updated, `fragment_to_atom_mapping` added to return tuple

### Toolchain (HACF)
- HACF updated to v0.5.1 (post-v0.5.0 skills: `hyper-contextualize`, `hyper-handoff`, `hyper-grill-docs`)
- Agent instruction files (AGENTS.md, CLAUDE.md, GEMINI.md) framing-banner aligned

## Known Limitations (v0.2.1)
- TiCat1/3/4 round-trip bonding inference fails (SINGLE bonds instead of AROMATIC — de-aromatization trade-off)
- `SMILESToXYZ` in `translator.py` is incomplete (dummy atoms, not real SMILES parsing) — TD-003
- `XYZToSMILES.convert()` defined twice — second shadows first — TD-001

## Key Files for Next Session
- `spec/worklog/NOTES.md` — session-persistent state for the stereo arc; read first
- `src/oinsmiles/generation/oin_parser.py` — Direct Parser implementation (blockers above)
- `spec/audit/DirectParser_IntegrationAudit_20260506.md` — full analysis of 5 blockers
- `spec/compiled/SuperPRD_Stereo*.md` — per-feature SuperPRDs (the monolithic `SuperPRD.md` no longer exists)
- `CHANGELOG.md` `[0.3.0]` — the released-but-unpushed 0.3.0 entries
