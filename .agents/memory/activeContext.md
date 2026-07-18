
# Active Context

## Purpose
Captures the current state of OIN-SMILES development. Updated after significant task completions. Read first to understand where to pick up.

## Current State (as of 2026-07-18)

**Release lineage now runs through v0.4.2** — see `CHANGELOG.md` for the authoritative per-version
detail (the block below is retained as history; it predates v0.3.6 and is superseded).

- **v0.3.6 / v0.3.7** — tmCAT/tmPHOTO round-trip fix waves; v0.3.7 pushed to `origin/main`
  (`a0a8b513`), FF-path pass ≈89% (`--quick` sampled).
- **v0.4.0** — performance wave (P0–P11), byte-identical geometry; pushed (`5538b722`, tag `v0.4.0`).
- **v0.4.1** — round-trip tooling, real `SMILESToXYZ`, CI/mypy hardening; pushed (`c7edeeb6`, tag `v0.4.1`).
- **v0.4.2 — round-trip ACCURACY wave. MERGED to LOCAL `main` (`e6febd16`), UNPUSHED.** `pyproject`
  bumped to 0.4.2; `v0.4.2` tag deferred to push. Phases S1/S3/S5/S6a/S6b/S7 (+docs).
  - Full tmCAT/tmPHOTO sweep **complete: 25,197 molecules** — 88.4% `--quick` round-trip pass /
    95.8% accuracy-clean (`docs/ACCURACY_v0.4.1.md`, a screening floor, not a headline).
  - Quick-mode A/B vs a matched v0.4.1 control over all 2,917 failures: **+107 code-attributable
    fixes, 0 deterministic regressions** (`docs/ACCURACY_v0.4.2.md`; tool `tools/ab_compare.py`).
  - New tooling: `tools/ab_compare.py` (two-dir A/B), `tools/milestone_report.py` (backlog snapshots).

**Not yet pushed** (user directive: additional docs/testing first). `README.md` is owned by a concurrent
session and intentionally left untouched here.

## Current State (as of 2026-07-09)

### Release Status
- **v0.2.0** — Released 2026-03-07: Molassembler backend, P/N stereocenter encoding, CLI, OIN v3.6
- **v0.3.0** — In `CHANGELOG.md [0.3.0] - 2026-07-04` (pyproject synced): OIN v3.7 descriptor-free metal token; stereo round-trip arc (winding preservation, haptic-face control, Zone-A P `[P@]`/`[P@@]` encoding + square-planar enforcement); eta-ring canonicalization; bidentate incompatible-bite → DG routing; plus the prior v0.2.1 work (CLI fix, P/N fixtures, TiCat ETKDG fix, Direct Parser audits). Not yet pushed/tagged.
- **v0.3.1** — In `CHANGELOG.md [0.3.1] - 2026-07-05`: MACE MLIP optimizer support, FF convergence knobs, Distinct RMSD error codes, and stability fixes for Ir drift, BDPP/BDNN, FeCO5/FeH2CO4. (Note: `pyproject.toml` was NOT bumped for 0.3.1 — reconciled to 0.3.2 below.)
- **v0.3.2** — In `CHANGELOG.md [0.3.2] - 2026-07-06` (pyproject bumped 0.3.0 → 0.3.2): geometry-code-aware conformer selection (fixes stochastic BDNN square-plane; BDNN 5/5, full round-trip 25/25); eta RMSD recovery for TiCat1–4 (coord-sphere element-key drop + robust RMSD; TiCat2 Cp aromaticity); wider deduplicated UFF pool; generator energy-None sort crash fix; **eta-ligand winding rac/meso diastereomer fix + generalization to N rings** (TiCat3/4 swap; TiCat5/6 fixtures). Committed on `feature/metallogen-3d-generator` (`04ad753`, `5d6260f`, `118b82c`, `7f0e880`), not pushed.
- **v0.3.3** — In `CHANGELOG.md [0.3.3] - 2026-07-06` (pyproject 0.3.2 → 0.3.3): **metallogen is now the default `OIN3DGenerator`/CLI engine** (was legacy Molassembler) with **g-xTB (`optimizer="xtb"`) as the default optimizer** — a subprocess `xtb --gxtb --opt` wrapper that falls back to FF if the binary is absent. MACE (`mace-omol-*`, fails loudly w/o torch+weights) and legacy Molassembler (`engine="legacy"`, reference for Zone-A P enforcement) are opt-in. Added `oin2xyz --engine/--optimizer`, `tools/install_gxtb.sh`, `run_optimization_grid.py`. `verify_roundtrip.py --optimizer` default → `xtb`. Legacy-specific real-generation unit tests pinned to `engine="legacy"` to keep the fast suite deterministic + heavy-optimizer-free. History squashed to `10e18a1` (consolidated) + `0e264d7` (removed `spec/`) and **pushed** to `origin/feature/metallogen-3d-generator`; the g-xTB layer landed as a follow-on local commit.
- **v0.3.4** — In `CHANGELOG.md [0.3.4] - 2026-07-07` (pyproject 0.3.3 → 0.3.4): default g-xTB optimizer crash fixed (`ASEOptimizer.optimize()` read energy from a non-existent calculator on the subprocess path); dead frag-based charge/BO paths rerouted through the PuLP solver; **MACE/PyTorch moved to an opt-in `--extra mace`** (default install is lightweight FF + g-xTB, no torch); configurable g-xTB `timeout` (default 300s) with FF fallback; `ruff` + CI added (repo lint-clean); packaging metadata + `docs/OPTIMIZERS.md`; repo clutter removed.
- **v0.3.5** — In `CHANGELOG.md [0.3.5] - 2026-07-08` (pyproject 0.3.4 → 0.3.5): the tmCAT/tmPHOTO round-trip effort. **Structure-level canonical comparator** `oin/compare.py::canonical_roundtrip_key` (metal+geo, sorted RDKit-canonical fragment multiset, winding multiset — collapses notation drift; harness now compares by key, `a333617`); **carbene/dative-amine m-SMILES H-count fix** (bare-C NHC + dative N–H, ACAWOR/ABESAD, `58913bb`); **backbone P/S/Si stereocentre recovery** (`build_contract_mol` stamps `_OIN_CIPCode` on backbone P; Si/S perceive-then-flip, ABOPOY, `494629c`); **C=C E/Z preserved end-to-end** (encoder carries stereo atoms + `/`\`\` dirs, generator enforces on embed for geometrically-free bonds, `fb8505f`/`0079ac7`); **CN-8 square-antiprismatic `SQA` template + quinoid `_dearomatize_stuck_rings`** (AFEPIM; `KeyError(BondType.AROMATIC)` fix, `6c50286`); **encoder Track A1 canonical symmetric-donor `{slot}`** (`3b9b6ae`); **η³-allyl double-bond-loss fix** (radical-clear in `_flatten_template`; ABAZEK/ABETIK/ABETOQ/ACALOI/AGOVOK, `86aec45`/`b2270f3`); **Zone-A chiral P donor lone-pair stereo recovery** (primes `_OIN_CIPCode_LP` from `rdCIPLabeler`, ACUWUT, `65255a1`); **`--quick` crash fix** (`max_attempts` filtered before `TMCOptimizer`) + **optimizer rename `xtb` → `g-xtb`** (both spellings accepted, `f3cab8f`); dataset harness enhancements + `recalculate_oin_smiles.py`/`rebuild_summary.py` (`1123423`). Merged to `main` via squash-PR (follow-on to PR #2 which brought MetalloGen + v0.3.4).
- **v0.2.2 (Direct Parser bugfixes)** — still planned/deferred (5 P0/P1 blockers, see below); superseded numbering-wise by 0.3.0 but the work itself is untouched.

### Active Sprint — v0.3.6 tmCAT/tmPHOTO round-trip fix wave (S1–S6)
Six parallel sessions, each in its own git worktree, each owning a **disjoint set of files** and one defect class from the tmCAT/tmPHOTO baseline. Handoffs live in `spec/handoffs/v0.3.6/S1-S6.md` (**gitignored**, main-checkout only; each has a `▶ START HERE` bootstrap that creates its own worktree). Shared protocol + file-ownership matrix: `spec/handoffs/v0.3.6/README.md`.

| Session | Class | Owns | Status |
|---|---|---|---|
| S1 | donor-H | `metallogen_adapter.py` → `convert_parsed_to_msmiles` | **DONE** — `ac3a689` on `main` |
| S2 | eta-diene | `metallogen_adapter.py` → `_flatten_template`, `build_contract_mol`, `_oin_fragment_templates` | **DONE** — `bbe567e` on `main` |
| S3 | aromatic-perception | `utils/xyz2mol.py`, `generator3d/process.py` | **DONE** — `d96fd03` on `main` |
| S4 | eta-winding | `utils/oin_aligner.py`, `oin/compare.py` | **DONE** — `882cefb` on `main` |
| S5 | metrics | `tests/integration/rmsd_utils.py`, `tools/*` | **DONE** — `2e6e8a2` on `main` |
| S6 | stereo | `generator3d/{ligand,embed,chem,clean_geometry,__init__}.py`, `core/translator.py` | **DONE** — `9d026e2` on `main` |

**Suite after the full wave:** `discover tests/unit` **350 OK** (skip=5 on rdkit 2026.03.3 and 2025.09.3 — both blessed); `verify_xyz_to_oin.py` **27/27** (S6 fixed the VOacac2 E/Z); lint clean. Note the handoffs' quoted baseline "233 OK / 4 skipped" is **stale** — pristine `main` in a fresh `uv sync` worktree was already 245 OK / 5 skipped before S1. Measure your own baseline before claiming a regression.

**S6's result, and the four things its handoff got wrong:**
- **`atom_stereo` (34 rows) contains ZERO P/N stereocentres** — 31 carbon, 2 sulfur, 1 silicon. The handoff's whole mission ("preserve P/N donor @-stereo") did not apply. The Zone-A P work (`65255a1`) already holds.
- **The real hole:** `generator3d` carried no atom chirality at all (`get_ace_mol_from_rd_mol` copies only Z / charge / bond order), and `embed.py` seeded distance geometry from an **unseeded** `random.randint`. Every sp3 centre was an independent coin flip *per run*. Fixed by mirroring v0.3.5's E/Z carry-through for atoms, plus `generate_3d_structures(seed=42)`. Pristine A/B: the `@@` enantiomer came back mirrored 3/3 runs.
- **`EZ_bond_stereo` + VOacac2 are one ENCODER bug**, not a generator bug: chelate rings close through the metal, and **RDKit's SSSR ignores DATIVE bonds**, so a ring-locked C=C looks acyclic and gets a spurious `/` marker. Ring test now runs on a copy with DATIVE→SINGLE.
- **`no_conformers` (58 rows) is NOT a `--quick` artifact.** All 58 ran under `--quick` during a concurrent sweep, but a 6-row serial retry on pristine gave **2 contention flakes, 4 genuine failures** (`FIXYER`, `ZIHGEE`, `DAHXOB`, `IREPAX`). A failing row burns all 250 embed attempts (`ZIHGEE`: 1696s).
- **Verification trap:** `get_tmc_mol(path, 0)` defaults `with_stereo=False`; the `@` tags come from `CIPAssigner.assign_all()`, which only `XYZToSMILES.convert()` calls. Measured that way, 22/34 rows "pass"; through `convert()`, **32/34 fail**. Always verify stereo through `convert()`.
- **Latent bug found, NOT fixed:** `_apply_double_bond_stereo` force-sets a carried bond back to DOUBLE without adjusting formal charges. Narrowing `ligand.py`'s near-donor filter to the (correct) chelate-ring predicate enforces AFECIZ's free dangling imine, PuLP wants that C=N single with a charged N, and `SanitizeMol` then rejects a 4-valent N — **every `ff_clean` raises** and AFECIZ goes 553s → full attempt budget. The widening was reverted; `AFECIZ`/`XIZXAG` still fail that bond. Fix the charge handling, then narrow. See `tests/unit/test_chelate_locked_ez.py`.
- **P1 of the v0.4.0-perf wave is superseded in part.** S6 already threaded the seed, and the default is now `seed=42` (deterministic), not P1's "without `--seed`, behavior unchanged". Banner added to `spec/handoffs/v0.4.0-perf/P1-seed-hygiene.md`.
- **Not measured:** the 39-row dataset round-trip check, and `no_conformers` triage beyond the 6-row sample.

**Before you delete your session branch, tag it** — `git tag -a archive/<slug> ... && git branch -D <branch>`. A squash-merge means your commits are never ancestors of `main`, so `-D` silently discards the granular history (S3 and S4 both lost theirs; S4's was recovered from a dangling commit before `gc` ran). Rescued as `archive/s3-aromatic-perception` and `archive/s4-eta-winding`. Rule now in `AGENTS.md` → **Version Control**, along with the "N commits not in main is expected after a squash" check and the `git stash push` trap.

**The wave is complete (S1–S6), cut as `v0.3.6`, and `main` is now PUSHED and level with `origin/main`.** A rebase after S6 landed **rewrote every SHA** — `e317999`, `9be2d08`, `d96fd03`, `882cefb`, `2e6e8a2` and friends no longer resolve. S6's code lives at `9d026e2`; its granular history is preserved as tag `archive/s6-stereo` (`694a38a`, on an orphaned line, kept alive only by the tag). Remaining wave-end task: **re-run the dataset on current `main` and regenerate `CASE_REGISTRY.md`** — the existing registry measures pre-S1 code and must not be quoted.

**S4's result, and the corrections a later session would otherwise re-derive:**
- **Winding is meaningless iff a fragment automorphism reverses the eta group's cyclic order** — because turning a ring over is a *proper* rotation that leaves the metal alone. Not "all ring atoms in one symmetry class", which the S4 handoff proposed: mesitylene (2 classes) and an arm-substituted Cp* (4 classes) are both orientation-free, and that weaker test would have left them broken.
- **`compare.py` was deliberately NOT touched.** `test_roundtrip_canonical_key.py::test_eta_winding_flip` pins `>` != `<` on a Cp* at the comparator level (its `N`/`[NH]` difference collapses, so it is a pure winding guard). The encoder-side fix keeps it green; the handoff's suggested "belt-and-braces comparator collapse" would break it.
- **`_fragment_mol_for_canonicalization` returns None on a 4-coordinate neutral boron** (BPh4-, and NEFNER's CH2B(ArF)3 arm) because `SanitizeMol` rejects the valence. RC2 then silently fell through to the *geometric* heading atom, which tracks the embedding. That silent None — not the winding sign — was why those headings wandered.
- **RDKit's SMILES writer ignores `SetNumExplicitHs` on an inherited atom** and prefers its own implicit-valence guess, so two graph-identical phenyls serialize differently and an automorphism test built on an edited copy fails *open*. `_orientation_symmetry_graph` therefore rebuilds the graph from scratch and carries the H count in the **isotope** field (H+1), which neither `CanonicalRankAtoms` nor the writer may drop.
- **A surviving winding difference is now a REAL diastereomer difference.** `QOFTOU_comp_0` is the only such case; MetalloGen builds its ansa-bis(indenyl) zirconocene as rac or meso non-deterministically across runs. Routed to S6 via `tools/triage_overrides.json` (the durable mechanism — `CASE_REGISTRY.md` is gitignored derived data, so hand-editing it is ephemeral).
- **`XIVJEU_comp_0` now passes**; the registry's rac/meso reading of it and of QOFTOU came from the stale `d950f2a` reports. Don't trust registry evidence without an A/B against current `main`.
- 55-case S4 sweep on `882cefb`: **44 success, 0 residual winding-notation failures.** The other 10 are pre-existing defects the fix newly *exposes* by getting past the string check (the harness checks string → RMSD → atom count in that order): 6 RMSD sentinels (`996.0`/`999.0` are literal sentinels in `rmsd_utils.py`, not RMSDs — S5), 2 carborane `@`-stereo (wontfix), 1 E/Z direction (S6), 1 generation atom-count (S1).

**Traps for the next session (learned the hard way in S2):**
- **The shared `tmCAT-tmPHOTO_xyz_dataset/20260707-results/CASE_REGISTRY.md` is STALE** — its reports predate S1/S2, so it still lists 79 `eta_diene_localization` rows that now pass. Re-running `classify_failures.py` does **not** fix this (it re-derives from the same old reports); the dataset must be re-run on `bbe567e` first. Don't pick cases from it blindly. Also: `classify_failures.py --output-dir` **writes** to the dir you point it at — never aim it at the shared results.
- **S3's worktree is branched from `ac3a689`** and must rebase onto `bbe567e`. Three of its cases (AYOVUH, SEJPEE, TESFIH) already re-encode byte-identically after S2.
- **Never run two `tools/test_dataset_roundtrip.py` sweeps concurrently** (or alongside the unit suite) — it fabricates `MetalloGen failed to generate any conformers` errors. `status: pending_g-xtb` is *not* a failure. Conformer choice is stochastic (`TBP` vs `SPY`). Always A/B a suspected regression against pristine code before believing it.
- Two zero-byte reports (`GAHKIK_comp_0.json`, `OVEQAI_comp_0.json`) sit in the shared reports dir, left by a killed continuous runner.

NOTE: the `spec/` tree (worklog/process/compiled) was **removed from the repo** in `0e264d7` and `spec/process/` + `spec/handoffs/` are gitignored — session notes do not live in-repo; rely on git history + this file + the auto-memory. Deferred, non-blocking follow-ups: Zone-A N encoding (needs Option-C out-of-band marker); a real compatible-bite bidentate 3D fixture; DG-path set-based enforcement limitation; extending Zone-A P stereo enforcement to the metallogen engine (currently legacy-Molassembler-only); and the Direct Parser bugfixes below. (The TiCat1/3 `[Ti_TET]`↔`[Ti_TPY]` eta string drift was **closed** in `118b82c`; the eta winding rac/meso diastereomer swap in `7f0e880`.)

## Direct Parser — Deferred to v0.2.2
**Audit doc**: `spec/audit/DirectParser_IntegrationAudit_20260506.md`

Integration is blocked by 5 issues in `src/oinsmiles/generation/oin_parser.py`:
1. **(P0 Blocker)** Fragment rank ↔ atom index mapping missing → "bond to self" errors
2. **(P1)** Polydentate ligand connectivity not handled
3. **(P0 Blocker)** No permutation/isomerism selection (cis vs. trans, fac vs. mer)
4. **(P1)** Eta bond translation to atom indices broken
5. **(P2)** Missing test coverage for direct parser

**Production pipeline (current)**: Uses `OINParser.parse()` + `MetalloGenAdapter.generate()` — all integration tests pass. (The legacy Molassembler backend was removed once MetalloGen became the sole engine.)

## Recent Completions

### S5 — RMSD mapping sentinels + hard watchdog (2026-07-09, `2e6e8a2`)
**All three of the handoff's stated hypotheses were wrong.** `996`/`999` are sentinels, yes — but `999` is *"no metal found"*, not "hard failure inside the metric" (the exception sentinel is `995`), and **all 32 `999` rows are Y or Sc**: `Y` (39) and `Sc` (21) were the only two metals missing from the hand-copied `METAL_ATOMIC_NUMBERS` set in `rmsd_utils.py`. `core/constants.py::TRANSITION_METALS_NUM` has exported the correct list all along, with a docstring warning not to duplicate it (TD-005). The handoff's "Y complexes fail 96% — eta-heavy chemistry" was a hardcoded-set omission. There is also **no eta-centroid reduction** anywhere in the module, though the handoff *and the module's own docstring* both claimed it.

The `996`/`997` rows are a **distance-cutoff artifact on the input side only**: `mol_orig` comes from `MolFromXYZFile` (no bonds), so a covalent-radius cutoff stands in for connectivity while the generated side uses real bonds. It errs in **both** directions, so no radius tuning fixes it — `DAPZIF`'s real apical Pd–N at **2.57 Å** falls outside a **2.54 Å** cutoff, while `ROJXIY` admits a non-donor C at 2.19 Å. Fix: select the input sphere to match the *bonded* sphere's composition (k nearest heavy atoms per element), rejecting when the k-th exceeds `r_cov(M) + r_cov(el) + CEILING_TOL`. The **ceiling, not the string gate, is the safety net** — `verify_ir_complexes.py`/`compare_dg_strategies.py` call the metric with no string gate.

**The decisive structural fact:** `canonical_roundtrip_key` runs *before* the RMSD (`test_dataset_roundtrip.py`), so every sentinel row had already proved OIN(in) == OIN(gen). The chemistry was correct; only the metric failed. 9/9 exemplars now pass with real RMSD 0.25–0.75 (DAPZIF, CAWYOR, ROJXIY, RUBNUZ, ABETIK, NEZWEU, YICXIS, BUJPIH, XEMSAK).

**Zero-regression is provable and was measured**, not assumed: where the old cutoff admitted exactly k atoms of an element, the k *nearest* atoms of that element are the same set. 44 passing molecules across 24 metals, each generated **once** and scored by both old and new metric (re-running the harness cannot isolate a metric change — conformer choice is stochastic): **44/44 identical**.

Taxonomy: `calculate_tmc_rmsd_detailed()` → `(rmsd, None)` | `(None, reason)`; `calculate_tmc_rmsd()` stays a float wrapper for the three integration scripts. Harness emits `RMSD mapping failed at <tier>: <reason>`; `classify_failures.py` bins `rmsd_mapping_failed` **and still reads the legacy `>=900` sentinels** (all 2607 historical reports reclassify unchanged). Note it *already* had an `rmsd_sentinel` bin — what was missing was an honest error string from the harness.

**Watchdog — the handoff was wrong about the hang site too.** `--mol-timeout` was a `signal.alarm`, which cannot interrupt native code; but it also only wrapped *generation*, while **`UGUHAH_comp_0` wedges inside `XYZToSMILES.convert()`** (>150 s, measured), which ran in the parent. Both passes now run in a `spawn` subprocess SIGKILLed on expiry, **with the encode inside the child**; a pass-1 child ships its encode back immediately so a later kill still leaves `smiles_1`. Untimed runs stay in-process (no spawn cost). Verified: 45 s cap kills UGUHAH mid-encode, run continues, 99 s total instead of wedging. Guarded by `tests/unit/test_roundtrip_watchdog.py` (fakes, not real children — spawning under `python -m unittest` makes each child re-run discovery).

**Adjudication (no code change).** The 4 `high_rmsd` rows are FF-fallback *sampling* noise at the 1.0 Å gate on 3–4-atom spheres: `AMUKAV` (worst, 1.1686) **passes at 0.9976** once given a 5-conformer tier; `NUPQOG`/`OWAHEB`/`REQCOI` sit at 1.004–1.022. Don't chase them until an `xtb` binary exists. The 4 `geometry_or_fragment_change` rows fail the *string* gate with `rmsd: null` and **never reach the metric**, so the class was rerouted off S5 to `generator-geometry (unassigned)`: `PUVVEJ` now passes, `FAKZAU` flips between tiers (TBP/SPY are near-degenerate at a 5-coordinate centre), while `CUBCAE` (`[Ni_SPL]`→`[Ni_TET]`) and `NODLAW` (`[Ti_TPL]`→`[Ti_TPY]`) are reproducible generator defects.

**Registry corrections a later session would otherwise re-derive:** the handoff's evidence pack mis-assigns three molecules — `SIMNUZ` is `winding_flip` (**S4**), `UGUHAH` is `macrocycle_perception` (**S3**), `YICXEO` is `interrupted` (rerun-only; the real S5 row is `YICXIS`). S5 owned **70** rows (62 sentinel + 4 high_rmsd + 4 geom), not the "41" the handoff quoted. Also: **`UGUHAH_comp_0.xyz` exists in BOTH `photo/` and `cat/`**, and reports are keyed by basename — the two overwrite each other in `individual_reports/`. Pre-existing harness collision, not fixed here.

### S3 — aromatic/charge perception in the encoder (2026-07-09, `d96fd03`)
**Three of the five handoff hypotheses were wrong.** The biggest bucket was not quinoid chemistry: `RWMol.AddBond(u,v,type)` copies the bond TYPE but creates the bond with `IsAromatic=False`, while `AddAtom` **does** copy the atom flag. `build_contract_mol` kekulizes in place (`Chem.Kekulize` keeps flags on), so the fragment rebuild in `get_oin_string` dropped the ring's only aromatic evidence; `OINSanitizer` then upgraded only the SINGLE ring bonds, emitting the unparseable `Cc1c=c(C)…=c(C)c=1` → `RAW:` token in `compare.py` → key never matches. Fix = normalize an aromatic-flagged bond to `BondType.AROMATIC`; copying the flag alone is **not** enough (Kekulé orders perturb `CanonicalRankAtoms` → Ferrocene/TiCat1 winding flipped).

Other corrections: `fix_equivalent_Os` (where all 17 kekulize tracebacks point) rewrites **ZERO** matches — it is merely the first full sanitize; the real cause is `AC2mol` drawing `P=c` at a wrong ligand charge. `get_tmc_mol` **never returns None** (it raises) — the 8 `xyz2mol_none_crash` rows are `lig_checks` calling `len(res_mols)`, which spends the `ResonanceMolSupplier` cursor, then iterating it. Also: blanket `SetFormalCharge(0)` in `get_oin_string` destroyed nitro → unparseable `N(O)=O`; charges must be restored **in pairs** (a lone `[O+]` rewrites bound CO `C{0}#O` → `C{0}#[O+]`, seen on 3/12 controls).

New shared `utils/aromaticity.py` (`stuck_ring_atoms`, `clear_ring_aromaticity`, `dearomatize_stuck_rings` moved out of `process.py`, `kekulize_safe_sanitize`, `OINEncodeError(ValueError)`).

Results: garbled-`c=` 8/8 round-trip (clears the stale ABERIK 2-iminopyridine cohort); nitro 8/8; `none_crash` all 8 encode. All 17 kekulize crashers encode with zero bare tracebacks — **but do not round-trip**: failures move downstream to atom-count/string mismatches because MetalloGen cannot reproduce their `[CH]` radicals and `=P` ylides. Regression shard: 250 random previously-passing molecules, **248 byte-identical, 0 errors** (the 2 diffs are cross-rdkit E/Z canonicalization, reproducible on unmodified `main`).

**`macrocycle_perception` is NOT aromatic perception** → `docs/KNOWN_LIMITATIONS.md`. A porphyrin is a dianion; `get_tmc_mol` sees two pyrrole N at −1 and, sanitizing with the metal attached, marks the 18-π macrocycle aromatic. OIN carries no formal charge, so MetalloGen builds all four N neutral and the contract mol fails to sanitize (`Explicit valence for atom # 7 N, 4`) — no valence model remains for aromaticity to run on. Two encoder-side fixes were written and reverted after proving both run *after* that point. It belongs to the donor-charge layer (S1-adjacent). Forward-encode stability holds (3 identical encodes).

> Cross-session: S3's fix does **not** shrink S2's `eta_diene` or S4's `winding` buckets — sampled rows re-encode byte-identically before and after (they were never garbled). Their defects are genuinely separate.

### S2 — η-alkene/diene bond-order localization (2026-07-08 → 07-09, `bbe567e`)
`build_contract_mol` recovered a generated ligand's bond orders by substructure-matching `_flatten_template(t)` into the generated fragment. **Both sides are heavy-atom, all-single connectivity graphs of equal size, so that match is an automorphism search** and RDKit returns an arbitrary one; the template's bond orders were copied onto whatever edges it picked. COD's flattened 8-ring has |Aut| = 16 and only 4 maps keep the C=C on the metal-bound carbons — hence `[CH2]=[CH2]` on the backbone (GASBIN/PENGAT), a double bond on a methyl (ABIRIO `C{3<}(=[CH3])`), an alkyne migrated onto a para-ethyl (PIJCAO). Explains **78/78** rows of the class.

> The S2 handoff's stated root cause (an allyl-style *match failure* where nothing transfers, à la v0.3.5) was **WRONG** — the match succeeds and lands wrong. Distinguishing evidence is free: nothing-transferred ⇒ no `=` at all; wrong-automorphism ⇒ the *right number* of double bonds in the *wrong places*. Two earlier handoffs also guessed wrong here.

Fix (all inside `metallogen_adapter.py`):
- `_flatten_template(t, donor_slots)` stamps `_oinSlot` on **every** atom (RDKit's `SubstructMatchParameters.atomProperties` treats an *absent* property as a non-match). Colour by **OIN coordination slot, not donor/non-donor** — a porphyrin's four N are all donors, so a binary colour leaves 8 macrocycle rotations legal and re-picking one shifts a slot label (BOQPIG regressed exactly this way, deterministically 3/3).
- `_generated_donor_slots` assigns generated donors to slots **globally** (bitmask DP, each slot taking the template's donor count). A haptic ligand straddles its slot vector over a wide arc, so per-atom nearest-vector mis-groups an η³-allyl terminus on a chelate (FIKXIJ, an allyl-phosphine).
- `_transfer_score` breaks intra-slot ties against the embedded 3D geometry (the true C=C is ~0.1 Å shorter; MetalloGen embedded from an m-SMILES that still had the true bond orders). The **legacy map is preferred** when slot-valid and within `SCORE_TOL` — formal charges and `_CIPCode` stereo ride along on the same map, so this is a *repair*, not a re-pick.
- Fast path for templates with no non-aromatic multiple bond (bond-order-invariant, and exactly the tBu/phenyl ligands with |Aut| up to ~3e4); `MATCH_MAX = 512`; donor-count guard on fragment↔template pairing; unanchored fallback pass so the change can never do worse than what it replaced.

Result: **`eta_diene_localization` 78 → 0**, 72/78 of the bucket pass. Residuals belong to other sessions (AYOVUH→S3/S4, BIKRIX→S6, FIKXIJ→S5, SEJPEE/TESFIH→S3/S5, CISDOZ `N(O)=O` vs `N(=O)O` is a **comparator** bug in `oin/compare.py`→S4 — RDKit canonicalizes both to the same molecule). No perf regression (porphyrin 46.4s vs 46.3s pristine). Guarded by `tests/unit/test_contract_mol_diene_transfer.py` (8 tests); regression floor `test_contract_mol_allyl_transfer.py` stays green.

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
