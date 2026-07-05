# HANDOFF: η-ligand round-trip recovery (XYZ→OIN→XYZ→OIN)

> **Purpose of this file.** Self-contained resumption doc for the effort opened
> 2026-07-04: fixing the 6 round-trip failures surfaced by
> `uv run bash tests/run_verification.sh`, extending the test suites, and
> reaching accurate lossless round-trip for ALL transition-metal complexes.
> Read this first, then `spec/worklog/NOTES.md`. This effort is **new** and runs
> alongside the (now-complete, zero-xfail) v3.7/stereo backlog — do not conflate.
>
> **Status at handoff:** DIAGNOSIS COMPLETE + ROADMAP DRAFTED. No code changed
> yet. Both code-trace investigations and the architect roadmap are captured
> below. Next action = review roadmap, then execute WS-0 → WS-6 in order.
>
> **UPDATE 2026-07-04 — PHASE 0 LANDED (WS-0/1/2).** Policy Q1–Q4 resolved
> (Q4=both string+RMSD; **Q1=ENCODE phase / V3.8 bump** — not accept+document;
> Q2=template; Q3=emit-with-warning). `TASK-40/41/42-*.md` materialized +
> executed. Round-trip **19→20/25** — **Ferrocene now PASSES** (E1 fixed).
> WS-2's fix differs from §3's proposal: staged `SanitizeMol` fails (kekulize
> raises on charge-less Cp; `SetAromaticity` leaves SINGLE bond types), so the
> real fix restores `AROMATIC` bond type on intra-ring SINGLE bonds
> (`IsInRing`-guarded). Suites green: 55 / 127 (skip=3, xfail=0) / 25. Not
> committed (staged for review). See NOTES.md 2026-07-04 Log. **Next = Phase 1
> WS-3 (MiniPRD-D).**

---

## 1. The problem (measured 2026-07-04)

Run: `uv run bash tests/run_verification.sh` → artifacts in
`verification_artifacts_20260704_064522/`.

- **XYZ→OIN (one-way): 25/25 PASS** (`summary_integration.json`,
  `verify_xyz_to_oin.py`). The encoder is sound on *real* input geometry.
- **Unified round trip XYZ→OIN→XYZ→OIN: 19/25 PASS, 6 FAIL**
  (`summary_roundtrip.json`, `verify_roundtrip.py`). Every failure is an
  η-bound (haptic / Cp / indenyl) ligand.

| # | Complex | Failing check(s) | RMSD | Notes |
|---|---|---|---|---|
| Ex4  | Ferrocene | ~~string only (`[cH]-[cH]`)~~ **FIXED 2026-07-04 (WS-2)** | **0.9773** (<1.0) | now PASSES round-trip; E1 flipped it as predicted |
| Ex16 | TiCp2Me2 | string (`[cH]-[cH]`) + RMSD | **1.6752** | E1 fixes string; stays red on RMSD |
| Ex19 | TiCat1 (ansa Me2Si(Cp)2TiMe2) | string (step2 garbage `[Ti_PBP]`, shredded frags) + RMSD | **1.6012** | needs G1 (string) + geometry |
| Ex20 | TiCat2 (CGC Me2Si(Cp)(NtBu)TiMe2) | string (kekulized `C1[CH]=[CH]…`) + RMSD | **999** (sphere mismatch) | fell to molassembler DG; needs G4 |
| Ex21 | TiCat3 (Me2Si(indenyl)2TiMe2, TPY) | **crash** — no step2.oin | — | E3 crash; no RMSD baseline yet |
| Ex22 | TiCat4 (Me2Si(indenyl)2TiMe2, TPY) | **crash** — no step2.oin | — | E3 crash; no RMSD baseline yet |

Crash message (TiCat3/4): `xyz2mol failed: cannot unpack non-iterable NoneType object`.

**Round-trip harness semantics** (`tests/integration/verify_roundtrip.py`):
- Step 2 (re-encode) calls `get_oin_string(gen_result.mol)` **directly** on the
  generator's bonded mol (lines 264-279); only if that raises does it fall back
  to `convert(gen_xyz_path)` (line 283) which re-perceives from coordinates.
  This split is why pass-1 and pass-2 strings differ.
- `normalize_oin_for_comparison` (lines 46-92) strips winding markers
  (`{0>}`→`{0}`), metal stereo descriptors, `[OH2]`→`O`, and canonicalizes slot
  numbering — so **only bond-order / aromaticity diffs survive** to fail the
  string check.
- RMSD threshold **< 1.0** PASS (lines 376-384). **999.0** is a sentinel
  (`tests/integration/rmsd_utils.py:84-92`) for coordination-sphere
  extraction/count/element mismatch, not a distance.

---

## 2. Root causes — CONFIRMED by two code traces

### Encoder side (XYZ→OIN)
- **E1 — serialize-without-sanitize.** `OINSanitizer.generate_robust_smiles`
  (`src/oinsmiles/utils/oin_aligner.py:24-75`) never runs
  `SanitizeMol`/`Kekulize`/aromaticity perception; it serializes whatever atom
  flags + bond types are present. The generator's mol has aromatic-flagged Cp
  atoms joined by **SINGLE** bonds → RDKit emits `[cH]-[cH]`. Fragment copy in
  `get_oin_string` (`xyz2mol.py:899,918`) copies flags verbatim. **Primary cause
  of the Ferrocene/TiCp2Me2 string mismatch.**
- **E2 — distortion-fragile hapticity pruning.** `xyz2mol.py:574-599` drops
  ring donors when `|d(M,Ci) − d(M,Cj)| ≥ 0.4 Å`. Slightly tilted/puckered
  *generated* Cp rings lose donors (η5→η3) → wrong charge/aromaticity → kekulized
  ring output (TiCat2) or total perception failure.
- **E3 — bare `None` → unpack crash.** `get_tmc_mol` returns bare `None`
  (`xyz2mol.py:542-543`) when `get_lig_mol` fails; `translator.py:28` does
  `tmc_mol, xyz_coords = get_tmc_mol(...)` → "cannot unpack non-iterable
  NoneType" (TiCat3/4). Needs a real diagnostic error.
- **E4 — distance-cutoff connectivity.** `xyz2mol_local.py:1035-1052`
  (`Rcov+Rcov+0.5`, over-valence trims longest bond) is also distortion-sensitive
  (background contributor, not a direct fix target yet).

### Generator side (OIN→XYZ, `src/oinsmiles/generation/molassembler_adapter.py`)
- **G1 — `mol=None` from multi-eta.** `_stitch_multi_eta_fragment` (line 591;
  any single fragment binding ≥2 slot groups) returns `mol=None`
  unconditionally (line 1133) → `has_all_mols=False` → **no bond topology**
  emitted for ansa/bis-indenyl → step-2 re-perceives bonds from distorted raw
  coordinates. Highest-leverage generator fix.
- **G2 — quality gate bypassed.** Multi-eta fragments are excluded from
  `eta_frag_ranges` (lines 2092-2095), so the only geometry gate (inter-fragment
  clash check `final_min < 1.45 → return None → DG fallback`, lines 2236-2291)
  never runs for TiCat1/3/4 → clashing garbage geometry (<1 Å atom pairs)
  emitted unconditionally and stamped `Template-generated from OIN (…)`.
- **G3 — de-aromatize without re-aromatize.** `_embed_fragment` converts
  aromatic→SINGLE for ETKDG (lines 719-728, documented in
  `docs/ETKDG_AROMATIC_FIX.md`) and never restores it; combined with G1 no valid
  topology survives.
- **G4 — mixed η/σ fragment unsupported.** TiCat2 has Cp ring + N donor in ONE
  fragment. `_stitch_multi_eta_fragment` treats the single N as a "ring", hits
  `ring_radius < 1e-6 → return None` (line 828) → molassembler DG fallback →
  bad geometry → RMSD 999.
- **G5 — indenyl benzo drag.** `_bfs_atoms` (line 686) pulls the fused benzo
  ring along rigidly while the SVD plane/centroid use only the 5 binding
  carbons → distortion (TiCat3/4).
- **G6 — unconstrained rotamer phase.** Ferrocene/TiCp2Me2 use the single-slot
  path (`_stitch_eta_fragment`, line 1136); geometry is decent but the ring's
  rotational phase about the metal→centroid axis is chosen by a clash-avoidance
  sweep (lines 2246-2272), **not** by any encoded signal. OIN encodes no rotamer
  phase → some of this RMSD may be **format-inherent** (see Policy Q1). For
  symmetric rings the winding marker is deliberately a no-op (line 1401-1402).

### Feasibility findings from the architect pass (corrections/sharpenings)
1. **E1 flips only Ferrocene by itself.** TiCp2Me2 (RMSD 1.675) and TiCat1
   (1.601) still fail on RMSD after E1; TiCat3/4 have no RMSD baseline until the
   crash is removed and a real mol is produced.
2. **G1 has a hidden caller contract.** `_assemble_combined_mol`
   (`molassembler_adapter.py:556-588`) wires dative bonds at `frag_start + bidx`;
   the multi-eta caller passes *fragment-SMILES* atom indices. The stitch
   function's output atom order is a permutation of `mol_h`
   (ring0+H, ring1+H, Si, Me-C+3H, Me-C+3H), so returning a real mol **requires
   also returning remapped binding indices** — the existing `heavy_atom_map`
   (lines 1050-1067) already computes most of this map.
3. **Bridge substituents are hardcoded as 2 methyls** (`_place_tetrahedral_methyls`,
   line 939), placed geometrically, not read from `mol_h`. The G1 mol
   reconstruction must map them back to `mol_h`'s bridge substituents (identical
   methyls → arbitrary pairing is safe) and must **degrade gracefully to
   `mol=None` for non-SiMe2 bridges** (never emit a wrong mol).
4. **Sanitizing in the encoder is safe for OIN bookkeeping.** Slot mapping uses
   explicitly-passed `binding_indices_in_ligand` + a post-hoc
   `MolFromSmiles`/substructure match; RDKit atom props (`__origIdx`,
   `_OIN_CIPCode`) survive `SanitizeMol`. The Cp-anion kekulization failure needs
   a **staged fallback**: full sanitize → `SANITIZE_ALL ^ KEKULIZE` +
   `SetAromaticity` → current no-op.

---

## 3. Roadmap (dependency-ordered; adapt, don't treat as final)

Goal state at every landing point: `verify_roundtrip.py` improving toward 25/25
with these suites green — `discover tests` (55 OK), `discover tests/unit`
(124 OK, 3 skips, **0 xfail**), `verify_xyz_to_oin.py` (25/25). Each workstream
ends with a `NOTES.md` log entry + status-table row. MiniPRDs archive to
`spec/archive/`. Scoped `git add` only; **no push** (standing instruction);
metal is always `fragments[0]`.

### Phase 0 — Harness + cheap encoder robustness (de-risks everything)

- **WS-0 · TASK-40 — per-complex round-trip filter + failure artifacts.**
  *TASK-file, Haiku, ~0.25 session.* Add `--only SUBSTRING` to
  `verify_roundtrip.py` (filter `get_examples()`); persist intermediates in the
  per-example `except` path (TiCat3/4 currently leave no step2.oin). Accept:
  `--only Ferrocene` runs 1 example; full run byte-identical to today. Risk: none.

- **WS-1 · TASK-41 — E3 real errors instead of bare `None`.** *TASK-file, Haiku,
  ~0.25 session.* Replace `return None` at `xyz2mol.py:542-543` with a raised
  `ValueError` naming the offending ligand fragment; audit for other bare-None
  returns. `translator.py:27-30` already wraps → no caller change (grep confirms
  only caller). Accept: TiCat3/4 log a descriptive message (still FAIL — flips
  nothing); all suites green. Risk: none (already-failing paths only).

- **WS-2 · TASK-42 — E1 sanitize/aromaticity-normalize before serialize.**
  *TASK-file, **Sonnet**, with a mandatory Phase-0 probe (TASK-30 pattern),
  ~1 session.* In `OINSanitizer.generate_robust_smiles`
  (`oin_aligner.py:24-75`), before the H-locking loop, on the RWMol copy attempt
  staged: (a) `SanitizeMol`; (b) on kekulize failure `SANITIZE_ALL ^
  SANITIZE_KEKULIZE` + `SetAromaticity`; (c) fall back to current no-op. Probe
  decides whether to also normalize bond types on the fragment copy in
  `get_oin_string` (`xyz2mol.py:918`). **Strict acceptance:**
  1. **Diff gate** — pass-1 OIN strings for all 25 complexes byte-identical
     before/after (protects re-pinned haptic goldens
     `test_haptic_face_golden_match`, `test_winding_inertness.py`). Any diff =
     regression → stop, escalate to MiniPRD.
  2. Ferrocene round-trip → PASS; TiCp2Me2 string check passes (stays FAIL on RMSD).
  3. New unit test (§4) green; all suites green.
  Risk: charge-less Cp fragments hit kekulize failure → fallback (b) is
  load-bearing; unit test must pin it.

### Phase 1 — Generator topology (highest leverage for ansa/indenyl trio)

- **WS-3 · MiniPRD-D — G1+G3: `_stitch_multi_eta_fragment` returns a real bonded
  mol.** *Full HACF chain; Fable/Opus design, Sonnet execute, ~2 sessions.*
  Build the return mol from `mol_h` with aromatic bond types restored from `mol`
  (the sanitized parse at line 650 — fixes G3 too); `RenumberAtoms` into the
  emission order (ring0, ring1, Si, Me…) using `heavy_atom_map`; **coverage
  guard** → `mol=None` if atoms don't fully cover or bridge isn't exactly 2 CH3;
  **return remapped binding indices** and update the caller (lines 2085-2097) so
  `_assemble_combined_mol` wires dative bonds correctly. Accept: TiCat1/3/4
  `gen_result.mol is not None` with correct topology; step-2 goes through
  `get_oin_string`; **normalized string identity passes** for TiCat1/3/4 (needs
  WS-2). RMSD may still fail — DoD is topology + string, **not** RMSD. All suites
  green; Ferrocene stays green. Risks: every `GeneratedStructure.mol` consumer
  now gets non-None for multi-eta (rmsd sphere extraction, Zone-A-P loop at
  line 2306 — inert for these P-less complexes but red-team must confirm);
  conformer/`all_pos` index alignment after renumber.

### Phase 2 — Generator geometry quality (RMSD)

- **WS-4 · MiniPRD-D2 — G2: multi-eta ring phase + quality gate.** *MiniPRD;
  Fable/Opus design, Sonnet execute, ~1.5 sessions. Depends on WS-3.* Two ordered
  pieces: (1) **bridge-aware ring phase** — after plane alignment in Phase 5
  (lines 772-864), rotate each ring about its slot axis so its ipso carbon points
  toward the other ring's slot (the ansa bridge physically determines the phase;
  no encoded signal needed); preserve Phase-3 winding bookkeeping (measure
  circulation *after* the rotation). (2) **run the clash gate unconditionally**
  (currently `if eta_frag_ranges:` at 2236) so multi-eta output is never emitted
  ungated. Do NOT add multi-eta to the free-rotor sweep (2246-2272) — the bridge
  supersedes it. Accept: TiCat1 flips fully (string+RMSD); TiCat3/4 RMSD measured
  (ideally <1.0 → else WS-6); full-25 diff shows only the 6 targets change;
  `test_winding_inertness.py` goldens hold. Risk: gate can route a still-clashing
  complex to DG (worse for haptics) → placement lands first; gate DoD includes
  "does not fire on TiCat1/3/4 after piece 1".

- **WS-5 · MiniPRD-E — G4: mixed η/σ fragment support (TiCat2 CGC).** *MiniPRD;
  recommend a short `/hyper-consult-cto` first (SPL-P Option-E precedent);
  Fable/Opus design, Sonnet execute, ~1.5 sessions. Depends on WS-3.* Generalize
  slot groups in `_stitch_multi_eta_fragment`: a group with `len(bidxs)==1` is a
  σ-donor — place it at `slot_unit * _bond_length`, skip plane-SVD/winding, BFS
  its substituents (tBu) rigidly, let the bridge (Si between Cp-ipso and N)
  proceed. Guard σ-groups out of `signed_circulation`/`_eta_ring_is_symmetric`.
  Accept: TiCat2 template path taken (no DG subprocess), RMSD <1.0, string
  identity (kekulized diff gone via bonded-mol re-encode); all suites + 25/25
  encode green. Alternative (curated DG with haptic constraints) should be
  red-teamed but is expected to lose (DG is the proven-bad path here).

### Phase 3 — Residuals + policy

- **WS-6 · TASK-50 — G5 indenyl diagnostic (then conditional fix).** *Diagnostic
  TASK (TASK-30 pattern), Sonnet, ~0.5 session; fix sized after. Depends on
  WS-3+WS-4.* Measure whether the benzo-ring rigid drag actually costs RMSD after
  Phases 1-2 (coordination-sphere RMSD only scores donor atoms, so benzo
  distortion may be invisible to the gate). Output: a written MIX-style
  classification → cheap fix or escalate. TiCat3/4 may already pass by this point.

- **WS-7 (implied) — G6 rotamer-phase policy for TiCp2Me2/Ferrocene.** Resolve
  Policy Q1 (below). If "accept + document", write it into README/format spec and
  relax the RMSD assertion or mark these as known-limit; if "encode phase", that's
  a format-version bump (V3.8) and a much larger MiniPRD.

---

## 4. Test-suite extensions (fold into the relevant workstream)

- **Per-root-cause unit pins**, following the diagnostic-round-trip pattern in
  `tests/unit/test_stereo_roundtrip_diagnostics.py`:
  - E1: a unit test that feeds `generate_robust_smiles` a mol with aromatic-atom
    + SINGLE-bond Cp state and asserts no `-` in output / consistent aromaticity
    (pins the fallback (b) branch).
  - E3: assert `get_tmc_mol` raises a descriptive `ValueError` (not `TypeError`)
    on the TiCat3 generated XYZ.
  - G1: assert `OIN3DGenerator.generate(TiCat1 OIN).mol is not None` with expected
    dative-bond count.
  - G4: assert TiCat2 takes the template path (no DG) and RMSD <1.0.
- **Fast unit-level round-trip** for one small eta complex (Ferrocene) so the
  loop doesn't require the full `run_verification.sh`.
- Make `verify_roundtrip.py --only` (WS-0) the standard iteration entry point.

---

## 5. Open policy questions (human decides)

- **Q1 — Rotamer phase (G6).** Is the ring rotational-phase RMSD for
  Ferrocene/TiCp2Me2 a **format limitation to accept + document**, or should OIN
  grow a phase encoding (format-version bump, large MiniPRD)? Cheapest path:
  accept + document, relax/annotate the RMSD assertion for symmetric-ring free
  rotors. **This gates whether TiCp2Me2 can ever hit RMSD <1.0 without a format
  change.**
- **Q2 — TiCat2 mixed η/σ.** Template support (WS-5 recommended) vs a curated
  DG path with haptic constraints? Recommendation: template. Confirm before WS-5.
- **Q3 — Clash-gate posture (WS-4).** When the unconditional gate fires on a
  multi-eta complex, is DG fallback (proven-worse for haptics) acceptable, or
  should it be emit-with-warning? Recommendation: emit-with-warning; gate should
  not fire on TiCat1/3/4 after bridge-aware placement lands.
- **Q4 — Definition of "lossless".** Is byte-identical normalized OIN string the
  bar, or also RMSD <1.0? The 6 failures split differently under each (Ferrocene
  passes RMSD; TiCp2Me2 passes string). Recommend: **both** required for "lossless",
  which makes Q1 load-bearing.

---

## 6. Suite invariants (must stay green at every landing point)

```
uv run python -m unittest discover tests            # 55 OK
uv run python -m unittest discover tests/unit       # 124 OK, skipped=3, xfail=0
uv run python tests/integration/verify_xyz_to_oin.py # 25/25
uv run bash tests/run_verification.sh               # target: roundtrip 25/25
```
Root discovery does NOT recurse into `tests/unit` (no `__init__.py` — see
AGENTS.md); run both. Package manager `uv`. Python ≥3.10. Branch is ~20 commits
ahead of origin, **not pushed** (standing instruction). Watch the uncommitted
ruff-adoption tree caveat in NOTES.md if it's still open; use `git commit
--no-verify` only while the HACF pre-commit hook is red, and say so.

---

## 7. Session model / configuration guidance

- **Design / MiniPRD authoring (WS-3, WS-4, WS-5):** Opus 4.8 or Fable 5, higher
  reasoning effort. Use the HACF chain: `/hyper-architect` → `/hyper-redteam` →
  `/hyper-resolve` → `/hyper-execute` → `/hyper-audit`. Consider
  `/hyper-consult-cto` for the WS-5 template-vs-DG decision (SPL-P Option-E
  precedent).
- **TASK-file execution (WS-0, WS-1, WS-2, WS-6):** Sonnet 5 for WS-2/WS-6
  (needs RDKit judgment); Haiku 4.5 fine for WS-0/WS-1 (mechanical, fully
  specified). Each TASK-NN file must be self-contained (exact file+line edits,
  acceptance commands) per the NOTES.md process.
- **Diagnostics / probes:** any tier; the probe *step* inside WS-2 is the
  load-bearing part — do not skip it.
- Keep one workstream per session/context to avoid adversarial-agent
  cross-contamination (CLAUDE.md mandate).

---

## 8. Where things stand / immediate next action

1. Diagnosis: **DONE** (this doc).
2. Roadmap: **DRAFTED** (this doc, §3). Not yet turned into `TASK-40..42` /
   `TASK-50` files or `MiniPRD-D/D2/E` specs.
3. Code: **UNCHANGED.**

**Next action:** decide Policy Q1–Q4 (at least Q4 and Q1), then materialize
Phase 0 as `spec/worklog/TASK-40/41/42-*.md` and execute WS-0 → WS-2. Phase 1+
graduate to MiniPRDs via the HACF chain when their phase begins.

The exploratory plan-mode file (superseded by this doc) lives at
`~/.claude/plans/i-just-ran-a-magical-lightning.md`.
