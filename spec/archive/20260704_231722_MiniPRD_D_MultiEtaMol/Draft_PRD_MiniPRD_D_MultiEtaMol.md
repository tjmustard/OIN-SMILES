# Draft PRD — MiniPRD-D: Multi-Eta Fragment Bonded-Mol Reconstruction (G1 + G3)

## Metadata
- **Project Name**: OIN-SMILES — η-ligand round-trip recovery, Phase 1 · WS-3
- **Version**: 0.1.0 (Draft)
- **Status**: Draft
- **Owner**: Architect Agent (seeded from `spec/worklog/ROUNDTRIP-eta-recovery-handoff.md` §3 WS-3 + §2 G1/G3)
- **Target Node**: `atom_molassembler_adapter` (`src/oinsmiles/generation/molassembler_adapter.py`)
- **Effort/Model guidance**: Fable/Opus design (this doc) → Sonnet execute (per handoff §7)

## 1. Introduction & Goals

### 1.1 Problem Statement
On the OIN→XYZ generation side, `_stitch_multi_eta_fragment`
(`molassembler_adapter.py:591`) **unconditionally returns `mol=None`**
(`:1133`, `# mol=None for now (XYZ is the deliverable)`). Any single fragment
that binds a metal at ≥2 distinct slot directions — i.e. an **ansa /
bis-indenyl metallocene** (TiCat1, TiCat3, TiCat4) — flows through this path.
The `None` propagates to the caller (`:2096-2099`): `has_all_mols` is set
`False`, so `_template_generate` never builds a combined RDKit mol
(`:2298`), and `GeneratedStructure.mol` is `None`.

Downstream, the round-trip harness re-encodes the *generated* structure. When
`gen_result.mol is None` it falls back from the direct bonded-mol re-encode
`get_oin_string(gen_result.mol, …)` to `convert(gen_xyz_path)`
(`verify_roundtrip.py:273-294`), which **re-perceives bonds from the raw,
distorted generated coordinates**. This is the failure surface:

- **TiCat1** → coordinate re-perception yields a shredded/garbage step-2 OIN
  string (wrong metal token, mangled fragments) → normalized string check FAILS.
- **TiCat3 / TiCat4** → coordinate re-perception of the indenyl geometry **fails
  outright**, surfacing (post-WS-1) as `ValueError: get_lig_mol failed for
  ligand fragment …` — no step-2 OIN produced at all.

Compounding this is **G3**: `_embed_fragment` de-aromatizes rings for ETKDG
(`:719-728`, aromatic bonds → `SINGLE`, atom aromatic flags cleared) and never
restores them; the fragment's aromatic identity is lost even in the coordinates
that *are* emitted.

**Measured ground truth (2026-07-04, this session — do not re-derive):**
- `OIN3DGenerator.generate("[Ti_TET].C[Si](C)(c{0}1[cH]{0}[cH]{0<}[cH]{0}[cH]{0}1)c{1}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1.[CH3]{2}.[CH3]{3}").mol`
  → **`None`** (G1 reproduced). Generated XYZ has **36 atoms**.
- Parsed fragments: `[Ti]` (rank 0, metal), ansa
  `C[Si](C)(c1[cH][cH][cH][cH]1)c1[cH][cH][cH][cH]1` (frag 1), `[CH3]` (frag 2),
  `[CH3]` (frag 3). `winding_by_slot = {0:'<', 1:None, 2:None, 3:None}`.
- The ansa fragment's `MolFromSmiles(…, sanitize=False)` parse retains **10
  aromatic ring atoms + 10 aromatic bonds** — confirming `mol` (`:650`) is a
  reliable aromatic source for G3.

This is the highest-leverage remaining round-trip fix. Post-Phase-0 the suite
sits at **20/25** round-trip; the 5 residual failures are TiCp2Me2 + TiCat1–4,
all Phase 1+. WS-3 unblocks the topology + string for the ansa/indenyl trio.

### 1.2 Solution Overview
Make `_stitch_multi_eta_fragment` **reconstruct and return a real bonded RDKit
mol** (plus remapped binding indices) whenever it can do so *provably
correctly*, and update its single caller to wire dative bonds at the correct
atom indices. All work is confined to `atom_molassembler_adapter`.

Reconstruction recipe (all inputs already exist at the return site):
1. **Build from `mol_h`** (`:742`) — the ETKDG-embedded fragment already carries
   the complete, correct bond graph (rings, Si–ipso, Si–methyl, C–H) and all 27
   atoms for the SiMe2 bis-Cp case.
2. **Restore full aromatic state (G3, Decision D-1):** copy `IsAromatic` atom
   flags **and** `AROMATIC` bond types from the sanitize-false parse `mol`
   (`:650`) onto the `mol_h` copy. `mol` and `mol_h` share heavy-atom indices
   (both parse the same `frag_smiles`; `AddHs` appends), so the copy is a direct
   per-bond/per-atom transfer.
3. **`RenumberAtoms` into emission order** (ring0 heavy+H, ring1 heavy+H, Si,
   Me-C+3H, Me-C+3H) using the existing `heavy_atom_map` (`:1050-1067`) extended
   to a total `output→mol_h` bijection (ring H atoms follow the same
   `output_idx + local_i` formula already used for heavies; the 2 identical CH3
   groups pair arbitrarily-but-safely).
4. **Attach the conformer** from the already-assembled `all_positions` (emission
   order) so `frag_mol[i] ↔ all_positions[i]`, satisfying
   `_assemble_combined_mol`'s conformer contract (`:582-585`).
5. **Return remapped binding indices** (emission-space) so the caller passes
   *those* to `fragment_mol_parts`, and `_assemble_combined_mol` (`:577-580`)
   wires each metal→binding DATIVE bond at the correct atom.

**Coverage guard + graceful degrade (Decision handoff §2 #3):** produce a mol
**only** when the reconstruction provably covers the whole fragment and the
bridge is exactly SiMe2; otherwise **withhold** the mol (return
`(positions, symbols, None, decisions, [])`, keeping today's XYZ-only behavior)
and emit a `logger.warning`. **Never emit a wrong mol.**

### 1.3 Target Audience
Internal — the OIN→XYZ generator (`atom_molassembler_adapter` →
`atom_oin3d_generator`) and the round-trip re-encoder (`get_oin_string`,
`atom_xyz2mol`). No public API change.

## 2. Confidence Mandate
**Confidence Score**: 8/10.
Root causes are pinned to exact file:lines and empirically reproduced; all three
design decisions are resolved (§5.1); the reconstruction arithmetic is verified
(27 emission atoms = 27 `mol_h` atoms for TiCat1). Residual uncertainty is
execution-local and settled by the acceptance gate, not more human input: (a)
whether `bond.IsInRing()` (relied on by the WS-2 re-aromatizer and by the
guard's ring/methyl classification) needs an explicit ring-perception call
(`Chem.FastFindRings` / partial sanitize) on the reconstructed mol; (b) the
exact methyl-pairing/BFS mechanics for indenyl (TiCat3/4), where `_bfs_atoms`
also pulls the fused benzo ring (topologically fine, part of `mol_h`); (c)
confirming the reconstructed mol flows through `get_oin_string` without a
behavior-changing sanitize. A **Phase-0 probe task** (TASK-30 pattern) de-risks
(a)–(c) before any production edit.

**Clarifying Questions** (all RESOLVED with the user this session — see §5.1):
- [x] G3 scope → **full** aromatic restoration (atom flags + bond types) from `mol`. (D-1)
- [x] Degrade posture → soft-degrade to `mol=None` + `logger.warning` (never `warnings.warn`/`OINStereoWarning`). (D-2)
- [x] Acceptance model → enforced `tests/unit` invariant pins + integration confirmation; **no candidate artifacts**. (D-3)

## 3. Scope

### 3.1 In-Scope
- **G1** — `_stitch_multi_eta_fragment` returns a real bonded mol for
  SiMe2-bridged **bis-Cp** and **bis-indenyl** multi-eta fragments
  (TiCat1 / TiCat3 / TiCat4), plus **remapped emission-space binding indices**.
- **G3** — full aromatic restoration (atom `IsAromatic` + `AROMATIC` bond types)
  on the reconstructed mol, copied from the `mol` parse (`:650`).
- **Caller wiring** — update the multi-eta branch (`:2085-2099`) to unpack the
  new 5-tuple and pass the remapped binding indices to `fragment_mol_parts`.
- **Coverage guard** — exactly 2 slot groups; bridge atom with exactly 2
  non-ipso substituents each a pure `CH3`; total `output→mol_h` bijection;
  `etkdg_ok` (real `mol_h` present). Any miss → withhold mol.
- **Graceful degrade** — soft return (`mol=None`, XYZ preserved) + a
  `logger.warning` naming the reason.
- **Tests (D-3)** — new enforced `tests/unit` pins (mol-not-None, dative-bond
  count, aromatic ring bonds, conformer/atom-count alignment, normalized
  `OIN1==OIN2`) + integration confirmation via `verify_roundtrip.py`.

### 3.2 Out-of-Scope (explicit — later workstreams)
- **WS-4 (MiniPRD-D2, G2):** ring-phase / bridge-aware rotation and the
  unconditional clash gate. **RMSD is NOT a DoD gate here.** DoD is topology +
  string only.
- **WS-5 (MiniPRD-E, G4):** mixed η/σ single-fragment support (TiCat2 CGC).
- **WS-6 (TASK-50, G5):** indenyl benzo-ring rigid-drag geometry.
- **WS-7 (G6):** rotamer-phase encoding (V3.8) for TiCp2Me2/Ferrocene.
- The **single-eta** path (`_stitch_eta_fragment`), the **DG fallback**, and the
  **encoder** side (`oin_aligner.py`, `perception_tmc.py`) — read-only here.
- Any change to winding computation (`signed_circulation`, `_determine_winding`)
  or the Phase-3 haptic-face correction block (`:1069-1130`).

## 4. User Stories (Atomic)
| ID | User Story | Acceptance Criteria | Priority |
| :--- | :--- | :--- | :--- |
| US-001 | As the generator, I want `_stitch_multi_eta_fragment` to return a real bonded mol for a SiMe2 bis-ring fragment so bond topology is emitted instead of `mol=None`. | 1. TiCat1/3/4 `generate(oin).mol is not None`.<br>2. Reconstructed fragment mol built from `mol_h`; atom count == `mol_h` count.<br>3. Ferrocene (single-eta) byte-identical (untouched path). | High |
| US-002 | As the generator, I want full aromatic state restored on the reconstructed mol so a re-encode serializes aromatic rings. | 1. Ring bonds that are aromatic in `mol` (`:650`) are `AROMATIC` in `gen.mol`.<br>2. Ring atoms carry `IsAromatic=True`.<br>3. WS-2 re-aromatizer (`oin_aligner.py:50-58`) is a **no-op** on `gen.mol` (belt-and-suspenders, not a dependency). | High |
| US-003 | As `_assemble_combined_mol`, I want emission-space binding indices so dative bonds land on the correct atoms. | 1. Function returns remapped binding indices (a permutation of `mol_h` via `heavy_atom_map`).<br>2. Caller passes them to `fragment_mol_parts`.<br>3. `gen.mol` has exactly **12** DATIVE bonds on Ti for TiCat1/3/4 (10 ring-C + 2 methyl-C via the normal path). | Critical |
| US-004 | As a safety reviewer, I want a mol produced ONLY when provably correct. | 1. Non-SiMe2 bridge / >2 rings / coverage mismatch / no `mol_h` → `mol=None` (XYZ still emitted), never a wrong mol.<br>2. Degrade emits a `logger.warning` (not `warnings.warn`).<br>3. Bare `None` (full placement abort → DG) is never returned by the degrade path. | Critical |
| US-005 | As a maintainer, I want the TiCat1/3/4 round-trip string to close. | 1. Step-2 re-encode routes through `get_oin_string` (not the `convert()` fallback).<br>2. `normalize(OIN1)==normalize(OIN2)` for TiCat1/3/4.<br>3. TiCat3/4 no longer raise (crash removed via the topology path). | High |
| US-006 | As a consumer of `GeneratedStructure.mol`, I want the now-non-None mol to be safe for every reader. | 1. RMSD coordination-sphere extraction (`verify_roundtrip.py:333/381`, now fed `gen.mol`) does not newly crash; sphere neighbor count == 12.<br>2. Zone-A-P enforcement loop (`:2306`) stays inert (P-less fragments).<br>3. RMSD value is *measured, not gated*. | Medium |

## 5. Technical Specifications (The Blueprint)

### 5.1 Architecture & Resolved Trade-offs

**Data flow (unchanged shape):** `OIN3DGenerator.generate` (`engine.py`) →
`_template_generate` (`:1989`) → per-fragment loop (`:2058`) →
`_stitch_multi_eta_fragment` (`:591`) → `_assemble_combined_mol` (`:556`) →
`GeneratedStructure(xyz, mol, …)`. The re-encoder consumes `.mol` via
`get_oin_string` (`perception_tmc.py:751`). WS-3 changes only what
`_stitch_multi_eta_fragment` returns and how the caller consumes it.

**Current mechanism / root causes (`molassembler_adapter.py`):**
- **G1:** the final `return` (`:1133`) hard-codes the mol slot to `None`; the
  caller sets `has_all_mols=False` (`:2099`); no combined mol is built.
- **G3:** `_embed_fragment` (`:719-728`) sets every aromatic/1.5-order bond to
  `SINGLE` and clears every atom's `IsAromatic` flag before ETKDG, and never
  restores them.
- **Hidden caller contract:** `_assemble_combined_mol` (`:577-580`) adds
  `AddBond(0, frag_start + bidx, DATIVE)` — it indexes **`frag_mol`'s own
  atoms** by `bidx`. The caller currently passes `all_binding_idxs`
  (fragment-SMILES indices, `:2088`), valid only if `frag_mol` is in
  fragment-SMILES order. The reconstructed mol is in **emission** order (a
  permutation), so emission-space binding indices are mandatory.

**Reconstruction correctness facts (verified this session):**
- `mol_h` (27 atoms for TiCat1) already holds the correct bond graph; only
  aromatic typing is lost. `mol` (`:650`) holds the aromatic typing and shares
  heavy-atom indices with `mol_h`.
- `heavy_atom_map` (`:1050-1067`) already maps ring-heavy + Si `mol_h` indices to
  emission indices; ring **H** follow `output_idx + local_i`; the 2 CH3 groups
  are the only atoms NOT in the map (placed geometrically, `:939-1019`) and pair
  1:1 to `mol_h`'s two methyls.
- Binding atoms are ring carbons (heavy) → all present in `heavy_atom_map`, so
  remapped binding indices = `[heavy_atom_map[b] for b in all_binding_idxs]`.

**Resolved Trade-offs Log:**

- **D-1 (G3 scope — the WS-2 subsumption question).** *Issue:* WS-2's
  re-encoder (`oin_aligner.py:50-58`) already restores `AROMATIC` bond type on
  in-ring `SINGLE` bonds — does the reconstruction still need G3? *Analysis:*
  WS-2 fires **only when both endpoint atoms are already `IsAromatic`**, and
  `_embed_fragment` **clears** those atom flags. A verbatim `mol_h` copy would
  therefore serialize the Cp as aliphatic `C1CCCC1` and FAIL the string check —
  so restoring **atom flags** is mandatory regardless of WS-2. *Options:* (A)
  full restoration (atom flags + bond types) from `mol`; (B) atom-flags-only,
  lean on WS-2 for bond types. *Resolution (user-confirmed):* **(A) full
  restoration.** One extra `SetBondType` in the same loop that sets the flag;
  makes `gen.mol` self-consistently aromatic and valid for **every** consumer
  (RMSD sphere, SDF writers, molassembler) independent of the re-encoder, and
  makes WS-2 a provable no-op on multi-eta rather than a hidden dependency. (B)
  leaves `gen.mol` internally inconsistent (aromatic atoms, SINGLE bonds).

- **D-2 (degrade posture).** *Issue:* what happens when the coverage guard
  cannot produce a provably-correct mol. *Options:* (a) silent `mol=None`
  (status quo); (b) `warnings.warn`; (c) `logger.warning`. *Resolution
  (user-confirmed):* **(c) `logger.warning`** naming the reason, plus **soft
  degrade** (return positions + `mol=None`, never bare `None`). `logger.*`
  cannot trip the `-W error::OINStereoWarning` gate in
  `test_zone_a_p_genenforce.py` (the exact coupling class that bit TASK-31/32);
  `warnings.warn` risks it. For TiCat1/3/4 the branch never fires (SiMe2
  succeeds), so it is pure future-observability with zero current-suite effect.

- **D-3 (acceptance model).** *Issue:* how to gate the DoD (topology + string,
  not RMSD) and whether any output is a Candidate Artifact. *Analysis:* topology
  (12 dative bonds, aromatic ring bonds, atom counts) is deterministic given the
  OIN; the OIN(2) string acceptance is the self-checking identity
  `normalize(OIN1)==normalize(OIN2)` (OIN1 recomputed each run → no golden
  file). *Resolution (user-confirmed):* **enforced `tests/unit` invariant pins +
  integration confirmation; nothing routes to the Candidate Artifact protocol.**

### 5.2 System Graph Blast Radius
Nodes in `spec/compiled/architecture.yml`:
- **Modified:** `atom_molassembler_adapter` — the ONLY node whose source
  changes (function body + its single caller `_template_generate`).
- **Behavior shifts (not source-modified here, but now exercised for
  multi-eta):**
  - `atom_generated_structure` — `.mol` becomes non-None for SiMe2 multi-eta.
  - `atom_oin3d_generator` — surfaces the populated mol via `GeneratedStructure`.
  - `atom_xyz2mol` — `get_oin_string` becomes the active re-encode path for
    TiCat1/3/4 (was `convert()` coordinate re-perception).
- **Untouched by construction** (gate every new branch on "SiMe2 bis-ring
  multi-eta" + coverage bijection): the single-eta path, the DG worker
  (`_molassembler_worker`), Zone-A-P enforcement (`_verify_zone_a_p`, inert for
  P-less fragments), the encoder nodes (`atom_oin_aligner`, `atom_oin_sanitizer`,
  `atom_cip_assigner`), and `generation.molassembler_adapter.haptic_face_correction`
  (the `:1069-1130` block is read, not changed).

### 5.3 Execution Checklist (MiniPRDs)
- [ ] `spec/compiled/MiniPRD_MultiEtaMol.md` (to be compiled by `/hyper-resolve`)

### 5.4 API Contracts / Schema
Internal to `molassembler_adapter.py` (no public API change):

- **`_stitch_multi_eta_fragment` return type**: `4-tuple → 5-tuple`.
  - Was: `tuple[np.ndarray, list[str], "Chem.Mol | None", list] | None`
    = `(positions, symbols, mol, decisions)`.
  - Now: `tuple[np.ndarray, list[str], "Chem.Mol | None", list, list[int]] | None`
    = `(positions, symbols, mol, decisions, binding_out_idxs)`.
  - `binding_out_idxs` = binding-atom indices in **emission/output** space
    (same length/order as the caller's `all_binding_idxs`), derived via
    `heavy_atom_map`. On the degrade path (`mol is None`) it is `[]` (caller
    ignores it). Plain tuple (matches the sibling `_stitch_eta_fragment` 4-tuple
    / `_stitch_fragment` 3-tuple conventions — no NamedTuple).
  - Early/guard `return None` paths (`:647,652,683,778,829,852` and the
    not-multi-eta returns) are **unchanged** (whole-function `None` → caller
    `return None` → DG fallback).
- **Caller (`_template_generate`, `:2085-2099`)**: unpack the 5-tuple; use
  `binding_out_idxs` in `fragment_mol_parts.append((frag_mol, binding_out_idxs))`
  when `frag_mol is not None`; else `has_all_mols=False` as today.

### 5.5 Dependencies
- RDKit (`Chem.RenumberAtoms`, `RWMol`, bond/atom `SetBondType`/`SetIsAromatic`,
  `Chem.Conformer`, and — if needed — `Chem.FastFindRings`) — already a hard dep.
- numpy / scipy — already deps. No new libraries.

## 6. Negative Constraints (The "Do NOTs")
- **DO NOT** return a bonded mol unless the coverage guard passes fully
  (SiMe2 bridge + total `output→mol_h` bijection + `etkdg_ok`). Fail closed to
  `mol=None`. **Never emit a wrong mol.**
- **DO NOT** return bare `None` on the degrade path — always return
  `(positions, symbols, None, decisions, [])` so the XYZ is still emitted.
- **DO NOT** use `warnings.warn` / `OINStereoWarning` for the degrade signal;
  use the module `logger` only.
- **DO NOT** change winding computation (`signed_circulation`,
  `_determine_winding`) or the Phase-3 haptic-face correction block
  (`:1069-1130`). Reconstruct topology only; leave `all_positions` (post-Phase-3)
  as the geometry source.
- **DO NOT** add the multi-eta fragment to `eta_frag_ranges` (`:2092-2095`
  comment stands — the free-rotor sweep must not touch it).
- **DO NOT** modify the single-eta path, the DG worker, or the encoder side.
- **DO NOT** shift any existing golden: the **25 pass-1 encode strings**
  (`verify_xyz_to_oin.py`), `test_regression_stability`, the winding-inertness
  goldens, and Ferrocene's round-trip must stay byte-identical.
- **DO NOT** pursue RMSD < 1.0 as a gate — that is WS-4. RMSD is measured only.

## 7. Risks & Mitigation
- **R1 — miswired dative bonds (CRITICAL).** If the returned binding indices are
  not remapped to emission space, `_assemble_combined_mol` bonds the metal to the
  wrong atoms → corrupt re-encode string **and** corrupt RMSD coordination sphere
  (silent, since RMSD now trusts `gen.mol`). → **Mitigation:** remap via
  `heavy_atom_map`; assert exactly 12 DATIVE bonds and sphere neighbor-count == 12
  (US-003, US-006).
- **R2 — coverage-guard false positive.** A near-SiMe2 bridge (e.g. one methyl
  replaced) could slip through and produce a mis-mapped mol. → **Mitigation:**
  strict predicate (bridge has exactly 2 substituents beyond the 2 ipso, each a
  C with 3 H and degree-1 heavy) AND a total-bijection check on atom count; any
  miss → `mol=None`.
- **R3 — ring perception / conformer alignment after `RenumberAtoms`.**
  `bond.IsInRing()` (WS-2 + guard classification) needs valid ring info; the
  conformer must align to `all_positions` after renumber. → **Mitigation:**
  Phase-0 probe decides whether `Chem.FastFindRings`/partial sanitize is needed;
  acceptance asserts `mol.GetNumConformers()==1` and atom-count == len(all_positions).
- **R4 — every `GeneratedStructure.mol` consumer now non-None for multi-eta.**
  RMSD sphere extraction and the Zone-A-P loop both read `gen.mol`. → **Mitigation
  / Red-Team focus (§9):** RMSD must-not-crash gate (US-006); Zone-A-P confirmed
  inert (P-less), but red team must verify no `normal_frag_meta` entry is created
  for the ansa fragment.
- **R5 — `mol` (`:650`) not aromatic if `frag_smiles` lacked lowercase `c`.**
  Then the G3 copy restores nothing and the string comes out aliphatic. →
  **Mitigation:** probe confirmed the ansa fragment parses with 10 aromatic
  atoms/bonds; the restore loop is a fail-safe copy (no aromatic bonds in `mol` →
  no change), and the string acceptance would *detect* (not silently pass) any
  regression.
- **R6 — indenyl benzo atoms (TiCat3/4).** `_bfs_atoms` pulls the fused benzo
  ring into the ring set; it must be fully covered by the bijection. →
  **Mitigation:** the bijection/atom-count check covers all ring atoms; benzo
  *geometry* distortion (G5) is explicitly WS-6, out of scope — topology only.

## 8. Success Metrics
- TiCat1/3/4: `generate(oin).mol is not None`, **12 DATIVE bonds** on Ti,
  ring bonds `AROMATIC`, `GetNumConformers()==1`, atom count == len(all_positions).
- TiCat1/3/4: step-2 re-encode routes through `get_oin_string`;
  `normalize(OIN1)==normalize(OIN2)` passes; TiCat3/4 no longer raise.
- `verify_roundtrip.py`: TiCat1/3/4 run to completion, **string check PASS**,
  RMSD *measured but not gated*.
- Ferrocene round-trip stays PASS (single-eta path byte-identical).
- Suites green: `discover tests` → **55 OK**; `discover tests/unit` → **127 + new
  pins, skipped=3, expected failures=0**; `verify_xyz_to_oin.py` → **25/25**.
- The 25 pass-1 encode strings byte-identical (generator-only change, but a
  shared-path diff guard).

## 9. Notes for /hyper-redteam (next phase — seed focus)
Per the WS-3 brief, stress these first:
1. **`gen.mol` consumer sweep.** Every reader now receives non-None for
   multi-eta: **RMSD coord-sphere extraction** (`verify_roundtrip.py:333/381` →
   `rmsd_utils._extract_coordination_sphere(use_bonds=True)`), the **Zone-A-P
   enforcement loop** (`:2306`, expected inert — confirm no `normal_frag_meta`
   entry is ever created for the ansa fragment), and any MOL/SDF writer.
2. **Conformer / `all_pos` index alignment after `RenumberAtoms`** — prove
   `frag_mol[i]` position == `all_positions[i]` for the whole fragment, including
   the geometrically-placed Si + methyls.
3. **Emission-space binding-index correctness** — the R1 hazard; verify the
   `heavy_atom_map` remap against `_assemble_combined_mol`'s `frag_start + bidx`
   contract.
4. **Coverage-guard completeness** — enumerate what must route to `mol=None`
   (non-SiMe2 bridge, >2 rings, ETKDG-analytic-fallback with no `mol_h`, atom
   -count mismatch) and confirm each degrades softly, not to bare `None`.
5. **Ring-info dependency** — does `bond.IsInRing()` hold on the reconstructed
   mol without an explicit ring-perception call (Phase-0 probe item).
