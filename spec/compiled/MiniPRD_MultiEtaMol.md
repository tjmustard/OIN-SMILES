# MiniPRD: Multi-Eta Fragment Bonded-Mol Reconstruction (G1 + G3)

**Hypergraph Node ID:** `atom_molassembler_adapter`
**Parent Node:** `atom_oin3d_generator` (module `mod_generation`)
**Workstream:** WS-3 of the η-ligand round-trip recovery effort (`spec/worklog/ROUNDTRIP-eta-recovery-handoff.md` §3, Phase 1)
**Compiled by:** `/hyper-resolve` from `spec/active/Draft_PRD_MiniPRD_D_MultiEtaMol.md` + `spec/active/RedTeam_Report.md`
**Status:** Compiled — ready for `/hyper-execute` (Sonnet execute, per handoff §7)

---

## 0. Resolution Ledger (decisions that hardened the draft)

Red Team raised 10 findings (F1–F10) + missing NFRs. All are resolved:

| Ref | Decision (user-confirmed this session) |
| :-- | :-- |
| **Q1 — Indenyl scope** | **Probe-gated, all three.** A Phase-0 probe MUST run on an **indenyl** fragment (not only Cp) before ANY production edit. If green, TiCat1/3/4 all close here with **per-family acceptance gates**. If the indenyl probe is red, escalate (do not silently ship indenyl). |
| **Q2 — Failure gating** | **Full enforcement.** F1, F2, F5 become **hard unit pins** + a harness `mol_reencode_failed` metric that fails the run. No silent-failure surface may pass as green. |
| **F4** | Add a pin: the Zone-A-P loop (`:2306`, now non-empty via the two `[CH3]` metas) runs over the methyl metas against the combined mol containing the aromatic ansa fragment, returns no mismatch, does not raise. |
| **F6** | Assertion strengthened: the **metal atom** carries exactly **12 DATIVE bonds** (not just 12 anywhere), for all of TiCat1/3/4. |
| **F7** | `sorted()` the heavy-atom sets before building **both** emission positions and the renumber/bijection (kept in lockstep). |
| **F8** | The degrade `logger.warning` states **both** the mol *and* the bridge coordinates are unreliable for a non-SiMe2 bridge. |
| **F9** | A 999 RMSD sentinel on TiCat1/3/4 is a **flagged observation**, never folded into "measured not gated"; assert gen-side sphere count == 12. |
| **F10** | Downgrade the ~12 `print(…, file=sys.stderr)` in the touched function to `logger.debug`; **no new `print()`**. |
| **NFR-a** | Unit baseline frozen at **127** (post-Phase-0; the 124→127 delta is `tests/unit/test_xyz2mol_errors.py` + `tests/unit/test_oin_sanitizer_aromaticity.py`). WS-3 pins add on top; enumerate names. |
| **NFR-b** | Add a unit test asserting the documented **5-tuple** return shape (guards a future second caller). |
| **D-1** (carried) | **Full** aromatic restoration (atom `IsAromatic` + `AROMATIC` bond types) from `mol` (`:650`). Makes WS-2 (`oin_aligner.py:50-58`) a provable no-op. |
| **D-2** (carried) | Soft degrade: return `(positions, symbols, None, decisions, [])` + `logger.warning` (never `warnings.warn`/`OINStereoWarning`, never bare `None`). |
| **D-3** (carried) | Enforced `tests/unit` invariant pins + integration confirmation. **Nothing routes to the Candidate Artifact protocol** (all outputs deterministic: `randomSeed=42`, F7 `sorted()`, self-checking `normalize(OIN1)==normalize(OIN2)`). |

---

## 1. The Confidence Mandate

**Confidence Score: split by ligand family (per Q1).**
- **TiCat1 (Cp): 9/10.** Reconstruction arithmetic empirically verified (27 emission
  atoms = 27 `mol_h` atoms; `mol` at `:650` parses with 10 aromatic atoms/bonds);
  all decisions resolved.
- **TiCat3/4 (indenyl): ~6/10 until the Phase-0 probe runs on an indenyl fragment.**
  Benzo-ring atoms enlarge the bijection, `_bfs_atoms` pulls the fused ring, and the
  SMILES-parse aromatic model may disagree with the geometry-perceived model
  (string-mismatch risk, W-1/W-11). **The probe raises indenyl to ≥9 or escalates.**

**Agent instruction:** Do not begin the production edit until **Task 0 (Phase-0 probe,
including an indenyl fragment)** is green. The probe is the load-bearing de-risk step;
skipping it forfeits the confidence basis for TiCat3/4.

---

## 2. Atomic User Stories

* **US-001 (High):** As the generator, I want `_stitch_multi_eta_fragment` to return a
  real bonded mol for a SiMe2 bis-ring fragment so bond topology is emitted instead of
  `mol=None`.
  *AC:* (1) TiCat1/3/4 `generate(oin).mol is not None`. (2) Reconstructed mol built from
  `mol_h`; `mol.GetNumAtoms() == mol_h.GetNumAtoms()`. (3) Ferrocene (single-eta)
  byte-identical (untouched path).

* **US-002 (High):** As the generator, I want full aromatic state restored so a re-encode
  serializes aromatic rings.
  *AC:* (1) Ring bonds aromatic in `mol` (`:650`) are `AROMATIC` in `gen.mol`.
  (2) Ring atoms carry `IsAromatic=True`. (3) WS-2 re-aromatizer (`oin_aligner.py:50-58`)
  is a **provable no-op** on `gen.mol` (feed `gen.mol` ring bonds through the WS-2 loop;
  assert no bond changes).

* **US-003 (Critical):** As `_assemble_combined_mol`, I want emission-space binding
  indices so dative bonds land on the correct atoms.
  *AC:* (1) Function returns `binding_out_idxs` = `[heavy_atom_map[b] for b in
  all_binding_idxs]` (fragment-local emission indices). (2) Caller passes them to
  `fragment_mol_parts`. (3) **[F6] The metal atom carries exactly 12 DATIVE bonds** for
  TiCat1/3/4 (10 ring-C + 2 methyl-C).

* **US-004 (Critical):** As a safety reviewer, I want a mol produced ONLY when provably
  correct.
  *AC:* (1) Non-SiMe2 bridge / >2 slot groups / coverage mismatch / no `mol_h` →
  `mol=None` (XYZ still emitted), never a wrong mol. (2) Degrade emits `logger.warning`
  (not `warnings.warn`). (3) Bare `None` (full placement abort → DG) is never returned by
  the degrade path. **[F8] (4)** For a non-SiMe2 bridge the warning states the XYZ
  bridge coordinates are also unreliable.

* **US-005 (High):** As a maintainer, I want the TiCat1/3/4 round-trip string to close.
  *AC:* (1) Step-2 re-encode routes through `get_oin_string` (not the `convert()`
  fallback). (2) `normalize(OIN1)==normalize(OIN2)` for TiCat1/3/4. (3) TiCat3/4 no
  longer raise.

* **US-006 (Medium):** As a consumer of `GeneratedStructure.mol`, I want the now-non-None
  mol safe for every reader.
  *AC:* (1) RMSD coord-sphere extraction (`verify_roundtrip.py:333/381`) does not newly
  crash; **gen-side** sphere neighbor count == 12. (2) Zone-A-P loop (`:2306`) stays inert
  (see US-009). (3) **[F9]** RMSD value measured; a 999 sentinel is *flagged*, not gated.

* **US-007 (Critical) [F1]:** As the acceptance gate, I want the `get_oin_string` re-encode
  path proven, not assumed.
  *AC:* (1) A unit pin calls `get_oin_string(gen.mol, coords)` directly on a reconstructed
  **TiCat1 and TiCat3** mol and asserts it **does not raise**. (2) Output contains
  lowercase aromatic ring atoms (no `[cH]-[cH]` single-bonded serialization).
  (3) `verify_roundtrip.py` records `mol_reencode_failed` when it falls back to
  `convert()`; a fallback on TiCat1/3/4 **fails the run**.

* **US-008 (Critical) [F2]:** As the re-encoder, I want the emission bijection proven for
  *coordinates*, not just dative bonds.
  *AC:* (1) `gen.mol.GetNumAtoms() == len(all_pos)`. (2) A known ring carbon's `gen.mol`
  conformer position equals its `all_pos` entry (identity coord mapping holds, since
  `get_oin_string` has no `__origIdx` and falls back to identity, `xyz2mol.py:819-825`).
  (3) `renumber_order` is a **complete permutation** of `range(mol_h.GetNumAtoms())`
  (set-equality guard); any miss → `mol=None` + warning.

* **US-009 (Medium) [F4]:** As a safety reviewer, I want the newly-live Zone-A-P path
  proven inert.
  *AC:* For TiCat1, `combined_mol` is non-None → the loop at `:2306` runs over the two
  `[CH3]` metas; `_verify_zone_a_p` returns no mismatch and does not raise on a combined
  mol containing the reconstructed aromatic ansa fragment.

---

## 3. Implementation Plan (Task List)

**Confined to `atom_molassembler_adapter` (function body + its single caller
`_template_generate`), plus the enforced test files.** The re-encoder and `oin_aligner.py`
WS-2 fix are **read-only** (WS-2 is a committed precondition).

- [ ] **Task 0 — Phase-0 probe (BLOCKING, TASK-30 pattern, includes an indenyl
  fragment).** Standalone script; **no production edit yet.** For a Cp ansa
  (`C[Si](C)(c{0}1[cH]{0}[cH]{0<}[cH]{0}[cH]{0}1)c{1}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1`) **and**
  an indenyl ansa (TiCat3 OIN), determine empirically:
  (a) whether `bond.IsInRing()` holds on the reconstructed mol without an explicit
  ring-perception call (`Chem.FastFindRings` / partial sanitize) — check it as WS-2 sees
  it (post `get_oin_string` fragment rebuild), not only on `gen.mol`;
  (b) the exact methyl-pairing/BFS mechanics for indenyl (does `_bfs_atoms` pull the fused
  benzo ring, and is the total bijection still complete?);
  (c) whether `get_oin_string(reconstructed_mol, coords)` runs to completion **without a
  behavior-changing kekulize/sanitize** (the F1 make-or-break), for **both** Cp and
  indenyl. **Exit criterion:** indenyl passes (a)–(c) → proceed; else escalate.

- [ ] **Task 1 — `sorted()` the heavy sets [F7].** In `_stitch_multi_eta_fragment`, sort
  each `ring_heavy_idxs` set before it drives emission order AND before it drives the
  bijection, so atom order is reproducible and the two structures stay in lockstep. (Apply
  in both the Phase-5 loop and the `heavy_atom_map` loop consistently.)

- [ ] **Task 2 — Reconstruct the bonded mol from `mol_h`.** Copy `mol_h` (`:742`), then
  **[D-1 / G3]** restore aromatic state: for every heavy–heavy bond that is `AROMATIC` in
  `mol` (`:650`), set the corresponding bond `AROMATIC` + set both atoms `IsAromatic=True`
  on the copy (same heavy-atom indices; `AddHs` appended H, so indices align).

- [ ] **Task 3 — Build the TWO distinct index structures (do not conflate):**
  - `binding_out_idxs = [heavy_atom_map[b] for b in all_binding_idxs]` — ring-C only,
    fragment-local emission indices (US-003).
  - `renumber_order` — the **total** `output→mol_h` permutation: ring heavy (sorted) +
    ring H (`output_idx + local_i`, same formula as heavies) + Si (`si_atom_idx`) + the two
    methyls (identify `mol_h`'s two Si-bonded CH3 carbons and their 3 H each; pair 1:1 to
    the emission methyl block C,H,H,H,C,H,H,H — safe because both are pure CH3 on Si).
  - **[F2 / US-008]** Guard: `set(renumber_order) == set(range(mol_h.GetNumAtoms()))` and
    `len(renumber_order) == mol_h.GetNumAtoms()`; any miss → `mol=None` + warning.

- [ ] **Task 4 — Renumber + attach conformer.** `Chem.RenumberAtoms(mol_copy,
  renumber_order)` → emission order. Attach a single conformer from the assembled
  **post-Phase-3** `all_positions` so `frag_mol[i] ↔ all_positions[i]`. (Ring-perception
  call from Task 0(a) applied here if the probe found it necessary.)

- [ ] **Task 5 — Coverage guard [US-004].** Produce the mol ONLY when: exactly 2 slot
  groups; bridge atom has exactly 2 non-ipso substituents, each a C with 3 H and degree-1
  heavy (pure CH3); the total bijection (Task 3) is complete; `etkdg_ok` (real `mol_h`).
  Any miss → degrade.

- [ ] **Task 6 — Degrade path [D-2 / F8].** On any guard miss return
  `(positions, symbols, None, decisions, [])` (never bare `None`) and
  `logger.warning(...)` naming the reason; for a non-SiMe2 bridge, the message states the
  XYZ bridge coordinates are also unreliable (two force-placed generic methyls).

- [ ] **Task 7 — Return the 5-tuple.** `_stitch_multi_eta_fragment` returns
  `(positions, symbols, mol, decisions, binding_out_idxs)`. All early/guard `return None`
  paths (`:647,652,683,778,829,852`) are unchanged (whole-function `None` → caller `return
  None` → DG fallback).

- [ ] **Task 8 — Caller wiring (`_template_generate`, `:2085-2099`).** Unpack the 5-tuple;
  when `frag_mol is not None`, `fragment_mol_parts.append((frag_mol, binding_out_idxs))`;
  else `has_all_mols=False` as today. Do **not** add the multi-eta fragment to
  `eta_frag_ranges` (`:2092-2095` comment stands).

- [ ] **Task 9 — Harden `_assemble_combined_mol` [F5].** Before the `try`, validate
  `len(all_pos) == combined_rw.GetNumAtoms()` (after builds) and each fragment's
  `natoms == positions`; in the `except`, `logger.warning`/`logger.error` naming the
  failure (no more silent `return None`).

- [ ] **Task 10 — Observability [F10].** Downgrade the `[DEBUG]` `print(…, file=sys.stderr)`
  calls in `_stitch_multi_eta_fragment` to `logger.debug`. No new `print()`.

- [ ] **Task 11 — Tests (enforced pins, §5) + harness metric.** Add the unit pins for
  US-001…US-009; add `mol_reencode_failed` to `verify_roundtrip.py` metrics and fail the
  run when TiCat1/3/4 fall back to `convert()`.

---

## 4. The Negative Space (Constraints)

* **DO NOT** return a bonded mol unless the coverage guard AND the total-permutation guard
  pass fully. Fail closed to `mol=None`. **Never emit a wrong mol.**
* **DO NOT** return bare `None` on the degrade path — always
  `(positions, symbols, None, decisions, [])` so the XYZ is still emitted.
* **DO NOT** use `warnings.warn` / `OINStereoWarning` for the degrade signal; use the
  module `logger` only (it cannot trip `-W error::OINStereoWarning` in
  `test_zone_a_p_genenforce.py`).
* **DO NOT** rely on `set` iteration order for atom emission — `sorted()` the heavy sets
  (F7).
* **DO NOT** add new `print()` calls; use `logger.debug`/`logger.warning` (F10).
* **DO NOT** let a `get_oin_string` failure pass silently — the harness must record
  `mol_reencode_failed` and fail the run for TiCat1/3/4 (F1).
* **DO NOT** conflate `binding_out_idxs` (ring-C, from `heavy_atom_map`) with
  `renumber_order` (total N-atom permutation) — they are distinct structures.
* **DO NOT** change winding computation (`signed_circulation`, `_determine_winding`) or the
  Phase-3 haptic-face correction block (`:1069-1130`). Reconstruct topology only; leave
  post-Phase-3 `all_positions` as the geometry source.
* **DO NOT** add the multi-eta fragment to `eta_frag_ranges` (`:2092-2095`).
* **DO NOT** modify the single-eta path, the DG worker, or the encoder side
  (`oin_aligner.py`, `xyz2mol.py` are read-only; WS-2 is a committed precondition).
* **DO NOT** shift any existing golden: the **25 pass-1 encode strings**
  (`verify_xyz_to_oin.py`), `test_regression_stability`, the winding-inertness goldens,
  and Ferrocene's round-trip must stay byte-identical (compare **after** WS-2 is merged, so
  WS-2's diff is not misattributed to WS-3).
* **DO NOT** pursue RMSD < 1.0 as a gate — that is WS-4. RMSD is measured only; a 999 is a
  flagged observation.
* **DO NOT** begin the production edit before Task 0 (Phase-0 probe, indenyl included) is
  green.

---

## 5. Integration Tests & Verification

**Enforced `tests/unit` pins (deterministic):**
* **T-US001:** `OIN3DGenerator.generate(<TiCat1 OIN>).mol is not None`; `mol.GetNumAtoms()
  == mol_h.GetNumAtoms()` (27 for TiCat1). Repeat for TiCat3, TiCat4.
* **T-US002:** every ring bond aromatic in `mol` is `AROMATIC` in `gen.mol`; ring atoms
  `IsAromatic=True`; the WS-2 loop (`oin_aligner.py:50-58`) changes **zero** ring bonds on
  `gen.mol`.
* **T-US003/F6:** the metal atom in `gen.mol` has **exactly 12 DATIVE bonds** (TiCat1/3/4).
* **T-US004/F8:** a non-SiMe2 bridge fixture → `gen.mol is None`, XYZ still emitted,
  `logger.warning` fired naming the reason AND flagging XYZ bridge coords; **no bare
  `None`**.
* **T-US007/F1 (Critical):** `get_oin_string(gen.mol, coords)` on reconstructed **TiCat1
  and TiCat3** does not raise and emits lowercase aromatic ring atoms (no `-` between ring
  atoms).
* **T-US008/F2 (Critical):** `gen.mol.GetNumAtoms() == len(all_pos)`; a known ring-C's
  `gen.mol` conformer position `== all_pos[that_index]`; `renumber_order` is a complete
  permutation.
* **T-US009/F4:** for TiCat1, `_verify_zone_a_p(combined_mol, …)` over the methyl metas
  returns no mismatch and does not raise.
* **T-NFR-b:** `_stitch_multi_eta_fragment(...)` returns a 5-tuple of the documented shape.
* **T-Ferrocene:** single-eta round-trip byte-identical (regression guard).

**Integration confirmation (`verify_roundtrip.py`):**
* **Test 1 (Deterministic):** TiCat1 / TiCat3 / TiCat4 → run to completion, step-2 routes
  through `get_oin_string` (metric `mol_reencode_failed == False`),
  `normalize(OIN1)==normalize(OIN2)` **PASS**, TiCat3/4 no longer raise. RMSD measured
  (not gated); a 999 sentinel is reported as a flagged observation.
* **Test 2 (Non-deterministic routing):** **None.** All outputs are deterministic
  (`randomSeed=42`, F7 `sorted()`, self-checking string identity). **No Candidate Artifact
  routing is triggered** (D-3).

**Suite invariants (green at landing):**
```
uv run python -m unittest discover tests            # 55 OK
uv run python -m unittest discover tests/unit       # 127 (baseline, post-Phase-0) + new pins; skipped=3, xfail=0
uv run python tests/integration/verify_xyz_to_oin.py # 25/25 (pass-1 strings byte-identical)
uv run bash tests/run_verification.sh               # roundtrip: TiCat1/3/4 string PASS (RMSD measured)
```

---

## 6. System Graph Blast Radius (`architecture.yml`)

* **Modified (source):** `atom_molassembler_adapter` — `_stitch_multi_eta_fragment` body,
  its caller `_template_generate`, and `_assemble_combined_mol` (F5 logging/validation).
* **Newly-exercised paths (same node, dead for these complexes today):** the Zone-A-P loop
  (`:2306`, via the two `[CH3]` metas) and the normal-path methyl mol-assembly — both first
  run for TiCat1 once `combined_mol` is non-None (US-009).
* **Behavior shift (not source-modified):** `atom_generated_structure` (`.mol` non-None for
  SiMe2 multi-eta); `atom_oin3d_generator` (surfaces it); `atom_xyz2mol` — `get_oin_string`
  becomes the active, **newly-exercised** re-encode path for TiCat1/3/4 (heavy-atom-only
  fragment rebuild on an aromatic-never-kekulized ring; proven by T-US007).
* **Untouched by construction:** single-eta path, DG worker, encoder nodes
  (`atom_oin_aligner`, `atom_oin_sanitizer`, `atom_cip_assigner`), and
  `generation.molassembler_adapter.haptic_face_correction` (`:1069-1130` read, not changed).

---

## 7. Out-of-Scope (later workstreams)
- **WS-4 (G2):** ring-phase / bridge-aware rotation + unconditional clash gate. RMSD is not
  a DoD gate here.
- **WS-5 (G4):** mixed η/σ single-fragment (TiCat2 CGC).
- **WS-6 (G5):** indenyl benzo-ring rigid-drag geometry (topology only here; benzo
  *geometry* distortion is WS-6).
- **WS-7 (G6):** rotamer-phase encoding (V3.8).
