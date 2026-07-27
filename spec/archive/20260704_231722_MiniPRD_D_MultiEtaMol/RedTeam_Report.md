# RedTeam_Report — MiniPRD-D: Multi-Eta Fragment Bonded-Mol Reconstruction (G1 + G3)

- **Target PRD:** `spec/active/Draft_PRD_MiniPRD_D_MultiEtaMol.md` (v0.1.0 Draft)
- **Blast-radius source:** `spec/compiled/architecture.yml`
- **Method:** Every `file:line` claim in the Draft PRD was cross-checked against the
  live source (`molassembler_adapter.py`, `perception_tmc.py`, `oin_aligner.py`,
  `verify_roundtrip.py`, `rmsd_utils.py`). Findings below are code-verified, not
  prose-only. Verdicts: **CONFIRMED** = I reproduced the claim in source;
  **HAZARD** = a real gap the PRD under-weights or misses.
- **Scope discipline (CRITICAL RULE 1):** No new product features proposed. All
  findings target technical execution, edge cases, resilience, and missing NFRs of
  the *already-scoped* WS-3 work.

---

## 0. Verification Ledger (what I checked against the code)

| PRD claim | Verdict | Evidence |
| :-- | :-- | :-- |
| `_stitch_multi_eta_fragment` returns `mol=None` unconditionally | **CONFIRMED** | `molassembler_adapter.py:1133` `return …, None, decisions` |
| Caller passes fragment-SMILES indices `all_binding_idxs`; unpacks 4-tuple | **CONFIRMED** | `:2085` 4-tuple unpack, `:2088` `all_binding_idxs`, `:2097` append |
| `_assemble_combined_mol` bonds metal at `frag_start + bidx` (frag-local) | **CONFIRMED** | `:578-580` |
| `_embed_fragment` de-aromatizes and never restores | **CONFIRMED** | `:722-728` |
| `mol` (`:650`) shares heavy-atom indices with `mol_h` (AddHs appends) | **CONFIRMED** | `:650` sanitize-false parse; `:732` `AddHs` (default appends H) |
| `heavy_atom_map` covers ring-heavy + Si, **omits ring H and both methyls** | **CONFIRMED** | `:1055-1066` — `if …GetAtomicNum() != 1` skips H; methyls placed geometrically `:1021,1044-1047` (never entered into the map) |
| WS-2 re-aromatizer fires only when both endpoints already `IsAromatic` | **CONFIRMED** | `oin_aligner.py:50-58` (`GetBeginAtom().GetIsAromatic() and …EndAtom()`) |
| 27 emission atoms = 27 `mol_h` atoms (TiCat1 Cp) | **CONFIRMED (Cp only)** | arithmetic: 2×(5C+4H)+Si+2×(C+3H)=27; indenyl **not** verified |
| Blast-radius node ids exist in `architecture.yml` | **CONFIRMED** | `atom_molassembler_adapter:290`, `atom_generated_structure:334`, `atom_oin3d_generator:244`, `atom_xyz2mol:499`, `haptic_face_correction:265` all present |
| Zone-A-P loop creates no `normal_frag_meta` for the ansa fragment | **CONFIRMED — but incomplete** | multi-eta branch `:2078-2100` never touches `normal_frag_meta`; **however the two Ti-methyls do** via the normal path `:2198`, so the loop at `:2306` now runs (see F4) |
| Step-2 re-encode uses `get_oin_string`, else falls back to `convert()` | **CONFIRMED — fallback is SILENT** | `verify_roundtrip.py:288` then `except …:292` prints only a `Note:` and reverts to `convert()` (see F1) |
| `get_oin_string` maps atoms→coords by `__origIdx`, else identity | **CONFIRMED — reconstructed mol has no `__origIdx` → identity** | `perception_tmc.py:819-825`; identity mapping makes the emission bijection load-bearing for *coordinates* (see F2) |
| Coordination sphere uses bonds for gen, distance for orig | **CONFIRMED** | `rmsd_utils.py:62` (`use_bonds=False` orig), `:77` (`use_bonds=True` gen); DATIVE counts as a neighbor `:195` (see F10) |

**Net:** the PRD's factual spine is accurate and unusually well-pinned. The
weaknesses are not wrong facts — they are **silent-failure surfaces**, an
**empirical base that covers only TiCat1 (Cp) while the DoD gates TiCat3/4
(indenyl)**, and **newly-exercised assembly paths** the blast-radius glossed.

---

## 1. Introduction & Goals — Analysis

### Clarifying Questions
1. **Silent-fallback masking.** The harness catches any `get_oin_string`
   exception and reverts to `convert()` (`verify_roundtrip.py:289-292`), printing
   only `Note: get_oin_string failed …`. If the reconstructed aromatic mol trips
   `get_oin_string`, the string check runs against the *same coordinate
   re-perception that fails today* — i.e. the test fails **exactly as now** but
   the run *looks* like the mol path was taken. **What assertion proves the
   `get_oin_string` path actually executed and succeeded, rather than silently
   falling back?** (US-005.1 has no such enforced check.)
2. **Which `frag_smiles` reaches the function at re-encode?** §1.1's measured
   ground truth uses the ansa SMILES *with lowercase aromatic `c`*. Is that the
   string the generator receives at runtime, or can `parsed_oin.fragments[i]`
   arrive already-kekulized/uppercased for some producers? R5 hinges on lowercase
   `c` surviving into `mol` (`:650`).
3. **Is `all_pos[0]` guaranteed to be the metal?** The whole coordinate/atom
   identity chain (F2) assumes combined-mol atom 0 = metal = `all_pos[0]`.
   Confirmed structurally in `_assemble_combined_mol:569` (metal seeded first),
   but the caller's `all_pos` seeding is outside the read window — please pin it.

### What-If Scenarios
- **W-1 (indenyl aromatic-model divergence).** OIN1 is produced by perceiving
  bonds from *real geometry* (`convert(input_xyz)`), OIN2 by copying aromatic
  state from a *SMILES parse* (`:650`). For plain Cp these two aromatic models
  agree. For **indenyl (TiCat3/4)** the fused benzo system can be perceived
  differently by the two paths → different lowercase/kekulé pattern → the
  `normalize(OIN1)==normalize(OIN2)` identity FAILS even though topology is
  "correct." The DoD lists TiCat3/4 but the empirical proof covers only Cp.
- **W-2 (kekulize on serialize).** The reconstructed ring is an aromatic-flagged,
  never-kekulized, explicit-H neutral C5 ring. If *any* step inside
  `get_oin_string` → serialize calls a strict `SanitizeMol`/`Kekulize` (the
  `:801` comment "Sanitize will adjust Zone A" implies a downstream sanitize
  exists beyond the read window), it raises `KekulizeException` → silent
  `convert()` fallback (W-1's failure, dressed as success).
- **W-3 (non-SiMe2 XYZ still wrong).** For a `SiMePh` / `SiEt2` bridge the guard
  withholds the mol — but `all_positions` was **still** built by force-placing two
  *generic* methyls (`_place_tetrahedral_methyls` always emits 2 C + 6 H,
  `:1021`). The "keep today's XYZ-only behavior" degrade therefore emits an XYZ
  with fabricated bridge substituents. Not a WS-3 regression, but the PRD frames
  the degraded XYZ as acceptable when it is silently wrong.

### Points for Improvement
- Add an **enforced unit pin** that calls `get_oin_string(gen.mol, coords)`
  directly on a reconstructed TiCat1 **and** TiCat3 mol and asserts (a) it does
  **not** raise, (b) the output contains lowercase aromatic ring atoms. This
  converts F1 from "integration hope" to "unit fact."
- Make the harness fallback **loud**: on `get_oin_string` exception, record a
  `mol_reencode_failed=True` metric so a silent revert can't pass review.
- The degrade `logger.warning` (D-2) should state that **both the mol and the
  bridge coordinates** are unreliable for a non-SiMe2 bridge, not just the mol.

---

## 2. Confidence Mandate — Analysis

### Clarifying Questions
1. The **8/10** score rests on "reconstruction arithmetic is verified (27 = 27
   for TiCat1)." **What is the verified atom count and aromatic-atom count for an
   indenyl fragment (TiCat3/4)?** Until that exists, confidence for the
   *indenyl half of the DoD* is materially lower than 8/10.
2. Residual uncertainty (c) — "flows through `get_oin_string` without a
   behavior-changing sanitize" — is listed as a minor probe item. Given F1/W-2,
   **should (c) be re-ranked as the primary make-or-break risk** rather than a
   footnote?

### What-If Scenarios
- **W-4 (probe on the wrong molecule).** The Phase-0 probe (TASK-30 pattern) is
  specified generically. If it is run only on TiCat1 (Cp), it will pass (a)-(c)
  and green-light production — yet TiCat3/4 (the crash cases, the actual
  motivation) remain unproven. The probe passing on Cp is **not** evidence for
  indenyl.

### Points for Improvement
- Split the confidence score: **9/10 for TiCat1 (Cp)**, **~6/10 for TiCat3/4
  (indenyl)** until the probe is run on an indenyl fragment. Make "probe includes
  ≥1 indenyl fragment" an explicit Phase-0 exit criterion.

---

## 3. Scope — Analysis

### Clarifying Questions
1. In-scope lists "remapped emission-space binding indices … derived via
   `heavy_atom_map`." But `heavy_atom_map` **omits ring H and both methyls**
   (`:1055-1066`). The binding indices (ring carbons) are covered, yet the
   **`RenumberAtoms` total permutation is a different, larger structure** the map
   does not yet contain. Is building the *total* `output→mol_h` bijection
   (ring H via `output_idx+local_i`, methyls via pairing) understood to be
   **new code**, not a reuse of `heavy_atom_map`?
2. Out-of-scope explicitly excludes the encoder side (`oin_aligner.py`,
   `perception_tmc.py`) as "read-only." But WS-3's success **depends on** WS-2's
   re-aromatizer already being merged (`oin_aligner.py:50-58`). Is WS-2 a hard
   *precondition* dependency (must be committed/verified first), and is that
   ordering pinned anywhere the executor will see?

### What-If Scenarios
- **W-5 (map conflation off-by-one).** Because the PRD describes one map
  (`heavy_atom_map`) doing two jobs, an executor may build the `RenumberAtoms`
  order by inverting `heavy_atom_map` alone — which is missing 14 atoms (all H) +
  8 methyl atoms for TiCat1. `Chem.RenumberAtoms` requires a **complete**
  permutation of all N atoms; a partial list raises or silently corrupts. This is
  the most likely concrete execution bug.

### Points for Improvement
- Rename the two derived structures in the spec so they can't be conflated:
  `binding_out_idxs` (ring-C only, from `heavy_atom_map`) vs `renumber_order`
  (total N-atom permutation, separately constructed). State that
  `len(renumber_order) == mol_h.GetNumAtoms()` and every index appears once — and
  make that a guard, not a comment.
- Add WS-2-committed as an explicit **precondition** row in §5.3.

---

## 4. User Stories — Analysis

### Clarifying Questions
1. **US-003.3 ("exactly 12 DATIVE bonds")** — `_assemble_combined_mol:579`
   guards `if global_bidx < GetNumAtoms()` and **silently skips** an out-of-range
   index. So a remap off-by-one yields **11** bonds with *no error*. Is the
   "assert 12 DATIVE" pin run for **all three** of TiCat1/3/4, and does it also
   assert the 12 are **on the metal atom** (not just 12 DATIVE anywhere)?
2. **US-002.3 (WS-2 is a no-op on `gen.mol`)** — this is only true if G3 restores
   **bond types to AROMATIC** (then `oin_aligner.py:51` `== SINGLE` is false →
   skip). Confirmed self-consistent with D-1(A). But if a future edit reverts G3
   to atom-flags-only (option B), WS-2 silently becomes load-bearing again. Should
   US-002.3 be an **enforced** assertion (feed `gen.mol` to
   `generate_robust_smiles`, assert byte-identical in/out on ring bonds)?
3. **US-006.1 ("sphere neighbor count == 12")** — the orig sphere is
   distance-derived (`rmsd_utils.py:62`, `use_bonds=False`; `MolFromXYZFile` has
   no bonds). Is the assertion on the **gen** sphere (bond-derived, deterministic
   12) or on **both**? The distance-derived orig sphere can return ≠12 → 999
   sentinel even when gen is a clean 12.

### What-If Scenarios
- **W-6 (US-001.3 Ferrocene "byte-identical").** Ferrocene is single-eta
  (`_stitch_eta_fragment`), a different function — untouched **source**. But
  Ferrocene has only one ligand fragment + no methyls; TiCat1 additionally routes
  its two `[CH3]` through the *normal* path whose combined-mol contribution is
  **newly exercised** (F4). "Ferrocene byte-identical" does not cover that new
  interaction — a separate TiCat1-specific assembly pin is needed.

### Points for Improvement
- Strengthen US-003.3 to: "the **metal atom** has exactly 12 DATIVE bonds, and
  `combined_mol.GetNumBonds()` increased by exactly the expected amount." The
  silent-skip at `:579` otherwise hides under-wiring.
- Add US-00X: "for TiCat1, the two `[CH3]` fragments each contribute a mol whose
  atom count equals its emitted position count" — pin the newly-live normal-path
  alignment (F4), not only the ansa fragment.

---

## 5. Technical Specifications — Analysis

### Clarifying Questions
1. **Determinism of emission order.** `ring_all_idxs = _ring_atoms_with_H(mol_h,
   ring_heavy_idxs)` starts from `list(heavy_idxs)` where `heavy_idxs` is a
   **`set`** (`_bfs_atoms` returns a set, `:696`; `_ring_atoms_with_H:758`
   `result = list(heavy_idxs)`). Emission order → XYZ atom order → `get_oin_string`
   identity coord mapping (F2) all inherit **set iteration order**. For `int`
   indices CPython is effectively stable, but this is an *implicit* contract.
   **Is the atom order guaranteed reproducible across runs/interpreters, or should
   the heavy sets be `sorted()` before emission?**
2. **`mol` vs `mol_h` H divergence at the ipso carbon.** `mol` (`:650`) is the
   aromatic parse (ipso `c` bonded to Si → 0 H). `mol_h` de-aromatizes to SINGLE
   *after* `SANITIZE_PROPERTIES` (`:715`), then `AddHs` (`:732`). Does the ipso
   carbon acquire an unexpected explicit H in `mol_h` (making the count ≠ 27)?
   The "27 = 27" measurement implies not — but this is exactly the fragile,
   version-sensitive perception the Phase-0 probe must lock for **indenyl**, whose
   ring-fusion carbons have a different H pattern.
3. **§5.4 return-type change is unguarded at the type level.** The 4→5-tuple
   change means any *other* caller of `_stitch_multi_eta_fragment` breaks on
   unpack. Grep confirms one caller (`:2078`), but is there a test that would
   catch a second caller added later? (No static typing enforces the arity.)

### What-If Scenarios
- **W-7 (conformer/atom-count cascade).** `_assemble_combined_mol` builds one
  conformer sized to `combined_rw.GetNumAtoms()` and fills it from `all_pos`
  (`:582-584`). The invariant is **per-fragment** `frag_mol.GetNumAtoms() ==
  len(frag_positions)`. If the reconstructed ansa mol has 27 atoms but
  `frag_positions` has 27 in a *different order* (renumber ≠ emission), every
  **downstream** fragment's atoms read the wrong coordinates too — the corruption
  is not local to the ansa fragment.
- **W-8 (`_assemble_combined_mol` swallows the failure).** `:587-588`
  `except Exception: return None`. A reconstruction that **passes the coverage
  guard** but produces a mol `CombineMols`/`Conformer` rejects degrades to
  `combined_mol = None` with **no log** — bypassing the D-2 `logger.warning`
  precisely when a subtle bug exists. The "never a silent failure" intent is
  defeated one layer below the guard.

### Points for Improvement
- **`sorted()` the heavy-atom sets** before building both `ring_pos` emission and
  the renumber order. It is free and removes the entire determinism class (F7,
  W-1's ordering component).
- Add a `logger.warning`/`logger.error` inside `_assemble_combined_mol`'s
  `except` (or validate `len(all_pos) == combined_rw.GetNumAtoms()` and
  per-fragment counts *before* the try) so W-8 cannot fail silently.
- Add an explicit **invariant assertion** in the reconstruction:
  `renumber_order` is a permutation of `range(mol_h.GetNumAtoms())` (set-equality
  check) — fail closed to `mol=None` + warning if not.

---

## 5.2 System Graph Blast Radius — Analysis (`architecture.yml`)

### Clarifying Questions
1. The PRD says `atom_molassembler_adapter` is "the ONLY node whose source
   changes." Confirmed for *source*. But **behaviorally**, enabling `combined_mol`
   for TiCat1 lights up code inside the same node that was **dead for this
   complex**: the Zone-A-P enforcement loop (`:2306`, now non-empty
   `normal_frag_meta` from the two methyls) and the normal-path methyl
   mol-assembly. Should the blast radius record these as **newly-exercised paths
   within the node**, not just "inert"?
2. `generation.molassembler_adapter.haptic_face_correction` (node `:265`) is
   declared "read, not changed." But the Phase-3 rotation (`:1095-1118`) mutates
   `all_positions` **after** `heavy_atom_map` is built. Is it confirmed that the
   *conformer attached to the reconstructed mol* uses the **post-rotation**
   `all_positions` (so `frag_mol[i]` ↔ rotated `all_pos[i]`), and that the
   topological map is rotation-invariant? (It is index-based, so yes — but this
   should be an explicit acceptance note, since a future refactor that snapshots
   positions pre-rotation would silently desync geometry from topology.)

### What-If Scenarios
- **W-9 (`atom_xyz2mol` re-encode path swap is under-tested).** The blast radius
  correctly flags that `get_oin_string` becomes the active re-encode path for
  TiCat1/3/4. But `get_oin_string` rebuilds each fragment **heavy-atom-only** with
  H counts derived by counting explicit-H neighbors (`perception_tmc.py:882-893,911-912`)
  and copies bond types verbatim (`:925`). This path has **never** been exercised
  on a mol whose ring bonds are AROMATIC-but-never-kekulized. It is a genuine new
  code path in a node the PRD marks "behavior shift, not source-modified."

### Points for Improvement
- Reclassify `atom_xyz2mol` and the Zone-A-P block from "inert / behavior shift"
  to **"newly-exercised, needs a dedicated pin."** The distinction matters for
  where the executor spends test effort.

---

## 6. Negative Constraints — Analysis

### Clarifying Questions
1. "DO NOT return a bonded mol unless the coverage guard passes fully." The guard
   lives in `_stitch_multi_eta_fragment`, but the **actual mol emission** happens
   two layers up in `_assemble_combined_mol`, which can *independently* return
   `None` (W-8). Does "never emit a wrong mol" also require that a guard-passing
   reconstruction that fails assembly is **logged**, not silently dropped?
2. "DO NOT shift any existing golden … 25 pass-1 encode strings byte-identical."
   WS-3 is generator-only, but it depends on WS-2 (encoder). Is the byte-identical
   check run **after** WS-2 is merged (so the baseline already includes WS-2's
   effect), or against a pre-WS-2 baseline (which would show WS-2's diff and be
   misattributed to WS-3)?

### What-If Scenarios
- **W-10 (observability constraint vs. existing stderr spam).** The degrade
  posture mandates `logger.warning` and forbids `warnings.warn`. Sound. But
  `_stitch_multi_eta_fragment` already emits ~12 raw `print(…, file=sys.stderr)`
  per call (`:623,640,643,751,790,822,928,931,1027,1067`). Mixing new
  `logger.warning` with legacy stderr prints yields inconsistent, un-filterable
  observability — and stderr spam on every generate could mask the one warning
  that matters.

### Points for Improvement
- Add a negative constraint: "the reconstruction/guard path must not add new
  `print()` calls; use the module `logger`." Opportunistically downgrade the
  existing `[DEBUG]` prints to `logger.debug` while the function is open (they are
  already `file=sys.stderr`, so this is behavior-preserving for real output).

---

## 7. Risks & Mitigation — Analysis

### Clarifying Questions
1. **R1 is under-scoped.** The PRD frames miswired binding indices as corrupting
   "the re-encode string and the RMSD coordination sphere." Verified deeper:
   `get_oin_string` maps atoms→coordinates by **identity** for the reconstructed
   mol (no `__origIdx`, `perception_tmc.py:819-825`), so a bad emission bijection also
   feeds **wrong coordinates into PAI alignment** (`:773`) and slot geometry —
   corrupting OIN2 *even if* the dative bonds happened to land right. Should R1's
   blast be expanded to "coordinate mapping," not just "dative bonds + sphere"?
2. **R3 (ring perception after `RenumberAtoms`)** — the mitigation defers to the
   Phase-0 probe. But `bond.IsInRing()` is *also* relied on by the WS-2
   re-aromatizer (`oin_aligner.py:53`) on the re-encoder side. If the reconstructed
   mol reaches WS-2 without valid ring info, WS-2's `IsInRing()` guard mis-fires.
   Is the probe checking `IsInRing()` on the mol **as WS-2 sees it** (post
   `get_oin_string` fragment rebuild), not only on `gen.mol`?

### What-If Scenarios
- **W-11 (R5 fail-safe is not fail-detectable for indenyl).** R5 says "no
  aromatic bonds in `mol` → no change → string acceptance detects it." True for
  Cp (all-or-nothing). For **indenyl**, a *partial* aromatic-perception mismatch
  (5-ring aromatic, benzo perceived differently) yields a mol that is **partly**
  aromatic — plausible enough to serialize, wrong enough to mismatch, and **not**
  obviously detectable as "restored nothing." The fail-safe assumes a binary
  outcome the indenyl case doesn't have.

### Points for Improvement
- Add **R7 — silent assembly degrade** (W-8): `_assemble_combined_mol` returns
  `None` on any exception with no log; mitigate by pre-validating counts and
  logging in the `except`.
- Add **R8 — silent re-encode fallback** (F1): harness reverts to `convert()` on
  `get_oin_string` failure; mitigate with a unit-level direct assertion and a
  harness `mol_reencode_failed` metric.
- Rewrite R1's mitigation to assert, in addition to "12 DATIVE bonds," that a
  **known ring carbon's `gen.mol` coordinate equals its `all_pos` entry**
  (proves the identity coord mapping, the real R1 surface).

---

## 8. Success Metrics — Analysis

### Clarifying Questions
1. "`discover tests/unit` → 127 + new pins, skipped=3, **expected failures=0**."
   The handoff/NOTES baseline is **124** unit OK (`ROUNDTRIP-eta-recovery-handoff.md`
   §6). The PRD asserts **127**. **Which three tests account for 124→127, and are
   they WS-3's new pins or pre-existing drift?** A moving baseline hides
   regressions.
2. The metric "step-2 re-encode routes through `get_oin_string` (not the
   `convert()` fallback)" is **not machine-checkable** in the current harness
   (the fallback is a caught exception + print). How is this metric *measured* at
   acceptance — by grepping stdout for the `Note:` line? That is fragile.

### What-If Scenarios
- **W-12 (RMSD 999 masquerades as "measured not gated").** If the bond-derived
  gen sphere (clean 12) and the distance-derived orig sphere (fuzzy) disagree on
  count/element, `calculate_tmc_rmsd` returns the **999 sentinel**
  (`rmsd_utils.py:92`), not a distance. Because "RMSD is measured, not gated," a
  999 passes WS-3's DoD — but it also means US-006.1's "sphere neighbor count ==
  12" is only half-true (gen side), and a genuine sphere-extraction regression on
  the newly-bonded mol would be **invisible** behind "measured not gated."

### Points for Improvement
- Pin the exact unit count with the enumerated test names (not "127 + new pins"),
  so the executor and auditor share one baseline.
- Make "routes through `get_oin_string`" checkable: assert at unit level that
  `get_oin_string(gen.mol, coords)` returns a string **and** does not raise —
  independent of the integration harness's silent fallback.
- Distinguish, in the success table, **999 (extraction failure)** from a real
  high RMSD; a 999 on TiCat1/3/4 post-WS-3 should be a *flagged observation*, not
  folded into "measured."

---

## 9. Notes for /hyper-redteam (seed items) — Disposition

The PRD's own five seed items are answered here, with verdicts:

1. **`gen.mol` consumer sweep** → **Confirmed + extended.** RMSD sphere
   (`rmsd_utils.py:195`, DATIVE counts as a neighbor — good). Zone-A-P loop
   (`:2306`) — the **ansa** fragment creates no `normal_frag_meta`, but the two
   **`[CH3]` fragments do** (`:2198`), so the loop **runs** for TiCat1 once
   `combined_mol` is non-None. MOL/SDF writers (`verify_roundtrip.py:373-376`) now
   receive the reconstructed mol — confirm `MolToMolFile` tolerates
   aromatic-never-kekulized rings. **See F4.**
2. **Conformer/`all_pos` alignment after `RenumberAtoms`** → **Primary hazard
   (W-7).** The invariant is per-fragment `natoms == len(positions)`, and a
   mismatch corrupts *downstream* fragments too. Pin it.
3. **Emission-space binding-index correctness** → **Confirmed, broader than
   stated (F2).** It also drives the identity coordinate mapping in
   `get_oin_string`, not only dative bonds.
4. **Coverage-guard completeness** → **Confirmed the mol degrades softly**; note
   the **XYZ is still wrong** for non-SiMe2 bridges (W-3), and the *assembly*
   layer can degrade silently (W-8).
5. **Ring-info dependency (`IsInRing`)** → **Confirmed dual-sided:** relied on by
   the reconstruction guard **and** by WS-2 on the re-encoder side
   (`oin_aligner.py:53`). Probe both.

---

## 10. Missing NFRs / Unknown Unknowns

| NFR | Gap | Recommendation |
| :-- | :-- | :-- |
| **Determinism / reproducibility** | Emission (hence XYZ atom) order derives from `list(set(...))` (`:696,758`). Reproducible on CPython-with-int-keys, but an implicit contract. | `sorted()` the heavy sets before emission + bijection. Cheap; eliminates the class. |
| **Failure observability** | `_assemble_combined_mol` (`:587`) and the harness (`:289`) both swallow exceptions and degrade silently, *below* the guard that is supposed to guarantee "never silent." | Log in both `except` blocks; add a `mol_reencode_failed` harness metric. |
| **Logging consistency** | ~12 `print(file=sys.stderr)` in the touched function vs. the new `logger.warning`. | Downgrade `[DEBUG]` prints to `logger.debug`. |
| **Test-baseline integrity** | Success metric shifts unit baseline 124→127 without naming the three deltas. | Enumerate the exact tests; freeze the baseline in the MiniPRD. |
| **Empirical coverage** | All reconstruction proof is on TiCat1 (Cp); TiCat3/4 (indenyl) — the crash cases that motivate WS-3 — are unverified. | Phase-0 probe **must** include an indenyl fragment; split the confidence score. |
| **Acceptance falsifiability** | "routes through `get_oin_string`" and "999 vs high-RMSD" are not machine-checkable today. | Add unit-level direct assertions; treat 999 as a flagged observation. |
| **Return-type arity safety** | 4→5-tuple change is enforced only by the single caller not crashing. | A unit test that asserts the return arity/shape, so a future second caller is caught. |

---

## 11. Prioritized Findings (for /hyper-resolve triage)

| ID | Sev | Finding | Anchor | Fix posture |
| :-- | :-- | :-- | :-- | :-- |
| **F1** | **Critical** | Harness silently reverts to `convert()` on `get_oin_string` failure → a kekulize/sanitize failure on the reconstructed aromatic mol fails the test *exactly as today* but looks like the mol path ran. DoD US-005.1 has no enforced check. | `verify_roundtrip.py:289-292`; `perception_tmc.py:801` | Unit pin: `get_oin_string(gen.mol)` must not raise + emits aromatic; harness `mol_reencode_failed` metric. |
| **F2** | **Critical** | Emission-order bijection feeds `get_oin_string`'s **identity** atom→coord mapping (no `__origIdx`), so a bad permutation corrupts coordinates/PAI/slots, not just dative bonds. Broadens R1. | `perception_tmc.py:819-825,773`; `rmsd_utils.py:200` | Assert a known ring-C's `gen.mol` coord == its `all_pos` entry; total-permutation guard. |
| **F3** | **High** | "Provably correct" is verified for **Cp only**; TiCat3/4 (indenyl) atom count, aromatic-model agreement, and 12-bond count are unproven yet gated by US-001/003/005. | §1.1 measured GT (Cp); `_bfs_atoms:686` pulls benzo | Probe an indenyl fragment; split confidence; consider de-scoping TiCat3/4 to a follow-up if the probe is red. |
| **F4** | **High** | Enabling `combined_mol` for TiCat1 newly exercises the normal-path methyl assembly **and** the Zone-A-P loop (`:2306`, non-empty via the two `[CH3]` metas at `:2198`) — paths dead for this complex today. Blast radius says only "inert." | `:2078-2100` vs `:2197-2207,2306` | Dedicated TiCat1 assembly pin; verify `_verify_zone_a_p` tolerates the aromatic ansa atoms in `combined_mol`. |
| **F5** | **High** | `_assemble_combined_mol` swallows all exceptions → guard-passing-but-assembly-failing reconstruction degrades to `None` **with no log**, defeating D-2 "never silent." | `molassembler_adapter.py:587-588` | Pre-validate counts; log in `except`. |
| **F6** | **Med** | DATIVE under-wiring is silent: `:579` skips out-of-range indices rather than raising → 11 bonds, no error. Makes the "assert 12" pin load-bearing. | `:579` | Strengthen US-003.3 to "metal has exactly 12 DATIVE"; run for all three. |
| **F7** | **Med** | Determinism rests on `list(set(...))` atom ordering flowing into XYZ order + identity coord mapping. | `:696,758` | `sorted()` the heavy sets. |
| **F8** | **Med** | Coverage guard withholds the *mol* for non-SiMe2 bridges but the *XYZ* still contains two force-placed generic methyls — silently wrong coordinates. | `:1021,1044-1047` | Warning must flag the XYZ as unreliable, not only the mol. |
| **F9** | **Med** | RMSD sphere asymmetry (bond-based gen vs distance-based orig) can hit the 999 sentinel that "measured not gated" then hides. | `rmsd_utils.py:62,77,92` | Treat 999 as a flagged observation; assert count on both spheres. |
| **F10** | **Low** | Production `print(file=sys.stderr)` spam mixed with the new `logger.warning`. | `:623…1067` | `logger.debug`. |

**Revised overall confidence:** **8/10 for TiCat1 (Cp)** as the PRD claims;
**~6/10 for TiCat3/4 (indenyl)** until the Phase-0 probe is run on an indenyl
fragment and F1/F2 have enforced unit assertions. The design is sound and
exceptionally well-pinned; the residual risk is concentrated in **silent-failure
surfaces** (F1, F2, F5) and an **empirical base narrower than the DoD** (F3).

---

## Final Action
Run **`/hyper-resolve`** to triage these findings into the compiled
`spec/compiled/MiniPRD_MultiEtaMol.md`. Priority order for resolution:
**F1 → F2 → F3 → F4 → F5**, then fold F6–F10 into the acceptance pins and negative
constraints. In particular, `/hyper-resolve` should decide: (a) whether the
Phase-0 probe must include an indenyl fragment before any production edit
(recommended: **yes**), and (b) whether "routes through `get_oin_string`" and the
999-vs-RMSD distinction become **enforced** assertions rather than observations.
