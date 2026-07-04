# RedTeam Report — Eta-Ring Canonicalization (Fragment Order + Heading Atom)

**Target:** `spec/active/Draft_PRD.md`
**Blast-radius reference:** `spec/compiled/architecture.yml` (node `atom_oin_aligner`)
**Focus (per invocation):** (1) no golden shifts from the shared `base_sort_key`;
(2) the canonicalization must not mask a real reflection — winding sense stays out of the
fragment/heading identity except as a content-identical tiebreak, and the R2 skip stays
meaningful.

**Verdict:** The PRD's design decisions (D-RC1 canonical ring SMILES + winding tiebreak;
D-RC2 lowest canonical atom rank) are sound in principle, but the PRD **misidentifies the
code lever for RC1** and **assumes an unproven safety property for RC2**. Two findings are
correctness-blocking (RT-1, RT-4). Three are hazards the executor must neutralize
(RT-2, RT-3, RT-5). Full analysis below, per PRD section.

---

## §1 Problem Statement / §5.1 Root-Cause Mechanism — Analysis

### Clarifying Questions
- **RT-1 (BLOCKING) — RC1 lever is dead code.** The PRD (§1.2, §5.1, §5.4) says to "replace
  the arrival-order index `i` in `base_sort_key`/`chem_id`." But `base_sort_key`
  (`oin_aligner.py:250`) is stored as `"key"` at `:308` and **never read anywhere**
  (grep: `"key"` is write-only). The only *live* fragment identity is
  `chem_id = (mass, lig["smiles"])` at `:313`, consumed at `:447/:459/:462/:571`. Editing
  `base_sort_key` is a no-op. Does the resolution acknowledge that RC1 must be fixed at the
  fragment-**rank** assignment, not at `base_sort_key`?
- **RC1 is not a homogeneous-sort problem.** `chem_id` only regroups items *within an equal
  `chem_id` bucket* (`_permute_and_serialize:456-494`); the final string orders by
  `final_sorted_view.sort((rank, local_idx))` (`:497`), and `rank` originates as `"rank": i`
  (`:306`). Two **content-distinct** halide rings sit in **different** `chem_id` buckets and
  are therefore **never compared or reordered** by the homogeneous sort. So RC1 cannot be
  fixed by making `chem_id` canonical; it requires reassigning the ranks of the eta
  fragments themselves. Which code path owns that reassignment, and where is it gated?
- **RC2 heading vs. the pinned golden.** The golden puts the `<` marker on the O-bearing
  carbon (`Oc{0<}1...`), a position produced by the *old geometric* `best_idx` rule. Is the
  O-carbon actually the **lowest-`CanonicalRankAtoms`** atom of that ring? If canonical rank
  lands the heading on a *different* atom, the new rule will not reproduce the golden — the
  golden would then need re-pinning, which the PRD forbids (§6). The executor must confirm
  "min canonical rank == golden's marker atom" for BOTH rings before touching anything.

### What-If Scenarios
- **Global re-rank breaks cisplatin (catastrophic golden shift).** If RC1 is implemented as a
  global fragment re-rank by canonical SMILES, cisplatin's `[Cl]{0}.[Cl]{1}.N{2}.N{3}`
  inverts: canonical `N` (ASCII 0x4E) sorts before `[Cl]` (`[` = 0x5B), yielding
  `N{...}.N{...}.[Cl]{...}.[Cl]{...}`. Every monodentate golden with mixed ligand species
  would shift. The re-rank MUST be scoped to permute *only eta fragments among the rank slots
  they already occupy*, leaving all non-eta fragments and all cross-species interleaving
  byte-identical.
- **Two distinct rings occupy non-adjacent ranks.** In a mixed complex (e.g. an eta ring at
  rank 1 and another at rank 4 with monodentates between), a naive "sort the eta fragments"
  could pull rank-4's ring into rank-1's slot and displace the monodentates. The swap must
  preserve the *set* of rank positions the eta fragments hold, permuting only their
  assignment within that set.
- **Metal invariant.** Rank 0 is the metal by load-bearing invariant
  (`[[project_oin_invariants]]`; `generation/engine.py:158`). Any re-rank that can touch
  rank 0 corrupts every downstream consumer.

### Points for Improvement
- Rewrite §1.2/§5.1/§5.4 to name the **rank-assignment** step (the `"rank": i` origin at
  `_reduce_hapticity:306` and the `target_ranks = sorted(frag_groups.keys())` /
  `final_sorted_view` ordering) as the RC1 lever, and explicitly **retire the `base_sort_key`
  reference** (or note it is dead and to be removed).
- Add an explicit invariant: *"the RC1 canonicalization permutes only the rank labels of
  same-mass eta fragments among themselves; the multiset of rank positions is preserved; all
  non-eta fragments keep their exact rank."* Make this a hard test, not prose.

---

## §5.1 D-RC2 (Heading Rule) + §6 (Winding Frozen) — Analysis  *(SAFETY-CRITICAL)*

### Clarifying Questions
- **RT-4 (BLOCKING) — winding character start-invariance is UNPROVEN and possibly false.**
  The PRD's load-bearing safety claim is "changing the heading atom moves WHICH atom carries
  the marker, never the character." But `signed_circulation` (`oin/winding.py:44-51`) computes
  the sign from a **single edge**: `cross(v_star, v_next) · axis`, where
  `v_next = centered[(star_local_idx + 1) % n]` in **constituent/SMILES order**. The character
  therefore depends on the star atom *unless* (a) the ring is planar-convex AND (b)
  `constituent_indices` order equals the ring's cyclic-adjacency order (so every consecutive
  edge circulates with the same sign). `constituent_indices` is `sorted(local_idx)`
  (`:300`) — **sorted by index, with no guarantee it tracks ring adjacency.** If it doesn't,
  moving the heading to the lowest-canonical-rank atom can **flip the emitted `>`/`<`
  character** for identical geometry. That is not a pure relabeling — it silently changes
  stereo semantics and could fabricate or mask a face flip. This must be proven or defended
  before D-RC2 ships.
- Is `Chem.CanonicalRankAtoms` invoked with `breakTies=True` (order-invariant unique ranks)?
  With `breakTies=False`, symmetric-equivalent atoms share a rank and `min` is ambiguous —
  though the gate is "substituted/asymmetric," a partially-symmetric ring (two identical
  substituents) could still tie.

### What-If Scenarios
- **Non-planar / puckered generated ring.** ETKDG can produce a slightly non-planar Cp. If
  `signed_circulation`'s single-edge measure is start-sensitive on a puckered ring, the
  fixture (planar, hand-built) and the generated (puckered) structures could disagree on the
  character *after* the heading moves — reintroducing the very fixture-vs-generated
  non-determinism RC2 is meant to kill, or worse, flipping a marker.
- **Reflection-masking via heading choice.** If the character is start-dependent and the new
  heading rule is applied uniformly to a genuinely reflected ring, the rule could pick a star
  whose single-edge sign matches the *unreflected* reference — masking the reflection. This is
  exactly the failure mode `test_haptic_face_r2_geometric_fallback_never_auto_substituted`
  guards. The R2 skip must remain present and meaningful (it does under the PRD, but only if
  RC2 does not alter the character computation).
- **Tiebreak is safe, but confirm it.** The RC1 winding tiebreak only *orders* two
  content-identical rings; it never normalizes either ring's character (both are still emitted
  from geometry). Reordering cannot mask a reflection because per-ring characters are
  preserved. Confirmed non-issue — but resolution should state this explicitly so no future
  edit "simplifies" the tiebreak into a normalization.

### Points for Improvement
- **Add a start-invariance proof or a guard.** Either (a) prove that for the substituted-eta
  case `signed_circulation` is start-invariant given the actual `constituent_indices` ordering
  (and add a test asserting the character is identical for every choice of star on the fixture
  ring), or (b) if it is NOT start-invariant, treat the winding character as a property of the
  *ring traversal orientation* and recompute it independent of heading — so the heading choice
  provably cannot move the character. The MiniPRD must not ship D-RC2 on an unproven
  invariant.
- Add a **negative/reflection test** (can reuse the R2 machinery, kept as its own test): a
  synthetically reflected halide ring must still emit the flipped character under the new
  heading rule — proving RC2 is a relabeling, not a normalization.
- Pin `CanonicalRankAtoms(..., breakTies=True)` explicitly in §5.4.

---

## §5.4 API Contracts (mol availability & index provenance) — Analysis

### Clarifying Questions
- **RT-2 (HAZARD) — no RDKit mol at the site.** The `lig` dict carries only `smiles` and
  `binding_atoms` (`(global_idx, mass, coords, local_idx)`) — **no mol object**
  (`:188-189`, `:264`). Both D-RC1 (`MolToSmiles`) and D-RC2 (`CanonicalRankAtoms`) need a
  mol. Where does it come from — rebuilt via `Chem.MolFromSmiles(lig["smiles"])`, or threaded
  in from `xyz2mol`? If rebuilt, atom indices in the rebuilt mol will **not** match the
  `local_idx` provenance of `constituent_indices`, so mapping a canonical rank back to the
  correct `local_idx` (the value the serializer uses at `:605-611`) is a real
  index-aliasing bug surface.
- **RT-3 (HAZARD) — `lig["smiles"]` is the unstable input.** The PRD itself calls
  `lig["smiles"]` "heading/traversal-dependent" (§5.1). Building the RC1 canonical signature
  or the RC2 mol *from that same string* risks inheriting the instability unless the signature
  is a fresh `CanonicalRankAtoms`/`MolToSmiles` (which re-canonicalizes and is order-invariant)
  computed from a correctly-reconstructed mol. Confirm the signature is computed with a
  canonicalization that is provably invariant to `lig["smiles"]`'s atom order.

### What-If Scenarios
- **Bond-perception divergence (fixture vs generated).** The RC1 key is `MolToSmiles` of the
  ring. The fixture is hand-built; the generated structure is ETKDG + molassembler. If the two
  perceive ring bonds differently (aromatic `c` vs Kekulé `C=C`, or a differently-inferred
  substituent bond order), the canonical SMILES **differs** and the RC1 key fails to collapse
  the two — RC1 silently does nothing for the exact case it targets. The signature must be
  computed after a normalization identical to what already produces `lig["smiles"]` (so both
  sides pass through the same perception), or must key on a perception-robust invariant.
- **Local-index gaps.** If `MolFromSmiles` drops explicit Hs or reorders, `constituent_indices`
  (fragment-local) will not index the rebuilt mol's atoms. A wrong-but-plausible heading atom
  produces a wrong-but-plausible golden that passes no test yet corrupts real inputs.

### Points for Improvement
- Specify the **exact mol source** and a **round-trip-free index map** from
  `CanonicalRankAtoms` output to `constituent_indices` / `local_idx`. Prefer threading the
  already-built fragment mol (or an atom-map-preserving construction) over
  `MolFromSmiles(lig["smiles"])`.
- Add a **fail-safe**: if the canonical signature or rank cannot be computed / mapped, fall
  back to today's behavior for that fragment (never a silent partial reorder). The PRD's §5.4
  mentions this; make it a tested branch.

---

## §3 Scope / §7 Risks / §8 Success Metrics — Analysis

### Clarifying Questions
- **Completeness of the RC1 key.** For two rings that are **content-identical AND same
  winding**, the D-RC1 key (canonical SMILES + winding tiebreak) still ties → falls back to
  arrival order. Is residual arrival-order non-determinism acceptable here? (Argument: truly
  identical rings with identical faces produce byte-identical fragments regardless of order,
  so the output string is unchanged — but only if *every* per-atom field is identical.
  Confirm and document, or add a final deterministic tiebreak, e.g. lowest constituent
  global index.)
- §8 asserts "every pre-existing golden byte-identical" but the enumerated gate is the unit
  suite. Are the ansa-metallocene / multi-eta goldens (bridged rings) in the gate? Bridged
  multi-eta is the nearest neighbor to the changed code and most likely to regress.

### What-If Scenarios
- **Symmetric-eta regression via the shared path.** RC2 is gated to "substituted/asymmetric,"
  but the gate itself is new logic on the shared heading-selection loop
  (`_permute_and_serialize:508-597`). A mis-scoped gate could reroute ferrocene / plain-Cp
  (currently `SYMMETRIC_LIGANDS`-forced) through the new rank rule and shift the ferrocene
  golden. The `SYMMETRIC_LIGANDS` branch must remain the first-wins path for those SMILES.
- **`n < 3` degenerate default.** `signed_circulation` returns `">"` for `n < 3`. Ensure the
  RC1/RC2 gates never classify a monodentate or 2-atom binder as an "eta ring" and subject it
  to canonical-rank heading (it has no meaningful ring).

### Points for Improvement
- Add explicit hard tests, beyond the golden flip:
  1. **Non-eta order inertness:** cisplatin/transplatin/BDPP/BDNN/BINAP strings byte-identical.
  2. **Symmetric-eta inertness:** plain ferrocene / ansa-metallocene goldens byte-identical.
  3. **RC1 scoped-swap:** a synthetic complex with an eta ring at a non-1 rank keeps all
     non-eta ranks fixed.
  4. **RC2 start-invariance / reflection:** heading move does not change the character;
     a reflected ring still flips it (the R2 property, as a live assertion in addition to the
     kept skip).
- Elevate RT-4 (start-invariance) to a **precondition** of the MiniPRD: no code lands until it
  is proven or the character computation is made provably heading-independent.

---

## Blast-Radius Summary (architecture.yml)
- **Modified:** `atom_oin_aligner` only (confirmed — the fix is internal to `OINDiscreteAligner`).
- **Must stay byte-inert:** `atom_oin_writer` (consumes the canonical string), `atom_oin_sanitizer`,
  `atom_cip_assigner`, `atom_xyz2mol`, `oin.winding.signed_circulation` (FROZEN), entire
  generation side.
- **Sharpest regression risk:** the shared heading-selection loop and the fragment-rank
  ordering — both touch every fragment, so the scoping gates ("substituted/asymmetric eta"
  and "permute only eta ranks among themselves") are the load-bearing safety boundary.

---

## Findings Ledger (for `/hyper-resolve`)
| ID | Severity | Section | Summary |
|----|----------|---------|---------|
| RT-1 | BLOCKING | §5.1/§5.4 | RC1 lever wrong: `base_sort_key` is dead code; fix must reassign eta-fragment **ranks**, not `chem_id`; global re-rank breaks cisplatin. |
| RT-4 | BLOCKING | §5.1/§6 | Winding character start-invariance unproven; `signed_circulation` is single-edge and may flip on heading change → RC2 not provably a pure relabeling. |
| RT-2 | HAZARD | §5.4 | No mol in `lig` dict; canonical rank/SMILES needs a mol + a correct index map to `constituent_indices`. |
| RT-3 | HAZARD | §5.4 | Signature built from the unstable `lig["smiles"]`; must re-canonicalize order-invariantly and survive fixture-vs-generated bond-perception divergence. |
| RT-5 | HARDEN | §3/§7/§8 | Missing hard tests: non-eta inertness, symmetric-eta inertness, RC1 scoped-swap, RC2 start-invariance + live reflection assertion; RC1 key incomplete for identical-content/same-winding rings. |

**Next:** run `/hyper-resolve` to triage RT-1…RT-5 and compile the final SuperPRD + MiniPRD.
