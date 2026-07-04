# Red Team Report — MiniPRD-C: Zone-A P SPL Dummy-Metal Embed

**Target:** `spec/active/MiniPRD_ZoneA_P_SPL_DummyEmbed.md`
**Blast-radius reference:** `spec/compiled/architecture.yml`
**Seed:** `spec/worklog/SPL-P-enforcement-decision.md` (Option E)
**Analyst:** Red Team Agent (`/hyper-redteam`, 2026-07-03)
**Verdict:** The *decision* (Option E) is sound and well-argued. The *MiniPRD* has
**two factual defects that will actively mislead the executor** and **three
unresolved scoping/oracle gaps** that can produce silent, test-passing failures.
Recommend a `/hyper-resolve` pass before `/hyper-execute`.

---

## TOP FINDINGS (ranked)

1. **[BLOCKER — stale reference]** `PseudoAtomStrategy` and `PSEUDO_ATOMIC_NUM`
   **no longer exist in `src/`.** `grep -rn "PseudoAtomStrategy\|PSEUDO_ATOMIC_NUM" src/`
   returns nothing; `core/chirality.py:22` is now `from .constants import
   TRANSITION_METALS_NUM`. `architecture.yml:209` confirms it: *"PseudoAtomStrategy
   scaffolding removed, behaviour retained."* Yet the MiniPRD instructs the executor
   in **three places** (Task 2, Constraint §4, and the §1 context list — *"the
   dead-code device this MiniPRD revives"*) to revive/reuse it at `core/chirality.py:22`.
   The executor will follow a dead pointer.

2. **[BLOCKER — wrong reuse target]** The fallback reuse target, `_build_dummy_metal_copy`
   (`chirality.py:82`), **cannot be reused as-is for the attach.** It *converts an
   existing metal neighbour into a dummy* — its first act is `if metal_idx is None:
   raise ValueError("no metal neighbour found on P atom")`. In `_stitch_fragment` the
   fragment is **metal-free** (trivalent P) pre-ETKDG, so the function raises
   immediately. The MiniPRD's "reuse one wildcard-dummy device" story is broken on
   both ends: the named class is gone, and the working helper does the inverse
   operation (convert, not attach). A genuinely new `attach` helper is required, and
   the Confidence Mandate should reflect that this is new code, not a revival.

3. **[HIGH — dummy leak via bidentate branch]** Task 1 gates all new behaviour on
   `_zone_a_p_expected_labels(frag_smiles)` being non-empty — which is **also true
   for bidentate Zone-A-P fragments** (DIPAMP, BDPP; both have `tests/fixtures/`
   XYZ files: `Rh-RR-DIPAMP-Cl2.xyz`, `PdCl2-RR-BDPP.xyz`). But the attach (Task 2)
   and the strip (Task 3) are scoped to the **monodentate** orientation block
   (`len(binding_idxs) == 1`, lines 1466–1485), which `return`s at 1485. The
   bidentate branch (1487+) has its own `return` at ~1519 with **no strip**. If a
   dummy is attached for a bidentate Zone-A-P fragment and control reaches the
   bidentate branch, **the Z=0 atom leaks into `combined_mol` and the XYZ block** —
   the document's own #1 forbidden outcome. The gate must be co-conditioned on
   `len(binding_idxs) == 1`, and the bidentate Zone-A-P case explicitly routed to
   the honest-warning DG path (consistent with the TASK-30/31 deferral).

4. **[HIGH — oracle self-consistency in Test 1]** Test 1 asserts
   `_metal_present_cip_label` gives labels "each matching `_zone_a_p_expected_labels`."
   These are **two different CIP conventions**: `_zone_a_p_expected_labels` /
   `_lp_cip_label` use a **Z=0 dummy (lowest priority, lone-pair convention)**, while
   `_metal_present_cip_label` uses a **real metal via a dative bond**. A real Pt
   (Z=78) is *high* CIP priority, not low — so the metal-present R/S string can be the
   **opposite** of the lone-pair R/S string for the identical physical configuration.
   MiniPRD-B's `_verify_zone_a_p` is careful to use LP-convention on *both* sides
   (SuperPRD B1, "never cross-convention"). Test 1 as written mixes conventions and
   may assert a false equality (or accidentally pass by cancellation). It should
   compare like-for-like — `_lp_cip_label` on the dummy copy vs `_zone_a_p_expected_labels`,
   the exact pair `_verify_zone_a_p` already trusts — and treat `_metal_present_cip_label`
   only as a documented, separately-derived cross-check.

5. **[MEDIUM — untested strip in the re-embed loop]** Test 6 asserts the re-embed
   loop's attempt counter stays 0 for the acceptance fixtures. The strip logic in
   the loop's `_stitch_fragment` call (seed `42 + attempts*1009`, line 1946) is
   therefore **never exercised by any test**. A strip bug localized to the loop
   path would pass all six acceptance tests, then in production produce
   `len(new_positions) != end - start` → the `if` guard at 1958 silently skips the
   update → the loop can't converge → a *false* "could not be enforced" warning (or,
   worse, a leaked dummy). Add a test that forces the loop to fire *with* the dummy
   path active (the `_test_flip_chiral_idx` seam already exists for exactly this).

---

## SECTION-BY-SECTION ANALYSIS

### §1 The Confidence Mandate — Analysis
* **Clarifying Questions:**
  - Given Findings 1–2, what is the *actual* reuse surface? Is the executor expected
    to (a) resurrect the deleted `PseudoAtomStrategy` from git history, (b) generalize
    `_build_dummy_metal_copy` into an `attach`-or-`convert` helper, or (c) write a
    fresh `_attach_dummy_metal(mol, p_idx)`? These have materially different diffs and
    blast radii, and the document currently implies (a), which is impossible.
  - The seed note self-rates confidence 8/10, citing residual risk "entirely in the
    placement/strip plumbing." Does that 8 already account for the attach helper being
    net-new (not revived) code? If not, the honest score is lower.
* **What-If Scenarios:**
  - Executor greps for `PseudoAtomStrategy`, finds nothing, and either invents a new
    class with that name (scope creep) or stalls on a clarifying question that should
    have been pre-answered here.
* **Points for Improvement:**
  - Replace every `PseudoAtomStrategy` / `core/chirality.py:22` reference with the
    real target. State plainly: *"the wildcard-dummy device was deleted; this MiniPRD
    re-introduces an* attach *helper, symmetric to `_build_dummy_metal_copy`'s* convert
    *logic, not a revival of removed code."*

### §2 Atomic User Stories — Analysis
* **Clarifying Questions:**
  - **US-C1:** In what order are `SetNoImplicit` / `AddHs` / dummy-attach applied?
    The existing code does `Chem.AddHs(rw.GetMol())` (line 1424), which appends H
    atoms **after** all heavy atoms. If the dummy is attached before `AddHs`, the
    H's land *after* the dummy, so the dummy is **not** the highest index (contradicts
    Task 4). If attached after `AddHs`, `NoImplicit(True)` must already be set on P or
    `AddHs` will hang a phantom H on it, making P 5-coordinate garbage.
  - **US-C2:** What P→dummy bond distance does ETKDG assign to a `P–*` (Z=0) bond?
    Wildcards have no covalent radius / UFF parameters. Confirm the alignment uses the
    **normalized** P→dummy direction only (distance-invariant), and that
    `EmbedMolecule` does not return −1 on the wildcard-bearing mol.
* **What-If Scenarios:**
  - **US-C1 / US-C4:** `AddHs` is documented to preserve chiral tags, but adding a 4th
    explicit neighbour to a P that RDKit modelled as "3 explicit + implicit lone-pair"
    can reinterpret the tetrahedral parity if the neighbour is inserted at a different
    ordinal position than the phantom it replaces. If the parity flips on attach,
    `[P@]` embeds as the wrong enantiomer *silently* and the loop masks it as
    "could not enforce." The encode side (`_build_dummy_metal_copy`) never embeds, so
    "it already works on encode" does **not** cover embed-time parity.
  - **US-C2:** If the wildcard makes `EmbedMolecule` fail, `_stitch_fragment` returns
    `None` → DG fallback → honest warning but **zero enforcement** — i.e. the bug is
    "fixed" into a different unenforced path for the exact fixtures meant to prove it.
* **Points for Improvement:**
  - Pin the operation order explicitly in US-C1 (e.g. "SetNoImplicit(P) → attach dummy
    with SINGLE bond + NoImplicit(dummy) → AddHs → embed → align → strip dummy →
    strip is index-stable because dummy is re-stripped *before* return").
  - Add an assertion that the embedded conformer's P handedness matches the input tag
    *before* placement (cheap `_lp_cip_label` on the pre-strip embedded fragment),
    turning US-C4's "guaranteed structurally" into a checked invariant.

### §3 Implementation Plan — Analysis
* **Clarifying Questions:**
  - **Task 1 vs Task 3:** Is the gate `_zone_a_p_expected_labels(frag_smiles)` non-empty
    **AND** `len(binding_idxs) == 1`? As written, Task 1 says "gate ALL new behaviour"
    on the former alone, but Task 3's strip lives only in the monodentate branch. See
    Finding 3 — this must be resolved or bidentate Zone-A-P fragments leak the dummy.
  - **Task 4:** Does `_fresh_fragment_mol` (used by `_zone_a_p_expected_labels`) and
    `MolFromSmiles(...,sanitize=False)` in `_stitch_fragment` produce **identical
    heavy-atom ordering**? Both are parse-order and MiniPRD-B already relies on
    `atom_start + local_p_idx` (line 1925), so this is likely fine — but Task 4 should
    *assert* it, not assume it, since the dummy attach is the first thing that could
    perturb it.
* **What-If Scenarios:**
  - **Task 3 strip → XYZ writer:** a leaked Z=0 atom reaches line 1985
    (`f"{sym:<2} ..."` where `sym = mol.GetAtomWithIdx(i).GetSymbol()`). Z=0's symbol
    is `*` (or empty), producing an **invalid element token in the XYZ block**. This
    directly breaks **Test 4** (byte-stable round trip: `XYZToSMILES` re-parse of the
    generated XYZ), and any downstream `xyz2mol` consumer of `GeneratedStructure`.
    That is the concrete failure chain behind the "#1 concern."
  - **Task 5 loop:** covered in Finding 5.
* **Points for Improvement:**
  - Make Task 4 emit a hard `assert` (not a comment) that `len(returned_positions) ==
    heavy+H count with the dummy absent`, and that no atom in the returned mol has
    `GetAtomicNum() == 0`. Cheap, and it converts the whole "dummy must not leak"
    negative-space rule into an enforced postcondition at the single choke point.

### §4 The Negative Space — Analysis
* **Clarifying Questions:**
  - The constraints forbid mirror/improper transforms and dummy leakage well, but say
    nothing about the **bidentate Zone-A-P** case (Finding 3). Is the intended
    behaviour "bidentate Zone-A-P → DG fallback + honest warning, unchanged from
    MiniPRD-B"? If so, state it as an explicit DO-NOT (do not attach a dummy when
    `len(binding_idxs) > 1`).
* **What-If Scenarios:**
  - "DO NOT widen scope to Zone-A **N**" is honoured by the P-only filter in
    `_zone_a_p_expected_labels` (`GetAtomicNum() != 15: continue`, line 366) — good,
    structurally enforced. No gap here.
* **Points for Improvement:**
  - Add: *"DO NOT attach the dummy on the bidentate/higher path; those fragments keep
    the MiniPRD-B DG-fallback + `_warn_zone_a_p_fallback` behaviour verbatim."*

### §5 Integration Tests & Verification — Analysis
* **Clarifying Questions:**
  - **Test 1:** resolve the convention mismatch (Finding 4). Which oracle is
    normative — `_lp_cip_label` (matches `_verify_zone_a_p`) or `_metal_present_cip_label`?
  - **Test 4 (PAMP fixture):** the fixture is a "Candidate Artifact requiring
    human-in-the-loop verification." Until a human validates its absolute
    configuration, Test 4 cannot be a merge gate — is the acceptance gate
    conditional on that validation landing first?
* **What-If Scenarios:**
  - **Test 5 (inertness):** asserts "dummy branch not entered" via "a spy/count or a
    pre/post XYZ byte-diff." A byte-diff proves *output* inertness but not that the
    branch wasn't entered-and-reverted; a spy/count is stronger. For the guarantee to
    mean what §4 claims (structural inertness), prefer the branch-entry counter.
  - **Test 6:** false-confidence risk — see Finding 5. Passing Test 6 (loop never
    fires) is compatible with a broken strip *inside* the loop.
* **Points for Improvement:**
  - Add **Test 7 (loop-with-dummy):** use `_test_flip_chiral_idx` to force one
    mis-embed on an SPL Zone-A-P fixture so the loop fires ≥1 attempt, then assert
    (a) convergence, (b) no Z=0 atom in the final mol, (c) `len(positions)` stable.
    This is the only test that exercises the strip on the re-embed path.

### §6 Acceptance Gate — Analysis
* **Clarifying Questions:**
  - The single-sentence gate inherits Test 1's oracle ambiguity ("opposite, correct
    metal-present CIPs"). Once Finding 4 is resolved, restate the gate in the chosen
    convention so "correct" is unambiguous.
* **Points for Improvement:**
  - Add "…and no generated XYZ contains a Z=0/`*` atom" to the gate — it is the
    cheapest possible guard against the highest-severity failure mode and belongs in
    the one-sentence gate.

### §7 Next Steps & Blast Radius (`architecture.yml`) — Analysis
* **Clarifying Questions:**
  - The change is confined to `atom_molassembler_adapter` (status `clean`). Downstream
    `generation.molassembler_adapter.haptic_face_correction` and
    `atom_generated_structure` consume its outputs — is a leaked/renumbered atom in
    `combined_mol` capable of shifting the eta-ring indices those nodes rely on?
    (Low risk: `_stitch_fragment` is the monodentate/bidentate path, eta uses
    `_stitch_eta_fragment` — but confirm the `all_pos`/`all_frag_idxs` index space is
    untouched for eta fragments interleaved with a Zone-A-P monodentate fragment.)
* **What-If Scenarios:**
  - A Zone-A-P monodentate fragment placed *before* an eta fragment in `all_pos`: if
    the dummy strip renumbers incorrectly, `eta_frag_ranges` offsets (line 1764) and
    the winding correction shift — a cross-feature regression touching Phase-3
    haptic-face code that has nothing to do with this fix.
* **Points for Improvement:**
  - The "chain may be short" note in §7 is optimistic given Findings 1–4. Recommend
    the `/hyper-resolve` pass explicitly close: (a) reuse-target correction, (b)
    bidentate gate, (c) Test-1 oracle, (d) AddHs/attach ordering — before Sonnet
    executes.

---

## SUMMARY TABLE

| # | Severity | Finding | Fix locus |
|---|----------|---------|-----------|
| 1 | BLOCKER | `PseudoAtomStrategy` / `PSEUDO_ATOMIC_NUM` deleted from `src/`; MiniPRD says "revive" it | §1 context, Task 2, §4 |
| 2 | BLOCKER | `_build_dummy_metal_copy` requires a metal neighbour → can't attach to metal-free fragment | §1, Task 2 |
| 3 | HIGH | Bidentate Zone-A-P (DIPAMP/BDPP) satisfies the gate but has no strip → dummy leaks | Task 1/3, §4 |
| 4 | HIGH | Test 1 mixes lone-pair (Z=0) and metal-present (real Pt) CIP conventions | §5 Test 1, §6 |
| 5 | MEDIUM | Strip on the re-embed-loop path is untested (Test 6 asserts loop never fires) | §5 |
| 6 | MEDIUM | Leaked Z=0 atom writes an invalid element token to XYZ → breaks Test 4 & round trip | Task 3, §6 |
| 7 | LOW | AddHs-vs-dummy index ordering and embed-time parity unspecified | US-C1, Task 4 |

**Final Action:** Run `/hyper-resolve` to triage these findings — Findings 1–4 should
be closed in the MiniPRD text before `/hyper-execute` on Sonnet.
