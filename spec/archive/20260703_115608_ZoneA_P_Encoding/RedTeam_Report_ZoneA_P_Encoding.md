# Red Team Report — Zone-A P Stereocenter Encoding (Stereo Roadmap Phase 4)

- **Target:** `spec/active/Draft_PRD.md` (v0.1.0, 2026-07-03)
- **Blast-radius reference:** `spec/compiled/architecture.yml`
- **Analyst:** Red Team Agent (adversarial pass, code-grounded)
- **Date:** 2026-07-03
- **Verdict:** The core design (Option A, lone-pair convention, verify-and-flip)
  is sound and unusually well-evidenced. However the PRD contains **one internal
  contradiction (the §3.1.4 warning oracle compares labels across two conventions
  the PRD itself declares incomparable)**, **one dependency on machinery that does
  not exist (Phase-3 reflection)**, and **one correctness gap (fragment-wide
  mirroring corrupts co-resident stereocenters)**. These must be resolved in
  `/hyper-resolve` before MiniPRD compilation.

---

## §1 Introduction & Goals — Analysis

### Clarifying Questions
1. **Metal predicate ownership.** "Metal-bound P (P bonded to a metal atom)" —
   which predicate defines *metal*? `TRANSITION_METALS_NUM` lives in
   `utils/xyz2mol.py:30`; `core/chirality.py` has no metal list and `mod_core`
   does not currently import it. Will `CIPAssigner` import from `xyz2mol`,
   duplicate the list (a new TD-005-style duplication), or should a shared
   constants module finally be created? Note the project has a *known* stale
   `is_metal` bug (TASK-04 / `@SP1`); building a second ad-hoc metal check
   multiplies that failure surface.
2. **Polymetallic scope.** The problem statement and §5.1 say "the metal"
   (singular). What happens for bimetallic complexes or a μ₂-bridging phosphide
   (P bonded to **two** metals)? Is the claim "always exactly 3 fragment-local
   neighbours" actually an invariant, or just true for the monometallic fixtures
   on hand?
3. **P–H phosphines.** A secondary phosphine (R₂PH–M) is Zone-A with
   `total_degree == 3` where one neighbour is implicit H. Do the spike results
   (which used C/C/aryl substituents only) extend to P–H, where the implicit-H
   handling of `MolFromSmiles`/`AddHs` interacts with the trivalent tag?

### What-If Scenarios
1. **Bridging phosphide:** P bonded to Rh and Rh′. The dummy-metal copy
   ("keep only its bond to this P") makes P 5-coordinate, or — if only one metal
   is swapped — leaves a live transition metal in the CIP digraph that is absent
   from the fragment view, silently breaking the §5.1 equivalence. The label is
   stored, survives, and is **wrong** — a lossless-looking lie, strictly worse
   than today's honest drop.
2. **Metal identification miss:** a complex whose centre element is missing
   from (or mis-flagged by) the metal predicate → the P is treated as Zone-B,
   the metal-present `_OIN_CIPCode` is used to flip a trivalent tag → wrong
   convention emitted with full confidence.

### Points for Improvement
1. Pin the metal predicate to a single shared source of truth in the MiniPRD
   (explicit import path), and state the monometallic assumption as a **guard**:
   if a candidate P is bonded to ≥2 metal atoms, store nothing and
   `warnings.warn` (fall back to today's clearing behaviour).
2. Add "P bonded to exactly one metal atom" to the formal Zone-A definition used
   by both MiniPRDs, so encode and generation sides agree on eligibility.

---

## §2 Confidence Mandate — Analysis

### Clarifying Questions
1. **Q1 sharpened — which CIP implementation is authoritative?** The
   verify-and-flip in `recover()` recomputes CIP via the **legacy**
   `Chem.AssignStereochemistry` (`_CIPCode`), while the diagnostic oracle uses
   **`rdCIPLabeler`** (the modern, rules-complete implementation). These two are
   known to disagree on non-trivial cases (CIP rules 4b/5, like/unlike
   descriptors). If they disagree on a phosphine, the pipeline flips to satisfy
   legacy while the warning fires against rdCIPLabeler — a permanent,
   unactionable warning. Which one defines `_OIN_CIPCode_LP`?
2. **Q2 sharpened:** in `_template_generate` (`molassembler_adapter.py:977`),
   the assembled complex exists as `fragment_mol_parts` + positions only when
   `has_all_mols` is true; `GeneratedStructure.mol` is `Optional` and is `None`
   for eta-fallback cases. What is the specified behaviour of verify-and-reflect
   when **no assembled RDKit mol exists** (eta fallback, DG fallback)? Silent
   skip is a stereo hole; the PRD must name it.
3. **Q4 sharpened:** the spurious-tag gate relies on
   `AssignStereochemistry(cleanIt=True)` refusing a `_CIPCode` on
   graph-symmetric P in the dummy copy. Has this been spiked on the *actual*
   BDPP/BDNN dummy-copies (with a `*` fourth substituent), rather than inferred?
   A `*` neighbour changes the symmetry analysis input.

### What-If Scenarios
1. **Legacy/new CIP divergence on DIPAMP itself:** verify-and-flip settles on
   the legacy label; US-A4's rdCIPLabeler check disagrees → the "no warning on
   clean DIPAMP" acceptance criterion (US-A4.2) fails not because the code is
   wrong but because the PRD mixed oracles.

### Points for Improvement
1. Resolve Q1 by **decree, not test**: pick `rdCIPLabeler` (or legacy)
   end-to-end for every `_OIN_CIPCode_LP` computation, recomputation, and
   cross-check, and write the choice into the MiniPRD contract. Add the
   raw-parity unit test on top.
2. Add an explicit confidence deduction for the generation side: the encode
   side has three spikes behind it; the generation side's verify-and-reflect
   has **zero** (no reflection code exists yet — see §5 analysis).

---

## §3 Scope — Analysis

### Clarifying Questions
1. **Property precedence.** After the change, a metal-bound P may carry **both**
   `_OIN_CIPCode` (metal-present sense, set by the existing loop at
   `chirality.py:104-108`) and `_OIN_CIPCode_LP`. The current `recover()`
   branches on `total_deg` **first**. If a Zone-A P ever presents with
   `total_deg >= 4` in the fragment (implicit-H promotion during sanitize,
   protonated phosphine, bridging cases), the `elif stored_cip and total_deg >= 4`
   branch flips against the **metal-present** code — wrong convention. What is
   the exact branch order and property precedence? Should `assign_all()` stop
   setting `_OIN_CIPCode` on Zone-A P entirely?
2. **PseudoAtomStrategy deletion blast radius.** `architecture.yml` has
   `atom_chirality_recovery.edges.depends_on: [atom_pseudo_atom_strategy]` — the
   PRD deletes the node but does not mention the dangling edge. Also
   `tests/unit/test_stereo_roundtrip_diagnostics.py:181` references
   `PseudoAtomStrategy` in a docstring, and `recover()`'s final `else` branch
   comment cites it. Does US-A3's "no remaining import/reference" include
   docstrings/comments, and does the ≥4-neighbours-no-CIP **clearing behaviour**
   (currently labelled "PseudoAtomStrategy fallback") survive the deletion?
3. **Sanitization of the dummy copy.** C-3 makes `Chem.SanitizeMol` a hard
   precondition for CIP assignment. Deleting all metal–ligand bonds except the
   P bond can strand aromatic eta-ligand systems (a Cp ring that was
   charge-balanced by the metal) → `SanitizeMol` on the copy **throws**. Is the
   dummy copy sanitized, and what happens on failure?

### What-If Scenarios
1. **CpRh(PR₃)-type complex (eta ligand + chiral phosphine):** dummy-copy
   construction breaks the Cp aromaticity, sanitize raises, the exception
   propagates out of `assign_all()` per its documented contract → **the entire
   `convert()` fails for a complex family that round-trips today.** This is a
   regression far outside the stereo feature: the blast radius of an unguarded
   per-atom exception is the whole encode pipeline.
2. **Charge leakage:** the metal→Z=0 swap keeps the metal's formal charge
   (e.g. Rh⁺ → `[*+]`). A charged dummy can change sanitize outcomes and,
   in principle, isotope/charge-sensitive ranking. Nothing in the PRD zeroes it.

### Points for Improvement
1. Specify the dummy-copy recipe precisely: swap Z→0, **zero formal charge,
   clear isotope, drop all other metal bonds AND all other metals' bonds**, then
   sanitize inside a `try/except` where failure ⇒ store nothing + warn
   (degrades to today's behaviour instead of crashing).
2. Add the `architecture.yml` edge removal
   (`atom_chirality_recovery → atom_pseudo_atom_strategy`) to MiniPRD-A's
   deliverables so `/hyper-audit` reconciles cleanly.
3. State the branch order for the new `recover()` explicitly:
   `_OIN_CIPCode_LP` present ⇒ Zone-A verify-and-flip **regardless of degree**;
   only then fall through to the existing degree-keyed branches.

---

## §4 User Stories — Analysis

### Clarifying Questions
1. **US-A1.3 / US-A2 interplay:** for a symmetric P, `assign_all()` stores no
   `_LP` prop — but the **chiral tag** copied through fragmentation (set by the
   metal-present `AssignAtomChiralTagsFromStructure`, which tags any pyramidal
   geometry) is still on the fragment atom. The Zone-A-without-property branch
   clears it — good — but is there a test asserting the tag (not just the
   property) is absent from the emitted SMILES for BDPP/BDNN? Byte-identical
   goldens cover it implicitly; make it explicit.
2. **US-B1.3 testability:** "post-assembly verify-and-reflect corrects a
   mis-embedded pyramid" — ETKDG *reliably* builds the correct pyramid (spike 2),
   so under normal operation this branch never fires. How is a mis-embed
   **forced** in a test? Without an injection point (e.g. a test-only
   pre-mirrored fragment or a mocked embed), the reflect branch ships untested.
3. **US-B2 byte-stability scope:** byte-stable across *what*? Same process, same
   RDKit version, same seed? Canonical SMILES with stereo tokens has changed
   across RDKit releases; a golden asserting byte equality is implicitly pinning
   the RDKit version.

### What-If Scenarios
1. **Interdependent stereocenters in one pass:** DIPAMP's fragment contains
   **both** Zone-A P atoms. `recover()` runs `AssignStereochemistry` **once**
   before the loop (`chirality.py:145`); if the loop flips P₁'s tag, P₂'s
   `_CIPCode` (computed pre-flip) can be stale wherever CIP rules 4b/5 make one
   centre's descriptor depend on the other's. Result: P₂ flipped (or not
   flipped) against a stale label — an enantiomer-scrambling race entirely
   inside one function.
2. **Warning deduplication swallows the second P:** Python's default warning
   filter dedupes by (message, category, location). Two conflicting P atoms →
   identical warn call site → user sees **one** warning and fixes one atom.

### Points for Improvement
1. Require `recover()` to **recompute CIP after each flip** (or iterate to a
   fixed point with a bounded pass count) when a fragment holds >1 tagged P.
2. Specify a dedicated warning class (e.g. `OINStereoWarning(UserWarning)`),
   include the atom index **in the message string** (defeats dedup, aids
   filtering), and require tests to use `warnings.catch_warnings(record=True)`.
3. Add a US-B1 acceptance criterion for the *no-mol* generation paths (eta/DG
   fallback): defined behaviour (skip + warn) rather than undefined.

---

## §5 Technical Specifications — Analysis

### Clarifying Questions
1. **THE INTERNAL CONTRADICTION (blocking).** §3.1 item 4 and US-A4 specify the
   diagnostic as: run `rdCIPLabeler` CIP-from-3D **on the metal-present
   complex** and warn if it "conflicts with the stored `_OIN_CIPCode_LP`". But
   §5.1's opening paragraph states these two views **can legitimately differ**
   ("the R/S letter can legitimately differ between the two views of the same
   geometry") — that asymmetry is the entire reason the LP convention exists.
   As written, the warning fires on *correct* output whenever the metal-present
   and lone-pair labels differ. The cross-check must run `rdCIPLabeler` **on the
   dummy-metal copy** (same convention), or compare metal-present-vs-metal-present
   (`_OIN_CIPCode` vs `rdCIPLabeler` on the full mol) as a *separate* sanity
   check. Which is intended?
2. **"Phase-3-style reflection machinery" does not exist.** `grep -i
   "reflect\|mirror"` over `generation/molassembler_adapter.py` returns nothing,
   and Phase 3 (`Draft_PRD_StereoPhase3_HapticFace.md`) is itself an
   **unresolved draft** sitting in `spec/active/` in parallel. MiniPRD-B is
   therefore specified against machinery that is, today, vapor. Does MiniPRD-B
   (a) build its own reflection, (b) block on Phase 3 landing first, or
   (c) get re-scoped to share one reflection utility with Phase 3? Running both
   drafts through `/hyper-resolve` independently invites two divergent mirror
   implementations in the same file.
3. **Stale line reference:** `_molassembler_worker` is at
   `molassembler_adapter.py:1489`, not `:1418` as stated in §3.1 item 7 and the
   decision doc. Minor, but the MiniPRD will be executed by a cheaper model that
   takes line numbers literally.
4. **Downstream index mapping (verify-unchanged list is right, but name the
   mechanism):** after `recover()`, `xyz2mol.py:949` re-serializes and
   `xyz2mol.py:957` re-parses with `sanitize=False` to rebuild the
   fragment→SMILES atom map, and `_count_smiles_atoms_before`
   (`oin/inline.py:5`) counts atoms by token scanning. Both handle bracket atoms,
   so `P{0}` → `[P@]{0}` should be transparent — but has the token change been
   traced through `SLOT_REGEX` adjacency (a `{0}` that previously followed `P`
   now follows `]`)? One unit test on `parse_inline_string` with `[P@]{0}` and
   `[P@@]{1>}` (tag + winding co-occurrence) closes this.

### What-If Scenarios
1. **Fragment-wide mirroring corrupts co-resident stereocenters (blocking).**
   §3.1 item 6 fixes a wrong P pyramid by mirroring **the whole fragment**
   across the plane of P's three substituents. A mirror is an improper
   operation: it inverts **every** tetrahedral centre in the fragment. DIPAMP
   survives only because its backbone is achiral. Take a DIPAMP-like ligand
   with a stereogenic backbone carbon (BDPP-style backbone + P-stereogenic
   ends): reflecting to fix P silently enantiomerizes the carbons — the
   "correction" manufactures a diastereomer. Same failure inside one ligand
   with two P centres where only one mis-embeds: the mirror fixes one and
   breaks the other, and re-verification enters a flip-flop loop.
2. **Unbounded verify-reflect-replace loop:** mismatch → mirror → re-place →
   re-verify. Nothing in the PRD bounds this. A pathological fragment (or the
   diastereomer case above, which mirroring can *never* fix) loops until the
   60 s ProcessPoolExecutor timeout kills the whole generation, converting a
   stereo warning into a hard product failure.
3. **rdCIPLabeler pathological runtime:** the new CIP labeler has documented
   exponential corner cases on highly symmetric fused systems. Running it on
   every metal-present complex during `assign_all()` puts an unbounded
   computation on the hot encode path with no timeout.

### Points for Improvement
1. **Replace fragment mirroring with per-centre correction where possible:**
   invert only the P pyramid (swap two substituent positions locally, or
   re-embed the fragment from the corrected tag with a new seed) instead of a
   global mirror; permit the global mirror only when the fragment's *only*
   stereocentres are the mis-embedded P set and all of them disagree. Encode
   this decision table in MiniPRD-B.
2. **Bound the enforcement loop:** max 1 reflect + 1 re-embed attempt (NFR);
   on persistent mismatch, emit the structure with a hard `OINStereoWarning`
   rather than dying in the timeout.
3. Guard the rdCIPLabeler call (per-call try/except + optional time budget) and
   make the diagnostic skippable via a flag for batch users.
4. Fix the §5.2 blast radius list: add `oin/inline.py::parse_inline_string` /
   `_count_smiles_atoms_before` to the verify-unchanged audit set, and correct
   the `:1418` line reference.

---

## §6 Negative Constraints — Analysis

### Clarifying Questions
1. "**DO NOT** compute or store a fragment-local CIP for a *trivalent* P
   directly" — but `recover()`'s verify-and-flip **does** recompute trivalent
   CIP from the tag (per §3.1 item 2, backed by spike 1). The constraint as
   worded contradicts the mechanism unless it is scoped to *3D perception*
   ("do not derive a tag from 3D coordinates of a trivalent P"). Reword before a
   literal-minded executor reads it as banning the recover() recompute.
2. Is "DO NOT regress carbon `@/@@`" backed by a named test set (TASK-10) run
   in CI, or is it aspirational? US-A5.4 names it — good — but only for the
   encode MiniPRD; MiniPRD-B touches the adapter every carbon path also uses.

### What-If Scenarios
1. **Spurious-tag gate bypass via property survival:** atom properties survive
   `AddAtom` copies and pickling. If any *earlier* pipeline stage (or a stale
   test fixture mol) already set `_OIN_CIPCode_LP`, the gate reads presence, not
   provenance — a leftover property on a symmetric P emits a phantom tag.
   Low probability today, but the property name is now load-bearing API.

### Points for Improvement
1. Add a constraint pair: `assign_all()` must **clear any pre-existing
   `_OIN_CIPCode_LP`** before recomputing (idempotence), mirroring
   `cleanIt=True` semantics.
2. Add: **DO NOT** apply an improper (mirror) transform to any fragment
   containing stereocentres other than the target P set (feeds §5 improvement 1).

---

## §7 Risks & Mitigation — Analysis

### Clarifying Questions
1. **RISK-2's mitigation is "Red Team to pin the object."** Pinned:
   the assembled-complex RDKit mol in the adapter (`_template_generate`'s
   combined-mol assembly) is the right layer — it is the only place where both
   the metal bond and per-fragment atom provenance exist *before* positions are
   frozen into the XYZ string. Engine-level (`OIN3DGenerator`) only sees
   `GeneratedStructure`, whose `.mol` is `Optional` and post-hoc. But this
   forces the answer to §2-Q2: eta/DG-fallback paths (mol `None`) get **no
   enforcement** — the risk register must carry that residual risk explicitly.
2. RISK-5's "literature configuration cross-check" — cross-check against what,
   concretely? (R,R)-DIPAMP's P configuration is literature-documented; name
   the citation or CCDC refcode in the MiniPRD so HITL isn't a vibe check.

### What-If Scenarios
1. **Missing risk — encode-side crash regression** (the eta-ligand sanitize
   failure from §3 What-If 1). Every risk in the register is about stereo
   *quality*; none covers `assign_all()` newly **throwing** on complexes that
   convert fine today. That is the highest-severity unlisted risk.
2. **Missing risk — dual-CIP-implementation disagreement** (§2 above): flips
   satisfy legacy, warnings judge by rdCIPLabeler → unresolvable warning state.

### Points for Improvement
1. Add RISK-7: dummy-copy construction/sanitize failure ⇒ must degrade to
   store-nothing + warn, never propagate. Add a CpM(PR₃)-type fixture (even a
   synthetic one) as its regression test.
2. Add RISK-8: enforcement loop non-termination ⇒ bounded retries NFR.
3. Add RISK-9: fallback generation paths (mol=None) are stereo-unenforced ⇒
   documented residual risk + warning.

---

## §8 Success Metrics — Analysis

### Clarifying Questions
1. Are the generation-side metrics seeded/deterministic? "Assert on derived
   metal-present CIP" is the right oracle, but ETKDG across environments can
   still fail to embed at all; is a fixed seed part of the test contract?
2. What is the metric for the **warning channel** beyond "fires only on genuine
   conflict"? Genuine conflict is currently undefined because of the §5
   convention contradiction — the metric inherits that ambiguity.

### What-If Scenarios
1. **Green suite, wrong science:** every listed metric can pass while the
   bridging-P / polymetallic case emits confidently wrong labels, because no
   fixture exercises anything but monometallic Rh/Pd. Metrics measure the happy
   path only.

### Points for Improvement
1. Add a negative-capability metric: for any P where the LP label cannot be
   computed (multi-metal, sanitize failure, no CIP), the emitted OIN is
   **identical to today's output** and a warning names the atom — "graceful
   degradation is observable."
2. Add: `uv run` suite green **with `-W error::OINStereoWarning`** on the clean
   fixtures — turns "no spurious warnings" into an enforced gate instead of a
   sentence.

---

## §9 HITL Candidate Artifact — Analysis

### Clarifying Questions
1. The reviewer confirms (R,R) from a raw OIN string plus an XYZ file — by what
   means? An OIN `[P@]` token does not visually map to R/S without doing the CIP
   analysis in your head. Should the execution run also emit a depiction
   (mol block / image) or at least the per-atom rdCIPLabeler table it already
   promises, in a copy-pasteable form for the sign-off note?
2. Two copies of the fixture exist: `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz` and
   `tests/integration/Rh-RR-DIPAMP-Cl2.xyz` (both currently untracked). Which is
   canonical, and is the other deleted? A golden signed off against one file
   while tests read the other is a silent-drift channel.

### What-If Scenarios
1. **Sign-off against the wrong convention:** the reviewer confirms "(R,R)"
   thinking in metal-present CIP (the literature convention for DIPAMP
   complexes), while the string stores lone-pair-convention labels. If the two
   conventions differ for DIPAMP's P atoms, a *correct* string gets rejected —
   or a wrong one approved. The reviewer instructions never say which
   convention the `@`/`@@` encodes.

### Points for Improvement
1. Reviewer instructions must state explicitly: "the tag encodes the
   fragment-local (lone-pair) sense; the printed rdCIPLabeler table shows the
   metal-present sense; they may legitimately differ — confirm against the
   *metal-present* (R,R) column." Record the sign-off (who/when/verdict) in
   `spec/worklog/`.
2. Deduplicate the fixture before promotion; the golden's provenance line
   should name the exact fixture path + SHA.

---

## Consolidated Blocking Items (for `/hyper-resolve`)

| # | Severity | Item | PRD Ref |
|---|----------|------|---------|
| B1 | **Blocker** | Warning oracle compares across conventions the PRD declares incomparable — respecify against the dummy-metal copy | §3.1.4, US-A4, §5.1 |
| B2 | **Blocker** | Verify-and-reflect depends on nonexistent "Phase-3-style" reflection machinery; Phase 3 is an unresolved parallel draft — pin build/borrow/block | §3.1.6, §5.1 |
| B3 | **Blocker** | Global fragment mirror inverts co-resident stereocenters (diastereomer manufacture); needs per-centre correction decision table | §3.1.6 |
| B4 | **High** | Dummy-copy sanitize failure (eta ligands, charges) can crash `convert()` for complexes that work today — mandate degrade-to-clear + warn | §5.1, §7 |
| B5 | **High** | Single vs dual CIP implementation (legacy `AssignStereochemistry` vs `rdCIPLabeler`) must be pinned end-to-end | §2 Q1, §5.3 |
| B6 | **High** | Multi-P fragments: recompute CIP after each flip in `recover()` (stale-label race) | §3.1.2 |
| B7 | **Medium** | Polymetallic/bridging-P guard: ≥2 metal bonds ⇒ store nothing + warn | §1.1, §5.1 |
| B8 | **Medium** | Enforcement loop retry budget NFR; fallback paths (mol=None) explicitly unenforced + warned | §3.1.6, §7 |
| B9 | **Medium** | `recover()` branch precedence: `_LP` property beats degree check; `assign_all` idempotence (clear stale `_LP`) | §3.1.2, §6 |
| B10 | **Low** | Housekeeping: `architecture.yml` dangling edge to `atom_pseudo_atom_strategy`; stale `:1418` line ref (now `:1489`); fixture duplication (fixtures/ vs integration/); reword the "do not compute trivalent CIP" constraint | §3.1.3, §3.1.7, §6, §9 |

---

**Next step:** run `/hyper-resolve` to triage these findings and compile the
final SuperPRD + MiniPRDs.
