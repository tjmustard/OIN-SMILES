# Draft PRD — Zone-A P Stereocenter Encoding (Stereo Roadmap Phase 4)

## Metadata
- **Project Name**: OIN-SMILES — Metal-bound (Zone-A) P stereocenter encoding, Phase 4
- **Version**: 0.1.0 (Draft — pre `/hyper-redteam`)
- **Status**: Draft — ready for adversarial analysis
- **Owner**: Architect Agent (Fable) / Thomas Mustard
- **Roadmap ref**: `spec/worklog/ROADMAP-stereo.md` § Phase 4 (absorbs superseded Phase 2)
- **Decision ref**: `spec/worklog/PHASE4-decision.md` (CTO consult, 2026-07-03)
- **Design brief**: `spec/worklog/PHASE4-design-brief.md`
- **Depends on**: TASK-20 (diagnostic, DONE 2026-07-03); fixture `Rh-RR-DIPAMP-Cl2.xyz` (built 2026-07-03)
- **Child MiniPRDs**: `MiniPRD_ZoneA_P_Encode.md` (A), `MiniPRD_ZoneA_P_GenEnforce.md` (B)

---

## 1. Introduction & Goals

### 1.1 Problem Statement
OIN silently drops the stereochemistry of any phosphorus atom bonded directly to
the metal ("Zone A"). Verified 2026-07-03 (TASK-20) with `Rh-RR-DIPAMP-Cl2.xyz`
(both P atoms are genuine CIP stereocenters): `XYZToSMILES().convert()` emits
`[Rh_SPL].Cc1ccccc1P{0}(CCP{1}(c1ccccc1)c1ccccc1C)c1ccccc1.[Cl]{2}.[Cl]{3}` —
**no `@`/`@@` on either P.**

Mechanism: `ChiralityRecoveryUtility.recover()` (`core/chirality.py:155-158`)
unconditionally clears the chiral tag on any P/N atom with `total_degree < 4` in
the post-fragmentation ligand mol. A metal binder — with the metal excluded from
the fragment SMILES by OIN's construction — always has exactly 3 fragment-local
neighbours, so it is *always* Zone A and *always* stripped. The strip happens
BEFORE the branch that would consume the stored CIP code.

**The information is discarded, not lost.** `CIPAssigner` computes CIP from the
intact 3D structure (metal present) and attaches it to the atom; `recover()`
throws it away. A fix has ground truth to work from.

This breaks the project's core **lossless round-trip** promise for an entire,
important ligand family (chiral phosphines — DIPAMP, DiPAMP-type P-stereogenic
diphosphines used in asymmetric catalysis).

### 1.2 Solution Overview
Per `PHASE4-decision.md` (empirically validated by three RDKit spikes):

**Option A — lone-pair `[P@]`/`[P@@]` convention, phosphorus only.** Trivalent
`[P@]` is Daylight-legal, parses to a stable CIP label, survives sanitize /
AddHs / RemoveHs / 20 randomized renumber round-trips, and ETKDG builds the
correct pyramid from it with zero embed changes. No OIN grammar change and no
version bump — `@`/`@@` inside fragment SMILES is existing v3.x syntax; Zone-A P
tags are new *content*, not new *syntax*.

- **MiniPRD-A (encode side):** `CIPAssigner.assign_all()` stores a fragment-local
  (lone-pair-convention) CIP label `_OIN_CIPCode_LP` on each metal-bound P,
  computed while 3D is available via the **dummy-metal equivalence** (§5.1);
  `ChiralityRecoveryUtility.recover()` replaces the Zone-A *clear* with a
  *verify-and-flip* keyed on `_OIN_CIPCode_LP`. Delete the dead
  `PseudoAtomStrategy` (Option B, formally rejected).
- **MiniPRD-B (generation side):** parse already passes `@/@@` through untouched;
  ETKDG builds the pyramid; add a **post-assembly verify-and-reflect** on the
  fully-assembled complex (P is 4-coordinate → 3D perception valid) that mirrors
  the fragment if the re-derived label disagrees with the input tag.

**Nitrogen is explicitly OUT of scope** (§3.2): RDKit clears trivalent `[N@]` as
non-stereogenic (amine inversion), so an in-fragment tag on a metal-bound N
cannot survive any RDKit pass. Zone-A N is deferred to a future Option-C
out-of-band marker if a fixture and demand appear.

### 1.3 Target Audience
Internal: the XYZ→OIN and OIN→XYZ pipelines and their maintainers. Unblocks the
superseded Phase 2 (ETKDG stereo verification) by producing a P tag that
actually survives the fragment boundary. External: any conformant SMILES
consumer of OIN strings (accepts `[P@]` — standard Daylight).

---

## 2. Confidence Mandate
**Confidence Score**: 8/10 (Draft — pre Red Team)

Evidence base (RDKit spikes, `uv run`, 2026-07-03 — logged in `PHASE4-decision.md` §Evidence):
- String stability under 20 renumber round-trips: **proven**.
- ETKDG honors trivalent `[P@]`/`[P@@]` pyramid sense across 5 seeds: **proven**.
- 3D perception fails for trivalent P but works 4-coordinate: **proven** (shapes the design — all geometry↔tag work happens metal/dummy-present).
- Trivalent `[N@]` is cleared by RDKit: **proven** (scopes N out).

Open questions for Red Team:
- **Q1 (label parity):** does the dummy-metal LP-CIP computed encode-side equal
  the fragment-local CIP RDKit computes from the *trivalent* `[P@]` after metal
  removal, atom-order-independently? (The verify-and-flip in `recover()` is
  designed to absorb a mismatch, but the raw-parity assumption needs a unit test.)
- **Q2 (generation enforcement placement):** should verify-and-reflect live in
  the adapter (`_template_generate`) or one layer up in `OIN3DGenerator`? Which
  assembled-mol object has both the metal bond and a usable conformer?
- **Q3 (molassembler fallback):** does `masm.io.experimental.from_smiles` respect
  trivalent `[P@]`, or is an atom stereopermutator required in the worker?
- **Q4 (spurious-tag gate):** does "no `_OIN_CIPCode_LP` property ⇒ no tag"
  reliably keep symmetric phosphines (BDPP/BDNN two-identical-phenyl P) tag-free?

---

## 3. Scope

### 3.1 In-Scope
**MiniPRD-A — Encode side (`core/chirality.py`):**
1. `CIPAssigner.assign_all()`: for each metal-bound P (P bonded to a metal atom),
   compute the lone-pair-convention CIP via the dummy-metal equivalence (§5.1)
   while the 3D conformer is present; store as atom prop `_OIN_CIPCode_LP` ∈
   {`R`,`S`}. No label computable ⇒ store nothing (the spurious-tag gate).
2. `ChiralityRecoveryUtility.recover()`: replace the `total_deg < 4` **clear**
   with: if the atom is P and carries `_OIN_CIPCode_LP` → keep the chiral tag,
   recompute the fragment-local (trivalent) CIP, flip the tag on mismatch (the
   same verify-and-flip pattern already used for 4-coordinate zones). P/N without
   the property, and **all** N, keep today's clearing behaviour.
3. Delete `PseudoAtomStrategy` (class + `atom_pseudo_atom_strategy` node);
   remove `PSEUDO_ATOMIC_NUM` if unreferenced elsewhere.
4. **RDKit diagnostic warning:** during `CIPAssigner`, run `rdCIPLabeler` CIP-from-3D
   on the metal-present complex; if a metal-bound P's perceived label conflicts with
   the stored `_OIN_CIPCode_LP` (or perception errors/returns None where a tag was
   set), emit a `warnings.warn` naming the atom index — surfaced to the user per
   their HITL request.
5. Goldens: regenerate `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt` (see §4
   — Candidate Artifact, HITL sign-off required).

**MiniPRD-B — Generation side (`generation/molassembler_adapter.py`, maybe
`generation/engine.py`):**
6. Post-assembly **verify-and-reflect**: on the assembled complex (metal-bonded P
   → 4-coordinate → `AssignStereochemistryFrom3D` valid), recompute the
   lone-pair-convention label via the dummy-metal trick and compare to the input
   fragment tag; on mismatch, mirror the P-containing fragment across the plane of
   P's three substituents and re-place (Phase-3-style reflection machinery).
7. Molassembler fallback: investigate (Q3); if `from_smiles` drops trivalent
   `[P@]`, set an atom stereopermutator in `_molassembler_worker`
   (`molassembler_adapter.py:1418`).
8. Tests: enantiomer-discrimination (flip `@↔@@` in DIPAMP OIN → opposite
   metal-present CIP on regenerated structures); full XYZ→OIN→XYZ→OIN byte-stable
   round-trip on DIPAMP.

### 3.2 Out-of-Scope (explicit)
- **Nitrogen Zone-A stereo.** Trivalent `[N@]` is RDKit-cleared; deferred to a
  future Option-C out-of-band marker. `recover()` keeps clearing Zone-A N.
- **Option B (wildcard `*` pseudo-atom)** — formally rejected; scaffolding deleted.
- **Option C (OIN grammar / out-of-band marker)** — not built; reserved for N.
- **OIN version bump** — none; this is content-level within existing v3.x grammar.
- **Carbon `@/@@` paths** — already correct; must not regress.
- **Axial/atropisomeric P (BINAP)** — separate concern (`test_axial_chiral.py`),
  not a tetrahedral P stereocenter.
- **Winding / haptic face** — Phases 1/3, separate SuperPRDs.
- **Direct-parser nodes** (`parse_oin_direct`) — deferred to v0.2.2.
- **Legacy V2.4 sidecar / `oin/parser.py` path.**

---

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
| :-- | :-- | :-- | :-- |
| US-A1 | As the encode pipeline, I store a fragment-local CIP for each metal-bound P so it survives fragmentation. | 1. After `CIPAssigner.assign_all()` on `Rh-RR-DIPAMP-Cl2` mol, both P atoms carry `_OIN_CIPCode_LP` ∈ {R,S}.<br>2. Computed via dummy-metal copy (metal→Z=0, single P bond retained), 4-coordinate perception.<br>3. A symmetric phosphine P (two identical phenyls) gets NO `_OIN_CIPCode_LP`. | High |
| US-A2 | As the encode pipeline, `recover()` keeps and corrects the P tag instead of clearing it. | 1. `XYZToSMILES().convert()` on DIPAMP emits a `@`/`@@` on **both** `P{0}` and `P{1}`.<br>2. Zone-A P without `_OIN_CIPCode_LP` → cleared (unchanged).<br>3. All Zone-A N → cleared (unchanged).<br>4. Carbon `@/@@` unaffected. | High |
| US-A3 | As a maintainer, dead Option-B code is gone. | 1. `PseudoAtomStrategy` class and `atom_pseudo_atom_strategy` node removed.<br>2. No remaining import/reference. | Medium |
| US-A4 | As the user, I get a warning if RDKit's CIP-from-3D disagrees with the emitted tag. | 1. On the metal-present complex, `rdCIPLabeler` label conflict (or perception error) for a tagged P → `warnings.warn` naming the atom index.<br>2. No warning on the clean DIPAMP (R,R) case. | High |
| US-A5 | As the test suite, the DIPAMP golden is a verified Candidate Artifact. | 1. New OIN string written to `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt`.<br>2. NOT promoted to a `tests/fixtures/` trusted golden until HITL confirms (R,R).<br>3. BDPP/BDNN goldens byte-identical (negative control).<br>4. Carbon-chirality round-trips (TASK-10) stay green. | High |
| US-B1 | As the generation pipeline, the input P tag is preserved and enforced through embedding. | 1. Regenerating DIPAMP OIN yields a complex whose metal-present CIP on both P = (R,R).<br>2. Flipping `@↔@@` in the OIN yields opposite metal-present CIP.<br>3. Post-assembly verify-and-reflect corrects a mis-embedded pyramid. | High |
| US-B2 | As the round-trip contract, DIPAMP is lossless. | 1. XYZ→OIN→XYZ→OIN is byte-stable on the OIN string.<br>2. RDKit CIP-from-3D on the regenerated metal-present complex matches the original. | High |

---

## 5. Technical Specifications

### 5.1 Architecture — the lone-pair label convention (the correctness crux)
CIP labels are context-dependent: with the metal present it is usually P's
highest-priority substituent; in the fragment the lone pair is lowest. The R/S
letter can legitimately differ between the two views of the same geometry, so the
MiniPRD must pin ONE convention.

**The OIN string stores the fragment-local (lone-pair) CIP sense.** Computable
encode-side (3D still present) via the **dummy-metal equivalence**: a dummy atom
(Z=0) is also lowest priority and sits exactly where the metal was — i.e. where
the lone pair points. Therefore:

> lone-pair CIP of the fragment ≡ CIP of the 4-coordinate P in a copy of the full
> mol where the metal is swapped to Z=0, keeping only its bond to this P.

This makes the label 4-coordinate (RDKit 3D perception valid — spike 3) and
identical in priority-ordering to the trivalent fragment view (lowest-priority
substituent occupies the same spatial position either way).

Data flow (no topology change; new atom property + reused tag channel):
```
XYZ → get_tmc_mol() [full mol, 3D, metal present]
  → CIPAssigner.assign_all()
      for each metal-bound P:
        build dummy-metal copy → AssignAtomChiralTagsFromStructure + AssignStereochemistry
        store _OIN_CIPCode_LP on the real atom   [NEW]
        rdCIPLabeler cross-check → warn on conflict   [NEW]
  → fragmentation (metal excluded from fragment SMILES; tag copied by AddAtom)
  → OINSanitizer.generate_robust_smiles() → (smiles, fragment mol)
  → ChiralityRecoveryUtility.recover()
      Zone-A P WITH _OIN_CIPCode_LP: keep tag, recompute trivalent CIP, flip on mismatch  [CHANGED]
      Zone-A P/N WITHOUT it, and all N: clear (unchanged)
  → OIN string carries [P@]/[P@@]

OIN string → parse_inline_string [regex-only, @/@@ passes through — no change]
  → OINParser.parse → fragment SMILES with [P@]
  → _template_generate: MolFromSmiles(sanitize=False) → AddHs → ETKDG   [pyramid built from tag]
  → assembled complex [P now 4-coordinate]
  → verify-and-reflect: dummy-metal LP-CIP vs input tag; mirror on mismatch   [NEW]
```

### 5.2 System Graph Blast Radius
- **Modified:** `atom_cip_assigner` (`core/chirality.py`), `atom_chirality_recovery`
  (`core/chirality.py`), `atom_molassembler_adapter`
  (`generation/molassembler_adapter.py`); possibly `atom_oin3d_generator`
  (`generation/engine.py`) if verify-and-reflect sits at engine level.
- **Deleted:** `atom_pseudo_atom_strategy`.
- **Verify-unchanged (audit, not edit):** `atom_xyz2mol` (call site `perception_tmc.py:947`),
  `translator.py:35` (`CIPAssigner` call), `atom_oin_sanitizer` (property
  pass-through), `atom_oin_inline_handler` / `atom_oin_parser_gen` (regex parse
  passes `@/@@` through).
- **Untouched:** direct-parser nodes, `atom_oin_parser_oin`, legacy `oin/parser.py`.

### 5.3 API Contracts / Schema
```python
# src/oinsmiles/core/chirality.py

# NEW atom property (string): '_OIN_CIPCode_LP' ∈ {'R','S'} on metal-bound P only.
# Distinct from the existing '_OIN_CIPCode' (metal-present sense) so both views
# are inspectable; recover() keys the Zone-A branch on the _LP property.

class CIPAssigner:
    def assign_all(self, mol: Chem.Mol) -> Chem.Mol: ...
    # + computes _OIN_CIPCode_LP via dummy-metal copy; rdCIPLabeler cross-check → warnings.warn

class ChiralityRecoveryUtility:
    def recover(self, mol: Chem.Mol) -> Chem.Mol: ...
    # Zone-A P branch: keep+verify+flip when _OIN_CIPCode_LP present; else clear.

# DELETED: class PseudoAtomStrategy, const PSEUDO_ATOMIC_NUM (if unreferenced)
```
No new package dependencies (`rdCIPLabeler`, `warnings` are already available).

### 5.4 Dependencies
- RDKit (existing) — `rdCIPLabeler`, `AssignAtomChiralTagsFromStructure`,
  `AssignStereochemistryFrom3D`, `EmbedMolecule`.
- No new third-party dependencies.

---

## 6. Negative Constraints (The "Do NOTs")
- **DO NOT** emit `@`/`@@` on any nitrogen (Zone-A or otherwise) — trivalent
  `[N@]` is RDKit-cleared; keep clearing Zone-A N.
- **DO NOT** emit a tag on a non-stereogenic P (symmetric phosphine): no
  `_OIN_CIPCode_LP` ⇒ no tag. This is the spurious-tag gate; BDPP/BDNN goldens
  must stay byte-identical.
- **DO NOT** revive `PseudoAtomStrategy` / insert `*` wildcards into fragment SMILES.
- **DO NOT** bump the OIN format version or change the grammar. Content-level only.
- **DO NOT** regress carbon `@/@@` behaviour.
- **DO NOT** promote `Rh-RR-DIPAMP-Cl2_oin.txt` to a trusted fixture without HITL
  R,R sign-off (§4).
- **DO NOT** assert on raw 3D coordinates in generation tests — assert on derived
  metal-present CIP (deterministic oracle over non-deterministic ETKDG output).
- **DO NOT** compute or store a fragment-local CIP for a *trivalent* P directly
  (perception fails — spike 3); always use the 4-coordinate dummy-metal copy.

## 7. Risks & Mitigation
- **RISK-1 (label parity, Q1):** encode-side dummy-metal LP-CIP may not equal the
  fragment trivalent CIP under some atom orderings. → **Mitigation:** the
  `recover()` verify-and-flip recomputes and corrects; add a raw-parity unit test
  on the fragment mol BEFORE `recover()` to make any divergence visible, not silent.
- **RISK-2 (generation enforcement placement, Q2):** wrong layer lacks metal bond
  or conformer. → **Mitigation:** Red Team to pin the object; default to the
  assembled-complex mol in the adapter where both are present.
- **RISK-3 (molassembler fallback, Q3):** `from_smiles` may drop trivalent `[P@]`.
  → **Mitigation:** investigate in MiniPRD-B; atom stereopermutator fallback in
  the worker; molassembler is fallback-only so the primary ETKDG path is unaffected.
- **RISK-4 (spurious tags, Q4):** gate leaks a tag onto a symmetric P.
  → **Mitigation:** BDPP/BDNN byte-identical goldens as a hard regression gate.
- **RISK-5 (circular golden):** the DIPAMP golden encodes whatever the code emits.
  → **Mitigation:** Candidate Artifact protocol — HITL R,R sign-off + independent
  `rdCIPLabeler` oracle + literature configuration cross-check before promotion.
- **RISK-6 (non-RDKit consumers):** must accept `[P@]` on trivalent P.
  → **Mitigation:** standard Daylight SMILES; documented, not gated.

## 8. Success Metrics
- `XYZToSMILES().convert()` on `Rh-RR-DIPAMP-Cl2.xyz` emits `@`/`@@` on both P.
- Independent `rdCIPLabeler` CIP-from-3D on the original complex = (R,R);
  regenerated complex re-derives (R,R); `@↔@@` flip inverts both labels.
- XYZ→OIN→XYZ→OIN byte-stable on DIPAMP.
- `uv run python -m unittest discover tests` → OK; BDPP/BDNN + carbon-chirality
  tests unchanged; new encode + generation tests green.
- No spurious tags on symmetric phosphines; user-facing RDKit warning fires only
  on genuine CIP conflict.

---

## 9. HITL Candidate Artifact — reviewer instructions
When MiniPRD-A completes, the reviewer (Thomas) confirms the emitted stereochemistry:
- **File to check:** `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt` (the new OIN
  string) alongside the source geometry `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz`.
- **What to confirm:** both phosphorus stereocenters are **(R,R)** — i.e. the
  `@`/`@@` on `P{0}` and `P{1}` correspond to the R,R configuration of
  (R,R)-DIPAMP. The execution run will print the `rdCIPLabeler` CIP-from-3D labels
  for both P atoms and a WARNING if RDKit perceives any conflict.
- **If correct:** approve promotion of the string to a trusted golden.
- **If incorrect:** provide the corrected structure/configuration; MiniPRD-A
  re-runs against it. (The RDKit warning is advisory — human sign-off is authoritative.)
```
