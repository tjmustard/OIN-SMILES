# SuperPRD — Stereo Roadmap Phase 4: Zone-A P Stereocenter Encoding

## Metadata
- **Project Name**: OIN-SMILES — Metal-bound (Zone-A) P stereocenter encoding, Phase 4
- **Version**: 1.1.0 (Compiled — MiniPRD-C resolved via `/hyper-resolve`, 2026-07-03; adds SPL/bidentate dummy-embed, absorbs TASK-31)
- **Status**: Ready for `/hyper-execute`
- **Owner**: Architect Agent (Fable) / Thomas Mustard
- **Roadmap ref**: `spec/worklog/ROADMAP-stereo.md` § Phase 4 (absorbs superseded Phase 2)
- **Decision ref**: `spec/worklog/PHASE4-decision.md` (CTO consult, 2026-07-03)
- **Design brief**: `spec/worklog/PHASE4-design-brief.md`
- **Depends on**: TASK-20 (diagnostic, DONE 2026-07-03); fixture `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz`
- **Provenance**: Draft_PRD.md (v0.1.0) → RedTeam_Report_ZoneA_P_Encoding.md (2026-07-03) → this compilation. All ten Red Team findings (B1–B10) carry a documented decision — see §5.1 Resolved Trade-offs Log.
- **Child MiniPRDs**:
  - `spec/compiled/MiniPRD_ZoneA_P_Encode.md` (MiniPRD-A, encode side)
  - `spec/compiled/MiniPRD_ZoneA_P_GenEnforce.md` (MiniPRD-B, generation side)
  - `spec/compiled/MiniPRD_ZoneA_P_SPL_DummyEmbed.md` (MiniPRD-C, generation-side
    SPL/bidentate dummy-embed — fixes the one-sided-SPL bug MiniPRD-B warns on;
    seeded by `spec/worklog/SPL-P-enforcement-decision.md` Option E; Red Team +
    `/hyper-resolve` 2026-07-03; see §5.5 Resolved Trade-offs Log C1–C5)

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
the fragment SMILES by OIN's construction — has exactly 3 fragment-local
neighbours in the monometallic case, so it is always Zone A and always stripped.
The strip happens BEFORE the branch that would consume the stored CIP code.

**The information is discarded, not lost.** `CIPAssigner` computes CIP from the
intact 3D structure (metal present) and attaches it to the atom; `recover()`
throws it away. The fix has ground truth to work from.

This breaks the project's core **lossless round-trip** promise for an entire,
important ligand family (chiral phosphines — DIPAMP-type P-stereogenic
diphosphines used in asymmetric catalysis).

### 1.2 Solution Overview
Per `PHASE4-decision.md` (empirically validated by three RDKit spikes) and the
`/hyper-resolve` session of 2026-07-03:

**Option A — lone-pair `[P@]`/`[P@@]` convention, phosphorus only.** Trivalent
`[P@]` is Daylight-legal, parses to a stable CIP label, survives sanitize /
AddHs / RemoveHs / 20 randomized renumber round-trips, and ETKDG builds the
correct pyramid from it with zero embed changes. No OIN grammar change and no
version bump — `@`/`@@` inside fragment SMILES is existing v3.x syntax; Zone-A P
tags are new *content*, not new *syntax*.

- **MiniPRD-A (encode side):** `CIPAssigner.assign_all()` stores a fragment-local
  (lone-pair-convention) CIP label `_OIN_CIPCode_LP` on each **eligible** metal-bound
  P (bonded to exactly one metal — B7), computed while 3D is available via the
  **dummy-metal equivalence** (§5.1). All dummy-copy construction is guarded:
  any failure degrades to store-nothing + `OINStereoWarning` (B4).
  `ChiralityRecoveryUtility.recover()` replaces the Zone-A *clear* with a
  *verify-and-flip* keyed on `_OIN_CIPCode_LP`, checked **before** the degree
  branches (B9), recomputing CIP after each flip for multi-P fragments (B6).
  `rdCIPLabeler` is the authoritative CIP implementation end-to-end (B5); the
  diagnostic oracle runs on the **same dummy-metal copy** — never cross-convention
  (B1). Delete the dead `PseudoAtomStrategy` (Option B, formally rejected).
- **MiniPRD-B (generation side):** parse already passes `@/@@` through untouched;
  ETKDG builds the pyramid; add a **post-assembly verify-and-re-embed** on the
  fully-assembled complex (P is 4-coordinate → 3D perception valid): on label
  mismatch, re-run ETKDG on the fragment with a new seed, bounded at 3 attempts —
  **no mirror/improper transforms, ever** (B2/B3). Persistent mismatch emits the
  structure with a hard `OINStereoWarning`. Fallback paths with no assembled mol
  (eta/DG fallback) skip enforcement + warn (B8, residual RISK-9).

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
**Confidence Score**: 9/10 (post-resolution)

- Encode side: 9/10 — three spikes behind the tag mechanics; the dummy-copy
  recipe and CIP-implementation choice are now pinned (B4, B5); residual
  uncertainty is the raw-parity question (Q1), absorbed by verify-and-flip and
  made visible by a mandated raw-parity unit test.
- Generation side: 8/10 (explicit deduction per Red Team §2) — the re-embed
  enforcement path has zero spikes behind it; mitigated by the bounded-retry
  NFR and the fact that spike 2 shows ETKDG mis-embeds are rare, so the branch
  is a safety net, not a hot path. Q3 (molassembler `from_smiles` vs trivalent
  `[P@]`) remains an in-MiniPRD investigation task.

All four Draft-PRD open questions are resolved or converted into MiniPRD tasks:
- **Q1 (label parity)** → verify-and-flip absorbs divergence; raw-parity unit
  test mandated (RISK-1).
- **Q2 (enforcement placement)** → pinned to the assembled-complex mol in the
  adapter (`_template_generate`); engine level rejected (`GeneratedStructure.mol`
  is `Optional` and post-hoc).
- **Q3 (molassembler fallback)** → investigation task in MiniPRD-B; stereopermutator
  fallback in `_molassembler_worker` (`molassembler_adapter.py:1489`).
- **Q4 (spurious-tag gate)** → BDPP/BDNN byte-identical goldens as hard gate,
  plus an explicit emitted-SMILES tag-absence assertion and the `assign_all()`
  idempotence constraint (clear stale `_LP` before recompute).

---

## 3. Scope

### 3.1 In-Scope
**MiniPRD-A — Encode side (`core/chirality.py`):**
1. **Zone-A eligibility (formal, shared by both MiniPRDs — B7):** a P atom bonded
   to **exactly one** metal atom, where "metal" is decided by the single existing
   predicate source (`TRANSITION_METALS_NUM`, `utils/perception_tmc.py:20`) imported —
   never duplicated. P bonded to ≥2 metals (bridging phosphide): store nothing +
   `OINStereoWarning` → degrades to today's clearing behaviour.
2. `CIPAssigner.assign_all()`: for each eligible P, build the **dummy-metal copy**
   (§5.1 recipe: metal→Z=0, zero formal charge, clear isotope, drop all other
   metal bonds and all other metals' bonds, sanitize inside try/except), compute
   the lone-pair-convention CIP with **`rdCIPLabeler`** (B5), store as atom prop
   `_OIN_CIPCode_LP` ∈ {`R`,`S`}. Any failure (construction, sanitize, no label)
   ⇒ store nothing (+ warn on exception paths) — the spurious-tag gate + B4.
   `assign_all()` first **clears any pre-existing `_OIN_CIPCode_LP`** (idempotence, B9).
3. `ChiralityRecoveryUtility.recover()`: new branch order (B9) —
   `_OIN_CIPCode_LP` present ⇒ Zone-A verify-and-flip **regardless of degree**
   (keep tag, recompute fragment-local CIP via `rdCIPLabeler`, flip on mismatch);
   only then fall through to the existing degree-keyed branches. When a fragment
   holds >1 tagged P, **recompute CIP after each flip** (bounded fixed-point,
   max 2 full passes — B6). P/N without the property, and **all** N, keep today's
   clearing behaviour.
4. Delete `PseudoAtomStrategy` (class + `atom_pseudo_atom_strategy` node) and the
   `architecture.yml` dangling edge `atom_chirality_recovery → atom_pseudo_atom_strategy`
   (B10); remove `PSEUDO_ATOMIC_NUM` if unreferenced; scrub docstring/comment
   references (`tests/unit/test_stereo_roundtrip_diagnostics.py:181`, `recover()`
   final-else comment). The ≥4-neighbours-no-CIP clearing *behaviour* survives —
   only the naming/scaffolding goes.
5. **Diagnostic warning (B1-resolved):** during `CIPAssigner`, run `rdCIPLabeler`
   **on the dummy-metal copy** and warn (`OINStereoWarning`, atom index in the
   message string) if it conflicts with the stored `_OIN_CIPCode_LP` or perception
   errors where a tag was set. Metal-present labels are printed for HITL but
   **never compared cross-convention**. The call is guarded (try/except) and
   skippable via a flag for batch users.
6. Goldens: regenerate `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt` (§9 —
   Candidate Artifact, HITL sign-off required).

**MiniPRD-B — Generation side (`generation/molassembler_adapter.py`):**
7. Post-assembly **verify-and-re-embed** in `_template_generate`'s assembled-complex
   stage (the only layer holding both the metal bond and per-fragment atom
   provenance before positions freeze): recompute the lone-pair-convention label
   via the dummy-metal recipe (`rdCIPLabeler`) and compare to the input fragment
   tag; on mismatch, re-run ETKDG on the fragment with a new seed and re-place.
   **Max 3 attempts (B8); no mirror or improper transform ever (B2/B3).**
   Persistent mismatch ⇒ emit the structure + hard `OINStereoWarning` (never die
   in the 60 s ProcessPoolExecutor timeout). Paths with no assembled RDKit mol
   (eta fallback, DG fallback) skip enforcement + warn (residual RISK-9).
8. Molassembler fallback: investigate Q3; if `masm.io.experimental.from_smiles`
   drops trivalent `[P@]`, set an atom stereopermutator in `_molassembler_worker`
   (`molassembler_adapter.py:1489`).
9. Tests: enantiomer-discrimination (flip `@↔@@` in DIPAMP OIN → opposite
   metal-present CIP on regenerated structures, fixed ETKDG seed); full
   XYZ→OIN→XYZ→OIN byte-stable round-trip on DIPAMP (byte-stability scoped to
   the pinned RDKit version — `uv.lock`; an RDKit bump regenerates goldens via
   the Candidate Artifact protocol); a forced-mismatch test for the re-embed
   branch via a test-only injection point (mocked/pre-mirrored initial embed).

**Cross-cutting:**
10. New warning class `OINStereoWarning(UserWarning)`; atom index always in the
    message string (defeats warning dedup); tests use
    `warnings.catch_warnings(record=True)`; clean-fixture suite must pass with
    `-W error::OINStereoWarning`.
11. Unit test: `parse_inline_string` with `[P@]{0}` and `[P@@]{1>}` (bracket-atom
    token + slot + winding co-occurrence) — closes the `SLOT_REGEX`-adjacency
    question (B10).

### 3.2 Out-of-Scope (explicit)
- **Nitrogen Zone-A stereo.** Trivalent `[N@]` is RDKit-cleared; deferred to a
  future Option-C out-of-band marker. `recover()` keeps clearing Zone-A N.
- **Option B (wildcard `*` pseudo-atom)** — formally rejected; scaffolding deleted.
- **Option C (OIN grammar / out-of-band marker)** — not built; reserved for N.
- **OIN version bump** — none; content-level within existing v3.x grammar.
- **Carbon `@/@@` paths** — already correct; must not regress (TASK-10 set is
  the named gate, run for **both** MiniPRDs — MiniPRD-B touches the adapter every
  carbon path uses).
- **Mirror/reflection machinery** — rejected for this feature (B2/B3); if Phase 3
  builds a reflection utility, this feature does not consume it.
- **Axial/atropisomeric P (BINAP)** — separate concern (`test_axial_chiral.py`).
- **Winding / haptic face** — Phases 1/3, separate SuperPRDs.
- **Direct-parser nodes** (`parse_oin_direct`) — deferred to v0.2.2.
- **Legacy V2.4 sidecar / `oin/parser.py` path.**
- **Polymetallic/bridging-P label support** — guarded out (store nothing + warn),
  not implemented.

---

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
| :-- | :-- | :-- | :-- |
| US-A1 | As the encode pipeline, I store a fragment-local CIP for each eligible metal-bound P so it survives fragmentation. | 1. After `CIPAssigner.assign_all()` on `Rh-RR-DIPAMP-Cl2` mol, both P atoms carry `_OIN_CIPCode_LP` ∈ {R,S}, computed by `rdCIPLabeler` on the dummy-metal copy.<br>2. Dummy copy per §5.1 recipe (Z→0, charge 0, isotope cleared, single P–dummy bond retained).<br>3. A symmetric phosphine P (two identical phenyls) gets NO `_OIN_CIPCode_LP`.<br>4. P bonded to ≥2 metals: no property + `OINStereoWarning`.<br>5. Dummy-copy sanitize failure (CpM(PR₃) fixture): no property + warning, `convert()` completes. | High |
| US-A2 | As the encode pipeline, `recover()` keeps and corrects the P tag instead of clearing it. | 1. `XYZToSMILES().convert()` on DIPAMP emits `@`/`@@` on **both** `P{0}` and `P{1}`.<br>2. `_OIN_CIPCode_LP` branch runs before any degree check; recompute uses `rdCIPLabeler`.<br>3. Multi-P fragment: CIP recomputed after each flip (≤2 full passes).<br>4. Zone-A P without `_OIN_CIPCode_LP` → cleared; emitted BDPP/BDNN SMILES asserted tag-free (explicit, not only byte-golden).<br>5. All Zone-A N → cleared. 6. Carbon `@/@@` unaffected (TASK-10 green). | High |
| US-A3 | As a maintainer, dead Option-B code is gone. | 1. `PseudoAtomStrategy` class and `atom_pseudo_atom_strategy` node removed; `architecture.yml` dangling edge removed.<br>2. No remaining import/reference, including docstrings and comments.<br>3. The ≥4-no-CIP clearing behaviour is retained under a neutral name. | Medium |
| US-A4 | As the user, I get a warning only on genuine self-inconsistency. | 1. `rdCIPLabeler` re-run **on the dummy-metal copy** conflicting with stored `_OIN_CIPCode_LP` (or erroring where a tag was set) → `OINStereoWarning` naming the atom index in the message.<br>2. No warning on clean DIPAMP; clean-fixture suite passes with `-W error::OINStereoWarning`.<br>3. Diagnostic skippable via flag; metal-present labels printed, never compared. | High |
| US-A5 | As the test suite, the DIPAMP golden is a verified Candidate Artifact. | 1. New OIN string written to `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt`.<br>2. NOT promoted to `tests/fixtures/` trusted golden until HITL confirms (R,R) per §9 protocol.<br>3. BDPP/BDNN goldens byte-identical (negative control).<br>4. Carbon-chirality round-trips (TASK-10) stay green. | High |
| US-B1 | As the generation pipeline, the input P tag is preserved and enforced through embedding. | 1. Regenerating DIPAMP OIN yields metal-present CIP (R,R) on both P (fixed seed).<br>2. Flipping `@↔@@` yields opposite metal-present CIP.<br>3. Forced mis-embed (test injection) → corrected within ≤3 re-embed attempts, no mirror applied.<br>4. Persistent mismatch → structure emitted + `OINStereoWarning` (no timeout death).<br>5. mol=None paths (eta/DG fallback) → enforcement skipped + warning. | High |
| US-B2 | As the round-trip contract, DIPAMP is lossless. | 1. XYZ→OIN→XYZ→OIN byte-stable on the OIN string (pinned RDKit version).<br>2. `rdCIPLabeler` CIP-from-3D on the regenerated metal-present complex matches the original. | High |

---

## 5. Technical Specifications (The Blueprint)

### 5.1 Architecture & Resolved Trade-offs

**The lone-pair label convention (the correctness crux).** CIP labels are
context-dependent: with the metal present it is usually P's highest-priority
substituent; in the fragment the lone pair is lowest. The R/S letter can
legitimately differ between the two views of the same geometry, so exactly ONE
convention is pinned:

**The OIN string stores the fragment-local (lone-pair) CIP sense.** Computable
encode-side (3D still present) via the **dummy-metal equivalence**: a dummy atom
(Z=0) is also lowest priority and sits exactly where the metal was — i.e. where
the lone pair points. Therefore:

> lone-pair CIP of the fragment ≡ CIP of the 4-coordinate P in a copy of the full
> mol where the metal is swapped to Z=0, keeping only its bond to this P.

**Dummy-copy recipe (normative — B4):** swap metal→Z=0, **zero its formal
charge, clear its isotope**, drop all other metal–ligand bonds **and all other
metals' bonds**, then `Chem.SanitizeMol` inside try/except. Failure at any step
⇒ store nothing + `OINStereoWarning` (degrades to today's clearing behaviour;
`convert()` never newly crashes).

Data flow:
```
XYZ → get_tmc_mol() [full mol, 3D, metal present]
  → CIPAssigner.assign_all()
      clear stale _OIN_CIPCode_LP (idempotence)                       [NEW]
      for each P bonded to EXACTLY ONE metal (TRANSITION_METALS_NUM):
        build dummy-metal copy (recipe above, guarded)                [NEW]
        rdCIPLabeler on copy → store _OIN_CIPCode_LP on the real atom [NEW]
        rdCIPLabeler re-check on the SAME copy → warn on conflict     [NEW]
  → fragmentation (metal excluded; tag + properties copied by AddAtom)
  → OINSanitizer.generate_robust_smiles() → (smiles, fragment mol)
  → ChiralityRecoveryUtility.recover()
      _OIN_CIPCode_LP present ⇒ Zone-A verify-and-flip (any degree),
        rdCIPLabeler recompute; multi-P: recompute after each flip    [CHANGED]
      Zone-A P/N without it, and all N: clear (unchanged)
  → OIN string carries [P@]/[P@@]

OIN string → parse_inline_string [regex-only, @/@@ passes through — no change]
  → OINParser.parse → fragment SMILES with [P@]
  → _template_generate: MolFromSmiles(sanitize=False) → AddHs → ETKDG
  → assembled complex [P now 4-coordinate]
  → verify: dummy-metal LP-CIP (rdCIPLabeler) vs input tag
      mismatch → re-embed fragment, new seed, ≤3 attempts             [NEW]
      still mismatched → emit + OINStereoWarning                      [NEW]
      no assembled mol (eta/DG fallback) → skip + warn                [NEW]
```

**Resolved Trade-offs Log** (all decisions: `/hyper-resolve`, 2026-07-03, Thomas Mustard):

- **B1 — Warning oracle convention (Blocker).**
  *Issue:* Draft §3.1.4/US-A4 compared `rdCIPLabeler`-on-the-metal-present-complex
  against `_OIN_CIPCode_LP` — two conventions §5.1 declares legitimately
  divergent; the warning would fire on correct output.
  *Options:* (A) run the oracle on the same dummy-metal copy; (B) two separate
  same-convention checks; (C) drop the auto-warning.
  *Resolution:* **A** — dummy-copy oracle only. Apples-to-apples; warning means
  genuine self-inconsistency. Metal-present labels remain informational (printed
  for HITL, never compared).
- **B2/B3 — Mis-embed correction strategy (Blockers).**
  *Issue:* Draft's "Phase-3-style reflection" doesn't exist (Phase 3 is an
  unresolved parallel draft), and a whole-fragment mirror is an improper
  transform that inverts every co-resident stereocenter (diastereomer
  manufacture; flip-flop loops).
  *Options:* (A) bounded re-embed with new seed, no mirror; (B) build a
  decision-table mirror utility now; (C) block on Phase 3 and share.
  *Resolution:* **A** — re-embed only. Co-resident stereocenters safe by
  construction; zero Phase-3 dependency; cost is rare extra embeds.
- **B4 — Dummy-copy failure behaviour (High).**
  *Issue:* deleting metal–ligand bonds can strand eta-ligand aromaticity →
  `SanitizeMol` throws → unguarded exception crashes `convert()` for complexes
  that work today (worst unlisted regression).
  *Options:* (A) guard + degrade to store-nothing + warn; (B) fail loud.
  *Resolution:* **A** — degrade + warn, full normative recipe, CpM(PR₃)-type
  regression fixture added (RISK-7).
- **B5 — Authoritative CIP implementation (High).**
  *Issue:* legacy `Chem.AssignStereochemistry` (`_CIPCode`) and `rdCIPLabeler`
  disagree on CIP rules 4b/5; mixing them yields a permanent unactionable
  warning state.
  *Options:* (A) `rdCIPLabeler` end-to-end; (B) legacy end-to-end.
  *Resolution:* **A** — `rdCIPLabeler` for every `_OIN_CIPCode_LP` computation,
  recomputation, and cross-check. Legacy still runs for tag *perception*; its
  `_CIPCode` is never compared. Calls guarded per-call + skippable flag
  (pathological-runtime note, RISK-10).
- **B6–B10 + NFR package (Medium/Low).** *Resolution:* standard defaults
  approved as a group — recompute-per-flip (≤2 passes); exactly-one-metal
  eligibility guard with the predicate imported from `utils/perception_tmc.py`
  (no duplicate list); 3-attempt re-embed budget with mol=None paths
  skip+warn; `_LP`-before-degree branch precedence + `assign_all()` idempotence;
  `OINStereoWarning(UserWarning)` with atom index in-message and
  `-W error::OINStereoWarning` clean-fixture gate; housekeeping (dangling
  `architecture.yml` edge, `:1418`→`:1489` line ref, constraint reword,
  `parse_inline_string` bracket-adjacency test); byte-stability scoped to the
  pinned RDKit version.
- **HITL protocol (Red Team §9).** *Resolution:* amended protocol approved —
  see §9.

### 5.2 System Graph Blast Radius
The following nodes in `spec/compiled/architecture.yml` are affected:
- **Modified:** `atom_cip_assigner` (`core/chirality.py`), `atom_chirality_recovery`
  (`core/chirality.py`), `atom_molassembler_adapter`
  (`generation/molassembler_adapter.py`).
- **Deleted:** `atom_pseudo_atom_strategy` **and** the dangling edge
  `atom_chirality_recovery.edges.depends_on → atom_pseudo_atom_strategy`.
- **Verify-unchanged (audit, not edit):** `atom_xyz2mol` (call site
  `perception_tmc.py:947` UNVERIFIED; re-serialize/re-parse at `perception_tmc.py:949/957` UNVERIFIED),
  `translator.py:35` (`CIPAssigner` call), `atom_oin_sanitizer` (property
  pass-through), `atom_oin_inline_handler` — specifically
  `oin/inline.py::parse_inline_string` and `_count_smiles_atoms_before`
  (`oin/inline.py:5`) for `[P@]{0}` bracket-token adjacency —
  and `atom_oin_parser_gen` (regex parse passes `@/@@` through).
- **Untouched:** `generation/engine.py` (enforcement pinned to the adapter, not
  the engine — Q2), direct-parser nodes, `atom_oin_parser_oin`, legacy `oin/parser.py`.

### 5.3 Execution Checklist (MiniPRDs)
Execute in order — A establishes the property contract B consumes:
- [x] `spec/compiled/MiniPRD_ZoneA_P_Encode.md` (MiniPRD-A) — DONE
- [x] `spec/compiled/MiniPRD_ZoneA_P_GenEnforce.md` (MiniPRD-B; HITL golden
      sign-off per §9 may proceed in parallel with B) — DONE
- [ ] `spec/compiled/MiniPRD_ZoneA_P_SPL_DummyEmbed.md` (MiniPRD-C) — resolved,
      ready for `/hyper-execute`. Fixes the one-sided-SPL correctness bug
      MiniPRD-B honestly warns on; extends enforcement to all bidentate Zone-A-P
      including incompatible-bite (absorbs TASK-31). See §5.3.1.

### 5.3.1 Resolved Trade-offs Log — MiniPRD-C (`/hyper-resolve`, 2026-07-03, Thomas Mustard)

MiniPRD-C (`spec/compiled/MiniPRD_ZoneA_P_SPL_DummyEmbed.md`) was Red-Teamed
(`RedTeam_Report.md`, 7 findings: 2 blocker, 2 high, 2 medium, 1 low) and resolved
in this session. Full detail lives in MiniPRD-C §0; the SuperPRD-level record:

- **C1 — Attach helper (Blockers 1 & 2):** the draft's reuse targets are dead —
  `PseudoAtomStrategy`/`PSEUDO_ATOMIC_NUM` were deleted from `src/` (verified:
  `grep` empty; `chirality.py:22` now imports `TRANSITION_METALS_NUM`), and
  `_build_dummy_metal_copy` *raises* on a metal-free fragment and *converts* a
  metal rather than *attaching* a dummy. **Resolution:** a new dedicated
  `_attach_dummy_metal(mol, p_idx)` in `core/chirality.py` — net-new code (not a
  revival), mirroring only the SINGLE-bond + `NoImplicit` valence-fix recipe.
- **C2 — Denticity scope (Finding 3 + expansion):** the gate is also true for
  bidentate Zone-A-P, whose branch had no strip → dummy leak. **Resolution:**
  extend the dummy-embed to **all** Zone-A-P of any denticity, **including
  incompatible-bite chelates**. This **absorbs TASK-31** and **supersedes** commit
  `ee0b3f0`'s "route incompatible-bite Zone-A-P to DG without enforcement": DG
  still handles *placement* for incompatible bites, but the dummy-embed +
  `_verify_zone_a_p` enforcement now runs on that path too. The original
  MiniPRD-C "template-path only" / "DO NOT touch DG-fallback" constraints are
  formally overridden (MiniPRD-C §4 OVERRIDE notes).
- **C3 — Test-1 oracle (Finding 4):** never compare `_metal_present_cip_label`
  against `_zone_a_p_expected_labels` (divergent conventions). **Resolution:**
  enforcement asserts like-for-like `_lp_cip_label` == `_zone_a_p_expected_labels`
  (the pair `_verify_zone_a_p` trusts, B1); discrimination asserts opposite
  `_metal_present_cip_label` for `[P@]` vs `[P@@]`.
- **C4 — NFR/robustness defaults (Findings 5, 6, 7 + Test-5) — approved as a
  group:** hard `_stitch_fragment` postcondition (no Z=0 atom returned;
  `len(positions)` dummy-absent-stable) + acceptance-gate guard "no generated XYZ
  contains a Z=0/`*` atom"; pinned op order `SetNoImplicit(P)` → attach (SINGLE,
  `NoImplicit`) → `AddHs` → embed → align → strip-before-return + pre-placement
  `_lp_cip_label` parity assert (hard failure, not silent loop-mask); Test 7
  loop-with-dummy via `_test_flip_chiral_idx`; branch-entry-counter inertness
  proof (not byte-diff alone).
- **C5 — PAMP fixture gating (Candidate Artifact):** **hard merge gate
  immediately** — Test 4 (byte-stable `PtCl2-PAMP` round trip) blocks now,
  trusting the Avogadro fixture's absolute configuration as-is. *Accepted risk
  (RISK-C3):* round-trip stability pinned to an unverified reference; reviewer
  sign-off recorded post-hoc in `spec/worklog/`, non-blocking.

**Scope-change consequences to reconcile downstream:** flip the CONFIRMED-BUG row
and the TASK-31 row in `spec/worklog/NOTES.md` on MiniPRD-C completion; note that
`ee0b3f0`'s enforcement-free incompatible-bite routing is superseded for Zone-A-P.

### 5.4 API Contracts / Schema
```python
# src/oinsmiles/core/chirality.py

class OINStereoWarning(UserWarning):
    """All Phase-4 stereo diagnostics. Message ALWAYS embeds the atom index."""

# NEW atom property (string): '_OIN_CIPCode_LP' ∈ {'R','S'} — set ONLY on P
# bonded to exactly one metal, ONLY by rdCIPLabeler on the dummy-metal copy.
# Distinct from existing '_OIN_CIPCode' (metal-present sense, informational).
# assign_all() clears any pre-existing _OIN_CIPCode_LP before recomputing.

class CIPAssigner:
    def assign_all(self, mol: Chem.Mol) -> Chem.Mol: ...
    # + dummy-metal copy (normative recipe §5.1, guarded), rdCIPLabeler label,
    #   same-copy cross-check → OINStereoWarning; diagnostic skippable via flag.

class ChiralityRecoveryUtility:
    def recover(self, mol: Chem.Mol) -> Chem.Mol: ...
    # Branch order: _OIN_CIPCode_LP present ⇒ verify-and-flip (rdCIPLabeler),
    # any degree, recompute after each flip (≤2 passes); else existing branches.

# DELETED: class PseudoAtomStrategy, const PSEUDO_ATOMIC_NUM (if unreferenced)

# src/oinsmiles/generation/molassembler_adapter.py
# _template_generate: post-assembly verify (dummy-metal rdCIPLabeler vs input
# tag) → re-embed new seed ≤3 attempts → emit + OINStereoWarning on persistence;
# skip + warn when no assembled mol exists.
```
No new package dependencies (`rdCIPLabeler`, `warnings` already available).

### 5.5 Dependencies
- RDKit (existing, version pinned by `uv.lock`) — `rdCIPLabeler`,
  `AssignAtomChiralTagsFromStructure`, `AssignStereochemistryFrom3D`,
  `EmbedMolecule`.
- `utils/perception_tmc.py::TRANSITION_METALS_NUM` — the single metal predicate source.
- No new third-party dependencies.

---

## 6. Negative Constraints (The "Do NOTs")
- **DO NOT** emit `@`/`@@` on any nitrogen (Zone-A or otherwise); keep clearing Zone-A N.
- **DO NOT** emit a tag on a non-stereogenic P: no `_OIN_CIPCode_LP` ⇒ no tag.
  BDPP/BDNN goldens must stay byte-identical AND their emitted SMILES asserted tag-free.
- **DO NOT** store `_OIN_CIPCode_LP` on a P bonded to ≥2 metals — warn and skip.
- **DO NOT** duplicate the metal predicate — import `TRANSITION_METALS_NUM` from
  `utils/perception_tmc.py` (no second list; TD-005 lesson).
- **DO NOT** apply a mirror or any improper transform to a fragment — correction
  is re-embed-only (B2/B3).
- **DO NOT** compare CIP labels across conventions: the diagnostic oracle runs on
  the dummy-metal copy only; metal-present labels are print-only (B1).
- **DO NOT** use legacy `_CIPCode` values in any `_OIN_CIPCode_LP` computation,
  recomputation, or comparison — `rdCIPLabeler` end-to-end (B5).
- **DO NOT** let dummy-copy construction/sanitize exceptions escape
  `assign_all()` — degrade to store-nothing + warn (B4).
- **DO NOT** revive `PseudoAtomStrategy` / insert `*` wildcards into fragment SMILES.
- **DO NOT** bump the OIN format version or change the grammar. Content-level only.
- **DO NOT** regress carbon `@/@@` behaviour (TASK-10 gate, both MiniPRDs).
- **DO NOT** promote `Rh-RR-DIPAMP-Cl2_oin.txt` to a trusted fixture without HITL
  (R,R) sign-off per §9.
- **DO NOT** assert on raw 3D coordinates in generation tests — assert on derived
  metal-present CIP (deterministic oracle over non-deterministic ETKDG output).
- **DO NOT** derive a chiral tag from **3D perception of a trivalent P**
  (perception fails — spike 3); geometry↔tag work always uses the 4-coordinate
  dummy-metal copy. (Graph-based CIP recompute from an existing tag on trivalent
  P — the `recover()` verify step — is fine and required.)

## 7. Risks & Mitigation
- **RISK-1 (label parity, Q1):** encode-side dummy-metal LP-CIP may not equal the
  fragment trivalent CIP under some atom orderings. → `recover()` verify-and-flip
  corrects; raw-parity unit test on the fragment mol BEFORE `recover()` makes
  divergence visible, not silent.
- **RISK-2 (enforcement placement, Q2 — RESOLVED):** pinned to the
  assembled-complex mol in `_template_generate`; engine level rejected.
- **RISK-3 (molassembler fallback, Q3):** `from_smiles` may drop trivalent
  `[P@]`. → investigate in MiniPRD-B; stereopermutator fallback in the worker;
  molassembler is fallback-only, primary ETKDG path unaffected.
- **RISK-4 (spurious tags, Q4):** gate leaks a tag onto a symmetric P.
  → BDPP/BDNN byte-identical goldens + explicit tag-absence assertion +
  `assign_all()` idempotence (property presence can't outlive provenance).
- **RISK-5 (circular golden):** the DIPAMP golden encodes whatever the code
  emits. → Candidate Artifact protocol (§9): HITL (R,R) sign-off + independent
  `rdCIPLabeler` oracle + literature configuration cross-check (citation/CCDC
  refcode named in MiniPRD-A) before promotion.
- **RISK-6 (non-RDKit consumers):** must accept `[P@]` on trivalent P.
  → standard Daylight SMILES; documented, not gated.
- **RISK-7 (encode-side crash regression — NEW, B4):** dummy-copy sanitize
  failure on eta-ligand/charged complexes. → guarded recipe, degrade to
  store-nothing + warn; CpM(PR₃)-type regression fixture proves `convert()`
  still completes.
- **RISK-8 (enforcement non-termination — NEW, B8):** verify/re-embed loop.
  → hard budget: ≤3 re-embed attempts, then emit + `OINStereoWarning`; never
  reaches the 60 s ProcessPoolExecutor timeout.
- **RISK-9 (unenforced fallback paths — NEW, B8):** eta/DG-fallback generation
  has no assembled mol ⇒ no stereo enforcement. → documented residual risk;
  skip + warn so the gap is observable, never silent.
- **RISK-10 (rdCIPLabeler pathological runtime — NEW, B5):** exponential corner
  cases on highly symmetric fused systems, on the hot encode path. → per-call
  try/except guard + diagnostic skippable via flag for batch users.

## 8. Success Metrics
- `XYZToSMILES().convert()` on `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz` emits
  `@`/`@@` on both P.
- Independent `rdCIPLabeler` CIP-from-3D on the original complex = (R,R);
  regenerated complex re-derives (R,R); `@↔@@` flip inverts both labels
  (fixed ETKDG seed in tests).
- XYZ→OIN→XYZ→OIN byte-stable on DIPAMP (pinned RDKit version).
- `uv run python -m unittest discover tests` → OK; BDPP/BDNN + carbon-chirality
  (TASK-10) tests unchanged; new encode + generation tests green.
- **Clean-fixture suite passes with `-W error::OINStereoWarning`** — "no
  spurious warnings" is an enforced gate, not a sentence.
- **Negative-capability metric (graceful degradation is observable):** for any P
  where the LP label cannot be computed (multi-metal, sanitize failure, no CIP),
  the emitted OIN is byte-identical to today's output and an `OINStereoWarning`
  names the atom.

---

## 9. HITL Candidate Artifact — reviewer instructions (amended per Red Team §9)
When MiniPRD-A completes, the reviewer (Thomas) confirms the emitted stereochemistry:
- **File to check:** `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt` alongside
  the **canonical** source geometry `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz`.
  The duplicate `tests/integration/Rh-RR-DIPAMP-Cl2.xyz` is **deleted** before
  promotion; the golden's provenance line names the exact fixture path + SHA.
- **Convention statement (read first):** the `@`/`@@` in the OIN string encodes
  the **fragment-local (lone-pair) sense**. The printed `rdCIPLabeler` table
  shows the **metal-present sense**. They may legitimately differ. **Confirm
  (R,R) against the metal-present column** — that is the literature convention
  for (R,R)-DIPAMP complexes.
- **What the execution run provides:** a copy-pasteable per-atom `rdCIPLabeler`
  label table (both conventions, clearly headed) **and** a mol-block depiction
  of the complex; plus any `OINStereoWarning` output.
- **Literature cross-check:** MiniPRD-A names the (R,R)-DIPAMP citation/CCDC
  refcode the reviewer checks against — sign-off is not a vibe check.
- **If correct:** approve promotion to a trusted golden. **Record the sign-off
  (who/when/verdict) in `spec/worklog/`.**
- **If incorrect:** provide the corrected structure/configuration; MiniPRD-A
  re-runs against it. (Warnings are advisory — human sign-off is authoritative.)
