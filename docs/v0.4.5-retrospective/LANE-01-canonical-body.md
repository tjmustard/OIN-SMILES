# Lane 1 — Canonical ligand body (and canonical perception)

**What this lane was for:** make the ligand-body text inside an emitted OIN string a function of
the ligand's *graph* alone, so two presentations of one complex (rotated, renumbered, or
regenerated) write the same body bytes — the job `oin/compare.py` was already doing at
*comparison* time, moved upstream into the *emit* path.

---

## ELI5

An OIN string contains, for each ligand, a SMILES "body" — a text spelling of that ligand.
RDKit's spelling routine is already deterministic: give it the same molecule object and it always
writes the same text. The problem is one level up. The molecule object is *guessed* from the raw
3D coordinates (which atom pairs are bonded, and which of them are double bonds), and that guess
is not unique — a ring can come out as alternating single/double bonds ("Kekulé") or as an
aromatic ring, and a charged ligand can put its double bond on either of two equivalent places
("resonance forms"). So the same real molecule becomes two different molecule objects, and the
faithful speller writes two different — both correct — strings. This lane made the *guess*
canonical, not just the spelling: reparse each ligand body through a round trip
(text → molecule → text) so both spellings collapse onto one, and separately fix the two places
in the perception code where the guess depended on the order the atoms happened to appear in the
input file.

## The work, visually

```
        XYZ file  ──────────────────────────────────────────────┐
                                                                │
  ┌─────────────────── PERCEPTION (OIN_CANONICAL_PERCEPTION) ────┴────────────────┐
  │                                                                               │
  │   utils/xyz2mol_local.py::AC2BO            utils/xyz2mol.py::lig_checks        │
  │   "returns an arbitrary resonance form"    ranks candidates by                │
  │   ├─ _ordered_valences → itertools.product │   (most aromatic, fewest charges) │
  │   │     ⚠ walks per-atom lists in INPUT    │   and settled a tie by "whichever │
  │   │       index order                      │   came FIRST" = ResonanceMol-     │
  │   └─ get_UA_pairs → nx.max_weight_matching │   Supplier order = INPUT ORDER    │
  │         ⚠ result depends on edge insertion │                                   │
  │           order                            │  FIX: sort on                     │
  │                                            │  (-N_aromatic, N_pos+N_neg,       │
  │   FIX (by conjugation):                    │   canonical SMILES)               │
  │     relabel atoms into canonical order ──► perceive there ──► map BO matrix    │
  │     back.  Canonical order = the SMILES WRITE order (_smilesAtomOutputOrder),  │
  │     NOT CanonicalRankAtoms(breakTies=True).                                    │
  └───────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼   sanitized fragment mol + donor_indices
  ┌────────────────── EMIT (OIN_CANONICAL_BODY) — oin/canonical_body.py ──────────┐
  │  canonical_body_emit(mol, donor_indices)                                       │
  │                                                                                │
  │   pass k:  stamp map numbers 1..n on the donors  (identity CARRIED, not        │
  │            re-derived — GetSubstructMatch can hand back a wrong automorphism)  │
  │              │                                                                 │
  │              ▼  h_faithful_smiles(canonical=True)                              │
  │            "[N:1]CC[N:2]"  ──► compare._parse_fragment ──► reparsed mol        │
  │              │                                                                 │
  │              ├─ G1  composition: element counts, total H, total charge equal?  │
  │              │        (radical electrons deliberately EXCLUDED)          ✗→None│
  │              ├─ G1b aromatic-atom count did not DROP?                    ✗→None│
  │              ├─ G2  each map number unique, present, same element +           │
  │              │        same heavy degree as the donor it claims?           ✗→None│
  │              ├─ clear chelate-locked E/Z (dummy-Fe probe, rings through it)   │
  │              ├─ clear map numbers, then FINAL h_faithful_smiles write         │
  │              ▼                                                                 │
  │       smiles == previous pass?  ── no ──► loop (≤ _MAX_REPARSE_PASSES = 4)     │
  │              │ yes                       └─ still moving after 4 ⇒ OSCILLATOR  │
  │              ▼                              ⇒ return None (keep input body)    │
  │   (smiles, {donor_idx_in_mol: position_in_smiles}, reparsed_mol)               │
  └────────────────────────────────────────────────────────────────────────────────┘
                                     │  None ⇒ caller keeps its body for the
                                     │         WHOLE fragment (never partially)
                                     ▼
                     OINInlineHandler stamps {n} markers  ──►  OIN string
```

Legend: `✗→None` = guard trips, whole fragment reverts. `⚠` = the order-dependent seam that was
the actual defect. Everything inside the two boxes is behind a lever; both are default-**ON**
since v0.4.5.

## Initial assumptions and hypothesis

1. **The handoff's premise.** The v0.4.4 capstone bucket report labelled 500 molecules
   `rdkit_canonical` inside the 12.32% `key_equal` population, and those were taken to be ligand
   *body* drift — two serializations of one graph. Lane 1 was therefore scoped as the big lane:
   close ~500 molecules.
2. **The fix was believed mechanical.** `oin/compare.py::canonical_fragment_body` already folded
   exactly this drift at compare time — that is *why* the round-trip key was canonical while the
   emitted string was not. Promoting the same function to the emit path should be a plumbing job.
3. **The one thing known to be hard** was slot identity: `{n}` markers are placed from
   `_smilesAtomOutputOrder`, and a reparse changes that order, so the marker had to survive the
   round trip somehow.

## What was actually found

**Refuted — the 500 are mostly not body drift.** `tools/diagnose_body_drift.py` attributes each
row to exactly one root cause. Over the 500 `rdkit_canonical` rows:

| cause | count | share |
|---|---:|---:|
| `slot_or_order` — body multisets **identical**, only slots/order differ (⇒ Lane 2) | **396** | 79.20% |
| `reparse_fixable` — aromatic-vs-Kekulé / implicit-vs-explicit H (⇒ Lane 1) | **104** | 20.80% |
| `ez_chelate_locked`, `resonance_charge`, `formal_charge_placement`, `connectivity`, `unattributed` | 0 each | 0.00% |

Mechanism: `rdkit_canonical` is the **fallthrough** branch of `_key_equal_subclass`'s cascade. Its
three positive tests (`fragment_reorder`, `slot_renumber`, `winding_star_drift`) each require the
drift to be *exclusively* of one kind, so a pair mixing a renumbered slot with a reordered
fragment matches none of them. The headline was an upper bound **off by about 5×**. Over all 828
`key_equal` rows the body-drift population is **104 (12.56%)**, of which **101** would go
byte-exact from the body fix alone and 3 also need Lane 2.

**Refuted — the reparse fixes key-level defects.** `tools/canonicality_probe.py` holds the
molecular graph fixed and varies only presentation (random proper rotation, random atom
renumbering, both), so the correct answer is byte-**identical**. 250 dataset molecules (248
encoded), 3 trials × 3 transforms:

| arm | byte-stable | drifted | `rdkit_canonical` | `slot_renumber` | key broken |
|---|---:|---:|---:|---:|---:|
| levers off | 146 (58.87%) | 102 | 71 | 42 | **49** |
| `OIN_CANONICAL_BODY` | 152 (61.29%) | 96 | 65 | 39 | **49** |
| both levers | 156 (62.90%) | 92 | 58 | 40 | **41** |

Per molecule, against levers-off: the body lever is **6 fixed, 0 regressed, 0 key defects fixed,
0 broken**; both levers together are **10 fixed, 0 regressed, 9 key fixed, 1 broken (net −8)**.
An earlier reported read of the same shape was **drift 17 → 17** for the reparse alone
(`docs/V045_STATUS_2026-07-25.md`). So the reparse is strictly non-regressive and strictly
insufficient: **the perception lever is what repairs the key.**

**Confirmed — and this is the diagnosis that redirected the lane.** `rotate` alone causes **zero**
drift; the encoder was already orientation-invariant. **All** drift comes from renumbering. That
is not a disappointment but an explanation: the reparse folds two *serializations* of one graph,
whereas renumbering hands the serializer a genuinely *different graph*, and `MolToSmiles` is
faithful to whichever resonance form it receives.

**Confirmed — the canonical order is the SMILES write order, not `CanonicalRankAtoms`.** Measured
over 20 random renumberings of `CC(N)=NC`:

| source of order | invariant under renumbering? |
|---|---|
| `Chem.CanonicalRankAtoms(mol, breakTies=True)` | **no** — a different ranking 18 times in 20 |
| `Chem.MolToSmiles(mol)` canonical string | **yes** — one string every time |
| the graph induced by relabelling with `_smilesAtomOutputOrder` | **yes** — one graph every time |

Refined by independent reproduction: the *write order* is not invariant either (20 renumberings →
20 distinct orders). What is invariant is the **graph that order induces** — adjacency bytes plus
symbol vector: always exactly 1 distinct value. The residual freedom is exactly the automorphism
group, which preserves the adjacency matrix, which is all `_AC2BO_core` consumes; instrumenting
base and renumbered inputs shows byte-identical `(AC, atoms)` reaching the core.

**Confirmed — the fixed point exists for almost everything.** Over 6062 distinct capstone corpus
ligand bodies, **6056 converge in one pass**. The other **6 oscillate with period two**: RDKit
flips `@`/`@@` on adamantane-cage carbons across every parse/write cycle (a degenerate cage
"stereocentre" that is not one), so there is no fixed point at all. Those bail. The same
instability already affects `compare.canonical_fragment_body`, so it is pre-existing in the
comparison key, not introduced here.

**Measured cost of the aromaticity guard: zero.** 0 of the 713 distinct ligand bodies in the
capstone `rdkit_canonical` population lose aromaticity under the reparse, and all 500 rows still
converge.

**Regression fixture** `tests/fixtures/NAXDOI.xyz` — the smallest in-repo structure reproducing
the defect. Permuting its atom lines:

| levers | drift over 6 presentations | key stable |
|---|---:|---|
| none | 3/6 (`rdkit_canonical`) | **no** |
| `OIN_CANONICAL_BODY` only | 3/6 — unchanged | **no** |
| `OIN_CANONICAL_PERCEPTION` | **0/6** | yes |

## What was done

**New module `src/oinsmiles/oin/canonical_body.py`.**

- `canonical_body(body_smiles) -> str` — re-export of `compare.canonical_fragment_body` under an
  encoder-side name. Deliberately *not* a copy: one implementation, three consumers (comparison
  key, encoder, Lane 2's vertex colours), so they cannot drift apart. Returns `"RAW:<input>"` for
  an unparseable body.
- `canonical_body_emit(mol, donor_indices) -> (smiles, {donor_idx: position}, reparsed_mol) | None`
  — the encoder seam, hooked in `utils/xyz2mol.py` at the `lever_enabled("OIN_CANONICAL_BODY")`
  branch (~line 1746) that overwrites `sanitized_smiles` / `sanitized_mol` and sets
  `canonical_body_positions`.
- Helpers: `_reparse_once`, `_composition`, `_n_aromatic`, `_clear_chelate_locked_stereo`,
  `_heavy_degree`, `_output_order`. Constants `_MAX_REPARSE_PASSES = 4`, `_MAX_PROBE_DONORS = 10`.

**Design choices, and the alternatives rejected.**

- **Slot identity is carried through the reparse by ATOM MAP NUMBER, never re-derived.** The
  obvious alternative — find the donor again in the reparsed mol with `GetSubstructMatch` — is
  unsafe: on a near-symmetric ligand the matcher can return a *wrong automorphism*, e.g. the two
  ortho carbons of a cyclometalated aryl look identical to it, so the marker lands on a CH instead
  of the deprotonated X-type carbon (`c{N}` → `[cH]{N}`) and the coordination sphere is silently
  corrupted. Map numbers survive `MolToSmiles` → `MolFromSmiles`, so the donors are recovered by
  label, the labels are cleared, and only then is the final canonical SMILES written and its
  output order read.
- **All-or-nothing.** Every failure path returns `None` and the caller keeps the un-reparsed body
  for the **whole** fragment. A misplaced `{n}` is far worse than the notation drift being closed.
- **Three load-bearing guards.** (1) *Composition* — a map number forces brackets and a bracket
  changes implicit-H semantics (`n` in a five-ring means one implicit H; `[n:1]` means none), so
  element counts, total H and total charge must be unchanged. Radical electrons are excluded on
  purpose: RDKit re-derives them from valence on every parse, so comparing them would reject every
  carbene while telling us nothing. (2) *Donor identity* — the atom recovered by map number `k+1`
  must match in element and heavy degree, with no lost or duplicated label. (3) *Idempotence* —
  the emitted body must satisfy `canonical_body(body) == body`, i.e. be a canonical
  representative rather than one step along a walk.
- **Aromaticity guard (1b).** A metallo-porphyrin's macrocycle is perceived aromatic in the
  complex, but with the metal stripped its four N's carry no hydrogen and RDKit's default model
  calls the free base non-aromatic, so a naive reparse emitted `C1=C2C=CC(=N2)…` where the encoder
  had `c1c2nc(…`. Harmless at compare time (both sides lose it identically), a fidelity loss at
  emit time. The reparse now bails when the aromatic-atom count drops.
- **Chelate-locked E/Z clearing** mirrors `compare._chelate_locked_fragment_key` exactly: build a
  probe with a dummy Fe bonded to every donor, take the rings through it, clear `BondStereo` plus
  incident `BondDir` on the double bonds those rings lock. A double bond in *no* metal ring stays
  untouched, so real diastereomers still separate.
- **Perception, seam 1 (`lig_checks`).** Every consumer of its candidate list ranked resonance
  forms by *(most aromatic, fewest formal charges)* and settled ties by "whichever came first" —
  i.e. `ResonanceMolSupplier` enumeration order, a function of the input numbering. Sorting the
  candidate list on `(-N_aromatic, N_pos + N_neg, canonical SMILES)` fixes `_select_lig_mol`'s
  three accumulate loops and `_rescue_unusable_perception`'s `max` at once, changing no selection
  *logic* — only which member of an exact tie wins.
- **Perception, seam 2 (`AC2BO`).** Closed **by conjugation**: relabel the atoms into canonical
  order, perceive there, map the bond-order matrix back. Chosen over hardening the valence walk
  and the Kekulé matching separately, because it turns every index-order dependence inside the
  core into a function of one canonical labelling, in one place.
- **`get_tmc_mol` retry.** Reordering the valence walk can surface a different but equally *valid*
  Lewis structure, and "valid" to `AC2BO` is not "usable" once the dative bonds go on: on
  `AGUFEN.xyz` (a PPN counter-cation) the canonical order draws a `P=c` ylide with a pentavalent
  ipso carbon that passes the free-ligand check and raises `OINEncodeError` at assembly.
  `get_tmc_mol` now retries once under `suppress_canonical_perception()`. Molecules taking the
  retry stay order-dependent — a right answer that drifts beats a reproducible wrong one.

**v0.4.6 aftermath (consistency fix only, no accuracy claim).** `canonical_body_emit` had **two**
`MolToSmiles` writes, and the **final** one — whose output becomes the emitted body — silently
discarded the `OIN_H_FAITHFUL` repair applied upstream at `utils/xyz2mol.py`. Both now route
through `oin/hydrogen.py::h_faithful_smiles`. **This has no measured accuracy benefit:** an A/B
over the 45-molecule `Atom count mismatch` population gave **match 8 / mismatch 37 in both arms,
identical**. It is kept because two writes of one body should not disagree, and because it makes
the composition guard reachable in the useful direction (a fragment whose bare 0-H symbol used to
re-read one hydrogen heavier failed `_composition` and lost canonicalization entirely).

## Dead ends and refutations

| tried | what killed it |
|---|---|
| "the 500 `rdkit_canonical` molecules are body drift" | `tools/diagnose_body_drift.py`: 396/500 (79.2%) are `slot_or_order`; only 104 are reparse-fixable |
| "promoting the compare-layer reparse closes the gap" | canonicality probe: 6 byte-stability fixes, **0 key defects fixed**, drift 17→17 in the earlier read |
| re-derive the donor's new position with `GetSubstructMatch` | rejected on the measured `c{N}` → `[cH]{N}` failure mode: a wrong automorphism on a near-symmetric cyclometalated aryl puts the marker on a CH |
| build the canonical relabelling on `CanonicalRankAtoms(breakTies=True)` | not invariant — a different ranking in **18 of 20** renumberings of `CC(N)=NC`; `breakTies=True` settles symmetry ties on the *input index*. A single-fixture guard would have passed anyway (the Y2 lesson) |
| naive reparse with no aromaticity guard | broke `tests/unit/test_aromatic_reencode.py` on the Ni-porphyrin: `c1c2nc(…` → `C1=C2C=CC(=N2)…` |
| perception retry chosen by "keep whichever perception scored a higher total bond order" | silently re-imported the order dependence the lever removes (the input-order result *is* a function of the numbering). Caught only because it broke the `NAXDOI` invariance guard. The retry is now triggered by the canonical attempt **alone** |
| iterate the reparse until it stops moving, unconditionally | 6 of 6062 corpus bodies oscillate with period two (RDKit adamantyl `@`/`@@` flip) — no fixed point exists, so `_MAX_REPARSE_PASSES` exhaustion returns `None` |
| **(v0.4.6)** rescue Lane 6's P3 descriptor by copying the metal-locked chiral tag onto the reparsed donor | **MEASURED WRONG.** It does make P3 emit under `OIN_CANONICAL_BODY` and POJJOP passes — but setting a tag *after* the sanitize adds a stereocentre the canonical ranker did not know about, which moves the canonical **write order**, and `@`/`@@` is a parity relative to that order. On `RIFGUJ_comp_2` (three Cu-bound amines on one cyclohexane) the three ring-**carbon** tags then flip between a structure and its mirror. Geometry says they must not: `AssignStereochemistryFrom3D` + `rdCIPLabeler` label those carbons lowercase `s` — pseudo-asymmetric, a **relative** all-cis descriptor — identically for the structure and its reflection |
| **(v0.4.6)** "H-faithful in the canonical body will move the atom-count population" | A/B over the 45 `Atom count mismatch` molecules: **8/37 in both arms**. Kept as a consistency fix only |
| tests spelling "lever off" by **deleting** the env var | correct only while a lever defaults OFF; after promotion it silently means ON. **17 test failures in v0.4.5**, 6 more in v0.4.6. Now a lint |

## Where it landed

- **Levers:** `OIN_CANONICAL_BODY` and `OIN_CANONICAL_PERCEPTION`, both **default-ON** since
  v0.4.5 via `src/oinsmiles/oin/levers.py::_DEFAULT_ON`. Promotion evidence:
  `docs/PROMOTION_GATE_v0.4.5.md` — all six canonicality levers together took byte-stability under
  rotation/renumbering from **58.1% (173/298) to 69.6% (208/299)** and comparison-key instability
  from **60 molecules to 16** on a 300-molecule seed-42 sample; `rotate` drift 0 in both arms.
- **Code:** `src/oinsmiles/oin/canonical_body.py`; hooks in `src/oinsmiles/utils/xyz2mol.py`
  (`lig_checks` sort ~line 497, `AC2BO` conjugation ~line 850, `canonical_body_emit` ~line 1746);
  `src/oinsmiles/utils/xyz2mol_local.py::AC2BO`.
- **Guards:** `tests/unit/test_canonical_body.py` —
  `TestFlagOffIsByteIdentical::{test_goldens_unchanged, test_module_not_imported_when_flag_unset}`,
  `TestCanonicalBodyFunction::{test_kekule_and_aromatic_converge, test_is_idempotent_on_its_own_output, test_unparseable_body_gets_stable_raw_token}`,
  `TestCanonicalBodyEmitGuards::{test_returns_fixed_point_of_canonical_body, test_no_atom_map_residue_in_emitted_body, test_preexisting_map_number_bails, test_out_of_range_donor_bails, test_donor_positions_index_the_donor_element, test_oscillating_body_bails_rather_than_emitting_a_non_fixed_point, test_chelate_locked_double_bond_loses_its_ez_marker, test_double_bond_only_half_inside_the_metal_ring_is_untouched, test_pendant_double_bond_keeps_its_ez_marker}`,
  `TestFlagOnEncoderIntegration::{test_markers_stay_on_the_same_donors, test_isomer_is_preserved, test_facmer_stays_distinct, test_every_fixture_still_encodes}`,
  `TestCanonicalPerception::{test_atom_permutation_is_renumbering_invariant, test_naxdoi_drifts_under_renumbering_with_the_lever_off, test_naxdoi_is_renumbering_invariant_with_the_lever_on, test_lever_off_is_byte_identical, test_lever_on_preserves_the_isomer_on_the_goldens}`.
  Also `tests/unit/test_aromatic_reencode.py`, `tests/unit/test_regression_stability.py`,
  `tests/unit/test_locked_donor.py::TestRifgujRingCarbonsArePseudoAsymmetric`, and the trap lint
  `tests/unit/test_levers.py::TestNoTestUnsetsAPromotedLever::test_no_test_file_unsets_a_default_on_lever`.
- **Commits.** Branch `swimlane/v045-lane1`, tip `12569f03`, fully merged into `main`:
  `19d20042` (generator-free canonicality A/B instrument), `7b85e123` (frozen sweep-cohort
  builder), `20044883` (rotation/renumbering canonicality probe), `b7955bad` (Step-0 attribution +
  `OIN_CANONICAL_BODY`), `2ebb9935` (`OIN_CANONICAL_PERCEPTION`), `55002bf8` (results note +
  tie-break sentinel), `c616aa19` (the two safety fixes), `12569f03` (docs). Promotion:
  `1450b5ce`, release `0d165845`. v0.4.6 aftermath: `2aa728f5` (H-faithful writes + P3 negative
  result), `0bf35884` (record that the H-faithful fix has no measured benefit).
- **Recorded incompatibility:** `levers.py::_HELD_OFF["OIN_EMIT_LOCKED_DONOR"]` — P3 works only
  with `OIN_CANONICAL_BODY=0` today, which is how `tests/unit/test_locked_donor.py` runs.
- **Known residue:** `YOYBIY_comp_0` is the one molecule the perception lever makes worse (its
  canonical valence walk lands on an all-single-bond perception of a bis(pyridylamidine) ligand,
  `C1[CH][CH][CH][CH]N1` where other numberings give `c1ccccn1`); it is *usable*, so the
  `get_tmc_mol` retry does not fire. Net key effect of the lever is still −8.

## Open questions / for the next agent

1. **P3 under the shipped default.** The correct fix must preserve the metal-locked donor tag
   *without* perturbing the canonical ranking — either keep the donor bracketed through the
   sanitize, or re-derive parity from the parent geometry *after* the write order is fixed. **Next
   measurement:** with `OIN_EMIT_LOCKED_DONOR=1 OIN_CANONICAL_BODY=1`, assert POJJOP's enantiomer
   pair diverges **and** `RIFGUJ_comp_2`'s three ring-carbon tags are identical for structure and
   reflection (`test_locked_donor.py::TestMultiCentreDescriptor::test_flips_under_reflection`,
   `::test_three_locked_amines_all_invert_together`,
   `::TestRifgujRingCarbonsArePseudoAsymmetric`). Do not accept a fix that only passes POJJOP.
2. **`YOYBIY_comp_0`.** Decide whether the perception ranking should prefer the more-aromatic
   *assembled* perception rather than the more-aromatic free-ligand perception. **Next
   measurement:** re-run `tools/canonicality_probe.py` restricted to the amidine/amidinate
   population and count key-breaks, so a ranking change is judged on more than this one molecule.
3. **The 6 oscillators.** If RDKit ever stabilizes the adamantyl `@`/`@@` flip, the guard test
   flips to an equality and 6 more bodies canonicalize. Re-run
   `canonical_body(canonical_body(cage)) == canonical_body(cage)` on the pinned cage SMILES after
   any rdkit bump (currently pinned `==2025.9.3`).
4. **Stereo perception is still order-dependent** and is now the dominant residual with both
   levers on: `QUPWUT` `[S@]` → `[S@@]`, `OJOXAM` `[C@@H]` → `[C@H]`. Out of Lane 1's scope; owned
   by the stereo-stability lane (`OIN_STABLE_STEREO`).
5. **`OIN_H_FAITHFUL` next step is per-ATOM provenance, not another aggregate.** Two aggregate
   hypotheses are already refuted (`docs/V046_HFAITHFUL_FINDINGS.md`). Walk one molecule from each
   `dH` band (+1…+3 for 28/45; 0 for 4/45; the −14 / −16 / −36 losses) mapping parent atom →
   fragment atom → emitted token.
6. **`OIN_EMIT_AXIAL` must be re-measured under `OIN_CANONICAL_PERCEPTION=1`.** The Y2 cohort
   numbers backing its promotion evidence were taken with perception OFF, and perception feeds
   `_is_atropisomer_candidate`, whose steric wall keys off `not GetIsAromatic()` — measured on
   `YESKOZ`, hindered axes go 2 → 1 under canonical perception.
