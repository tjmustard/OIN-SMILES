# Lane 7 — Y3 research residuals

**What this lane was for:** close the three questions the Y3 injectivity wave left open (an
unsplittable 70-structure ambiguous residue, an unsound donor-swap probe, and a missing `@SP`
distinctness operator) and build the two Lane-5 metal-stereo fixtures that had been in the Y1
plan since the beginning and never existed — all as a **measurement-only** lane that changed no
encoder output.

Primary sources: `docs/INJECTIVITY_Y3_RESIDUALS.md`,
`docs/INJECTIVITY_Y3_UNKNOWN_UNKNOWNS.md`, `docs/INJECTIVITY_Y1_OVERVIEW.md`,
`docs/V045_STATUS_2026-07-25.md` §"Lane 7 — fixtures + instruments".

---

## ELI5

The project's only test of "is the 1D string a faithful hash of the 3D structure?" was to turn a
structure into a string, turn the string back into a structure, and check the string comes out
the same. That test can never catch the dangerous failure — two genuinely *different* molecules
being given the *same* string — so an earlier audit built a different instrument: take a
structure, make its **mirror image**, and see whether the encoder can still tell the two apart.
The problem is that a floppy molecule's mirror image is often just the same molecule after some
bonds have twisted, so the mirror test kept crying wolf; 70 structures were stuck in an "unknown"
pile because nobody could tell twisting apart from genuine handedness. This lane built the tool
that separates them (ask: *can you reach the mirror by rotating bonds, without breaking any?*),
built two new test molecules chosen so that a future descriptor cannot pass by accident, and
built two new ways to make a "twin" molecule that a plain mirror image cannot produce. Jargon:
*isomer* = same atoms, different arrangement; *conformer* = same arrangement, different twist;
*enantiomer* = mirror-image pair; *diastereomer* = a different arrangement that is not a mirror
image; *automorphism* = a relabelling of atoms that leaves the molecule looking identical.

---

## The work, visually

```
LANE 7 = 1 prerequisite + 3 tasks.  Nothing in src/ was touched.

┌─ PREREQUISITE (shipped FIRST, commit cc0dd117, because Lane 5 was blocked on it) ────────┐
│                                                                                          │
│  THE FIXTURE-RATIONALE DECISION TREE                                                     │
│                                                                                          │
│  "Validate a metal Δ/Λ (helicity) descriptor."                                           │
│            │                                                                             │
│            ├─ on fac-Ir(ppy)3 alone?                                                     │
│            │        │                                                                    │
│            │        └─ its 3 chelates are UNSYMMETRIC (C,N donors)                       │
│            │              │                                                              │
│            │              └─► a descriptor that secretly encodes fac/mer                 │
│            │                   INSTEAD of helicity STILL PASSES  ✗  the single-          │
│            │                   fixture trap (this is exactly how the Y2 wave             │
│            │                   shipped a reflection-INVARIANT axial token)               │
│            │                                                                             │
│            └─ add ZUMNEC (tris(catecholato)Mo, 37 atoms)                                 │
│                     │                                                                    │
│                     ├─ HOMOLEPTIC; the 2 O donors per chelate share ONE                  │
│                     │  CanonicalRankAtoms(breakTies=False) rank                          │
│                     ├─ ⇒ NO fac/mer distinction exists to lean on                        │
│                     ├─ ⇒ metal helicity is the SOLE stereogenic element                  │
│                     └─► a fac/mer-in-disguise descriptor FAILS here  ✓                   │
│                                                                                          │
│  "…and validate the @SP (square-planar) path."                                           │
│            │                                                                             │
│            └─ add JEGKOW (Rh(I) SP, donors N / P / C-carbonyl / I, 31 atoms)              │
│                     │                                                                    │
│                     ├─ mirror is NOT a distinct isomer (0.099 Å)  ← CORRECT, not weak    │
│                     │  a square-planar complex is PLANAR ⇒ its coordination plane IS     │
│                     │  a mirror plane ⇒ 4 different DONORS give DIASTEREOMERS             │
│                     └─► reflection is the WRONG operator for @SP; the right one is a      │
│                          DONOR SWAP  ──────────────────────────────┐                      │
└────────────────────────────────────────────────────────────────────┼──────────────────────┘
                                                                     │ (creates the need
                                                                     │  for Task C)
┌─ TASK A ─ torsion-aware configurational oracle ──────────────┐     │
│  tools/injectivity/torsion_oracle.py                         │     │
│                                                              │     │
│  OLD (oracle.py::geometric_chirality): RIGID fit             │     │
│      mirror ──[proper rotation × automorphisms]──► base      │     │
│      floppy achiral molecule ⇒ "chiral"  ✗                   │     │
│                                                              │     │
│  NEW: search the input's OWN TORSION ORBIT                   │     │
│      Θ = (θ1…θk) over rotatable, ACYCLIC, non-terminal bonds │     │
│      d_mirror = min over Θ, min over automorphisms,          │     │
│                 properRMSD( mirror(base), coords(Θ) )        │     │
│                                                              │     │
│      d_mirror ≤ 0.5 Å ⇒ CONFORMATIONAL (collapse is RIGHT)   │     │
│      d_mirror > 0.5 Å ⇒ CONFIGURATIONAL (collapse is a LOSS) │     │
│                                                              │     │
│      + PAIRED POSITIVE CONTROL, same molecule/optimiser/     │     │
│        budget, target reachable BY CONSTRUCTION:              │     │
│           ctrl ok & mirror ok      → conformational          │     │
│           ctrl ok & mirror not ok  → configurational         │     │
│           ctrl NOT ok              → inconclusive            │     │
│      ⇒ "no match in this budget" is EVIDENCE, not PROOF;     │     │
│        budget + threshold ride inside TorsionVerdict         │     │
└──────────────────────────────────────────────────────────────┘     │
                                                                     │
┌─ TASK C ─ the two twin operators that were never built ◄────────────┘
│  tools/injectivity/twin_operators.py
│
│      mirror_z (pre-existing)  : coordinate transform, trivially valid,
│                                 answers ONE question (whole-molecule enantiomerism)
│      swap_donor               : exchange two donors' coordination SITES  → DIASTEREOMERISM
│      invert_axial             : negate ONE biaryl dihedral                → ONE AXIS AT A TIME
│      invert_tetrahedral       : exchange two substituent BRANCHES         → one sp3 centre
│
│  these are STRUCTURAL EDITS, so they can produce nonsense:
│      edit ──► vdW clash count BEFORE/AFTER (generator3d/clash.py::vdw_clash_count)
│                       │
│                       ├─ clash_after − clash_before > 0  ⇒ probe_operator REFUSES to score it
│                       └─ ok ⇒ certify distinctness on the TORSION ORBIT, then encode both
└──────────────────────────────────────────────────────────────────────────────────────────┘

┌─ TASK B ─ the donor swap, done properly ──────────────────────────────────────────────────┐
│  tools/injectivity/positional_isomers.py                                                  │
│                                                                                           │
│  ✗ ORIGINAL (Y1): hand-write two OIN strings with donors in swapped slots, see if the      │
│    key folds them. It does. PROVES NOTHING — for a square-planar complex with two          │
│    identical ancillary ligands, swapping slots 0↔1 IS a reflection of an ACHIRAL           │
│    complex, so the two strings name the SAME molecule and folding them is CORRECT.         │
│                                                                                           │
│  ✓ REDONE, two independent geometry-driven lines:                                          │
│    Line 1  corpus pairs   : same constitution (canonical SMILES) + DIFFERENT trans-donor   │
│                             multiset ⇒ certified distinct WITHOUT any OIN string           │
│    Line 2  swap_donor     : real geometry, symmetry-equivalent donor pairs SKIPPED,        │
│                             trans swaps are the BUILT-IN NEGATIVE CONTROL                  │
│                                                                                           │
│  verdict vocabulary: REFUTED / CONFIRMED / UNDETERMINED (render() computes it)             │
└───────────────────────────────────────────────────────────────────────────────────────────┘

LEGEND
  ✓ / ✗          the choice taken / the choice rejected, with the reason attached
  ──►            "therefore"
  d_mirror       best RMSD from the input's torsion orbit onto its own mirror image
  d_control      same search, against a target that is reachable by construction
  Δ/Λ            metal-centred helicity (P1 blind spot)
  @OH / @SP      RDKit's octahedral / square-planar stereo descriptors
  P1 / P2 / P3   the three named blind spots: metal Δ/Λ, axial atropisomer, metal-bound 2° amine
  key            oin/compare.py's round-trip equivalence key (deliberately lossy)
```

---

## Initial assumptions and hypothesis

Going in, the lane held five beliefs. Four survived; the fifth was replaced by something
sharper.

1. **`fac-Ir(ppy)3` is not a sufficient fixture for metal Δ/Λ.** Its three chelates are
   unsymmetric (C,N donors), so a descriptor that encodes *fac/mer* rather than *helicity* would
   still appear to work on it. This is the same shape of trap the Y2 wave hit, where a
   reflection-**invariant** axial token passed every guard because the only fixture exercised the
   easy single-axis case. **Corollary carried forward as a rule: Lane 5 must NOT be validated on
   `fac-Ir(ppy)3` alone.**
2. **The Y3 ambiguous residue is unsplittable with a rigid instrument.** `oracle.py`'s
   `geometric_chirality` is a rigid superposition test, so on a flexible molecule the mirror
   reads as chiral because it is a different *conformer*, not a different isomer.
3. **The Y1 donor-swap probe was unsound**, not merely inconclusive: it hand-wrote OIN strings
   and never established that its two strings denoted different isomers.
4. **`@SP` metal stereo had no instrument at all**, because a mirror cannot produce a
   square-planar diastereomer and no other twin operator existed.
5. **The planned Task-A instrument was a conformer pool** — embed the molecule's own graph with
   `EmbedMultipleConfs` and ask whether any conformer superimposes on the mirror under a proper
   rotation. `INJECTIVITY_Y3_UNKNOWN_UNKNOWNS.md` §Status states exactly that. **This plan was
   discarded during implementation as unsound; see "Dead ends" below.**

Because the lane changed no encoder output, the acceptance bar was that the test-suite floor
had to hold **exactly**, not merely stay above a threshold — every new test is additive and no
existing count could move.

---

## What was actually found

### CONFIRMED — the two fixtures behave as their rationale requires

Both are real crystal geometries, selected by a **census of the 26,230-structure corpus** for
their archetype (single metal, ≤ 90 atoms).

| fixture | what | mirror per rigid oracle | metal descriptor | vdW clashes |
|---|---|---|---|---|
| `tests/fixtures/ZUMNEC.xyz` | tris(catecholato)Mo, 37 atoms | **distinct, 1.34 Å** | `@OH` permutation **11 → 9** (flips) | 0 |
| `tests/fixtures/JEGKOW.xyz` | Rh(I) square planar, donors N / P / C(carbonyl) / I, 31 atoms | **not distinct, 0.099 Å** | `@SP` permutation 2 (does **not** flip) | 0 |

- ZUMNEC is **the only homoleptic tris-bidentate in the corpus whose two donors per chelate are
  symmetry-equivalent** (`CanonicalRankAtoms(breakTies=False)` gives them a single rank).
- Its 1.34 Å mirror RMSD is **2.7× the 0.5 Å threshold** and **13–27× the achiral controls**
  (0.05–0.10 Å).
- It carries **no** axial and **no** bound-amine axis, so a Lane-5 pass on it is unambiguous.
- It **reproduces the P1 `key_blind` collapse today** — pinned by
  `test_metal_stereo_fixtures.py::TestZumnecTrisBidentate::test_currently_key_blind`.
- JEGKOW was reordered metal-first, **which leaves the emitted OIN byte-identical**.
- JEGKOW's mirror being *not* distinct is chemically correct, not a weak fixture, and
  `TestJegkowSquarePlanar::test_mirror_is_the_same_isomer` asserts the achirality deliberately.

### CONFIRMED — the torsion oracle is right on ten known answers, all ten

| structure | expected | `d_mirror` | `d_control` | verdict |
|---|---|---:|---:|---|
| `EDOQIZ` — unsubstituted biphenyl on linear Au | conformational | 0.076 | 0.075 | ✅ conformational |
| `WAVGOS` — bis-NHC on linear Au, propargyl arms | conformational | 0.048 | 0.025 | ✅ conformational |
| `PERPIO` — Ni alkynyl + phosphine | conformational | 0.292 | 0.115 | ✅ conformational |
| `PdCl2-R-BINAP` — hindered biaryl | configurational | 3.448 | 0.032 | ✅ configurational |
| `YESKOZ` — two hindered axes | configurational | 2.819 | 0.017 | ✅ configurational |
| `fac-Ir(ppy)3` — metal Δ/Λ | configurational | 2.673 | 0.000 | ✅ configurational (rigid) |
| `ZUMNEC` — metal Δ/Λ | configurational | 1.065 | 0.000 | ✅ configurational (rigid) |
| `POJJOP` — metal-bound amine | configurational | 1.707 | 0.009 | ✅ configurational |
| `CisPlatin` | achiral | 0.003 | — | ✅ rigid achiral |
| `JEGKOW` | achiral | 0.067 | — | ✅ rigid achiral |

The cross-validation named in the task brief — **EDOQIZ should be REACHABLE** (unsubstituted
biphenyl on linear Au, so not a real atropisomer) and **PdCl2-R-BINAP should NOT be reachable**
(hindered) — is rows 1 and 4 above, and both come out as required.

`WAVGOS`, `PERPIO` and `EDOQIZ` are the three cases Wave 3 triaged **by hand** and called
conformational; the tool reproduces all three, which is the closest thing to an external check
that exists here.

> ⚠ **Two different instruments produce two different numbers for the same molecule, and this is
> not a contradiction.** ZUMNEC reads **1.34 Å** in the fixture table (rigid oracle,
> `oracle.py::is_distinct_enantiomer`, H-explicit graph) and **1.065 Å** as `d_mirror` in the
> table above (torsion oracle, heavy atoms only, heavy-skeleton automorphisms). Likewise
> `fac-Ir(ppy)3` is quoted at 3.19 Å in `INJECTIVITY_Y1_OVERVIEW.md` and 2.673 Å here. Always
> name the instrument when quoting a mirror RMSD from this program.

### CONFIRMED — the operator line of Task B: the key separates every distinct swap it scored

| article | swaps | result |
|---|---|---|
| `PtMeNH3ClBr-Cis` (square planar, four different donors) | 6 total | the **2 trans swaps** correctly read `distinct=False` (a trans exchange is a 180° rotation of the whole complex — the built-in negative control); **all 4 cis swaps** are distinct isomers the encoder separates **at both the raw string and the key** |
| `FeH2(CO)4` (octahedral, added in `b315b929`) | 6 total | **all six** produce distinct isomers and the key separates **every one** — this is the cis/trans dihydride pair, i.e. the fac/mer half of the question |

### CONFIRMED — `invert_axial` settles the Y2 multi-axis question that Wave 2 could not

`YESKOZ` carries two symmetry-equivalent hindered axes of opposite sign; flipping **one** turns
the (+,−) diastereomer into (−,−), an edit no mirror can make.

| build | oracle | raw strings | key | verdict |
|---|---|---|---|---|
| default (`OIN_EMIT_AXIAL` off) | distinct | **identical** | equal | `encoder_blind` |
| `OIN_EMIT_AXIAL=1` | distinct | **differ** | equal | `key_blind` |

**The Y2 axial token works on a multi-axis molecule.** The Wave-2 cohort A/B could not establish
this because its multi-axis arm failed for a *generator* reason (torsions relaxed out of the
hindered window). This probe is generator-free, so that confound does not apply, and the residual
fold now lives in `compare.py`'s `_AXIAL_TOKEN_RE`, **not in the encoder**. → handed to Lane 4.

### CONFIRMED — Y3 context that reframes what a pass rate means

Carried in from `INJECTIVITY_Y3_UNKNOWN_UNKNOWNS.md`, because every Lane-7 conclusion is read
against it:

- Over **3917 molecules (2080 passing, 1837 failing)**, **77.8% of round-trip failures never
  test the notation**: generator timeout 1238 (**67.4%**), encoder refused the input 241
  (13.1%), generator produced nothing 191 (10.4%), canonicalization noise 123 (6.7%), output
  names a different isomer 44 (2.4%). The median failing molecule ran **300.3 s against a 300 s
  budget**.
- ⇒ **a round-trip pass rate is substantially a measure of generator throughput**, and compute
  buys apparent accuracy.
- The UU hunt over **299 structures** found **0 new confirmed blind spots** (a confirmed blind
  spot = the OIN collapses the twins *and* InChI separates them); **65 / 299 = 21.7%** of
  collapses trace to the already-named P1/P2/P3 axes; 130 (43.5%) are correct injectivity; 70
  (23.4%) are the ambiguous residue; 25 (8.4%) carry a CIP stereocentre; 9 (3.0%) are correct
  invariance.

### REFUTED — linkage isomerism as a blind spot

The key **distinguishes** both classic cases: thiocyanate N-bound vs S-bound (`N{0}=C=S` /
`S{0}=C=N`) and nitro N-bound vs nitrito O-bound. The vertex colouring carries donor-atom
identity, so linkage isomers do not collide. This Y1 target-map entry is closed.

### REFUTED — two beliefs the tool itself held while it was being built

See "Dead ends" — both made the tool *confidently wrong*, not merely imprecise.

### NOT RECORDED IN THE REPO — the Task-B corpus-pair line, and the residue split

⚠ Both of these are real gaps a future agent will hit.

- The Task-B **Line 1** tables and the final REFUTED/CONFIRMED/UNDETERMINED verdict are written
  to `results-injectivity-y3/positional_isomers.md`, which is **gitignored and absent** from this
  checkout. `docs/INJECTIVITY_Y3_RESIDUALS.md` §Part 2 points at that file and states **no
  verdict inline**. The only Line-1 fact that is durably recorded is the false positive the scan
  had to be hardened against (below).
- The **70-case residue split** (`tools/injectivity/config_split.py`) likewise writes to
  `results-injectivity-y3/`, and **no split counts appear in any tracked file**.
  `docs/V045_STATUS_2026-07-25.md` lists exactly this as outstanding: *"the 70-case UU residue
  split and its hand inspection, the Task-B corpus-pair line, and the exact suite count."*
- Consequently **`docs/INJECTIVITY_Y1_OVERVIEW.md` still records "donor swap undetermined"** in
  its Wave-3 roadmap line. The operator line (Line 2) is a clean *separates-everything* result
  as tabulated above; the corpus line was never published. Do not upgrade the Y1 overview's
  verdict without regenerating Line 1.

---

## What was done

### Prerequisite — the two fixtures (commit `cc0dd117`)

| artefact | detail |
|---|---|
| `tests/fixtures/ZUMNEC.xyz` | 39 lines; tris(catecholato)Mo, 37 atoms |
| `tests/fixtures/JEGKOW.xyz` | 33 lines; Rh(I) square planar, 31 atoms, **reordered metal-first** (emitted OIN byte-identical) |
| `tests/unit/test_metal_stereo_fixtures.py` | 166 lines, the guard module |
| `tools/injectivity/report.py` | both registered in `default_probe_set()` |

Guard classes and what each one is *for* — the asymmetry is deliberate and load-bearing:

- `TestFixtureGeometryIsClean::test_no_vdw_clashes` / `test_single_metal_center` — the oracle's
  verdict is meaningless on a bad geometry, so the geometry is gated first.
- `TestZumnecTrisBidentate::test_three_chelates_with_symmetry_equivalent_donors` — asserts three
  bidentate chelates, each with its two donors sharing **one** symmetry rank. This is the
  assertion that makes ZUMNEC a *guard* rather than just another fixture.
- `::test_mirror_is_a_distinct_isomer` — `rmsd > 1.0`.
- `::test_metal_descriptor_is_octahedral_and_flips` — shape `["OH"]`, and the permutation must
  differ between base and mirror.
- `::test_no_other_axis_is_implicated` — `rep.base.axial == []` and `rep.base.bound_amine == []`.
- `::test_currently_key_blind` — assert-current-behaviour: `VERDICT_KEY_BLIND`.
- `TestZumnecAspirational::test_metal_chirality_should_diverge_at_key` — `@unittest.expectedFailure`.
  This is the **twin** of the existing `fac-Ir(ppy)3` P1 aspirational test, so **Lane 5 must flip
  both**; fixing only the easy fixture leaves this one red.
- `TestJegkowSquarePlanar::test_mirror_is_the_same_isomer` — asserts `metal_flips` is **False**,
  the mirror is **not** oracle-distinct, `raw_equal` is True, and the verdict is
  `VERDICT_INVARIANT_OK`. Asserting the *opposite* here would be a category error.

Suite floor moved **605 OK / 3 skip / 3 xfail → 615 OK / 3 skip / 4 xfail** (+10 OK, +1 xfail).

**Rejected alternative:** validating Lane 5 on `fac-Ir(ppy)3` plus a synthetic/idealised
tris-chelate. Rejected because a hand-built geometry cannot certify "the only homoleptic
tris-bidentate in the corpus with symmetry-equivalent donors" — that claim requires the census,
and the census is what makes ZUMNEC a *guard against a fac/mer-in-disguise descriptor* rather
than merely a second example.

### Task A — the torsion-aware configurational oracle (commit `be6a6ae4`)

New: `tools/injectivity/torsion_oracle.py` (574 lines), `tools/injectivity/config_split.py`
(237 lines, the residue-split driver), `tests/unit/test_torsion_oracle.py` (226 lines); plus a
5-line change to `tools/injectivity/uu_hunt.py` so it persists the **FULL** residue rather than
the first 25, which is what makes all 70 splittable.

Public API: `configurational_verdict(xyz_path, *, charge, threshold, restarts, sweeps, grid,
seed, max_autos, mol, coords, name) -> TorsionVerdict`.

Defaults, all module constants so they can be quoted: `MATCH_THRESHOLD = 0.5` (Å — deliberately
**identical** to the rigid oracle's, so the two are directly comparable), `AUTO_DESCENT_K = 8`,
`N_RESTARTS = 12`, `N_SWEEPS = 4`, `COARSE_GRID = 12` (30° steps), `REFINE_STEPS = (10.0, 3.0)`,
`seed = 42`.

Key internals:

- `rotatable_torsions(mol, root)` — a bond qualifies when it is **acyclic**, both ends carry
  another **heavy** neighbour, and removing it disconnects the graph. Metal–ligand bonds *are*
  included: rotating a monodentate ligand about its donor bond is genuine conformational freedom.
- `_side_without(mol, n, *, cut, start)` — the fix for bug 1 below. **The cut removes the BOND,
  never either ATOM.**
- `dihedral_negating_theta(mol, coords, tors)` — **the most important starting point, and it is
  not a heuristic.** A reflection preserves every bond length and angle and **negates every
  signed dihedral**, so the torsion vector that negates all rotatable dihedrals is the physically
  correct guess; a purely conformational chirality is found there immediately. Both signs
  (`neg`, `-neg`) are seeded, so the search does not depend on `_rotate`'s sign convention.
- `_automorphism_perms(mol, heavy, max_autos)` — enumerated on the **heavy skeleton**
  (`heavy_skeleton()`), then filtered to permutations preserving each atom's hydrogen count and
  formal charge, so a `-CH2-` can never be matched onto a `-CH3`.
- `_batch_proper_rmsd(target, moved, perms)` — batched Kabsch with the **reflection forbidden**
  (`diag(1, 1, det)`), matching `oracle._kabsch_proper_rmsd` term for term.
- `_Scorer` — full automorphism set scores the *final* vector; only the best `AUTO_DESCENT_K = 8`
  are carried through the inner descent. More automorphisms can only *lower* an RMSD, so a subset
  biases the tool toward `configurational` — the direction that costs hand inspection rather than
  the direction that hides a blind spot.
- `_Scorer.residual_profile` — per-atom leftover deviation at the best vector, worst 6 first.
  This exists for the ring-pucker false positive (below) and is attached to every
  `configurational` verdict.
- `TorsionVerdict` carries `budget`, `threshold`, `d_control`, `n_autos_total`,
  `n_autos_descent`, `rigid_mirror_rmsd` and the note, so **a conclusion cannot be quoted without
  the search that produced it**. `budget = (restarts + 3) * max(1, sweeps) * max(1, len(tors)) *
  (grid + 4)`; the `+ 3` was added in `2f675157` because the three non-random starts (the
  undeformed structure and both signs of the negating seed) were being omitted, under-reporting
  the search by `3/(N+3)`.

Two documented modelling choices: **heavy atoms only** (terminal rotors move hydrogens without
changing configuration and would only add noise dimensions), and the automorphism-subset descent
described above.

Guards in `tests/unit/test_torsion_oracle.py`, in three layers:

1. the torsion model — `test_rotation_preserves_every_bond_length`,
   `test_dihedral_negating_seed_negates_dihedrals`,
   `test_batched_kabsch_matches_the_scalar_one`, `test_reflection_is_forbidden`;
2. the two bugs, each a regression test —
   `TestChelateBondsAreNotRotatable::test_chelate_locked_metal_bonds_are_excluded` (asserts
   `tors == []` for both `fac-Ir(ppy)3` and `ZUMNEC`),
   `::test_a_monodentate_metal_bond_is_still_rotatable` (the exclusion must not over-fire),
   `TestAutomorphismCompleteness::test_heavy_skeleton_finds_more_than_the_starved_path`,
   `::test_every_permutation_is_a_real_automorphism`;
3. verdicts on known answers — `TestVerdictsOnKnownAnswers::test_achiral_control_is_rigid_achiral`,
   `::test_square_planar_four_donors_is_rigid_achiral`,
   `::test_metal_delta_lambda_is_configurational` (also asserts `n_torsions == 0`),
   `::test_hindered_biaryl_is_configurational_with_a_passing_control` (the strong form:
   `d_mirror > 1.0` **and** `d_control < threshold`), `::test_verdict_carries_its_own_search_size`.

Suite delta: **+13 tests, all OK, ~16 s.**

### Task B — the donor swap, from real geometry (commits `df125292`, `b315b929`, `d42194ef`)

New: `tools/injectivity/positional_isomers.py` (428 lines at introduction; +58/−26 in
`d42194ef`).

- `TRANS_ANGLE_DEG = 150.0`; `TRANS_ANGLE_LADDER = (140.0, 150.0, 160.0)`; `MAX_ATOMS = 70`.
- `arrangement_signature(mol, coords, cutoff)` — the multiset of **trans donor pairs** as sorted
  element pairs. For a fixed constitution this is a **configurational invariant**: unchanged by
  any rotation of the complex or of any bond, and it differs between cis/trans and between
  fac/mer. **No OIN string enters the certification.**
- `constitution_key(mol)` — canonical SMILES of the perceived complex, stereo stripped.
- `_robustly_different(la, lb)` — three conditions, each removing one way a crystal geometry can
  fake an isomer difference: differ at **every** cutoff; trans-pair **count** matches at every
  cutoff; the count is stable across the ladder within each structure.
- `scan(dataset, limit, checkpoint)` — pass 1 is rdkit-free formula bucketing; pass 2 perceives
  only the shortlist and is **checkpointed to a JSONL sidecar and resumable**.
- `operator_pairs(paths)` — the Line-2 driver, over `PtMeNH3ClBr-Cis`, `CisPlatin`, `JEGKOW`,
  `TransPlatin`, and `FeH2(CO)4`.
- `render(...)` — computes the verdict: `CONFIRMED` if anything distinct is folded, `REFUTED` if
  distinct pairs were scored and none folded, `UNDETERMINED` if nothing was scorable.

**Why the resumability change was needed and is worth keeping** (`d42194ef`): pass 2 is hours of
perception, and the whole result used to be written only at the very end. **Two runs were lost to
a harness task timeout — the positional scan at 2200/3116 and `uu_hunt` at 150/300** — roughly
1.5 hours of a loaded machine for nothing. Failures are recorded too, otherwise a resumed run
retries every unperceivable structure; a torn final line from a hard kill is skipped rather than
fatal. Verified: a second `--limit 300` run resumes all 11 perceptions instead of redoing them.
Long sweeps in this lane are launched with
`systemd-run --user -p OOMPolicy=continue -p MemoryMax=…`, per the project sweep protocol, so
they are detached from the harness process tree.

**Rejected alternative:** keeping the hand-written string pair and simply adding more of them.
Rejected because the defect is not sample size — it is that *nothing independent certified the
two strings denoted different isomers*. Adding strings multiplies an unsound comparison.

### Task C — the twin operators (commit `df125292`)

New: `tools/injectivity/twin_operators.py` (566 lines), `tests/unit/test_twin_operators.py`
(240 lines); `tools/injectivity/report.py` gains `--operators` (+108 lines), wiring all three
into the curated probe table.

> ⚠ **Naming discrepancy, flagged.** The task brief and the commit subject both call the second
> operator `invert_stereocenter`. There is **no function of that name.** It ships in two
> flavours — `invert_axial(mol, coords, axis_index=0)` and
> `invert_tetrahedral(mol, coords, center)` — and `invert_stereocenter` is only a section header
> in the module. Grep for the two real names.

- `EditedTwin` — carries `operator`, `detail`, `coords`, `clash_before`, `clash_after`, `error`,
  `freeze_bonds`, and a `geometry_ok` property. `CLASH_TOLERANCE = 0`: a twin is rejected when the
  edit introduces vdW clashes the original did not have.
- `swap_donor(mol, coords, donor_a, donor_b)` — each ligand is rotated about the metal by the
  minimal rotation carrying its own donor direction onto the other's, with the metal–donor
  distance preserved. Refuses a **chelate** donor (not detachable by cutting its own metal bond)
  and refuses two donors on the same ligand.
- `enumerate_donor_swaps(mol)` — **skips symmetry-equivalent donor pairs**, because exchanging
  those is the identity on the isomer. This is the trap the hand-written probe fell into.
- `invert_axial(...)` — negates the signed dihedral of **one** hindered axis via
  `detect_axial_axes`, and records `freeze_bonds = ((ax.a1, ax.a2),)`.
- `invert_tetrahedral(...)` — requires **all substituents symmetry-distinct** (else "not a
  stereocentre"), requires ≥ 2 detachable branches, generates **every** eligible exchange, ranks
  by introduced clash, and relaxes the winner in torsion space (`relax_torsions`) before the gate
  rules.
- `_exchange_branches(...)` — rotate **then translate**, an isometry.
- `torsion_orbit_distance(...)` — distinctness against an arbitrary target (what a *diastereomer*
  comparison needs), with the same paired positive control; `freeze_bonds` removes named bonds
  from the search.
- `probe_operator(...)` — if `not twin.geometry_ok`, returns with `oracle_distinct=None` and
  `key_equal=None`. **A rejected twin is never scored.** Otherwise it certifies distinctness,
  encodes both through `XYZToSMILES().convert()`, and classifies via `twin_collision.classify`.

Guards in `tests/unit/test_twin_operators.py`, in three groups mirroring the three risks:

- rigidity — `TestSwapDonorIsRigid::test_intra_ligand_bond_lengths_are_untouched`,
  `::test_chelate_donors_are_refused`, `::test_symmetry_equivalent_donors_are_not_enumerated`;
- filtering — `TestClashGateFilters::test_some_swap_is_rejected_by_the_clash_gate` (the gate must
  bite **and** must not reject everything), `::test_a_rejected_twin_is_never_scored`;
- honest distinctness in **both** directions —
  `TestDistinctnessIsHonestBothWays::test_trans_swap_is_the_same_isomer`,
  `::test_cis_swap_is_a_distinct_isomer_the_encoder_separates`;
- `TestInvertAxial::test_flips_exactly_one_axis`,
  `::test_the_flipped_axis_is_frozen_during_certification`,
  `::test_default_build_collapses_the_single_axis_flip` (skipped if `OIN_EMIT_AXIAL` is set);
- `TestInvertTetrahedral::test_non_stereocentre_is_refused`,
  `::test_ring_locked_centre_is_refused`, `::test_branch_internal_geometry_survives`,
  `::test_it_actually_inverts_the_centre` (the P3 oracle's signed tetrahedral volume must change
  sign).

Suite delta: **+14 tests, all OK, ~41 s.**

⚠ `TestInvertAxial.setUp` pins `OIN_CANONICAL_PERCEPTION=0`. Reason recorded in the test:
`OIN_CANONICAL_PERCEPTION` (default-ON since v0.4.5) reads more of the porphyrin as aromatic, so
one YESKOZ meso wall stops satisfying `_is_atropisomer_candidate`'s `not GetIsAromatic()` gate and
the hindered count measures **1 instead of 2**. With one axis there is nothing to isolate. Neither
count changes any emitted string — both YESKOZ axes are non-stereogenic and the token is empty
either way.

---

## Dead ends and refutations

### 1. The planned instrument — a free conformer pool — is UNSOUND here, and was discarded

⚠ **This is a direct discrepancy with the task brief.** The brief describes Task A as *"Added a
conformer pool of the molecule's own graph and asked whether ANY conformer superimposes on the
mirror under a PROPER rotation."* That is the **plan** as written in
`INJECTIVITY_Y3_UNKNOWN_UNKNOWNS.md` §Status. It is **not what shipped**, and the module docstring
plus commit `be6a6ae4` both explain why in the same words:

> `EmbedMultipleConfs` on the graph would generate **both handednesses** of exactly the axes under
> test — metal Δ/Λ and atropisomerism are **not carried in the graph** — so the pool would fit the
> mirror **every time** and the test would answer `conformational` **unconditionally**.

The delivered instrument is a **torsion-space search from the actual input structure**: rotating a
dihedral is a continuous deformation that cannot change any configuration, so the torsion orbit of
the input *is* the conformational orbit of its configuration. The "N conformers / threshold must be
stated" discipline from the brief survives intact, just expressed as
`budget = (restarts + 3) * sweeps * len(tors) * (grid + 4)` plus `threshold`, both carried inside
`TorsionVerdict`.

### 2. Cutting the ATOM instead of the BOND — a confidently wrong answer, not a near miss

Rotatability was tested by asking whether the graph disconnects when an **atom** is blocked. **A
metal is a cut vertex**, so blocking it detaches every ligand from every other and each
metal–donor bond of a **chelate** looked rotatable. Result: `fac-Ir(ppy)3`'s Δ/Λ mirror was
"reached" at **`d_mirror = 0.042 Å`** by swinging whole ligands off their own chelate rings.
Fixed in `_side_without`; regressed by
`test_torsion_oracle.py::TestChelateBondsAreNotRotatable`.

### 3. Automorphism starvation on the H-explicit graph — and it is STILL LIVE in `oracle.py`

Automorphisms were enumerated on the H-explicit graph, where methyl rotations consume the whole
`maxMatches` budget on permutations that leave every heavy atom fixed. On `EDOQIZ` (two
*tert*-butyls) **4000 full-graph matches collapse to 6 distinct heavy images, while the heavy
skeleton has 864.** Since more automorphisms can only *lower* an RMSD, starvation **inflates every
number** — which is why a freely rotating biphenyl read as configurational.

> ⚠ **`tools/injectivity/oracle.py` has the same starvation and it was deliberately LEFT IN
> PLACE.** It changes **no curated fixture's verdict** (`fac-Ir(ppy)3`, BINAP, POJJOP and
> CisPlatin all stay on the same side of the threshold), so nothing published in `BASELINE.md` §3
> moves. But it inflates the rigid mirror RMSD on methyl-rich species — **`EDOQIZ` reads 2.55 Å
> starved and 0.50 Å complete** — so **the rigid oracle over-reports chirality at dataset scale**,
> and `uu_hunt`'s "chiral" gate inherits that. Correcting it would change published per-fixture
> RMSDs, so it was flagged for the release owner rather than changed inside a measurement lane.
> The residue split is unaffected: `config_split` uses the complete enumeration and re-decides
> every case.
>
> **This latent defect then cost real time in Lane 8.** `oracle.py` reported `ROGYAO_comp_0` as
> *"distinct, mirror RMSD 2.586 Å, ENCODER-BLIND (total)"*. The molecule is **achiral** — the
> torsion oracle superimposes its mirror at **0.423 Å over 4 automorphisms** — and the false
> positive came from the 4000-automorphism cap being starved by four methyls. Two Lane-8
> assertions had the chemistry backwards because of it. See `LANE-08-stable-stereo-renumbering.md`.

### 4. Rotate-then-SCALE for a branch exchange — a similarity, not an isometry

Branches were moved by rotating then **scaling** by the distance ratio. Swapping an N–H branch
against an N–Pt branch that way **compressed every bond in the moved branch by half**. Fixed to
rotate-then-translate in `_exchange_branches`; regressed by
`TestInvertTetrahedral::test_branch_internal_geometry_survives`.

### 5. `invert_tetrahedral` accepted non-stereocentres — the hand-written mistake in a new costume

It "inverted" a metal-bound ammine by exchanging one of three **equivalent** hydrogens for the
metal, producing the **identical molecule** and reporting success. Now requires all substituents
symmetry-distinct (the same test `config_oracle.bound_amine_centers` applies). Regressed by
`TestInvertTetrahedral::test_non_stereocentre_is_refused`.

### 6. A free torsion orbit un-does an axial flip — every axial twin would score "not distinct"

Atropisomerism is conformational **plus a rotational barrier**, and a barrier is exactly what a
geometric torsion search does not model. So the distinctness search **must freeze the edited
bond**; otherwise the orbit rotates the flip straight back. Hence `EditedTwin.freeze_bonds` and
`torsion_orbit_distance(..., freeze_bonds=...)`. Regressed by
`::test_the_flipped_axis_is_frozen_during_certification`.

### 7. `BOCYEA_comp_0` vs `BOCYEA_comp_1` — the first Line-1 "pair" was a false positive

Two **crystallographically independent copies of the same compound**, separated only by one
donor–M–donor angle sitting either side of the 150° trans cutoff. Hardening applied, all three
parts recorded in `_robustly_different` and `TRANS_ANGLE_LADDER`: require the arrangement to
differ at **every** cutoff in 140/150/160°, require the trans-pair **count** to match (a differing
count means a straddled boundary or a different coordination number, neither of which is a
positional isomerism), and **exclude same-refcode pairs**.

### 8. `invert_tetrahedral`'s measured scope limit — recorded, not papered over

The inversion itself is correct: POJJOP's signed tetrahedral volume flips sign. But **every
stereocentre in the fixture set is locked** — inside a chelate ring (BDPP, DPDME) or bound to the
metal (POJJOP) — and for all of them the rigid exchange drives a substituent into the coordination
sphere. Torsion relaxation reduces the damage (**POJJOP 3 → 2 clashes, DPDME 10 → 4**) but does not
clear it. Reaching those would need **bond-angle relaxation, which a generator-free instrument must
not do.** For a locked centre the instrument remains `mirror_z` or `invert_axial`; this operator is
for stereocentres on a freely rotating pendant.

### 9. Two things a `conformational` / `configurational` verdict does NOT mean

- **`conformational` is geometric, not energetic.** An isolable atropisomer whose axis is a free
  (non-ring-locked) torsion is "reachable by rotating bonds" and reads `conformational`. That is
  the strict geometric truth and the chemically incomplete one. It does **not** corrupt the residue
  split, because the barrier question is triaged **upstream**: `uu_hunt` routes any structure with
  a stereogenic hindered axis into `known_axis` (P2) *before* the residue is formed, using
  `oin.axial.detect_axial_axes`' hindrance heuristic. The two instruments are complementary and
  **neither should be used for the other's question.**
- **`configurational` can over-report on a RING PUCKER** (recorded in `4f0a81fc`). The torsion
  orbit deliberately excludes ring bonds — rotating one would tear the ring open — but a saturated
  ring has *conformational* freedom the orbit therefore cannot express. The λ/δ pucker of an
  ethylenediamine-type chelate is the textbook case: its mirror is a different pucker of the same
  isomer, no torsion vector reaches it, and the tool reports `configurational` even though
  collapsing it is correct. **This is what the per-atom residual profile is for**: a residual
  concentrated on the saturated atoms of a chelate ring is the signature of a pucker false
  positive; one localised on a metal's donor set or on a biaryl is a genuine configurational
  difference. **Every `configurational` survivor is hand-checked against this before it is called a
  candidate blind spot.**

---

## Where it landed

**Branch `swimlane/v045-lane7`, tip `d42194ef`. Fully merged** — `main` is 164 commits ahead of
the tip and 0 behind, so `git log main..swimlane/v045-lane7` is empty. Merged via
`28d0870d` (→ `trial/v045-integration`) and `4d92d828` (→ `release/v0.4.5`).

Commits, in order:

| commit | subject |
|---|---|
| `cc0dd117` | `test(lane7): Δ/Λ tris-bidentate + 4-donor square-planar fixtures for Lane 5` |
| `be6a6ae4` | `tools(lane7): torsion-aware configurational oracle + residue splitter` |
| `df125292` | `tools(lane7): swap_donor + invert_stereocenter twin operators` |
| `dd6aa6dc` | `docs(lane7): Y3 residuals write-up -- fixtures, torsion oracle, donor swap, operators` |
| `b315b929` | `tools(lane7): add the octahedral donor-swap probe` |
| `4f0a81fc` | `docs(lane7): record the ring-pucker over-report the torsion oracle can make` |
| `2f675157` | `tools(lane7): count the seeded starts in the reported search budget` |
| `d42194ef` | `tools(lane7): make the corpus scan resumable` |

**Levers: none.** No env var was introduced and **no file under `src/` was touched**, so
levers-OFF bytes are unchanged by construction. That is why the suite floor had to hold exactly:
605 OK / 3 skip / 3 xfail → 615 OK / 3 skip / 4 xfail after the fixtures, then +13 and +14 from
Tasks A and C. ⚠ The **final combined suite count was never recorded** —
`docs/V045_STATUS_2026-07-25.md` lists "the exact suite count" as outstanding for this lane.

Files, final state:

| path | status |
|---|---|
| `tests/fixtures/ZUMNEC.xyz` | new |
| `tests/fixtures/JEGKOW.xyz` | new |
| `tests/unit/test_metal_stereo_fixtures.py` | new — 6 tests + 1 `@expectedFailure` |
| `tools/injectivity/torsion_oracle.py` | new |
| `tools/injectivity/config_split.py` | new — residue-split driver |
| `tools/injectivity/twin_operators.py` | new |
| `tools/injectivity/positional_isomers.py` | new |
| `tests/unit/test_torsion_oracle.py` | new — 13 tests |
| `tests/unit/test_twin_operators.py` | new — 14 tests |
| `tools/injectivity/report.py` | extended — `default_probe_set()` + `--operators` |
| `tools/injectivity/uu_hunt.py` | +5 lines — persist the FULL residue |
| `docs/INJECTIVITY_Y3_RESIDUALS.md` | new — 242 lines |

Reproduce (from the worktree root, main `.venv`; rdkit is pinned, do **not** `uv sync`):

```bash
PYTHONPATH=$PWD/src python -m tools.injectivity.report --probes --operators
PYTHONPATH=$PWD/src python -m tools.injectivity.torsion_oracle tests/fixtures/PdCl2-R-BINAP.xyz
PYTHONPATH=$PWD/src python -m tools.injectivity.uu_hunt --n 300      # writes the residue
PYTHONPATH=$PWD/src python -m tools.injectivity.config_split         # splits it
PYTHONPATH=$PWD/src python -m tools.injectivity.positional_isomers --scan
```

Output lands in `results-injectivity-y3/` — **gitignored, regenerable, seeded at 42**, and
**absent from this checkout**. Never write sweep output to `/tmp`.

Downstream state as of v0.4.6: Lane 5 shipped a chelate-aware Δ/Λ descriptor
(`oin/metal_config.py`) validated on **both** fixtures — it emits and inverts under reflection on
ZUMNEC and correctly emits nothing for JEGKOW — behind `OIN_EMIT_METAL_CONFIG`, which is
**held opt-in** in `oin/levers.py::_HELD_OFF`. So
`TestZumnecAspirational::test_metal_chirality_should_diverge_at_key` is still an
`@expectedFailure` on the default path.

---

## Open questions / for the next agent

1. **Regenerate the two missing result sets and publish their numbers.** `config_split` (the
   70-case residue split) and `positional_isomers --scan` (Task-B Line 1) both write only to the
   gitignored `results-injectivity-y3/`, and **neither's counts exist in any tracked file**. Until
   Line 1 is regenerated, `INJECTIVITY_Y1_OVERVIEW.md`'s "donor swap undetermined" is the
   authoritative verdict and should not be upgraded on the strength of the operator line alone.
2. **Every `configurational` survivor still needs hand inspection**, specifically against the
   residual profile, for the ring-pucker false positive. `configurational` is a *candidate* blind
   spot, never a finding.
3. **`oracle.py`'s automorphism starvation is a live, known, unfixed defect** with a release-owner
   decision attached (fixing it moves published per-fixture RMSDs). It has already produced one
   documented wrong answer that cost Lane 8 time (`ROGYAO_comp_0`). Decide it or re-flag it; do
   not use `oracle.py` for new work — use `torsion_oracle.py`.
4. **`invert_tetrahedral` cannot reach a locked stereocentre**, and closing that needs bond-angle
   relaxation, which a generator-free instrument must not do. If P3 work needs a single-centre
   twin on a locked centre, that is a design question, not a tuning question.
5. **The `OIN_EMIT_AXIAL` promotion evidence predates `OIN_CANONICAL_PERCEPTION` going default-ON**,
   and perception feeds `_is_atropisomer_candidate`. Measured on YESKOZ, hindered axes go **2 → 1**
   under canonical perception. `levers.py::_HELD_OFF` requires both Y2 cohorts to be re-measured
   with perception ON before promotion. Lane 7's `TestInvertAxial` sidesteps this by pinning
   perception OFF — that pin is a workaround, not a resolution.
6. **When the axial token is promoted, `compare.py`'s `_AXIAL_TOKEN_RE` fold must be removed in
   the same commit.** Lane 7 established that the residual collapse now lives in the key, not the
   encoder; leaving the fold in place makes the round trip structurally unable to verify what the
   token encodes.
7. **Lane-5 validation rule, restated because it is the whole point of the fixture work:** a metal
   Δ/Λ descriptor validated on `fac-Ir(ppy)3` alone has not been validated. It must also flip
   `TestZumnecAspirational`.
