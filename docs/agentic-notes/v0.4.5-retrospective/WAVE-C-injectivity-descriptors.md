# Wave C — the injectivity descriptors (Lanes 5 and 6): metal Δ/Λ and the metal-locked amine

**Purpose of the wave:** close the two remaining Y1 injectivity blind spots — **P1**, metal-centred
Δ/Λ helicity, and **P3**, the stereogenic secondary amine that is stereogenic *only because a metal
occupies its fourth position* — by **emitting** a descriptor the encoder currently discards. Both
lanes were scheduled after [Wave B](WAVE-B-canonical-slots.md) because both were believed to need
Lane 2's exported `canonical_slot_permutation()` to give the descriptor a reproducible frame.

Neither lane went as scheduled. **Lane 6 landed in v0.4.5 and consumed nothing from Lane 2. Lane 5
was never started in v0.4.5** and was built in v0.4.6, where it took **four formulations, three of
which were refuted by measurement.**

---

## ELI5

A round trip can only prove the encoder is *not lossy* if two genuinely different molecules produce
two different strings. Two cases where they did not: a molecule shaped like a three-bladed propeller
around the metal comes in left-handed and right-handed forms (Δ and Λ), and the string could not tell
them apart; and a nitrogen atom bonded to a metal is a left/right handed centre only because the metal
plugs its fourth slot — strip the metal, as the encoder does when it writes the ligand out, and the
chemistry software throws the handedness away as a nitrogen that can flip freely. This wave built two
extra descriptors to recover each one from the original coordinates. Both ship switched off, for the
same reason: a descriptor that *adds* information to the string means the 3D generator must be able to
rebuild what the string now claims, and until it can, switching the descriptor on turns a quiet wrong
answer into a loud failure. The most instructive result is that handedness around a metal is **not** a
property of where the donor atoms sit — it is a property of **which donors are paired into which
ligand** — and three plausible descriptors were built and thrown away before that became clear.

## The wave, visually

```
   ══ WAVE B: LANE 2 canonical_slots.py ══
      exports canonical_slot_permutation() / canonical_slot_map()
      ‼ and MEASURES: 0 of 150 molecules emit a metal `@` tag at all
                    │
        ┌───────────┴────────────────────────────┐
        │ believed dependency                    │ believed dependency
        ▼                                        ▼
  ┌───────────────────────────┐         ┌──────────────────────────────┐
  │ LANE 5 — P1 metal Δ/Λ     │         │ LANE 6 — P3 locked amine     │
  │                           │         │  OIN_EMIT_LOCKED_DONOR       │
  │  ✖ NOT STARTED in v0.4.5  │         │  ✔ LANDED in v0.4.5          │
  │    the 0/150 measurement  │         │  ◆ Lane 2 dependency was     │
  │    rescoped it from       │         │    MEASURED NOT TO EXIST     │
  │    "un-fold a collapse"   │         │    (nothing from             │
  │    to "CREATE a           │         │     canonical_slots consumed)│
  │    descriptor"            │         │                              │
  └────────────┬──────────────┘         └──────────────┬───────────────┘
               │ deferred to v0.4.6                    │
               ▼                                       ▼
  ╔═══════════════════════════════╗       ╔════════════════════════════════════╗
  ║ v0.4.6 — FOUR FORMULATIONS    ║       ║ SHIPPED BUT UNUSABLE IN THE        ║
  ║                               ║       ║ DEFAULT CONFIGURATION              ║
  ║ 1 signed volume in canonical  ║       ║                                    ║
  ║   slot order                  ║       ║ OIN_CANONICAL_BODY (default-ON)    ║
  ║   ✖ REFUTED — homoleptic ties ║       ║   reparses the ligand body …       ║
  ║     resolve by an ODD perm    ║       ║   … sanitizing a METAL-FREE        ║
  ║     ⇒ 1 → −1 on 2 of 4        ║       ║      fragment CLEARS the [N@]      ║
  ║       shuffles                ║       ║   ⇒ the tag is stamped, then       ║
  ║                               ║       ║     DISCARDED                      ║
  ║ 2 magnitude threshold on a    ║       ║   ⇒ test_locked_donor.py runs with ║
  ║   permutation-invariant       ║       ║     OIN_CANONICAL_BODY=0           ║
  ║   chirality index             ║       ║                                    ║
  ║   ✖ REFUTED — achiral JEGKOW  ║       ║ v0.4.6: THE OBVIOUS FIX IS         ║
  ║     −3.287e−04 vs chiral      ║       ║ MEASURABLY WRONG                   ║
  ║     ZUMNEC −4.807e−04         ║       ║   copy the tag onto the reparsed   ║
  ║     (1.5×) ⇒ NO threshold     ║       ║   donor → P3 emits, POJJOP passes  ║
  ║                               ║       ║   … but a tag set AFTER sanitize   ║
  ║ 3 unconstrained mirror-       ║       ║   MOVES the canonical WRITE ORDER, ║
  ║   superposition               ║       ║   and @/@@ is a parity relative    ║
  ║   ✖ REFUTED — calls chiral    ║       ║   to that order                    ║
  ║     ZUMNEC ACHIRAL: the perm  ║       ║   ⇒ RIFGUJ_comp_2's three ring     ║
  ║     search RE-PAIRS donors    ║       ║     CARBON tags flip between a     ║
  ║     into different ligands    ║       ║     structure and its MIRROR       ║
  ║                               ║       ║   ⇒ geometry says they must NOT:   ║
  ║ 4 CONSTRAIN perms to PRESERVE ║       ║     pseudo-asymmetric, lowercase   ║
  ║   CHELATE MEMBERSHIP          ║       ║     `s`, RELATIVE all-cis          ║
  ║   ✔ WORKS                     ║       ║   ⇒ REVERTED + permanent guard     ║
  ║     ZUMNEC 1.3752 Å           ║       ║                                    ║
  ║     JEGKOW 0.0582 Å           ║       ║ single-centre POJJOP COULD NOT     ║
  ║     = 24× margin              ║       ║ catch it; the MULTI-centre RIFGUJ  ║
  ║                               ║       ║ guards did                         ║
  ║ LESSON: Δ/Λ is a property of  ║       ╚════════════════════════════════════╝
  ║ CHELATE CONNECTIVITY, not of  ║
  ║ donor POSITIONS               ║        Fixtures both lanes are validated on
  ╚═══════════════╤═══════════════╝        came from LANE 7 (Wave A):
                  │                          ZUMNEC.xyz — helicity is its SOLE
                  ▼                            stereogenic element
   descriptor → emit (OIN_EMIT_METAL_CONFIG,   JEGKOW.xyz — 4 different donors,
                      default OFF)               square planar, correctly achiral
              → key FOLD (_METAL_CONFIG_TOKEN_RE,
                          un-fold obligation recorded)
              → GENERATOR reproduction in accept_fn
                 ◄── MUST exist SEPARATELY: the key FOLDS the token, so a key
                     match says NOTHING about helicity, and accepting on it
                     alone returns the WRONG ENANTIOMER while reporting success

Legend  ✔ works / landed    ✖ refuted or not started    ◆ dependency refuted
        ‼ the measurement that rescoped a lane
        Both wave levers ship default-OFF. Neither was promoted.
```

## Initial assumptions and hypothesis

1. **Both lanes depend on Lane 2.** A chiral descriptor needs a reproducible ordering of the metal's
   donors; Lane 2's canonical slot permutation was to be that ordering. Wave C was therefore
   scheduled last.
2. **Both lanes are "un-folding a collapse."** The Y1 audit had recorded P1 as `key_blind` —
   `fac-Ir(ppy)₃` and its enantiomer produce different *raw* strings but only by non-reproducible slot
   renumbering, which the round-trip key deliberately folds — and P3 as `encoder_blind`. Both were
   confirmed **recoverable from the input 3D** by the Y2 feasibility wave: the encoder discards a
   signal it has, it does not lack one.
3. **P3's fix is a carve-out at `core/chirality.py:722-727`,** narrowing the Zone-A predicate that
   drops the chiral tag of any P/N with `total_degree < 4`.
4. **The phosphorus case is the harder half of P3** and probably has a wider blast radius than the
   amine case.
5. **P1's descriptor is a signed volume** over the donors, taken in canonical slot order.
6. **Both levers stay opt-in for v0.4.5.** This was a confirmed product call, not a discovery:
   a lever that **adds information** to the string obliges the generator to reproduce what it emits,
   so promoting one converts a *silent false positive* into a *loud false negative*. That is the right
   direction of error but it is a separate decision from determinism, and determinism came first.

## What was actually found

### Refuted — the wave's own scheduling premise

**Lane 2 was not a dependency for Lane 6.** Nothing from `canonical_slots` is consumed. The
reproducible ordering came from a completely different insight, and it also refutes
[Lane 8](LANE-08-stable-stereo-renumbering.md)'s approach *for this case*: a chiral tag is a
**parity** for an atom with 3–4 real bonds (so bond-order dependence cancels — reorder the bonds and
both the parity and the writer's reading of it flip together) but an **absolute label** for an atom
with 2 real bonds + implicit H (so it does not). Measured directly, by holding one molecule fixed and
handing the SMILES writer two different bond orders at the stereocentre:

| fragment-atom class | does the emitted `@`/`@@` depend on bond order? |
|---|---|
| 3 or 4 real bonds | **yes** — the tag is a parity relative to bond order |
| 2 real bonds + implicit H | **no** — the same tag emits the same character either way |

A metal-locked donor lands in **both** classes: a PR₃ phosphine keeps three heavy neighbours once the
metal bond is cut, a secondary amine or phosphine keeps only two. So there are two frames. The
2-bond frame **drops the metal** — the neighbour the fragment is missing, standing in for the lone
pair, i.e. the one RDKit's own rule would place opposite the other three — and orders the remaining
three by `CanonicalRankAtoms(breakTies=False)` **symmetry-class rank**, which the stereogenicity gate
has already guaranteed are distinct, so **no index-dependent tie-break is ever consulted**. Measured
on POJJOP: signed volume **+9.405 across 6 random renumberings, −9.405 for the z-mirror**.

**Lane 5 was not started at all in v0.4.5.** The decisive input was [Lane 2](LANE-02-canonical-slots.md)'s
measurement that **0 of 150 molecules emit a metal `@` tag**. That converts the lane from "un-fold a
collapse the encoder already half-records" into "**create** a descriptor from nothing" — a larger
piece of work with a different risk profile than the lane budget assumed, and one that needs the Δ/Λ
tris-bidentate fixture [Lane 7](LANE-07-research-residuals.md) had only just built. Recorded as
deliberate deferral in `docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md` "Known gaps" §4 and in the release commit
message, not as an omission.

### Refuted — three of Lane 5's four formulations (v0.4.6)

| # | formulation | measured refutation |
|---:|---|---|
| 1 | **signed volume in canonical slot order** | Fails on **homoleptic** complexes. The canonical ordering *ties*, and some tie resolutions differ from others by an **odd** permutation — which flips the sign of a signed volume. Measured: the descriptor goes **1 → −1 on 2 of 4** donor shuffles, under pure renumbering. |
| 2 | **magnitude threshold on a permutation-invariant chirality index** | The index (commit `c72fbc1b`) genuinely solved the ordering blocker — 6 random donor permutations give an identical value where the ordered descriptor flipped. But the claim that "achirality falls out of the index rather than needing a planarity test" rested on two *synthetic* controls returning exactly `+0.000e+00`. On real structures: chiral **ZUMNEC −4.807e−04** vs achiral **JEGKOW −3.287e−04** — the **same order of magnitude, 1.5× apart**. Crystallographic pucker in an achiral complex produces a chirality index comparable to genuine helicity, so **no threshold separates them** and any `_CHIRALITY_EPS` is arbitrary. The exact-zero cancellation is a property of **idealized** coordinates only. |
| 3 | **unconstrained mirror-superposition** (mirror the donor set, ask whether any permutation + proper rotation superimposes it) | Right move for the *decision*, wrong answer: it calls chiral **ZUMNEC achiral**. As a bare point set ZUMNEC *is* achiral — six oxygens at octahedral vertices admit improper operations. A permutation search over **unlabelled** points is free to **re-pair the donors into different ligands**, which no physical operation can do, so it always finds a "symmetry" that is not chemically available. |
| 4 | **constrain the permutation search to preserve CHELATE MEMBERSHIP** (colour each donor by the ligand it belongs to) | **Works.** Best mirror-superposition RMSD: **ZUMNEC 1.3752 Å** vs **JEGKOW 0.0582 Å** — a **24× margin**. A mirror that has to re-pair chelates is correctly rejected. |

**The lesson, stated so nobody re-derives it: Δ/Λ helicity is a property of CHELATE CONNECTIVITY, not
of donor positions.** Reflecting a Δ complex yields Λ *only because* the reflection cannot be undone
while keeping the chelate pairing intact. This also retro-explains formulation 2: the chirality
index's non-zero reading for ZUMNEC was **residual crystallographic distortion** — the same magnitude
as achiral JEGKOW's pucker. **It was never detecting helicity at all.** Three of its four proven
properties (rotation invariance, reflection inversion, permutation invariance) were real and remain
real; they were properties of a quantity that does not mean what the lane needs. It is the one
refutation of the release that invalidated the descriptor's **input** rather than its method.

### Confirmed — P3 works, and is unusable in the shipped configuration

`POJJOP` (square-planar Pd whose sole stereocentre is a Pd-bound secondary amine) and its mirror
encoded **byte-identically**:

```
base   [Pd_SPL].Cc1ccc([NH]{0}Cc2ccccn{1}2)cc1.[Cl]{2}.[Cl]{3}
mirror [Pd_SPL].Cc1ccc([NH]{0}Cc2ccccn{1}2)cc1.[Cl]{2}.[Cl]{3}   ← identical
```

Behind `OIN_EMIT_LOCKED_DONOR` they now encode as `[N@@H]{0}` and `[N@H]{0}`, with
`mirror == swap(base)`. And it is **not usable in the shipped default**: `OIN_CANONICAL_BODY`
(default-ON since v0.4.5) reparses the ligand body, and sanitizing a **metal-free** fragment clears
the `[N@]` on a 2-degree amine — RDKit sees a freely inverting amine, which is *the exact behaviour
the descriptor exists to work around*. With both levers on, the tag is stamped and then discarded, so
`tests/unit/test_locked_donor.py` runs with **`OIN_CANONICAL_BODY=0`**. Recorded in
`src/oinsmiles/oin/levers.py::_HELD_OFF["OIN_EMIT_LOCKED_DONOR"]`.

### Confirmed — the obvious fix for that is measurably wrong (tried in v0.4.6)

Four lines. The correspondence is already available: `canonical_body._reparse_once`'s **Guard 2**
proves `donors[k] ↔ new_donors[k]` (same element, same heavy degree). Copying the chiral tag onto the
reparsed donor **does** make P3 emit under `OIN_CANONICAL_BODY`, and **POJJOP passes.**

It is still wrong. Setting a chiral tag **after** the sanitize introduces a stereocentre the canonical
ranker did not account for, which moves the canonical **write order** — and `@`/`@@` is a **parity
relative to that order**, not an absolute label. On `RIFGUJ_comp_2` (three Cu-bound amines on one
cyclohexane) the three ring-**carbon** tags then flip between a structure and its mirror. The geometry
says they must not: `AssignStereochemistryFrom3D` + `rdCIPLabeler` label those carbons lowercase
**`s`** — pseudo-asymmetric, a *relative* (all-cis) descriptor — and they read `s` identically for the
structure **and** its reflection:

```
atom   6  base=s  mirror=s   same
atom   8  base=s  mirror=s   same
atom  10  base=s  mirror=s   same
```

**Single-centre POJJOP could not have caught this.** The multi-centre `RIFGUJ_comp_2` guards did.
This is the Y2 lesson intact: *a measurement that only exercises the easy case will confirm a wrong
belief.* Reverted, with the mechanism written into `canonical_body.py::_reparse_once`,
`levers.py::_HELD_OFF` and the test module, plus a **lever-independent** guard
`tests/unit/test_locked_donor.py::TestRifgujRingCarbonsArePseudoAsymmetric` so the four-line "obvious
fix" cannot be silently re-attempted. **v0.4.5's decision to defer this rather than rush it was
therefore correct.**

### Refuted — three things Lane 6 was told, and one population estimate

1. **Narrowing the Zone-A predicate would have achieved nothing.** Two lines before the Zone-A branch,
   `recover()` calls `Chem.AssignStereochemistry(cleanIt=True, force=True)`, which clears a trivalent
   nitrogen **regardless of what the Zone-A branch then does**. Exempting the atom from the Zone-A
   clear exempts it from nothing. Moreover the unconditional clear is **load-bearing** — it is what
   keeps a genuinely invertible amine and a symmetric phosphine from acquiring spurious sp3 handedness
   — so it is left unconditional on purpose. The fix is a **restore that runs after the whole loop**,
   writing only atoms left `CHI_UNSPECIFIED`, so blast radius is exactly "metal-locked donors that
   would otherwise carry no tag at all" — which is precisely the gap
   [Lane 8](LANE-08-stable-stereo-renumbering.md) could not reach, since `stable_stereo.py:112-118`
   only *corrects tags that already exist.*
2. **The trivalent-P gap does not exist as described.** Measured over 400 corpus molecules: **all 7
   eligible metal-locked P donors, across all 3 molecules carrying them, already had a chiral tag**
   from the existing Zone-A lone-pair path — the restore was a no-op for every one. Worse for the
   briefing, the two molecules cited as Lane 8's phosphorus residuals are **not stereocentres at
   all**: `FEQFIS_comp_0`'s metal-bound P is aromatic P(N)(O)(O)Au with only **3** distinct symmetry
   classes (the two oxygens are equivalent), and `CEBVIR_comp_0`'s N donors are pyridine-type with 3
   neighbours. Emitting for either would be the over-sensitivity failure, not a fix.
   **Do not cite this lane as having closed a phosphorus defect.**
3. **Lane 2 was not a dependency** — see above.
4. **Y1's 2.85% (171/6000) blast-radius figure is an over-count**, because its geometric pre-filter
   ("N within 2.6 Å of a transition metal with exactly 1 H and 2 C neighbours") never tested
   symmetry-distinctness. Re-measured on a 400-molecule seed-42 sample (395 encoded): **5.3% carry at
   least one eligible metal-locked donor** (4.6% with an eligible N, 0.8% with an eligible P; 27
   eligible N centres / 7 eligible P centres) and **4.6% of strings change** with the lever ON. The
   two figures count different populations — the Y1 number over-counts the *stereogenic* motif, the
   lane numbers measure the gated one — so quote them with their provenance. And `CULGOF`, `ATEFIP`,
   `MOKCEV` were on the Y1 "blind" list and are **not stereocentres at all** (Rh-bound morpholine with
   symmetry-equivalent ring α-carbons; Pd-bound diethylamine plus a primary `[NH2]`; Ir-bound
   macrocyclic triamine with equivalent ring carbons).

### Confirmed — why the generator check must exist *separately* from the key

`compare.py` **folds** both new tokens: `_AXIAL_TOKEN_RE` at `src/oinsmiles/oin/compare.py:104`, and
`_METAL_CONFIG_TOKEN_RE = re.compile(r"\s*\|mc:[+\-]\|")` at `:117`. So a **key match says nothing
about helicity**, and a generator that accepted a conformer on key equality alone would hand back the
**wrong enantiomer while reporting success**. Lane 5's generator step is therefore an explicit
reproduction check in `accept_fn`, not a key comparison. The same reasoning is recorded for P3: the
key folds the descriptor too (`POJJOP` and its mirror still share a key with the lever ON), so
promoting that lever cannot move `facmer_divergent` or any other harness count — and equally, **the
harness cannot confirm the descriptor. Only the raw string can.**

## What was done

### Lane 6 — metal-locked N/P donor chirality, Y1 blind spot P3 → [LANE-06-locked-donor-amine-P3.md](LANE-06-locked-donor-amine-P3.md)

**Lever:** `OIN_EMIT_LOCKED_DONOR`, default **OFF**. (The promotion-gate document
`docs/agentic-notes/v0.4.5/PROMOTION_GATE_v0.4.5.md` §3 lists this lever under its charter name
`OIN_EMIT_BOUND_AMINE`; the shipped name in `levers.py::_HELD_OFF` and `locked_donor.py::ENV_LEVER`
is `OIN_EMIT_LOCKED_DONOR`. Same lever.)

**Code:** `src/oinsmiles/oin/locked_donor.py` (new) · `src/oinsmiles/core/chirality.py` (the restore
after `recover()`'s main loop) · `src/oinsmiles/utils/perception_tmc.py` (stamps the property during the
fragment rebuild) · `src/oinsmiles/oin/inline.py`. Guards: `tests/unit/test_locked_donor.py`.
Design document: `docs/agentic-notes/v0.4.5/LOCKED_DONOR_v0.4.5.md`.

**The notation.** Measured on rdkit 2025.09.3, the standard `[N@]`/`[N@@]` tag *is* usable, provided
it is applied after the last sanitizing step: `MolToSmiles` of a manually-set tag emits `[N@@H]`, but
`AssignStereochemistry(cleanIt=True)` clears it, and so does a sanitizing `MolFromSmiles`. Hence the
mechanism stamps an atom **property** during the fragment rebuild and converts it to a tag as
`recover()`'s final action — properties survive sanitization, tags on nitrogen do not. One further
place had to be taught this: `oin/inline.py` re-parses each fragment SMILES with `sanitize=False` to
attach `{slot}` markers, and the following `MolToSmiles` re-runs `assignStereochemistry(cleanIt=True)`
whenever the mol carries no `_StereochemDone` — deleting the descriptor at the last possible moment.
Marking perception done there fixes it, **gated twice** (lever on *and* the fragment actually contains
`[N@`) so no other fragment's serialization changes.

**The stereogenicity gate** (`plan_locked_donors`), four conditions, all narrow on purpose:
(1) N or P, **not aromatic**; (2) bonded to **exactly one** transition metal; (3) **exactly four**
neighbours in the metal-present mol — the metal plus three, which is what "the metal locks the fourth
position" means; (4) those four in **four distinct symmetry classes**
(`CanonicalRankAtoms(breakTies=False)`). Condition 4 is the over-sensitivity guard and it does real
work: `CisPlatin`, `PtMeNH₃ClBr` (ammine: M,H,H,H) and `Cis-PtCl₂(en)` (primary amine: M,H,H,C) show
no eligibility and no string change; `JUCCUH` is byte-identical.

**Measured encoder behaviour, lever ON** (`mirror` is the z-reflection):

| molecule | lever OFF | lever ON |
|---|---|---|
| `POJJOP` (1 locked amine, sole stereocentre) | mirror byte-identical | `[N@@H]{0}` vs `[N@H]{0}`, `mirror == swap(base)` |
| `RIFGUJ_comp_2` (3 Cu-bound amines) | mirror byte-identical | all three invert, nothing else moves |
| `TEGFET`, `KANYUW`, `LARYEI` | already diverged via a second stereocentre | still diverge, now with explicit amine descriptors |
| `JUCCUH` | pendant `[N@@H]` already captured | **byte-identical** — no eligible donor |

Three-property test: invariant under input renumbering (POJJOP whole-string 3/3; RIFGUJ descriptors
4/4), invariant under proper rotation (whole-string 3/3 both), flips under reflection
(`mirror == swap(base)`; on RIFGUJ, `swap` of only the `[N@]` tags).

Two scoping notes, both forced by measurement rather than convenience: `mirror == swap(EVERY tag)` is
**false** on RIFGUJ, correctly — its 1,3,5-trisubstituted cyclohexane ring carbons encode a *relative*
(all-cis) arrangement that does not invert under reflection, verified with the lever **off** so it is
not something the lane introduced. And full-string renumbering stability cannot be asserted on RIFGUJ:
with the lever off it drifts in 3 of 3 renumberings (ring-carbon tags flip, one loses its tag, and the
`{2}`/`{3}` slots swap) — that is Lane 8's stereo class plus Lane 2's slot drift, and asserting it
here would import two other lanes' open defects into this lane's guard.

**Honest limitation, stated in the design document:** for nitrogen the descriptor is
*encoder-authoritative but not RDKit-round-trippable* — a sanitizing re-parse of the OIN drops it, so
the generator cannot rebuild it and the round trip reports a **loud false negative** rather than a
silent collision. That is the same trade the axial lever makes and the reason this lever is
default-OFF. An out-of-band token (`|amine:0+|`, following `|ax:±|`) would survive re-parse and was
**rejected**: it needs its own atom identity, sign convention and multi-token ordering — three fresh
chances to re-run the Y2 "sorted by sign, therefore reflection-invariant" mistake — and it would give
phosphorus a *second* representation on top of the `[P@]` the existing Zone-A lone-pair path already
emits.

**Self-contradictory acceptance criteria, resolved not papered over:** the `@expectedFailure`
`test_metal_bound_amine_should_diverge_in_raw_string` runs with the lever **OFF**, so flipping it
would directly contradict "levers-OFF byte-identical". Left in place with the flip condition
documented. ([Lane 2](LANE-02-canonical-slots.md) hit the same shape of contradiction over
`get_input_order_key`.)

### Lane 5 — metal Δ/Λ helicity, Y1 blind spot P1 → [LANE-05-metal-delta-lambda-P1.md](LANE-05-metal-delta-lambda-P1.md)

**Not started in v0.4.5. Built in v0.4.6.** Lever `OIN_EMIT_METAL_CONFIG`, default **OFF**, emitting
a trailing `|mc:+|` / `|mc:-|` sidecar.

**Code:** `src/oinsmiles/oin/metal_config.py` (new, ~429 lines) · `src/oinsmiles/utils/perception_tmc.py`
(the emit hook) · `src/oinsmiles/oin/compare.py` (`_METAL_CONFIG_TOKEN_RE`, the fold) ·
`src/oinsmiles/generation/metallogen_adapter.py` (the generator reproduction check).
Guards: `tests/unit/test_metal_config.py` (16 tests). Functions of record:
`chirality_index`, `metal_config_sign`, `metal_config_token`, `_kabsch_proper_rmsd`, `is_achiral`,
`metal_config_sign_symmetry`, `_admissible_permutations`, `is_achiral_chelate_aware`,
`metal_config_token_chelate`, `token_for_mol`, `parse_metal_config_token`.

**The four-stage pipeline, and why each stage exists:**

1. **descriptor** — the chelate-aware symmetry test above. Its **input** is perception's
   metal-incident bonds, not a distance cutoff; that was the fix for the donor-set half of the
   ordering blocker.
2. **emit** behind `OIN_EMIT_METAL_CONFIG`, default OFF. Verified: **ZUMNEC → `|mc:-|`, its mirror →
   `|mc:+|`** (it inverts, as a chirality descriptor must); **JEGKOW and CisPlatin → nothing**
   (achiral, no false positive).
3. **key fold** — `_METAL_CONFIG_TOKEN_RE` strips the token before comparison, **with the un-fold
   obligation recorded** in `compare.py` and in `levers.py::_HELD_OFF`: the fold must be removed in
   the same commit that promotes the lever, for the same reason `_AXIAL_TOKEN_RE` must (folding makes
   the round trip structurally unable to verify what the token encodes).
4. **generator reproduction in `accept_fn`** — the generator now *requires* an accepted conformer to
   reproduce the requested helicity. **This must be a separate check precisely because stage 3 folds
   the token**: a key match says nothing about helicity, so accepting on the key alone would return
   the wrong enantiomer while reporting success.

Note the consequence recorded in `_HELD_OFF`: because `compare.py`'s key does not yet *know* the
token, **lever-ON round trips will report mismatches until it does.**

**Lane 5's status by property**, as measured:

| property | state |
|---|---|
| permutation / relabelling invariance of the index | SOLVED (measured over 8 permutations) |
| invariance under proper rotation | PROVEN |
| inversion under reflection | PROVEN (exact negation) |
| donor SET determination | SOLVED — perception's metal-incident bonds, not a distance ratio |
| point-set achirality test | BUILT and correct as such |
| **detects Δ/Λ helicity** | **YES, with formulation 4** — chelate-membership-constrained mirror superposition, ZUMNEC 1.3752 Å vs JEGKOW 0.0582 Å |
| wired to emit | **YES, behind a default-OFF lever** |

## Dead ends and refutations

| tried / believed | what killed it |
|---|---|
| "Lane 5 and Lane 6 both need Lane 2's `canonical_slot_permutation()`" | Lane 6 consumes nothing from `canonical_slots`; the reproducible frame came from the parity-vs-absolute-label distinction instead |
| narrow the Zone-A `total_degree < 4` predicate at `core/chirality.py:722-727` | `recover()` calls `AssignStereochemistry(cleanIt=True, force=True)` **two lines earlier**, which clears trivalent N regardless. The carve-out had to be a restore *after* the whole loop |
| use `bound_amine_centers`' `CanonicalRankAtoms(breakTies=True)` ordering as the emission mechanism | measured **not** invariant under renumbering: 2–11 distinct rank vectors over 20 renumberings. Fine as an oracle on one fixed molecule, useless for emission |
| "use the fragment's own bond order and let the canonical writer normalise it" (Lane 8's argument, correct *there*) | **drifted under renumbering in 2 of 3 trials on POJJOP.** A metal-locked donor can be a 2-real-bond atom, where the tag is an **absolute label**, not a parity — so the bond-order dependence does not cancel |
| an out-of-band `\|amine:0+\|` token instead of `[N@]` | rejected: needs its own atom identity, sign convention and multi-token ordering, i.e. three fresh chances to repeat the Y2 sign-sort mistake, and it duplicates phosphorus' existing `[P@]` representation |
| "the trivalent-P case is the harder half of P3" | **all 7 eligible P donors already carried a tag**; the restore was a no-op for every one. The two cited residuals are not stereocentres |
| Y1's 2.85% P3 blast radius | over-count — the geometric pre-filter never tested symmetry-distinctness. Measured 5.3% eligible / 4.6% of strings change |
| **copy the metal-locked tag onto the reparsed donor to rescue P3 under `OIN_CANONICAL_BODY`** (v0.4.6, four lines, POJJOP passes) | **MEASURED WRONG.** Setting a tag after the sanitize moves the canonical **write order**, and `@`/`@@` is a parity relative to it: `RIFGUJ_comp_2`'s three ring-**carbon** tags flip between a structure and its mirror, where geometry says they must not (pseudo-asymmetric, lowercase `s`, base=s / mirror=s at atoms 6, 8, 10). Reverted; permanently guarded by `TestRifgujRingCarbonsArePseudoAsymmetric` |
| Lane 5 formulation 1 — signed volume in canonical slot order | homoleptic complexes tie, and some tie resolutions differ by an **odd** permutation ⇒ **1 → −1 on 2 of 4** shuffles under pure renumbering |
| Lane 5 formulation 2 — magnitude threshold on a permutation-invariant chirality index | achiral **JEGKOW −3.287e−04** vs chiral **ZUMNEC −4.807e−04**: same order of magnitude, 1.5× apart. **No threshold separates them.** The exact-zero cancellation holds only for **idealized** coordinates (ideal square: `+0.000e+00`) — two clean synthetic controls were read as evidence of a general property |
| Lane 5 formulation 3 — unconstrained mirror-superposition | calls chiral **ZUMNEC achiral**. A permutation search over unlabelled points may **re-pair donors into different ligands**, which no physical operation can do |
| "achirality falls out of the chirality index, so no planarity/symmetry test is needed" | refuted by the two real structures above. **`chirality_index` was measuring residual crystallographic distortion, not helicity** — three of its four proven properties were real properties of the wrong quantity |
| accept a generated conformer on `canonical_roundtrip_key` equality when a `\|mc:\|` token is present | the key **folds** the token, so a match says nothing about helicity — it would return the wrong enantiomer while reporting success. Hence a separate `accept_fn` reproduction check |

**One residual recorded rather than resolved:** `ABIFAV_comp_0` (Pd, two adamantanecarboxylate O
donors and two *N*-isopropylcyclohexylamine donors) emits **both** amine descriptors with **opposite**
signs — `{1}` as `[N@@H]`, `{3}` as `[N@H]` — and yet base and mirror stay byte-identical with the
lever ON. Its two amine ligands are constitutionally identical with opposite configurations, so
internal compensation would explain it: reflecting the molecule exchanges which ligand is which and
the descriptor pair is reflection-invariant. `LARYEI_comp_0` shows the same pattern at its two
ethylenediamine nitrogens, and there the string *does* differ because its backbone carbons invert.
**What is not established** is whether the fold in `ABIFAV` is purely that (R,S) compensation or
whether the reflection-sensitivity of *slot assignment* also contributes. Distinguishing them needs a
conformer-independent isomer oracle, which the lane did not build.

## Where it landed

### Lane 6 (v0.4.5)

- **Branch** `swimlane/v045-lane6`, tip **`63ba9090`**, 5 commits: `bac039ef` (recover metal-locked
  N/P donor chirality behind `OIN_EMIT_LOCKED_DONOR`) · `b6e2310a` (three-property + over-sensitivity
  guards) · `9bc8b025` (the P3 write-up, the key-fold guard, the xfail correction) · `d64a4154`
  (measured prevalence; the trivalent-P gap refuted) · `63ba9090` (the measured acceptance table).
- **Merged into `release/v0.4.5`** as **`df7417a9`**, then into local `main` in **`0d165845`**, tag
  **`v0.4.5`**. Not pushed.
- **Suite on the branch: 620 OK / 3 skip / 3 xfail** against a 605/3/3 baseline — delta **+15**,
  exactly `tests/unit/test_locked_donor.py`, 0 regressions. Ruff `check` + `format` clean on
  `src/ tools/ tests/`.
- **Lever `OIN_EMIT_LOCKED_DONOR` is default-OFF and stays so**, with the reason in
  `src/oinsmiles/oin/levers.py::_HELD_OFF`: the information-adding trade, **plus** the recorded
  incompatibility with default-ON `OIN_CANONICAL_BODY`, **plus** the warning that the obvious fix was
  tried in v0.4.6 and measured wrong.
- **New defect found in passing, unrelated:** `tests/fixtures/BENVOG_comp_0.xyz` **hangs
  `get_tmc_mol` at baseline** (>120 s wall, 5 s CPU — blocked, not compute-bound). It killed a
  whole-fixture survey. Consistent with `BENVOG` also appearing in the `encode_fail` lane's unresolved
  `resonance_timeout` cohort and in `valsearch`'s over-cap set — likely one root cause, and it is in
  the repo's own fixture set.

### Lane 5 (v0.4.6)

- **Branch** `swimlane/v046-hfaithful`. Commits, in order:
  `26f504e3` (descriptor — sound, with two measured blockers named) ·
  `c72fbc1b` (permutation-invariant chirality index — solves the homoleptic ordering blocker) ·
  `13f2b999` (docs: the magnitude threshold is REFUTED) ·
  `c172a57c` (the descriptor's INPUT is wrong — donor positions cannot express Δ/Λ) ·
  `848af5c6` (**Δ/Λ DETECTED** — chelate-aware symmetry test, both fixtures correct) ·
  `5b363057` (put the descriptor under the suite) ·
  `e8485ebf` (wire it to the emit path behind `OIN_EMIT_METAL_CONFIG`) ·
  `820b92f5` (teach the round-trip key to FOLD `|mc:|`, un-fold obligation recorded) ·
  `27089512` (**generator now REPRODUCES the requested Δ/Λ** — lane complete end to end).
- **Released** in **`d799de1f`** (`release(v0.4.6): boron cage promotion + Lane 5 metal Delta/Lambda
  (P1)`), local `main` only, not pushed.
- **Suite:** 838 tests OK / 3 skipped / 4 xfail on the branch before the boron promotion; the
  promotion left exactly one failure out of 840, which was the boron blast-radius leak (see
  [Wave D](WAVE-D-integrate-promote-release.md)).
- **Default path verified byte-identical on all 61 fixtures.**
- **Both fixtures come from [Lane 7](LANE-07-research-residuals.md)**: `tests/fixtures/ZUMNEC.xyz`
  and `tests/fixtures/JEGKOW.xyz`. Without ZUMNEC, formulations 1–3 would each have looked like they
  worked.

### Documentation of record

`docs/agentic-notes/v0.4.5/LOCKED_DONOR_v0.4.5.md` (Lane 6) · `docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md` (Lane 5's four
formulations, and the P3 restoration negative result) ·
`docs/KNOWN_LIMITATIONS.md` "Encoder injectivity blind spots (Y1 audit)" ·
`docs/agentic-notes/injectivity/INJECTIVITY_Y1_P1_METAL.md`, `docs/agentic-notes/injectivity/INJECTIVITY_Y1_P3_AMINE.md` ·
`src/oinsmiles/oin/levers.py::_HELD_OFF` — the single best place to read *why* neither lever is on.

## Open questions / for the next agent

1. **P3 under the shipped default is the top open item, and the fix is constrained, not open-ended.**
   A correct fix must preserve the metal-locked donor tag **without perturbing the canonical
   ranking**: either keep the donor **bracketed** through the sanitize, or **re-derive parity from the
   parent geometry once the write order is fixed**. Do **not** re-attempt the four-line tag copy.
   **The acceptance measurement, exactly:** with `OIN_EMIT_LOCKED_DONOR=1 OIN_CANONICAL_BODY=1`,
   assert POJJOP's enantiomer pair diverges **and** `RIFGUJ_comp_2`'s three ring-carbon tags are
   identical for structure and reflection
   (`test_locked_donor.py::TestMultiCentreDescriptor::test_flips_under_reflection`,
   `::test_three_locked_amines_all_invert_together`, `::TestRifgujRingCarbonsArePseudoAsymmetric`).
   **A fix that only passes POJJOP has not been tested.**
2. **Promoting either lever is a package, not a flag flip.** For `OIN_EMIT_METAL_CONFIG`: remove the
   `_METAL_CONFIG_TOKEN_RE` fold **in the same commit**, and take a **corpus population measurement**
   first — nobody knows how many molecules would emit `|mc:|`. For `OIN_EMIT_LOCKED_DONOR`: the
   blocker is generator-side (RDKit drops `[N@]` on the parse the generator performs, so
   selection-by-descriptor would need the descriptor carried out of band into the parsed structure
   first), plus the `OIN_CANONICAL_BODY` incompatibility in item 1.
3. **The generator side of Δ/Λ is verified on two fixtures only.** `accept_fn` reproduction was shown
   for ZUMNEC and its mirror. Two fixtures is exactly the sample size that produced the release's
   repeated wrong answers — extend to the corpus population from item 2 before trusting it.
4. **`ABIFAV_comp_0` needs a conformer-independent isomer oracle** to separate genuine (R,S) internal
   compensation from slot-assignment reflection-sensitivity. [Lane 7](LANE-07-research-residuals.md)'s
   `tools/injectivity/torsion_oracle.py` is the closest existing instrument and was built after Lane 6
   closed — try it before building anything new.
5. **`tests/fixtures/BENVOG_comp_0.xyz` hangs `get_tmc_mol`** and appears in three lanes' residues.
   Fix it or remove it from the fixture set; a hanging fixture silently truncates any whole-fixture
   survey.
6. **Do not re-derive the Δ/Λ lesson.** If a future formulation reasons about donor **positions**,
   it is formulation 1, 2 or 3 again. Δ/Λ is a property of **chelate connectivity**; the ligand
   partition has to be threaded in from the caller.
