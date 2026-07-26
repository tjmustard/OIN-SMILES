# v0.4.6 — H-faithful canonical body, and two refuted hypotheses

Branch `swimlane/v046-hfaithful`, off tagged v0.4.5. **Nothing here is merged.**

The v0.4.5 close-out ranked four items toward the accuracy goal. Two of the four turned out to be
wrong as stated, and finding that out is the main result of this branch.

## 1 · H-faithful canonical body — CORRECT, but NO measured benefit

**The mechanism was real.** `canonical_body_emit` had two writes, both plain `Chem.MolToSmiles`:
the intermediate that feeds the reparse, and — the one that matters — the **final emit whose
output becomes the body**. So `xyz2mol.py:1710` computed an H-faithful string and `:1736`
overwrote it. That is exactly how `OIN_CANONICAL_BODY` "undid" `OIN_H_FAITHFUL`, as recorded in
`levers.py::_HELD_OFF`. Both writes now route through `h_faithful_smiles`.

**The benefit was not.** A/B over the 45-molecule `Atom count mismatch` population, comparing the
OIN-implied atom count against the input XYZ (generator-free):

| arm | match | mismatch |
|---|---|---|
| `OIN_H_FAITHFUL=0`, canonical body ON | 8 | 37 |
| `OIN_H_FAITHFUL=1`, canonical body ON | **8** | **37** |

Identical. `h_faithful_smiles` guarantees only that a string re-reads with the hydrogens it was
*written* with, so whatever moves the count is not that.

### Where the count actually diverges — and what is NOT yet known

Three measurements, and two hypotheses killed along the way. Recorded in order because both wrong
turns were over-reads of small samples.

| measurement | result |
|---|---|
| `perceived_H` (parent mol) vs `input_H` (XYZ) | **agrees in 36/45** — perception is right |
| `oin_H` (implied by the string) vs `input_H` | diverges in **41/45** |
| `dH` vs count of BARE donor atoms in the string | **matches in only 4/45** |

So the divergence appears **between the perceived parent and the emitted string**, not in
perception and not in write/read fidelity.

⚠ **Hypothesis 1, WRONG:** "it lives in perception (`get_lig_mol` / `AC2BO`)". Inferred from the
flat A/B without checking where the H entered. Refuted by `perceived_H == input_H` in 36/45.

⚠ **Hypothesis 2, WRONG:** "one implicit H per severed metal–donor bond — a bare `N{1}` gains the
hydrogen the metal bond was spending". It is a real effect and it reads convincingly off examples
like `CUDBOU` (`N{1}c1ccccc1N{2}`), but it does not survive the corpus: `dH` equals the bare-donor
count in 4/45, and `dH` is **bidirectional**, spanning −36 to +14. A single donor-cut rule cannot
produce negative `dH`.

**What the distribution says instead:** the class is heterogeneous, with at least two mechanisms.
28 of 45 sit at `dH` = +1…+3 (consistent with something small and donor-adjacent), 4 are `dH` = 0
(the mismatch is not hydrogen at all), and three are large losses (−14, −16, −36) which no
donor-cut story explains — those need their own attribution, probably haptic/eta bodies or the
`[CH]`-radical writing.

**Next step is per-ATOM attribution, not another aggregate.** Walk one molecule from each `dH`
band, mapping parent atom → fragment atom → emitted token, and record which atom's H changed and
at which step. Aggregates have now produced two plausible-and-wrong answers; only per-atom
provenance will settle it.

**Why the change is kept anyway:** it removes a genuine inconsistency — one lever silently
undoing another — and it is a prerequisite for `OIN_H_FAITHFUL` ever mattering. It is **not** an
accuracy win and must not be counted as one.

**Verified:** byte-identical on all 61 fixtures against unmodified `main` (`#DONE 61` sentinel on
both arms, outputs diff-clean), and the full suite is 838 tests OK / 3 skipped / 4 xfail.

## 2 · Restoring the locked-donor tag through the reparse — REFUTED

Four lines. The correspondence is already available: `_reparse_once`'s Guard 2 proves
`donors[k] <-> new_donors[k]` (same element, same heavy degree). It makes P3 emit under
`OIN_CANONICAL_BODY`, and **POJJOP passes**.

It is still wrong. Setting a chiral tag *after* the sanitize introduces a stereocentre the
canonical ranker did not account for, which moves the canonical **write order** — and `@`/`@@` is
a parity relative to that order, not an absolute label. On `RIFGUJ_comp_2` (three Cu-bound amines
on one cyclohexane) the three ring-**carbon** tags then flip between a structure and its mirror.

The geometry says they must not. `AssignStereochemistryFrom3D` + `rdCIPLabeler` label those
carbons lowercase **`s`** — pseudo-asymmetric, a *relative* (all-cis) descriptor — and they read
`s` identically for the structure **and** its reflection:

```
atom   6  base=s  mirror=s   same
atom   8  base=s  mirror=s   same
atom  10  base=s  mirror=s   same
```

So the restoration silently rewrote stereochemistry that must not move. The lane's multi-centre
mirror guards caught it; single-centre POJJOP could not — the Y2 lesson intact. **v0.4.5's
decision to defer this rather than rush it was therefore correct.**

Reverted, with the mechanism written into `canonical_body.py::_reparse_once`,
`levers.py::_HELD_OFF` and the test module, plus a new lever-independent guard
`test_locked_donor.py::TestRifgujRingCarbonsArePseudoAsymmetric` so the four-line "obvious fix"
cannot be silently re-attempted.

A correct fix must preserve the tag **without perturbing the ranking**: keep the donor bracketed
through the sanitize, or re-derive parity from the parent geometry once the write order is fixed.

## 3 · The timeout bucket masks defects — REFUTES "the gap is mostly compute"

24 molecules (seed 42, stratified 12 `UFF_1` + 12 `g-xTB_1`) that all hit the 300 s wall in full
mode, re-run on the cheaper `--quick` path:

| outcome | n |
|---|---|
| SUCCESS | 6 (25%) |
| String mismatch | 6 (25%) |
| Atom count mismatch | 6 (25%) |
| MetalloGen failed | 6 (25%) |
| **timed out again** | **0** |

Only ~25% is compute-limited. More compute buys ~44 of 936 molecules (~4.7%), not the 174 the
timeout count implies.

⚠ The probe's `elapsed_s` values are unusable — it ran alongside the 6-shard 5k sweep at load
21–33, and `tools/v045_state.sh` says wall-clock is meaningless above ~12. Only the pass/fail
outcomes survive, and only because none of the 18 failures was a timeout.

## Revised ranking (the previous one is superseded)

1. **`OIN_BORON_CAGE`** — the only item with a measured accuracy gain: 0/36 → **34/36** encodes on
   the boron population, 0.2–4.2 s each. Costs 14 silent false passes becoming honest failures.
   Needs a full sweep, not a mid-release flip.
2. **Perception-side hydrogen** — the atom-count class, now correctly located in `get_lig_mol` /
   `AC2BO` rather than in serialization. Unscoped.
3. **MetalloGen generation failures** (~80 molecules). Unscoped.
4. **String mismatch** (~55 including timeout-masked). Unscoped.
5. **Compute** for the ~44 genuinely timeout-limited. Buys the least of any option.

## On the 100% target

Not reachable from here, and not only for effort reasons. `xyz2AC_obabel` can perceive a
genuinely **different graph** between two conformers near the covalent-radius + 0.45 Å cutoff
(`xyz2mol_local.py`), which bounds what any canonicalization can achieve — this is stated as
out-of-scope in the v0.4.5 plan and remains true. The <30 s target is much closer: already 83.6%
of succeeding molecules, median 6.9 s, p90 50 s.

---

# The <30 s tail: driver identified and quantified

Previously this target was measured (83.6% of succeeding molecules ≤30 s, median 6.9 s, max
277 s) but its cause was never attributed. It is **eta-ring winding**.

Measured over the 634 succeeding molecules of the re-baseline, by presence of an eta winding
marker (`{n>}` / `{n<}`) in the emitted OIN:

| band | n | eta-winding present |
|---|---|---|
| ≤30 s | 530 | **19.8%** |
| >30 s | 104 | **63.5%** |
| 30–60 s | 51 | 58.8% |
| 60–120 s | 30 | **76.7%** |
| >120 s | 23 | 56.5% |

A **3.2× enrichment** in the tail. This corroborates the earlier v0.4.4 finding that eta is 23% of
molecules but 35.6% of CPU, and it locates the cost precisely.

⚠ Deliberately a **string-presence count, not a timing measurement**. A 6-shard sweep was running
at load 21–33 and `tools/v045_state.sh` states wall-clock is meaningless above ~12 — the same trap
that already cost the timing half of the quick-mode probe. Marker presence is identical on an idle
and a loaded machine, so this result survives the contention.

Size is NOT the driver: the twelve slowest molecules are 73–139 atoms (35–66 heavy) against a
6.9 s median at 2–3× smaller. The relationship is not one of scale.

## The mechanism, and the fix candidate

`metallogen_adapter.py`, in the pool-sizing block:

```python
needs_winding = bool(_eta_winding_multiset(getattr(parsed, "original_oin", None)))
base_pool = ETA_SELECT_POOL if needs_winding else DEFAULT_SELECT_POOL
pool_n = max(self.ensemble_size, base_pool)
...
if needs_winding:
    uff_pool_size = max(uff_pool_size, 2 * pool_n)
```

An eta molecule pays a **wider select pool AND a 2× wider UFF pre-pool, committed UP FRONT**,
because "winding-aware selection can only pick a winding that exists in the pool".

The candidate fix is **incremental widening**: start at `DEFAULT_SELECT_POOL` and widen only if no
conformer with an acceptable winding was found, rather than paying 2× for every eta molecule
regardless. `OIN_EARLY_EXIT` (default-ON since v0.4.4) already short-circuits on the first
accepted conformer, so molecules whose winding is sampled early would stop paying for a pool they
never walk.

**NOT implemented, and the honest reason:** the payoff depends on a number this probe cannot
obtain under load — how many candidates an eta molecule typically walks before its winding
matches. If matches usually come early, incremental widening is a large win; if they usually come
late, it changes nothing and adds a retry path. That distribution needs an idle machine, so it is
the first measurement to take once the 5k sweep finishes. Guessing which regime holds is exactly
the error that produced four refuted hypotheses earlier in this release.

---

# ⚠ CAUTION ON THE BORON PROMOTION — blast radius is wider than boron

Promoting `OIN_BORON_CAGE` left exactly one suite failure out of 840, and it is not a stale
golden. `test_canonical_body::test_unparseable_body_gets_stable_raw_token`:

```
AssertionError: 'C#O' != 'RAW:C#O'
```

**Carbon monoxide now parses under boron mode.** `compare.py::_parse_fragment` gates its
valence-check-free rungs (`_NO_VALENCE`, `_NO_VALENCE_NO_KEKULIZE`) on `OIN_BORON_CAGE`, and those
rungs apply to **every fragment**, not only to boron cages. `C#O` fails the valence check and
nothing else, so with valence checking skipped it succeeds and never reaches the `RAW:` fallback.

That matters because **CO is one of the most common ligands in transition-metal chemistry**. The
promotion therefore changes fragment parsing for a large population that has nothing to do with
boron — far beyond the 34 molecules the promotion was justified on.

## What is NOT yet established

The 61-fixture byte-identity check in this branch was run **before** the promotion, with
`OIN_BORON_CAGE` OFF. There is **no byte-identity evidence for non-boron molecules under
boron-ON**. Two questions are open and both need measuring before this is merged:

1. does any non-boron fixture's emitted OIN change when the lever is on?
2. for a CO-containing complex, is the permissively parsed `C#O` mol *correct*, or is it an
   over-valent carbon that the `RAW:` fallback existed to keep out of the string?

## The narrower fix, if the answer to either is bad

Scope the valence-free rungs to fragments that actually contain a boron cage — the detector
already exists (`_is_electron_deficient_cluster`, and `_has_boron_cage` at the `xyz2mol` call
sites) — instead of enabling them globally whenever the lever is set. That keeps the measured
0/36 → 34/36 gain while removing the collateral change to CO and every other over-valent
fragment.

Recorded rather than fixed because the deciding measurement (a boron-ON vs boron-OFF byte-identity
A/B over the fixture set and a CO-containing cohort) has not been run, and this release already
contains four hypotheses that looked obvious and were refuted by exactly the measurement that was
skipped.

## Merge readiness for the boron promotion — what actually gates it

Two things were conflated. Separating them:

**NOT a gate: the byte-identity A/B.** I deferred it as "needs an idle machine once the 5k sweep
frees the cores". That was wrong. String equality is **deterministic** — the same input produces
the same OIN whether the box is idle or at load 30. Only *wall-clock* is unusable under load.
Deferring a load-independent measurement on load-dependent grounds cost real time, and it is the
mirror image of the earlier mistake of running a timing probe *during* the sweep.

**IS a gate: the merge itself.** The 5k sweep runs one subprocess per molecule, each importing
from the main checkout's `src`. Merging this branch mid-run would mean the early molecules were
measured under v0.4.5 and the later ones under v0.4.6 — a mixed-config sweep, the exact asymmetry
that manufactured v0.4.4's 11 phantom regressions. So the merge waits for the sweep regardless of
what the A/B says.

Correct sequencing, therefore:

1. run the byte-identity A/B **now** (done — load-independent);
2. let the 5k sweep finish and publish the clean v0.4.5 absolute number;
3. merge, then re-sweep to diff the boron gain against that baseline on identical molecules.

Step 2 is not optional cost: without a clean v0.4.5 number there is nothing to diff the boron
promotion against, and "34 molecules now encode" would be an unanchored claim rather than a
measured delta.

## The scoping is VERIFIED — 0 boron-free fragments affected, n=1194

The full-encoder byte-identity run was killed mid-flight (6 of 61 fixtures, no `#DONE` sentinel —
the sentinel is what made "incomplete" immediately distinguishable from "agreement"). It was
replaced with a **better** instrument rather than restarted.

My change touches exactly one predicate, `_parse_fragment`'s cage rung, so exercising that
directly is both cheaper and more targeted than comparing whole OIN strings downstream. Harvested
every distinct fragment body the corpus actually emits (from `smiles_1`/`smiles_2` across the
936-molecule re-baseline, slot markers and `_GEO` stripped) and compared the parse result with the
lever ON and OFF:

| | |
|---|---|
| distinct emitted fragment bodies | **1,194** |
| fragments whose parse result differs ON vs OFF | **56** |
| …of those, containing boron | **56** — the intended scope |
| …of those, boron-FREE | **0** — the leak is closed |

All 56 go `None → parsed` (e.g. an 8-atom `[BH]1B[BH][B@H]2...` cage, a 46-atom carborane
thioether), which is the recovery the lever exists for. `C#O` and every other over-valent
non-boron fragment is untouched.

Better than the instrument it replaced on every axis: n=1194 rather than 61, it tests the changed
code rather than a downstream proxy, it runs in seconds rather than ~80 minutes, and it is
load-independent so the running sweep cannot corrupt it.

---

# The eta "incremental widening" fix is WRONG — refuted BEFORE implementing it

An earlier section named incremental pool widening as the fix for the eta-driven >30 s tail: start
at `DEFAULT_SELECT_POOL` and widen only on failure, instead of committing a 2× UFF pre-pool up
front for every eta molecule. I was about to implement it. **It would have achieved nothing.**

`generator3d/__init__.py`'s fill loop already does this. From `_try_accept`'s own docstring:

> Returns the accepted ``Molecule`` (the one appended to ``successful_mols``) or ``None`` -- the
> SL1 early-exit hook keys off this to test the fresh conformer against ``accept_fn``.

`accept_fn` is consulted **per conformer, during the fill**, and `OIN_EARLY_EXIT` (default-ON since
v0.4.4) short-circuits the moment one is accepted. So a wide `uff_pool_size` is a *ceiling*, not a
cost: a molecule whose winding appears on attempt 3 stops at attempt 3 whether the ceiling was 10
or 40. Incremental widening would re-implement, with added retry machinery, a short-circuit that
already exists.

**So the eta cost is not pool bookkeeping — it is a LOW ACCEPTANCE RATE.** Eta molecules are
expensive because the embed seldom produces the requested ring face: nothing accepts, so the loop
runs toward `reject_budget` / `embed_time_budget` instead of exiting early. The surrounding code
says exactly this about the analogous E/Z case — AFECIZ "went from 565s unconstrained to >27 min,
all of it spent rejecting", which is why `reject_budget` exists at all.

**What would actually help, and why it is not a quick win:** raise the *probability* of generating
the right winding — constrained embedding, or seeding the ring face directly — not the pool size.
That is construction rather than selection, and this project already carries **three recorded
negative results for construction over selection**. So it is a research question with prior
evidence against the obvious approach, not an afternoon's work.

This is the sixth hypothesis refuted in this release and the only one refuted *before* being built.
Reading the code the fix would have touched cost minutes; implementing and then measuring it would
have cost hours for a null result. Worth keeping as a habit: **when a fix is named from a mechanism
rather than from a measurement, read the code it would replace first — the short-circuit you are
about to add may already be there.**

**The boron promotion now has no unmeasured risk.** What remains before it can be called an
accuracy delta is arithmetic, not safety: the 5k sweep has to publish a clean v0.4.5 baseline to
diff against on identical molecules.

---

# Lane 5, second measurement: the MAGNITUDE THRESHOLD is refuted

The permutation-invariant chirality index (commit `c72fbc1b`) genuinely solved the ordering
blocker — 6 random donor permutations give an identical value where the ordered descriptor flipped
1 → −1. That part stands.

But I also claimed "achirality falls out of the index rather than needing a planarity test", on the
strength of two synthetic controls returning **exactly** `+0.000e+00` (a perfect square, an ideal
octahedron). That claim is **wrong on real structures.** Re-measured with donors taken from
perception rather than a distance cutoff — which fixes the donor-set half of blocker 2:

| fixture | index |
|---|---|
| ZUMNEC — genuinely chiral Δ/Λ | **−4.807e-04** |
| JEGKOW — achiral square planar, only puckered | **−3.287e-04** |
| ideal square (synthetic) | +0.000e+00 |

**The two real structures are the same order of magnitude, 1.5× apart.** Crystallographic pucker in
an *achiral* complex produces a chirality index comparable to genuine helicity, so **no threshold
separates them** and any `_CHIRALITY_EPS` is arbitrary. The exact-zero cancellation is a property
of idealized coordinates only, and reading two clean synthetic controls as evidence of a general
property is the same mistake as every other refutation in this release: **a measurement that only
exercises the easy case confirmed a wrong belief.**

## What Lane 5 now needs

Not a better threshold. A **symmetry test**: does an improper operation map the complex onto itself
within tolerance? That is point-group detection with a distortion tolerance, which is a real piece
of work and the honest remaining scope of the lane.

## Lane 5 status, precisely

| property | state |
|---|---|
| permutation / relabelling invariance | **SOLVED** (index form, measured over 8 permutations) |
| invariance under proper rotation | **PROVEN** |
| inversion under reflection | **PROVEN** (exact negation) |
| donor SET determination | **SOLVED** — take perception's metal-incident bonds, not a distance ratio |
| chiral vs achiral DECISION | **OPEN** — needs symmetry detection, not a magnitude threshold |
| wired to the emit path | **NO** — and it must not be until the decision above is sound |
