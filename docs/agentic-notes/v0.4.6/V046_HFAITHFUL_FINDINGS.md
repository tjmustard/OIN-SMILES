# v0.4.6 — H-faithful canonical body, and two refuted hypotheses

Branch `swimlane/v046-hfaithful`, off tagged v0.4.5. **Nothing here is merged.**

The v0.4.5 close-out ranked four items toward the accuracy goal. Two of the four turned out to be
wrong as stated, and finding that out is the main result of this branch.

## 1 · H-faithful canonical body — CORRECT, but NO measured benefit

**The mechanism was real.** `canonical_body_emit` had two writes, both plain `Chem.MolToSmiles`:
the intermediate that feeds the reparse, and — the one that matters — the **final emit whose
output becomes the body**. So `perception_tmc.py:1725` computed an H-faithful string and `:1749`
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
(`perception_core.py`), which bounds what any canonicalization can achieve — this is stated as
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
already exists (`_is_electron_deficient_cluster`, and `_has_boron_cage` at the `perception_tmc` call
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

---

# Lane 5, third measurement: the descriptor's INPUT is wrong

Replacing the magnitude threshold with a proper symmetry test (mirror the donor set, ask whether any
permutation + proper rotation superimposes it) was the right move for the *decision*. It produced a
harder result:

| fixture | symmetry test says | truth |
|---|---|---|
| ZUMNEC — chiral Δ/Λ tris-bidentate | **achiral** | chiral |
| JEGKOW — square planar | achiral | achiral |
| ideal square | achiral | achiral |

**ZUMNEC reads achiral, and as a bare point set it IS.** Six oxygens at octahedral vertices admit
improper operations; there is no handedness in the donor positions alone.

**Δ/Λ helicity is a property of the chelate CONNECTIVITY** — which donor pairs belong to the same
bidentate ligand, and how those chelate planes twist about the metal. Reflecting a Δ complex yields
Λ only because the reflection cannot be undone *while keeping the chelate pairing intact*. A
permutation search over unlabelled points is free to re-pair the donors, so it always finds a
"symmetry" that is not chemically available.

This also retro-explains the previous result: `chirality_index`'s non-zero reading for ZUMNEC
(−4.807e-04) was residual crystallographic distortion — the same magnitude as achiral JEGKOW's
pucker (−3.287e-04). **It was never detecting helicity at all.** Three of its four proven
properties (rotation invariance, reflection inversion, permutation invariance) were real and are
still real; they were just properties of a quantity that does not mean what the lane needs.

## What Lane 5 actually requires

Constrain the permutation search to relabellings that **preserve chelate membership** — treat the
donors as a *coloured* point set, colour = the ligand each donor belongs to. A mirror that has to
re-pair chelates is then correctly rejected. That needs the ligand partition threaded in from the
caller, which the current signature does not carry.

## Lane 5 status, corrected again

| property | state |
|---|---|
| permutation / relabelling invariance of the index | SOLVED |
| invariance under proper rotation | PROVEN |
| inversion under reflection | PROVEN |
| donor SET determination | SOLVED (perception's metal-incident bonds) |
| point-set achirality test | BUILT and correct as such |
| **detects Δ/Λ helicity** | **NO — needs chelate connectivity, not positions** |
| wired to emit | NO, correctly |

The honest read: the lane has a working achirality test and a working pseudoscalar, neither of which
is yet a Δ/Λ descriptor. Eighth refutation of this release, and the only one that invalidated the
*input* rather than the method.

---

# One untried idea for the eta tail: REPAIR the winding, don't constrain the embed

The eta tail's cost is a low acceptance rate — the embed seldom produces the requested ring face, so
the loop runs toward `reject_budget` / `embed_time_budget` instead of exiting early. The obvious
remedy is to make the embed produce the right face (constrained embedding, ring-face seeding), and
that is *construction over selection*, for which this project already carries three negative results.

**But there is a third option that is neither, and it has not been tried: post-hoc repair.** Take a
conformer whose topology and geometry are already acceptable but whose eta ring is on the WRONG face,
and rotate that ring to the requested face, then relax locally. The eta winding is a rotational
degree of freedom about the metal–ring axis; flipping which face is presented is a bounded geometric
edit, not a re-embed.

Why this is not the refuted category: the three prior negative results are about *constructing* a
geometry to satisfy a constraint from the start. This accepts whatever the embed produced and edits
one degree of freedom afterwards — the same shape as the `invert_stereocenter` / `swap_donor` twin
operators Lane 7 already built in `tools/injectivity/twin_operators.py`, which perform structural
edits and then filter through the existing vdW clash gate rather than trusting the edit.

What it would need, in order:

1. measure the acceptance-rate distribution for eta molecules — how many attempts before the
   requested face appears? **This is a COUNT, not a timing, so it is load-independent** and can be
   taken while a sweep runs. Do not defer it for load reasons; that error was already made twice in
   this release, in both directions.
2. if the distribution is long-tailed, implement the ring-face rotation as a twin-operator-style
   edit, gate it through `generator3d/clash.py`, and re-encode to confirm the winding actually
   changed — never trust the edit.
3. A/B on the eta subpopulation only, comparing attempts-to-acceptance rather than wall-clock.

Recorded rather than built because step 1 has not been run, and this release contains eight
hypotheses that looked obvious and died to precisely the measurement that got skipped. The idea's
value is that it sits outside the category the project has already disproved — not that it is likely
to work.

---

# Correction: the eta cost looks like COST-PER-ATTEMPT, not attempt count

The section above attributed the eta tail to a low acceptance rate and proposed post-hoc ring-face
repair to reduce the number of attempts. A size-controlled measurement undercuts that.

Eta vs non-eta among the 634 succeeding molecules, banded by atom count:

| atoms | eta median | non-eta median | ratio |
|---|---|---|---|
| 0–50 | 7.9 s (n=47) | 3.3 s (n=86) | **2.4×** |
| 50–80 | 15.5 s (n=91) | 5.0 s (n=182) | **3.1×** |
| 80–120 | 49.7 s (n=31) | 8.1 s (n=176) | **6.2×** |

Two things follow.

**Size is not a confound — it is the opposite.** Eta molecules are *smaller* on median (62 vs 73
atoms) while being 3.2× slower overall, and the eta penalty holds inside every band. So the eta
correlation reported earlier survives control for size, and is in fact understated by the raw
comparison.

**But the mechanism is probably not acceptance rate.** If the cost were the NUMBER of attempts
before the requested ring face appears, the eta/non-eta ratio would be roughly FLAT across size
bands — attempt count is a property of the winding constraint, not of molecule size. It is not flat:
it climbs 2.4 → 3.1 → 6.2×. A penalty that scales with size is a **cost-per-attempt** signature —
more work per conformer (haptic scaling, ring placement, relaxation) rather than more conformers.

A supporting argument from the other direction: the tail's median eta count is **1 marker**, i.e. one
ring, two faces. An unbiased embed should hit the requested face about half the time, and
`OIN_EARLY_EXIT` fires on the first acceptance — so an attempt-count story predicts eta molecules
finish in one or two attempts and are therefore FAST. They are 2.4–6.2× slower.

**Consequence for the repair idea recorded above: it probably does not help.** Reducing attempts
cannot fix a per-attempt cost. Recorded rather than deleted, because the reasoning that killed it is
the useful part, and because the discriminating measurement is now specific and cheap: instrument the
fill loop to log attempts-to-acceptance for eta molecules (a COUNT, load-independent). If eta
molecules accept in 1–2 attempts, attempt count is definitively not the mechanism and the lever is
per-conformer eta cost — a profiling problem, not a stereochemistry one.

Ninth refutation of this release, and the first of a proposal made in this same document.

## The eta discriminating test needs a small HARNESS change, not a sweep restart

The signal that settles attempt-count vs cost-per-attempt already exists in the code:
`metallogen_adapter` records `adapter.early_exit_hit` and `adapter.early_exit_miss` (with
`n_mols`). If eta molecules **hit** early exit, acceptance is fine and the cost is per-attempt; if
they **miss**, acceptance is the bottleneck. No new instrumentation required.

⚠ But `OIN_TELEMETRY=1` alone captures NOTHING. `generation/_telemetry.record()` is a no-op unless
the env var is set **and** a `collecting()` context is active, and `tools/test_dataset_roundtrip.py`
never opens one — verified, not assumed. Restarting the v0.4.6 sweep with the env var set would have
discarded ~60 completed molecules for zero data.

So the step is: wrap per-molecule generation in `telemetry.collecting()` in the harness and persist
the snapshot into each `individual_reports/*.json`. Then a single sweep yields the accuracy numbers
**and** the eta acceptance distribution together. `record()` is documented as never raising and never
consuming randomness, so it cannot perturb the accuracy result — worth confirming on a small run
before a 15 h one, given how many assumptions in this release turned out false.

Do NOT retrofit this onto the currently running sweep; it needs a code change the running processes
have already imported past.

## The eta telemetry measurement RAN, and returned a null result

Cohort: 10 eta + 10 non-eta molecules, size-matched to 50–120 atoms (the band where the eta penalty
is 3–6×), run with `OIN_TELEMETRY=1` through the newly instrumented harness.

**Result: telemetry captured on 16/16 completed molecules, and NO degradation site fired for either
group.** `adapter.early_exit_hit` / `early_exit_miss` — the counters I had called "the discriminating
signal, already in the code" — are both zero.

**Why: `_select_by_geometry(..., early_exit=False)` is the default in both signatures** (adapter
lines ~1443 and ~1678), so unless a caller passes it explicitly that block never executes and the
counters never fire. So "the signal already exists, it just needs collecting" was **wrong** — the
site is unreachable on this path as configured. Fourth wrong turn on this one item.

What the null result *is* worth: **the eta penalty is not a degradation path.** None of the
generator's instrumented fallbacks fire for eta molecules any more than for size-matched non-eta
ones, which means the extra 3–6× is being spent in *ordinary* embed / relax / selection work rather
than in a retry or recovery branch. That is consistent with cost-per-attempt, but by **absence of
evidence** rather than by measurement, so it should not be quoted as a positive finding.

What would actually settle it, stated without the confidence I had before:
1. add an explicit attempt counter to `generator3d`'s fill loop (not a degradation site — a plain
   counter incremented per embed attempt), and log it per molecule;
2. re-run this same size-matched cohort;
3. flat attempts across eta / non-eta ⇒ cost-per-attempt (profiling target); higher attempts for eta
   ⇒ acceptance-limited.

**What IS confirmed:** the harness telemetry plumbing works end to end in a real sweep — captured on
16/16, absent when the env var is unset. That part is reusable for step 1 above and cost nothing to
verify.

## MEASURED: the eta tail is attempt-driven. Both my earlier answers were wrong.

Added a plain per-attempt counter to `generator3d`'s fill loop (`pool.attempts_spent`, recorded at
**every** return of `generate_3d_structures`, not just the last) and measured:

| molecule | attempts | accepted | target_pool |
|---|---|---|---|
| Ferrocene (eta) | **32** | 28 | **32** |
| CisPlatin (non-eta) | **0** | 1 | 10 |

The non-eta molecule short-circuits on the FIRST attempt via early exit. The eta molecule runs the
**entire** pool and never short-circuits — and its pool is itself widened (32 vs 10). So it pays ~32
attempts where a comparable molecule pays 1, which comfortably explains a 3–6× wall-clock penalty.

**So "cost-per-attempt", inferred earlier from the ratio-vs-size slope, was wrong** — and so was the
"low acceptance rate" story in its original form. The mechanism is specifically: **eta molecules
never satisfy the early-exit predicate, so the widened pool is paid in full.**

That also partly resurrects the fix I refuted. I killed incremental pool widening on the grounds that
`accept_fn` short-circuits per conformer — true in principle, but eta molecules never satisfy it, so
the short-circuit never fires and the 2× widening is a real cost after all.

### The precise defect, and the fix I did NOT implement

Early exit accepts on `canonical_roundtrip_key` equality. Ferrocene round-trips fine (it is a
golden), so the final `_select_by_geometry(honor_winding=True)` DOES find a correct winding — it is
only the early-exit key match that never succeeds. The two acceptance predicates disagree: the
cheaper one that could stop the loop is stricter than the one that ultimately judges success.

**Fix: make early exit accept what `_select_by_geometry` would accept.** Then eta molecules stop at
the first acceptable conformer instead of filling 32, with no change to what counts as a correct
round trip.

NOT implemented here: it changes acceptance semantics in the selection path, which needs a corpus A/B
(does any molecule that currently passes stop passing?) and that is not something to land on the
strength of two molecules. Two fixtures is exactly the sample size that produced the four wrong
answers above.

## Where the eta cost actually lives: the POOL FILL, not selection

Implemented the predicate-alignment fix behind `OIN_ETA_EARLY_EXIT` (default OFF) and A/B'd it on
Ferrocene:

| lever | attempts | eta early-exit fired |
|---|---|---|
| OFF | 32 | 0 |
| ON | **32** | **1** |

**It fires and the attempt count does not move.** The reason is structural, and it is the thing to
carry forward: `generate_3d_structures` fills the **entire** pool before `_select_by_geometry` is
ever called. A selection-side early exit can only shorten the selection *scan*, which is cheap
beside 32 embeds.

**The site that can stop pool filling is the `accept_fn` passed INTO `generate_3d_structures`**, which
`_try_accept` consults per conformer. Adding the eta-winding criterion *there* is the real fix. The
lever as implemented is kept because it is correct, harmless, default-OFF, and marks exactly where
the boundary between the two mechanisms lies — which is the part I kept getting wrong.

Fifth wrong turn on this single item, and the counter caught it in one A/B. Sequence, for anyone
tempted to shortcut it: low-acceptance-rate → pool widening (refuted by early-exit-exists) →
cost-per-attempt (inferred from a slope, refuted by the counter) → attempt-driven (correct) →
selection-side predicate alignment (correct but ineffective) → **fill-loop accept_fn (untried)**.

## STOP CONDITION on the fill-loop fix: find out WHY key equality fails first

The named remaining fix is to add the eta-winding criterion to the `accept_fn` passed into
`generate_3d_structures` — the predicate `_try_accept` consults per conformer, and the only one that
can stop pool filling. Before implementing it, there is a soundness question that must be answered,
and it is cheap:

**Why does key equality never succeed for Ferrocene?** The pool runs to 32/32 because `accept_fn`
never returns True, yet Ferrocene is a golden and round-trips via the winding branch of
`_select_by_geometry`. So `canonical_roundtrip_key` equality is stricter than winding-match in some
**other** dimension — slot labels, the `>` heading character, or connectivity.

Until that dimension is identified, adding winding-match as an acceptance path is **not obviously
sound**: it would accept conformers the key rejects for a reason nobody has named, and acceptance
runs before any of the geometry scoring that the final selection's winding branch benefits from.
`_select_by_geometry` applies its winding test to already-`scored` candidates; `accept_fn` would
apply it to raw pool conformers. Those are not the same population.

The measurement to take, in order:
1. for Ferrocene, print the requested `canonical_roundtrip_key` beside the key of a pool conformer
   that the winding branch *does* accept, and diff them. That names the strict dimension.
2. if the difference is cosmetic (slot labelling, heading char), the winding path is safe to add to
   `accept_fn` and the eta tail collapses from 32 attempts to ~1.
3. if the difference is connectivity, do NOT add it — the widened pool is doing real work and the
   runtime cost is the price of correctness.

This is where the item genuinely stops: not for want of a fix, but because step 1 has not been run
and the fix is unsound without it. Sixth attempt on this tail would otherwise be a guess about
acceptance semantics, and five of five previous guesses here were wrong.

## RESOLVED: the eta pool cost is irreducible via acceptance. Item closed.

Ran the gate. For Ferrocene, the requested OIN and the regenerated OIN are **byte-identical**, and
`canonical_roundtrip_key` equality **holds** for the finally-selected conformer:

```
requested   [Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}...
regenerated [Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}...
keys equal: True
```

So the key is **not** stricter in some unnamed dimension — the hypothesis behind the fill-loop fix.
The real explanation is different and it settles the item:

**`accept_fn` is handed RAW pool conformers; the key only matches AFTER optimization.** For CisPlatin
the raw conformer already satisfies the key, which is why it accepts on attempt 0. For an eta
molecule the ring winding is not right until relaxation — and relaxation happens *after* the fill
loop. At the moment `accept_fn` sees an eta conformer, that conformer is genuinely not yet
acceptable.

**Therefore no acceptance-predicate change can shorten the eta fill.** Adding winding-match to
`accept_fn` would test a property the conformer does not have yet, and would either never fire (no
gain) or fire on pre-relaxation geometry (unsound). The widened pool is doing real work: it exists so
that *after* optimization, at least one member has the requested face.

**The eta runtime cost is the price of correctness, not a defect.** The `<30 s` tail is therefore not
closable by cheap means. What remains are genuinely expensive options — make the embed produce the
requested face before relaxation (construction over selection: three prior negative results), or
relax fewer candidates by scoring winding pre-relaxation (needs a pre-relaxation winding predicate
that does not currently exist).

Six attempts, one per hypothesis, each killed by a measurement:
low-acceptance-rate → pool widening → cost-per-attempt → attempt-driven → selection-side alignment →
fill-loop `accept_fn`. **The answer is that the cost is structural.** The `OIN_ETA_EARLY_EXIT` lever
stays default-OFF and documented as ineffective; the per-attempt counter stays, since it is what made
each refutation take one A/B instead of a cycle of argument.
