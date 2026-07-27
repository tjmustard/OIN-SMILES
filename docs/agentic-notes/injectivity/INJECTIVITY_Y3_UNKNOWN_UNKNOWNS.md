# Injectivity Audit (Y3) — the unknown unknowns and the missed successes

**Status:** measurement wave — no encoder code changed. Parents:
`INJECTIVITY_Y1_OVERVIEW.md` (the confusion-matrix frame), `INJECTIVITY_Y2_FEASIBILITY.md`
(fix feasibility). Reproducers: `tools/injectivity/{missed_success_audit,uu_hunt}.py`.

Waves 1–2 worked the **false positive** cell: a lossy OIN that still passes the round trip.
Wave 3 closes the frame by attacking the other two open questions — what a round-trip
*failure* actually means, and whether there are blind spots nobody has named.

---

## Part 1 — The missed-success (false-negative) audit

### The question

A round-trip FAIL is normally read as "something is wrong." With what, though? The program's
starting calibration was that such a failure is *usually the generator* — MetalloGen not yet
fully accurate, or the job timing out — not the notation. That is testable, and the v0.4.4
regression sweep already holds the evidence, so no re-run was needed.

Every failing row was attributed to exactly one cause, grouped by what it says about the
**notation**.

### The result: 77.8% of failures never test the notation

Over 3917 molecules (2080 passing, 1837 failing):

| cause | count | share of failures | informative about the OIN? |
|---|---:|---:|---|
| generator timeout | 1238 | 67.4 % | no — untested |
| encoder refused the input | 241 | 13.1 % | encoder coverage |
| generator produced nothing | 191 | 10.4 % | no — untested |
| canonicalization noise | 123 | 6.7 % | test too strict |
| output names a different isomer | 44 | 2.4 % | ambiguous |

The median failing molecule ran **300.3 s against a 300 s budget**. Two thirds of the failure
mass is a stopwatch expiring.

### What follows

**The calibration is confirmed**, and it has sharper consequences than "the generator is
imperfect":

1. **The headline round-trip pass-rate is substantially a measure of generator throughput**,
   not of notation quality. Reading a pass-rate movement as a change in losslessness is a
   category error unless the failure mix is held fixed.
2. **Compute buys apparent accuracy.** With most failures sitting at the wall clock, raising
   the time budget converts failures into measurements without touching the encoder. Any A/B
   that also changes runtime moves the pass-rate for reasons unrelated to the change under
   test — precisely the config-asymmetry artefact behind the v0.4.4 regression sweep, where
   all 11 apparent "regressions" turned out to be 300 s timeouts.
3. **Confirmed evidence of a lossy OIN in the failure population: zero.** Not because none
   exists — Y1 proved lossy encodings do — but because *this instrument cannot see them*.
   The 44 divergent-isomer cases are deliberately **not** scored against the notation: a
   divergent isomer is equally consistent with the generator having built the wrong thing and
   with the OIN having licensed it, and the sweep cannot separate those.

The last point is the reason the whole program exists. Losslessness must be measured by
generator-free collision probes, as in Waves 1–2. The round-trip pass-rate does not answer
that question and never did.

Reproduce: `PYTHONPATH=$PWD/src python -m tools.injectivity.missed_success_audit`
(reads the sweep's `bucket_report.json`; no generation, so it is cheap and deterministic).
Output: `results-injectivity-y2/missed_success_audit.{json,md}`.

---

## Part 2 — The unknown-unknown hunt

### The method, and the trap it has to avoid

Same generator-free mirror twin as Wave 1, but triaged so only the *unexplained* residue
surfaces. For each structure: mirror it, encode both, and ask whether the encoder separates
them. If it does not, ask whether a **known** cause explains the chirality — an RDKit CIP
stereocentre, or one of the P1/P2/P3 axes from the configurational oracle.

The trap is the one Wave 1 already hit: the geometric oracle is a **rigid** superposition
test, so on flexible molecules it over-reports chirality — a floppy achiral molecule's mirror
is a non-superimposable *conformer* and reads as chiral. Wave 1 refused to publish a rate for
this reason, and simply sorting by a rigidity proxy is not enough.

### The discriminator: InChI

**InChI is configuration-based and conformation-independent.** That makes it a sharp knife
here. If InChI separates the two twins while the OIN does not, the encoder has lost something
standard cheminformatics keeps — with *no* conformational confound available to explain it
away. That is a rigorous finding, not a candidate.

So the residue splits in two, and the two halves are not equally strong:

- **CONFIRMED blind spot** — collapse, and InChI separates the twins.
- **Ambiguous residue** — collapse, and InChI agrees. Two very different things live here:
  - *conformational chirality*, where the mirror is a different conformer of the same isomer
    and collapsing it is **correct** (that is conformer invariance, which the encoder has
    already been shown to have);
  - *jointly-blind axes*, because InChI is itself blind to metal Δ/Λ and to atropisomerism —
    so P1/P2/P3 sit in this bucket, and a genuinely new axis that InChI also misses would
    hide here too.

### The result: no new rigorous blind spot in 299 structures

| verdict | count | share | meaning |
|---|---:|---:|---|
| chiral, encoder separates | 130 | 43.5 % | correct injectivity |
| collapse, named P1/P2/P3 axis | 65 | 21.7 % | known blind spot |
| collapse, InChI agrees | 70 | 23.4 % | ambiguous residue |
| collapse, CIP stereocentre present | 25 | 8.4 % | modelled stereo type |
| achiral, encoder collapses | 9 | 3.0 % | correct invariance |
| **collapse, InChI SEPARATES them** | **0** | **0.0 %** | **confirmed blind spot** |

Two things are worth reading off this.

**No new axis was found.** Across 299 structures there is not a single case where InChI keeps
a configurational difference the OIN loses. On this sample the encoder is at least as
discriminating as standard cheminformatics. That is a negative result for the sample size,
not a proof of injectivity — but it is the strongest statement this instrument can make, and
it means the named axes are not the tip of a large unnamed iceberg.

**The named axes dominate the real collapses.** 65 of 299 (21.7 %) collapse for a P1/P2/P3
reason. The Wave-1 targeting was pointed at the right things: the blind spots that exist are
overwhelmingly the ones already catalogued, which is why fixing P2 was worth the effort.

### Manual triage of the ambiguous residue

The top-ranked rigid candidates were opened by hand, and they are conformational:

| structure | what it is | verdict |
|---|---|---|
| `EDOQIZ` | unsubstituted biphenyl phosphine on linear Au | freely-rotating biaryl; the stereogenicity gate already refuses to call this axis stereogenic — **correct collapse** |
| `WAVGOS` | bis-NHC on linear Au, propargyl arms | flexible arms; InChI agrees — conformational |
| `PERPIO` | Ni with alkynyl + phosphine | InChI agrees; mirror RMSD 1.07 Å, low — conformational |

All three have **identical InChI** for base and mirror. The encoder is right to collapse them.

### Status

Separating *conformational* from *jointly-blind* inside the ambiguous residue needs a
torsion-aware configurational test rather than a rigid one — enumerate conformers and ask
whether any of them superimposes on the mirror by a proper rotation. That is the remaining
instrument, and it is what would let this hunt scale past manual triage.

Reproduce: `PYTHONPATH=$PWD/src python -m tools.injectivity.uu_hunt --n 300`.
Output: `results-injectivity-y2/uu_hunt.{json,md}`.

---

## Part 3 — The plan's other UU candidates

The Y1 target map listed two further suspects, both flagged as "key blind". They are
key-level claims, so the key — a pure string function — is the right instrument to test them
with directly.

### Linkage isomerism — REFUTED (not a blind spot)

A ligand bound through a different donor atom is a genuinely different compound whatever the
coordination slots do, so the comparison is sound. The key **distinguishes** both classic
cases:

| pair | key |
|---|---|
| thiocyanate N-bound vs S-bound (`N{0}=C=S` / `S{0}=C=N`) | differs |
| nitro N-bound vs nitrito O-bound | differs |

The vertex colouring carries donor-atom identity, so linkage isomers do not collide. This
target-map entry is refuted.

### Symmetric / asymmetric donor swap — INCONCLUSIVE, and a methodological warning

An obvious next probe is to swap which coordination slot each donor of an unsymmetric chelate
occupies, and see whether the key folds it. Done naively — by hand-writing the two OIN strings
— the key does fold them, which looks like a fresh blind spot.

**That conclusion would be wrong, and it is worth recording why.** For a square-planar complex
with two identical ancillary ligands, swapping slots 0↔1 is related to the original by a
reflection, and a square-planar complex is planar and therefore achiral — so the reflected
arrangement is the *same* isomer, and the key is right to fold it. The probe never established
that its two strings denote different isomers. (The octahedral variant was worse: it placed a
short bidentate across *trans* slots, which is geometrically impossible.)

This is exactly the trap Y1 documented in the pre-existing
`test_isomer_divergence.py::test_metal_stereo_raw_only` — **hand-fed strings prove nothing
about isomerism unless something independent certifies the two are distinct isomers.** The
rule the audit adopted for encoder claims applies to key claims too: drive the comparison from
real geometry, with the oracle establishing distinctness. Doing that for donor swap needs real
linkage/positional isomer pairs from the corpus, and is the natural next probe.

Status: **undetermined**, not refuted and not confirmed.

## Where the confusion matrix stands now

| | Round-trip PASS | Round-trip FAIL |
|---|---|---|
| **OIN truly correct** | True Positive | **False Negative — measured (Part 1): 77.8 % of failures never test the notation** |
| **OIN wrong / lossy** | **False Positive — measured (Waves 1–2): 3 named axes, P2 now fixed** | True Negative |

Both off-diagonal cells now have instruments. Neither is measurable by the round-trip test
itself, which was the thesis this program set out to establish.
