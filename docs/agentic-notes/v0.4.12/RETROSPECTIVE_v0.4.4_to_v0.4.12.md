# Retrospective — has OIN-SMILES gotten more accurate, or faster, since v0.4.4?

<!-- artifact-url: https://claude.ai/code/artifact/9f3a0c17-8f85-4ee6-922b-df599a51978e -->

**Scope:** v0.4.4 (2026-07-23) through v0.4.13 (2026-07-28). Ten releases.
**Baseline commit:** `main` @ `9fbdf678`, tag `v0.4.13`, `pyproject` 0.4.13.
⚠ `main` moves under this project — re-read the tip before trusting any SHA here.
**Method:** committed evidence only. No sweep, benchmark or A/B was run for this report.
**Refresh it with:** `/release-retrospective <version>` — see [Appendix A](#appendix-a--how-to-update-this).
**Published page:** <https://claude.ai/code/artifact/9f3a0c17-8f85-4ee6-922b-df599a51978e> · GitHub Pages source on the `gh-pages` branch.

---

## ELI5

Think of the project as a machine that turns a 3D molecule into a short piece of text, and back
again. It works if the text you get out the second time is *exactly* the text you put in.

For two releases the machine genuinely got better. Then someone checked the scoreboard and found
it was cheating: when the machine was asked "did you rebuild the molecule right?", it had been
grading its own homework — it looked at what it *meant* to build instead of what it *actually*
built. Once that was fixed, the score dropped from 83 out of 100 to 72 out of 100. **The machine
didn't get worse. The scoreboard got honest.**

Then five releases in a row changed nothing about how the machine works. That sounds bad, but it
isn't sloppiness — it's the opposite. Each one built the thing it was supposed to build, measured
it properly, and in most cases discovered the idea was wrong and shelved it. One built a fix worth
8 points and then found it worked by *throwing away* information — quietly turning left-handed
molecules into right-handed ones. It was switched off.

**The sixth one finally moved the score up: 72 out of 100 to 76.** It did that by taking the
switched-off fix, adding the safety check that tells the good half from the harmful half, and
proving on two separate samples that the harmful half was gone. It also cancelled the expensive
re-measurement it was supposed to run — because re-running a machine that makes slightly different
choices each time would have added noise to the very number it was trying to read, and there was a
cheaper check that gave an exact answer.

So: **more accurate than v0.4.4, yes. More accurate than v0.4.6 — now yes, and honestly measured
this time. And much harder to fool.**

---

## ELI20

Two separate questions, two separate answers.

### Accuracy

`byte_exact` — the share of molecules whose generated 3D structure re-encodes to a byte-identical
OIN string — improved materially across v0.4.4, v0.4.5 and v0.4.6, stopped moving for five
releases, and **moved up for the first time at v0.4.13**.

| release | default-path accuracy change | key number |
|---|---|---|
| **v0.4.4** | **Real gain** — `OIN_EARLY_EXIT` + OIN-direct assembly promoted | worst-cohort byte-exact **44.7% → 60.5%** |
| **v0.4.5** | **Real gain** — six canonicality levers promoted together | byte-stability under rotation/renumbering **58.1% → 69.6%**; key instability **60 → 16** molecules |
| **v0.4.6** | **Real, small** — `OIN_BORON_CAGE` promoted | boron cages **0/36 → 34/36** encoding |
| **v0.4.7** | **Zero, by design** — five lanes, four refuted their own premise | bucket report byte-for-byte unchanged |
| **v0.4.8** | **Zero behaviour change. The measurement changed.** | `byte_exact` **82.80% → 72.46%** |
| **v0.4.9** | **Zero** — confirmed live | 296/300 agreement, **0/300** encoder drift |
| **v0.4.10** | **Zero, gated** | ARM 1 **62/62**, ARM 2 **90/90** byte-identical |
| **v0.4.11** | **Zero** — built a +7.86-point fix and refuted it | fold collapses enantiomers in **221 of its own 393 gains** |
| **v0.4.12** | **Zero** — parity filter works, ships default-OFF | mirror audit **19 → 0** collapses; surviving gain **+3.42 pts** |
| **v0.4.13** | 🎉 **Real gain — the first upward move in the project's history** — fold + parity veto promoted together | `byte_exact` **72.46% → 75.88%** (+3.42 pts, 171 molecules, **0** in a bad direction) |

The v0.4.8 drop is the one that needs explaining, and it is not a regression. The harness had been
scoring a round trip with `get_oin_string(gen_result.mol, coords)` — the *generator's own bond
graph*. That is circular: `gen_result.mol` is precisely the artifact that would have to be wrong
for the test to fail. `FIYHUT_comp_0` ships both cyclopentadienyl rings 0.85 Å off the iron — ten
bonded carbons down to zero — and scored a byte-exact pass. Re-deriving the verdict from the
generated coordinates alone moved 613 molecules down and 30 up. **72.46% and 82.80% describe the
same underlying performance; only one of them is true.**

The v0.4.13 rise is the mirror image and deserves the same scrutiny, because **this project has
already raised `byte_exact` once by deleting information**. The donor fold moves 393 molecules in
one direction with the comparison key untouched — and collapses enantiomers in 221 of them, which
neither `byte_exact` nor the key can see. What makes +3.42 trustworthy is not the transition
matrix: it is that the reflection-parity veto zeroes the collapse on **two independent draws**
(19 → 0 cat-only, 33 → 0 mixed) while the achiral population stays **unmoved** in both, and that
every collapsing molecule is accounted for as a *separated pair* rather than a decline
(73 + 19 = 92; 83 + 33 = 116). **171 survivors was also reached from the opposite direction** —
v0.4.11 bounded the safe set from above at "at most ~172" by counting damage.

### Speed

Three real default-path wins since v0.4.4, and **no current corpus-wide number**.

| release | default-path speed change |
|---|---|
| **v0.4.4** | `OIN_EARLY_EXIT` default-ON — **~5×** on the worst cohort |
| **v0.4.5** | Removed a duplicate **48–57 s** re-encode running once per rejected conformer (unconditional). First encoder profile ever taken: `AC2BO` is **99.8%** of a slow encode |
| **v0.4.6–v0.4.9** | **None on the default path.** The two big candidates were held OFF on measurement, not caution |
| **v0.4.10** | Deleted a discarded `.index()` scan that eigendecomposed a Coulomb matrix per candidate — **on by default, no lever**: `CAHQEJ_comp_0` **−32.9%**, `FOSNEI_comp_0` **+0.3% (nil)** |
| **v0.4.11–v0.4.13** | None claimed. v0.4.13's promotion is **generator-neutral by measurement** — 9669 corpus strings, 1019 moved by the fold, **0 comparison keys changed** — so `accept_fn` returns bit-identical conformers and runtime cannot have moved |

The headline speed figure — **994/5000 = 19.88% of molecules over 30 s, median 7.19 s** — comes
from the v0.4.6 sweep and **has not been re-measured since**. Every speed number after it is a
per-molecule A/B, and every one of those is **bimodal**: the same change measures −32.9% on one
molecule and nothing on the next one over. Goal B (`max(elapsed_s) < 30 s`) is **not delivered**.

### What actually improved most

Neither number. What improved most is the project's ability to tell a real result from a fake one.
Between v0.4.4 and v0.4.13 the unit suite went **551 → 993**, and the instruments added in that
window — the honest re-score, the corpus encoder-identity gate, the two-arm byte-identity gate, the
coordination-integrity check, the mirror audit — are what caught the 10.34-point inflation, the
59%-false-positive cohort, the dead gate arm, and the enantiomer collapse. Five releases that moved
no points each removed a way of being wrong.

**v0.4.13 is the return on that.** The +3.42 it banks is the fix v0.4.11 built and refused; it
could only be shipped once an instrument existed that could see the damage. The same release then
caught **four instruments printing plausible nothing** — including one that reported the *refuted*
number under the correct heading, and a byte-identity gate whose PASS was meaningless because
**0 of its 62 fixtures were among the molecules the change moved**. The lesson generalises past
this project: *a gate that cannot see your change is not evidence about your change* — state its
coverage of the moved population, not its verdict.

---

## The two questions, answered directly

> ### Is it more accurate?
>
> **Yes — and as of v0.4.13, measurably so on the honest metric: 72.46% → 75.88%.**
>
> The long answer is the interesting one. From v0.4.6 to v0.4.12 the default path did not change a
> single answer, and the reported figure *fell* from 82.80% to 72.46% at v0.4.8 because the
> measurement was corrected, not because the software regressed. **v0.4.13 is the first release in
> the window to move the number up**, and it moved it by shipping a fix that had been built,
> measured, refuted and shelved two releases earlier — once the veto that separates its safe half
> from its harmful half could be proven on two independent draws.
>
> **The 24.12 points that remain are now decomposed by mechanism, not just by bucket.** Of the 767
> genuine failures, **280 are the generator returning a structure with ligands off the metal** —
> a one-site guard, not a capability limit.

> ### Is it faster?
>
> **On the molecules that hit the specific costs removed, yes — up to −86.7%.** Corpus-wide,
> **unknown since v0.4.6**, because no full sweep has been run since — and v0.4.13 deliberately
> declined to run one. Three default-path optimisations landed (v0.4.4, v0.4.5, v0.4.10); the
> largest measured single win (`VAFMIA_comp_0`, 81.89 s → 10.87 s) sits behind a lever that ships
> **off**. Goal B (`max(elapsed_s) < 30 s`) is **not delivered**.

---

## How to read the numbers in this document

Three traps, all of which this project has already fallen into once and documented:

1. **The bucket reports are not a time series.** They are drawn from four different cohorts —
   6719 (v0.4.2 capstone), 3917 (v0.4.4 regression), 936 (v0.4.5 rebaseline), 5000 (seed-42
   sweep). Every accuracy figure below carries its N and its cohort. **Do not draw a line between
   two of them.** The only genuine like-for-like comparison in this whole document is
   scored-vs-honest at v0.4.8, which classifies *the same 5000 reports* two ways.

2. **`metrics.elapsed_s` is nested and is a SUM.** Read from the top level it silently yields `0`.
   It also accumulates up to three separately SIGKILLed harness attempts, so the retired headline
   "max 759.9 s against a 300 s budget" was arithmetic on a sum — all 4658 single-attempt rows in
   the 5k sweep finish within **0.2 s** of their cap. The harness enforces to ε ≈ 0.2 s.
   *(`v0.4.9/ELAPSED_S_IS_A_SUM_v0.4.9.md`)*

3. **Suite count is a rigour proxy, not an accuracy metric.** 551 → 988 tests means the project
   can detect more kinds of wrongness. It does not mean the notation got better.

---

## Per-release detail

Each section: what shipped · accuracy · speed · what was refuted · what it cost.

---

### v0.4.4 — accuracy + the instrument to measure it

**What shipped.** Six parallel worktree swimlanes (SL0–SL5). The centrepiece is a **fac/mer-aware
canonical round-trip key** (`oin/compare.py::canonical_roundtrip_key`) that separates a genuine
*fac*↔*mer* miss from benign slot relabeling, plus `tools/roundtrip_bucket_report.py`, which
decomposes every round trip into `byte_exact` / `key_equal` / `facmer_divergent` / `structural` /
`hard_fail` / `encode_fail`. Two levers earned promotion; three were measured and kept opt-in.

**Accuracy: real gain.** On the worst-cohort A/B, `OIN_EARLY_EXIT` — accept the first conformer
that independently re-encodes to the key, instead of exhausting the pool — took byte-exact
**44.7% → 60.5%** and key-match **55.3% → 73.7%**, with zero regressions. OIN-direct assembly
(build the internal `MetalComplex` straight from the parsed OIN, retiring the winding-lossy
m-SMILES bridge to a fallback) was proven **accuracy-neutral twice** and promoted for what it
unblocks: metal `@SPn` chirality now survives into 3D generation. A separate 3917-molecule
regression sweep measured **11 regressions / 1092 fixes**, and all 11 regressions are 300 s
timeouts — **zero correctness regressions**.

**Speed: real gain.** `OIN_EARLY_EXIT` is ~**5×** faster on the worst cohort. Encoder-side, the
`AC2BO` valence-combination sort was capped and `ResonanceMolSupplier` enumeration moved into a
forked, `RLIMIT_CPU`-bounded child, recovering previously-unencodable molecules (`BENVOG`,
`HUCNAU`) without changing any currently-encodable OIN.

**Refuted.** Rigid η-winding construction, difficulty-ordered Kabsch greedy placement, and a
stretched-bond acceptance metric were all measured net-negative or neutral and kept OFF — a third
independent confirmation that **selection beats construction** in this pipeline.

**Cohort caution.** The `byte_exact 81.19%` figure often associated with v0.4.4 is the *v0.4.4 key
applied to v0.4.2-generated data* (6719 molecules, `results-v0.4.4-sl4/`). It measures the new
classifier, not new generation.

Suite: **488 → 551 OK**.

---

### v0.4.5 — the emitted string becomes canonical

**What shipped.** Six canonicality levers promoted to default-ON *together*, via a new
single-source registry (`src/oinsmiles/oin/levers.py`): `OIN_CANONICAL_BODY`,
`OIN_CANONICAL_PERCEPTION`, `OIN_CANONICAL_SLOTS`, `OIN_CANONICAL_ETA_WINDING`,
`OIN_STABLE_METAL_AC`, `OIN_STABLE_STEREO`. What made them safe to promote as a block: each
**repairs a renumbered presentation without rewriting the canonical answer**, which is why the
corpus shows no churn.

**Accuracy: real gain, and the important number is not `byte_exact`.** On a generator-free probe
(300 molecules, seed 42, graph held fixed, varying only proper rotation and atom renumbering — so
the correct answer is byte-identical by construction):

| arm | byte-stable | comparison **key** broken |
|---|---|---|
| all levers OFF | 173/298 — 58.1% | 60 — 20.1% |
| all levers ON | 208/299 — **69.6%** | **16 — 5.4%** |

The key-instability figure matters most: the comparison key *is* the harness's acceptance
predicate and the basis of every accuracy number the project reports. It went from unstable on
**1 molecule in 5** to **1 in 19**. On the 936-molecule re-baseline, **145 of 436** previously
failing molecules were fixed (33.3%), with all 11 apparent regressions being 300 s timeouts
against an 1800 s baseline ⇒ **zero correctness regressions**.

**Speed: real gain, unconditional.** The perf lane found `OIN_EARLY_EXIT`'s confirming step
(`_reencode_oin`, a full `XYZToSMILES().convert()` measured at **48–57 s per call** on a 59-atom
eta molecule) running **twice** per conformer — once in `accept_fn` during pool fill, then again,
redundantly, in `_select_by_geometry_impl`'s own re-scan over the same pool. Anything surviving to
`_select_by_geometry` is provably already a confirmed non-match, so the second test recomputed a
known answer. Memoised, no lever. The encspeed lane took the project's **first ever encoder
profile**: `perception_core.AC2BO` is **99.8%** of a slow encode — not "dominant", essentially all
of it — and it is neither an eta nor a size phenomenon (a 39-atom non-eta molecule profiles the
same and is *slower*).

**Refuted / corrected.** The **inverted CIP goldens** for `PdCl2-RR-BDPP` and `PdCl2-RR-BDNN` had
been wrong for four months: the test that "verified" them ran `rdCIPLabeler` on a SMILES reparsed
from the encoder's *own* output, and `rdCIPLabeler` converts a parity tag into a label without
checking it — so an inverted tag was self-consistent and passed. Ground truth from
`AssignStereochemistryFrom3D` is (R,R). Separately, an unplanned Lane 8 found **13% of molecules
emitted different absolute stereochemistry under pure atom renumbering**; no prior instrument
could see it.

**Known limitation shipped honestly.** P3 (metal-bound 2° amine) is unusable in the shipped
default — `OIN_CANONICAL_BODY`'s reparse clears the `[N@]` it stamps, and the obvious fix is
measurably wrong.

Suite: **837 OK**.

---

### v0.4.6 — boron, Δ/Λ, and the first crack in the metric

**What shipped.** `OIN_BORON_CAGE` promoted to default-ON; a complete metal-centred Δ/Λ descriptor
pipeline (`oin/metal_config.py`) behind `OIN_EMIT_METAL_CONFIG` (default OFF); and — decisively for
everything that follows — `OIN_INDEP_SCORE` **added** (default OFF), recording the honest round-trip
verdict beside the scored one.

**Accuracy: real, small, and correctly priced.** On the 936-molecule re-baseline, 34 of the 36
`XYZToSMILES failed` molecules are electron-deficient boron clusters, and the lever takes that
population from **0/36 encoding to 34/36** at 0.2–4.2 s each. ⚠ It also moves **14 molecules from
scored-passing to failing** — which is *correct*: they passed while describing the wrong graph
(`VEJXOZ` invents a C=B double bond). It trades 14 silent false positives for 14 loud honest
failures, so a headline pass rate can move either way.

**Speed: a measured regression, accepted.** Promoting the boron lever moved most of that class
from failing **instantly** to failing **slowly**. Of 33 measured, only **2 generate** a 3D
structure, **25 burn the entire generation cap** producing nothing, 6 fail instantly — roughly
**2.1 CPU-h wasted per 5000-molecule sweep**. `XIQKOY_comp_0` is the two-point proof: lever OFF it
fails in 0.87 s; lever ON it encodes a correct B10 cage and then runs past 340 s. Still the right
call — a right-graph loud failure beats a wrong-graph silent pass — but the cost is now on record.

**Refuted — eight hypotheses, each with the measurement that killed it.** Recorded in
`V046_HFAITHFUL_FINDINGS.md`: P3 tag restoration through the reparse; "the accuracy gap is mostly
compute" (75% of timeouts hide real failures, not latent passes); a donor-cut rule for hydrogen;
eta incremental pool widening; and three successive Δ/Λ formulations.

**The finding that set up v0.4.8.** The 5000-molecule seed-42 sweep landed here:
`byte_exact` **4140/5000 = 82.80%**, round-trip 93.2%, median 7.2 s, **994 molecules over 30 s**.
Alongside it, the false-positive measurement: **61/633 = 9.6%** overall, **28.1% of haptic**
molecules, plus 8/302 false negatives the other way — one root cause, net **~5.7 points
over-stated**. And `OIN_ACCEPT_SCORED`, which "recovered" 90 of 340 failures, was held OFF because
`report["coordination"]` showed those 90 to be **60 DEGRADED / 21 BOUNDARY / 9 INTACT** — a claimed
+90 is really **+9**.

Suite: **838 OK**.

---

### v0.4.7 — a release of measured negatives

**What shipped.** Five swimlanes. **Four refuted the thing they were built to ship, and the fifth
refuted its own sibling.** Every lever default-OFF.

**Accuracy: zero, and the evidence that it is zero is the deliverable.** The
`results-v0.4.6-sweep` bucket report is **byte-for-byte unchanged** — `byte_exact` 4140 / 82.80%,
`> 30 s` 994 / 19.88%.

**Speed: zero on the default path.** The encode floor was characterised, not reduced.

**Refuted.**
- **`OIN_ACCEPT_SCORED`: DO NOT PROMOTE**, superseding this lever's own earlier "promote" reading.
  The gate that recommended it was **circular** — `passed` is computed with the same predicate the
  lever accepts on, so "18/22 both arms, zero regressions" could not detect what dropping the
  strict step costs. Measured with a genuinely independent arm: **15/20 → 7/20, 8 regressions,
  0 fixes**, one-way. Six of the eight lose haptic coordination outright, with the metal geometry
  tag degrading in lockstep (`[Ru_TET]`→`[Ru_TPL]`). **Byte-identical notation, changed geometry** —
  the cost is invisible to the very metric that would police it.
- **`OIN_BORON_GEN_FASTFAIL`'s discriminator was refuted before it could do harm.** `ULODUU_comp_0`
  is `[Zr_TET]` and *does* generate, in 61.8 s — but the lane's sweep capped at 30 s, so a 61.8 s
  success read there as a cap-burner. **The class boundary moves with the compute you give it.**
  Cross-tabulated, every geometry with a success also has failures, so geometry separates nothing —
  and it was the last discriminator standing after hapticity, size and denticity.
- **The encode floor is three cost regimes, not one.** `XIRMER`'s 20-minute encode is **3 forked
  resonance timeouts at 95.8% of the fork budget each** — the CPU-limit backstop firing, not slow
  arithmetic, so tuning the solver would never have touched it.

**What it built that mattered later.** `OIN_ATTACH_CHECK` (7 of 8 known regressions separated,
0 false positives, 7–81 ms against the 48–57 s strict test it replaces) — and it **never reads a
bond object** as evidence of attachment, because a detached ligand keeps its bond. Plus the
**two-arm byte-identity gate** (`tools/gate_v047.sh`), the first instrument that could tell a
notation change from a compute change.

Suite: **909 OK**.

---

### v0.4.8 — the honest number

**What shipped.** `OIN_INDEP_SCORE` promoted to default-ON. **No encoder, generator or notation
behaviour changed.** The measurement changed.

**Accuracy: the reported figure falls 10.34 points, on purpose.** Same 5000 molecules, same
conformers, same key, same `status` gate — the single variable is which round-trip string the
verdict reads:

| bucket | scored | % | honest | % | delta |
|---|---:|---:|---:|---:|---:|
| `byte_exact` | 4140 | 82.80% | **3623** | **72.46%** | **−517** |
| `key_equal` | 520 | 10.40% | 610 | 12.20% | +90 |
| `facmer_divergent` | 1 | 0.02% | 16 | 0.32% | +15 |
| `structural` | 9 | 0.18% | **417** | **8.34%** | **+408** |
| `hard_fail` | 315 | 6.30% | 319 | 6.38% | +4 |
| `encode_fail` | 15 | 0.30% | 15 | 0.30% | ±0 |

643 molecules moved bucket: **613 degraded, 30 improved** — the correction runs in both directions,
which is what distinguishes it from a pessimism knob. At corpus scale **36.7% of haptic
`byte_exact` passes were false**, against 6.7% non-haptic.

**Speed: unchanged, and stated as such.** Nothing was re-generated, so `metrics.elapsed_s` is
untouched and `> 30 s` stays 994 / 19.88%.

**Why the drop is trustworthy.** An encoder that mis-perceived generated geometries would
manufacture exactly this result, so the honest arm was cross-checked against `report["coordination"]`,
which reads distances only and consults no bond graph:

- **control** (3595 `byte_exact` in both arms): **1.3% flagged** — inside its 3.7% false-alarm band
- **moved** (428 pass → `structural`/`facmer_divergent`): **64.7% flagged — a 50× enrichment**
- and the two agree on **mechanism**, not just count: `contacts_lost ÷ sites_lost` lands on the
  integers 1–7, and 95.2% of the ratio-above-1 population carries a haptic token — because an η⁵-Cp
  is five *contacts* but one coordination *site*.

**Refuted.** The recorded suspicion that a ~19× `structural` inflation on re-encode was an artifact
of the method. It is real, and two thirds of it is independently confirmed. Also refuted: the cost
argument that had held the lever off. It priced the second encode at "0.4–1.5 s/molecule"; the
measured figure is **0.33 s/molecule** and the whole corpus re-scores in **334 s** against the
~55 CPU-h a live re-sweep costs. **There was never a cost case for scoring dishonestly.**

**A second gate found load-bearing.** The `Atom count mismatch` check is not a third error
direction — it is the only instrument that catches 18 molecules at corpus scale, hydrogen-only in
18/18, of which **8 re-encode byte-identically to their input**. No string comparison of any kind
can separate a structure carrying two extra hydrogens from the original. `XAKCAP_comp_0` defeats
four instruments at once: scored string PASS, honest string PASS (the same string), key EQUAL,
`coordination` INTACT, atom count 61 ≠ 63.

**Scope, deliberately narrow.** This changed what is *reported*, not what is *accepted*.
`accept_fn` is untouched. Scoring the acceptance ladder honestly would move runtime *and* the
failure mix in the same release that re-baselines the number, making both unmeasurable.

Suite: **920 OK**.

---

### v0.4.9 — speed becomes measurable, and the release refutes its own charter

**What shipped.** `OIN_ENFORCE_BUDGET` (default OFF), `BudgetExhaustedError`, and a **frozen,
stratified 328-molecule runtime benchmark** reproducing to **0.28%** (277.01 s vs 277.79 s,
byte-identical rows) — the noise floor every later runtime claim must clear.

**Accuracy: zero.** Confirmed by a live 300-molecule run: **296/300** per-molecule agreement with
the offline re-score, **0/300** encoder drift.

**Speed: nothing promoted, but the goal became expressible for the first time.**
`OIN3DGenerator(timeout=)` was confirmed **advisory** — 60 s requested, 60.0–172.8 s spent, worst
case **2.9×** — because the deadline is checked only at the top of the embed attempt loop and an
in-flight attempt always runs to completion.

**Refuted — its own chartered justification.** The charter rested on one number: *"759.9 s against
a 300 s budget."* That is arithmetic on a sum. Split by `tier_passed`, all **4658** single-attempt
rows finish within **0.2 s** of their 300 s cap. The advisory-timeout defect is real and two direct
probes measure it; **the corpus number was simply never evidence for it.**

**Refuted — the first design of its own lever.** Threading the deadline into `get_embedding` alone
— the function two profiles had indicted — measured **ε = +48.4 s on a 30 s budget: it changed
almost nothing.** There is **no single sink**: `FOSNEI` is 75% `get_embedding`; `CAHQEJ` is 77%
`get_embedding` + `numpy.linalg.eig`; `VAFMIA` is **99% `chirality._reparse_cip_label_once`**. A
bound threaded into whichever function profiled expensive last is not a bound.

**Refuted — the inherited cohort.** The v0.4.7 gate's 100-molecule slow cohort is **59% false
positive** when cross-tabulated against the honest baseline (59 of 100 are `byte→FAIL`), because it
selected on the *scored* verdict — and a slow cohort is made of exactly the haptic molecules the
scored verdict over-states most.

**Measured, not delivered.** With the bound on, **ε = +32.8 s** on a 30 s budget and the **same 11
molecules** exceed it in every arm. The bound compresses the tail (2.63× → 2.09×); it does not
remove it. `max(elapsed_s) < 30 s` is **not delivered**.

**The finding that reordered the roadmap.** `structural` was 9 molecules / 0.18 points when scored
dishonestly. Honestly it is **417 / 8.34 points — the second-largest block in the gap**, ahead of
`hard_fail`. Also found and deliberately left for the next release: `get_embedding`'s discarded
`.index()` call, worth 22% of an eta generation.

**Fixed, quietly important.** `--mol-timeout` never reached the generator (the harness hardcoded
it, so a bound could not be A/B-ed through the harness at all), and the v0.4.7 gate silently
resolved to an arbitrary interpreter — once an unrelated project's venv with rdkit 2025.09.2
against the pinned 2025.9.3. **A byte-identity gate on a different rdkit reports MISMATCHes that
read as code regressions.**

Suite: **930 OK**.

---

### v0.4.10 — cost per attempt, and the arbiter was broken first

**What shipped.** Two performance changes and a gate repair.

**The first thing measured was the gate itself, and the gate was dead.** `tools/gate_v047.sh arm1`
had been exiting 1 **before comparing anything** since `dd51a515`: `ULODUU_comp_0.xyz` was added as
a fixture and the golden was never extended, so `EXPECTED_FIXTURE_COUNT = 61` hard-refused every
run. **v0.4.9 froze a 328-molecule benchmark and measured a 0.28% noise floor while the encoder arm
of its own gate was non-runnable.** Behind the refusal sat a second, real drift: the v0.4.7
`xyz2mol` → `perception_tmc` rename was behaviour-neutral but not *string*-neutral, and ARM 1
hashes error strings on purpose. The other **60 rows were byte-identical**, so the encoder itself
had not moved; both rows re-frozen at 62.

> **A gate that fails before it compares is indistinguishable from a gate that is merely
> inconvenient to run, and it silently stops covering everything else it was watching.**

**Accuracy: zero, by construction.** ARM 1 **62/62** and ARM 2 **90/90** byte-identical, for both
lanes, on all four gate runs, plus identical generated-structure fingerprints across every A/B.

**Speed: one default-path win, one shelved.**

| molecule | class | `.index()` deleted (default-ON) | CIP memo (default-OFF) |
|---|---|---|---|
| `CAHQEJ_comp_0` | eta, `[Ni_TPL]`, 2 haptic | **−32.9%** | −2.4% |
| `VAFMIA_comp_0` | `[Cu_LIN]`, adamantyl NHC | — | **−86.7%** (81.89 s → 10.87 s) |
| `FOSNEI_comp_0` | non-eta, boron cage | **+0.3% (nil)** | — |

`get_embedding`'s outer loop called `alternative_ace_mol_list.index(alternative_ace_mol)` and
**threw the result away** — not a no-op, because `list.index` compares via `Molecule.__eq__` →
`is_same_molecule` → `get_c_eig_list` → `numpy.linalg.eig`. It eigendecomposed a Coulomb matrix per
candidate per outer iteration. Ships **on by default with no lever**: gating provably dead code
would permanently ship `if not lever: <discarded computation>`.

**The bimodality is attributed, not merely observed.** `CAHQEJ` makes **99** such comparisons
costing **38.52 s (38.8% of generation)**; `FOSNEI` makes **3**, costing 0.03 s. The `FOSNEI` null
is explained, not excused as noise — and the two independent measurements agree on size (38.8%
attribution vs 32.9% quiet A/B).

**Corrected before release.** A first A/B taken while four gate processes competed for the box
(load average **35**) reported −50.2% and +9.6%. Quiet-box re-measurement gives **−32.9%** and
**+0.3%** — the gain was over-stated by 17 points and the null was buried in 30% within-arm spread.
**Byte-identity gates are load-immune and can be parallelised freely; wall-clock is neither.**

**Refuted.** The chartered Lane 1 (SVD in `_finalize_positions`) — its premise does not reproduce:
all 165 `update_xyz` calls across 3 molecules converge in **exactly 1 iteration** against a cap of
30, and the charter's "1189 SVD calls / 15.6 s on HIDCIH" measures 74 calls / 1.75 s. Chartered
Lane 2 (per-attempt redundancy) was already fixed by an earlier release. **Executing the charter as
written would have spent the release on its third-best target.**

**Not done, deliberately.** `OIN_MEMO_CIP_REPARSE` is **not promoted**. The charter permits
same-release promotion only if byte-identity holds on the whole benchmark; this release ran the
**fast band — 90 of 328**. Fast-band-only evidence is not that. So the single largest measured
speed win in the project's history ships **off**.

Suite: **946 OK**.

---

### v0.4.11 — the fix that worked, and had to be refused

**What shipped.** `OIN_CANONICAL_DONOR_FOLD`, **default OFF**, plus
`tools/mirror_audit_donor_fold.py` — the instrument that caught the problem.

**Accuracy: FLAT — and this is the most instructive release in the set.** The lane built the fix
v0.4.5 Lane 2 had specified in writing for the **largest single block in the gap** (496 molecules /
9.92 points). It works, exactly as designed:

- **393 molecules** move `key_equal/slot_renumber → byte_exact`, **none in any other direction**
- `facmer_divergent` holds at 16
- the comparison key moves on **0 of 992** strings
- both gate arms stay byte-identical with the lever off
- **+7.86 `byte_exact` points**

**It also collapses enantiomers.** A corpus mirror audit on a *uniform* 250-molecule draw found
**19 (7.6%)** whose mirror image encodes identically once the fold is on. Run directly on the 393
molecules the fold claims as wins: **221 of 393 (56.2%)** collapse. An independent geometric oracle
(`tools/injectivity/oracle.py`, sharing no machinery with the encoder) confirms **18 of 19**
uniform-draw collapses and **26 of a 30-sample** of the 221 are genuinely chiral. **More than half
the gain is the damage.**

**Why the safety argument failed.** It assumed two donors in one `breakTies=False` symmetry class
are interchangeable. `CanonicalRankAtoms` computes the symmetry of the **isolated ligand graph**,
but those donors sit at distinct vertices whose relation to the *other* ligands is chirality-bearing
— so the vertex permutation their exchange induces can be **improper**. *A fragment's automorphism
says nothing about the parity of the vertex permutation it induces.* v0.4.5's restriction to proper
rotations was not conservatism; it was the load-bearing correctness condition.

> ### 🔴 The finding that outranks the points
> **`byte_exact` can be raised by deleting information, and the comparison key will agree.** Both
> are blind to reflection, because `_parse_vertex_colors` folds that axis deliberately. **A
> one-directional transition matrix is not evidence of safety for anything touching
> canonicalization.** Mirror-audit every future canonicality lever, on a uniform draw, before
> quoting its points.

**Speed: no change claimed.**

**Two cohorts, not one.** Uniform 19/250 (7.6%) and the runtime-stratified cohort 31/300 (10.3%) —
the damage is not an artifact of one sample. The named Δ/Λ fixtures (`ZUMNEC`, `fac-Ir(ppy)₃`)
**both pass with the lever on**: neither carries the vulnerable motif, and **fixtures alone could
never have caught this.**

**A methodology error recorded rather than buried.** The stratified audit was read mid-run as "0
regressions in the first 200" and briefly written up as a *"wrong stratum"* finding. The tool prints
a verdict only every 50th molecule, so four clean progress lines had been mistaken for 200 clean
molecules; the completed run reports 31 collapses. **A partial run is not a result, and a progress
line is not a tally.**

**Also measured.** All 496 `slot_renumber` molecules classified for the first time:
`same_vcolor_identical` **496/496**; `diff_occupancy`, `diff_geometry`, `diff_colors` and
`postpass_BUG_diverges` all **0** — the charter's headline risk refuted in the release's favour.
And **90 of the 118 `distinct_donors_LOCAL` are frozen-resonance-form artifacts** (acac binds
through two equivalent oxygens but is written ketone/enol), which is a ligand-*body* gap ⇒ v0.4.14
is re-sized from 2.28 to **~4.08 points**.

Suite: **957 OK**.

---

### v0.4.12 — the parity filter

**Status: RELEASED on `main` at `675e5425`, `pyproject` 0.4.12 — but no `v0.4.12` tag exists yet.**

**What shipped.** `OIN_FOLD_PARITY_VETO` (default OFF) — declines the donor fold on any
molecule where folding would make the structure's mirror encode identically. And
`OIN_ETA_ACCEPT_EXIT` (default OFF) — the eta winding criterion relocated from
`_select_by_geometry_impl` into `accept_fn`, the only site consulted *during* pool filling.

**Accuracy: FLAT, deliberately — both levers ship OFF.** But the filter works: a uniform
250-molecule mirror audit goes **19 → 0** `REGRESSION_raw_collapsed`, with those 19 moving into
`distinct_both_arms` (73 → 92) and `achiral_or_preexisting_fold` **unmoved at 157** — the direct
evidence that the veto did not simply refuse to fold everything. **The surviving gain is +3.42
points** (171 of 393). v0.4.11 bounded the safe set at *"at most ~172 (~3.44 pts)"* by counting
collapses; this filter counts survivors through the shipped predicate and lands on **171**. Two
independent routes agreeing to one molecule.

**And it holds at a second, disjoint seed — which is the part that matters.** v0.4.11's recorded
methodology lesson was that one draw cannot separate a corpus property from a sampling artifact.
Seed 11 answers it: the **baseline independently reproduces 19/250** on a disjoint sample, so the
**7.6% collapse rate is a property of the corpus**, not of seed 7. The veto zeroes both seeds while
`achiral_or_preexisting_fold` stays untouched at **157 and 141** respectively — so it is separating
the right set rather than declining to fold. *That the achiral population is unmoved is the
second-most-important number here: a veto without that left conjunct would have swallowed all 157
as well, still reported `REGRESSION_raw_collapsed = 0`, and looked like a total success while
destroying the fold's entire benefit.* Lane 1's gate is now fully met; the ~55 CPU-h re-sweep and
the re-frozen goldens remain as v0.4.13 preconditions.

**Speed: no runtime claim in the release.** `OIN_ETA_ACCEPT_EXIT`'s timing A/B was **stopped rather
than banked** at load 26 — *never interleave timing runs with gate runs* — so its predicted tail
reduction shipped **unverified**.

> ⚠ **Post-release, not part of v0.4.12.** A lane commit on `swimlane/v0412-eta-accept`
> (`223c7244`, unmerged at the time of writing) re-ran that A/B on a quiet box against a
> **re-derived** 12-molecule population: median **61.8 s → 18.82 s (−69.5%)**, total −47.4%,
> `>30 s` **11 → 5**, with the same bimodal signature as v0.4.10 (`RIRYOJ` 119.04 s → 4.4 s, a 27×
> win with a byte-identical string; two molecules nil). **The runtime thesis is confirmed — and the
> lever still stays off.** A fifth gate arm, the metal-configuration check v0.4.11 paid for, caught
> `KIHHUG_comp_0` returning a *different Δ/Λ descriptor* with a **byte-identical emitted string**;
> four of the five gates call that molecule untouched. The charter makes any Δ/Λ divergence
> blocking. *A large real speedup that costs something four of five instruments cannot see is
> precisely the trade this release exists to detect.*

**Refuted before it was built.** Reflection parity is **not a property of the emitted string**: a
donor swap is a transposition fixing every other vertex, so the obvious `det > 0` test on the
polyhedron rejects *every* swap and degenerates the fold to the identity. The veto had to run on
the pristine conformer in `get_oin_string` instead.

**Two silent defects caught only by disbelieving a clean result.** `tmc_mol`'s atom order is **not**
the coordinate order (`__origIdx`) — zipping positionally encoded `BIWDIV` as `[Co_TBP]` with
invented bonds. And the mirror was encoded with the fold *inherited*, which disarmed the achiral
guard and made the self-check decline on **18 of 18** movers **while all three fixture tests
passed** — because declining to fold also separates a mirror pair.

> **Ask what a broken version of the instrument would print. If it is the same thing, you have
> measured nothing.**

**Also corrected.** `OIN_ETA_EARLY_EXIT`'s promotion gate is **void, not unrun** — it runs
downstream of a fully-filled pool and cannot reduce embed count by construction (its own in-code
A/B already recorded Ferrocene at 32 attempts both arms). The chartered v0.4.6 accept-gap cohort is
**stale on 8 of 8**; the real eta target population is **405 molecules whose key never matches, 378
of them > 30 s**. And the carry-forward licence is now **measured**: re-encoding the frozen corpus's
inputs *and* stored generated structures reproduces the v0.4.8 strings on all 393 movers (0 drift),
so v0.4.9–v0.4.11's byte-identity claims are confirmed rather than inherited.

Suite: **988 OK**.

---

### v0.4.13 — the fold ships, and the release cancels its own 55 CPU-h precondition

**Status: RELEASED, tag `v0.4.13`, `main` @ `9fbdf678`, `pyproject` 0.4.13.**

**What shipped.** **PROMOTED to default-ON, together:** `OIN_CANONICAL_DONOR_FOLD` and
`OIN_FOLD_PARITY_VETO`. They are one promotion with two names — the fold is the change that pays,
the veto is the condition under which it is safe, and
`test_levers::TestDonorFoldAndParityVetoAreCoupled` pins them so neither can be demoted alone.
**Added but default-OFF:** `OIN_PREFILTER_ADVISORY`, which makes the cheap acceptance prefilter
advisory instead of dispositive.

**Accuracy: REAL GAIN — `byte_exact` 3623 → 3794, 72.46% → 75.88%, +3.42 points.** 171 molecules
move `key_equal/slot_renumber → byte_exact`; **0** move in a bad direction; `facmer_divergent`
unchanged at 16. This is the first upward movement in the window this report covers.

**The safety gate ran on two independent draws, and the accounting is what carries it.**

| | cat/ (v0.4.12's draw) | cat+photo (mixed) |
|---|---:|---:|
| `REGRESSION_raw_collapsed`, veto OFF → promoted | **19 → 0** | **33 → 0** |
| `achiral_or_preexisting_fold` | 157 → 157 | 134 → 134 |
| `distinct_both_arms` | 73 → 92 | 83 → 116 |
| accounting | **73 + 19 = 92** ✓ | **83 + 33 = 116** ✓ |

A veto that works and a veto that declines everything both report zero collapses — v0.4.12's own
first implementation declined on **18 of 18** while all three fixture tests passed. Here every
collapsing molecule lands in `distinct_both_arms` and **nothing else moves**, so the zero is bought
by *separating*, not by *abstaining*. The cat/ arm reproduces v0.4.12's published table
(157 / 73 / 1) line for line; the mixed draw reads 13.2% against cat-only's 7.6%, which is a second
population rather than a contradiction.

**Refuted: the 55 CPU-h re-sweep this release was chartered to run.** v0.4.12 made it a promotion
precondition. It was not run, and running it would have produced **weaker** evidence:
`results-v0.4.8-honest` — the 72.46% baseline itself — is an offline **re-score**, not a sweep, and
**no generator sweep has run since v0.4.6**. Re-running a *stochastic* generator to A/B an
*encoder-side* change contaminates the signal with run-to-run variation and is not even
like-for-like with the table it is compared against.

What licensed the offline route was measured rather than assumed. `accept_fn` accepts by comparing
the round-trip **key**, so the question is whether the fold ever moves one:

| `tools/fold_key_invariance.py`, whole corpus | |
|---|---:|
| strings compared | **9669** |
| strings the fold MOVED | **1019** |
| strings whose KEY changed | **0** |

The 1019 is the load-bearing half — *a lever that never fired would also print 0 key changes.*
Verdict `GENERATOR_NEUTRAL`, reproducing v0.4.11's "0 of 992" at ~10× the scale.

**Two n = 1 classes finally sized, and one re-mechanises a release four rungs out.** Over the 767
genuine failures: **GAVSED (`DETACHED`) 280**, **MEDZUR (`INTACT`) 99**, `BOUNDARY` 53,
`NO_STRUCTURE` 335 — against a `byte_exact` control of **1.32%** `DETACHED` vs 24.11% on the
failing side (**18.2× enrichment**, `UNKNOWN` = 0). **`structural` is 266/417 = 63.8% `DETACHED`**:
those molecules did not fail because the generator could not assemble them — it assembled something
and *returned* it with ligands off the metal, because `_select_by_geometry`'s fallback ranking is
not attachment-aware. A one-site **return-path guard**, not a capability floor. Independently
cross-checked: `missed_success_audit` partitions the same 767 by cause and its largest bucket is
**437 molecules it labels "ambiguous — generator OR notation"**; this split resolves **280 of them
to generator**.

**Four instruments caught printing plausible nothing.** The mirror audit died silently on a dataset
a branch switch had deleted (26,232 files; `git status` stays clean). `fold_transition_sim.py`
excluded **all 393 movers** — its `--dataset` default is relative and the run was from a worktree —
and printed **+7.86**, the *refuted* bare-fold figure, under a heading saying "veto", exit 0.
`run_sweep.sh` still carried the sibling-glob venv trap `gate_v047.sh` was hardened against in
v0.4.9, live in the one place it costs 55 CPU-h. And **ARM 1 passed byte-identically before *and*
after the promotion because 0 of its 62 fixtures were fold-movers** — a PASS that means "no
regression", not "the change works". Each now refuses or reports its denominator.

> **A gate that cannot see your change is not evidence about your change. State its coverage of the
> moved population, not its verdict.**

**What it cost.** ARM 2's v0.4.9 golden had **11 of 325 rows** genuinely move; they were re-run
individually and patched rather than bulk-accepting a regenerated manifest, and `MANIFEST_SHA256`
was recomputed (arm2 does not check it, so a stale one would have gone unseen). Lane 1's corpus
prevalence is **n = 1 and therefore unquotable** — the lever, telemetry and harness are built and
the defect is confirmed on `AROHIA_comp_0`, but the measurement needs a quiet machine and is handed
to v0.4.14. The 48 `byte_exact` molecules that read `DETACHED` remain unexplained. One merge commit
shipped without a trailer (`--no-edit`; the hook only rewrites, never inserts).

Suite: **993 OK**.

---

## What is not known

Stated explicitly, because the gaps are as decision-relevant as the numbers.

1. **Corpus speed since v0.4.6 is unmeasured.** `994/5000 = 19.88% over 30 s, median 7.19 s` is a
   v0.4.6-era figure. v0.4.10's default-ON deletion has never been run at corpus scale. The cheap
   way to close this is the frozen 328-molecule v0.4.9 benchmark (~1–2 CPU-h), not a 55 CPU-h
   sweep.
2. **Corpus accuracy still rests on the v0.4.6 generator run — and v0.4.13 argues that is correct,
   not a gap.** No generator sweep has run since v0.4.6; 75.88% is an offline re-score of those
   stored structures, exactly as 72.46% was. The re-score is **exact** for an encoder-side change
   that moves no comparison key, and v0.4.13 measured that condition (0 of 9669). What remains
   genuinely unmeasured is anything that would change what the *generator returns* — so the next
   release that touches acceptance or selection owes a real sweep, and the key-invariance check is
   how it finds out. Note the drift control passed on all 393 movers, which turns three releases'
   byte-identity *claims* into measurements.
3. **`OIN_MEMO_CIP_REPARSE`'s promotion gate has not been run** — the full 328-molecule cohort with
   the lever on, ~10 CPU-h sharded 6-way. Until it is, the project's largest measured single speed
   win stays off.
4. **The four bucket-report cohorts are not comparable**, so there is no honest like-for-like
   accuracy line from v0.4.4 to today. Building one requires re-running the seed-42 5000-molecule
   corpus on current `main`.
5. **`OIN_ETA_ACCEPT_EXIT` has no runtime claim at all** — its A/B was stopped, not banked.
6. **`PREFILTER_VETO`'s corpus prevalence is n = 1.** v0.4.13 built the lever, the telemetry and
   the harness, and confirmed the defect on `AROHIA_comp_0` — 2 conformers the cheap prefilter
   rejects that the strict test accepts. Whether that is 2 molecules or 200 corpus-wide is
   **unmeasured**, and the release declined to quote n = 1 as a prevalence.
7. **The 48 `byte_exact` molecules that read `DETACHED` are unexplained.** Either the notation does
   not express the lost metal contact, or `coordination.intact` is over-sensitive at that
   tolerance — the 1593-molecule `BOUNDARY` band says the second is possible. Neither was tested.
8. **The MEDZUR class — 99 molecules with attachment intact and re-perception still disagreeing —
   has no mechanism.** It is now sized but no better understood than when it was n = 1.

---

## Where the remaining 24.12 points are

From `docs/agentic-notes/ROADMAP_100_100.md`, honest baseline, **re-derived at v0.4.13**. The
previous copy of this table read 27.54 and is superseded: the promotion consumed 171 molecules.

| block | n | pts | nature | owning release |
|---|---:|---:|---|---|
| **`structural`** | **417** | **8.34** | 🔴 **re-mechanised: 266 (63.8%) are `DETACHED`** — an unguarded *return* path, not a capability floor | v0.4.17 — **candidate to pull forward; ~5.32 pts is one site** |
| `hard_fail` | 319 | 6.38 | compute; **315/319 produce no structure at all** | — |
| `key_equal` → `slot_renumber` | **325** | **6.50** | canonicality, encoder-side. **Was 496; the fold took 171.** 90 of the residue are frozen resonance forms (ligand **body**) | v0.4.14 |
| `key_equal` → `rdkit_canonical` | 114 | 2.28 | canonicality, encoder-side | v0.4.14 (**⇒ ~4.08 with the 90**) |
| `facmer_divergent` | 16 | 0.32 | wrong isomer | v0.4.15 |
| `encode_fail` | 15 | 0.30 | encoder coverage | v0.4.18 |
| **sum** | **1206** | **24.12** | | |

**The failure side also has a second decomposition now — by mechanism rather than by bucket.**
Over the same 767 genuine failures: **280 `DETACHED`** (the generator returned a structure with
ligands off the metal), **99 `INTACT`** (attachment fine, re-perception still disagrees), 53
`BOUNDARY`, 335 `NO_STRUCTURE`. Two independently-written tools partition that set and agree to
within ~5 molecules.

⚠ **The open sequencing question is now sharper, not softer.** `structural`'s 266 `DETACHED`
molecules are worth up to **5.32 points at one code site**, and they sit at v0.4.17 — behind
v0.4.15, which is explicitly chartered as *"knowledge, not points"*. Whether that should be
reordered is recorded as an **open decision for the project owner** in
`LADDER DECISION 2026-07-27 (v0.4.13)`; it was deliberately not applied by the session that found
it, because the last such reorder was the owner's call.

**The other structural fact:** the two goals are one goal. Of the 340 failures in the v0.4.6 sweep,
**78.8% never test the notation** — 240 are generator timeouts and 28 produced nothing. The
notation-attributable gap is **57/5000 ≈ 1.1%**: *where the notation is actually exercised, it is
~98.9% correct.* A per-molecule 30 s cap recovers **37.84 CPU-h** and costs **251 passes
(5.02 points)**; 93.1% of honest passes already finish under 30 s.

---

## Sources

| what | where |
|---|---|
| Per-release narrative and figures | `CHANGELOG.md` §§ [0.4.4]–[0.4.13] |
| **v0.4.13's frozen numbers, tracked and public** | **`measurements/v0.4.13-honest/`** (13 files: both mirror draws, the transition record naming all 171, the generator-neutrality proof, the attachment split) |
| The promotion, the cancelled sweep, the four blind instruments | `docs/agentic-notes/v0.4.13/PROMOTION_AND_CLASSES_v0.4.13.md` |
| MEDZUR / GAVSED sizing + the `missed_success_audit` cross-check | `docs/agentic-notes/v0.4.13/LANE-02-attach-classes.md` |
| `PREFILTER_VETO` — confirmed defect, unmeasured prevalence | `docs/agentic-notes/v0.4.13/LANE-01-prefilter-advisory.md` |
| Honest re-baseline + transition matrix | `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/bucket_report_both.md` |
| Cohort bucket reports | `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.{4-sl4,4-regression,5-rebaseline,6-sweep}/` |
| v0.4.5 promotion evidence | `docs/agentic-notes/v0.4.5/PROMOTION_GATE_v0.4.5.md`, `ENCODER_PERF_v0.4.5.md`, `PERF_v0.4.5.md` |
| Speed A/Bs, bimodality attribution | `docs/agentic-notes/v0.4.10/SPEED_v0.4.10.md` |
| Cost regimes, benchmark provenance | `docs/agentic-notes/v0.4.9/BUDGET_BOUND_v0.4.9.md`, `RUNTIME_BENCHMARK_v0.4.9.md` |
| The `elapsed_s`-is-a-sum refutation | `docs/agentic-notes/v0.4.9/ELAPSED_S_IS_A_SUM_v0.4.9.md` |
| v0.4.8 honest baseline + atom-count gate | `docs/agentic-notes/v0.4.8/HONEST_BASELINE_v0.4.8.md`, `ATOM_COUNT_GATE_v0.4.8.md` |
| Enantiomer collapse / parity filter | `docs/agentic-notes/v0.4.11/LANE-02-donor-fold.md`, `v0.4.12/LANE-01-parity-veto.md` |
| Shipped configuration, every held-off lever's reason | `src/oinsmiles/oin/levers.py` (`_DEFAULT_ON`, `_HELD_OFF`) |
| Gap decomposition and ladder | `docs/agentic-notes/ROADMAP_100_100.md` |

---

## Appendix A — how to update this

Refresh at any release close-out with:

```
/release-retrospective v0.4.13      # append one release (the common case)
/release-retrospective --full       # re-derive every figure from source
```

The command is self-contained at `.claude/commands/release-retrospective.md` — it carries the
sources, the section skeleton, the standing caveats and the guardrails, so a cold session needs
no other context.

### The three artifacts it keeps in sync

| file | what it is |
|---|---|
| `docs/agentic-notes/v0.4.12/RETROSPECTIVE_v0.4.4_to_v0.4.12.md` | this document — the full record |
| `docs/agentic-notes/v0.4.12/retrospective.page.html` | the page **body only** (no `<!doctype>`/`<html>`/`<head>`/`<body>`) |
| `tools/build_retrospective_page.py` | wraps that body into a standalone `index.html` for GitHub Pages |

One body source, two destinations — the Artifact publisher and the `gh-pages` branch — so the
hosted page and the shared link cannot drift. The build script **refuses** a source containing
document-level tags, because that failure would otherwise be invisible until someone opened the
published page.

### Where each file lives — settled, do not re-litigate

| branch | contents |
|---|---|
| **`main`** | **all source** — the report, the page body, the build script, the command |
| **`gh-pages`** | **only the generated `index.html` + `.nojekyll`.** Nothing else, ever. |

The source cannot live anywhere but `main`: a slash command only exists on the branch you have
checked out, so `/release-retrospective` is simply *absent* from a side branch or from `gh-pages`.
The report also cites lane docs, `levers.py` and `CHANGELOG.md`, all on `main` — splitting them
leaves the next updater without one side.

`gh-pages` is an **orphan** branch, no shared history with `main`. That is the point: a session
rewriting `main` cannot disturb the published site. **Never merge `main` into it**, and never put
source on it — everything there is publicly served.

```bash
git worktree add ../oin-ghpages gh-pages
python3 tools/build_retrospective_page.py \
    --source docs/agentic-notes/v0.4.12/retrospective.page.html \
    --out ../oin-ghpages/index.html
git -C ../oin-ghpages commit -am "site: retrospective through <version>"
git -C ../oin-ghpages push origin gh-pages     # authorised; `main` is not
```

### Work in a worktree — but not on a side branch

A `docs/release-retrospective` branch was tried once and was **wrong**. It conflated *isolation
from another session's checkout* — which a **worktree** solves — with *needing separate history*,
which was never true. It was merged to `main` and deleted.

The primary checkout is often on another session's swimlane, and `main` gets rewritten underneath
you; both happened while this report was being written. So: check out **`main`** in a dedicated
worktree (`git worktree add ../oin-retrospective main` — no `-b`), and after **every** commit
re-read `git log --oneline -1` to confirm your commit survived. If it did not, recover it with
`git reflog` + `git cherry-pick`.
