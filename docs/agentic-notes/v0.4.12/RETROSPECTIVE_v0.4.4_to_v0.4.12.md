# Retrospective — has OIN-SMILES gotten more accurate, or faster, since v0.4.4?

<!-- artifact-url: https://claude.ai/code/artifact/9f3a0c17-8f85-4ee6-922b-df599a51978e -->

**Scope:** v0.4.4 (2026-07-23) through v0.4.18 (2026-09-23). Fifteen releases.
**Baseline commit:** `main` @ tag `v0.4.18`, `pyproject` 0.4.18.
⚠ `main` moves under this project — re-read the tip before trusting any SHA here.
**Method:** committed evidence, **plus one 5000-molecule generator sweep run for the v0.4.16
refresh** — the first since v0.4.6, and the reason several figures below moved. v0.4.17 and
v0.4.18 each ran their own sweeps of record as a condition of promoting a default; this refresh
reads them, re-derives the v0.4.18 headline and every runtime row from the frozen
`measurements/` tree, and runs nothing.
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

Then the seventh release moved the score up again — and this is where it gets uncomfortable. It
reported +1.6, and the tool it used **could not have reported a loss even if there had been one**.
It measured by re-scoring molecules the machine had built *months earlier*, so anything it broke in
today's machine was invisible by construction. Checked properly, it had broken 7 molecules while
fixing 78. The same flaw was then found in the release before it: that one's +3.4 is really +2.4.

So the project finally did the expensive thing it had been avoiding — rebuilt all 5000 molecules
from scratch and re-graded them. **The answer is 77 out of 100.** The previously advertised 77.3 was
almost right, but by accident: one error made the gains look too big and another made the baseline
look too small, and they happened to cancel.

**And the machine got much faster without anyone claiming credit for it** — molecules taking over
30 seconds went from 994 to 678, and the typical one from 7.2 seconds to 4.0.

The twelfth release then went looking for the next four points in two obvious-looking places — and
**found that neither could be reached the way everyone assumed.** One built a check that makes sure
the rebuilt molecule hasn't dropped a piece; it changed the answer 51 times and made it right
**zero** times, because *every* candidate the machine had built was broken in the same way. The
other found the machine sometimes builds a perfect **mirror image** — correct in every respect
except handedness, like a left glove where a right one was asked for — and discovered the machine
only ever builds *one* candidate for those, and it is the mirror. **You cannot pick a better answer
out of a bag that only contains the wrong one.**

That release found 48 molecules it *could* fix — but almost none of them were in the place it was
aiming. And taking them makes the machine four times slower on the molecules it touches, so the fix
is built, measured, and deliberately left switched off until someone decides whether that trade is
worth it.

The thirteenth release tried to make that trade affordable. The slowness has an obvious-looking
cause — when the machine can't find a perfect answer it keeps searching to the end of its budget,
and **94% of that searching is done on molecules that never find one**. So: put a limit on how far
past a good-enough answer it keeps looking, and find the sweet spot.

**There is no sweet spot.** Measured at every possible limit, keeping 79% of the benefit costs 89%
of the slowdown. The cost and the benefit rise together, because the *first* few extra attempts are
the expensive ones, not the last. So the trade is exactly as bad as it was, it is now priced, and
nobody has to wonder again.

That release also finally explained a group of 187 molecules that had been filed under "nobody
knows why" for three releases. **In 82% of them the machine built the right skeleton and the
*reader* disagreed about the details** — which bonds are double, where the hydrogens sit. That
matters because the next release had been planned on the assumption they were build failures, and
they mostly are not.

**The useful conclusion: about two-thirds of what's left can't be fixed by choosing better. The
machine has to build differently — and a further chunk may not be a building problem at all, but a
reading one.** That is a much more specific to-do list than "get to 100".

The fourteenth release did something nobody had done: before choosing what to fix, it gave **every
one of the 5000 molecules exactly one reason for failing**, using rulers that never ask the machine
to grade itself. That census found the biggest block was filed under the wrong name — 252 molecules
everyone thought the machine *could not build* were really the *text* being written two different
ways for the same shape. And it found that the earlier "safety check" against left/right mix-ups
had been guarding against a problem the text was creating for itself. Fix the text-writing rule
properly and the check has nothing left to do. **Score: 77 to 83 — the biggest single jump in the
window.** The census also started a stricter second score, where a rebuilt molecule must pass an
independent inspection and not just re-read to the same text: **69 to 75** on that one.

The fifteenth release then took on the biggest thing the census *did* leave with the machine: about
340 molecules where a ring of carbons that should sit *on* the metal was built floating just off
it. Everyone had assumed the machine couldn't find the right answer among its candidates. A ruler
said otherwise: **every candidate was built 0.7 Å too far out**, because two lines of borrowed code
aimed the ring at the edge of the bonding range and then stopped it from stepping closer. Even
ferrocene — the textbook example — was being built wrong, and no test had ever looked. Fixing the
two lines, then teaching the machine to notice when a ring never arrived and try once more, took
the score **83 to 87 on the loose count and 75 to 78 on the strict one**. Each step was predicted
on all 1146 affected molecules before the full run, and the full run landed on the prediction. The
release also found that one of the project's own guard-rail tests **cannot see this whole class of
defect** — it grades the machine on its own bond list — and 36% of the rows it passes fail the
honest check. That is left open, on purpose, with the evidence attached.

So: **more accurate than v0.4.4, yes. More accurate than v0.4.6 — yes, measured by rebuilding
everything, and now graded twice: once by the machine's own re-reading and once by an independent
ruler. And much harder to fool.**

---

## ELI20

Two separate questions, two separate answers.

### Accuracy

`byte_exact` — the share of molecules whose generated 3D structure re-encodes to a byte-identical
OIN string — improved materially across v0.4.4, v0.4.5 and v0.4.6, stopped moving for five
releases, **moved up for the first time at v0.4.13**, and moved twice more, by larger steps, at
v0.4.17 and v0.4.18. **From v0.4.17 the project carries two numbers**: *self-consistent* is
`byte_exact` under the honest score (`E(x) == E(G(E(x)))`); *VERIFIED* is a pass whose string
also parses back to the input's graph *and* whose generated structure a neutral ruler finds
graph-isomorphic to the input with no stereo element mirrored or changed. Every figure before
v0.4.17 in this document is the first kind.

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
| **v0.4.13** | **Real gain, but SMALLER THAN PUBLISHED** — fold + parity veto promoted together | published `byte_exact` **72.46% → 75.88%** (+3.42); **measured ~+2.42** — ~150 real gains against **~30 losses** its instrument could not see |
| **v0.4.14** | **Real gain — `OIN_RESONANCE_DONOR_FOLD` promoted** | **+1.42 pts measured end-to-end** (78 gains, **7 losses**, n=182 of 182). Also re-filed **5.94 pts** off the encoder ladder without moving a molecule |
| **the sweep** | **First real generator sweep since v0.4.6** | **`byte_exact` 3858/5000 = 77.16%** honest (86.88% scored — a **9.72-pt** gap) |
| **v0.4.15** | **Zero, by design — both lanes measured, both ship OFF** | L1 **0 gains of 289** (51 structures moved, 0 improved); L2 **+48 available** (+0.96 pts, 0 losses) but held at **4.00×** runtime |
| **v0.4.16** | **Zero, by design — the release prices a decision rather than moving a number** | bounding the search is **refuted as a rescue**: keeping 79% of the +48 costs **89%** of the runtime penalty. Separately, the 187 "unexplained" molecules are **82% perception, not construction** |
| **v0.4.17** | **Real gain — `OIN_EXACT_DONOR_FOLD` promoted; the largest step in the window, and the first measured end-to-end *before* promotion** | sweep: self-consistent **77.16% → 82.72%**, VERIFIED **69.24% → 74.60%**. The lever isolated on its complete 504-molecule mover set: **+273 / +265**, 3 / 4 losses. `slot_renumber` **252 → 5**. The census that preceded it moved 5.04 pts from "the generator cannot build it" to "the encoder writes it two ways" |
| **v0.4.18** | **Real gain — three generator-side levers promoted in two steps** | release sweep: self-consistent **82.72% → 86.94%**, VERIFIED **74.60% → 78.40%**. Each step predicted by an A/B over all **1,146 η-bound** molecules at a **noise floor of zero** and landed to within the 300 s budget boundary. `DETACHED` was a **distance** (+0.73 Å). **39 + 5 verified passes lost and accepted** |

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

**v0.4.14 then found the flaw in that reasoning, and it applies to both releases.** The licence for
measuring offline was `tools/fold_key_invariance.py` reading **0 comparison keys changed** — argued
to mean the generator returns the same conformers either way. It does not. `accept_fn` decides by
key, so key-invariance bounds **acceptance**; but the generator's input is the OIN **string**, so a
slot relabeling changes `ParsedOIN`, the CoordMap, and **the pool itself**. And because an offline
re-score holds the generated structure *fixed*, for a molecule that already round-trips it can only
print "still fine" — **it reports `bad_direction = 0` whether or not losses exist.** Measured
end-to-end, v0.4.13 cost **~30 losses** against ~150 real gains (**~+2.42, not +3.42**) and v0.4.14
cost **7** against 78 (**+1.42, not +1.56**). Neither release was wrong to skip a sweep — v0.4.14
showed the affected population is *derivable*, so an exact A/B is cheaper — but both were wrong
about *why* they could.

**v0.4.17 and v0.4.18 are the first releases whose gains were measured the way that lesson
demands, and the first that get the same scrutiny in the other direction.** Both ran a real
harness A/B over the *complete* affected population — 504 movers derived from coordinates, then
all 1,146 η-bound molecules — with both arms generated, both numbers reported, and a full sweep
owed and paid at each promotion. The sweeps landed on the predictions: 82.72 against 82.62
predicted, 86.10 against 86.04, 86.94 against 87.00, every residual a row at the 300 s budget. But
v0.4.17 also retired evidence this document had been carrying as settled: the veto's three
"oracle-confirmed enantiomer" fixtures — the molecules v0.4.12 built it on and v0.4.13 promoted it
against — turned out to encode to their *mirror's* string 8, 7 and 6 times in 12 renumberings
under the shipped encoder. **One renumbering per fixture had been drawn. A coin, tossed three
times.** v0.4.11's "221 of 393 collapse enantiomers" and v0.4.12's "19 → 0" were reading it. The
fold was over-folding, the veto was patching the over-fold, and the fix needed neither.

### Speed

Three real default-path wins since v0.4.4 through v0.4.14, then two more that were measured
corpus-wide and like-for-like — **and the corpus figure is current: every sweep since v0.4.14 has
run in one configuration.**

| release | default-path speed change |
|---|---|
| **v0.4.4** | `OIN_EARLY_EXIT` default-ON — **~5×** on the worst cohort |
| **v0.4.5** | Removed a duplicate **48–57 s** re-encode running once per rejected conformer (unconditional). First encoder profile ever taken: `AC2BO` is **99.8%** of a slow encode |
| **v0.4.6–v0.4.9** | **None on the default path.** The two big candidates were held OFF on measurement, not caution |
| **v0.4.10** | Deleted a discarded `.index()` scan that eigendecomposed a Coulomb matrix per candidate — **on by default, no lever**: `CAHQEJ_comp_0` **−32.9%**, `FOSNEI_comp_0` **+0.3% (nil)** |
| **v0.4.11–v0.4.13** | None claimed. v0.4.13's promotion was argued **generator-neutral by measurement** — 0 of 9669 keys moved. ⚠ v0.4.14 showed that argument bounds *acceptance*, not *embedding*: the generator consumes the **string**, so runtime **can** move. It was never measured either way |
| **v0.4.14** | None claimed on the default path. Measured on the 182 molecules its lever touches: **+57.4% total CPU**, `> 30 s` 29 → 36, worst molecule **+511.5 s** — a real cost, on a small population |
| **v0.4.15** | **None on the default path — both levers ship OFF.** L1 measured **free** (1.01×, 0 gains). L2 measured **4.00×** on its 365 molecules, `> 30 s` **30 → 122** there (~678 → ~770 corpus-wide): the reason its +0.96 pts is *held*, not taken |
| **v0.4.16** | **None, and none available — that is the release's result.** The full curve over all 365: bound 3 holds `> 30 s` to **52** but keeps only 19 of 48; bound 12 keeps 38 and costs `> 30 s` **104**. The frontier is close to linear, so a limit moves *along* v0.4.15's trade instead of improving it |
| **v0.4.17** | **Real, corpus-wide, and like-for-like with the v0.4.14 sweep.** The parity veto's two extra encodes had been paid *per conformer* — `accept_fn` re-encodes every conformer in the pool — and the exact fold needs neither. Σ `elapsed_s` **38.7 → 33.2 h**, median **4.01 → 3.32 s**, `> 30 s` **678 → 579**. Over the lever's own 504 movers: Σ 6.19 → 2.93 h, `> 30 s` 110 → 58, and 16 molecules that used to be killed at 300 s because the string asked for an unbuildable arrangement now finish in seconds |
| **v0.4.18** | **Real, small, and reported both ways.** Σ **33.18 → 31.95 h**, `> 30 s` **579 → 521**, `> 300 s` 180 → 173, max 698 → 628 s. ⚠ The median moved the other way, **3.32 → 3.53 s**: the tail came down, the middle did not. On the η cohort the pair's A/B read Σ 14.03 → 11.77 h and the retarget's read flat (11.63 → 11.81 h); why the corpus median rose was not measured |

**The headline speed figure is no longer stale.** The v0.4.6-era `994/5000 = 19.88% over 30 s,
median 7.19 s` stood unmeasured for eight releases. The sweeps of record since read:

| corpus runtime, N = 5000, from the **nested** `metrics.elapsed_s` | v0.4.6 sweep | v0.4.14 sweep | v0.4.17 sweep | **v0.4.18 release sweep** |
|---|---:|---:|---:|---:|
| `> 30 s` | 994 (19.88%) | 678 (13.56%) | 579 (11.58%) | **521 (10.42%)** |
| median | 7.19 s | 4.01 s | 3.32 s | **3.53 s** |
| Σ `elapsed_s` | — | 38.7 h | 33.18 h | **31.95 h** |
| max | 759.9 s | 728.8 s | 697.9 s | **628.5 s** |

**The −6.3 points from v0.4.6 to v0.4.14 nobody claimed.** No release in that window asserted a
corpus speed win; three landed per-molecule optimisations and the rest measured none. ⚠ **That
step is not strictly like-for-like:** the v0.4.14 run capped BLAS threads to 1
(`OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS=1`) and the v0.4.6 run did not — a deliberate choice,
because `OIN3DGenerator(timeout=)` is *advisory* and CPU starvation shrinks the embed pool, which
would have biased **accuracy**. Treat that direction as solid and the magnitude as approximate.
**The three columns from v0.4.14 on *are* like-for-like** — same cohort, seed, 6 shards, 300 s
budget, BLAS = 1 — and the v0.4.17 and v0.4.18 rows were re-derived from the frozen per-molecule
extracts for this refresh. Goal B (`max(elapsed_s) < 30 s`) is still **not delivered**.

### What actually improved most

Neither number. What improved most is the project's ability to tell a real result from a fake one.
Between v0.4.4 and v0.4.18 the unit suite went **551 → 1076**, and the instruments added in that
window — the honest re-score, the corpus encoder-identity gate, the two-arm byte-identity gate, the
coordination-integrity check, the mirror audit, and from v0.4.17 the four-instrument census and its
neutral ruler — are what caught the 10.34-point inflation, the 59%-false-positive cohort, the dead
gate arm, the enantiomer collapse, and then the fact that the "enantiomer collapse" had been
measured on a coin. Six releases that moved no points each removed a way of being wrong.

**v0.4.15 added the sharpest one yet, and it is a rule rather than a tool.** Its six A/B arms first
returned a flawless *"0 gains, 0 losses"* over 1107 molecule-pairs — because a `sys.path.insert`
inside the measuring tool silently overrode `PYTHONPATH` and both arms ran identical code. The
standing rule *"ask what a broken version would print"* does not catch that: **here the broken and
the working version print the same thing.** What caught it was a **second** instrument built to
verify the *mechanism* (does the lever fire?) rather than the *outcome* (did the metric move?).
Two instruments disagreeing is the signal; either alone reads clean.

**v0.4.13 is the return on that, and v0.4.14 is the bill.** The gain v0.4.13 banks is the fix
v0.4.11 built and refused, shippable only once an instrument existed that could see the damage — but
v0.4.14 then showed that *instrument itself* could not see losses, and corrected both releases'
headlines downward. The pattern is now four deep: v0.4.8 corrected the scoring predicate, v0.4.11
corrected the fold, v0.4.14 corrected the *measurement method*, and each correction was found by
building the instrument the previous one lacked. v0.4.13 caught **four instruments printing plausible
nothing** — including one that reported the *refuted*
number under the correct heading, and a byte-identity gate whose PASS was meaningless because
**0 of its 62 fixtures were among the molecules the change moved**. The lesson generalises past
this project: *a gate that cannot see your change is not evidence about your change* — state its
coverage of the moved population, not its verdict.

**v0.4.16 hit that same rule twice more, and added a harder case.** Its confirmation arm was killed
partway through a population that happened to be *sorted*, so every completed row was a molecule the
change keeps: **30/30 agreement, confirming nothing**, because an arm containing only cases your
mechanism accepts cannot discriminate. And its classifier was a **normalizer** — the component that
decides whether two things are "the same" — which is the one place "ask what a broken version would
print" fails hardest: a string-level heavy-atom comparison disagreed with a canonical one on **57 of
109** molecules, a coin flip, and printed a clean table that **agreed with the roadmap's existing
assumption**. It was caught only because a *read example* contradicted it. **When a normalizer
decides your headline, validate it against an independent canonical comparison** — not against your
reading of a handful of cases, and least of all against your expectations.

**v0.4.17 added the rule that makes a headline mean something: carry a second number the lever
cannot teach.** `byte_exact` is self-agreement — `E(x) == E(G(E(x)))` — and a canonicality lever
can raise it by making the second encode agree with whatever was built. The census's VERIFIED
predicate judges the *structure* with a ruler that never imports the encoder; a lever that moves
both together is doing what it claims, and one that moves only the first is not. It also found that
**a gate re-frozen from a predicted mover list had carried 15 wrong rows for four releases**,
because the list was a sample frame and the gate runs in full only at a release. **A mover list is
a sample frame, not a verdict** — hence `/refreeze-goldens`, from a full run, every promotion.

**v0.4.18 then used that full run to find what the gate itself is blind to.** ARM 2's round-trip
field is computed through the generator's own bond graph — the pre-v0.4.8 *scored* path — so
whether a ring sits 2.1 Å or 2.9 Å from its metal never enters. **138 of 198 η rows built a
different structure and the gated field moved on 0.** A construction regression passes ARM 2, and
did, for every release that shipped ferrocene detached. The release added an honest observation
column beside the gated one and left the decision — re-freeze ~135 passing rows as honest fails, or
change what the gate means — with the owner and the evidence. The general form: **a gate that
scores an artifact by the artifact's own account of itself cannot see the class of defect where the
account is wrong**, and this project has now met it in the harness (v0.4.8) and in its own
regression gate (v0.4.18).

---

## The two questions, answered directly

> ### Is it more accurate?
>
> **Yes — 86.94% self-consistent and 78.40% VERIFIED at v0.4.18, measured by a full sweep at the
> shipped defaults, up from 77.16% / 69.24% at v0.4.14 on the same cohort in the same
> configuration.** Two releases account for all of it: v0.4.17 (+5.56 / +5.36, the encoder's donor
> fold) and v0.4.18 (+4.22 / +3.80, the generator's η target). Each was predicted by an A/B over
> its complete affected population before the sweep, and each sweep landed on the prediction to
> within the 300 s budget boundary. The VERIFIED number is the one to quote: it is the first in this
> document that a lever cannot raise by teaching the re-encode to agree with a wrong structure.
>
> The long answer is the interesting one. From v0.4.6 to v0.4.12 the default path did not change a
> single answer, and the reported figure *fell* from 82.80% to 72.46% at v0.4.8 because the
> measurement was corrected, not because the software regressed. **v0.4.13 is the first release in
> the window to move the number up**, and it did so by shipping a fix built, measured, refuted and
> shelved two releases earlier. v0.4.14 moved it again — and, more consequentially, showed that both
> releases' headlines had been measured with an instrument that **could not report a loss**.
> Corrected: v0.4.13 is **~+2.42**, not +3.42; v0.4.14 is **+1.42**, not +1.56.
>
> **A 5000-molecule sweep run for this refresh settles the absolute at 77.16%** — the first real
> generator sweep since v0.4.6, ending an eight-release chain of offline re-scores. The previously
> published 77.30% was only 0.14 off, but by **two errors cancelling**: the chain over-stated the
> lever deltas and under-stated the base.
>
> **v0.4.15 then tested both of those decompositions by building them — and refuted both.**
> "A one-site guard, not a capability limit" was wrong: the return-path guard **changed the
> returned structure 51 times across 289 molecules and improved it zero times.** The mirror class
> was wrong in the other direction: it is real, but **1 of 201 is recoverable** because the pool
> holds exactly one key-matching conformer and it *is* the mirror.
>
> **So the gap is now partitioned by measured reachability, not by bucket name:** 0.96 points are
> reachable today (behind a lever held off on cost), **10.04 points are proven unreachable by any
> selection predicate**, and 5.54 produce no structure at all.
> **15.36 of 22.84 points (67%) require the generator to BUILD something different** — not to
> choose differently, and not to emit differently.
>
> **v0.4.16 closed the last uncharacterised block and split it three ways.** The 3.74 points that
> had never been looked at are now: **2.82 PERCEPTION** — the generated heavy-atom skeleton is
> already correct and the re-reading disagrees about bond orders, aromaticity or hydrogens —
> **0.50 construction**, and **0.12 stereo inversion**. That does not move `byte_exact`, but it
> moves the *schedule*: a release had been sized to take all of it as build work.
>
> It also priced the one reachable point and found the price fixed. **Bounding the search is
> refuted as a rescue**: keeping 79% of the +48 costs 89% of the runtime penalty, because the
> *early* extra attempts are the expensive ones. The +0.96 is still available, still costs ~1.8
> points of the speed goal, and is now a standing decision rather than an open engineering task.
>
> **v0.4.17 then found that the partition above had the biggest block under the wrong owner.** A
> census gave every one of the 5000 molecules exactly one fault with instruments that never route
> through the encoder to judge the generator or the reverse. The 201 "enantiomers the generator
> cannot build" were **9 molecules**; the 252 `slot_renumber` were the *encoder* writing an achiral
> molecule two ways, because the donor fold swapped symmetry *buckets* rather than ligand
> *automorphisms* — the same over-fold that had collapsed enantiomers at v0.4.11 and made a veto
> necessary at v0.4.12. Folding over true automorphisms needs neither: **+273 self-consistent /
> +265 VERIFIED** on the complete 504-molecule mover set, 3 / 4 losses, and a sweep at
> **82.72% / 74.60%** against 82.62 / 74.54 predicted.
>
> **v0.4.18 took the block the census left with the generator and found it was a distance.** The
> 343 `DETACHED` molecules — 86% η-bound — had every candidate in the pool built with the ring a
> median **+0.73 Å** too far from the metal, because upstream's cleaner aimed an η group at 1.1–1.2
> × the covalent sum (the encoder's bonding cutoff is ~1.20 ×) and its clash guard reverted the
> first step inward on Fe, Co and Ni. Ferrocene itself was built at Fe–C 2.85 Å against a real
> 2.05 Å. Two promotions — the target + exemption pair, then a per-molecule retarget for the rings
> the pair left stalled — and the release sweep reads **86.94% / 78.40%**, +4.22 / +3.80 over
> v0.4.17, with 39 + 5 verified passes lost and accepted. **The gap is now 13.06 pts
> self-consistent and 21.60 VERIFIED, every molecule attributed, and the three largest blocks —
> construction 3.94, perception 3.52, serializer 2.96 — are within a point of each other.**

> ### Is it faster?
>
> **Yes, corpus-wide, and since v0.4.14 the sweeps are like-for-like — but v0.4.15 found where the
> next accuracy point would cost speed.** `> 30 s` went **994 → 678 → 579 → 521** across the
> v0.4.6, v0.4.14, v0.4.17 and v0.4.18 sweeps (19.88% → **10.42%**), Σ `elapsed_s` 38.7 → 33.18
> → 31.95 h from v0.4.14 on, and max 759.9 → 628.5 s. The first step nobody claimed: three
> default-path optimisations landed (v0.4.4, v0.4.5, v0.4.10) and every release to v0.4.14
> measured none. The next two were measured: v0.4.17 stopped paying the parity veto's two extra
> encodes per conformer (median 4.01 → 3.32 s), and v0.4.18's attached rings finish inside the
> budget more often (`> 300 s` 180 → 173). ⚠ v0.4.18's median rose 3.32 → 3.53 s while its tail
> fell; the cause was not measured. ⚠ The v0.4.6 → v0.4.14 step crosses a thread-configuration
> change and its magnitude is approximate. The largest single per-molecule win (`VAFMIA_comp_0`,
> 81.89 s → 10.87 s) still sits behind a lever that ships **off**. Goal B (`max(elapsed_s) < 30 s`)
> is **not delivered** — 521 molecules are over it.
>
> **v0.4.15 made the trade explicit for the first time.** Its Lane 2 recovers +0.96 accuracy points
> and costs **4.00×** runtime on the molecules it touches — `> 30 s` 30 → 122 there, ~678 → ~770
> corpus-wide. Roughly **+1 point of Goal A against ~+1.8 points of Goal B.** The two goals have
> been stated together since v0.4.4; this is the first release where a shipped lever had to choose
> between them, and it is why the lever is held rather than taken.
>
> **v0.4.16 tried to remove that choice and could not.** The obvious move is to cap how far past a
> good-enough answer the generator keeps searching — 94% of the extra cost is spent on molecules
> that never find a better one. Measured at every bound over all 365: **there is no favourable
> setting.** Keeping 79% of the accuracy costs 89% of the slowdown; the two rise together, because
> the *first* extra conformers are the expensive ones (each is a full embed plus a full
> re-perception), not the last. **The trade is exactly as bad as v0.4.15 found it — and now it is
> priced, so no future release need re-open it as an engineering question.**

---

## How to read the numbers in this document

Four traps, all of which this project has already fallen into once and documented:

1. **The bucket reports are not a time series.** They are drawn from five different cohorts —
   6719 (v0.4.2 capstone), 3917 (v0.4.4 regression), 936 (v0.4.5 rebaseline), 5000 (v0.4.6 seed-42
   sweep), and 5000 (the **v0.4.14 baseline sweep**, which the v0.4.17, v0.4.18 and v0.4.18
   *release* sweeps re-use). Every accuracy figure below carries its N and its cohort. **Do not
   draw a line between two of them.**
   Two comparisons *are* genuine like-for-like: scored-vs-honest at v0.4.8 and again at v0.4.14,
   each classifying *the same 5000 reports* two ways (**10.34** and **9.72** points respectively).
   ⚠ The v0.4.6 and v0.4.14 sweeps share a cohort and a seed but **not** a code version or a thread
   configuration — closer than anything before them, and still not a series.
   ✅ **From v0.4.14 on the sweeps of record *are* a series**: one cohort, one seed, 6 shards,
   `--mol-timeout 300`, BLAS threads = 1, and the only variable is the code under test. 77.16 →
   82.72 → 86.10 → 86.94 self-consistent (69.24 → 74.60 → 77.66 → 78.40 VERIFIED) is the first
   line in this document that may legitimately be drawn — and each step was independently
   predicted by an A/B on its affected population before the sweep confirmed it.

   ✅ **Every figure in that list is now independently checkable, which it was not before.** All
   five cohorts are frozen as per-molecule extracts in `measurements/`, and each re-derives its
   published headline exactly through the real classifier rather than a re-implementation —
   81.19 / 44.91 / 60.26 / 82.80 / 72.46 / 77.16, six for six. The v0.4.18 release headline
   re-derives from its frozen ruler verdicts, parse-back and honest bucket report alone
   (5000 / 4347 / 3920, checked for this refresh), and the v0.4.17 and v0.4.18 runtime rows
   re-derive from their extracts. **The sample definitions are frozen too**
   (`measurements/cohorts/`): seed, N, the dedup priority, and every molecule name. Until
   v0.4.16 those manifests were untracked, so the corpus under every headline this project has
   published lived one `rm -rf` from being unreproducible — while the *lane* populations had been
   frozen since v0.4.15 under the rule "a rate without its sample is not reproducible". The rule
   was right; it had just never been pointed at the corpus.

2. **`metrics.elapsed_s` is nested and is a SUM.** Read from the top level it silently yields `0`.
   It also accumulates up to three separately SIGKILLed harness attempts, so the retired headline
   "max 759.9 s against a 300 s budget" was arithmetic on a sum — all 4658 single-attempt rows in
   the 5k sweep finish within **0.2 s** of their cap. The harness enforces to ε ≈ 0.2 s.
   *(`v0.4.9/ELAPSED_S_IS_A_SUM_v0.4.9.md`)*

3. **Suite count is a rigour proxy, not an accuracy metric.** 551 → 1076 tests means the project
   can detect more kinds of wrongness. It does not mean the notation got better. And a green suite
   is only as good as what it looks at: **1,067 of 1,067 tests passed while ferrocene was being
   built with its rings 0.8 Å off the iron**, because no test measured that distance.

4. **From v0.4.17 there are two numbers, and `byte_exact` is the weaker one.** *Self-consistent*
   (`byte_exact` under the honest score) is self-agreement: `E(x) == E(G(E(x)))`. *VERIFIED*
   additionally requires that the string parses back to the input's graph and that a neutral ruler
   — one that never imports the encoder — finds the generated structure graph-isomorphic to the
   input with no stereo element mirrored or changed. The gap between them (406 passes at v0.4.17,
   427 at v0.4.18) is accuracy the metric *asserts*, not accuracy the project *has*. A lever must
   move both, or it may be teaching the second encode to agree with a wrong structure. Every figure
   in this document before v0.4.17 is self-consistent only; the 77.16% baseline is 69.24% VERIFIED.

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

**Accuracy: REAL GAIN, BUT SMALLER THAN PUBLISHED. Published +3.42 (72.46% → 75.88%); measured
~+2.42 (→ ~74.88%).** The release reported 171 molecules moving `key_equal/slot_renumber →
byte_exact` with **0** in a bad direction. **The 0 is the artifact**: it came from an offline
re-score, which holds the generated structure fixed and therefore *cannot express a loss* — for a
molecule that already round-trips it can only ever print "still fine".

Measured end-to-end by v0.4.14 (`tools/generator_ab_honest.py` — real generation both arms, scored
by re-perceiving the written XYZ): gains **22/25 = 88%** (seed 17) ⇒ ~150 real; losses **6/40 = 15%**
(seed 13) over a **197-molecule at-risk population** — movers that were *already* `byte_exact`, which
nobody sampled — ⇒ **~30 losses**. Still the first upward movement in the window, and still real.

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

> 🔴 **v0.4.14 FOUND THE HOLE IN THIS ARGUMENT.** `accept_fn` decides by key, so key-invariance
> bounds **acceptance** — but the generator's input is the OIN **string**, so a slot relabeling
> changes `ParsedOIN`, the CoordMap, and **the pool itself**. Key-invariance says nothing about
> **embedding**. Skipping the sweep was still the right call — v0.4.14 showed the affected population
> is *derivable*, making an exact A/B far cheaper than 55 CPU-h — but the stated reason was wrong,
> and the headline it licensed was over-stated by ~1 point.

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

### v0.4.14 — the release that corrected its own instrument, and the sweep that ended the chain

**Status: RELEASED, tag `v0.4.14`, `pyproject` 0.4.14.**

**What shipped.** **PROMOTED to default-ON:** `OIN_RESONANCE_DONOR_FOLD` — it widens the donor
fold's equivalence test to the fragment's *constitutional skeleton* (bond orders, aromatic flags,
charges and hydrogens erased; connectivity, element and chiral tag kept), so acac written
ketone/enol, carboxylates and sulfonates stop reading as two inequivalent donors. Ester `-O-`/`=O`,
ether/ketone and amide N/O still do not merge. **Nothing else promoted**; `OIN_PREFILTER_ADVISORY`
stays off.

**Accuracy: REAL GAIN — +1.42 points, 78 gains against 7 losses, net +71 molecules.** Measured
**end-to-end over all 182 molecules the lever can affect** (`n = 182 of 182`, nothing sampled). The
release *first* reported **+1.56 with `bad_direction = 0`**, from an offline re-score, and correcting
that is the release's real content.

**Refuted — its own measurement method, and the previous release's.** The licence for measuring
offline was `tools/fold_key_invariance.py` reading **0 comparison keys changed**. That bounds
**acceptance** — `accept_fn` decides by key — but the generator's input is the OIN **string**, so a
slot relabeling changes `ParsedOIN`, the CoordMap and **the pool itself**. And an offline re-score
holds the generated structure *fixed*, so for a molecule that already round-trips it can only print
"still fine":

> **An offline re-score reports `bad_direction = 0` whether or not losses exist. Ask what a broken
> version of an instrument would print — here, the same thing a working one does.**

The population nobody had sampled is the one that can only lose: **movers that were already
`byte_exact`** — 39 for v0.4.14, **197 for v0.4.13**. Both releases' headlines came down:
v0.4.14 **+1.56 → +1.42**, v0.4.13 **+3.42 → ~+2.42**.

**And no sweep was needed to find that out.** The affected population is **derivable**: a molecule
whose `encode(input)` is byte-identical in both arms hands the generator the same string, and
generation is **seeded**, so it is unchanged *by construction*. `tools/lever_string_movers.py`
derives it from coordinates — **93** of 5000 move the input string, 182 including the generated side
— and an A/B over that set is *exact*. ⚠ Derive it from **coordinates, never from a frozen sweep's
stored strings**: that route read **179** and was wrong in both directions (89 do not move under
today's encoder; **3 that do were missing**), because a stored `smiles_1` was emitted by the
v0.4.8-era encoder.

**Two buckets whose NAMES were hypotheses.** `key_equal` is documented as *"benign canonicalization
— the win reclaimed"*; **183 of its 361 members (50.7%) are the generator building the ENANTIOMER**,
invisible because `compare._parse_vertex_colors` folds reflection deliberately *and `accept_fn`
decides by that key* (222/222 classified, 0 excluded). And `rdkit_canonical` is **80.7% η-set
denticity drift** — the generated ring slips and fewer carbons fall inside the bonding cutoff — not
canonicality at all. **5.94 points re-filed off the encoder ladder without a single molecule
changing bucket**, leaving it **~1.28 points** of reachable work.

**Lane 2 (measurement only): `PREFILTER_VETO`'s prevalence, finally sized.** v0.4.13 could measure
only `n = 1`. Over a stratified 50-molecule sample (seed 7; 49 measured, 1 hard timeout, **5
`INSTRUMENT_DEAD` excluded — live denominator 44**): the cheap prefilter vetoed **261** conformers
and the strict test disagreed on **4** — **1.5% of vetoes, 4.5% of molecules**. `AROHIA_comp_0`'s
`0/48`-vs-`16/48` is an **outlier**. It stays off — but **its latency objection is REFUTED**: it is
measured *faster* (−498 s of which is 2 molecules exiting early; **+19.9 s** across the 42 where it
did nothing). Anyone citing cost to keep it off is citing a refuted number.

**What it cost.** 7 molecules that round-tripped no longer do. They are **not** the encoder emitting
a wrong string — it is canonical in both arms; the **generator**, handed an equally valid but
differently-*labelled* input, builds a worse structure. *The generator's output depends on the slot
labeling of its input, and it should not.* Every future canonicalization lever pays this toll, which
is now the strongest argument for fixing that before spending more releases on canonicality. Also:
arm2 goldens re-frozen (7 of 325, 1 of 100, fields 1–6 only — splicing whole rows silently breaks
`--band`); ARM 1 PASS byte-identical at **0/62** coverage.

Suite: **1007 OK**.

---

### The baseline sweep — run for this refresh, and it ends an eight-release chain

**Not a release.** A 5000-molecule generator sweep at the shipped v0.4.14 defaults, run because
every `byte_exact` figure since v0.4.6 was an offline re-score of v0.4.6-era structures while
v0.4.9–v0.4.14 all changed what the generator builds. **The deltas were sound; the absolute had
drifted.**

| | v0.4.6 sweep | **v0.4.14 sweep** |
|---|---:|---:|
| `byte_exact` (honest) | — | **3858 / 5000 = 77.16%** |
| `byte_exact` (scored) | — | 4344 = 86.88% |
| scored − honest | 10.34 pts (at v0.4.8) | **9.72 pts** |
| `> 30 s` | 994 (19.88%) | **678 (13.56%)** |
| median | 7.19 s | **4.01 s** |

**The published 77.30% was 0.14 points off — by two errors cancelling.** The offline chain
over-stated the lever deltas (losses invisible) *and* under-stated the absolute (v0.4.6-era
structures are worse than what the current generator builds). That is luck, not method.

**Runtime improved by 6.3 points and nobody claimed it.** No release in the window asserted a corpus
speed win. ⚠ BLAS threads were capped to 1 here and not in the v0.4.6 run — deliberate, because
`OIN3DGenerator(timeout=)` is *advisory*, so CPU starvation shrinks the embed pool and would have
biased **accuracy**. Direction solid, magnitude approximate.

⚠ **Preserved as an extract, not whole.** The run is 261 MB and gitignored; the bucket JSONs are
2 MB against a 512 KB per-file cap. `measurements/v0.4.14-sweep/per_molecule_extract.tsv` carries
5000 rows (bucket under *both* scorings, subclass, nested `elapsed_s`) and was **verified to
re-derive 3858 / 77.16% / median 4.01 s / >30 s 678 exactly**.

---

### v0.4.15 — two lanes built to spec, and the data refuted both targets

**What shipped.** `OIN_ATTACH_RETURN` + `OIN_ATTACH_RETURN_STRICT` (Lane 1) and
`OIN_ACCEPT_STRING_EXACT` (Lane 2), all three **added and default-OFF**. **Nothing was promoted, so
nothing changed for a user.** The default path is byte-identical: 1039 tests pass with every lever
unset.

**Accuracy: Zero, by design — and the zero is two different findings.**

*Lane 1* wired v0.4.7's coordinate-only attachment predicate into the **return** path, closing the
gap `OIN_ATTACH_CHECK`'s own lever entry names as its residual ("the check guards ACCEPTANCE, not
RETURN"). Measured over the 289 molecules where the guard's own predicate fires: **0 gains,
0 losses, 51 structures changed, runtime 1.01×.** The lever fires hard — telemetry shows **10–16
winding-matching conformers per molecule and every one detached** — promotes a better-attached
conformer 51 times, and improves the answer zero times. **Attachment is necessary, not sufficient**,
which had been a single-molecule anecdote (MEDZUR) and is now a measured class.

*Lane 2* replaced the charter's approach outright. The chartered fix — use `metal_config` as an
acceptance predicate — **cannot work**: acceptance needs a reference handedness and the generator's
only input is the OIN string, which is why the helicity branch already in the adapter is dead code.
Measured replacement: over the 183 known `MIRROR_MATCH` molecules, **normalized strings differ
183/183 while the comparison keys agree 183/183** — the handedness survives normalization and is
folded only by the key. Result: **+48 molecules (+0.96 pts), 0 losses over 565.**

**The decomposition is the finding.** Of the 48 gains, **43 are `rdkit_canonical` and 5 are
`slot_renumber`** — **28.7%** recovery on the non-enantiomer 164 against **0.5%** (1 of 201) on the
enantiomers. **57× apart.** Everything the lane recovers is *outside* the class it was aimed at;
scoped as chartered, the release would have shipped **+1 molecule**.

**Speed: none on the default path.** Lane 1 is free (1.01×). Lane 2 costs **4.00×** on its 365
molecules, `> 30 s` **30 → 122** there — which is why its +0.96 points is held rather than taken.

**Refuted.**
- 🔴 **The charter's own headline claim.** It opened, in three files, with "the absolute baseline is
  NOT trustworthy — v0.4.13 really measures ~74.88%". The sweep reads **77.16%**: over-stated by
  **0.14** points, not 2.4. What the offline chain drifted was bucket *composition*
  (`structural` +67, `hard_fail` −53), which an offline re-score cannot see because it holds the
  structure fixed.
- 🔴 **v0.4.14's write-off of `rdkit_canonical`.** It was re-filed as "80.7% η-set denticity drift,
  not reachable by canonicalization, do not re-propose a string fix" and removed from the ladder.
  Lane 2 recovered **43 of its 113 — 38.1%, the highest rate of any block measured.** Both
  statements are true: v0.4.14's was about the *encoder* and stands; the block is reachable from the
  *generator*. **Reachability is a property of a mechanism, not of a block** — and the project made
  that error in both directions in consecutive releases.
- 🔴 **Emitting `|mc:±|` as v0.4.16's plan.** `pool.accept_incumbent_recorded = 1`: exactly one
  key-matching conformer exists and it is the mirror. Emitting would make the generator responsible
  for a handedness it does not build — 201 silent wrong answers become 201 loud failures for no
  gain.

**What it cost.** 🔴 **All six A/B arms were run twice, because the first six measured the wrong
tree.** `generator_ab_honest.py` does `sys.path.insert(0, <its own>/../src)`, which **overrides
`PYTHONPATH`** — so every arm imported `main`'s `oinsmiles`, where neither lever exists, both sides
of each A/B ran identical code, and the result was a flawless *"0 gains, 0 losses, output moved 0"*
across 1107 molecule-pairs. It was believed for an hour. It was caught only because a **second**
instrument (`selection_pool_probe.py`, written to prove a lever *fires* before believing any null)
contradicted the telemetry.

**The standing rule "ask what a broken version would print" is necessary and not sufficient** — here
the broken and working versions print the same thing. The defence that worked was a second
instrument measuring the *mechanism* rather than the *outcome*. Third occurrence of this family
after v0.4.9's sibling venv and v0.4.13's sibling glob; now guarded by a hard refusal rather than a
note. Also caught before publication: the harvester would have **silently dropped every arm JSON**
(no pattern matched `lane1_pop_*.json`) and then, once fixed, would have let the *void* arms
overwrite the real ones — they share basenames.

**No full sweep, deliberately.** Both levers ship OFF, so the shipped default is unchanged and the
v0.4.14 baseline sweep *is* v0.4.15's measurement. Re-running 38.7 CPU-h to reproduce an unchanged
default would be waste. **If `OIN_ACCEPT_STRING_EXACT` is promoted, a sweep becomes mandatory.**
No combined arm, for the same reason: with Lane 1 off, "both" is not a shipping configuration.

Suite: **1039 OK**.

---

### v0.4.16 — the trade gets a price, and the price does not come down

**What shipped.** `OIN_STRING_EXACT_BOUND` **added and default-unset**, plus `lever_int()` beside
`lever_enabled()` for levers where `0` is a meaningful value. **Nothing was promoted, so nothing
changed for a user.** 11 levers default-ON, 17 held OFF — unchanged from v0.4.15.

**Accuracy: Zero, by design — the release exists to price a decision, not to move a number.**

v0.4.15 left **+48 molecules (+0.96 pts, 0 losses)** behind a **4.00×** runtime cost, on the theory
that the cost was a reclaimable tail: the lever declines to *stop* the pool, and the **317
molecules that never gain consume 93.9% of its bill**. v0.4.16 built the bound that would reclaim
it and measured the whole curve over all 365 molecules.

**The two halves of the promotion bar are mutually unreachable:**

| | needs | and then gives |
|---|---|---|
| keep ≥ 75% of the 48 | bound **≥ 12** | `> 30 s` = **104** against a limit of 52 |
| `> 30 s` ≤ 52 | bound **≤ 3** | recovered **19** of 48 |

**The shape is the finding, not the threshold.** To keep 79% of the gain you pay 89% of the runtime
penalty — the frontier is close to linear, so bounding moves *along* the v0.4.15 trade rather than
improving it. **The charter's mechanism hypothesis is refuted**: the tail was never wasted
(`early_hit` already stopped both fill loops), and the cost is that each extra conformer is a full
embed **plus** a full re-perception, so the *early* extras are the expensive ones.

**Speed: none on the default path**, and none available. That is the point of the release.

**The second lane characterised a block the roadmap had carried for three releases as
*"nobody knows why"*.** All **187** molecules classified, 0 unaccounted, every class read against
an example. Over the 172 `structural` INTACT+BOUNDARY:

**PERCEPTION 141 (82.0%) · CONSTRUCTION 25 (14.5%) · STEREO 6 (3.5%)**

The heavy-atom graph is *already correct* in 82% of them; what differs is bond orders, aromaticity,
charge or hydrogen count. `facmer_divergent` is **100% `ARRANGEMENT_ONLY`** — the one block on the
board whose name turned out to be a measurement rather than a hypothesis.

**Refuted.**
- 🔴 **That bounding could rescue the string-exact trade.** It cannot. The +0.96 points remains
  available and still costs ~1.8 points of the runtime half; it is now a **standing owner decision,
  not a lane**, because no further engineering makes it cheaper.
- 🔴 **That `structural` is one block with one mechanism.** It is three: `DETACHED` 301 (6.02,
  construction), INTACT+BOUNDARY 172 (3.44, **82% perception**), `NO_STRUCTURE` 11 (0.22). v0.4.17
  had been sized to take all of it as construction work and now declines the middle third.
- 🔴 **The charter's own framing of `BOUNDARY`**, which called it "the attachment call is inside the
  tolerance band" — i.e. treated it as a *cause*. Over the molecules that round-trip **perfectly**,
  BOUNDARY is **35.4%** and INTACT is **63.0%**. Neither discriminates; BOUNDARY is the modal state
  of a *passing* molecule.

**What it cost. 🔴 Three instruments produced a plausible, precise, wrong answer before anything
agreed.**

1. **A normalizer wrong 57 times in 109 — a coin flip.** The heavy-atom comparison began as string
   normalization; SMILES ring-closure digits and atom ordering are arbitrary labels, so an
   identical graph written two ways read as different. It printed `SKELETON` at **74/172 (43%)**
   where the canonical comparison reads **4/172 (2.3%)** — **a factor of 18**, and *both* look like
   a finished measurement. **The broken version confirmed the roadmap's existing assumption**, which
   is exactly why it would not have been questioned. Caught only because a **read example
   disagreed with the instrument**; neither alone sufficed, since the eyeball read was too shallow
   on a 300-character macrocycle and the tool was confident.
2. **A derived runtime curve that understated every bounded row**, by ignoring the work after the
   pool loop. The error grew with the bound — worst exactly where the decision sits — and biased
   toward *promoting*.
3. **`pgrep -f` matching its own `bash -c` body**, a trap already in this project's notes, reported
   three live workers for an arm that had been dead for minutes; the same mistake in a different
   form (`ps | grep` matching the watcher's own command line) then cost a confirmation arm, silently
   SIGTERMed at 30 of 48 without writing its output.

**The method result is the transferable one.** A knee curve is **arithmetic, not a parameter
sweep**: the accepted fallback conformer is returned whatever the pool does afterwards, so bounding
at *N* changes the answer only for molecules whose hit lies beyond *N*. Recording each molecule's
minimum viable bound makes both curves derivable from **one** run — the charter had budgeted 1–2
hours *per point*. That is legitimate only because it was checked rather than assumed, at four
points the derivation could not fake: the recovered ceiling **48** against v0.4.15's independently
measured **48** (exact); both runtime endpoints within **5%** of the frozen arm; a live
bound-0 gate reading **0 gains, 0 losses, 0 output moved**; and a live bound-12 arm agreeing
**48 of 48, in both directions**.

⚠ **The first version of that confirmation arm was worthless and looked perfect.** The population
was sorted, so when it was killed at 30 of 48 every completed row happened to be a molecule the
bound *keeps* — 30/30 agreement, confirming nothing about the bound. An arm containing only cases
your mechanism accepts cannot discriminate. The 10 discriminating molecules were re-run to close it.

**Durability, and a self-inflicted error in it.** The release froze the corpus of record as a
per-molecule extract — 5000 rows in 0.26 MB against 268 MB for the directory — and a post-tag
audit then found two larger gaps: the **cohort manifests had never been tracked** (the sample
under every corpus figure this project has published), and only the *current* baseline had an
extract, so v0.4.6-sweep and v0.4.8-honest — the two the entire 10.34-point correction rests
on — were summary-only. Both closed; seven sweeps now frozen and each re-derives its headline
exactly.

🔴 **The extract tool leaked local paths into a public tree, and it was the same mistake twice
in one session.** The cohort manifests were scrubbed by hand *because* `measurements/` is
public — and the tool written minutes later had no scrub. **1421 rows** carried an absolute
source path, not through any field that looks path-shaped but through `error`, because a
Python traceback embeds them. `harvest_measurements.py` has refused on this since v0.4.12; a
tool writing into the same tree by a different door has to apply the same guard, or the guard
only covers the door nobody uses. Fixed at the tool, all extracts regenerated, whole-tree
audit clean — but one file reached the remote first, because a sibling session pushed `main`
mid-session.

**No full sweep, deliberately** — the shipped default did not change, the same reason v0.4.15
skipped its own.

Suite: **1050 OK**.

---

### v0.4.17 — the defect was the fold, and the release starts carrying two numbers

**Status: RELEASED, tag `v0.4.17` (2026-09-19), `pyproject` 0.4.17.**

**What shipped.** **PROMOTED to default-ON:** `OIN_EXACT_DONOR_FOLD` (owner decision 2026-09-18).
The donor fold's candidate set becomes rotations × the slot permutations induced by *real* ligand
automorphisms (`canonical_slots._donor_automorphism_permutations` —
`GetSubstructMatches(useChirality=True)` on a pendant-stripped, class-pinned skeleton), and
`fold_parity.resolve` returns that fold directly: no veto, no mirror encode. `OIN_FOLD_PARITY_VETO`
stays ON and becomes **inert, not wrong** — `OIN_EXACT_DONOR_FOLD=0` brings the bucket fold back
and the veto is load-bearing again at once; the coupling test still pins that path. **12 levers
default-ON, 17 held OFF. Nothing else promoted.** Also added: the fault-attribution census
(`tools/census/`, four instruments and the join), `tools/v0417/`, and `/refreeze-goldens`.

**Accuracy: REAL GAIN — the largest step in the window, and the first measured end-to-end *before*
it was promoted.** Two numbers from here on:

| `results-v0.4.17-sweep`, N = 5000, shipped defaults | v0.4.14 sweep of record | **v0.4.17** | A/B predicted |
|---|---:|---:|---:|
| **self-consistent** — `byte_exact` under the honest score, `E(x) == E(G(E(x)))` | 3858 = 77.16% | **4136 = 82.72%** | 82.62% |
| **VERIFIED** — the string parses back to the input's graph **and** a neutral ruler finds the generated structure graph-isomorphic to it with no stereo element mirrored or changed | 3462 = 69.24% | **3730 = 74.60%** | 74.54% |
| passes the metric cannot verify | 396 | 406 | |
| `byte_exact` · `key_equal` · `structural` · `hard_fail` | 3858 · 365 · 484 · 266 | 4136 · 122 · 479 · 236 | |

*Self-consistent* is what every earlier section of this document calls `byte_exact`. *VERIFIED* is
new, and it is the one to trust: a canonicality lever can raise the first by teaching the second
encode to agree with a wrong structure, and only the second says it did not. Here both moved
together, **+5.56 / +5.36 pts**.

**The lever, isolated.** The real sweep harness over the lever's *complete* mover set — 504
molecules, derived from coordinates: 392 whose input string moves ∪ 286 whose generated structure's
re-encode moves; every other molecule gets the same two strings in both arms and is unchanged by
construction — OFF and ON arms run *simultaneously* so load was common to both: self-consistent
**+273 (276 gains / 3 losses)**, VERIFIED **+265 (269 / 4)**; dead-lever check 392; noise floor
3/504. The sweep then reproduced the ON arm's verdict on **504 / 504** movers and left 4,494 of
4,496 non-movers unchanged, so +278 = +273 (the lever) + 3 (the A/B's own noise-floor flips) + 2
(non-mover timeouts that finished). **Book the lever at +5.46 / +5.30; the rest is run-to-run.**
`slot_renumber` **252 → 5**. All four losses are generator-side label dependence — an equally valid,
differently *labelled* input and a worse structure — the dependence v0.4.14 measured at 7 in 182,
still open.

**The census came first, and it re-pointed the release before a line of it was written.**
Owner-approved: every one of the 5,000 molecules of the v0.4.14 sweep given exactly one fault,
`UNATTRIBUTED 0`, by four instruments none of which routes through the encoder to judge the
generator or through the generator to judge the encoder. **22.84 pts = G 10.54 · E1 8.54 · E2 2.96
· P 0.60 · data 0.20.** It moved `slot_renumber` — 252 molecules / 5.04 pts, filed for three
releases as enantiomers the generator could not build — onto the **encoder**: the molecules are
achiral by the ruler and `E(mirror x) ≠ E(x)`. `G_HANDEDNESS` is **9 molecules**. And the string
is wrong for **59% of `hard_fail`** (158 of 266 carry a string that does not describe the input),
so the "compute floor" was 2.08 pts, not 5.32. It also confirmed v0.4.16's perception estimate by
an independent path: `E2_P_*` 148 / 2.96 against the 141 / 2.82 estimated.

**Refuted.**
- 🔴 **The handoff's own prescription — an achirality test for the parity veto.** The veto was not
  the defect. The bucket fold permutes each symmetry class of a fragment *independently*; a
  ligand's C2 axis moves all of its classes at once, so a one-class swap is not a symmetry of the
  ligand — it relabels the complex as a different arrangement (the mirror image, for a cis-α
  tetradentate). That over-fold is why v0.4.11 collapsed enantiomers, why v0.4.12 built a veto,
  and why the veto fired on one hand of 568 achiral molecules. Folding over true automorphisms
  needs no veto and no mirror encode: 550 of the 568 pairs unified, positive control 5,792 / 5,792,
  offline == live on 2,443 / 2,443, key-level drift 255 → 255 with membership changed on 0 — the
  lever touches labels and nothing else.
- 🔴 **The three "oracle-confirmed enantiomer" fixtures the veto was built on (v0.4.12) and
  promoted against (v0.4.13) were a coin.** Under the *shipped* encoder a renumbered `x` already
  landed on its **mirror's** string 8, 7 and 6 times in 12 for `BIWDIV`, `CIHVAT`, `OJEKET`.
  `TestPresentationInvariance` drew one renumbering per fixture — one toss, three times lucky. The
  oracle that confirmed them chiral conflates conformation with configuration; the census ruler
  calls all three achiral. This retires the evidence v0.4.13 was promoted on: v0.4.11's "221 of
  393 collapse enantiomers" and the veto's "19 → 0" were reading this coin. Replaced by a pinned
  refutation and an exact-fold property test.
- 🔴 **"v0.4.17 is construction"** — every earlier section of this document says so. The census
  made it an encoder-canonicality release and sent the construction block to v0.4.18.
- 🔴 **"`hard_fail` is the compute floor."** 65 shipped strings described an arrangement the input
  is not — the bucket fold's lex-min can land outside the input's own orbit — often one that cannot
  be built. 19 of them were 300 s timeouts and 16 now run (`ROLYIB` 300 s → 3 s). The generator
  had been searching for the impossible until it was killed.

**Speed: a real default-path win, unclaimed until measured.** The veto's two extra encodes were
being paid *per conformer* — `accept_fn` re-encodes every conformer in the pool. Over the 504
movers: Σ `elapsed_s` 6.19 → 2.93 h, `> 30 s` 110 → 58. Corpus-wide, like-for-like with the v0.4.14
sweep: Σ **38.7 → 33.2 h**, median **4.01 → 3.32 s**, `> 30 s` **678 → 579 (11.58%)**, `> 300 s`
214 → 180, max 729 → 698 s.

**What it cost.**
- **17 of the 36 pairs the veto "protected" are unified — owner decision: ACCEPT.** Their strings
  differ only in which of two *automorphic* atoms carries which integer (`PEBVEZ` is tmeda written
  forwards and backwards); the chirality is real but lives where the string does not look — a
  metal-bound N–H / C–H centre, the Y1 P3 blind spot. 13 of the 17 already flipped under
  renumbering. Census "chiral & same" is 163, not 146, from here on; the repair is to encode the
  centre.
- **18 pairs the census called achiral stay split.** They are chiral by how a ligand wraps (Hf/Zr
  bis(phenolate)-bis(ether) cis-α, helical pentadentates); the ruler matches donors by class and
  shares the bucket fold's blind spot. **A ruler and the thing it measures can share a blind spot
  without sharing code** — by sharing a method.
- 🔴 **15 of 425 ARM 2 golden rows had been wrong since v0.4.13, and nothing said so.** Earlier
  re-freezes re-ran only *predicted* movers, and ARM 2 runs in full only at a release. The full
  gate run re-froze **51 of 425** (40 moved by this lever, 15 stale — 14 restored by
  `OIN_CANONICAL_DONOR_FOLD=0` alone, which dates them); this release's own mover list would have
  missed **19 of the 51 (37%)**. **A mover list is a sample frame, not a verdict** — hence
  `/refreeze-goldens`, from a FULL run, every promotion. ARM 2's predicate is also **circular**
  (byte identity through the generator's own graph): none of the 51 went `in == out → in ≠ out`
  and 21 went the other way, which says the gate lost nothing, not that accuracy moved.
- **"ARM 2's field 3 is a fresh stochastic generation"** — written into `levers.py` and into this
  document's v0.4.14 section — is **false**. Generation is seeded; the full gate reproduced all 368
  untouched rows six shards wide. Marked, not deleted.
- The gate scores a SIGKILLed row as a string mismatch (`EQEROI`, `MUKGUW` fail on a loaded box
  whatever the golden says); the freeze tool copied a `.gz` unscrubbed and would have kept 1 of 20
  staged files. Both caught before the public tree.
- 9 of 1,059 tests broke at promotion; 7 were the coin above, 2 bookkeeping.

Suite: **1064 OK**.

---

### v0.4.18 — DETACHED was a distance

**Status: RELEASED, tag `v0.4.18` (2026-09-23), `pyproject` 0.4.18.**

**What shipped.** **PROMOTED to default-ON, generator-side, in two steps:**
`OIN_ETA_COVALENT_TARGET` + `OIN_VDW_EXEMPT_BINDING` together (owner decision 2026-09-20), then
`OIN_ETA_RETARGET` (owner delegation 2026-09-23). All three live in `generator3d/`
(`clean_geometry.py::ff_clean`, `__init__.py::_try_accept`); none touches a string the encoder
emits — `smiles_1` moved on **0 of 5,000** in every sweep. **Added and held OFF, measured:**
`OIN_ETA_TARGET_UNSCALED`. **15 levers default-ON, 18 held OFF.** Also: `tools/v0418/`, an honest
observation column on ARM 2, and `tests/unit/test_eta_construction_levers.py`, which pins ferrocene
attached under shipped defaults and detached at `"0"` so the defect stays reproducible.

**Accuracy: REAL GAIN — +4.22 self-consistent / +3.80 VERIFIED, both predicted before they were
swept.**

| N = 5000, shipped defaults | v0.4.17 sweep | after the pair (`results-v0.4.18-sweep`) | **release sweep** |
|---|---:|---:|---:|
| self-consistent | 4136 = 82.72% | 4305 = 86.10% (predicted 86.04) | **4347 = 86.94%** (predicted 87.00) |
| VERIFIED | 3730 = 74.60% | 3883 = 77.66% (predicted 77.60) | **3920 = 78.40%** (predicted 78.42) |
| `structural` · `key_equal` · `hard_fail` | 479 · 122 · 236 | 320 · 129 · 228 | **284 · 119 · 230** |

Each promotion was predicted by a harness A/B over **all 1,146 η-bound molecules** — the ones
that pass today are in the cohort because they are where a loss comes from — and each sweep landed
on its prediction to within the 300 s budget boundary: 7 verdicts differed after the pair and 11
after the retarget, **every one a row at the 300 s budget**. `tools/v0418/sweep_vs_ab.py` checked
the two claims each projection stood on rather than repeating them: **3,713 / 3,713** non-η
structures byte-identical to the previous sweep (a non-η molecule cannot reach any of the three
levers), and **1,068 / 1,068** then **1,070 / 1,070** η structures byte-identical to the A/B's ON
arm (generation is deterministic at this load). The pair's OFF arm reproduced the v0.4.17 sweep
byte-for-byte on 1,063 / 1,063 structures — **a noise floor of zero**, so every transition is the
levers. ⚠ The retarget's A/B used the sweep of record's own rows as its OFF arm, so its "noise
floor" line prints 0 by construction; its levers were shown to fire (158 of 1,067 structures
differ) and its OFF-lever probe matched the sweep 103 / 103.

**The census re-run came first — twice.** The v0.4.17 table was stale for every mover, so before a
lane was chosen it was re-run on the v0.4.17 sweep: **17.28 / 25.40 pts = G 9.62 · E2 3.50 · E1
3.38 · P 0.58 · data 0.20**; noise floor of attribution 14 non-movers = 0.28 pts. E1 had fallen
8.54 → 3.38; the residual non-canonical strings (L1d) were worth **0.24 pts** and were demoted;
**`DETACHED` was 343 molecules / 6.86 pts, 86% of it η-bound**, piano stools (`TET` 163, `TPL`
55) — 39% of everything that failed and the largest block by a factor of 3. v0.4.15 had closed
the *selection* door on it (0 of 289: every conformer in the pool was detached), so the first
instrument was a ruler, not a hypothesis.

**The ruler.** A mapping-free distance audit — ligands matched by formula, contact cutoff = the
encoder's own `S + 0.45 Å`; control: the input against itself reads 0.000 Å — found the generator
puts η carbons a median **+0.73 Å** too far from the metal in the DETACHED class, +0.17 Å even in
verified passes; Fe +0.84, Co +0.78, Ni +0.71. Two mechanisms, both in upstream's `ff_clean`:
(1) an η group's target distance was multiplied by 1.1 / 1.2 ("elongation of haptic interaction")
while real M–C(η) over 1,378 groups is 1.008 × the covalent sum and the encoder's contact cutoff is
~1.20 × — the ring was aimed at the edge of the bonding range; (2) the FF scan's vdW guard perceives
bonds by distance, and for a small metal the bond threshold lies *inside* the clash threshold
(Fe–C 2.70 vs 2.81 Å; Co, Ni, Cu, Pt, Zn, Au), so a ring stepping inward registered five new
"clashes", the step was reverted and the scan stopped at 2.87 Å whatever the target. **Ferrocene
itself was being built detached — Fe–C 2.85 Å against a real 2.05 Å, honest round trip failing —
and no test looked**: 1,067 of 1,067 tests were green when the default flipped.

**The second lever is a selector that judges from the string.** The pair's A/B found its losses
(39 verified: Ru 11 · Zr 8 · Ir 6 · Rh 5, metals with no dead zone) had one mechanism: the pool
sweeps a global scale 0.8 … 1.2, first acceptable wins, and the 0.8 conformer used to carry an η
target of 0.8 × 1.2 = 0.96 S — right *by accident*. `OIN_ETA_TARGET_UNSCALED` (target = S at every
scale) was run as a full third arm and is **dominated** (+179 / −59 verified against +189 / −39) —
but the two configurations break *different* molecules, and a per-molecule oracle over the three
arms is **+222 verified with no losses**. The FF scan had never checked that a binding group
*reached* its target (a stalled ring reads as converged), so `ff_clean` now counts the declared
binding atoms it left outside the contact cutoff, and an η conformer that failed to clean or left a
donor out is cleaned once more from the same embedding with the target at S; the one with more
donors in is kept. One conformer per attempt, so the shipped seed sequence is untouched and a
conformer whose donors arrived is never touched. An offline simulation over three structures per
molecule with known verdicts — one that *can* express a loss — bounded it first at +45 / −1; the
harness A/B measured **+49 / −4 self-consistent, +43 / −5 VERIFIED**, booked at the
budget-excluded row as **+0.80 / +0.70 pts**.

**Speed: a default-path improvement, small, and reported both ways.** Corpus-wide, like-for-like
with the v0.4.17 sweep: Σ `elapsed_s` **33.18 → 31.95 h**, `> 30 s` **579 → 521 (10.42%)**,
`> 300 s` 180 → 173, max 698 → 628 s. ⚠ The median moved the other way, **3.32 → 3.53 s**. On the
η cohort the pair's A/B read Σ 14.03 → 11.77 h and the retarget's read flat (11.63 → 11.81 h); why
the corpus median rose was not measured. Goal B (`max(elapsed_s) < 30 s`) is **not delivered**.

**Refuted.**
- 🔴 **The rigid placer as the fix.** Where it places, it is perfect (`FEPZOQ` +0.02 Å in 1 s), but
  it failed to build 5 of 18 and **broke 5 of 6 controls**.
- 🔴 **`OIN_ETA_TARGET_UNSCALED` as a default.** It attaches more rings (182 against 168 of 295
  intact) and breaks more passes. Recorded as a negative *as a default*; it is one half of the
  retarget's choice.
- 🔴 **"Selection is exhausted" (v0.4.15) as a statement about the block.** It was a statement
  about the *pool*: every conformer in it was detached. With the targets fixed the pool can hold
  attached ones, and a selector then recovered +43 verified. Reachability is a property of a
  mechanism — the rule this document already carries, met again from the other side.
- 🔴 **Surplus contacts as a selection criterion.** Penalising contacts the string does not declare
  cost 20 losses in simulation — BOUNDARY contacts are the modal state of a *passing* molecule
  (v0.4.16).

**What it cost.**
- **39, then 5, verified passes lost and accepted** (the pair: Ru 11 · Zr 8 · Ir 6 · Rh 5; the
  retarget: Ru 3 · Pt 1 · Mo 1). 18 of the 39 and 3 of the 5 still round-trip byte-exact and fail
  only the ruler — which is why the selector must never be the round trip. The retarget brought
  13 of the pair's losses back.
- 🔴 **ARM 2 cannot see what the generator builds.** Its `smiles_2` is
  `get_oin_string(result.mol, coords)` — through the generator's own bond graph, the pre-v0.4.8
  *scored* path. In the full gate run **138 of 198 η rows built a different structure** than under
  v0.4.17 and the gated field moved on **0**; at the second promotion 41 of 198 and 0 again. It
  gates the encoder and the generator's graph + stereo; **a construction regression passes it, and
  did for every release that shipped ferrocene detached.** The release adds an honest observation
  column (sha256 of an independent `XYZToSMILES().convert` of the generated coordinates) beside the
  gated field: of 401 measured rows, **135 of the 378 the gate passes fail the harness's own honest
  predicate — 36%** — the scored-vs-honest gap v0.4.8 measured on the whole cohort, on cohorts
  chosen for being hard. Validated on the 370 rows whose gate structure is byte-identical to the
  release sweep's: the column agrees with the sweep's honest verdict 370 / 370. (The lane doc and the
  golden comment block first quoted 134 of 377 — a transcription error against the table's own
  401-row denominator, corrected at the source on 2026-09-23.) Gating on the column would re-freeze ~135 rows as
  honest fails or drop them — a change to what ARM 2 *means*, left to the owner with the evidence
  frozen. 0 of 425 rows owed at either promotion.
- **Found on the way, parked:** σ bonds in verified passes are a median **0.31 Å too short** — the
  pool's global 0.8 scale wins first. Invisible to both accuracy numbers; visible to anything that
  reads geometry.
- Two full sweeps (6 h 24 min + 6 h 29 min wall), four η-cohort A/B arms, two full ARM 2 gate
  runs. The carry-forward licence was voided twice in one release; both debts paid.

**The census of record, on the release sweep.** FAIL 653 = **13.06 pts** · VERIFIED gap
**21.60 pts**, `UNATTRIBUTED 0`, and **the noise floor of attribution is zero** — all 4,157
non-movers keep their fault, because the release sweep is byte-identical to the v0.4.17 sweep
outside the molecules the levers touch. G 481 → 278 (5.56 pts); `DETACHED` 334 → 142; 177 movers
went `G_CONSTRUCTION → NONE` (verified, not merely passing) and 18 went the other way. The three
largest blocks are now within a point of each other: **G construction 3.94 · E2 perception 3.52 ·
E1 serializer 2.96.**

Suite: **1076 OK**.

---

## What is not known

Stated explicitly, because the gaps are as decision-relevant as the numbers.

1. ~~**Corpus speed since v0.4.6 is unmeasured.**~~ **CLOSED by this refresh's sweep:** `> 30 s`
   **678/5000 = 13.56%**, median **4.01 s**, max **728.8 s**. ⚠ Not strictly like-for-like — BLAS
   threads capped to 1 here, not in the v0.4.6 run. **What is still unknown is attribution of
   the first step:** no release claimed the v0.4.6 → v0.4.14 win, so a −6.3-point improvement has
   no owner, and bisecting it across v0.4.7–v0.4.14 has not been attempted. The two steps since
   are owned: v0.4.17 (678 → 579, the veto's per-conformer encodes) and v0.4.18 (579 → 521,
   attached rings).
2. ~~**Corpus accuracy rests on the v0.4.6 generator run.**~~ **CLOSED: 77.16%, measured.** And
   v0.4.13's argument that the re-score was exact is **refuted** — key-invariance bounds
   *acceptance*, not *embedding*, because the generator consumes the OIN string. The re-score chain
   over-stated the deltas *and* under-stated the base; they cancelled to within 0.14 points, which
   is luck. **The standing rule replacing it:** any lever that relabels slots must be A/B'd over the
   coordinate-derived affected population, never re-scored over frozen structures.
3. **`OIN_MEMO_CIP_REPARSE`'s promotion gate has not been run** — the full 328-molecule cohort with
   the lever on, ~10 CPU-h sharded 6-way. Until it is, the project's largest measured single speed
   win stays off.
4. **The five bucket-report cohorts are still not comparable**, so there is no honest like-for-like
   accuracy line from v0.4.4 to today. The v0.4.6 and v0.4.14 sweeps share a cohort and seed but not
   a code version or a thread configuration — closer than anything before, still not a time series.
   **Partly CLOSED from v0.4.14 on:** the v0.4.14, v0.4.17, v0.4.18 and v0.4.18-release sweeps
   are one configuration and differ only in the code, so 77.16 → 82.72 → 86.10 → 86.94 *is* a
   series. Nothing before v0.4.14 joins it.
5. **`OIN_ETA_ACCEPT_EXIT` has no runtime claim at all** — its A/B was stopped, not banked.
6. ~~**`PREFILTER_VETO`'s corpus prevalence is n = 1.**~~ **CLOSED by v0.4.14: 4/261 vetoes = 1.5%,
   2/44 molecules = 4.5%** (stratified sample, seed 7, 5 `INSTRUMENT_DEAD` excluded). AROHIA is an
   outlier. Its latency objection is **refuted** — measured *faster*. **What is still unknown is the
   recovery count**, because the only one measured came from a stochastic A/B; a non-A/B count needs
   `probe_accept_gap.py` over the 26 molecules where the prefilter actually vetoes.
7. ~~**The `byte_exact` molecules that read `DETACHED` are unexplained.**~~ **CLOSED by v0.4.15:
   the second explanation was right.** Of the 52 on the fresh sweep, **51 hold every claimed
   coordination site** under the guard's own predicate — they lose 1–3 *light or ambiguous* donors
   (H, Si, B, F) with the actual count often unchanged or higher. `coordination_report` asks "did
   the donor set change"; `ligands_attached` asks "did a site go empty". **Two different tests, and
   the bucket's verdict is not the guard's verdict.**
8. ~~**The MEDZUR class still has no mechanism.**~~ **CLOSED by v0.4.16, and the answer re-points
   the ladder.** All 187 classified, 0 unaccounted, every class read against an example. Over the
   172 `structural` INTACT+BOUNDARY: **PERCEPTION 141 (82.0%) · CONSTRUCTION 25 (14.5%) · STEREO 6
   (3.5%)** — the heavy-atom graph is *already correct* in 82% and what differs is bond orders,
   aromaticity, charge or hydrogen count. `facmer_divergent` is **100% `ARRANGEMENT_ONLY`**, i.e.
   its bucket name was a measurement rather than a hypothesis — the only one on the board that was.
   ~~⚠ **What is still unknown is whose fault the perception is**: a mis-assignment on a faithful
   geometry (fixable by perception, up to 2.82 pts) and a correct assignment on a distorted one
   (construction, ~0 by that route) are indistinguishable in these strings.~~ **CLOSED by the
   v0.4.17 census, by construction of its rule order:** a molecule is filed `E2_P_*` only after
   the neutral ruler has passed its generated structure and chirality, so the perception block is
   perception *on a correct structure* by definition — 148 / 2.96 pts at v0.4.14, **176 / 3.52 on
   the v0.4.18 release sweep** (`E2_P_FRAGILE` 105: strings that change under renumbering or
   0.02 Å noise, 92 of them non-η, Ni 24 / Zn 22; `E2_P_OTHER` 71). It is an encoder lane, and the
   census of record recommends it, or the serializer's η-rank swap, for v0.4.19.
9. ~~**The mechanism splits are one cohort behind the bucket table.**~~ **CLOSED by the census of
   record:** every fault is attributed on the release sweep itself, 5,000 / 5,000, `UNATTRIBUTED 0`,
   with a control that reproduces the frozen v0.4.14 table byte-for-byte and a non-mover change
   count of **0** as the noise floor. There is no longer an older corpus behind any figure in the
   gap table.
10. **Why runtime improved from v0.4.6 to v0.4.14 is unexplained.** Median more than halved with
    no release claiming it. Candidates include v0.4.10's default-ON deletion finally being seen at
    corpus scale and the thread caps, but nothing separates them. The v0.4.17 and v0.4.18 steps
    are explained (item 1). ⚠ **New:** v0.4.18's corpus median *rose* 3.32 → 3.53 s while its
    tail fell (`> 30 s` 579 → 521); not measured.
11. ~~**Whether `OIN_ACCEPT_STRING_EXACT`'s 4.00× cost can be bought down is unmeasured.**~~
    **CLOSED by v0.4.16: it cannot.** The knee experiment was run over all 365 molecules and there
    is no favourable bound — keeping 79% of the +48 costs **89%** of the runtime penalty, because
    the *early* extra conformers are the expensive ones (each is a full embed plus a full
    re-perception), not a wasted tail. **The +0.96 points is still available, still costs ~1.8
    points of the speed goal, and is now a standing owner decision rather than an open engineering
    question.**
12. ~~**Whether the generator can build the correct enantiomer *at all* is unmeasured.** v0.4.15
    showed the pool holds exactly one key-matching conformer and it is the mirror — but it also
    found **`TAYDUV_comp_0`, 1 of 201**. So the correct handedness is **rare, not categorically
    absent**, and "the generator cannot build it" is not yet established.~~ **CLOSED by the v0.4.17
    census — the question was mis-sized by a factor of ~22.** The 201 "enantiomers" were the
    encoder writing achiral molecules two ways; `G_HANDEDNESS` — the generator building the mirror
    of a molecule whose string *does* separate the hands — is **9 molecules** at v0.4.14 and **4**
    on the release sweep. Forcing the pool wide was never run and no longer needs to be.
13. **`SPHERE_DIFF` — 183 verified-gap molecules, 133 of them passes — is G or the ruler's
    tolerance, and nobody knows which.** The ruler disagrees about the coordination sphere on a
    structure the round trip accepts; the C1 squeeze control puts the ruler's own floor at ~1.16%.
    It is the largest single block G owns on the VERIFIED side and has been carried as "a question
    for the L2 lane" since the census. Still open after v0.4.18.
14. **The generator's label dependence has no lane.** An equally valid, differently *labelled*
    input string and a worse structure: 7 in 182 at v0.4.14, all four of v0.4.17's losses. Every
    canonicality lever pays it until the CoordMap stops caring which of two automorphic donors is
    called `{2}`. Named in three releases, sized in two, owned by none.
15. **σ bonds in verified passes are a median 0.31 Å too short**, because the pool's global 0.8
    scale wins first and every σ bond is 0.8 × its length. Invisible to both accuracy numbers (the
    atoms stay in band); visible to anything that reads geometry. Parked with the owner's
    generation lane; **no accuracy figure in this document says anything about geometric
    fidelity beyond the contact cutoff.**
16. **What ARM 2 should mean.** Its gated field cannot see a construction regression (138 of 198 η
    rows changed, field moved on 0), and its new honest column says **135 of the 378 rows it
    passes fail the harness's own predicate**. Gating on the column re-freezes ~135 rows as honest fails or drops
    them; not gating leaves the release gate blind to the class of defect v0.4.18 just repaired.
    Owner decision, evidence frozen. Also owner's: the gate scores a SIGKILLed row as a string
    mismatch (`EQEROI`, `MUKGUW`).
17. **The residual `DETACHED` is two populations and only one has a mechanism.** η-bound 103
    (Ru 37, Ti 12, Zr 11) — mostly the ~30 molecules of oracle headroom where a per-atom arrival
    count ties or points the wrong way, plus rings the retarget did not reach; and **non-η 48**
    (Y 7, Pt 6, Ir 4, Pd 4, Zn 4), a σ-donor that never arrives, which no v0.4.18 lever touched
    and no ruler has yet measured.
18. **Stereo the string does not carry: 154 passes / 3.08 verified pts** (`E1_NOT_ENCODED` 118 +
    `E1_NONINJECTIVE` 36 on the release sweep). Invisible to the self-consistent number by
    construction; the Y1 P3 blind spot (a metal-bound N–H / C–H centre with its tag cleared) lives
    here, and v0.4.17's accepted 17 unified pairs moved into it. The repair is to encode the
    centre; no release owns it.

---

## Where the remaining 13.06 / 21.60 points are

**From the census of record on the v0.4.18 release sweep**
(`measurements/v0.4.18-release-sweep/v0418_rcensus_*`, `docs/agentic-notes/v0.4.18/CENSUS_OF_THE_RELEASE_SWEEP.md`),
N = 5000, honest scoring, every molecule given exactly one fault by four instruments none of which
routes through the encoder to judge the generator or the reverse, `UNATTRIBUTED 0`. Supersedes the
22.84-point bucket table this section carried from v0.4.14 through v0.4.16 — that table partitioned
by *bucket name*, this one by *measured fault*, and the census's first finding was that two of the
old table's largest blocks were under the wrong owner.

**The FAIL side — 653 molecules = 13.06 pts** (v0.4.17: 864 / 17.28):

| owner | fault | v0.4.17 | **v0.4.18** | pts | what it is |
|---|---|---:|---:|---:|---|
| **G** generator | `G_CONSTRUCTION` | 386 | **197** | **3.94** | generated graph ≠ input: **`DETACHED` 142** (η 103 / non-η 48), `SPHERE_DIFF` 50, `LIGAND_DIFF` 5 |
| | `G_NOTHING` | 76 | 69 | 1.38 | no structure: `TIMEOUT` 61, `NOCONF` 8 — the compute floor |
| | `G_HCOUNT` · `G_DIASTEREOMER` · `G_HANDEDNESS` | 19 | 12 | 0.24 | 3 · 5 · 4 |
| | **G total** | **481** | **278** | **5.56** | |
| **E2** perception | `E2_P_FRAGILE` | 101 | 105 | 2.10 | the input's own string changes under renumbering / 0.02 Å noise; 92 non-η, Ni 24 / Zn 22 |
| | `E2_P_OTHER` | 74 | 71 | 1.42 | bond orders / labels on a structure the ruler passes |
| | **E2 total** | **175** | **176** | **3.52** | |
| **E1** encoder | `E1_HCOUNT` | 79 | 77 | 1.54 | the string's H count ≠ the input's; 72 of 84 sit in `hard_fail` |
| | `E1_GRAPH` | 77 | 71 | 1.42 | the string's graph ≠ the input's: `LIGAND_DIFF` 39, `SPHERE_DIFF` 32 |
| | `E1_NONCANONICAL` | 12 | 11 | 0.22 | one molecule, two strings — was 266 / 5.32 before v0.4.17 |
| | `E1_NONINJECTIVE` | 1 | 1 | 0.02 | |
| | **E1 total** | **169** | **160** | **3.20** | |
| **P** perception on the input | `P_DETACHED` 18 · `P_E1_COVERAGE` 11 | 29 | 29 | 0.58 | |
| **data** | `DATA_MULTI` | 10 | 10 | 0.20 | two molecules in one file |
| | **sum** | **864** | **653** | **13.06** ✓ | |

**The VERIFIED side adds 427 passes that carry a fault** (v0.4.17: 406) — 1,080 molecules =
**21.60 pts**: `G_CONSTRUCTION` 205 (`SPHERE_DIFF` 133, `LIGAND_DIFF` 63, `DETACHED` 9),
`E1_NOT_ENCODED` 118, `E1_GRAPH` 61, `E1_NONINJECTIVE` 36, `E1_HCOUNT` 7. **By owner, verified:
G 483 / 9.66 · E1 382 / 7.64 · E2 176 / 3.52 · P 29 / 0.58 · data 10 / 0.20.**

**The totals that matter:**

| | FAIL pts | VERIFIED pts | |
|---|---:|---:|---|
| **construction** — the generator builds the wrong graph | 3.94 | 9.66 (with `G_NOTHING` and the rest of G) | `DETACHED` 142 is two populations (item 17); `SPHERE_DIFF` 183 verified may be the ruler's tolerance (item 13) |
| **perception** — a correct structure re-read differently | 3.52 | 3.52 | encoder lane; `E2_P_FRAGILE` is an encoder *stability* defect — one molecule, two hashes |
| **serializer** — the string was wrong and the generator built it | 2.96 | 4.32 | `E1_HCOUNT` + `E1_GRAPH`; 129 of the 148 sit in `hard_fail`. Mapped: three defects, one in string-writing code |
| **stereo the string does not carry** | — | 3.08 | `E1_NOT_ENCODED` 118 + `E1_NONINJECTIVE` 36; Y1 P3 |
| **produce nothing** — the compute floor, with encode timeouts | 1.60 | 1.60 | `G_NOTHING` 69 + `P_E1_COVERAGE` 11 = 80 molecules, the ceiling of an always-works fallback |
| canonicality residue | 0.22 | 0.22 | the residual of v0.4.17's 102; a canonicality defect, not an accuracy lane |
| other generator (H count, diastereomer, handedness) + non-injective | 0.26 | 0.26 | `G_HANDEDNESS` is **4** |
| data + input perception | 0.56 | 0.56 | |
| **sum** | **13.06** ✓ | | VERIFIED rows overlap the FAIL rows and do not sum |

⚪ **The three largest blocks are within a point of each other** — construction 3.94, perception
3.52, serializer 2.96 — which has not been true before. Through v0.4.16 this document said **"15.36
of 22.84 points (67%) require the generator to BUILD something different"**; that sentence was
the most consequential number here, and the census cut it in half twice. First by re-filing 5.04
pts of "enantiomers the generator cannot build" onto the encoder (v0.4.17), then by finding the
construction block was a distance and closing 6.86 → 2.84 pts of it in one release (v0.4.18). What
G still owns on the FAIL side is **5.56 pts**, and 1.38 of that produces nothing at all.

🔴 **The old table's two largest blocks were under the wrong owner, and the ownership column had
been wrong in both directions before.** `slot_renumber` 252 / 5.04 was filed as construction for
three releases; it was the encoder. `structural`/`DETACHED` 301 / 6.02 was filed as "proven NOT
reachable by selection" — true, and the reason was that every conformer in the pool was detached,
which is a construction defect, not a limit. This document has now recorded the same lesson at
v0.4.14 (`rdkit_canonical`), v0.4.15 (both lanes), v0.4.17 (`slot_renumber`) and v0.4.18
(`DETACHED`): **a block's name is a hypothesis, and reachability is a property of a mechanism.** The
census exists so that the next lane is chosen from a measured fault rather than from a bucket.

**The ladder, as the census of record leaves it.** Both v0.4.18 promotions moved what the
generator builds, so the v0.4.17 table was stale for 843 molecules and was re-run on the release
sweep before a lane was chosen; the noise floor of that attribution is **zero** (4,157 / 4,157
non-movers keep their fault). The recommendation to the owner for **v0.4.19** is an *encoder* lane
— E2 perception fragility (105 / 2.10, non-η, Ni/Zn) or the serializer's η-rank swap
(`spec/handoffs/v0.4.19/SERIALIZER_LANE_MAP.md`) — because neither needs a generator sweep to
measure; the residual `DETACHED` needs a different selector than the arrival count and is smaller
per molecule of effort now. The rule from the v0.4.17 re-attribution stands: **do not choose a
lane smaller than ~0.3 pts from this table**, the size of the attribution noise it measured
(0.28 pts) before the release sweep brought it to zero.

**The other structural fact, refreshed.** The two goals are one goal, but less so than they were.
On the v0.4.6 sweep 78.8% of failures never tested the notation; on the release sweep the
molecules that produce nothing are **69 of 653 failures (10.6%)**, and the census showed the
string was wrong for **59% of `hard_fail`** at v0.4.14 — the generator working on the wrong
molecule, not a throughput problem. A 30 s cap would still cost passes: 521 molecules are over it
and 4,347 pass, so at most **~10.4%** of the cohort is on the wrong side of Goal B, against 19.88%
at v0.4.6. The v0.4.6-era figures (37.84 CPU-h recoverable, 251 passes / 5.02 pts at risk) are
retired; the equivalent numbers on the release sweep have not been derived.

---

## Sources

| what | where |
|---|---|
| Per-release narrative and figures | `CHANGELOG.md` §§ [0.4.4]–[0.4.18] |
| **v0.4.18 — the release sweep of record, its census, and the ARM 2 full run with the honest column** | **`measurements/v0.4.18-release-sweep/`** (28 files: `v0418_rsweep_RUN.md`, the honest bucket report + ruler verdicts + parse-back that re-derive 4347 / 3920, `sweep_vs_ab`, the census table and its byte-for-byte control on the v0.4.14 record, the reattribution diff, and every ARM 2 row with field 9) |
| **v0.4.18 — the first promotion's sweep (86.10 / 77.66) and its ARM 2 full run** | **`measurements/v0.4.18-sweep/`** |
| **v0.4.18 — L2: the distance audit, every probe, and all three η A/B arms** | **`measurements/v0.4.18-l2/`** (`v0418_l2_eta_distance_audit`, `eta_ab_report{,_on3}`, per-molecule rows of off / on / on3) · `docs/agentic-notes/v0.4.18/L2_ETA_DETACHED.md` |
| **v0.4.18 — the selection lane: the offline bound and the retarget arm** | **`measurements/v0.4.18-selection/`** · `docs/agentic-notes/v0.4.18/SELECTION_ETA_RETARGET.md` |
| **v0.4.18 — the census re-run on the v0.4.17 sweep (17.28 / 25.40), before a lane was chosen** | **`measurements/v0.4.18-census/`** · `docs/agentic-notes/v0.4.18/CENSUS_OF_THE_v0.4.17_SWEEP.md` |
| **v0.4.18 — the census of record on the release sweep (13.06 / 21.60)** | `docs/agentic-notes/v0.4.18/CENSUS_OF_THE_RELEASE_SWEEP.md` |
| **v0.4.17 — the sweep of record (82.72 / 74.60), frozen** | **`measurements/v0.4.17-sweep/`** (`v0417_sweep_RUN.md`, `bucket_report_honest.md`, ruler verdicts, parse-back, `two_numbers.txt`, the per-molecule extract) |
| **v0.4.17 — the exact-fold lane: the offline audit, the generator A/B, the ARM 2 re-freeze rows and the stale-cause trace** | **`measurements/v0.4.17/`** (23 files) · `docs/agentic-notes/v0.4.17/L1_EXACT_DONOR_FOLD.md`, `ARM2_REFREEZE.md` |
| **v0.4.17 — the census: every molecule of the v0.4.14 sweep given one fault** | **`measurements/v0.4.17-census/`** · `docs/agentic-notes/v0.4.17/CENSUS.md`, `CENSUS_C1_generator_vs_input.md`, `CENSUS_C2_encoder_selfconsistency.md`, `CENSUS_C3_string_sufficiency.md` |
| **v0.4.16 — the priced trade, with the full recovered-vs-bound curve** | **`docs/agentic-notes/v0.4.16/LANE-price-the-string-exact-trade.md`** |
| **v0.4.16 — the 187 characterised, with a read example per class** | **`docs/agentic-notes/v0.4.16/LANE-characterise-the-unmeasured.md`** |
| **v0.4.16 — the two method results, and the three instruments that lied first** | **`docs/agentic-notes/v0.4.16/METHOD_one_run_beats_a_parameter_sweep_v0.4.16.md`** |
| **v0.4.16's frozen evidence** | **`measurements/v0.4.16/`** (12 files: the per-molecule `min_bound` ordinals behind the curve, the 48/48 live confirmation with each row's provenance, the classification with every class's membership, and every population's list) |
| **Every sweep in this document, per molecule** | **`measurements/*/sweep_extract_*.jsonl.gz`** — ten sweeps, **66,769 molecules in 3.1 MB**. The seven through v0.4.14 each verified to re-derive their published headline exactly; the three since (v0.4.17, v0.4.18, v0.4.18-release) re-derive every runtime row in this document (checked for this refresh), and the two v0.4.18 headlines re-derive from the frozen ruler verdicts + parse-back + honest bucket report beside them (v0.4.17's honest bucket report is frozen as a table, not JSON, so its 4136 / 3730 is read from its `two_numbers.txt`, not re-derived). Covers all five cohorts in the figure above plus v0.4.0 (25,197) and v0.4.8-honest. ⚠ No geometries: a mirror audit, clash re-score or distance audit still needs the original directory |
| **The sample definitions** | **`measurements/cohorts/`** — seed, N, dedup priority and every molecule name for all four cohorts. Untracked until v0.4.16, and 1033 basenames exist in both dataset subdirs, so the dedup priority is load-bearing: a cohort rebuilt with a different one is a different cohort that looks identical |
| **The measured reachability map — the partition this document's gap table is built on** | **`docs/agentic-notes/v0.4.15/REACHABILITY_MAP_v0.4.15.md`** |
| **v0.4.15's six A/B arms + the five frozen populations** | **`measurements/v0.4.15/`** (13 files: every arm's per-molecule verdict and every population's membership, so each rate is reproducible) |
| v0.4.15 lanes, both refutations, and the void-arm postmortem | `docs/agentic-notes/v0.4.15/LANE-attach-return.md`, `LANE-enantiomer-accept.md`, `BASELINE_SWEEP_CORRECTIONS_v0.4.15.md` |
| **The v0.4.14 baseline sweep — the absolute, and the per-molecule extract that re-derives it** | **`measurements/v0.4.14-sweep/`** (6 files: the authoritative honest table, the scored-vs-honest comparison, `RUN.md`, and `per_molecule_extract.tsv`) |
| **v0.4.14's frozen instruments** | **`measurements/v0.4.14/`** (18 files, incl. the full 182-molecule A/B, v0.4.13's corrected at-risk measurement, and every sample's seed + membership) |
| The generator-neutrality hole, and both corrected headlines | `docs/agentic-notes/v0.4.14/GENERATOR_NEUTRALITY_HAS_A_HOLE_v0.4.14.md` |
| `key_equal` is 50.7% enantiomers | `docs/agentic-notes/v0.4.14/VETO_RESIDUE_OWNERSHIP_v0.4.14.md` |
| `rdkit_canonical` is η-set drift; the residue re-scoped | `docs/agentic-notes/v0.4.14/RESIDUE_RESCOPED_v0.4.14.md` |
| v0.4.14 lanes + predicted-vs-actual | `docs/agentic-notes/v0.4.14/LANE-01-resonance-fold.md`, `LANE-02-prefilter-prevalence.md`, `CLOSEOUT_v0.4.14.md` |
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
