# v0.4.17 L1 — the donor fold swaps look-alikes, not symmetries

**Branch** `research/v0417-encoder-canonicality` (worktree `../oin-v0417`, off `research/census`).
**Lever** `OIN_EXACT_DONOR_FOLD`, **default ON since 2026-09-18** (owner decision). **Status:** built, unit-tested, gated offline
(exact), live over the full cohort, and A/B'd through the generator on its complete mover set
(§ L1b: **+273 self-consistent / +265 verified, 3 / 4 losses**). **Promoted, and swept: 82.72% self-consistent /
74.60% verified** (§ The sweep) — the new baseline of record.

## Headline

The census (C2) located the `E1_NONCANONICAL` defect — 266 failures, 5.32 pts — in
`fold_parity.resolve`'s left conjunct and prescribed a geometric achirality test. The defect is
one level further down, and fixing it there needs **no achirality test, no veto and no mirror
encode**:

> `_donor_swap_permutations` permutes every symmetry class of a fragment **independently**. A
> ligand automorphism moves all of its classes **at once**. A one-class swap is therefore not a
> symmetry of the ligand; the labeling it produces describes a *different arrangement*.

For a linear tetradentate `t1–i1–i2–t2` wrapped cis-α, exchanging the terminal pair alone yields
the mirror image. That single over-fold is why the fold collapsed enantiomers in v0.4.11, why
v0.4.12 needed a veto, and why that veto — asked to tell an enantiomer pair from an achiral
molecule using strings alone — fires on one hand of 568 achiral molecules.

Fold over **true fragment automorphisms** instead
(`canonical_slots._donor_automorphism_permutations`) and every candidate labeling describes the
molecule in hand. Then a chiral complex and its mirror have disjoint candidate sets (nothing to
veto), and an achiral complex's mirror labeling is `g·L·a` for a rotation `g` and an automorphism
`a` — already inside its own candidate set (one string, no test).

| | shipped (bucket fold + veto) | exact fold, no veto |
|---|---|---|
| the 568 achiral pairs a hand is vetoed on: `E(mirror x) == E(x)` | 0 | **550** |
| `E1_NONCANONICAL` failures whose pair is unified | 0 / 266 | **252 / 266** (5.04 pts upper bound) |
| the 36 ruler-chiral pairs the veto protects: still split | 36 | **19 by construction**, 17 unified — see § "The 17" |
| extra encodes on a fold-active molecule | 2 (self-check + mirror) | **0** |

## The instrument, and what a broken one would print

`tools/v0417/autofold_audit.py`. The slot post-pass is a pure function of the inline string and
its candidate set is an orbit, so applying a fold to the census's stored rotation-only strings
(`veto_probe*.jsonl`, arm `fold_off`) gives byte-for-byte what the encoder emits. It calls the
`src/` implementation with the lever forced — the first draft carried its own copy of the
automorphism search and was rewired before any number here was taken.

**Positive control, read first:** the SHIPPED bucket fold applied to `fold_off` must reproduce the
stored `veto_off` string. **5,792 / 5,792** (x and mirror, all four populations). A path that did
not reach the encoder's code, or a stale string, would print less than that.

**The controls that must not move:** `achiral & same` 150/150 unified under both folds;
`chiral & differs` sample 0/150 unified under both. A fold that unified everything, or nothing,
would fail one of them.

## Results (offline, exact)

`E(mirror x) == E(x)`, bucket fold × exact fold:

| population (C2's, by the ruler) | n | both | **bucket only** | exact only | neither |
|---|---:|---:|---:|---:|---:|
| achiral & differs (the defect) | 670 | **550** | **18** | 0 | 102 |
| achiral & same (control) | 150 | 150 | 0 | 0 | 0 |
| chiral & differs (control sample) | 150 | 0 | 0 | 0 | 150 |
| chiral & differs (ALL) | 1,623 | 17 | 19 | 0 | 1,587 |
| — of which the 36 protected pairs | 36 | 17 | 19 | 0 | 0 |

"exact only" is 0 everywhere: the exact candidate set is a subset of the bucket one, as it must be.

By census fault, the 550: `E1_NONCANONICAL` **252**, `NONE` (passes at risk) 195, `G_CONSTRUCTION`
55, `G_NOTHING` 15, `E1_GRAPH` 13, `E1_HCOUNT` 12, other 8. The 266 `E1_NONCANONICAL` failures
split **252 unified · 5 kept split · 9 in the 102**.

## The 18 — a correction to the census, not a miss

Eighteen pairs the ruler calls achiral, the bucket fold unifies, and the exact fold keeps apart.
Read one by one they are **chiral by how a ligand wraps**:

- six Hf/Zr/Y bis(phenolate)–bis(ether) `[OOOO]` catalysts (`VOLGOV`, `VOLHAI`, `VOLHIQ`,
  `VOLPUK`, `YUXREP`, `ZAYWOO`): phenolates trans, ethers cis — the textbook C2-symmetric cis-α
  tetradentate, chiral and used as such. x and mirror differ by the **ether pair alone**;
- helical pentadentates (`GEFGUW` PNNNP, `VIMQAJ` SNNNS, `ULIDIE`), a tetrahedral bis(N,O)
  (`AMAROX`), three boron-capped clathrochelates labelled `OCT` (`FOSNEI`, `HUTCIK`, `HUTCOQ`).

The ruler's sphere test (`neutral_graph.sphere_fit`) matches donors by **symmetry class and
ligand id**, not by automorphism, so inside one ligand it allows exactly the one-class swap the
bucket fold allows. **The ruler and the bucket fold share a blind spot**, which is why C2's
"568/568 agree with the veto off" could not see these. It is not in the ruler's stated list and
should be added: *chirality carried by chelate connectivity within one ligand reads as achiral.*

Consequences for the census table: five `E1_NONCANONICAL` rows (`FEDYUI`, `JOWBOP`, `TIZLAQ`,
`ULIDIE`, `VIMQAJ`) are chiral molecules whose failure is something else (most plausibly
handedness), and C2's 670 "non-canonical" is at most 652. The table is frozen; the correction
is recorded here, not applied.

## The 17 — the open question, stated rather than buried

The owner's gate was "the 36 protected pairs must stay split". Nineteen do, by construction. The
other seventeen are unified, and I am not going to claim that is free. What they are:

- **the two strings name the same molecule by any reading.** `PEBVEZ`:
  `CN{0}(C)CCN{1}(C)C` vs `CN{1}(C)CCN{0}(C)C` — tmeda written forwards and backwards. `JOCGOZ`:
  the two pyrazolyl arms of one tripod. `CUKPUV`: the two amido N of a dimethyl-cyclen. The only
  difference is which of two **automorphic** atoms carries which integer, and that names no
  property of a molecule;
- **the chirality is real but lives where the string does not look**: a metal-bound N–H or C–H
  stereocentre whose tag is cleared (`JOCGOZ`, `LISNUX`, `GIMKIY`, `PEBVEZ`, `SIDNOK/SIDNUQ` — the
  Y1 audit's P3 blind spot), an η² alkene face (`FOJZOU`), or a conformational twist of a sphere
  that is achiral as a configuration (`ODIHOY` bis-nacnac D2d, `YEGDIA` neocuproine, `CUKPUV`);
- **the split was not stable.** With the SHIPPED encoder **13 of the 17 already change string
  under three random renumberings** (C2's own `e_selfconsistency.jsonl`). A label that atom order
  flips was not protecting anything.

So these move census "chiral & same" from 146 to 163. The honest reading is that they were always
in that class and the veto was hiding it behind presentation noise; the fix is to encode the
centre (Y1 P3), not to keep the label.

> **DECIDED by the owner, 2026-09-18: ACCEPT the 17.** They are counted as NON-INJECTIVE
> (146 → 163) from here on, and the repair is the Y1 P3 lane — encode the metal-bound N–H / C–H
> centre — not a veto kept alive for these shapes.

## Side finding — the shipped string can describe a different arrangement

If the bucket fold's lex-min lands on a non-automorphic candidate and the mirror test does not
veto it, the encoder emits a labeling **outside the input's own orbit**: a string for an isomer
the input is not (often a geometrically impossible one). Among the 2,443 audited molecules the
shipped `E(x)` is outside x's orbit for **33** — and they are `hard_fail` **13 times (39% vs a
6.5% base rate)**, `G_NOTHING` 12, `G_CONSTRUCTION` 7. The generator was handed a complex that
cannot be built. The exact fold cannot emit such a string. The census put these under **G**; the
fault is **E1**. The un-audited remainder (~2,550, mostly achiral & same) is covered by the live
run below.

## Live, full cohort

`tools/census/e_selfconsistency.py` re-run unchanged with `OIN_EXACT_DONOR_FOLD=1` — same seeds,
same eleven presentations per molecule as C2, so every column is comparable line by line.
Output `results-v0.4.17-exactfold/e_selfconsistency_exact.jsonl`.

`#DONE 5000`, 6,909 s wall, stderr empty, `INSTRUMENT` 0, `MISSING` 0. Report:
`tools/v0417/exact_fold_live_report.py` (reads both files, encodes nothing).

**Read this row first — it is what licenses the offline audit.** For all 2,443 audited
molecules the string the live encoder emitted equals the audit's prediction, for x **and** for
its mirror: **2,443/2,443 and 2,443/2,443**. Molecules where the live string equals the SHIPPED
one although the audit predicted a change — what a lever that never reached the encoder would
print: **0**.

| | shipped | exact fold | |
|---|---:|---:|---|
| base encoded (of 5,000) | 4,987 | 4,989 | two fewer 300 s timeouts — no veto encodes |
| again / rewrite / rotate not byte-stable | 0 | 0 | determinism floor unchanged |
| **achiral & `E(mirror) ≠ E(x)`** — NON-CANONICAL | **670** | **120** | −550, the audit's number exactly |
| chiral & `E(mirror) == E(x)` — NON-INJECTIVE | 146 | 163 | +17, the pairs of § "The 17" |
| **renumber drift at slot level** (3 renumberings) | **494** | **33** | −461; C2 attributed 464 to the veto |
| renumber drift at KEY level (perception) | 255 | 255 | membership changed on **0** molecules |
| `E1_NONCANONICAL` failures now mirror-stable | 0 / 266 | **252 / 266** | |
| CPU-h for the 55,000 encodes | 17.76 | 13.29 | ⚠ not like-for-like load (10 vs 9 workers + a unit suite); direction only |

The key-level row is the control that matters: a lever that touched perception would move it,
and it moved on no molecule. The 120 that remain are the 102 (§ "Not done") plus the 18.

**Movers — the generator A/B population, derived from coordinates.** `E(x)` changes for
**392 / 4,987**: `key_equal` 166, `byte_exact` 132, `structural` 54, `hard_fail` 40; by census
fault `E1_NONCANONICAL` 158, `NONE` 117, `G_CONSTRUCTION` 48, `G_NOTHING` 25, other 44. List:
`results-v0.4.17-exactfold/exact_fold_movers.txt`. This is the INPUT side only. A molecule whose
input string does not move can still change verdict because the re-encode of its generated
structure does — 252 − 158 = at least 94 `E1_NONCANONICAL` failures are of that kind — so the A/B
population is this list ∪ the generated-side movers (`tools/lever_string_movers.py`).

**What this run does not say.** It is an encoder-only run. It cannot show a generator-side loss
among the 132 passing movers, which now hand the generator a differently-labelled string
(v0.4.14 measured 7 such losses in 182). That is L1b.

## L1b — the generator A/B: +273 self-consistent, +265 verified, 3 and 4 losses

**Population 504, derived from coordinates.** Input-side movers 392 (`E(x)` changes) ∪
generated-side movers 286 (`E` of the stored generated structure changes;
`tools/v0417/generated_side_movers.py`, 4,748 structures × 2 arms, OFF arm == sweep of record on
**4,743/4,743**); 174 are in both. 249 of the 266 `E1_NONCANONICAL` failures are inside. Every
other molecule has the same input string and the same re-encode in both arms, so it is unchanged
**by construction** and the A/B over the 504 is exact, not a sample.

**Instrument.** The sweep harness itself (`tools/test_dataset_roundtrip.py --mol-timeout 300`,
BLAS=1, g-xTB) — not `generator_ab_honest.py`, which keeps one bit per arm and discards the
structure the verified number is computed from. OFF and ON arms ran **at the same time**, three
shards each (`tools/v0417/run_generator_ab.sh`), so whatever the box did it did to both. Both
arms: lane tree, commit `de4a9d2c`, `xtb_available=True`, 504/504 reports, no traceback.

| read in this order | |
|---|---|
| **dead-lever check** — `smiles_1` differs between arms | **392**, exactly the input movers |
| **noise floor** — OFF arm vs sweep of record, pass/fail flips | **3 / 504** (all `hard_fail` → pass: a quieter box) |

| | pass OFF | pass ON | gains | losses | net | headline (n = 5,000) |
|---|---:|---:|---:|---:|---:|---|
| **self-consistent** (`byte_exact`, honest) | 135 | 408 | 276 | **3** | **+273** | 77.16% → **82.62%** (+5.46) |
| **VERIFIED** (census predicate, ruler on the arm's own structures) | 119 | 384 | 269 | **4** | **+265** | 69.24% → **74.54%** (+5.30) |

The two columns move together, which is the point of carrying both: this lever is not teaching
the second encode to agree with a wrong structure. 223 of the 269 verified gains are
`E1_NONCANONICAL`; 225 of the 249 in the cohort now pass (still failing: 18 `structural`,
5 `key_equal`, 1 `hard_fail`).

Bucket transitions: `key_equal → byte_exact` 244 · `structural → byte_exact` 18 ·
**`hard_fail → byte_exact` 13** · `hard_fail → structural` 6 · `encode_fail → byte_exact` 1.
Gains by side: 184 input-moved, 92 generated-side only.

**The losses, each read.** All four are generator-side — an equally valid, differently
*labelled* input string and a worse structure, the dependence v0.4.14 measured (7 in 182) and
which is still open:

| molecule | what the ON arm built |
|---|---|
| `AQENIV` Co(dmgH)₂(py)(R) | both axial ligands detached; re-read as `[Co_SPL]` (`SPHERE_DIFF`) |
| `KOBNUN` Ni bis(thiosemicarbazone) | a short M···N contact re-read as η² `N{1>}N{1}` |
| `QALWEH` Cp\*V(O)(silsesquioxane) | Cp\* re-read η⁴ — one ring carbon lost its slot |
| `VEKWEQ` Os PNP pincer (verified only) | still `byte_exact`; the ruler reads `LIGAND_DIFF` on the structure |

3 of 132 passing input movers = 2.3%. The exact fold did not emit a wrong string for any of them.

**What the offline prediction got wrong, in both directions.** Holding the generated structure
fixed predicted 247 gains and 8 losses. Live: the generated-side-only molecules matched the
prediction on **111 / 112** (exact by construction, one nondeterministic). Of the **8 predicted
losses, 7 did not happen** — handed the new string, the generator built a structure that
re-encodes to it. And 35 molecules predicted to stay failing (or with no stored structure at
all) passed. An offline re-score cannot express a loss *or a recovery* on an input mover; v0.4.14
recorded the first half of that sentence.

**Runtime — the veto was being paid per conformer.** `accept_fn` re-encodes every conformer in
the pool, and each re-encode of a fold-active molecule ran the veto's two extra encodes.

| over the 504 | OFF | ON |
|---|---:|---:|
| Σ `elapsed_s` | 6.19 h | **2.93 h** |
| median | 7.0 s | 3.6 s |
| > 30 s | 110 | **58** |
| > 300 s | 30 | **9** |
| `hard_fail` | 36 | **19** |

**The `hard_fail` recoveries are the side finding, measured over the whole cohort.** The shipped
`E(x)` lies outside the input's own orbit — a labeling for an arrangement the input is not — for
**65 of the 392 input movers** (the 33 of § "Side finding" were the audited subset). In the A/B:

| shipped string | n | OFF pass → ON pass | OFF `hard_fail` → ON `hard_fail` |
|---|---:|---|---|
| **outside** the input's orbit | 65 | 24 → **49** | **19 → 3** |
| inside it (relabeled only) | 327 | 111 → 267 | 15 → 14 |

19 molecules leave `hard_fail` (13 to `byte_exact`, 6 to `structural`); **all 19 were 300 s
generation timeouts in the OFF arm and 16 of the 19 carried an out-of-orbit string.** The
generator was asked for an arrangement that cannot be built and searched until it was killed:
`ROLYIB` 300 s → **3 s**, `TEGVOR` 300 s → 4 s, `ONOXOG` 300 s → 5 s, `DAQDAB` 565 s → 8 s, all
now `byte_exact`. The census filed these as `G_NOTHING`, "the compute floor". For these sixteen
the floor was the encoder asking for the impossible — and in-orbit `hard_fail` barely moves
(15 → 14), so this is the string, not the saved veto time.

**Decided afterwards by the owner: promote** (§ Promotion). The A/B says what the sweep it owes
should read — 82.62 / 74.54, ±0.06.

## Promotion — default-ON, and what the suite said when the default flipped

**Owner decision 2026-09-18: promote.** `OIN_EXACT_DONOR_FOLD` joins `levers._DEFAULT_ON`; its
evidence block sits above that set in the file's own convention. `OIN_FOLD_PARITY_VETO` stays ON
and becomes **inert, not wrong** — `OIN_EXACT_DONOR_FOLD=0` brings the bucket fold and its veto
back together, and the coupling invariant stays pinned for that path.

Flipping the default failed **9 of 1,059** tests. Two were bookkeeping (`held off` → `promoted`).
The other seven were the finding.

**The veto's own three "oracle-confirmed enantiomer" fixtures collapse under the exact fold** —
`BIWDIV`, `CIHVAT`, `OJEKET` (`tests/fixtures/fold_parity/`), the molecules v0.4.12 built the veto
on and v0.4.13 promoted it against. That had to be understood before anything launched:

| fixture | the two "enantiomer" strings differ in | renumber x 12×, SHIPPED encoder | exact fold |
|---|---|---|---|
| `BIWDIV` Co(pincer)₂ | the two amide N of ONE pincer, `{4}`↔`{5}` | lands on the **mirror's** string **8 / 12** | one string, 12 / 12, both hands |
| `CIHVAT` Mo(phen)(CO)₂Br(allyl) | the two N of the phenanthroline, `{0}`↔`{2}` | **7 / 12** | one string, 12 / 12 |
| `OJEKET` Zn(tripod)(ONO₂) | the two thione S of one tripod, `{1}`↔`{2}` | **6 / 12** | one string, 12 / 12 |

The shipped encoder never separated these enantiomers. It emitted one of two strings per
structure, chosen by input atom order, and `TestTheVetoSeparatesConfirmedEnantiomers` compared
**one presentation of each hand**. `TestPresentationInvariance` — the test `fold_parity`'s
docstring cited as "asserted, not assumed" — drew **one** renumbering per fixture: one toss of
what is measured above as a fair coin, three times lucky. The oracle that confirmed them chiral
(`tools/injectivity/oracle.py`) is a whole-molecule rigid superposition whose own docstring says
it "conflates conformation with configuration"; the census ruler, which reads configuration
only, calls all three achiral. What differs between the hands is a twist.

This is the class the owner accepted as "the 17", met again in the project's own gate fixtures.
It is recorded rather than quietly re-pinned because it retires the evidence v0.4.13 was
promoted on: v0.4.11's "collapses enantiomers in 221 of 393 gains" and the veto's "19 → 0" mirror
audit were both reading this coin.

**What changed in the suite, and why each change is not a re-pin:**

- `test_fold_parity.py` — every veto test now pins `OIN_EXACT_DONOR_FOLD="0"` (written, never
  unset: the veto path only exists there). `TestPresentationInvariance` is **replaced** by
  `TestTheVetoPathIsNotPresentationInvariant`, which pins a measured renumbering per fixture that
  sends x onto its mirror's string, and `TestTheExactFoldGivesTheseOneStableString` (both hands
  + four renumberings → one string; outcome `exact_fold`). A test that asserted a false invariant
  and passed on one draw is gone; the refutation is what is pinned.
- `test_exact_donor_fold.py` — `TestLeverIsHeldOff` → `TestLeverIsPromoted` (unset means ON, `=0`
  disables); `TestOffIsByteIdentical` → `TestUnsetIsTheExactFold` (unset ≡ `"1"`; `"0"` still
  selects the bucket fold, shown on a fixture where the two differ).
- `test_chelate_locked_ez.py` — `VOacac2`'s expected string re-pinned. The v0.4.16 string was one
  of **two** the encoder emitted for that one file (10 renumberings: 8 / 2) and differed from
  its mirror's; VO(acac)₂ is C2v and now has one. Each acac still spans two cis basal vertices.
- `fold_parity.py` — the PRESENTATION-INVARIANCE paragraph is kept and marked false, with the
  measurement. No behaviour change.

**Gate ARM 1** (62 fixtures): lever `=0` reproduces the v0.4.16 golden **byte-identically**
(`#DONE 62`), so the OFF path is untouched; lever ON moves **2 of 62** rows, `KAXVOX` (6 / 4 under
renumbering before, 10 / 10 now) and `VOacac2`. Re-frozen, reason recorded inside the golden.
**ARM 2** goldens: 6 of 100 (`v047`) and 30 of 325 (`v049`) rows are *predicted* movers; their
re-freeze needs real generation on a quiet box and was **owed after the sweep**, not run beside it.
**PAID 2026-09-19, from a FULL gate run: 51 of 425 rows, not 36** — 40 moved by this lever, 15 were
already stale (since v0.4.13's donor fold), and the predicted list would have missed 19 of the 51.
See `ARM2_REFREEZE.md`.

**The sweep.** `tools/v0417/launch_sweep.sh` → the project's own `tools/run_sweep.sh`, the
sweep-of-record configuration (same cohort, 6 shards, `--mol-timeout 300`, BLAS=1), **no `OIN_*`
variable set** — the run tests the shipped default, and `run_config.json`'s lever block must be
empty. Launched from the commit that carries this section (`814abff3`); result in § "The sweep". Read with
`tools/v0417/post_sweep.sh` → `sweep_two_numbers.py`, whose predicate was validated **before the
data existed**: over the sweep of record it reproduces the census exactly, 5,000 / 3,858 / 3,462.
Expected: **82.62% self-consistent, 74.54% verified, ± 0.06.**

## The sweep — 82.72% self-consistent, 74.60% verified

`results-v0.4.17-sweep`, n = 5,000, launched from `814abff3`, shipped defaults — `run_config.json`
lever block **empty**, no `OIN_*` in a live shard's environment — sweep-of-record configuration
(6 shards, `--mol-timeout 300`, BLAS=1, g-xTB). All six shards exited 0; `#DONE 5000`; 6 h 42 min.
Read by `tools/v0417/post_sweep.sh`; the reader's control ran first and reproduced the census on
the sweep of record, 5,000 / 3,858 / 3,462.

| | sweep of record (v0.4.14) | **v0.4.17** | A/B predicted |
|---|---:|---:|---:|
| **self-consistent** (`byte_exact`, honest) | 3,858 = 77.16% | **4,136 = 82.72%** | 82.62% |
| **VERIFIED** (string describes the input AND the ruler passes the structure) | 3,462 = 69.24% | **3,730 = 74.60%** | 74.54% |
| passes the metric cannot verify | 396 | 406 | |

**It decomposes exactly, which is the claim.** The A/B said every molecule outside its 504 movers
is unchanged by construction. Measured:

- on the 504 movers the sweep's pass/fail verdict equals the A/B's ON arm on **504 / 504**
  (408 pass in both) — the seeded generator reproduced itself in a different run, on a different
  day's load;
- of the 4,496 non-movers, **4,494 are unchanged**; the other 2 are `hard_fail` timeouts that
  finished this time;
- so self-consistent **+278 = +273 (the lever) + 3 (the A/B's own noise-floor flips against the
  record) + 2 (non-mover timeouts)**; verified +268 = +265 + 3. Losses: 3 and 4, all movers, all
  the generator-side ones of § L1b. **The lever is worth +5.46 / +5.30 points; +0.10 / +0.06 is
  run-to-run**, in the direction of fewer timeouts, and should not be booked as accuracy.
- `smiles_1` equals the exact-fold encoder run's string on **4,988 / 4,988**;
- re-deriving parse-back from the NEW strings moves **1 of 5,000** string-fault verdicts — "a
  slot relabeling cannot change these" was an assumption in § L1b; it is now a measurement.

| bucket | record | v0.4.17 | |
|---|---:|---:|---:|
| `byte_exact` | 3,858 | **4,136** | +278 |
| `key_equal` | 365 | 122 | −243 — **`slot_renumber` 252 → 5**; `rdkit_canonical` 113 → 117 |
| `structural` | 484 | 479 | −5 |
| `hard_fail` | 266 | 236 | −30 |
| `facmer_divergent` | 15 | 16 | +1 |
| `encode_fail` | 12 | 11 | −1 |

| runtime (`metrics.elapsed_s`, nested, a SUM) | record | v0.4.17 |
|---|---:|---:|
| Σ | 38.7 h | **33.2 h** |
| median | 4.01 s | 3.32 s |
| > 30 s | 678 (13.56%) | **579 (11.58%)** |
| > 300 s | 214 | 180 |
| max | 729 s | 698 s |

**This is the new baseline of record.** v0.4.14's carry-forward licence is void; nothing measured
against it may be quoted beside 82.72 / 74.60. The gap is now **17.28 pts self-consistent,
25.40 verified**, and the census's fault partition no longer sums to it. By census fault, among
the molecules that FAILED the sweep of record: `E1_NONCANONICAL` 266 → **41 still failing** (225
now pass; the 41 are 29 `structural`, 10 `key_equal`, 2 other — the 102-pair residual, the 5
mis-attributed wrap-chiral rows, and generator-side failures the label defect had been hiding);
`G_NOTHING` 104 → 86 (18 now pass — the unbuildable strings). 864 failures remain; read against
the OLD table they are `G_CONSTRUCTION` 380, `G_NOTHING` 86, `E1_HCOUNT` 79, `E1_GRAPH` 77,
`E2_P_FRAGILE` 76, `E2_P_OTHER` 62, `E1_NONCANONICAL` 41, other 63 — indicative only, because a
mover's old fault describes a string it no longer has. The next lane should be chosen from a
re-attribution of THIS sweep, not from that table.

⚠ Not frozen. `results-v0.4.17-sweep/` and `results-v0.4.17-exactfold/` are gitignored;
`/freeze-measurements` writes to the PUBLIC `measurements/` tree on `main` and is the owner's call.

## What was built

| file | change |
|---|---|
| `src/oinsmiles/oin/canonical_slots.py` | `_automorphism_query`, `_donor_automorphism_permutations`; `canonical_slot_relabeling` selects it under the lever |
| `src/oinsmiles/oin/fold_parity.py` | `resolve` returns the exact fold with outcome `exact_fold`, no veto, no mirror |
| `src/oinsmiles/oin/levers.py` | `_HELD_OFF["OIN_EXACT_DONOR_FOLD"]` with the evidence and the open question |
| `tests/unit/test_exact_donor_fold.py` | 9 tests on six real cohort pairs, incl. the defect pinned (bucket collapses `VOLGOV`, exact does not) |
| `tools/v0417/autofold_audit.py` | the offline gate |
| `tools/v0417/launch_sweep.sh`, `post_sweep.sh`, `sweep_two_numbers.py` | the promotion's sweep: launcher wrapper, reader, and the two-number predicate (control: 3,858 / 3,462 on the record) |
| `tools/v0417/exact_fold_live_report.py`, `automorphism_extension_check.py` | live-vs-offline report; the pendant-stripping premise |
| `tools/v0417/generated_side_movers.py`, `run_generator_ab.sh`, `post_generator_ab.sh`, `generator_ab_report.py` | L1b: population, the two harness arms, post-processing, the two-number report |

Two implementation facts worth keeping:

- **Radicals.** `CanonicalRankAtoms` ignores radical electrons; the substructure matcher does
  not. `LAMTAX`'s corrole has skeleton classes A~B, C~D and *zero* non-trivial self-matches until
  the `[CH][CH]` radicals are zeroed. Five pairs moved on that one line.
- **Pendant explosion.** Eight t-Bu groups are 6⁸ automorphisms acting identically on the donors.
  Non-donor leaves are stripped to a fixpoint and each surviving atom's full-graph symmetry class
  is pinned as its isotope, so the stripped graph cannot gain a symmetry. Over the cap the
  fragment contributes the identity — a missed fold, the safe direction.

## Not done, on purpose

- **The 102.** Not one mechanism: 63 differ in slot digits only (42 are `SPY` — a cis-dioxo /
  planar-tridentate sphere that is Cs in reality and formally chiral once one oxo is called
  "apical"; a template artifact), 15 fragment order, 11 η winding char, 13 carry inverted `@`
  tags (chiral by the string's own account). The first three groups are where C2's achirality
  bit is genuinely needed — and it must be **automorphism-based**, or it will collapse the 18.
- **The generator's label dependence.** L1b's four losses are all of this kind. It is a
  generator lane (the CoordMap should not care which of two automorphic donors is called `{2}`),
  and every future canonicality lever pays it until it is closed.
- ~~**ARM 2 golden re-freeze** (36 mover rows) and `/freeze-measurements` — after the sweep.~~
  Both DONE (`ARM2_REFREEZE.md`; `measurements/v0.4.17{,-sweep}/`). Left open by the re-freeze: the
  gate scores a SIGKILLed row as a string mismatch (`EQEROI`, `MUKGUW` fail on a loaded box), and
  field 8 of the v0.4.9 golden is a stale observation column.

## Traps this lane added

- A ruler and the thing it measures can share a blind spot **without sharing code** — here, by
  sharing a *method* (match by class). Reading the 18 disagreements is what found it; the counts
  alone said "568/568".
- "Must stay split" is only a gate if the split is a property of the molecule. Check it under
  renumbering before defending it.
- `GetSubstructMatches` on a `Mol` query reads radicals, charges and isotopes;
  `CanonicalRankAtoms` does not read radicals. Two "symmetry" notions in one function will
  disagree on exactly the molecules perception mangled.
