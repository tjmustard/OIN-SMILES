# Roadmap — `byte_exact = 100%` and `max(elapsed_s) < 30 s`

**Git-durable counterpart of `spec/handoffs/roadmap-100-100/README.md`**, which is gitignored along
with the rest of `spec/handoffs/`. The handoff tree carries the executable charters; this file
carries the ladder, so it survives a clone. **Where the two disagree, re-measure — do not pick one.**

**Authored:** 2026-07-26, against `main` @ `0c21729c` (197 commits ahead of `origin/main`, 0 pushed).
⚠ `main` moved **three times** while this was being written. Re-read the tip.
**Corpus of record:** `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.6-sweep`, N = 5000, seed-42,
`--mol-timeout 300`, PASS 1 + PASS 2 complete.

---

## Why a ladder instead of one push

Both goals are one goal, and neither is currently measurable.

1. **The headline is ~5.7 points over-stated.** The harness scores with
   `get_oin_string(gen_result.mol, coords)` — the generator's own bond graph — so it asserts bonds
   the geometry does not support. FP **61/633 = 9.6%** overall, **48/171 = 28.1% haptic**; FN 8/302 =
   2.6% the other way. `FIYHUT_comp_0` ships both Cp rings 0.85 Å off the iron (10 bonded carbons →
   0) and scores a pass. *(`docs/agentic-notes/v0.4.6/METRIC_FALSE_POSITIVES.md`)*
2. **78.8% of failures never test the notation.** Of 340 failures, 240 are generator timeouts and 28
   produced nothing. The notation-attributable gap is **57/5000 ≈ 1.1%** — where the notation is
   exercised it is ~98.9% correct. *(`docs/agentic-notes/v0.4.6/SWEEP_v0.4.6_5K.md` §2)*
3. **Speed alone does not buy accuracy honestly.** `OIN_ACCEPT_SCORED` "recovered" 90 of those 340;
   `report["coordination"]` showed **60 (66.7%) had ligands off the metal**, 21 boundary, 9 genuine.
   +90 → 95.0% is really **+9 → 93.4%**. *(`docs/agentic-notes/v0.4.6/SWEEP_v0.4.6_5K.md` §3)*
4. **The budget is not a budget.** `OIN3DGenerator(timeout=)` is advisory: 60 s requested,
   60.0–172.8 s spent. `p100 < 30 s` is not currently expressible.

---

## Definitions of done

| goal | metric | today | target |
|---|---|---|---|
| **A — accuracy** | `byte_exact`, scored by **independent** re-perception of the generated XYZ | **77.16%** *(v0.4.14 baseline sweep, n=5000 — the first REAL generator sweep since v0.4.6)* | **100%** |
| **B — speed** | per-molecule wall-clock against an **enforced** budget | **678/5000 = 13.56%** over 30 s; median **4.01 s**; max 728.8 s | **< 30 s, p100** |

⚠ Two live traps in one field. `metrics.elapsed_s` is **nested** — read from the top level it
silently yields `0` — **and it is a SUM** over up to three separately SIGKILLed harness attempts.
The old "max 759.9 s against a 300 s budget" headline was the second trap: all 4658 single-attempt
rows finish within **0.2 s** of their cap. See `v0.4.9/ELAPSED_S_IS_A_SUM_v0.4.9.md`.

## The gap — `100 − 77.16 = 22.84` points (v0.4.14 baseline sweep, n=5000)

⚠ **This table is post-v0.4.14's promotion.** `OIN_RESONANCE_DONOR_FOLD` shipped default-ON for a
net of **+71 molecules / +1.42 points** (78 gains, 7 losses). The previous copy read `100 − 75.88 = 24.12`.

🔴 **The counts below are the OFFLINE ones and are approximately, not exactly, right.** v0.4.14's
offline re-score reported +78/0 losses; measured end-to-end it is ~74 gains against ~6 losses,
because an offline re-score cannot express a loss at all. **The next release owes a sweep** — see
`v0.4.14/GENERATOR_NEUTRALITY_HAS_A_HOLE_v0.4.14.md`. Source:
`measurements/v0.4.14/` · `docs/agentic-notes/v0.4.14/LANE-01-resonance-fold.md`.

🔴 **v0.4.14 also RE-FILED two blocks that had been on the encoder ladder for four releases.**
Neither number moved; who owns them did. See `VETO_RESIDUE_OWNERSHIP_v0.4.14.md` and
`RESIDUE_RESCOPED_v0.4.14.md`.

| block | n | pts | nature | release |
|---|---:|---:|---|---|
| **`structural`** | **484** | **9.68** | 🔴 **NOT ONE BLOCK, AND NOT ONE MECHANISM — split measured across v0.4.15 and v0.4.16.** `DETACHED` **301** (6.02): a return-path guard fires on 10–16 conformers per molecule, promotes a better-attached one 51 times and improves the answer **0** times — attachment necessary, NOT sufficient, so **construction**. `INTACT`+`BOUNDARY` **172** (3.44): 🔴 **82% PERCEPTION** — the heavy-atom graph is already correct and bond orders / aromaticity / charge / H differ (`PERCEPTION` 122 + `GEOM_CODE` 19); only **25 mol / 0.50 pts** is construction and **6 / 0.12** is stereo inversion. `NO_STRUCTURE` **11** (0.22) | `DETACHED` → **v0.4.17 CONSTRUCTION**; the 172 → **a PERCEPTION lane, candidate v0.4.19**; `NO_STRUCTURE` → v0.4.18 |
| `hard_fail` | 319 | 6.38 | compute; **315/319 produce no structure at all** | — |
| `key_equal` → `slot_renumber` | **252** | **5.04** | 🔴 **SPLIT MEASURED v0.4.15: the non-enantiomer part is REACHABLE (28.7%, +48 mol behind `OIN_ACCEPT_STRING_EXACT`); the 201 enantiomers are NOT (0.5%) — the pool holds one key-matching conformer and it is the mirror** | reachable half **available now**; enantiomers → v0.4.16 CONSTRUCTION |
| `key_equal` → `rdkit_canonical` | 114 | 2.28 | 🔴 **NOT canonicality.** 92/114 (80.7%) is **η-set denticity drift** — a perception/geometry problem no string change can reach | **re-filed off the encoder ladder** |
| `facmer_divergent` | 16 | 0.32 | ✅ **CHARACTERISED v0.4.16: 100% `ARRANGEMENT_ONLY`** — metal and ligand bodies agree, only the arrangement signature differs, i.e. exactly what the bucket name says. **"wrong isomer" was right**, which is worth recording: this is the one block whose name was a measurement rather than a hypothesis. 11 of 15 are also `DETACHED`, so v0.4.17 must move them first | v0.4.18 |
| `encode_fail` | 15 | 0.30 | encoder coverage | v0.4.18 (opportunistic) |
| **sum** | **1128** | **22.56** | | |

⚠ **This table sums to 22.56, the headline says ~22.75, and the ~0.19 difference is real rather
than a rounding slip.** The per-block counts come from the OFFLINE re-score, which credited 78
gains and could not express a loss; end-to-end there are ~74 gains and ~6 losses, so ~10 molecules
sit in a bucket this table assigns wrongly. Which buckets they land in is **not known without a
sweep** — that is exactly what the sweep would settle. The table is kept at the offline counts
rather than silently fudged, because a fabricated redistribution would be worse than a stated
inconsistency.

### The re-filing, stated as owners

| | n | pts | owner |
|---|---:|---:|---|
| generator built the enantiomer (`MIRROR_MATCH`, 183/222 of the veto-reverted) | 183 | **3.66** | **generator** — `byte_exact` failing is CORRECT here; recovering it would delete the fact |
| η-set denticity drift (`rdkit_canonical`) | 114 | **2.28** | **perception / generator** |
| unfoldable same-enantiomer labeling (`NOT_A_MIRROR`) | 39 | 0.78 | encoder — needs a **parity-aware** canonicalization |
| resonance residue after v0.4.14 | 25 | 0.50 | encoder (20 slot swaps the skeleton keeps apart, 5 never fire) |

**5.94 points move off the encoder ladder.** `byte_exact` is unchanged by this; the schedule is not.

⚠ **`key_equal` is documented as "benign canonicalization — the win reclaimed". For 183 of its
361 remaining members (50.7%) that is false** — they are enantiomer pairs, invisible because
`compare._parse_vertex_colors` folds reflection deliberately and `accept_fn` decides by that key.
Third time a headline has rested on a metric that folds the axis under test (v0.4.8, v0.4.11,
v0.4.14). **A bucket name that asserts a cause is a hypothesis, not a measurement.**

### The other decomposition v0.4.13 added: by ATTACHMENT, not by bucket

Over the **767 genuine failures** (non-`byte_exact`, excluding `key_equal`'s benign
canonicalization). Control: `byte_exact` reads **1.32%** `DETACHED` against **24.11%** on the
failing side — **18.2× enrichment**, `UNKNOWN` = 0, `#DONE 5000`.

| class | n | what it is |
|---|---:|---|
| **GAVSED** — `DETACHED` | **280** | returned with ligands off the metal; `_select_by_geometry`'s fallback is not attachment-aware |
| **MEDZUR** — `INTACT` | **99** | attachment fine, independent re-perception still disagrees — **still unexplained** |
| `BOUNDARY` | 53 | the attachment call itself is inside the tolerance band |
| `NO_STRUCTURE` | 335 | nothing generated — no returned conformer to be wrong |

Both classes were handed forward on **n = 1** for two releases.

### 🔴 The honest number REORDERED the ladder

`structural` was **9 molecules / 0.18 points** under the scored metric and was parked at v0.4.16 as
*"knowledge, not points"*. Honestly scored it is **417 molecules / 8.34 points — the second-largest
block in the gap**, ahead of `hard_fail`. It is not a rounding change: 350 of them came directly
from `byte_exact`, molecules whose generated structure re-perceives with different coordination
than the OIN claims.

**v0.4.16 is therefore mis-sized and mis-ordered in the ladder below, and the "wrong isomer" label
is doing too much work** — a molecule that loses a haptic contact is not an isomer problem. A
release that owns 8.34 points cannot sit behind two canonicality releases worth 12.20 combined
without someone having decided that deliberately. **Re-sequence before planning v0.4.11.**

#### ⚠ Still not re-sequenced at v0.4.10's close-out — and now one release more urgent

v0.4.10 was byte-identical by construction, so it moved **no** point of the gap and the table above
is unchanged. That means the re-sequencing decision v0.4.9 asked for was simply not taken, and the
ladder still spends its **next** release on its **smallest** block:

| block | pts | currently scheduled |
|---|---:|---|
| `key_equal` → `slot_renumber` | **9.92** | v0.4.14 |
| `structural` | **8.34** | v0.4.16 |
| `hard_fail` | 6.38 | v0.4.9–13 |
| `rdkit_canonical` | 2.28 | v0.4.15 |
| `facmer_divergent` | 0.32 | v0.4.16 |
| **`encode_fail`** | **0.30** | **v0.4.11 — next** |

**18.26 of the 27.54 points sit at v0.4.14 and v0.4.16**, behind a release worth 0.30.

There is a real argument for the current order — v0.4.12 (honest acceptance) and v0.4.13 (harness
false negatives) are *measurement* prerequisites, and moving a 9.92-point release ahead of them
risks attributing a harness artefact to a canonicality fix. That argument covers v0.4.12/13. **It
does not cover v0.4.11**, whose 0.30 points depend on nothing and unblock nothing.

**Recommendation, for the human to accept or reject** — this reorders eight planned releases and
the target definition is the user's call, so it is written down rather than applied:

1. Keep v0.4.12 and v0.4.13 where they are; they are prerequisites, not point-scorers.
2. **Re-point v0.4.11 at `slot_renumber` (9.92 pts)** — encoder-side canonicality, the largest single
   block, and the area where v0.4.5 already built six working levers, so it has the shortest path
   from diagnosis to change.
3. Demote **encode floor R3** to an opportunistic slot. It is genuinely cheap and this session hit
   its cost first-hand — ARM 1 takes ~15 minutes almost entirely because several fixtures route
   through `_resonance_candidates_isolated` at a **120 CPU-s budget per fork** — but 0.30 points is
   a rounding error against a 27.54-point gap, and "it also speeds up our own gate" is a developer
   convenience, not a roadmap argument.
4. Re-label `structural` before scheduling it. 417 molecules whose generated structure re-perceives
   with **different coordination than the OIN claims** is a generator-capability problem sitting
   under an "isomer" heading; it likely belongs with v0.4.18, not beside `facmer_divergent`.

Runtime, from the primary reports: overall n = 5000, median **7.19 s**, `> 30 s` **994 (19.88%)**;
**eta** n = 1146, median 24.08 s, `> 30 s` 528 (46.1%) — **53.1% of the whole tail**; non-eta
n = 3854, median 6.18 s, `> 30 s` 466 (12.1%).

### What the tail actually costs, and what bounding it costs (v0.4.9)

| band | n | CPU-h | % of sweep | honest passes |
|---|---:|---:|---:|---:|
| `< 30 s` | 4006 | 8.67 | 15.8% | 3380 (84.4%) |
| `30–60 s` | 292 | 3.45 | 6.3% | 136 (46.6%) |
| `60–300 s` | 411 | 15.36 | 28.0% | 112 (27.3%) |
| **`≥ 300 s`** | **291** | **27.31** | **49.8%** | **3 (1.0%)** |
| total | 5000 | **54.8** | | 3631\* |

**Half the sweep's compute buys three molecules.** And **93.1% of all honest passes already finish
under 30 s**, so the two goals are far more separable than thesis 3 below assumed:

| per-molecule cap | CPU-h recovered | passes lost | `byte_exact` cost |
|---:|---:|---:|---:|
| 300 s | 3.06 | 3 | 0.06 pts |
| 120 s | 21.33 | 52 | 1.04 pts |
| 60 s | 30.97 | 115 | 2.30 pts |
| **30 s** | **37.84** | **251** | **5.02 pts** |

A 5k sweep under a 30 s bound costs **~17 CPU-h instead of 55**, which makes live sweeping an
affordable instrument again rather than something v0.4.8 had to route around with an offline
re-score.

\* 3631 counts string equality; the frozen report says **3623**. The eight-molecule difference is
the **atom-count gate** — molecules that encode byte-identically to their input while having lost
atoms. A string comparison cannot see them. Use `tools/roundtrip_bucket_report.py`, never an ad-hoc
`honest_class.endswith("->byte")`.

---

## Ordering thesis

1. **Instrument before engine.** The headline goes **down once, on purpose**, inside its own release
   boundary — otherwise it is indistinguishable from a regression, a confound this project has been
   caught by twice.
2. **Bound before optimize.** A p100 target against an unbounded budget is not a target.
3. **Compute is the accuracy work** for the failure side — but only when honestly scored, or it
   manufactures phantom passes at a measured rate of 2 in 3.
   ⚠ **Weakened by v0.4.9's measurement.** More compute is *worth less* than this thesis implies:
   the `≥ 300 s` band spends half the sweep's CPU for a **1.0%** honest pass rate, and 93.1% of all
   passes already land under 30 s. Compute buys the tail very little. Re-read thesis 3 as "compute
   is *some* of the accuracy work, and the cheapest part of it is already bought."
4. **Canonicality is the larger axis and is generator-independent** — 520 molecules / 10.40 points vs
   340 failures / 6.80 points.
5. **Information-adding tokens go last.** P1/P2/P3 convert silent collapses into loud failures;
   promoting them lowers the headline until the generator reproduces them.

---

## The ladder

| ver | theme | headline prediction |
|---|---|---|
| **v0.4.7** | close-out: land the 5 finished swimlanes; tag `v0.4.7`; bump `pyproject` 0.4.6 → 0.4.7 | flat |
| **v0.4.8** | **the honest number** — promote `OIN_INDEP_SCORE` to *the* score; 5k honest sweep; + the atom-count gate | **DOWN ~5.7 pts, planned** |
| **v0.4.9** ✅ | **speed becomes measurable** — `OIN_ENFORCE_BUDGET` (default OFF); 328-molecule stratified benchmark frozen. **Refuted its own premise:** `elapsed_s` is a sum, the harness already enforces to ε ≈ 0.2 s | pass flat ✓; benchmark reproduces to **0.28%** |
| **v0.4.10** ✅ | **cost per attempt**, byte-identical. Deleted the discarded `.index()` scan (**`CAHQEJ` −32.9%**, `FOSNEI` +0.3% nil) and added `OIN_MEMO_CIP_REPARSE` (**`VAFMIA` −86.7%**, `CAHQEJ` −2.4%), both **byte-identical on ARM 1 and ARM 2**. ⚠ **Found its own arbiter broken:** ARM 1 had been non-runnable since `dd51a515` | pass FLAT ✓ by construction; **every speed number bimodal by molecule** |
| **v0.4.11** ✅ | **`slot_renumber`** (496 / 9.92 pts) — built the within-fragment donor fold v0.4.5 specified. 🔴 **REFUTED IT.** It works (+7.86 pts, 393 molecules, one direction) and **collapses enantiomers: 221 of its own 393 gains (56.2%), 19/250 on a uniform draw**, 18 oracle-confirmed chiral. Ships **default OFF**. Lane 1 classified all 496: `diff_occupancy` **0**, and 90 of the 118 `distinct_donors_LOCAL` are frozen-resonance-form artifacts ⇒ **v0.4.14, not here** | **FLAT** — the lever that would raise it is unsafe |
| **v0.4.12** | **reflection parity + honest acceptance** — L1 `OIN_FOLD_PARITY_VETO`, the filter v0.4.11's close-out specified; L2 `OIN_ETA_ACCEPT_EXIT`, the winding criterion moved from selection into `accept_fn`, the only site that can stop pool filling. ⚠ `OIN_ETA_EARLY_EXIT` is **runtime-inert as sited** and its promotion gate is **void, not unrun** | default path **FLAT** (both levers OFF); L1 lever-ON gain and L2 tail move measured separately — see `docs/agentic-notes/v0.4.12/` |
| **v0.4.13** | **the donor fold ships + harness false negatives** — promotes `OIN_CANONICAL_DONOR_FOLD` + `OIN_FOLD_PARITY_VETO` to default-ON (owner's call; voids the carry-forward licence, so a full re-sweep and re-frozen goldens are owed). L1 `PREFILTER_VETO` acceptance-side prevalence; L2 the MEDZUR/GAVSED split — **both classes were n = 1 and are now 99 and 280** | **UP ~3.42 pts** from the promotion (72.46% → ~75.88%); Lane 1/2 are measurements, not headline movers |
| **v0.4.14** ✅ | **the residue, re-scoped.** Promoted `OIN_RESONANCE_DONOR_FOLD` (skeleton-level donor equivalence: acac ketone/enol, carboxylate, sulfonate) — **+1.56 pts, 78 molecules, 0 bad**, mirror audit **0 regressions at 179/179 = 100% coverage**. 🔴 **And REFUTED both halves of its own charter:** `rdkit_canonical` is **80.7% η-set denticity drift**, not canonicality; and **183 of the 222 veto-reverted molecules are the generator building the ENANTIOMER**, so `byte_exact` failing on them is CORRECT. **5.94 pts re-filed off the encoder ladder** | predicted **+3.5–4.1**, actual **+1.56** — the miss IS the finding (114 of the chartered 204 were never reachable) |
| **v0.4.15** ✅ | 🔴 **BOTH LANES MEASURED; ONE REFUTED, ONE RE-AIMED BY ITS OWN DATA.** Baseline re-established by the first real generator sweep since v0.4.6: **77.16%** — the charter's "over-stated base, really ~74.88%" alarm is **REFUTED** (0.14 pts, not 2.4). **L1 `OIN_ATTACH_RETURN`: 0 gains of 289, 51 structures moved, 0 losses, ratio 1.01 — REFUTED as an accuracy lever.** The guard promotes a better-attached conformer and re-perception still disagrees: attachment necessary, NOT sufficient, at scale. **L2 `OIN_ACCEPT_STRING_EXACT`: +48 mol / +0.96 pts, 0 losses — but 28.7% on the non-enantiomer 164 vs 0.5% on the 201 enantiomers (57×).** Everything it recovers is OUTSIDE the class the charter aimed it at. Both ship **default OFF**; L2's default is an open owner decision (4.00× runtime, `>30 s` 30→122 on its population) | predicted **+2 to +7**, actual **0.00 at the shipped default** (+0.96 available behind L2). The miss IS the finding |
| **v0.4.16** | 🔴 **RE-POINTED AGAIN, 2026-07-30 — the row below it already claimed the construction work.** The previous text said "this release should be CONSTRUCTION: make the embed able to produce both handednesses at all", which is **v0.4.17's L1** verbatim. Emitting `|mc:±|` stays refuted (the pool holds exactly ONE key-matching conformer and it is the mirror; 1 of 201 recovered), but the successor is not construction — it is **take the measured point and characterise the unknown**. **L1 resolves `OIN_ACCEPT_STRING_EXACT`**, whose +48 molecules / +0.96 pts and 0 losses are already measured and blocked only on a 4.00× runtime trade nobody had priced; the cost is bounded by `OIN_STRING_EXACT_BOUND`, and the knee is DERIVED from one instrumented run rather than a run per candidate bound. **L2 characterises the 3.74 unmeasured points** — because the project's own rule forbids aiming a lane at a block whose mechanism is a guess, and both later releases were sized against guesses about it | **UP 0.0–0.96** — the ceiling is known before the release starts; the open question is the runtime price |
| **v0.4.17** | 🔴 **RE-POINTED 2026-07-30 — CONSTRUCTION: build what selection cannot find.** The old "structural is no longer this release's subject, v0.4.15 L1 owns it" is superseded: **L1 was REFUTED**, so structural returns as a construction problem. Owns the **10.04 pts v0.4.15 proved unreachable by any selection predicate** — L1 the **201 enantiomers** (4.02; first deliverable is *can the embed produce both handednesses at all*, not a fix — `TAYDUV_comp_0` proves it is rare, not absent), L2 the **301 `DETACHED`** (6.02; 87.2% haptic, whole Cp/arene groups lost, 10–16 pool conformers all detached). ⚠ Cheapest decisive experiment: does `clean_geometry` **relaxation** detach rings that embedded correctly? | honest — **may go DOWN**; a documented limitation is an acceptable outcome |
| **v0.4.18** | 🔴 **RE-POINTED 2026-07-30 — THE FLOOR, stated with evidence.** `hard_fail` 266 (5.32, of which **262 produce NOTHING**) + `NO_STRUCTURE` 11 + `facmer` 15 + `encode_fail` 12, then the encoder ladder's last **1.28 pts** (39 `NOT_A_MIRROR` + 25 resonance residue) and the 57 notation molecules. ⚠ `OIN3DGenerator(timeout=)` is ADVISORY and `elapsed_s` is a SUM — **separate budget-exhaustion from genuine incapacity or the floor is over-stated**. **Product: the honest achievable ceiling** = 100 − produce-nothing floor − construction limit − notation-undecided, each with its evidence | **FLAT or DOWN** — the product is a defensible limit, not points |

### LADDER DECISION 2026-07-30 (v0.4.15) — **MEASURED, awaiting owner sign-off**

```
Baseline  77.16% byte_exact (n=5000, first REAL generator sweep since v0.4.6).
          The v0.4.15 charter's "over-stated base, really ~74.88%" alarm is REFUTED:
          the published 77.30% was over-stated by 0.14 pts, not 2.4. What the offline
          re-score chain drifted was bucket COMPOSITION (structural +67, hard_fail -53),
          not the headline.

Result:   v0.4.15 ships FLAT at 77.16%. Both levers land default-OFF.

          L1 OIN_ATTACH_RETURN        0 gains / 289   REFUTED as an accuracy lever
          L2 OIN_ACCEPT_STRING_EXACT +48 gains / 365  works, default is an owner call

🔴 THE TWO FINDINGS THAT RE-AIM THE LADDER:

  1. SELECTION IS EXHAUSTED FOR structural/DETACHED. L1 fires hard -- 10-16
     winding-matching conformers per molecule, EVERY ONE detached -- promotes a
     better-attached conformer 51 times, and improves the answer ZERO times.
     Attachment is necessary, not sufficient. The 6.02 pts in that bucket are NOT
     reachable by a return-path predicate.

  2. THE ENANTIOMER CLASS IS A CONSTRUCTION PROBLEM, NOT AN ACCEPTANCE ONE.
     pool.accept_incumbent_recorded = 1: exactly ONE key-matching conformer in the
     whole pool, and it is the mirror. L2 recovers 28.7% of the non-enantiomer 164
     and 0.5% of the 201 enantiomers -- 57x apart.

     ⇒ v0.4.16's chartered plan to EMIT |mc:±| is now the wrong next move. Emitting
       makes the generator responsible for reproducing a handedness it does not build,
       turning 201 silent wrong answers into 201 loud failures for no gain. v0.4.16
       should be CONSTRUCTION: make the embed able to produce both handednesses at all.

⚠ AND THE SCOPE CALL IS WHAT SAVED THE RELEASE. Scoped to the chartered 201 enantiomers,
  L2 would have shipped +1 molecule. The owner-accepted widening to all 365 key_equal
  (2026-07-29) is the entire lane. The charter's framing -- "enantiomers are the target,
  slot_renumber is a later lane" -- was measured backwards.

OPEN OWNER DECISION: promote OIN_ACCEPT_STRING_EXACT?
  FOR      +48 molecules = +0.96 pts, ZERO losses across 565 molecules
  AGAINST  4.00x runtime on its population; >30 s 30 -> 122 there, ~678 -> ~770 corpus-wide
  The target is byte_exact 100% AND max(elapsed_s) < 30 s, so this is ~+1 pt of one half
  against ~+1.8 pts of the other. If promoted, a full 5k sweep becomes MANDATORY.
```

### LADDER DECISION 2026-07-28 (v0.4.14) — **ACCEPTED by the project owner**

```
Change:   v0.4.15  was  the 57 notation molecules (knowledge, not points)
                   is   THE GENERATOR OWNS THE GAP -- two lanes:
                          L1  structural / DETACHED      266 mol, up to 5.32 pts
                          L2  generator enantiomer class 183 mol,       3.66 pts
          v0.4.16  information-adding tokens P1/P2/P3   (unchanged theme, moved up)
          v0.4.17  generator capability floor           (structural REMOVED -- L1 owns it)
          v0.4.18  the 57 notation molecules + the encoder ladder's last 1.28 pts

Because: it resolves BOTH standing open decisions at once -- v0.4.13's (should structural's
         266 DETACHED displace v0.4.15?) and v0.4.14's (should the enantiomer class?). Both
         are generator lanes; taking them together is what makes the ordering coherent
         rather than arbitrary. Combined they are 8.98 of the remaining 22.70 points.

         v0.4.14 established the premise: the ENCODER ladder has 1.28 points of reachable
         work left in total. Four releases were scheduled against blocks that are mostly
         generator work. Continuing down the encoder ladder would be optimising the
         smallest remaining block -- the exact mistake v0.4.9/v0.4.10 made.

⚠ ACCEPTED WITH A KNOWN CONFOUND. Two headline movers in one release means a reader cannot
  attribute a move to either lane without reading both apart. v0.4.13 accepted the same
  shape and flagged it. Mitigation is MANDATORY, not optional: each lane behind its own
  default-OFF lever, THREE arms measured (L1-only, L2-only, both), and separate commits so
  a bisect can separate them even if the release doc does not.

Not changed, and why: v0.4.16's P1/P2/P3 keeps its theme -- v0.4.15 L2 *uses* the Delta/Lambda
         descriptor as an acceptance predicate without emitting it, so emitting stays a
         separate, larger change. v0.4.17 keeps the capability floor but LOSES `structural`,
         which is now v0.4.15 L1's.

Accepted 2026-07-28 by the project owner. Recorded by the v0.4.14 close-out session.
```

#### The superseded proposal, kept for the record

```
Change:   v0.4.15  was  the 57 notation molecules (knowledge, not points)
          proposed  is  the GENERATOR ENANTIOMER class -- 183 mol / 3.66 pts

Because:  v0.4.14 measured that 183 of the 222 molecules the parity veto reverts are
          MIRROR_MATCH: the round trip encodes byte-for-byte as the input's own
          mirror (tools/veto_residue_chirality.py, 222/222 classified, 0 excluded).
          They sat in `key_equal/slot_renumber` -- a bucket documented as "benign
          canonicalization" -- because compare._parse_vertex_colors folds reflection
          deliberately. accept_fn decides by that same key, so the GENERATOR accepted
          a mirror-image structure and the harness scored it as a string difference.

          At 3.66 pts it is the largest single mechanism left outside `structural`,
          and unlike `structural` its mechanism is already named.

⚠ NOT AN ENCODER LANE, and that is the whole point of raising it. No canonicalization
  change should ever "recover" these 183 -- doing so would raise byte_exact by
  deleting the fact that the generator produced the wrong enantiomer. That is
  v0.4.11's refuted move with a different lever.

⚠ COMPETES WITH the still-open v0.4.13 decision below (`structural`'s 266 DETACHED,
  up to 5.32 pts, also one site). Both are now candidates to displace v0.4.15, and
  they are not mutually exclusive -- but they are both GENERATOR lanes, and the
  encoder ladder has 1.28 pts of reachable work left in total (39 NOT_A_MIRROR +
  25 resonance residue). The honest reading is that the ENCODER ladder is nearly
  exhausted and the remaining distance to 100% is mostly generator work.

Owner decision required. Recorded 2026-07-28 by the v0.4.14 session; NOT self-approved.
```

### LADDER DECISION 2026-07-27 — accepted, by the project owner

```
Change:   v0.4.11  was  encode floor R3 (15 mol / 0.30 pts)
                   is   slot_renumber   (496 mol / 9.92 pts)
          v0.4.14-v0.4.18 shift up one; `structural` pulled out of v0.4.16 and re-filed
          at v0.4.17, re-labelled as generator capability; encode floor R3 folded into
          v0.4.18 as opportunistic work.

Because:  the ladder pointed its NEXT release at the SMALLEST block in the gap while
          18.26 of 27.54 points sat four and six releases out. 0.30 points is a rounding
          error against a 27.54-point gap, and encode floor R3 depends on nothing and
          unblocks nothing, so no prerequisite argument protects its position. Its real
          merit -- ARM 1 of the byte-identity gate takes ~15 min, almost all of it
          fixtures hitting `_resonance_candidates_isolated` at 120 CPU-s per fork -- is a
          developer-convenience argument, not a roadmap-points one.

Not changed, and why:
          v0.4.12 (honest acceptance) and v0.4.13 (harness false negatives) stay early.
          They are MEASUREMENT PREREQUISITES: shipping a 9-point canonicality fix before
          the harness's false-negative rate is known risks attributing an artefact to the
          fix. That argument covers 12 and 13. It does not cover a 0.30-point release.

Also corrected: the ladder carried STALE counts -- `slot_renumber (459 / 9.18)` and
          `rdkit_canonical (61 / 1.22)` predate the v0.4.8 honest baseline. The honest
          figures are 496 / 9.92 and 114 / 2.28, reproduced 2026-07-27 from
          `results-v0.4.8-honest/bucket_report_honest.md`.
```

### LADDER DECISION 2026-07-28 — accepted, by the project owner

```
Change:   v0.4.12  was  honest acceptance (2 lanes)
                   is   reflection parity (L1) + honest acceptance (L2), one lane each

Because:  v0.4.11 refuted its own fix but left the block REACHABLE and named the
          parity filter as the concrete next step. 7.54 pts is 27% of the remaining
          27.54, the oracle is already built (mirror_audit_donor_fold.py +
          injectivity/oracle.py), and the defect is pinned by
          TestDonorFoldCollapsesEnantiomers. Deferring it means the largest live
          block waits behind a release worth 0 points by design.

Not changed, and why:
          v0.4.13 (harness false negatives) stays where it is. PREFILTER_VETO's
          SCORING half was already closed by OIN_INDEP_SCORE in v0.4.8 -- measured
          on the frozen corpus, cheap-fails-but-independent-passes is 28/5000, and
          the honest metric already counts those correctly. What remains is purely
          an ACCEPTANCE-side defect, and it belongs there with the MEDZUR class.
          The v0.4.12 sketch's "PREFILTER_VETO may belong here" is DECLINED.

Override recorded: this breaks the roadmap's own "one theme, one or two lanes"
          rule. The release therefore carries TWO headline predictions, not one,
          and CLOSEOUT §3 diffs both separately.
```

### LADDER DECISION 2026-07-27 (v0.4.13) — two decisions, two different deciders

> ⚠ Attribution matters here and the two halves are NOT the same kind of call.
> **The promotion was accepted by the project owner**, against a stated confound, in answer to a
> direct question. **The `structural` re-sizing was decided by the session** on the strength of a
> measurement taken this release; it is recorded so the owner can overturn it cheaply, and the
> reordering it implies is deliberately left open rather than applied.

```
Change:   the DONOR FOLD PROMOTION.  [DECIDED BY: project owner] OIN_CANONICAL_DONOR_FOLD + OIN_FOLD_PARITY_VETO
          go DEFAULT-ON in v0.4.13. byte_exact 72.46% -> ~75.88% (+3.42 pts, 171
          molecules). This VOIDS the carry-forward licence: v0.4.13 owes a full
          ~55 CPU-h re-sweep and re-frozen ARM 1 / ARM 2 goldens.

Because:  v0.4.12 satisfied precondition 1 (a second uniform mirror-audit seed
          reading 0, with its baseline independently reproducing 19/250 on a
          disjoint draw). The remaining two preconditions are COMPUTE and
          BOOKKEEPING, not science. The owner elected to spend them here rather
          than defer to v0.4.14.

Confound recorded, raised BEFORE the decision and accepted: v0.4.13 is also the
          release that first measures the harness's false-negative rate, so a
          reader cannot attribute a headline move to the fold rather than to the
          harness without reading the two apart. Mitigation: the sweep that
          produces the headline runs with BOTH Lane levers OFF, so the +3.42 is
          attributable to the promotion alone, and Lane 1/Lane 2 are reported as
          measurements rather than as headline movers.

Change:   v0.4.17  [DECIDED BY: session, from this release's measurement]
                   was  `structural` (417 / 8.34 pts), "generator capability
                        floor ... bounded by what the generator can assemble"
                   is   RE-SIZED and RE-MECHANISED. 266 of its 417 molecules
                        (63.8%) are a RETURN-PATH GUARD DEFECT, not a capability
                        limit, and are worth up to 5.32 pts on their own.

Because:  v0.4.13 Lane 2 classified the attachment state of every stored generated
          structure (tools/attach_class_audit.py, #DONE 5000, 0 UNKNOWN).
          `structural` reads 266 DETACHED / 48 BOUNDARY / 98 INTACT / 5 none.
          The 266 did not fail because the generator could not assemble them --
          it assembled something and RETURNED it with ligands off the metal,
          because `_select_by_geometry`'s fallback ranking is not attachment-aware.
          The guard already exists for ACCEPTANCE (OIN_ATTACH_CHECK); applying it
          to RETURN is one site. v0.4.12 handed this forward by name as "the
          GAVSED class ... closing it changes arm A's behaviour and needs its own
          gate."
          Control: byte_exact reads 1.32% DETACHED against 24.11% on the failing
          side -- 18.2x enrichment -- so detachment discriminates and the 266 is
          evidence rather than an artifact of the predicate.

Not changed, and why:
          v0.4.13's own theme is CONFIRMED, not refuted. Its sketch carried a
          deletion clause -- "if both turn out to be single-molecule curiosities,
          this release should be deleted from the ladder" -- and it does NOT fire:
          GAVSED goes n=1 -> 280 and MEDZUR n=1 -> 99 on the genuine-failure cut.
          v0.4.14 (`rdkit_canonical` + winding residue) keeps its slot. The
          promotion consumes 171 of `slot_renumber` but leaves the 90
          frozen-resonance-form molecules v0.4.11 re-filed there.

Open for the owner: whether the 266-molecule return-path guard should DISPLACE
          v0.4.15 ("the 57 notation molecules", explicitly "knowledge, not
          points"). That would trade a 0-point release for a measured 5.32-point
          one, two rungs earlier. Not applied here -- it reorders three releases
          and the last such reorder was taken by the owner, not by a session.
```

#### What the v0.4.12 re-baseline established before any code was written

- 🔴 **`OIN_ETA_EARLY_EXIT` is runtime-inert AS SITED, and its promotion gate is VOID rather
  than unrun.** It lives in `_select_by_geometry_impl`, which runs *after*
  `generate_3d_structures` has filled the whole pool; the only site that can stop pool filling
  is `accept_fn`. Its own in-code A/B already said so (Ferrocene: fires, attempts 32 → 32).
  The sketch's Lane 2 — "run the corpus A/B its promotion gate demands" — would have spent
  real CPU confirming a documented structural null. **Overridden:** the criterion moved to
  `accept_fn` as `OIN_ETA_ACCEPT_EXIT`.
- **The v0.4.7 stored A/B arms no longer exist** (`spec/handoffs/v0.4.7/runs/` is gone), so
  every G1–G4 number must be re-*measured*, not re-scored.
- **The v0.4.6 accept-gap cohort classification is STALE.** All 8 of its `CHEAP_ONLY`/`GAP`
  molecules now satisfy the key inside `accept_fn` (telemetry: `adapter.early_exit_hit` fires
  in both arms, `adapter.eta_accept_*` never). The lever's real target population had to be
  re-derived from the frozen sweep: **405 eta molecules whose key does not match, 378 of them
  over 30 s.**
- **Three releases' byte-identity claims hold.** Re-encoding the frozen corpus's inputs and
  stored generated structures with today's encoder reproduces the v0.4.8 strings exactly on
  every molecule tested, so v0.4.9/v0.4.10/v0.4.11's "default path unchanged" is confirmed
  rather than assumed.



v0.4.8–v0.4.10 are written as full charter sets. v0.4.11–v0.4.18 are **SKETCH** and are promoted to
full one at a time by the close-out ritual below.

---

## Self-propagation

Every release ends by running `spec/handoffs/roadmap-100-100/CLOSEOUT.md`:

1. land the lanes · 2. re-measure and freeze the bucket report · 3. **diff predicted vs actual — a
miss is a deliverable** · 4. write the git-durable `docs/agentic-notes/<release>/*.md` ·
**4b. freeze the NUMBERS into `measurements/<release>/` via `tools/harvest_measurements.py`** ·
5. tag, bump, CHANGELOG ·
**6. generate `v0.4.<N+1>/` and promote it from SKETCH to full using the numbers just measured** ·
7. prune worktrees.

> Step 4b exists because steps 1–7 previously made only the *prose* durable.
> `tmCAT-tmPHOTO_xyz_dataset/` is gitignored in full, so six consecutive releases shipped
> write-ups whose underlying numbers lived on one disk — and **v0.4.11's mirror audit, the
> measurement that chartered v0.4.12, was lost outright.**

Step 6 runs `spec/handoffs/roadmap-100-100/PROMPT-next-release.md`, which requires the next session
to re-derive the gap decomposition, **re-check whether the ladder is still right**, and say
explicitly what the new numbers refuted. **v0.4.19 and beyond never need to be written by hand.**

---

## Standing traps

- **A sample that only exercises the common case confirms whatever you already believed.** On
  2026-07-26 a 10-molecule boron sample read "0/10 generate" and was used to propose a blanket
  fast-fail; the full 33 read **2/33** — `RAWJEG` (LIN, 2 slots) and `ULODUU` (TET, 4 slots) do
  assemble — and **the proposal is refuted**, because a blanket gate costs 2 real passes and no clean
  discriminator exists.
- **Circular measurement.** If the arm that scores an A/B uses the predicate the change accepts on,
  it cannot detect what the change costs. `OIN_ACCEPT_SCORED` read "18/22 both arms, zero
  regressions"; the honest arm read **15/20 → 7/20, 8 regressions, 0 fixes**.
- **Canonicalization can destroy what it canonicalizes.** Sorting symmetry-equivalent axes by sign
  made the axial token reflection-invariant and silently killed the chirality it encoded; every guard
  written against the single-axis fixture passed. **Prove the descriptor still flips for the mirror,
  on a multi-center case.**
- 🔴 **`byte_exact` CAN BE RAISED BY DELETING INFORMATION — and the comparison key will agree.**
  v0.4.11's donor fold moved **393 molecules into `byte_exact`, 0 in any other direction**, changed
  the key on **0 of 992** strings, and passed both gate arms — while **collapsing enantiomers in 221
  of those same 393 gains (56.2%)**. The metric and the key are *both* blind to reflection, because
  `compare._parse_vertex_colors` folds that axis deliberately. A one-directional transition matrix
  is therefore **not** evidence of safety for anything touching canonicalization.
  **Mirror-audit every canonicality lever (`tools/mirror_audit_donor_fold.py`) before quoting its
  points**, and confirm chirality independently with `tools/injectivity/oracle.py`. Two independent
  cohorts found it at comparable rates (uniform 19/250 = 7.6%, stratified 31/300 = 10.3%).
- **A partial run is not a result, and a progress line is not a tally.** During v0.4.11 the
  stratified audit was read mid-run as "0 regressions in the first 200" and written up as evidence
  that the cohort was the *wrong sample*. The tool prints a verdict only every 50th molecule, so
  four clean progress lines had been mistaken for 200 clean molecules; run to completion the same
  cohort reported **31 collapses**. The "wrong stratum" conclusion was an artifact of reading a
  sampled log as a census. **Wait for the summary line, or count with `--verbose`.**
- **A fragment's automorphism says nothing about the parity of the vertex permutation it induces.**
  Two donors interchangeable in the *isolated ligand graph* can sit at vertices whose exchange is an
  **improper** operation on the coordination sphere. This is why folding is restricted to proper
  rotations; narrowing a fold's scope *within a fragment* is not a substitute for a parity check.
- **`OIN_EMIT_AXIAL`'s promoting evidence is stale** — measured with `OIN_CANONICAL_PERCEPTION` OFF,
  now default-ON; YESKOZ's hindered-axis count goes 2 → 1 under it.
- **P3's obvious fix is measured-wrong** — stamping `[N@]` after the sanitize moves the canonical
  write order, and `@`/`@@` is a parity relative to that order; RIFGUJ's pseudo-asymmetric ring
  carbons then flip between a structure and its mirror.
- **vdW clash count is not a sufficient quality proxy** — `POVPIA_comp_0` goes 16 clashes → 0 while
  re-perception shows a detached hydrogen and a C–N read as C=N.
- **"Unset means off" is a trap.** `os.environ.get("OIN_EMIT_AXIAL")` returns `"0"`, a non-empty
  string, and *enables* the lever. Always go through `oin/levers.py::lever_enabled`.
- **The shared checkout moves under you** — one `.git/index` across ~33 worktrees.
- 🔴 **`metrics.elapsed_s` is a SUM, not a duration** (v0.4.9). Up to three separately SIGKILLed
  harness attempts are added into one field, so the ceiling is `3 × mol_timeout`. The "759.9 s
  against a 300 s budget" figure that chartered a whole release is `300 + 300 + 160`; split by
  `tier_passed`, all 4658 single-attempt rows finish within **0.2 s** of their cap. This is the
  **second** trap in the same field, and it survived one audit because everyone was busy
  remembering the first one (it is nested).
- 🔴 **A string-equality count is not a pass count** (v0.4.9). Comparing `smiles_1` to the honest
  round-trip string gives 3631 passes; the bucket report gives 3623, because it applies the
  `status` gate first. The eight are exactly the **atom-count gate** population — byte-identical
  strings with atoms missing. Use `tools/roundtrip_bucket_report.py`, never an ad-hoc
  `honest_class.endswith("->byte")`.
- 🔴 **A gate can silently run on the wrong interpreter** (v0.4.9). `gate_v047.sh` globbed
  `$(dirname $REPO)/*/.venv/bin/python` because worktrees have no `.venv`; from a worktree it
  selected an unrelated project's venv, then one with rdkit **2025.09.2** against the pinned
  **2025.9.3**. A byte-identity gate on a different rdkit reports MISMATCHes that read as code
  regressions. Resolve via `git rev-parse --git-common-dir` and refuse on version drift.
- 🔴 **`grep -c '[p]rocess_name'` matches your own shell** (v0.4.9). A wait-loop built on it never
  exits — twice in one session, once in a leftover watcher from an earlier session that had been
  spinning for hours. Match on something the watcher's own command line cannot contain.
- **Profile before bounding.** v0.4.9's charter named the unbounded CBC solve and the 48–57 s
  `accept_fn` re-encode as the cost sinks. Measured: **2.1%** and **0.8%**. The sink is
  `embed.get_embedding` (61.5 s of *self* time in an 82.4 s generation). Bounding either suspect
  would have measured as "no change" and shipped nothing while looking like a fix.

---

## Sources

`docs/agentic-notes/v0.4.6/SWEEP_v0.4.6_5K.md` · `docs/agentic-notes/v0.4.6/METRIC_FALSE_POSITIVES.md` · `docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md` §10 ·
`docs/agentic-notes/v0.4.5-retrospective/` · `src/oinsmiles/oin/levers.py` ·
and on `swimlane/v047-*` until v0.4.7 lands: `docs/agentic-notes/v0.4.7/ACCEPT_SCORED_v0.4.7.md`,
`docs/agentic-notes/v0.4.7/ATTACH_CHECK_v0.4.7.md`, `docs/agentic-notes/v0.4.7/ENCODE_FLOOR_v0.4.7.md`, `docs/agentic-notes/v0.4.7/COHORT_v0.4.7.md`,
`docs/agentic-notes/v0.4.7/BORON_GEN_CEILING_v0.4.7.md`
