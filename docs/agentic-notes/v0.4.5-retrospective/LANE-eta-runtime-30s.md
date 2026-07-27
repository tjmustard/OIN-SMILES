# Lane: eta runtime — the `<30 s per job` goal

**Status:** mechanism FOUND, fix BUILT and MEASURED, lever **default OFF**.

> ⚠ **§3.2's "pass-rate: neutral" and §4's "promote once the sweep confirms" are SUPERSEDED.**
> A parallel v0.4.7 lane (`swimlane/v047-promote`, `docs/agentic-notes/v0.4.7/ACCEPT_SCORED_v0.4.7.md`) built on these
> commits and added the arms this cohort A/B lacked. Its **G2 gate FAILS**, and it fails on a flaw
> in *my* measurement: `passed` here is computed with `get_oin_string(gen.mol, coords)` — the very
> predicate the lever accepts on — so it is **circular** and structurally cannot detect what
> dropping the independent confirm costs. With a genuinely independent arm: **indep 15/20 → 7/20,
> 8 regressions, 0 fixes**, and the 8 are not cosmetic (6 lose haptic coordination outright, 1
> reassigns the donor atom, 1 detaches a hydrogen). Its G3 gate does pass — the emitted string is
> byte-identical 20/20 — so the honest one-line trade is **byte-identical notation, changed
> geometry**. Read that document, not this section, for the promotion decision. Everything in
> §2 (the mechanism), §2.4 (`timeout` is advisory), §2.5 (`PREFILTER_VETO`), §3.3/§3.4 (quality,
> and the dead metric) stands.
**Branch:** `research/eta-relax-invariant-proxy` (`e11a4cfe`, `15b7b7ff`, + this report). **UNLANDED
at time of writing.**
**Supersedes:** commit `41d2f52e` ("the eta pool cost is structural, not a defect. Item closed.")
— that verdict was **wrong**, and this document explains why it was wrong and how it was caught.

---

## ELI5

Imagine a factory that has to produce one acceptable widget. It builds a widget, holds it up to a
very strict inspector, and the inspector says no. So it builds another. And another — forty-eight
times, for an hour and a half.

Meanwhile the *customer* has a much simpler test than the inspector does. And it turns out the
**very first widget** already passed the customer's test. Every widget after that was built for
nothing.

That's the eta runtime problem. The pool-filling loop asks a stricter question than the question
that actually decides whether the round trip succeeded. Ask the *right* question and the median
molecule finishes in 3.6 s instead of 16 s, and the count of molecules over 30 s drops from 10 to 3.

The twist: the strict inspector was not being silly. It is the only test in the whole pipeline
that doesn't take the generator's word for anything. So removing it makes us *faster* and makes
our reported score *self-consistent* — but it does not make the chemistry more correct. Both facts
are in this report, because only reporting the first one would be dishonest.

---

## The shape of the work

```
   GOAL: "would be great to also do so with each job running less than 30 seconds"
                                  │
        ┌─────────────────────────┴──────────────────────────┐
        │  PRIOR ATTEMPTS (all reasoning, no profiling)       │
        │  0a7d1a33  "cost-per-attempt, not attempt count"    │ ← WRONG
        │  860c1988  eta early-exit lever: fires, ineffective │ ← true, unexplained
        │  41d2f52e  "RESOLVED: structural, not a defect"     │ ← WRONG, closed the item
        └─────────────────────────┬──────────────────────────┘
                                  │  reopened by actually profiling
                                  ▼
   STEP 1 · cProfile one mid-range eta molecule  (tools/profile_eta.py)
        HIDCIH_comp_1: generate 96.2s
          get_embedding      47.0s / 33 attempts
          ff_clean           34.5s / 33 attempts
          accept_fn          14.2s / 26 calls ... and returns False EVERY time
                                  │
                          "...on a molecule the sweep scores as a SUCCESS"
                                  ▼
   STEP 2 · read the two predicates side by side — they are NOT the same test
        score     : key(get_oin_string(gen.mol, coords))     == key(oin_in)
        accept_fn : cheap prefilter (can only REJECT)
                    AND key(XYZToSMILES().convert(xyz))      == key(oin_in)   ← extra
                                  │
                                  ▼
   STEP 3 · force a full pool fill, record BOTH verdicts per conformer
            (tools/probe_accept_gap.py)
        HIDCIH: 48 conformers │ cheap 46 matches, FIRST AT INDEX 0 (1.66s)
                              │ strict 2 matches, first at INDEX 25 (49.4s)
                                  │
                        n=1. Second molecule CONTRADICTED it (YIYGAP: both at 0).
                                  ▼
   STEP 4 · 20-molecule stratified cohort  (tools/summarize_accept_gap.py)
        CHEAP_ONLY     4 (20%)  scored-success exists, strict NEVER matches
        GAP            8 (40%)  cheap fires earlier than strict
        PREFILTER_VETO 1 ( 5%)  ← runs OPPOSITE; possible ACCURACY defect
        NO_GAP         5 (25%)  all fast already; lever irrelevant
        DEAD           2 (10%)  nothing ever reached the predicate
        17/17 with a scored-success hit it at CONFORMER 0.
        83% of observed wall-clock came AFTER that hit.
                                  ▼
   STEP 5 · build OIN_ACCEPT_SCORED, A/B on 3 arms  (tools/ab_accept_scored.py)
        runtime ✓   pass-rate ✓   structure quality ✓ (and BETTER, not worse)
                                  ▼
   STEP 6 · held OFF. Gate text says CORPUS; 22 molecules is a cohort.
            Recommendation recorded: PROMOTE after a lever-ON sweep.
```

---

## 1. Initial assumptions and hypothesis

### 1.1 What the project believed going in

Three prior conclusions were on record, and **two of them were mine and wrong**:

| commit | claim | verdict now |
|---|---|---|
| `0a7d1a33` | the eta tail is **cost-per-attempt**, not attempt count | **wrong** — HIDCIH spends 33 attempts vs JAKRUL's 3, at comparable per-attempt cost |
| `860c1988` | the eta early-exit lever **fires but is measured ineffective**; "the cost is the pool fill" | **true, and now explained** — it fires, then `accept_fn` rejects the conformer anyway |
| `41d2f52e` | **RESOLVED — structural, not a defect.** `accept_fn` sees raw pool conformers; the key only matches post-relaxation, which happens after the fill loop, so no acceptance predicate can shorten it | **wrong on the premise.** `_try_accept` calls `cleaner.clean_geometry` *before* filing, so `accept_fn` sees a **cleaned** conformer. The stage attribution was invented, not read. |

`levers.py`'s `OIN_ETA_EARLY_EXIT` entry had actually recorded the right instinct all along — *"the
cheap early-exit predicate is STRICTER than the one that decides success, so eta molecules pay the
full widened pool for nothing"* — but nobody had identified **which** part of the predicate was the
strict part, so the follow-up went to the wrong layer.

### 1.2 The hypothesis this lane started from

> The eta cost is the conformer pool filling to `target_pool` when an acceptable conformer already
> exists inside it. If so, the acceptance predicate — not the embed or the force field — is the
> lever.

### 1.3 What I expected to be wrong about

Explicitly, before measuring: **that bypassing `_select_by_geometry`'s clash-first ranking would
degrade structure quality.** `clash.VDW_ACCEPTANCE_ENABLED` defaults ON and selection sorts by
clash count first, so accepting the *first* scoring conformer instead of the *best-ranked* one
should hand back worse geometry. This is why the A/B measures quality as a first-class arm. **The
concern was refuted** — see §3.3.

---

## 2. What was found out

### 2.1 The two predicates are different tests, and one is a subset of the other

```
harness score  :  canonical_roundtrip_key(oin_in)
                    == canonical_roundtrip_key(get_oin_string(gen_result.mol, coords))

accept_fn      :  _reencode_oin_fast(cmol) key match      ← same family as the score
                  AND canonical_roundtrip_key(XYZToSMILES().convert(xyz)) == target
                                                           ← an INDEPENDENT re-perception
                                                             the score never asks for
```

`tools/test_dataset_roundtrip.py` scores with `get_oin_string(gen_result.mol, coords)` — the
generator's own perceived mol. `accept_fn` additionally writes the geometry to XYZ and re-perceives
**everything from coordinates alone**. That second test is strictly harder, so `accept_fn` ⊂
harness-success. Conformers in the difference are scored successes that acceptance discards.

### 2.2 The gap is large, and it is at conformer 0

`tools/probe_accept_gap.py` monkeypatches the predicate to record both verdicts and always return
False, forcing the pool to fill so every conformer is observed. HIDCIH_comp_1:

| | matches | first at | wall-clock of first |
|---|---|---|---|
| cheap (what the score uses) | **46 / 48** | **index 0** | **1.66 s** |
| strict (what acceptance uses) | 2 / 48 | index 25 | 49.4 s |

44 conformers were scored-successes that acceptance rejected. The unpatched run spends 96 s
reaching the one conformer step 2 accepts.

### 2.3 It is heterogeneous — and n=2 was enough to prove it

The second molecule tested (YIYGAP_comp_0) showed **no gap at all**: both predicates fire at
conformer 0, total 1.9 s. Had I stopped at HIDCIH I would have shipped a universal claim from a
single fixture — the exact failure mode recorded in [[mirror-audit-confirms-coupled-false-positives]].
The 20-molecule stratified cohort (`docs/agentic-notes/v0.4.6/eta_accept_gap_cohort.md`):

| case | n | meaning | does the lever help? |
|---|---|---|---|
| `CHEAP_ONLY` | 4 (20%) | scored-success exists; strict **never** matches in the whole pool | **yes, maximally** |
| `GAP` | 8 (40%) | cheap fires earlier than strict | **yes** |
| `PREFILTER_VETO` | 1 (5%) | strict matches 16/48, cheap matches 0 | **no** — see §2.5 |
| `NO_GAP` | 5 (25%) | both fire at index 0; all already 1.5–6.6 s | irrelevant |
| `DEAD` | 2 (10%) | nothing ever reached the predicate | no |

Two invariants held across every molecule with a scored-success:

1. **17 / 17 hit it at conformer 0.** Not "early" — *first*.
2. **83 %** of observed wall-clock (923 s of 1114 s) was spent *after* that hit.

### 2.4 A second, independent finding: `timeout` is not a bound

Measured while timing boron molecules: 60 s requested, **60.7–137.9 s spent** (GOHWOQ 2.3× over).
`embed_time_budget=self.timeout` bounds the embed *attempt loop*, not the OIN-direct assembly
around it. Only the harness's per-molecule SIGKILL subprocess actually enforces a budget.
Consequences: any timing taken without that watchdog **understates** the tail, and my own first
A/B script inherited the flaw and sat on one molecule for 30+ minutes. `tools/ab_accept_scored.py`
now imposes its own hard cap, applied identically to both arms.

### 2.5 `PREFILTER_VETO` — the case that runs opposite, and is not a latency bug

AROHIA_comp_0: the cheap test matched **0/48**, the strict test matched **16/48**. Because
`_reencode_key_matches` returns `False` on a cheap mismatch *before* the strict test runs,
production accepts **neither** — 16 independently verifiable conformers are unreachable today. This
is the one case in the cohort that is a candidate **accuracy** defect rather than a runtime one, and
`OIN_ACCEPT_SCORED` cannot help it (it makes the cheap test *more* decisive, not less). Left open,
named, and not conflated with the runtime story.

---

## 3. What was done

### 3.1 `OIN_ACCEPT_SCORED`

`_reencode_key_matches` gains `independent_confirm=True`. When the lever is on, a cheap match with
no stretched bond accepts immediately, skipping step 2. `fast is None` (perception failed) still
falls through to step 2 rather than accepting blind — an unperceivable conformer is not a scored
success. Default OFF, so with the lever unset every call is byte-identical to before.

### 3.2 The A/B, three arms, two runs

Both runs use the same 22-molecule cohort, stratified across the capstone elapsed distribution
plus 6 non-eta controls, one subprocess per molecule (the lever is read at predicate-construction
time, so one process cannot host both arms).

| | run 1 · no hard cap | | run 2 · hard cap + honest clash metric | |
|---|---|---|---|---|
| | **A** default | **B** scored | **A** default | **B** scored |
| pass | 18/22 | **18/22** | 16/22 | 18/22 |
| median | 16.01 s | **3.63 s** | 13.87 s | **5.54 s** |
| total | 1980.7 s | **626.0 s** | 1202.9 s | **351.3 s** |
| > 30 s | 10 | **3** | 8 | **2** |
| clash total | *not measured* | *not measured* | 16 over 1/17 | **2 over 2/19** |
| severe clashes | *not measured* | *not measured* | **7** | **0** |
| worst_overlap min | — | — | **0.4344** | 0.7283 |
| regressions | — | **none** | — | **none** |

**Read them together, not separately.** Run 1 has honest pass counts but a dead quality metric
(§3.4). Run 2 has an honest quality metric but its arm-A pass count is truncated by the hard cap —
its "2 fixes" (GAVSED, QIDKUL) are cap artifacts, both having passed uncapped in run 1 at 302.5 s
and 390.3 s. The defensible combined statement:

- **Pass-rate: neutral.** 18/18, zero regressions, on the uncapped comparison.
- **Runtime: median 3–4× faster; molecules over 30 s cut from 10 to 3.**
- **Structure quality: not degraded — improved on this cohort.**

### 3.3 The quality concern was refuted

My a priori expectation (§1.3) was that arm B would hand back worse geometry. It did not:

| molecule | arm A | arm B |
|---|---|---|
| POVPIA_comp_0 | clash **16**, severe 7, worst **0.4344** | clash **0**, worst 0.75 |
| RATPEK_comp_0 | clash 0, worst 0.7599 | clash 1, worst 0.7461 |
| DAKGON_comp_0 | clash 0, worst 0.8164 | clash 1, worst 0.7283 |

Arm B trades two mild single clashes for the elimination of a badly clashing structure — 7 severe
clashes to **0**. This is consistent with a comment already in `_select_by_geometry`: minimising
clash "has donors splayed to the edge of the gate and re-perceives as detached." Clash-first
selection is not reliably better geometry, and on POVPIA it was substantially worse.

### 3.4 A defect found in my own promotion gate

The first A/B built its quality arm on `clash.mol_clash_count(gen_result.mol)`. That function
duck-types on `mol.atom_list` and **returns 0 on `AttributeError`** — and `gen_result.mol` is a
bare `rdkit.Chem.rdchem.Mol`. So it returned 0 for **all 44 measurements in both arms**, and the
gate would have certified the quality arm without measuring it. The whole reason that arm existed
was that this lever bypasses clash ranking.

Caught by asking why a column was uniformly zero. Fixed by computing
`vdw_clash_count(positions, atomic_numbers)` from the returned coordinates, and by recording the
**continuous** `worst_overlap` beside the thresholded count — a real reading looks like `0.7502`, a
dead one looks like `0`. Recorded as a general rule in `[[degenerate-metric-certifies-nothing]]`.

### 3.5 Boron, found in passing

The `DEAD` case XIQKOY_comp_0 took 309.7 s here against **1.5 s** in the capstone. Cause: the
v0.4.6 `OIN_BORON_CAGE` promotion. Lever OFF, the encoder emits the cage as a disconnected fragment
and generation raises `UncoordinatedFragmentError` in 0.01 s (total 0.87 s). Lever ON, the encode is
a *correct* fully-coordinated B₁₀ cage and generation runs past 340 s. Sampling 10 of the 34 boron
molecules: **0/10 produce a 3D structure**, 6/10 burn the whole cap. So the promotion moved this
class from failing *instantly* to failing *slowly* — roughly 2.8 CPU-hours per full sweep for zero
extra passes. Full write-up in `docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md` §10; that document's headline "34 now
encode and round-trip" was also corrected, because all nine of its 34/34 checks are encoder-side
and none invokes the generator.

---

## 4. Where this landed

**Built, committed to `research/eta-relax-invariant-proxy`, lever default OFF.**

Not promoted, deliberately. The gate text in `levers.py` says **corpus** A/B; 22 stratified
molecules is a cohort. Promoting on n=22 is precisely the mistake this lane exists to correct —
two fixtures produced four wrong answers about this same tail, twice by my own hand. The corpus
A/B is a full sweep with the lever ON, and it cannot run while a default-path 5 k sweep is in
flight. **Recommendation on the evidence: promote once that sweep confirms.**

### Against the `<30 s` goal

| | molecules > 30 s |
|---|---|
| cohort, default | 10 / 22 |
| cohort, `OIN_ACCEPT_SCORED=1` | **3 / 22** |

And the three residuals each have a **different named cause**, none of which this lever addresses:

1. **AROHIA** — `PREFILTER_VETO` (§2.5). An accuracy question, not a latency one.
2. **QIDKUL** — genuine cost-per-attempt; 390 s → 166 s, still over.
3. **XIQKOY** — boron generator ceiling (§3.5). Assembling a polyhedral borane cage from m-SMILES
   is an open `generator3d` problem.

For context on the population this targets: **78 % of eta molecules exceed 30 s** in the capstone
(1221/1565), median eta success 54.5 s. Eta is the right lane; this lever is the largest single
runtime lever found so far; and it is not sufficient on its own.

### Open, in priority order

1. Corpus A/B with the lever ON → promote (or refute at scale).
2. `PREFILTER_VETO`: why does the cheap prefilter veto conformers the authoritative test accepts?
3. Boron generation: cage assembly, or a fast-fail so the class stops consuming the budget.
4. Cost-per-attempt profiling for the QIDKUL class — `_finalize_positions` → `ic.update_xyz` →
   `pinv`/`svd` was 15.6 s over **1189 SVD calls** on HIDCIH and is the obvious next target.

---

## 5. Reproduce

```bash
cd <worktree>; V=/path/to/main/.venv/bin/python; export PYTHONPATH=$PWD/src

# 1. where does one eta molecule's time actually go?
$V tools/profile_eta.py <dataset>/cat/HIDCIH_comp_1.xyz --top 22 --out prof.txt

# 2. the two verdicts, per conformer, full pool fill
$V tools/probe_accept_gap.py <dataset>/cat/HIDCIH_comp_1.xyz --json gap.json
#    add --stop-after-both for a cohort-affordable run

# 3. the population split
$V tools/summarize_accept_gap.py <dir-of-gap-json> --md cohort.md

# 4. the three-arm A/B (hard cap, because `timeout` is advisory -- see §2.4)
$V tools/ab_accept_scored.py --cohort cohort.json --out ab.json --workers 2 \
     --timeout 300 --hard-cap 330

# 5. boron: same molecule, one lever, two outcomes
for v in 0 1; do OIN_BORON_CAGE=$v GEN_CAP=60 $V tools/boron_gen_time.py XIQKOY; done
```

**Runtime caveat that applies to every number above:** they were taken on a box concurrently
running a 5,000-molecule sweep (load ~34–42 on 12 cores). Read them as **ratios within a run**.
The pass-rate, clash, and case-classification results are load-independent; the seconds are not.
