# v0.4.15 — what the first real generator sweep since v0.4.6 corrected

**Release:** v0.4.15 · **Measured:** 2026-07-29 · **Source:**
`results-v0.4.14-sweep`, N=5000, 6 shards, BLAS threads capped to 1, 38.7 CPU-hours / ~7.7 h wall.
Frozen: `measurements/v0.4.14-sweep/` and `measurements/v0.4.15/`.

This note records the corrections. Four of them change a number the v0.4.15 charter was written
against, and two of them change what a lane should *build*.

Completeness was verified by **report count, not exit status** — `individual_reports` = 5000,
plus `SHARDS_DONE reports=5000`, `SUMMARY_DONE`, `LANE2_REDERIVE_DONE`. That matters here: this
sweep's `PROGRESS` contains a `SUMMARY_DONE` sentinel echoed *after* `rebuild_summary.py` had
already failed on a missing flag, so the sentinel reports success independent of the step it
follows.

---

## 1. 🔴 The charter's loudest claim is REFUTED: the absolute was never over-stated

The v0.4.15 charter opens with it, in three separate files:

> **THE ABSOLUTE BASELINE IS NOT TRUSTWORTHY.** `byte_exact` is published as 77.30%, but no real
> generator sweep has run since v0.4.6. […] v0.4.13's published 75.88% measures ~74.88%; v0.4.14's
> 77.30% is a *measured delta on an over-stated base*.

**The sweep reads 77.16%.** The published figure was over-stated by **0.14 points.**

The offline re-score chain was accurate **in aggregate**. What it got wrong is the bucket
**composition**:

| bucket | charter (offline) | **SWEEP** | Δ |
|---|---:|---:|---:|
| `byte_exact` | 3865 / 77.30% | **3858 / 77.16%** | −7 |
| `structural` | 417 / 8.34% | **484 / 9.68%** | **+67** |
| `hard_fail` | 319 / 6.38% | **266 / 5.32%** | **−53** |
| `key_equal` | 361 | **365 / 7.30%** | +4 |
| `facmer_divergent` | 16 | **15** | −1 |
| `encode_fail` | 15 | **12** | −3 |

Sums exactly to 5000, and the failing blocks sum to 22.84 = `100 − 77.16`.

`structural` +67 against `hard_fail` −53 is the substantive movement: **~53 molecules that used to
produce nothing now produce something that is not good enough.** That is the generator improving
across v0.4.9–v0.4.14 and re-filing its own failures — and an offline re-score holds the structure
fixed, so it is structurally incapable of seeing it.

### The corrected lesson, which is narrower than the one that was recorded

> *An offline re-score chain over stale structures drifts WHICH bucket a molecule lands in, not the
> headline total.*

This does **not** retract the rule that a generator-changing lever must be measured by real
generation — both v0.4.15 lanes change what the generator returns, and both are measured that way.
It retracts only the claim that the absolute had drifted by ~2.4 points, which two releases spent
real caution on and which turned out to be 0.14.

⚠ **A trap this file nearly fell into.** `roundtrip_bucket_report.py` defaults to
`--score scored`, which reads **86.88%** on this same sweep — the circular predicate
`OIN_INDEP_SCORE` replaced at a measured 9.6% false-positive rate. **The scored/honest gap on this
sweep is 9.72 points.** Quoting the default is the v0.4.8 mistake; `--score honest` is required.

## 2. Both lanes are BIGGER than chartered

| lane | charter | **sweep** |
|---|---|---|
| 1 — `structural`/`DETACHED` | 266 / 5.32 pts | **301 of 484 / 6.02 pts** (62.2%) |
| 2 — `key_equal` enantiomer | 183 / 3.66 pts | **201 of 242 / 4.02 pts** (83.1%) |
| release total | 449 / 8.98 | **502 / 10.04** |

`hard_fail` is **262 of 266 `NO_STRUCTURE`** — a genuine capability floor, correctly out of scope.

## 3. 🔴 The two attachment predicates are DIFFERENT TESTS — Lane 1's exposure is 1, not 52

The charter points Lane 1 at `structural`/`DETACHED` and never states what the change could
**cost**. 52 of the 3858 `byte_exact` molecules also read `DETACHED`, so an attachment-preferring
return path could reorder the pool away from a conformer that currently passes. Taken at face
value that is −1.04 points of exposure.

`tools/attach_return_preflight.py` (new; output frozen at
`measurements/v0.4.15/attach_return_preflight.json`) evaluates the **guard's own** predicate — a
claimed coordination *site* holding no atom, coordinate-only, with the OIN's distinct slot count
as the claim:

| population | n | `SITE_LOST` | |
|---|---:|---:|---|
| `structural`/`DETACHED` | 301 | **289 (96.0%)** | Lane 1's target, 5.78 pts |
| `byte_exact`/`DETACHED` | 52 | **1 (1.9%)** | Lane 1's exposure, **0.02 pts** |
| `byte_exact`/`INTACT` | 250 | **0 (0.0%)** | control |

**50× separation, and the exposure is 52× smaller than the bound.** The reason:

- **`coordination_report`** (what `attach_class_audit.py` reports) compares the **input** donor set
  against the **generated** donor set — *did the donor set change.*
- **`attach_check.ligands_attached`** (what the guard uses) asks whether a claimed **site went
  empty**, on the generated geometry alone.

A molecule can trip the first and pass the second: swapping which hydride the metal binds changes
the donor set without emptying a site. The `coordination.reason` strings show that is exactly what
the two populations are:

```
the 52 :  1-3 LIGHT/AMBIGUOUS donors, count often unchanged or HIGHER
          Ag: lost {H:1} (4 -> 4)    Zr: lost {Si:1} (13 -> 14)    Y: lost {F:1} (7 -> 12)

the 301:  WHOLE MULTI-CARBON HAPTIC GROUPS -- a Cp/arene off the metal wholesale
          13x Ni: lost {C:6} (7 -> 1)     11x Fe: lost {C:5} (8 -> 3)
           7x Fe: lost {C:8} (8 -> 0)      9x Ru: lost {C:6} (9 -> 3)
```

⚠ **The `INTACT` control arm is what makes 1.9% a reading instead of a broken claim count.** A
probe that miscounts the claim fires on molecules with no attachment problem at all; this one fires
on **0 of 250**. Without that arm, "1/52" and "the probe is wrong" are indistinguishable — which is
the shape of instrument this project has now caught printing plausible nothing six times.

**Generalisation worth carrying:** *a bucket label and a guard's predicate are two different
measurements even when they share a name.* v0.4.14 recorded "a bucket name that asserts a cause is
a hypothesis"; this is the same error one level down — the bucket's *verdict* is also not the
guard's verdict.

## 4. 🔴 Lane 2 as chartered cannot work — and the replacement was chosen by measurement

The charter's Lane 2: use `oin/metal_config.py` as an acceptance predicate, do not emit the token.

**An acceptance test needs a reference handedness, and the generator's only input is the OIN
string.** `_select_by_geometry_impl` already contains a helicity-aware branch, and it is dead by
construction — `parse_metal_config_token(parsed.original_oin)` returns `None` for the entire corpus
because `OIN_EMIT_METAL_CONFIG` is off, and emitting the token is v0.4.16 by the charter's own
scoping. So the chartered lane has nothing to compare against.

Measured over the 183 prior `MIRROR_MATCH` molecules, re-read from this sweep's own
`smiles_1` / `smiles_2_indep`:

```
norm_differs = 183/183      key_same = 183/183      (norm_same 0, key_differs 0)
```

**The handedness survives `normalize_oin_for_comparison` and is folded only by the key.** The
mechanism is visible in the strings:

```
AFADOC_comp_0  input  ...c2O{5})c(O{4})c...      generated  ...c2O{4})c(O{5})c...
AGAVIQ_comp_0  input  ...P{3}(...OP{4}...)...    generated  ...P{4}(...OP{3}...)...
```

A **transposition of two same-coloured donors** — an odd permutation, i.e. a reflection.
`_parse_vertex_colors` colours every donor atom of a ligand with that ligand's *whole* canonical
body, so the two vertices are interchangeable and `_polyhedron_signature` (lexmin over the proper
rotation group) cannot distinguish them. The normalized string keeps absolute slot numbers, so it
can.

⇒ **Lane 2's predicate is normalized-string equality** (`OIN_ACCEPT_STRING_EXACT`). It invents no
descriptor, emits nothing, and separates 183/183 by construction. It is also exactly the charter's
own stated rule: *a predicate that decides what to accept must be sensitive to every axis the
answer is graded on* — `byte_exact` is graded on `smiles_1 == smiles_2`, so acceptance must compare
strings.

**Normalized, not raw.** Raw carries the metal `@OH`/`@SP` labels the encoder documents as
atom-order-dependent and irreproducible, so raw equality would reject nearly everything.
`normalize_oin_for_comparison` strips precisely those and keeps slots and winding — strictly
between the key (too loose) and raw (too strict).

**Rejected alternative:** a narrow "is this difference a reflection?" test on the string. Reflection
parity is not a property of the emitted string — a donor swap is a transposition fixing every other
vertex and does not preserve the Gram matrix, so the obvious `det > 0` test rejects *every* swap.
Refuted in v0.4.12; not re-opened.

## 5. 🔴 The lanes are NOT independent, contrary to both lane docs

Both lane docs state "Depends on: nothing. Blocks: nothing." **Mechanically false.**

An accepted conformer is returned as the **sole pool member** (`return [early_hit]`). So a stricter
acceptance predicate means more molecules exhaust the pool and fall through to
`_select_by_geometry_impl`'s fallback return — **which is the exact unguarded site Lane 1 exists to
fix.** Lane 2 makes Lane 1 load-bearing.

Two consequences, both acted on:

1. **Lane 1 lands first.**
2. **Lane 2 cannot simply reject the incumbent.** Rejecting a key-equal conformer would let the
   energy-sorted pool hand back a *different* mol — a real regression. Hence
   `generator3d.ACCEPT_INCUMBENT`: a third `accept_fn` verdict meaning "acceptable to fall back on,
   not to stop for". The pool keeps filling; if nothing string-exact appears, the first incumbent is
   returned, which is byte-identically what the pre-lever run returned. The lever therefore costs
   latency, not accuracy.

Both lanes use the same shape: **prefer better, never return worse.** The charter framed both as
rejection filters and then correctly warned that rejection converts silent wrong answers into loud
failures and lowers the headline. That trade is real but avoidable, and it is now a *separate*
lever in each lane (`OIN_ATTACH_RETURN_STRICT`) rather than the default behaviour.

## 6. Runtime moved — and it is NOT like-for-like

Recomputed from the **nested** `metrics.elapsed_s` (present on 5000/5000; a top-level read yields
`None`), which is also a **SUM** over up to three separately SIGKILLed attempts:

| | v0.4.8 era | **SWEEP** |
|---|---:|---:|
| `> 30 s` | 994 / 19.88% | **678 / 13.56%** |
| median | 7.19 s | **4.01 s** |
| max | 759.9 s | **728.8 s** |

⚠ **This sweep capped BLAS threads to 1 and the v0.4.6 sweep did not.** Direction is solid,
magnitude is approximate. The runtime improvement is **not claimed as a win** — accuracy is what
the run was for.

## 7. An open question, stated rather than explained away

`veto_outcome_audit.py` on the fresh structures: the parity veto reverted **all 242** folds it could
measure — `kept = 0`, `decided_against` 242/242 — against **171 kept** on the older corpus.

Either the freshly generated structures are systematically more chiral, or the veto behaves
differently against them. **Lane 2's predicate does not depend on this** (it is a direct three-way
encode comparison, not a veto reading), so it is not a blocker. But it is unexplained, and a
`kept = 0` reading is exactly what a broken outcome classifier would also print.

## 8. 🔴 IN-FLIGHT: Lane 2's lever fires correctly and recovers nothing on the first two molecules

The lever was smoke-tested before the population arms were launched, on the "ask what a broken
version would print" rule. It printed the *right* nothing, and the telemetry says why.

A/B over 3 known `MIRROR_MATCH` molecules, real generation, honest scoring:

```
byte_exact  OFF=False ON=False   3     REAL gains 0   REAL losses 0
input string moved: 0                  GENERATED output moved: 0
```

Telemetry on two of them, `OIN_TELEMETRY=1 OIN_ACCEPT_STRING_EXACT=1`:

| site | AFADOC_comp_0 | AGAVIQ_comp_0 |
|---|---:|---:|
| `adapter.string_exact_incumbent` | 2 | 5 |
| `pool.accept_incumbent_recorded` | **1** | **1** |
| `pool.accept_incumbent_returned` | 1 | 1 |
| `adapter.string_exact_early_exit_incumbent` | 1 | 1 |
| result | `key_equal`, not `norm_equal` | same |

**The wiring is live and the predicate fires.** `early_exit` defaults to on
(`os.environ.get("OIN_EARLY_EXIT", "1") != "0"`), so `accept_fn` is built on the default path;
`accept_incumbent_recorded` proves the sentinel reached the pool filler, and
`accept_incumbent_returned` proves the retain-incumbent path executed. The null is **not** a
wiring failure — which is the distinction the whole smoke test existed to draw.

🔴 **`accept_incumbent_recorded = 1` is the finding.** Exactly **one** conformer in the entire pool
carried the requested key, and it was the mirror. The pool does not hold an alternative
handedness, so **no acceptance predicate can fix these molecules — only construction can.** That
is Lane 2's chartered Q4 ("is the correct handedness anywhere in its pool?") answered *negative*
for these two, and the charter was right to demand it before a fix:

> Shipping a filter that rejects the mirror without knowing the pool holds an alternative converts
> 183 silent wrong answers into 183 loud failures and lowers the headline.

The retain-incumbent design means that does **not** happen here — the incumbent comes back and
nothing regresses — but it also means the lever costs latency and buys nothing on this sample.

**Two molecules are not a rate.** The population arms (201 `MIRROR_MATCH`, 365 `key_equal`, plus a
200-molecule `byte_exact` control) are what settle it. If the rate holds, **Lane 2 is refuted as a
`byte_exact` lever** and its result is a generator-capability finding for a later release rather
than an acceptance bug — a legitimate outcome, and the fourth time this project has ended a lane
by refuting its own plan.

⚠ Note for whoever reads the arms: both lanes hinge on the **same** unmeasured quantity — *does the
pool contain a better conformer at all?* Lane 1 needs an attached one, Lane 2 needs a
correctly-handed one. A near-zero result in either is a statement about construction, not about
the predicate.

## 9. Instruments added this release

- **`tools/attach_return_preflight.py`** — the guard's predicate vs the bucket's verdict, with a
  mandatory control arm. §3.
- **`harvest_measurements.py`** — `ALLOW` extended to `attach_*.json` and `string_exact_*.json`.
  ⚠ The pre-existing entry was the narrow literal `attach_class_audit.json`, which would have
  **silently dropped** `attach_return_preflight.json` — the same too-narrow-allowlist failure the
  `mirror_arm*` comment in that file records one release earlier. An allowlist fails exactly like a
  broken instrument: it prints a plausible total.
