# The valence search — v0.4.5 (`swimlane/v045-valsearch`)

Owner: valsearch swimlane, based on `main` @ `ebd6aabc`. The git-durable record; the
gitignored `spec/handoffs/v0.4.5/PROGRESS-valsearch.md` has the working narrative.

## Why this lane existed

`docs/ENCODER_PERF_v0.4.5.md` (the `encspeed` lane) profiled the encoder for the first time
and found `AC2BO` is **99.8%** of a slow encode. It then memoized the charge-*independent*
half of the search value-identically (1.62–2.16x) and stopped, because what remained —
`BO_is_OK` + `charge_is_OK` + `get_BO` over ~60 000 candidates, ~13–15 s — is exactly the part
that differs per charge/carbene ladder arm and so cannot be memoized. It named that the floor
and declined to attack it, correctly, because every remaining lever **changes perceived bond
orders**.

It also recorded the observation this lane exists to test: for these ligands the search
**fails**. Each ladder arm exceeds `_VALENCE_COMBO_CAP` (500 000), grinds all
`_VALENCE_FALLBACK_TRIES` = 20 000 candidates, finds **none valid**, and returns `best_BO` — a
fallback guess. *"The 20 000-candidate search is buying a guess. Whether that guess is worth
15 s per ladder arm is a product question this lane surfaces rather than answers."*

This is the answer.

## Base correction, stated up front

**`swimlane/v045-encspeed` is not merged into `main`** — `git merge-base --is-ancestor
3b15e26c main` is false. This lane is based on `main` @ `ebd6aabc`, so its baseline has **no
memoization**. The "~13–15 s per arm" figure is a *post-memo* residual; on `main` the same
path is ~35–50 s per arm. Every speedup below is measured against `main` and therefore
**overlaps the encspeed win rather than composing with it** — both remove work from the same
loop. Do not add the two ratios.

## Headline

**The 20 000-candidate search is not merely buying a guess — on every molecule measured it is
buying the *same* guess it already had after a few hundred candidates.** `best_BO` reaches its
final value early and is then overwritten hundreds of times with an identical matrix, because
the update predicate is `BO.sum() >= best_BO.sum()` (note `>=`, so ties overwrite) and the
maximum attainable sum is found near the front of the lazy product order.

## Q1 — how often does the search actually fail?

The over-cap decision is a function of the adjacency matrix alone, so it can be measured
**without running the search**: `tools/valsearch_scan.py` does `get_basic_mol` +
`MetalDisconnector` + `GetAdjacencyMatrix` + `possible_valences` at ~0.04 s/molecule instead of
~50 s. That is what made the population question affordable while the release sweep saturated
the host.

| sample | molecules | ligand frags | frags over cap | mols with ≥1 over cap |
|---|---|---|---|---|
| random, seed 1 | 100 | 315 | 0 (0.0%) | **0 (0.0%)** |
| random, seed 7 | 1992 | 6160 | 4 (0.1%) | **4 (0.2%)** |
| slow/failing cohort | 313 | 1005 | 4 (0.4%) | **4 (1.3%)** |

The slow cohort is `spec/handoffs/v0.4.5/hard_fail_worklists.json` (all five buckets) plus the
7 `resonance_timeout` molecules from `docs/ENCODE_FAIL_v0.4.5.md`.

**The rate is not the finding — the enrichment is.** The 4 over-cap molecules in the slow
cohort are `KESWUB`, `BENVOG`, `HICLAG`, `HOHKUL`: **4 of the 7** molecules
`docs/ENCODE_FAIL_v0.4.5.md` lists as *"unresolved this session — genuine encoder-side cost,
contended-machine time budget exhausted"*. That is **57% of that cohort**. The other three
(`FAQYUU`, `KEMTED`, `NAKLET`) are **not** over-cap, so their cost has a different cause and
this lane's lever cannot help them.

So: **~0.2% of the corpus, and a majority of the encoder's known unresolved residue.** A
general-purpose speedup this is not. A targeted fix for the encoder's hardest cohort it may be.

Positive control, because a scanner that reports 0 everywhere may simply be broken: the scan
flags both molecules `encspeed` measured, each at combo_size **1 259 712** > 500 000. It also
shows they share the **same 37-atom ligand** — that lane's "two molecules" are one ligand
measured twice, which is worth knowing before generalising from it.

## Q2 — does the fallback answer depend on how long you search?

### Whole-encode A/B, `QIDKUL_comp_0`

All arms in one process on one input file; budgets run 20000/5000/1000/200 and then **20000
again** as a determinism self-check.

| `_VALENCE_FALLBACK_TRIES` | wall | candidates | `best_BO` sha | OIN sha |
|---|---|---|---|---|
| 20 000 | 148.10 s | 60 002 | `1ea2ed6bafe7` | `0d428d9dfc56` |
| 5 000 | 27.29 s | 15 002 | `1ea2ed6bafe7` | `0d428d9dfc56` |
| 1 000 | 3.50 s | 3 002 | `1ea2ed6bafe7` | `0d428d9dfc56` |
| 200 | **0.92 s** | 602 | `1ea2ed6bafe7` | `0d428d9dfc56` |
| 20 000 (repeat) | 148.13 s | 60 002 | `1ea2ed6bafe7` | `0d428d9dfc56` |

`REPEAT-OK` — the two 20 000 arms agree to 0.03 s and to the sha, so the harness is
deterministic and the other rows are not noise. The OIN sha `0d428d9dfc56` **also matches the
value `docs/ENCODER_PERF_v0.4.5.md` recorded independently** for this molecule, so the encode
is confirmed correct across lanes rather than merely self-consistent.

`best_BO` identical, OIN identical, **161x** less wall clock.

### Whole-encode A/B, `QIDKIZ_comp_0`

| `_VALENCE_FALLBACK_TRIES` | wall | candidates | `best_BO` sha | OIN sha |
|---|---|---|---|---|
| 20 000 | 631.88 s | 60 001 | `1ea2ed6bafe7` | `796221298c6c` |
| 5 000 | 70.03 s | 15 001 | `1ea2ed6bafe7` | `796221298c6c` |
| 1 000 | 14.64 s | 3 001 | `1ea2ed6bafe7` | `796221298c6c` |
| 200 | **2.50 s** | 601 | `1ea2ed6bafe7` | `796221298c6c` |
| 20 000 (repeat) | 334.78 s | 60 001 | `1ea2ed6bafe7` | `796221298c6c` |

`REPEAT-OK`, and the OIN sha `796221298c6c` again matches `docs/ENCODER_PERF_v0.4.5.md`
independently. Note the two 20 000 arms took 631.88 s and 334.78 s for **bit-identical work** —
a 1.9x spread from host contention alone. The shas are the claim; the seconds are weather.

### Why — the mechanism, measured on the ligand directly

The 37-atom over-cap ligand, `charge = -8`, carbenes on:

| tries | candidates | `best_BO` reassignments | `found_valid` | `BO.sum()` | `BO == AC` |
|---|---|---|---|---|---|
| 200 | 200 | **7** | **0** | 94 | False |
| 20 000 | 20 000 | **485** | **0** | 94 | False |

Three things at once:

* **`found_valid = 0`** — the search genuinely fails, confirming the `encspeed` reading. No
  valid Lewis structure exists among 20 000 candidates, so the returned value is a guess.
* **The guess is not the input.** `BO != AC` and `BO.sum()` is 94 vs `AC.sum()` 78, so
  `best_BO` did move off its `AC.copy()` initialisation — the search is not a complete no-op.
* **But it stops improving almost immediately.** The maximum attainable sum, 94, is reached
  inside the first 200 candidates. The 478 further "improvements" are all `>=` **ties that
  write back an identical matrix**. That is why the answer is budget-invariant, and it is an
  empirical property, not a guarantee: a tie carrying a *different* matrix of equal sum would
  overwrite, and nothing in the code prevents that.

### Where the answer saturates — the number that sets the safe budget

Same ligand, same arm, budget swept finely. This is the table that matters, because it shows
the search is **not** useless — it is useful for about 100 candidates and useless for the
remaining 19 900.

| tries | candidates | `best_BO` reassignments | `BO.sum()` | `best_BO` sha | wall |
|---|---|---|---|---|---|
| 1 | 1 | 1 | 92 | `53205ad7836d` | 0.00 s |
| 2 | 2 | 2 | 92 | `53205ad7836d` | 0.00 s |
| 5 | 5 | 3 | 92 | `53205ad7836d` | 0.00 s |
| 10 | 10 | 3 | 92 | `53205ad7836d` | 0.01 s |
| 25 | 25 | 3 | 92 | `53205ad7836d` | 0.02 s |
| 50 | 50 | 4 | **94** | `facbe31d0ae9` | 0.23 s |
| **100** | 100 | 7 | 94 | **`9d234c9103a8`** | 0.33 s |
| 200 | 200 | 7 | 94 | `9d234c9103a8` | 1.20 s |
| 500 | 500 | 15 | 94 | `9d234c9103a8` | 5.04 s |
| 1 000 | 1 000 | 26 | 94 | `9d234c9103a8` | 4.95 s |
| 5 000 | 5 000 | 116 | 94 | `9d234c9103a8` | 109.00 s |
| 20 000 | 20 000 | 485 | 94 | `9d234c9103a8` | 231.53 s |

Read it in three bands:

* **1–25 candidates: genuinely worse.** Sum 92, a different matrix. Cutting this far would
  silently degrade perception. Anyone tempted to "just set it to 10" is wrong, and this row is
  the evidence.
* **50 candidates: right sum, wrong matrix.** `facbe31d0ae9` attains sum 94 but is not the
  answer the full search returns. The sum alone is not a sufficient stopping criterion.
* **100 candidates onward: converged.** `9d234c9103a8` from 100 through 20 000 — **0.5% of the
  budget produces 100% of the answer.** The 485 reassignments at 20 000 versus 7 at 100 are all
  ties rewriting the same matrix.

So the honest shape of the result is not "the grinding is pointless" but **"the grinding
saturates at ~0.5% of the budget"**. 200 is the smallest round number with margin over the
observed convergence point, and it is the value this lane recommends.

Note also that the two 20 000 timings taken minutes apart under different contention were
100.01 s and 231.53 s. That 2.3x spread on identical work is why every claim here is anchored
to counters and shas, with seconds as context only.

## Q3 — is `max_weight_matching` the wrong algorithm?

**No, and the premise of the question is wrong in an instructive way.**

The premise is confirmed as a unit test, not just a note: the graph `get_UA_pairs` builds
carries **no weight attributes** (`test_the_graph_get_ua_pairs_builds_carries_no_weights`), so
`nx.max_weight_matching` really is solving maximum *cardinality* matching.

But the graphs are **tiny**. Capturing the 200 real graphs from a production run of the over-cap
ligand gives 14–17 nodes and 13–17 edges. **The cost was never per-call complexity — it was
call count: 60 001 invocations on 15-node graphs.** "O(n³)-ish is too slow" is the wrong
diagnosis; a faster asymptotic matcher has almost nothing to bite on.

Interleaved microbenchmark, 3 repetitions over those 200 production graphs, arms alternating so
contention averages across them:

| matcher | wall (3 reps) | speedup | total cardinality | matching sha |
|---|---|---|---|---|
| `nx` — `nx.max_weight_matching` | 1.027 s | 1.00x | 1917 | `6a86cfa95a17` |
| `maxcard` — same, `maxcardinality=True` | 5.469 s | **0.19x** | 1917 | `6a86cfa95a17` |
| `greedy` — `nx.maximal_matching` | 0.059 s | **17.46x** | **1901** | `10d0add385b9` |

And the resulting perception:

| matcher | `BO.sum()` | `best_BO` sha | `found_valid` |
|---|---|---|---|
| `nx` | 94 | `9d234c9103a8` | 0 |
| `maxcard` | 94 | `9d234c9103a8` | 0 |
| `greedy` | 94 | **`e0a7c536aaa9`** | 0 |

* **`maxcard` is rejected** — 5.3x *slower* for a bit-identical answer. Asking networkx's Blossom
  implementation for the max-cardinality variant costs more, it does not save.
* **`greedy` is rejected on fidelity.** It is 17x faster but *maximal*, not *maximum*: 1901 edges
  against 1917, and it returns a **different `best_BO`** — the same total bond order (94) via a
  different Kekulé structure. That is precisely the "silently perceives a different molecule"
  trade this release has repeatedly caught, and it buys a component cost that Q2 has already
  eliminated.

**Q3 is moot once Q2 lands, and that is the real answer.** Matching was ~28 s of a 49 s encode
only because it ran 60 001 times; at 200 candidates per arm it runs ~600 times and costs ~0.3 s.
The budget cut and the matcher swap attack the *same* seconds — and the budget cut takes them
without changing a single perceived bond order, while the matcher swap changes them. **The
budget cut strictly dominates.** The `OIN_VALENCE_MATCHER` lever is kept only so this
measurement is reproducible; nothing should ship behind it.

## Gating

Both levers are **default OFF** and `main`'s output is byte-unchanged:

* `OIN_VALENCE_FALLBACK_TRIES` — unset is the historical 20 000. Garbage or non-positive values
  fall back to the default rather than silently collapsing every over-cap perception.
* `OIN_VALENCE_MATCHER` — unset is `nx.max_weight_matching`.

A structural point that makes the blast radius provable rather than sampled: `_fallback_tries()`
is read **only inside `if over_cap:`**. Sub-cap ligands — which Q1 measures as **99.8%** of the
corpus — cannot reach it, so they are byte-identical *by construction*, not by A/B. What needs
measuring is only the over-cap population, and that population is small enough to test
exhaustively rather than sampled.

## Instrument defects found (this release keeps producing them)

* **My own scan under-reported its denominator.** Results were keyed by `path.stem`, and
  `cat/` and `photo/` share refcodes, so 8 of a 2000-molecule sample silently collapsed into
  1992 — shrinking the denominator of the very rate the tool exists to report. Now keyed
  `<subdir>/<stem>`, and it warns if any key still collapses.
* **`docs/ENCODER_PERF_v0.4.5.md`'s two molecules are one ligand.** Both `QIDKUL_comp_0` and
  `QIDKIZ_comp_0` reduce to the same 37-atom over-cap fragment at combo_size 1 259 712. The
  profile agreeing across "two molecules" is therefore weaker evidence of generality than it
  reads as — it is one ligand measured twice.
