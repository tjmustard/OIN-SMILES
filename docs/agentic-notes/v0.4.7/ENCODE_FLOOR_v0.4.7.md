# The encode floor — v0.4.7 (`swimlane/v047-encfloor`)

Owner: L3-encfloor, based on `release/v0.4.7` @ `49f8d69b`. The gitignored
`spec/handoffs/v0.4.7/PROGRESS-encfloor.md` carries the working narrative; this file is
the git-durable record.

## Why this lane existed

`OIN_ACCEPT_SCORED` collapses most of the generation-side slow tail. The molecules it
cannot reach are the ones whose **bare encode** is the floor: `QIDKUL_comp_0` goes
390 s → 90 s and is still over target because `XYZToSMILES().convert()` on its input
alone costs 46–71 s. No generation-side lever can touch that.

Wall clock is unusable here — this host ran a 5 000-molecule sweep plus three sibling
lanes throughout, load 35–51 on 12 cores. Every claim below is either a **counter**
(`AC2BO_STATS`, call counts, candidates examined), a **sha**, or a **within-one-process
ratio**. Absolute seconds are marked ADVISORY and are never the basis of a conclusion.

---

## 1. Headline: the inherited cost model is right about one regime and inverted about another

`docs/agentic-notes/v0.4.5/ENCODER_PERF_v0.4.5.md` states that `xyz2mol_local.AC2BO` is **99.8 %** of a slow
encode — "not dominant, essentially all of it". Re-measured with the previously
uninstrumented functions included, that is true for exactly one of **three** distinct
regimes, and is off by a factor of ~100 for another.

| regime | who actually pays | `_AC2BO_core` share of encode | example |
|---|---|---|---|
| **R1 — over-cap valence grind** | `nx.max_weight_matching` via `get_UA_pairs` | **99.5 %** | `QIDKIZ_comp_0` |
| **R2 — sub-cap, large `combo_size`, loop exits early** | `_ordered_valences` materialising a product the loop never reads | 37–91 % | `UDITAD_comp_0`, `NOCGAN_comp_0` |
| **R3 — sub-cap, `AC2BO` ~free, cost is the forked resonance** | `lig_checks` → `_resonance_candidates_isolated`, 120 CPU-s per fork | **0.2 – 3.8 %** | `XIRMER_comp_0`, `NAKLET_comp_0`, `HACYEQ_comp_0` |

The two slowest molecules measured in this lane are both **R3**:

| molecule | encode wall (ADVISORY) | `_AC2BO_core` | share | what the rest is |
|---|---:|---:|---:|---|
| `XIRMER_comp_0` | 1235.54 s | 46.69 s | **3.8 %** | **3 forked resonance timeouts = 95.8 %** (§4.1) |
| `NAKLET_comp_0` | 1160.92 s | 10.01 s | **0.9 %** | forked resonance (child observed via `ps`) |
| `HACYEQ_comp_0` | 54.89 s | 0.13 s | **0.2 %** | **1 forked resonance, `ok`, = 90.2 %** (§4) |

> **A slow encode is not evidence of an over-cap ligand** — the v0.4.5 lane already said
> this. Stronger version, measured here: *a slow encode is not evidence of `AC2BO` at all.*
> For the three molecules above, 96.2 % / 99.1 % / 99.8 % of the encode happens outside it.

Reproduce:

```
V=/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python
DS=<dataset>/tmCAT-tmPHOTO_xyz_dataset
PYTHONPATH=$PWD/src $V tools/encfloor_attribute.py --dataset $DS \
    --molecules XIRMER_comp_0,NAKLET_comp_0,HACYEQ_comp_0
```

### 1.1 The three functions nobody had ever measured — all negligible

The brief flagged that `_canonical_atom_permutation` (`xyz2mol_local.py:1204`) and
`_valence_search_is_truncated` (`:1251`) now run on **every** `AC2BO` call because
`OIN_CANONICAL_PERCEPTION` is default-ON, and that neither had appeared in any
measurement. They are now instrumented (`tools/perf_encode_profile.py`,
`tools/encfloor_attribute.py`). **They are not a cost.**

| function | worst case observed across every molecule measured | share of that encode |
|---|---|---|
| `_canonical_atom_permutation` | 0.031 s (`NOCGAN_comp_0`, 3 calls) | 0.1 % |
| `_valence_search_is_truncated` | 0.007 s (`NOCGAN_comp_0`, 3 calls) | 0.0 % |
| `possible_valences` | 0.006 s (`UDITAD_comp_0`, 7 calls) | 0.4 % |
| `valence_combo_size` | 0.000 s everywhere | 0.0 % |

Concern #1 in the lane brief — that the stale attribution might be hiding cost in the
canonicalising wrapper — is **REFUTED**. `MolToSmiles` on an all-single-bond graph is
cheap even at 220 atoms, and it is called once per `AC2BO`, not once per candidate.

---

## 2. R1 — the over-cap grind (`QIDKUL`/`QIDKIZ`): 99.8 % confirmed, and the leaf is the matching

`tools/valsearch_scan.py` classifies both targets as over-cap: they share one **37-atom**
ligand with `combo_size = 1 259 712` against `_VALENCE_COMBO_CAP = 500 000`.

`QIDKIZ_comp_0`, counters mode:

| | value |
|---|---|
| encode wall (ADVISORY) | 82.30 s |
| `AC2BO` inclusive | 81.89 s (**99.5 %**) |
| `get_UA_pairs` inclusive | 50.03 s (60.8 %), 180 002 calls |
| `nx.max_weight_matching` inclusive | **38.41 s (46.7 %)**, 26 667 calls |
| `get_BO` / `BO_is_OK` / `charge_is_OK` | 12.01 / 10.12 / 11.37 s |
| `_ordered_valences` | 1 call, **0.00 s** |
| `AC2BO_STATS.candidates` | 60 001 (60 000 over-cap) |
| `AC2BO_STATS.found_valid` | **0** |
| `AC2BO_STATS.over_cap_exhausted` | 3 |

So R1 is: three `AC2BO` calls each grind the full 20 000-candidate budget, **find nothing**,
and return `best_BO`. The single largest leaf is `nx.max_weight_matching`. Note this
means the *product* of R1 is a guess, not a perceived Lewis structure — 60 000 candidates
bought `found_valid = 0`.

**LEAD 1 cannot help R1**: `_ordered_valences` is only reached on the sub-cap branch, and
QIDKIZ's expensive fragment never takes it (1 call, 0.00 s — that call is the *other*,
trivial fragment).

---

## 3. R2 — LEAD 1: the sub-cap branch materialises a product it does not read

### 3.1 The population question, answered first

The brief warned that `_ordered_valences` is "bounded only by `_VALENCE_COMBO_CAP =
500 000`". In the real corpus it never comes close:

`$V tools/valsearch_scan.py --dataset $DS --n 500 --seed 47` — 499 molecules,
**1559 ligand fragments**:

| percentile | `combo_size` |
|---|---:|
| p50 | **3** |
| p75 | 8 |
| p90 | 24 |
| p95 | 72 |
| p99 | 432 |
| p99.9 | 11 664 |
| max | **20 736** |

`over _VALENCE_COMBO_CAP: 0/1559`. The cap is **24×** above the largest value the corpus
produced. So the unbounded-memory reading of `_ordered_valences` is theoretical; the real
question is dead work, not blow-up.

### 3.2 The dead work, measured exactly

`_ordered_valences` materialises `combo_size` candidates (twice, plus a dict over the
five-group product). The loop consuming it exits as soon as a valid Lewis structure is
found. **`combo_size − candidates_examined` is exact and load-independent.**

| molecule | encode s (ADV) | `_AC2BO_core` share | materialised | examined | **dead** | dead % | `_ordered_valences` % of encode |
|---|---:|---:|---:|---:|---:|---:|---:|
| `NOCGAN_comp_0` | 53.85 | 76.3 % | 45 568 | 24 833 | **20 735** | 45.5 % | **8.22 %** |
| `UDITAD_comp_0` | 1.50 | 37.5 % | 11 865 | 24 | **11 841** | **99.8 %** | **34.67 %** |
| `EYECAO_comp_0` | 6.43 | 84.0 % | 11 664 | 1 297 | 10 367 | 88.9 % | 10.36 % |
| `AKIKOV_comp_0` | 13.91 | 91.2 % | 5 184 | 2 309 | 2 875 | 55.5 % | 2.08 % |
| `HACYEQ_comp_0` | 54.89 | 0.2 % | 2 306 | 10 | 2 296 | **99.6 %** | 0.15 % |
| `VORREA_comp_0` | 3.50 | 58.2 % | 2 592 | 305 | 2 287 | 88.2 % | 6.11 % |
| `NAKLET_comp_0` | 1160.92 | 0.9 % | 1 538 | 1 537 | 1 | 0.1 % | 0.01 % |
| `XIRMER_comp_0` | 1235.54 | 3.8 % | 12 288 | 12 288 | **0** | 0.0 % | 0.06 % |
| `LEZWAO_comp_0` | 0.95 | 6.0 % | 24 | 24 | 0 | 0.0 % | 0.08 % |

The extreme single call: `NOCGAN_comp_0`'s third `AC2BO` materialised **20 736**
candidates and the loop consumed **1** — 1.956 s of that call's 1.972 s, **99.2 %**.
`UDITAD_comp_0`'s seventh: 11 664 materialised, **9** consumed, **94.2 %**.

### 3.3 The change

`xyz2mol_local.py` sub-cap branch: `_ordered_valences(...)` → `iter_ordered_valences(...)`.

Byte-identity does not rest on a code reading. `iter_ordered_valences` yields the identical
sequence, and that is an **equality already asserted by a test that predates this lane**
(`tests/unit/test_valence_order.py::test_lazy_order_is_identical_to_the_sorted_path`).
Given identical order, both loop exits are identical:

* early return — same first valid candidate, so the same `BO`;
* exhaustion — same candidates in the same order, so the `BO.sum() >= best_BO.sum()`
  tie-break (note `>=`: the **last** maximum wins) sees the same sequence and selects the
  same `best_BO`.

### 3.4 Extending the proof to the population the corpus actually emits

The existing equality test generates **1–7 atom** configurations with `combo_size ≤ 5000`.
Real ligand fragments are not that shape, so the equality was re-run on harvested corpus
fragments — 15 slow / high-`combo_size` molecules, **50 fragments** across both
`allow_carbenes` arms, including `KEMTED` (168 atoms) and `XIRMER` (108 atoms):

```
checked=50 mismatch=0
```

That harvest also characterised the real shape, which is what the new test now generates:

* fragments are 2–168 atoms, but only **1–12** atoms have more than one candidate valence
  — which is why a 168-atom fragment (`KEMTED`) has `combo_size = 16`;
* every per-atom option list has length **1, 2 or 3**, never more;
* in every fragment measured, the multi-option atoms were **C/N/O/P/S** — exactly the
  elements the heuristic groups.

That last point matters for what the test must cover. When every multi-option atom is
grouped, `order_idx` pins the whole candidate and `sorted()`'s lexicographic tie-break is
**never exercised**. So corpus shape and tie-break shape are different risks, and
`test_lazy_order_matches_on_CORPUS_SHAPED_configurations` generates both — its
`forced_tiebreak` arm puts the options on a non-grouped element (As/Cl/Se) at large `n`,
the case the corpus does not currently show but the code still has to get right.

`test_subcap_AC2BO_does_not_materialise_the_product` is a structural pin: it patches
`_ordered_valences` to raise and runs a sub-cap `_AC2BO_core`, so a silent revert to the
eager call fails the suite rather than quietly restoring the dead work.

### 3.5 Byte-identity A/B

Two arms, **sequential**, never in one interpreter, the single file swapped via
`git show HEAD:src/oinsmiles/utils/xyz2mol_local.py` with an EXIT trap (never `git stash` —
the index is shared across 30+ worktrees). The driver `grep`s the base file for the eager
call and aborts if it is absent, so a mis-swapped base cannot produce a passing A/B.

**Gate** — `tools/enc_byte_identity_ab.py`, the project's standard 24-molecule instrument,
which clears the `AC2BO` memo between molecules so a cross-molecule cache hit cannot be
what makes two revisions agree:

| arm | molecules | missing | `MANIFEST_SHA256` |
|---|---:|---|---|
| A (base) | 24 | `[]` | `13b9d5b9e9bc63dae3c554b087a0ecb033f797ce3fd363ed7a16aa212e2f596c` |
| B (mine) | 24 | `[]` | `13b9d5b9e9bc63dae3c554b087a0ecb033f797ce3fd363ed7a16aa212e2f596c` |

All 24 manifest lines are byte-identical, not merely the digest.

**Exercise set** — the gate's 24 are mostly low-`combo_size`, so they would pass without
the changed line doing anything different. These 8 are the molecules with **measured
non-zero dead work**, run through `tools/encfloor_attribute.py` on both arms:

| molecule | OIN sha equal | `AC2BO_STATS` equal | `_ordered_valences` calls A → B |
|---|---|---|---|
| `AGUFEN_comp_0` | yes | yes | 27 → **0** |
| `NOCGAN_comp_0` | yes | yes | 3 → **0** |
| `UDITAD_comp_0` | yes | yes | 7 → **0** |
| `EYECAO_comp_0` | yes | yes | 1 → **0** |
| `AKIKOV_comp_0` | yes | yes | 1 → **0** |
| `VORREA_comp_0` | yes | yes | 1 → **0** |
| `HACYEQ_comp_0` | yes | yes | 2 → **0** |
| `LEZWAO_comp_0` | yes (both raise) | yes | 12 → **0** |

`AC2BO_STATS` equality covers `candidates`, `matching_calls`, `found_valid`,
`ac2bo_calls`, `over_cap_calls`, `over_cap_candidates`, `over_cap_exhausted`,
`over_cap_best_bo_improved` — all exact and load-independent. **Identical candidate counts
are the stronger claim**: they say the two arms examined the same number of candidates, so
the shas do not agree by accident downstream of a different search.

`LEZWAO_comp_0` raises `OINEncodeError` on both arms with the **same message** (quinoid
ring at the same atom indices). An error is part of the contract; it must be the same error.

The `_ordered_valences` column is the proof the change took effect at all — 54 calls
across the 8 molecules became 0, and **50 407 candidate tuples that were materialised and
never examined are no longer built.**

Full suite after the change: `python -m unittest discover tests/unit` →
`Ran 859 tests` `OK (skipped=3, expected failures=5)`.

### 3.6 Honest scope

**This change does nothing for the molecules this lane was pointed at.** `QIDKUL`/`QIDKIZ`
are R1 (over-cap, never reaches the branch); `XIRMER` and `NAKLET` are R3 with **0** dead
candidates. It helps R2, which from the 500-molecule scan is roughly **2 % of molecules**
(10/499 have `max_combo_size ≥ 1000`; 3/499 have ≥ 4096) — and it costs the other 98 %
nothing, because at `combo_size = 3` there is nothing to materialise either way.

---

## 4. R3 — the real floor for the slow sub-cap cohort: the forked resonance

`XIRMER_comp_0` is the molecule the v0.4.5 lane flagged as "sub-cap, >35 minutes, forks a
child, cause never identified". Identified:

`lig_checks` (`xyz2mol.py`) calls `_resonance_candidates_isolated` whenever
`_resonance_needs_isolation` holds (≥ 50 heavy atoms, or ≥ 35 aromatic atoms). That forks
a child bounded by `resource.RLIMIT_CPU` at `_RESONANCE_CPU_BUDGET_S = 120` **CPU**
seconds, with `_RESONANCE_WALL_SAFETY_S = 900` s as a wall-clock backstop. `lig_checks` is
called once per successful `AC2mol` on the `_select_lig_mol` charge/carbene ladder **and**
again from `_rescue_unusable_perception`. So the floor for a large conjugated ligand whose
resonance enumeration does not terminate is

> **`k × 120` CPU seconds**, where `k` is the number of ladder rungs that reach `lig_checks`

— entirely independent of valence combinatorics.

Counter evidence (`tools/encfloor_attribute.py` now instruments `lig_checks` and
`_resonance_candidates_isolated`; `fork_wall` is inclusive of the child, measured in the
same process as the encode it is divided by):

| molecule | encode s (ADV) | `lig_checks` calls | forks | fork wall | **% of encode** | status |
|---|---:|---:|---:|---:|---:|---|
| `HACYEQ_comp_0` | 32.41 | 2 | **1** | 29.22 s | **90.2 %** | `{'ok': 1}` |
| `NOCGAN_comp_0` | 27.67 | 3 | 3 | 4.78 s | 17.3 % | `{'ok': 3}` |
| `AGUFEN_comp_0` | 1.56 | 15 | 0 | — | 0 % | — |
| `LEZWAO_comp_0` | 0.40 | 12 | 0 | — | 0 % | — |
| `UDITAD_comp_0` | 0.47 | 6 | 0 | — | 0 % | — |

Two things worth noticing. First, **`HACYEQ`'s fork returns `ok`, not `timeout`** — the
child is not hitting the 120 CPU-second wall, it simply takes 29 of the encode's 32
seconds to enumerate. So R3 is not only about hangs; a *completing* `ResonanceMolSupplier`
on a large conjugated ligand is itself the floor. Second, `lig_checks` is called **15**
times on `AGUFEN` and **12** on `LEZWAO` with zero forks — those ligands sit below
`_resonance_needs_isolation`'s thresholds, so the ladder's repetition is only expensive
once a fragment is large enough to be isolated.

### 4.1 `XIRMER_comp_0` — the cause, finally

Re-run with the fork counters:

```
FORKED RESONANCE: lig_checks=3 forks=3 status={'timeout': 3} wall=988.18s (95.8% of encode)
_AC2BO_core calls=3  sum(wall_core)=39.90s (3.9% of encode)
_ordered_valences total = 0.000s (0 calls)
```

**Three ladder rungs, three forks, all three `timeout`** — each burns the full 120 CPU-second
budget and returns nothing, so the caller degrades to the single perceived form each time.
988.18 s of a 1031.40 s encode: **95.8 %**. That is the answer to "sub-cap, >35 minutes,
forks a child, cause never identified".

And the three forks are doing the *same* work on the *same* fragment: this molecule has one
ligand fragment (`_canonical_atom_permutation` reports **1 distinct** AC across all three
`AC2BO` calls), and the rungs differ only in `(charge, allow_carbenes)` —
`(-2, True)`, `(-2, False)`, `(-4, True)`. The encode pays `3 × 120` CPU seconds to learn
the same thing three times.

**The obvious fix and its catch.** Caching "this fragment's resonance enumeration times
out" within one encode would cut XIRMER's floor by ~2/3. The catch is that `lig_checks`
receives the *perceived* `lig_mol`, whose bond orders differ per rung — so a key on the
fragment's connectivity is not obviously sound (a different resonance starting point could
in principle terminate where another does not), while a key on the fully-perceived mol may
simply never hit. **Which of those is true is a measurement nobody has made**, and it is
the single highest-value next step for the encode floor. It is a correctness question, not
a caching one, and it should not be shipped on the strength of the argument above.

Also worth noting from this run: post-change, XIRMER's `AC2BO_STATS` are
`candidates=12288 matching=4096 found_valid=0` — **identical** to the pre-change run, and
`ordered_len` now reports `-` (the lazy path), giving a byte-identity data point on a
17-minute molecule for free.

For `NAKLET_comp_0` the fork counters were not captured; the evidence there is `ps` showing
a forked child several minutes old plus `_AC2BO_core` at 0.9 % of its ~19-minute encode.

The CPU bound is deliberate and correct (it makes the *outcome* load-independent), but it
means R3's wall cost under contention is `k × 120` CPU-seconds stretched by the load
factor — at load 45 on 12 cores, a 120 CPU-second child can take ~8 wall minutes.

**This, not `AC2BO`, is where the remaining encode floor work is.** Candidate directions,
none measured by this lane and none of them free:

1. **Cache the `timeout` verdict within one encode** — measured to be worth ~2/3 of
   `XIRMER`'s entire encode (§4.1). Recommended next probe, with the soundness question
   spelled out there. Measure first whether the perceived `lig_mol` actually repeats across
   rungs; if it does, the key is sound and the win is free.
2. Memoize the enumeration *result* (not just the verdict) across ladder rungs — same
   distinctness question, larger payoff on `HACYEQ`-like molecules whose forks return `ok`.
3. Lower `_RESONANCE_CPU_BUDGET_S`. **Not byte-identical** — a child that would have
   completed at 119 s now returns `timeout` and the caller degrades to the single perceived
   form, which can change the emitted OIN. Would need its own fidelity A/B. Note `HACYEQ`'s
   fork returns `ok` after ~29 s, so a budget cut aimed at `XIRMER` would risk exactly the
   molecules that currently succeed.

---

## 5. LEAD 2 — ladder invariants ARE recomputed, and it does not matter

The redundancy is real and larger than the brief estimated. `AGUFEN_comp_0` issues **27**
`_AC2BO_core` calls, not the "13 worst case" — the ladder plus `_rescue_unusable_perception`
plus `get_tmc_mol`'s `suppress_canonical_perception()` retry (AGUFEN is the documented case
where canonical perception yields a pentavalent carbon):

| helper | calls | **distinct** results | wall |
|---|---:|---:|---:|
| `possible_valences` | 27 | **3** | 0.002 s |
| `valence_combo_size` | 27 | 2 | 0.000 s |
| `_canonical_atom_permutation` | 3 | **1** | 0.003 s |
| `_valence_search_is_truncated` | 3 | 1 | 0.000 s |

**Verdict: NEGATIVE — do not memoize.** The redundancy ratio is 9:1 but the total
addressable saving is under 0.05 s on any molecule measured, including a 1235-second one.
Memoizing would add a cache, a key, an invalidation question and a test, to remove
something that does not appear above the noise floor of a contended host. Recorded here so
the 9:1 ratio does not tempt a future lane into re-deriving the same non-result.

---

## 6. LEAD 3 — the `AC2BO` memo DOES hit across a renumbering now, but not where this lane's targets live

`docs/agentic-notes/v0.4.5/V045_STATUS_2026-07-25.md` concluded the memo "cannot hit on the generated side,
because its key is a permutation away". That measurement predates the canonicalising
wrapper: `AC2BO` now calls `_AC2BO_core` on the **permuted** matrix, and
`_ac2bo_memo_anchor` tags on the bytes of *that* matrix, so two encodes of the same graph
should share a tag.

Tested **generator-free**, with `tools/encfloor_memo_probe.py`. Renumbering the input XYZ
is the right instrument because it holds the graph *exactly* fixed while destroying atom
order — a generated conformer does not (its perceived graph matched the input only 16/36 =
44 %, and a sibling lane has since measured that under `OIN_ACCEPT_SCORED` six molecules'
generated structures no longer perceive the haptic ligand as coordinated at all). So this
is an **upper bound on the mechanism**, measured without spending generator time.

Two arms per molecule in the same process: `REUSE` encodes the original then the
renumbered copy without clearing the memo; `COLD` clears between. The payload is
`AC2BO_STATS['matching_calls']` on **encode 2** — exact and load-independent.

| molecule | cap | encode-2 matching, REUSE | encode-2 matching, COLD | verdict | OIN equal |
|---|---|---:|---:|---|---|
| `VORREA_comp_0` | sub | **0** | 305 | **HIT** | yes |
| `EYECAO_comp_0` | sub | **0** | 1 297 | **HIT** | yes |
| `AKIKOV_comp_0` | sub | **0** | 2 309 | **HIT** | yes |
| `NOCGAN_comp_0` | sub | **0** | 1 296 | **HIT** | yes |
| `QIDKIZ_comp_0` | **over** | 26 642 | 26 642 | **MISS** | yes |

**Sub-cap: the memo eliminates 100 % of `nx.max_weight_matching` on the second encode.**
Not "some" reuse — zero matching calls. The v0.4.5 conclusion is refuted for the sub-cap
population under canonical perception.

**Over-cap: it cannot hit, and this is structural.** `AC2BO` routes an over-cap fragment to
`plain()` on the **un-permuted** AC (`_valence_search_is_truncated` returns `True`, which
is deliberate — relabelling can move a valid assignment out of the truncated search
window). So the tag is a function of the input numbering and a renumbering guarantees a
miss. `QIDKIZ_comp_0` shows the underlying order-dependence directly: the renumbered copy
examines **59 239** candidates where the original examined **60 001**. The emitted OIN
still matched, but that is `best_BO` happening to serialize the same, not invariance.

Two further facts worth carrying forward: one capped over-cap search produced **~29 000**
memo entries (`entries=28715` on the COLD arm), confirming the brief's estimate against the
200 000-entry cap; and no eviction occurred at 59 386 entries on the REUSE arm.

**Consequence for this lane's targets: nothing.** `QIDKUL`/`QIDKIZ` are exactly the case
that cannot hit. The value of LEAD 3 is for the sub-cap round trip, and it is *already
delivered by code on `main`* — this lane's contribution is measuring that it works, and
identifying that extending it to over-cap requires giving the over-cap branch a
renumbering-invariant anchor, which is a correctness question (which candidates the
truncated window contains), not a caching one.

## 7. Ruled OUT

Inherited from v0.4.5 and **not** re-chased (each already has numbers): extending
`OIN_VALENCE_CHARGE_FILTER` to sub-cap; numpy-vectorising `charge_is_OK` / `get_UA` /
`check_sum`; `OIN_VALENCE_MATCHER=maxcard` or `greedy`;
`OIN_VALENCE_ORDERED_FALLBACK`; memoizing whole-`AC2BO` results; a `get_BO` memo.

Ruled out **by this lane**:

- **The canonicalising wrapper is a hidden cost.** Refuted — §1.1, ≤ 0.1 % everywhere.
- **`_ordered_valences` is a general encoder hot spot.** Refuted — §3.1: median
  `combo_size` is 3, and on the two slowest molecules measured it is 0.06 % and 0.01 %
  of the encode with **zero** dead candidates.
- **Memoizing the ladder's invariants (LEAD 2).** Refuted — §5, under 0.05 s addressable.
- **`AC2BO` is the encode floor.** Refuted as a general claim — §1. It is the floor for
  R1 only.

---

## 8. Tools added

| tool | what it answers |
|---|---|
| `tools/encfloor_attribute.py` | per-`AC2BO`-call attribution: `combo_size`, candidates actually consumed, matching calls, `_ordered_valences` share, plus the forked-resonance counters. Overhead is per-`AC2BO`-call, not per-candidate, so it is safe on a 20-minute molecule where `perf_encode_profile.py`'s per-candidate wrappers would distort. |
| `tools/encfloor_memo_probe.py` | generator-free upper bound on whether the `AC2BO` memo can hit across two differently-numbered encodes of one graph (§8). |
| `tools/perf_encode_profile.py` | extended with `_canonical_atom_permutation`, `_valence_search_is_truncated`, `_ordered_valences`, `possible_valences`, `valence_combo_size`, `_AC2BO_core`, and an `AC2BO_STATS` dump. |
