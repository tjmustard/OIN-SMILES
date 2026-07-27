# Lane: valence search (`swimlane/v045-valsearch` + `swimlane/v045-valorder`)

**What the lane was for:** `AC2BO` — the routine that turns an adjacency matrix into bond
orders — searches a Cartesian product over per-atom candidate valences, and on large
conjugated ligands that product is exponential and the encoder never returns. This lane
bounded that search, instrumented it, and answered whether the grinding was buying a better
answer or merely a slower one.

> **Scope note, read this first.** The retrospective brief lists
> `docs/agentic-notes/v0.4.5/VALENCE_ORDER_v0.4.5.md` as a source for the *other* lane in this pair. It is not:
> that document is about the order in which `AC2BO` **enumerates candidate valence
> assignments**, which is the same search this file covers, and its findings are folded in
> here. The lane documented in `LANE-valence-order.md` is a different order-dependence
> entirely — the atom-index order of the AC **valence-capping loop** in `xyz2AC_obabel`,
> fixed by `OIN_STABLE_METAL_AC`. Two lanes, two meanings of "order"; see the discrepancy
> note at the top of `LANE-valence-order.md`.

Primary sources: `docs/agentic-notes/v0.4.5/VALENCE_SEARCH_v0.4.5.md`, `docs/agentic-notes/v0.4.5/VALENCE_ORDER_v0.4.5.md`,
`docs/agentic-notes/v0.4.5/ENCODER_PERF_v0.4.5.md`, `docs/agentic-notes/v0.4.5/ENCODE_FAIL_v0.4.5.md`,
`src/oinsmiles/utils/xyz2mol_local.py`, `tests/unit/test_valence_search_levers.py`,
`tests/unit/test_valence_order.py`, `tools/valsearch_scan.py`,
`tools/valsearch_budget_ab.py`, `tools/valorder_probe.py`,
`tools/valorder_feasibility.py`, `tools/valorder_encode_ab.py`.

---

## ELI5

A 3D structure file gives you where the atoms are, not which bonds are double or triple. The
encoder works that out by guessing how many bonds each atom should have (its *valence*),
trying every combination of those guesses, and keeping the first combination that produces a
chemically legal molecule (a *Lewis structure*). The trouble is that combinations multiply: 37
atoms with a few choices each is 1.26 million combinations, and a 147-atom ligand is 53.7
million — so on big flat conjugated ligands the encoder would sit there forever and produce no
string at all. This lane put a hard ceiling on how much of that space gets searched, then
measured what the search was actually achieving in the part it did reach. The answer was
uncomfortable: for the hardest molecules it was examining 20,000 combinations, finding *zero*
legal ones, and returning a guess — and the reason was not that the ceiling was too low, but
that legal combinations are vanishingly rare and sit far past any reachable ceiling. Counting
them in closed form instead of sampling them turned one molecule's 124-second failed encode
into a 0.87-second correct one.

---

## The work, visually

```
                 XYZ coordinates
                        |
                        v
   xyz2AC_obabel  ->  AC (adjacency matrix)      atoms = [C, C, N, O, B, ...]
                        |
                        v
   possible_valences(AC_valence, atoms, allow_carbenes)      [xyz2mol_local.py:884]
        atom0: [4]     atom1: [4,2,3]     atom2: [3,2]   ...   atom36: [2,1]
                       \______ per-atom candidate valence lists (vll) ______/
                        |
                        v
   valence_combo_size(vll)  =  PRODUCT of the list lengths     [xyz2mol_local.py:920]
        cisplatin / ferrocene scale ligand ......... a few dozen
        QIDKUL / QIDKIZ  37-atom ligand ............ 1 259 712
        ZAZREZ          144-atom ligand ............ 8 503 056
        HICLAG          147-atom ligand ........... 53 747 712
                        |                            (computable WITHOUT running the
                        |                             search -> tools/valsearch_scan.py,
                        |                             ~0.04 s/molecule instead of ~50 s)
                        |
      +-----------------+------------------------------------------------+
      |                                                                  |
  combo <= _VALENCE_COMBO_CAP (500 000)              combo > _VALENCE_COMBO_CAP
  SUB-CAP  ~99.8% of corpus                          OVER-CAP  ~0.2% of corpus
      |                                                                  |
      v                                                                  v
  _ordered_valences(vll, atoms)                       candidate_source = one of:
  materialises the WHOLE product TWICE                  (default) itertools.product(*vll)
  (once to score, once to sort), groups                 [O] OIN_VALENCE_ORDERED_FALLBACK
  O -> N -> C -> P -> S, walks it.                          iter_ordered_valences(...)
  Valid structure found in 3-8 candidates.              [F] OIN_VALENCE_CHARGE_FILTER
      |                                                     iter_charge_feasible_valences(...)
      |                                                                  |
      |                                                                  v
      |                                          islice(candidate_source, _fallback_tries())
      |                                                    <-- [B] OIN_VALENCE_FALLBACK_TRIES
      |                                                        default _VALENCE_FALLBACK_TRIES
      |                                                        = 20 000
      |                                                                  |
      |     QIDKUL's 1 259 712-candidate space, drawn to scale in RANK:  |
      |                                                                  |
      |     |<-- 20 000 -->|                                             |
      |     [##############|.............................................]
      |      ^ the window the default searches                           |
      |                            ^ rank 209 858 = FIRST candidate that could
      |                              possibly be valid.  ONLY 16 EXIST in 1.26 M.
      |                              No prefix budget in ANY order reaches it.
      |                                                                  |
      +--------------------------+---------------------------------------+
                                 |
                                 v
        per candidate:  get_UA -> get_UA_pairs -> nx.max_weight_matching
                                -> get_BO -> BO_is_OK(BO, AC, charge, ...) ?
                                 |                        |
                    VALID -------+                        +------- NOT VALID
                      |                                              |
                      v                          if BO.sum() >= best_BO.sum()
              return BO   (found_valid++)          and valences_not_too_large
                                                   and charge_is_OK:
                                                        best_BO = BO.copy()
                                                   ^^ note ">=" : TIES OVERWRITE
                                                              |
                       window exhausted with nothing valid    v
                       ---------------------------------> return best_BO
                                                          = A FALLBACK GUESS
                                                            (over_cap_exhausted++)

LEGEND
  AC .................. adjacency matrix: 1 where two atoms are bonded, distance-derived
  BO .................. bond-order matrix: 1/2/3 per bond, what the search is looking for
  vll ................. valences_list_of_lists, the per-atom candidate valence lists
  sub-cap / over-cap .. below/above _VALENCE_COMBO_CAP; decided from (AC, atoms) alone
  best_BO ............. the highest-total-bond-order candidate seen; returned when the
                        search fails. NOT a validated Lewis structure.
  found_valid ......... counter: AC2BO calls that returned a structure passing BO_is_OK
  [B] [O] [F] ......... the three default-OFF env levers this lane shipped
  #### / .... ......... searched / unsearched region of the candidate space, by rank
```

---

## Initial assumptions and hypothesis

The lane inherited three beliefs, two from `docs/agentic-notes/v0.4.5/ENCODER_PERF_v0.4.5.md` (the `encspeed`
lane) and one it generated itself.

1. **"`AC2BO` is the encoder."** `encspeed` profiled the encoder for the first time and found
   `AC2BO` is **99.8%** of a slow encode. It memoized the charge-*independent* half
   value-identically (1.62–2.16x) and then stopped, because the remainder — `BO_is_OK` +
   `charge_is_OK` + `get_BO` over ~60,000 candidates, ~13–15 s — differs per charge/carbene
   ladder arm and cannot be memoized. It named that the floor and declined to attack it,
   correctly, because every remaining lever **changes perceived bond orders**.

2. **"The 20,000-candidate search is buying a guess."** `encspeed` observed that each ladder
   arm exceeds `_VALENCE_COMBO_CAP`, grinds all `_VALENCE_FALLBACK_TRIES = 20 000`
   candidates, finds none valid, and returns `best_BO`. Its closing question — *"Whether that
   guess is worth 15 s per ladder arm is a product question this lane surfaces rather than
   answers"* — is what this lane existed to answer.

3. **The ordering hypothesis (generated by the valsearch half, tested by the valorder half).**
   The over-cap branch does not merely cap the search; it *skips* `_ordered_valences`, the
   O/N/C/P/S grouping heuristic. Sub-cap ligands find a valid structure in 3–8 candidates
   *because of* that heuristic; over-cap ligands iterate the raw unsorted product and fail in
   20,000. So `found_valid = 0` might be a statement about search **order**, not about the
   ligand — and if so the right fix was applying the heuristic to a bounded prefix, which
   would be an *accuracy* result rather than a performance one.

A fourth, unstated assumption was live in the codebase and turned out to matter: that
`BO.sum()` — total bond order — tracks perception quality, since `best_BO` is selected by it.

**Base correction that governs every number below.** `swimlane/v045-encspeed` was **not**
merged into `main` when this lane forked (`git merge-base --is-ancestor 3b15e26c main` is
false). The lane is based on `main` @ `ebd6aabc`, so its baseline has **no memoization**: the
"~13–15 s per arm" figure is a *post-memo* residual, and on `main` the same path is ~35–50 s
per arm. Every speedup here **overlaps** the encspeed win rather than composing with it —
both remove work from the same loop. **Do not add the two ratios.**

---

## What was actually found

### CONFIRMED — the search fails, and the population is tiny but enriched

`tools/valsearch_scan.py` answers the population question without running the search, because
the over-cap decision is a function of the adjacency matrix alone (`get_basic_mol` +
`MetalDisconnector` + `GetAdjacencyMatrix` + `possible_valences`, ~0.04 s/molecule):

| sample | molecules | ligand frags | frags over cap | mols with ≥1 over cap |
|---|---|---|---|---|
| random, seed 1 | 100 | 315 | 0 (0.0%) | **0 (0.0%)** |
| random, seed 7 | 1992 | 6160 | 4 (0.1%) | **4 (0.2%)** |
| slow/failing cohort | 313 | 1005 | 4 (0.4%) | **4 (1.3%)** |

The slow cohort is `spec/handoffs/v0.4.5/hard_fail_worklists.json` (all five buckets) plus the
7 `resonance_timeout` molecules from `docs/agentic-notes/v0.4.5/ENCODE_FAIL_v0.4.5.md`.

**The rate is not the finding — the enrichment is.** The 4 over-cap molecules in the slow
cohort are `KESWUB`, `BENVOG`, `HICLAG`, `HOHKUL`: **4 of the 7** molecules
`docs/agentic-notes/v0.4.5/ENCODE_FAIL_v0.4.5.md` classifies as `resonance_timeout` and describes as *"unresolved
this session — genuine encoder-side cost, contended-machine time budget exhausted"*. That is
**57%** of that bucket. `FAQYUU`, `KEMTED` and `NAKLET` are **not** over-cap, so their cost
has a different cause and this lane's levers cannot help them.

Context for that cohort: `encode_fail` is 48 molecules — 11.0% of the 436-molecule gap on the
capstone arm (6,719 molecules) — where the encoder produced **no OIN string at all** for the
input. Its final re-triage is `boron_cluster` 34, `resonance_timeout` 7, `encodes_now` 3,
quinoid/ylide kekulize 3, `perception_charge_gap` 1. The exponential valence product is the
mechanism behind the timeout part of that cohort.

Positive control, because a scanner reporting 0 everywhere may simply be broken: the scan
flags both molecules `encspeed` profiled, each at `combo_size` **1 259 712** > 500 000.

### CONFIRMED — the budget is 99.5% waste on the class it was measured on

Whole-encode A/B, all arms in one process on one input file, budgets 20000/5000/1000/200 and
then **20000 again** as a determinism self-check (`tools/valsearch_budget_ab.py`):

`QIDKUL_comp_0`

| `_VALENCE_FALLBACK_TRIES` | wall | candidates | `best_BO` sha | OIN sha |
|---|---|---|---|---|
| 20 000 | 148.10 s | 60 002 | `1ea2ed6bafe7` | `0d428d9dfc56` |
| 5 000 | 27.29 s | 15 002 | `1ea2ed6bafe7` | `0d428d9dfc56` |
| 1 000 | 3.50 s | 3 002 | `1ea2ed6bafe7` | `0d428d9dfc56` |
| 200 | **0.92 s** | 602 | `1ea2ed6bafe7` | `0d428d9dfc56` |
| 20 000 (repeat) | 148.13 s | 60 002 | `1ea2ed6bafe7` | `0d428d9dfc56` |

`QIDKIZ_comp_0`

| `_VALENCE_FALLBACK_TRIES` | wall | candidates | `best_BO` sha | OIN sha |
|---|---|---|---|---|
| 20 000 | 631.88 s | 60 001 | `1ea2ed6bafe7` | `796221298c6c` |
| 5 000 | 70.03 s | 15 001 | `1ea2ed6bafe7` | `796221298c6c` |
| 1 000 | 14.64 s | 3 001 | `1ea2ed6bafe7` | `796221298c6c` |
| 200 | **2.50 s** | 601 | `1ea2ed6bafe7` | `796221298c6c` |
| 20 000 (repeat) | 334.78 s | 60 001 | `1ea2ed6bafe7` | `796221298c6c` |

`REPEAT-OK` on both. Identical `best_BO`, identical OIN, **161x** and **253x** less wall
clock. The OIN shas independently match the values `docs/agentic-notes/v0.4.5/ENCODER_PERF_v0.4.5.md` recorded for
these molecules, so the encode is confirmed across lanes rather than merely self-consistent.
Note the two 20,000 arms of `QIDKIZ` took 631.88 s and 334.78 s for **bit-identical work** — a
1.9x spread from host contention alone. **The shas are the claim; the seconds are weather.**

### CONFIRMED — where the answer saturates (the number that sets a safe budget)

Same over-cap 37-atom ligand, `charge = -8`, carbenes on, budget swept finely:

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

`best_BO` converges at **100 candidates — 0.5% of the budget produces 100% of the answer**.
The 485 reassignments at 20,000 versus 7 at 100 are all `>=` **ties rewriting an identical
matrix**, which is why the answer is budget-invariant. That is an empirical property, not a
guarantee: nothing in the code prevents a late tie from carrying a *different* equal-sum
matrix, and the 50-candidate row (`facbe31d0ae9`, right sum, wrong matrix) proves the
mechanism can.

### REFUTED — the ordering hypothesis

`iter_ordered_valences` reproduces `_ordered_valences`' order element for element in O(1)
memory, so the hypothesis could be tested at the shipped 20,000 budget rather than argued.
`_ordered_valences` sorts by `(order_idx, tuple)` where `order_idx` is the position in
`product(O_sums, N_sums, C_sums, P_sums, S_sums)` — **O varies slowest, S fastest**. In these
ligands the atoms *with a choice* are overwhelmingly O and N, so the heuristic defers changing
exactly the atoms that must change.

| molecule | ligand atoms | charge | space | `Q0 == charge` | **C1 (what the code enumerates)** | first feasible rank, raw → heuristic |
|---|---|---|---|---|---|---|
| `QIDKUL_comp_0` | 37 | -8 | 1 259 712 | **16** | **16** | 209 858 → 209 858 *(same order)* |
| `QIDKIZ_comp_0` | 37 | -8 | 1 259 712 | **16** | **16** | 209 858 → 209 858 *(same order)* |
| `LIYFAA_comp_0` | 92 | -10 | 1 679 616 | **0** | **0** | — |
| `HICLAG_comp_0` | 147 | -2 | 53 747 712 | 97 868 | 105 620 | **4 → 32 768** |
| `BENVOG_comp_0` | 148 | -2 | 26 873 856 | 74 108 | 80 548 | 14 → 224 |
| `ZAZREZ_comp_0` | 144 | -2 | 8 503 056 | 1 424 304 | 1 734 696 | 648 → 648 *(same order)* |
| `KESWUB_comp_0` | 188 | -2 | 6 718 464 | 41 020 | 45 358 | 14 → 896 |
| `HOHKUL_comp_0` | 220 | -2 | 26 873 856 | 74 108 | 80 548 | 14 → 224 |

**Where the two orders differ, the heuristic order is WORSE (5 of 8), and for 3 of 8 it is not
a different order at all** — when every multi-choice atom belongs to one element group and
that group is traversed in atom order, the heuristic sequence *is* the raw sequence. For
`QIDKUL`, `QIDKIZ` and `ZAZREZ` the hypothesis is not merely false, it is **vacuous**.

On the real code path at 20,000 tries:

| molecule | raw | heuristic order | charge filter |
|---|---|---|---|
| `QIDKUL` ligand | 20 000 cands, no valid, 46.0 s | 20 000 cands, no valid, 41.2 s | **VALID at candidate 10, 0.02 s** |
| `LIYFAA` ligand | timeout (19 312 cands / 90 s) | timeout (16 354 cands / 90 s) | 0 feasible — provably hopeless |
| `HICLAG` ligand | timeout (5 942 cands / 90 s) | timeout (6 153 cands / 90 s) | **VALID at candidate 1 129, 18.6 s** |
| `BENVOG` ligand | VALID at 15, 0.28 s | VALID at 225, 4.56 s | **VALID at 1, 0.01 s** |
| `ZAZREZ` ligand | VALID at 649, 8.10 s | VALID at 649, 8.80 s | **VALID at 1, 0.02 s** |
| `KESWUB` ligand | timeout (6 033 cands / 90 s) | timeout (6 188 cands / 90 s) | **VALID at candidate 1 216, 22.1 s** |
| `HOHKUL` ligand | timeout (3 356 cands / 90 s) | timeout (3 947 cands / 90 s) | **80 000 of 80 548 feasible examined, none valid**, 2 572 s |

`BENVOG` and `HOHKUL` share every number in the DP table, which is exactly what a cache-keying
bug looks like (the sibling lane shipped one), so it was checked: they are distinct ligands
(148 atoms `C76H52N4O16`, `AC.sum` 320, sha `fea3a4d18cac`; 220 atoms `C100H100N4O16`,
`AC.sum` 464, sha `4a37c79878f2`) with the **same choice-bearing census** — 4 N and 16 O with a
choice, everything else pinned. The DP depends only on that census and the target charge, so
identical profiles are expected, not a collision.

### CONFIRMED — the real mechanism is FEASIBILITY DENSITY, and it is countable

`AC2BO`'s acceptance predicate has a necessary condition that is **additive over atoms**,
which makes the whole space decidable without searching it. A valid return forces
`BO.sum(axis=1) == valences` exactly (either `UA` is empty and `BO = AC`, or
`valences_not_too_large` bounds every atom above by `valences` while `(BO - AC).sum() ==
sum(DU)` fixes the total, so every per-atom slack is zero). Therefore `charge_is_OK` evaluates
`get_atomic_charge` on the *candidate's own* valences — a pure per-atom function. Write
`Q0(valences) = Σ_i get_atomic_charge(z_i, v_i)`; its only correction is `Q += 2` per
trivalent single-bonded carbon, and only while running below the target. Two necessary
conditions follow:

* **C1** `Q0 ≤ charge` and `charge − Q0` even;
* **C2** `sum(valences) − sum(AC_valence)` even, since every added bond raises `BO.sum()` by 2.

Both additive, so a suffix DP over `(Q0, parity)` yields, for the **whole** space and without
running the search: how many candidates can possibly be valid, and the exact **rank** of the
first one in either lexicographic enumeration order (mixed radix, so exact far beyond 2⁶⁴).
`tools/valorder_feasibility.py` does this for all eight known over-cap molecules in **8.4 s of
CPU**, and the instrument was validated three ways rather than trusted:

* **brute force** — enumerating QIDKUL's full 1 259 712-candidate space by hand agrees exactly:
  `survivors=16 first_rank=209858`, in both orders;
* **the real code path** — every predicted rank was hit by production `AC2BO` on the nose:
  BENVOG predicted at rank 14, found as candidate **15** (1-based); ZAZREZ 648 → **649**;
  BENVOG under heuristic order 224 → **225**;
* **the sibling lane** — the probe reproduces `best_BO` sha `9d234c9103a8` and both
  whole-encode OIN shas (`0d428d9dfc56`, `796221298c6c`).

**Two counts, and the wider one is operative.** `Q0 == charge` is the strict condition; **C1**
also admits `Q0 < charge` with matching parity. `iter_charge_feasible_valences` enumerates
**C1** — it must, or the filter would not be necessary-only — so every "candidate *n*" index is
an index into the C1 sequence, and C1 is the count a budget has to cover.

### CONFIRMED — the encoder was emitting UNPARSEABLE strings, and the filter fixes them

This is the accuracy result, and it is what makes the lane more than a perf note. Both
molecules are metal–NHC complexes (Rh(COD)Cl and AuCl of a
bis(2,4-dinitrophenyl)imidazol-2-ylidene):

| | default (raw prefix, 20 000) | `OIN_VALENCE_CHARGE_FILTER=1` |
|---|---|---|
| `QIDKUL_comp_0` whole encode | 124.44 s, `found_valid=0` on all 3 over-cap calls | **0.87 s**, `found_valid` on all of them |
| emitted OIN | `…-n2c{0}n(…)cc2…` | `…N2C{0}N(…)C=C2…` |
| does the string re-read? | **3/4 fragments — the ligand fails `KekulizeException`** | **4/4 fragments** |
| `QIDKIZ_comp_0` whole encode | 129.53 s, `found_valid=0` | **0.32 s** |
| does the string re-read? | **2/3 fragments — same `KekulizeException`** | **3/3 fragments** |

**124.44 s → 0.87 s is ~143x, and the string goes from one that does not parse to one that
does.** (The budget lever reaches a comparable wall time on the same molecule — 148.10 s →
0.92 s — but it does *not* fix the parse: it returns the same guess faster. The two levers buy
different things and should not be conflated.)

Why the default's string does not re-read, exactly: the guess is a **charge-separated**
perception and the emitted OIN drops the charges. Taking the emitted ligand fragment with slot
markers stripped —

| the fragment | result |
|---|---|
| as emitted, `…-n2cn(…)cc2…` | **`KekulizeException`, unkekulized atoms 8 22 23** |
| with the ylide charges the perception actually had, `…-[n+]2[c-]n(…)cc2…` | **OK** |
| as an imidazolium cation instead, `…-n2c[nH+](…)cc2…` | `KekulizeException` |

An aromatic 5-ring whose carbene carbon has three heavy neighbours, no hydrogen and no charge
has no electron to give the ring. Restoring the metal bond that `{0}` denotes cannot fix it
either — that makes the carbon 4-connected, still not kekulizable. The filter's perception
needs no ring charges at all, so the question does not arise.

The rescued structures were verified rather than assumed (`HICLAG` and `KESWUB` produce **no**
string on the default path, so there is no before/after diff to show):

| | `HICLAG` (147) | `KESWUB` (188) | `BENVOG` (148) |
|---|---|---|---|
| `found_valid` at candidate | 1 129 | 1 216 | 1 |
| `BO.sum()` / `best_BO` sha | 392 / `7ff80f862e9d` | 450 / `3deb2077e1c0` | 398 / `f4e60eba2807` |
| formal charges sum to target `-2` | **yes** | **yes** | **yes** |
| atoms above max tabulated valence | none | none | none |
| `SanitizeMol` | OK | OK | OK |
| nonzero formal charges | 2 × `O⁻` | 2 × `O⁻` | 2 × `O⁻` |

### CONFIRMED — one molecule is provably impossible

**At charge −10 with carbenes allowed, `LIYFAA_comp_0`'s 92-atom ligand has ZERO candidates
that can satisfy the predicate — out of 1 679 616.** Not "none in the first 20,000": none at
all, at any rank, in any order. No budget, no ordering and no matcher can perceive it there.
The default path nevertheless spends its full budget discovering this (19,312 candidates in
90 s before the probe's cap, in both order arms).

That **relocates** the problem: `LIYFAA` is not a search failure, it is a **charge-proposal**
failure. The extended-Hückel proposal handed `AC2BO` a target the ligand cannot reach, and the
charge/carbene ladder's job is to move off it. The DP makes "is this charge even reachable?" a
millisecond question — a cheap gate the ladder does not currently have.

### REFUTED — `HOHKUL` (near-certainly), and the defect is elsewhere

**80 000 of `HOHKUL`'s 80 548 possibly-valid candidates were examined (2 572 s) and not one is
accepted** — and `best_BO` never left `AC.copy()`, so no candidate even qualified for the
fallback update. The last **548 (0.7%)** are unexamined; closing them costs a fresh 43-minute
run for 0.7% of one molecule, so it is stated rather than spent. `HOHKUL`'s defect is
therefore **not** the enumeration: it is that `get_BO`'s greedy iterated matching saturates
none of its feasible assignments.

### REFUTED — "`max_weight_matching` is the wrong algorithm"

The premise is confirmed as a unit test, not merely noted: the graph `get_UA_pairs` builds
carries **no weight attributes**
(`test_valence_search_levers.py::TestMatcherLever::test_the_graph_get_ua_pairs_builds_carries_no_weights`),
so `nx.max_weight_matching` really is solving maximum *cardinality* matching. But the graphs
are **tiny** — 200 real graphs captured from a production run of the over-cap ligand are 14–17
nodes and 13–17 edges. **The cost was never per-call complexity, it was call count: 60 001
invocations on 15-node graphs.**

Interleaved microbenchmark, 3 repetitions over those 200 production graphs, arms alternating:

| matcher | wall (3 reps) | speedup | total cardinality | matching sha | resulting `best_BO` sha |
|---|---|---|---|---|---|
| `nx` — `nx.max_weight_matching` | 1.027 s | 1.00x | 1917 | `6a86cfa95a17` | `9d234c9103a8` |
| `maxcard` — same, `maxcardinality=True` | 5.469 s | **0.19x** | 1917 | `6a86cfa95a17` | `9d234c9103a8` |
| `greedy` — `nx.maximal_matching` | 0.059 s | **17.46x** | **1901** | `10d0add385b9` | **`e0a7c536aaa9`** |

The cost model checks out bottom-up, which is what licenses the "call count, not complexity"
claim: one `nx.max_weight_matching` on a 16-node graph is **482.6 µs**;
`docs/agentic-notes/v0.4.5/ENCODER_PERF_v0.4.5.md` counted **60 001** calls; 60 001 × 482.6 µs = **29.0 s** against
that lane's independently measured `get_UA_pairs` inclusive time of **28.31 s**. Nothing is
left over for an asymptotic effect to explain.

**Q3 is moot once the budget question lands.** Matching was ~28 s of a 49 s encode only
because it ran 60,001 times; at 200 candidates per arm it runs ~600 times and costs ~0.3 s.
The budget cut and the matcher swap attack the *same* seconds, and the budget cut takes them
**without changing a single perceived bond order** while the matcher swap changes them. The
budget cut strictly dominates.

### CONFIRMED — sub-cap ligands are untouched, by construction and by measurement

Structural argument: `_fallback_tries()` and both enumeration levers are read **only inside
`if over_cap:`**, so sub-cap ligands — 99.8% of the corpus — cannot reach them and are
byte-identical *by construction*, not by sampling. Measured anyway:

| control | result |
|---|---|
| valsearch, 12 molecules, `tries` 20 000 vs 200 | OIN identical **12/12**, `best_BO` identical **12/12**, `REPEAT-OK` **12/12**, OIN changed **0**; every molecule examined only **3 to 8** candidates |
| valorder, 16 molecules, filter OFF vs ON | OIN identical **16/16**, candidates examined identical **16/16**, `REPEAT-OK` **16/16**, `over_cap_calls` **0** on every molecule and arm, OIN changed **0**; **4 to 6** candidates each |

The valorder set is 12 dataset molecules (a seed-3 sample confirmed sub-cap by
`tools/valsearch_scan.py`, 0 of 53 ligand fragments over cap) plus all **four goldens** —
CisPlatin, Ferrocene, fac-Ir(ppy)₃, PdCl₂-R-BINAP.

**The 13th molecule of the valsearch control is worth more than the 12 that passed.**
`XIRMER_comp_0` is *sub-cap* and did not finish a single encode in 35 minutes — and it forks a
child mid-encode, the signature of the v0.4.4 SL5 resonance wrapper. Since it cannot reach the
valence-search fallback at all, its cost is something else entirely. **A slow encode is not
evidence of an over-cap ligand**, and any future attempt to size "encoder slowness" by this
lane's mechanism would be measuring the wrong thing.

### Instrumentation cost, measured rather than waved through

`_maximum_matching` reads `os.environ` on the encoder's hottest loop, so it was benchmarked:
`os.environ.get` is **1.403 µs** against **482.6 µs** for one matching call — **0.291%**, i.e.
~84 ms added to a ~29 s loop. Per-candidate counter increments are two dict operations against
milliseconds of real work per candidate. Small, but stated, because "surely negligible" is how
default paths get slower.

---

## What was done

All source changes are in **`src/oinsmiles/utils/xyz2mol_local.py`** unless noted.

### Instrumentation (ungated, shipped)

| thing | where | why |
|---|---|---|
| `possible_valences(AC_valence, atoms, allow_carbenes=True)` | `:884` | lifted **verbatim** out of `AC2BO`, which now calls it, so the two cannot drift. Makes the over-cap population answerable in ~0.04 s/molecule instead of ~50 s |
| `valence_combo_size(vll, cap=_VALENCE_COMBO_CAP)` | `:920` | mirrors `AC2BO`'s own early break, so the returned number is exactly what the encoder compares against the cap: the true product at or below the cap, "greater than the cap" above it |
| `AC2BO_STATS` + `reset_ac2bo_stats()` | `:861`, `:878` | 13 deterministic counters. Wall clock is unusable on this host (release sweep holds load above 12), so **every claim in the lane docs is a counter or a sha**, never a second |
| `_HEURISTIC_ELEMENTS = (8, 7, 6, 15, 16)` | `:939` | the O/N/C/P/S grouping constant, named so `iter_ordered_valences` and `_ordered_valences` cannot drift apart silently |
| `iter_ordered_valences(vll, atoms)` | `:949` | `_ordered_valences`' order, lazily, in O(1) memory — six nested products over the element groups |
| `iter_charge_feasible_valences(vll, atoms, charge, ac_valence)` | `:1036` | the C1/C2 suffix DP as a generator over the raw product |
| `charge_filter_supported(atoms)` | used at `:1421` | all 30 transition metals have an `atomic_valence` entry (`[20]`) and **no** `atomic_valence_electrons` entry, so `get_atomic_charge` is a `KeyError` for them. The filter declines such fragments and takes the historical path rather than crashing *earlier* than the default would |

`_ordered_valences` itself (`:1098`) is a named function precisely so the sorted walk can be
tested, and the property that matters is that **the extracted sorter is a PERMUTATION of the
raw `itertools.product` — same members, different order**. That is pinned by
`tests/unit/test_encoder_robustness.py::TestAc2boCapIsByteIdentical::test_ordered_valences_matches_unsorted_content`.

### Constants

* **`_VALENCE_COMBO_CAP = 500_000`** (`:818`) — above this, `AC2BO` skips the materialising
  sort and takes the bounded lazy product. **It must stay ≥ 100 000.** A
  cisplatin/ferrocene-scale ligand's valence product is only a few dozen, so a materially
  lower cap would push **real** ligands onto the fallback branch and break byte-identity for
  molecules that are perceived correctly today. Pinned by
  `tests/unit/test_encoder_robustness.py::TestAc2boCapIsByteIdentical::test_cap_is_large`
  (`assertGreaterEqual(_VALENCE_COMBO_CAP, 100_000)`).
* **`_VALENCE_FALLBACK_TRIES = 20_000`** (`:826`) — how many candidates the over-cap fallback
  grinds before returning `best_BO`. Far smaller than the sort cap because each iteration runs
  the full BO/charge check.

### The three levers (all default OFF; `main`'s output byte-unchanged)

| lever | env name | what it changes | verdict |
|---|---|---|---|
| budget | `OIN_VALENCE_FALLBACK_TRIES` (`_FALLBACK_TRIES_ENV`) | how many candidates the bounded prefix contains | **recommended at 200**, not applied |
| order | `OIN_VALENCE_ORDERED_FALLBACK` (`_ORDERED_FALLBACK_ENV`) | prefix taken in `_ordered_valences`' order | measured **WORSE**; kept only so the refutation is reproducible. **Nothing should ship behind it** |
| filter | `OIN_VALENCE_CHARGE_FILTER` (`_CHARGE_FILTER_ENV`) | prefix contains only C1/C2-feasible candidates | **recommended ON**, not applied |
| (measurement only) | `OIN_VALENCE_MATCHER` (`_MATCHER_ENV`) | `nx` / `maxcard` / `greedy` | both alternatives **rejected**; kept only for reproducibility |

Dispatch is at `:1413-1448`: the filter wins if both enumeration levers are set (it subsumes
the question), `over_cap_filter_unsupported` records a declined fragment, and
`itertools.islice(candidate_source, _fallback_tries())` applies the budget to whichever source
was chosen.

`_fallback_tries()` (`:838`) hardens the parse deliberately: garbage or non-positive values log
a warning and **fall back to 20 000** rather than silently collapsing every over-cap perception
to one candidate. A typo in an env var must not change perception.

### Two design properties that make the filter strictly dominant

1. **It is a subsequence of the raw product, in the same relative order.** It skips candidates;
   it never reorders them. So *the first valid candidate found is the one an unbounded raw
   search would have found*. Confirmed empirically: on `BENVOG` and `ZAZREZ` all three arms
   return the **same `best_BO` sha** (`f4e60eba2807`, `fddae6872a0b`).
2. **When the feasible set is empty it falls back to the historical enumeration** (`:1429-1440`).
   `best_BO` is assembled from candidates the filter drops — its own charge test is on the
   *BO*, not on the candidate — so without this the lever could turn one guess into a
   *different* guess. With it, the lever's entire blast radius is "a guess becomes a real Lewis
   structure". `LIYFAA` is the real molecule this covers. The cost is that `LIYFAA` keeps
   grinding a budget the DP has already proved cannot succeed; short-circuiting that is a
   separate, provably safe optimisation, left as an explicit decision rather than folded in.

### Alternatives rejected, with the reason

* **`maxcardinality=True`** — 5.3x *slower* for a bit-identical answer. Asking networkx's
  Blossom implementation for the max-cardinality variant costs more; it does not save.
* **`nx.maximal_matching` (greedy)** — 17x faster but *maximal*, not *maximum*: 1901 edges
  against 1917, and a **different `best_BO`** (`e0a7c536aaa9`) — the same total bond order via a
  different Kekulé structure. That is exactly the "silently perceives a different molecule"
  trade this release repeatedly caught, and it buys seconds the budget question already
  removes.
* **Lowering `_VALENCE_COMBO_CAP`** — would move real ligands onto the fallback branch. See the
  guard test above.
* **Applying the ordering heuristic to a bounded prefix** — the lane's own leading hypothesis,
  implemented and then **refuted by measurement** (above).
* **Promoting either recommended lever inside the lane** — declined, following the Y2 axial
  precedent: a change to perceived bond orders is a product call, and a lane should not make it
  silently.

---

## Dead ends, refutations, and costs accepted

### `found_valid = 0` is a statement about FEASIBILITY DENSITY, not about the cap

This is the single most important correction in the lane, because the naive reading — "the cap
is too small, raise it" — is wrong in a way that would waste arbitrary amounts of compute.
`found_valid = 0` means *the searched window contained no acceptable candidate*, and the
window's size is only one of three reasons that can happen:

| reason | example | does a bigger budget help? |
|---|---|---|
| the feasible set is dense and reachable | `BENVOG` (rank 14), `ZAZREZ` (rank 648) | already succeeds; nothing to fix |
| the feasible set is sparse and starts past the window | `QIDKUL`/`QIDKIZ` (16 feasible of 1 259 712, first at rank **209 858**) | **no** — you cannot afford the prefix; only *filtering* reaches it |
| the feasible set is **empty** | `LIYFAA` — **0 of 1 679 616** | **never**, at any budget, in any order, with any matcher |

For `LIYFAA` no budget increase whatsoever would help, and that is not an inference — it is a
closed-form count validated by brute force on QIDKUL's smaller space. The distinction matters
because "raise the budget" and "the ligand is unperceivable at this charge" call for opposite
follow-ups: the first is a compute decision, the second relocates the molecule to the
charge-proposal ladder.

### ⚠ The inverted `BO.sum()` claim — the correction that must not be lost

`docs/agentic-notes/v0.4.5/VALENCE_SEARCH_v0.4.5.md`'s fine budget sweep recorded:

> **1–25 candidates: genuinely worse.** Sum 92, a different matrix. Cutting this far would
> silently degrade perception.

**That inference is INVERTED, and the sha proves it.** The sum-92 matrix flagged as "worse" is
`53205ad7836d` — which is exactly the matrix the charge filter returns with **`found_valid=1`**,
i.e. the **VALIDATED** structure. The sum-94 matrix (`9d234c9103a8`) that the full
20,000-candidate search returns is the **unvalidated guess**. Checked atom by atom on both arms
of QIDKUL's ligand:

| | raw, 20 000 tries | charge filter |
|---|---|---|
| `found_valid` / `over_cap_exhausted` | 0 / 1 | **1 / 0** |
| candidates examined | 20 000 | **10** |
| `BO.sum()` | 94 | 92 |
| formal charges sum to the target `-8` | yes | yes |
| atoms above their max tabulated valence | none | none |
| `SanitizeMol` | OK | OK |
| perceived form | aromatic imidazolium ylide, `[c-]` + `[n+]`, 10 formal charges | **neutral NHC carbene**, `[C]`, 8 formal charges |
| emitted OIN re-reads | **no** (`KekulizeException`) | **yes** |

Both forms are drawable in isolation. Only one satisfies the algorithm's own predicate, and
only one survives being written into an OIN and read back.

**`BO.sum()` is NOT a quality metric and must not be used as one.** It is the tie-break inside
`best_BO`'s update rule (`BO.sum() >= best_BO.sum()`) and nothing more. A future agent reading
the fine sweep table in isolation will draw the wrong conclusion; this paragraph is the fix.
Anyone tempted to *report* a bond-order sum as evidence of a better perception should report
`found_valid` and a re-parse check instead — which is what `tools/valorder_encode_ab.py` does.

A related nuance that looks like a charge bug and is not: `charge_is_OK` (and
`set_atomic_charges` with it) applies a `+2` correction per trivalent single-bonded carbon
while the running total is below the target. On the four validated structures the per-atom
charges hit the target **exactly**, so the correction never fires. On the 200-candidate
`best_BO` guesses captured for `HICLAG` and `KESWUB` they sum to `-4` against a target of `-2`
and only reach the target *through* the correction — which writes `+1` onto a carbon
`get_atomic_charge` calls `-1`. The net charge ends up right either way; the difference is that
the guess needs a heuristic patch on two carbons and the validated structure needs none.

### The molecules that benefit most are the ones whose fidelity cannot be verified

Every over-cap molecule, each arm its own subprocess under a 200 s cap:

| molecule | ligand atoms | `tries=200` | `tries=20 000` |
|---|---|---|---|
| `KESWUB_comp_0` | 188 | TIMEOUT | TIMEOUT |
| `BENVOG_comp_0` | 148 | TIMEOUT | TIMEOUT |
| `HICLAG_comp_0` | 147 | **145.74 s**, 401 cands, `found_valid=1`, OIN `02b9b4e59da5` | **TIMEOUT** |
| `HOHKUL_comp_0` | 220 | TIMEOUT | TIMEOUT |
| `LIYFAA_comp_0` | 92 | **34.43 s**, 801 cands, `found_valid=1`, OIN `c5c84be56a0c` | **TIMEOUT** |
| `ZAZREZ_comp_0` | 144 | TIMEOUT | TIMEOUT |

`HICLAG` and `LIYFAA` complete at 200 and do **not** complete at 20,000 — so they are the two
molecules the budget cut most clearly helps, and they are **exactly the two for which no
byte-identity check is possible**, because the 20,000 baseline never produces a string to
compare against. Their OIN shas are recorded above so a future run on a quiet host can settle
it.

So the byte-identity claim is precisely this and no wider: **0 OIN strings changed out of 14
molecules where both arms completed** — 12 sub-cap plus `QIDKUL_comp_0` and `QIDKIZ_comp_0`. Of
the remaining six over-cap molecules, two changed from *no string at all* to a string, and four
produce no comparison in either direction. **"0 changed" covers two of the over-cap class's
eight known members, not the class.**

### A ligand-size taxonomy that did not survive

`docs/agentic-notes/v0.4.5/VALENCE_SEARCH_v0.4.5.md` grouped the over-cap population by ligand atom count and
concluded the 148–220-atom class was bound by *per-candidate cost*. The valorder measurements
overturn the split: **`KESWUB` at 188 atoms is rescued and `BENVOG` at 148 was never failing on
this axis**, while the 37-atom pair is among the *hardest* by feasible density (16 candidates in
1.26 M). **Size predicts per-candidate cost; it does not predict whether the search can
succeed. Density does.** Also: per-candidate cost measured on the ligand alone is
**0.012–0.027 s**, not the ~0.36 s the valsearch doc reported for `HICLAG` — that figure was
whole-encode wall divided by candidates and carried the rest of the encode, so the inferred
"70x spread between ligand sizes" is really about **5x**.

### Costs accepted

* **Neither recommended lever was promoted.** The measured upside is concentrated in ~0.2% of
  the corpus; the risk is a silent perception change. Left as an explicit product call.
* **The saturation margin is 2x, measured on ONE ligand.** The fine sweep that located
  convergence at 100 candidates was affordable only on the 37-atom ligand. A second fine sweep
  on `LIYFAA`'s or `HICLAG`'s ligand should precede promoting `OIN_VALENCE_FALLBACK_TRIES=200`.
* **Stability of the fallback answer is empirical, not structural.** `BO.sum() >=
  best_BO.sum()` means ties overwrite, and the 50-candidate row proves the mechanism *can*
  return a different equal-sum answer.
* **No pre-change baseline suite run.** Both halves recorded the same honest gap: a full run
  cost 890–2256 s under load and a second one would have taken cores the release sweep needed.
  The evidence that nothing regressed is the arithmetic (605 + 13 = 618; 618 + 16 = 634)
  together with **identical skipped and expected-failure counts** — a new xfail or skip is how a
  regression hides in a green suite.
* **`_rescue_unusable_perception` sweeps up to eight further charges through `AC2mol`**, which on
  an over-cap ligand is 8 more full searches. Not triggered by anything measured here, and not
  investigated.

### Instrument defects found (this release kept producing them)

* **The scan under-reported its own denominator.** Results were keyed by `path.stem`, and
  `cat/` and `photo/` share refcodes, so 8 of a 2000-molecule sample silently collapsed into
  1992 — shrinking the denominator of the very rate the tool exists to report. Now keyed
  `<subdir>/<stem>`, with a warning if any key still collapses.
* **`docs/agentic-notes/v0.4.5/ENCODER_PERF_v0.4.5.md`'s "two molecules" are one ligand.** Both `QIDKUL_comp_0` and
  `QIDKIZ_comp_0` reduce to the **same** 37-atom over-cap fragment at combo_size 1 259 712. A
  profile agreeing across "two molecules" is therefore weaker evidence of generality than it
  reads as.
* **The re-parse check in `valorder_encode_ab.py` is naive about some OIN spellings** (`[cH]`
  with an explicit hydrogen, `C#O` carbonyls) and flags a few sub-cap molecules in **both** arms
  identically. It is evidence of a *difference between arms on one molecule*, not an absolute
  grammar check.

---

## Where it landed

**Both halves are merged.** `git log --oneline main..swimlane/v045-valsearch` and
`main..swimlane/v045-valorder` are both empty — the branches are ancestors of `main`, merged
into `release/v0.4.5` as `d075f0d6` (valsearch) and `e4661843` (valorder), released by
`1450b5ce` (`release(v0.4.5): integrate 16 lanes and PROMOTE the six canonicality levers`).

Lane commits, oldest first:

| commit | subject |
|---|---|
| `fb9aa226` | valsearch: instrument the valence search, and two default-OFF levers |
| `21077309` | docs(valsearch): the search saturates at 0.5% of its budget |
| `1a11b620` | docs(valsearch): three classes, not one — and the caveat that matters most |
| `1ad885b0` | docs(valsearch): the big-ligand 20 000 baseline — 2 of 6 rescued, 4 untouched |
| `a551d2f1` | docs(valsearch): suite 618 OK, the sub-cap control, and the verdict |
| `1497cc90` | docs(valsearch): fix a stale count in the sub-cap paragraph |
| `b9ab6821` | docs(v0.4.5): valence-search floor is reducible ~100x — and my framing was half wrong |
| `50355757` | valorder: two over-cap enumeration levers, and the instrument that decides between them |
| `0790946d` | valorder: the order hypothesis is refuted — and the predicate is countable, not samplable |
| `ba591882` | valorder: the filter must decline a fragment it cannot reason about |
| `2b3a84b3` | valorder: read the filter lever once, and verify the rescued structures atom by atom |
| `7fd10555` | docs(valorder): HOHKUL settled, the suite, and the BENVOG/HOHKUL cache-collision check |
| `6f807c6f` | docs(v0.4.5): valence-order refutes my hypothesis — and finds unparseable OIN strings |
| `e145da7a` | docs(valorder): the final tree runs 634 OK, and the number that matters is the skip/xfail pair |
| `7f54454e` | docs(valorder): HOHKUL is 99.3% refuted — and my own feasible-count column was the wrong count |

(valsearch based on `main` @ `ebd6aabc`; valorder based on `swimlane/v045-valsearch` @
`1497cc90`.)

**Final lever state — all four still default OFF as of the v0.4.5/v0.4.6 releases.** None is
registered in `src/oinsmiles/oin/levers.py::_DEFAULT_ON` (which holds `OIN_BORON_CAGE`,
`OIN_CANONICAL_BODY`, `OIN_CANONICAL_PERCEPTION`, `OIN_CANONICAL_SLOTS`,
`OIN_CANONICAL_ETA_WINDING`, `OIN_STABLE_METAL_AC`, `OIN_STABLE_STEREO`). Note that the
valence levers are read through a **local** `_lever_enabled` in `xyz2mol_local.py` with
semantics identical to the registry: `0/""/false/no/off` disable. That forward compatibility
was checked rather than hoped for — dropping the registry file in and importing from both the
package root and the submodule resolves to `oinsmiles.oin.levers` with **no circular import**
(`oinsmiles/oin/` is a namespace package), and unset still reads OFF because neither lever is
in `_DEFAULT_ON`.

**Guard tests.**

`tests/unit/test_valence_search_levers.py` (13 tests) — the point is not that the levers work,
it is that the **default path is unchanged**:

* `TestFallbackTriesLever::test_unset_env_is_the_historical_constant`
* `TestFallbackTriesLever::test_empty_string_is_treated_as_unset`
* `TestFallbackTriesLever::test_env_overrides`
* `TestFallbackTriesLever::test_garbage_and_nonpositive_fall_back_to_the_default`
* `TestMatcherLever::test_unset_env_is_nx_max_weight_matching`
* `TestMatcherLever::test_unknown_matcher_falls_back_to_nx`
* `TestMatcherLever::test_named_arms_are_reachable_and_return_matchings`
* `TestMatcherLever::test_the_graph_get_ua_pairs_builds_carries_no_weights`
* `TestPossibleValencesExtraction::test_matches_an_inline_replica_of_the_original_loop`
* `TestPossibleValencesExtraction::test_allow_carbenes_false_drops_the_divalent_carbon_option`
* `TestPossibleValencesExtraction::test_combo_size_short_circuits_at_the_cap_like_ac2bo_does`
* `TestStatsCounters::test_reset_zeroes_every_counter`
* `TestStatsCounters::test_a_sub_cap_perception_counts_candidates_and_finds_a_valid_structure`

`tests/unit/test_valence_order.py` (16 tests), with the claim each pins:

| claim | test |
|---|---|
| `iter_ordered_valences` == `_ordered_valences`, element for element | `TestLazyOrderedValences::test_lazy_order_is_identical_to_the_sorted_path` (400 random configs, ≥150 actually compared, with an explicit anti-vacuity assertion) |
| every element group plus a non-ascending atom is covered | `TestLazyOrderedValences::test_covers_the_element_groups_the_heuristic_names` |
| no exponential materialisation before the first yield | `TestLazyOrderedValences::test_first_candidate_comes_from_an_unmaterialisable_space` (200 carbons, `[4,2,3]` each) |
| the filter never drops a candidate the **real** predicate accepts | `TestChargeFilterIsNecessary::test_every_candidate_the_real_predicate_accepts_survives_the_filter` — brute force over 11 whole small spaces running real `BO_is_OK`/`get_BO`, with an anti-vacuity guard |
| output is a subsequence, never a reordering | `TestChargeFilterIsNecessary::test_output_is_a_subsequence_of_the_raw_product_in_the_same_order` |
| an unreachable charge yields nothing | `TestChargeFilterIsNecessary::test_an_unreachable_charge_yields_nothing_at_all` |
| the filter actually prunes (a keep-everything filter would pass everything above) | `TestChargeFilterIsNecessary::test_it_actually_prunes` |
| an infeasible charge reproduces the default answer exactly | `TestInfeasibleChargeFallsBackToTheHistoricalPath::test_an_infeasible_charge_reproduces_the_default_answer_exactly` |
| the metal gap the decline exists for is real | `TestChargeFilterDeclinesWhatItCannotReasonAbout::test_the_gap_this_guard_exists_for_is_real` |
| a metal-bearing fragment takes the historical path | `TestChargeFilterDeclinesWhatItCannotReasonAbout::test_a_metal_bearing_fragment_takes_the_historical_path` |
| the test case really is over-cap | `TestOverCapLeverDefaults::test_the_branch_under_test_is_actually_the_over_cap_one` |
| unset takes neither lever | `TestOverCapLeverDefaults::test_unset_environment_takes_neither_lever` |
| **`"0"` disables rather than enabling** | `TestOverCapLeverDefaults::test_zero_disables_rather_than_enabling` |
| each lever is reachable when enabled | `TestOverCapLeverDefaults::test_each_lever_is_reachable_when_enabled` |
| the filter wins when both are set | `TestOverCapLeverDefaults::test_the_filter_wins_when_both_are_set` |
| sub-cap perception cannot change under any lever combination | `TestSubCapIsUntouched::test_neither_lever_can_change_a_sub_cap_perception` |

**Suite state at each tip:** valsearch — `discover tests/unit` **618 OK** (3 skipped, 3
expected failures) in 2256 s at load ~35 = the documented 605 baseline + 13. valorder — **632
OK** at `0790946d` (890 s) and **634 OK** at `7fd10555` (975 s) = 618 + 14 and 618 + 16, with
**identical** skip/xfail counts in both. `tests/unit/test_regression_stability.py` is inside
every run and green on the default path. `uvx ruff@0.15.20 check` and `format --check` clean
across `src`, `tools`, `tests`.

**Tools shipped:** `tools/valsearch_scan.py`, `tools/valsearch_budget_ab.py`,
`tools/valorder_probe.py`, `tools/valorder_feasibility.py`, `tools/valorder_encode_ab.py`.

---

## Open questions / for the next agent

1. **Promote `OIN_VALENCE_CHARGE_FILTER`?** This is the highest-value open decision in the
   lane. For: sub-cap is byte-identical by construction (16/16 measured); where the default
   already succeeds the filter returns the **identical** `best_BO`; where it fails, the filter
   produces a structure satisfying the algorithm's own predicate instead of a guess, and on
   QIDKUL/QIDKIZ it converts an OIN that does **not** re-read into one that does; when nothing
   is feasible it falls back so it cannot turn one guess into a different guess. Against: it
   **changes emitted strings** for the over-cap population (2 measured, 4 more going from no
   string at all to a string), so any hardcoded expectation for those refcodes moves; the
   counterfactual for `HICLAG`/`KESWUB`/`HOHKUL` cannot be diffed; and the `get_BO` matching
   defect is untouched, so promotion must not be read as "the over-cap class is solved".
2. **Promote `OIN_VALENCE_FALLBACK_TRIES=200`?** Run a second fine budget sweep on `LIYFAA`'s
   or `HICLAG`'s ligand first — the 2x margin over the measured convergence point of 100 rests
   on one ligand. Note this lever is *dominated* by the filter on the molecules where both
   apply, so consider whether it is needed at all once the filter is on.
3. **Short-circuit the provably-infeasible grind.** The DP proves in milliseconds that
   `LIYFAA`'s budget cannot succeed, yet the fallback still spends it. Skipping it is provably
   safe *as a search decision* but changes the returned `best_BO`, so it needs its own explicit
   call.
4. **Add the reachability gate to the charge ladder.** `LIYFAA` relocates from "search
   failure" to "charge proposal failure". `iter_charge_feasible_valences`' DP makes "is this
   charge even reachable for this ligand?" a millisecond question — a cheap pre-check
   `_select_lig_mol`'s ladder does not have, and one that could stop the ladder wasting a full
   arm on an impossible target.
5. **`HOHKUL`'s real defect is `get_BO`.** 80 000 / 80 548 feasible candidates examined, none
   accepted, `best_BO` never leaving `AC.copy()`. That points at the greedy iterated matching
   inside `get_BO`, not at the enumeration. Unowned.
6. **Enumerate the over-cap population exactly.** The ~0.2% rate is estimated from a
   1992-molecule sample, not enumerated. A full static scan is cheap in CPU (~33 s per 2000
   molecules unloaded) but projected to ~3.6 h at load 40 and was killed to avoid starving the
   release sweep. Worth running on a quiet host to turn the estimate into an exact list.
7. **`KESWUB`/`BENVOG`/`HOHKUL` fork a child during the encode**, consistent with the v0.4.4
   SL5 resonance wrapper and their `resonance_timeout` classification. Their dominant cost may
   not be `AC2BO` at all. Untested — and `XIRMER_comp_0` (sub-cap, no encode in 35 minutes,
   also forks) says the same thing from the other direction.
8. **Do not re-run the matcher experiment.** `OIN_VALENCE_MATCHER` exists only so the
   refutation is reproducible. Both alternatives are rejected on measured evidence: `maxcard`
   is 5.3x slower for a bit-identical answer, `greedy` returns a different Kekulé structure.
9. **Re-measure the whole lane once `swimlane/v045-encspeed`'s memo is in the baseline.** The
   speedups here overlap it; the composed figure is unknown and is *not* the product of the two
   ratios.
