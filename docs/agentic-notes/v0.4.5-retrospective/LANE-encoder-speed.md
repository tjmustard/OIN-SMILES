# Lane — encoder speed

The first time anyone profiled the encoder: `AC2BO` is 99.8% of a slow encode, half of that half is
recomputation, and removing it is 1.62-2.16x with a byte-identical string.

## ELI5

Every performance effort in this project's history targeted the half that *builds* 3D coordinates
from a string. Nobody had ever measured the half that *reads* a structure and writes the string —
and it turned out a hard molecule spends 50-60 seconds there, essentially all of it inside one
function that tries to work out where the double bonds go. That function gets called five times on
the same molecule with slightly different assumptions about charge, and most of the work it does is
identical every time, because the expensive part does not depend on charge at all. Remembering that
part instead of redoing it cut the time roughly in half, and the emitted string is byte-for-byte
identical — verified on 24 molecules, in both directions, with the memory cleared between molecules
so a lucky cache hit could not be what made the two versions agree.

## The work, visually

```
  ENCODE, XYZ -> OIN, ON A SLOW MOLECULE   (QIDKUL_comp_0: eta, 59 atoms, 49.43 s)
  =================================================================================

   .xyz file
      |
      v
   XYZToSMILES().convert()
      |
      +-- xyz2AC / GetAdjacencyMatrix ......................... cheap
      |     (xyz2AC_obabel: 0 calls -- AC comes from GetAdjacencyMatrix)
      |
      +-- xyz2mol.py::_select_lig_mol  THE CHARGE/CARBENE LADDER
      |      up to 5x  AC2mol -> AC2BO  on the SAME adjacency matrix,
      |      varying only `charge` and `allow_carbenes`
      |      (_rescue_unusable_perception can add EIGHT more arms on an AC of its own)
      |
      |      # AC sha      charge  allow_carbenes   wall      candidates
      |      1 df3f6198      -1        True         0.00 s        1
      |      2 6fb274fd      -8        True        16.81 s   20 000   <-+
      |      3 6fb274fd      -8       FALSE        16.99 s   20 000     | same AC
      |      4 6fb274fd      -6        True        16.73 s   20 000   <-+ 2 and 4 differ
      |      5 a914ed00       0        True         0.00 s        1        ONLY in charge
      |
      |      per-atom valence product > _VALENCE_COMBO_CAP (500 000)
      |        => islice(product, _VALENCE_FALLBACK_TRIES = 20 000)
      |        => grinds exactly 20 000 candidates, finds NONE valid, returns best_BO
      |           i.e. perception FAILS and the 20 000-candidate search buys a GUESS
      |
      |      inside each arm:
      |        get_UA -> get_UA_pairs -> get_bonds -> nx.max_weight_matching   CHARGE-FREE
      |                                                                       (recomputes
      |                                                                        bit-identically)
      |        get_BO  --------------------------------------------------- reads valences only
      |        BO_is_OK / charge_is_OK ---------------------------------- THE ONLY charge readers
      |
      +-- get_oin_string()  (eta winding, ring-rotation canon, inline emission) ... 0.03 s
      +-- aligner _brute_force_symmetries 1 call .................................. 0.01 s
      +-- ResonanceMolSupplier ................................................... 0 CALLS

   ATTRIBUTION:  AC2BO inclusive = 49.31 s = 99.8 %   <-- not "dominant". essentially all of it.


  MEASURED REDUNDANCY (deterministic; load-independent)
  =====================================================
   key                              total calls   distinct   removable
   get_UA_pairs (running matching)      60 001      26 668     55.6 %
   get_bonds by (AC, UA)               180 005       4 007     97.8 %
   get_BO by (AC, valences)             60 002      26 669     55.6 %  <- memo DECLINED


  WHAT SHIPPED, AND WHAT IT BOUGHT
  ================================
   [1] bounded LRU memo, get_UA_pairs + get_bonds       max_weight_matching 60 001 -> 26 668
       coarse key (AC, tuple(UA), tuple(du > 1))        get_UA_pairs wall  28.31 -> 14.86 s
       6 AC slots / 200 000 entries; `slot["ac"] is AC` get_bonds wall      4.10 ->  0.46 s
       returns a FRESH MUTABLE (callers mutate it)
   [2] charge_is_OK, same logic minus interpreter tax   charge_is_OK wall  14.35 ->  3.59 s
       (BO == 1).sum(axis=1) instead of per-carbon      get_atomic_charge
       list().count(1); get_atomic_charge INLINED;        4 440 174 -> 132 calls (-99.997 %)
       .tolist() to unbox np.int64
   [3] get_UA / valences_not_too_large: .tolist(), one  BO_is_OK wall       8.45 ->  3.17 s
       subtraction, append bound outside the loop
   [4] AC2BO's elif evaluates charge_is_OK IN PLACE     peak RSS 157 MB
       (it sat behind two cheaper short-circuits)

   INTERLEAVED PAIRED A/B (BASE/NEW/BASE/NEW, load 22-38 -- the ratio is the claim)
     QIDKUL_comp_0 (eta, 59 at)   104.94 / 107.10 s  ->  66.02 / 64.66 s   = 1.62x
     QIDKIZ_comp_0 (non-eta, 39)  123.98 / 128.13 s  ->  54.09 / 62.38 s   = 2.16x
     all 8 runs byte-identical:  QIDKUL 0d428d9dfc56   QIDKIZ 796221298c6c


  THE FLOOR THAT REMAINS  (and why it cannot be memoized)
  =======================================================
       CHARGE-FREE HALF                     CHARGE-DEPENDENT HALF
       get_UA / get_UA_pairs / get_bonds    BO_is_OK + charge_is_OK + get_BO
       identical across ladder arms         DIFFERENT per arm, by construction
       ==> MEMOIZED (this lane)             ==> ~13-15 s, IRREDUCIBLE by memoization
                                               = the floor for over-cap ligands

  LEGEND
    ->      becomes / control flow            <-  annotation on the line above
    [n]     the four shipped changes          ==> conclusion
    CAPS    the load-bearing fact on that line
```

## Initial assumptions and hypothesis

1. **The lane existed because the generation-side perf lane found a wall it could not move.** That
   lane measured a single eta molecule's bare `XYZToSMILES().convert()` at **46-71 s**, and a round
   trip runs at least two encodes — so such a molecule's floor was **~90-140 s before the generator
   did anything**. No generation-side fix can reach a 30 s round trip through that.

2. **Nobody had ever profiled the encoder.** v0.4.0's P2-P11, v0.4.4's `OIN_EARLY_EXIT`, and
   v0.4.5's perf lane all targeted generation.

3. **Hypotheses going in, all of them plausible and all of them wrong** (each is refuted with a
   number below): the aligner is factorial in coordination number, so it is probably the aligner;
   eta canonicalisation is elaborate, so it is probably that; `ResonanceMolSupplier` needed a
   forked, CPU-bounded wrapper in v0.4.4 SL5, so it is probably that; the molecule is 59 atoms, so
   it is probably size; the slow molecule is eta, so it is probably eta.

4. **A standing methodology assumption:** **wall-clock timing is meaningless above roughly load 12
   on this machine** (`tools/v045_state.sh`). Two sibling lanes were running throughout this lane's
   work at load 2-38, so the lane's design commitment was that every headline number must be either
   a deterministic counter or a ratio from an interleaved paired A/B.

5. **A hard constraint, not an assumption: the encoder's output *is* the product.** Any speed change
   must be SHA-identical on the emitted string. This is stricter than the generation side, where a
   different-but-equivalent conformer is acceptable.

## What was actually found

### Confirmed

* **`xyz2mol_local.AC2BO` is 99.8% of a slow encode.** Not "dominant" — essentially all of it.

  | | `QIDKUL_comp_0` (eta, 59 atoms) | `QIDKIZ_comp_0` (non-eta, 39 atoms) |
  |---|---|---|
  | encode wall | 49.43 s | 57.58 s |
  | `AC2BO` inclusive | **49.31 s (99.8%)** | **57.47 s (99.8%)** |
  | `get_oin_string` (the entire OIN emission) | 0.03 s | 0.02 s |
  | aligner (`_brute_force_symmetries`) | 0.01 s | 0.01 s |
  | `ResonanceMolSupplier` | **0 calls** | **0 calls** |

* **It is not an eta phenomenon and not a size phenomenon.** The **non-eta 39-atom** molecule has
  the *same* profile and is *slower*. It is a property of the ligand's valence combinatorics.

* **The mechanism.** `xyz2mol.py::_select_lig_mol` is a charge/carbene ladder: it calls
  `AC2mol` -> `AC2BO` **up to five times on the same adjacency matrix**, varying only `charge` and
  `allow_carbenes`. For these ligands the per-atom valence product exceeds `_VALENCE_COMBO_CAP`
  (**500 000**), so each call takes the `itertools.islice(product, _VALENCE_FALLBACK_TRIES)` branch
  and grinds **exactly 20 000 candidate valence assignments, finds none valid**, and returns
  `best_BO`. (`_rescue_unusable_perception` can add **eight more arms** on an AC of its own.)

* **Calls 2 and 4 of the ladder differ *only* in charge**, and only `BO_is_OK` and `charge_is_OK`
  read `charge` — the candidate-generation half (`get_UA` -> `get_UA_pairs` -> `get_BO`) is
  charge-independent and recomputes bit-identical results. Measured redundancy: `get_UA_pairs`
  60 001 calls / 26 668 distinct (**55.6% removable**), `get_bonds` by `(AC, UA)` 180 005 / 4 007
  (**97.8%**), `get_BO` by `(AC, valences)` 60 002 / 26 669 (55.6%).

* **Deterministic counters after the fix** (`QIDKUL_comp_0`; exact, contention-robust):

  | counter | before | after | change |
  |---|---|---|---|
  | `nx.max_weight_matching` calls | 60 001 | **26 668** | -55.6% |
  | `get_atomic_charge` calls | 4 440 174 | **132** | -99.997% |
  | `get_bonds` calls | 180 005 | 146 672 | -18.5% |
  | `get_bonds` inclusive wall | 4.10 s | **0.46 s** | -89% |
  | `charge_is_OK` inclusive wall | 14.35 s | **3.59 s** | -75% |
  | `BO_is_OK` inclusive wall | 8.45 s | **3.17 s** | -62% |
  | `get_UA_pairs` inclusive wall | 28.31 s | **14.86 s** | -48% |
  | peak RSS | — | 157 MB | the memo is not a memory problem |

* **Wall clock, with conditions stated.** The host was **not quiet** — two sibling lanes ran
  throughout, load 2-22 — so absolute seconds are not clean, and the lane reported two measurements
  that are honest about it. Single runs: **49.43 s at load 2.4 (before)** vs **24.82 s at load 12.2
  (after)** — faster despite 5x the load, so the improvement is real and understated. Then the
  load-robust measurement, an **interleaved paired A/B** alternating BASE/NEW/BASE/NEW so contention
  averages across both arms rather than favouring whichever ran first, at load 22-38:

  | molecule | BASE | NEW | speedup |
  |---|---|---|---|
  | `QIDKUL_comp_0` (eta, 59 atoms) | 104.94 s, 107.10 s | 66.02 s, 64.66 s | **1.62x** |
  | `QIDKIZ_comp_0` (non-eta, 39 atoms) | 123.98 s, 128.13 s | 54.09 s, 62.38 s | **2.16x** |

  All eight runs re-confirmed byte-identity (`QIDKUL` sha `0d428d9dfc56`, `QIDKIZ` sha
  `796221298c6c` on every arm and repetition) — a third independent confirmation after the
  24-molecule manifest and the 3-convert warm-memo check.

* **Cross-encode reuse works, which is what matters for a round trip.** Three consecutive
  `convert()` calls on `QIDKUL_comp_0` in one process, load 17:

  | | wall | `max_weight_matching` calls | OIN sha |
  |---|---|---|---|
  | 1st | 34.30 s | 26 668 | `0d428d9dfc56e18e` |
  | 2nd | **14.79 s** | **0** | `0d428d9dfc56e18e` |
  | 3rd | 15.73 s | 0 | `0d428d9dfc56e18e` |

  Peak RSS 157 MB across all three. A round trip re-perceives the same ligand connectivity many
  times — the input, then every generated conformer the SL1 accept-first check re-encodes.

* **Byte-identity, 24 molecules, two-directional.** `tools/enc_byte_identity_ab.py` with
  `git show HEAD:… > …` in both directions (**never `git stash`** — it is shared across worktrees),
  and the memo cleared between molecules so a cross-molecule hit cannot be what makes the arms
  agree. Set: **4 goldens** (CisPlatin, TransPlatin, Ferrocene, POJJOP) **+ 8 fixtures**
  (Cis-PtCl2(en), fac/mer-Ir(ppy)3, PdCl2-R-BINAP, Zeise's salt, YESKOZ, FeCO5, CuCN2) **+ 12
  dataset molecules**, **7 of them eta, 0 errors**. Manifest sha256 on **both** arms:
  `9b679638dfeb3cf81abaf563e2b59cced3e3a5a6f3e0005fb983c9f57c3b014b`. The only differing line
  between the two runs is the `--label` echo.

* **Suite:** `discover tests/unit` **605 tests OK** (3 skipped, 3 expected failures) — matching the
  documented baseline — plus 8 new guards in `tests/unit/test_ac2bo_memo.py`;
  `tests/unit/test_regression_stability.py` green; `ruff check` / `format` clean.

* **The 30 s question, answered honestly.** The **encode** now fits under 30 s for these molecules:
  `QIDKUL_comp_0` is ~25 s at load 12, ~15 s warm. **The round trip still does not**, because it
  runs at least two encodes — cold + warm ~= 25 + 15 = **40 s of encode alone** before the generator
  starts. This lane does not change that and says so.

* **The class for which 30 s is unreachable, with the evidence.** Ligands whose per-atom valence
  product exceeds `_VALENCE_COMBO_CAP`. What remains after this lane is the **charge-dependent**
  half — `BO_is_OK` + `charge_is_OK` + `get_BO` over 60 000 candidates, **~13-15 s** — and it is
  *not* memoizable, because it is exactly the part that differs per arm. That is the floor.

* **The encoder's own honest framing, recorded rather than buried:** for these ligands perception
  **fails** (no valid Lewis structure among 20 000 candidates) and `best_BO` is a fallback. *The
  20 000-candidate search is buying a guess.* Whether that guess is worth 15 s per ladder arm is a
  product question this lane surfaced rather than answered — and the follow-on `valsearch` lane then
  answered it.

### Refuted

Each of these was ruled out **with a number**, which is why they should not be re-chased:

* **`ResonanceMolSupplier` is not on the encode hot path** — **0 calls** in a full `convert()` of
  either slow molecule. The v0.4.4 SL5 forked, CPU-bounded wrapper never fires here.
* **The aligner is not the encoder's cost.** `_map_to_template` 3 calls / 0.00 s,
  `_brute_force_symmetries` 1 call / 0.01 s. "Factorial in coordination number" is real in principle
  and **0.02%** of this encode in practice.
* **Eta canonicalisation is not the encoder's cost.** All of `get_oin_string` — eta winding,
  ring-rotation canonicalisation, inline emission — is **0.03 s**, 0.06% of the encode. **Eta
  molecules are slow to *generate*, not to encode.**
* **Repeated `SanitizeMol` / `MolToSmiles` on an unchanged mol is not the cost** — 19 and 9 calls,
  0.00 s total.
* **`xyz2AC_obabel` is not the cost** — no calls; the AC comes from `GetAdjacencyMatrix`.
* **Memoising whole `AC2BO` results wins nothing** — 5 calls, 5 distinct argument tuples.
* **`get_BO` memo declined** despite a 55.6% hit rate: caching 26 669 N x N integer arrays is
  **~370 MB for ~2.3 s**. The cost/benefit is wrong.
* **Encoder slowness is not eta-specific and not size-driven** — the non-eta 39-atom molecule is
  slower than the eta 59-atom one, with an identical profile.

## What was done

Lane branch `swimlane/v045-encspeed`, based on `main` @ `b23decb4`. Git-durable record:
`docs/agentic-notes/v0.4.5/ENCODER_PERF_v0.4.5.md`. As with the perf lane,
`git log --oneline main..swimlane/v045-encspeed` returns **nothing** — the branch is now fully
merged into local `main` (`git merge-base --is-ancestor 3b15e26c main` is true today), so its tip is
the merge-base. The lane's commits are:

```
git log --oneline b23decb4..swimlane/v045-encspeed
3b15e26c  perf(encode): paired A/B ratio 1.62x eta / 2.16x non-eta; snapshot the LRU walk
b2d7433b  docs(encspeed): the encoder profiled for the first time -- AC2BO is 99.8%
0a6ee660  perf(encode): halve a slow encode -- AC2BO was 99.8% of it, half of that redundant
ca462fbb  tools(encspeed): encode-side attribution probes -- AC2BO is 99.8% of an eta encode
```

All code changes are in `src/oinsmiles/utils/xyz2mol_local.py`, and all are **ungated** —
dead-work removal and memoization are the two classes this codebase has previously proved
byte-identical, and the byte-identity A/B covers 24 molecules.

### 1. A bounded LRU memo for `get_UA_pairs` and `get_bonds`

The soundness argument is **the coarse key**: `get_UA_pairs`' result is a function of
`(AC, tuple(UA), tuple(du > 1 for du in DU))`, because `DU` is read **only** through the `du > 1`
predicate that decides how many virtual matching nodes get allocated. Two calls agreeing on that key
are handed the identical edge list in the identical insertion order, and `nx.max_weight_matching` is
deterministic on identical input — so the memo can only make the same answer arrive sooner.

Three implementation details that are load-bearing, not decorative:

* **An LRU over adjacency *matrices*** — `_AC2BO_MEMO_SLOTS = 6`, `_AC2BO_MEMO_MAX = 200_000` —
  rather than a single slot, because a round trip re-perceives the same ligand connectivity many
  times (the input, then every generated conformer the SL1 accept-first check re-encodes), each a
  fresh array with identical contents.
* **It is read only for an array the cache holds a live reference to** (`slot["ac"] is AC`).
  `id()` is just a fast index; **the identity test is what makes a recycled `id` unable to alias a
  stale entry.**
* **A hit returns a fresh mutable object.** Callers mutate what they get — `get_UA_pairs` appends
  the virtual-node edges to `get_bonds`' list — so a shared object would poison every later hit.
  `tests/unit/test_ac2bo_memo.py` pins this by poisoning a result and checking the next hit is clean.

### 2. `charge_is_OK` — same logic, minus the interpreter overhead

It was the second-largest cost (14.4-16.7 s). Three sources, all removed:

* `list(BO[i, :]).count(1)` allocated a fresh 59-element list of boxed `np.int64` **per carbon, per
  call**. One `(BO == 1).sum(axis=1)` gives every atom's count at once.
* `get_atomic_charge` was one Python call per atom — **4 440 174 calls in a single encode**. Inlined
  as the identical `elif` ladder. Now 132 calls, from `BO2mol` only.
* `list(BO.sum(axis=1))` yields boxed `np.int64`, so every downstream `1 - v`, `v == 2` and `Q += q`
  went through numpy's scalar machinery. `.tolist()` unboxes once.

**The carbon corrections stay a Python loop on purpose:** `if ns == 3 and Q + 1 < charge` reads the
*running* total, so they cannot be reordered or vectorised without changing the answer. Its dead
`q_list` accumulator (built, never read) is dropped.

### 3. `get_UA` / `valences_not_too_large`

`.tolist()` over `list()`, the subtraction done once instead of twice, `append` bound outside the
loop.

### 4. `AC2BO`'s `elif` evaluates `charge_is_OK` in place

It previously ran eagerly on the line above, sitting behind two cheaper short-circuiting predicates
and computed even when `status` had already returned. Pure dead-work removal; its only side effect is
a DEBUG log line reachable solely on the `allow_carbenes=False` arm, so the emitted OIN is untouched.

### The tools, and exactly what each measures

**`tools/perf_encode_profile.py`** — the attribution instrument, two modes on one molecule.
`--mode counters` monkeypatches ~25 named hot-spot candidates (`COUNTER_NAMES`) to count **and**
accumulate inclusive wall time: RDKit `SanitizeMol` / `MolToSmiles` / ring perception /
`GetSubstructMatches` / `AssignStereochemistry` / `rdCIPLabeler` / `CanonicalRankAtoms`,
`ResonanceMolSupplier`, all of `xyz2mol_local`'s `AC2BO` / `get_UA_pairs` / `get_BO` / `BO_is_OK` /
`charge_is_OK` / `get_UA` / `get_bonds` / `get_atomic_charge` / `xyz2AC_obabel`, `nx.max_weight_matching`,
the driver's `get_oin_string` / `get_tmc_mol`, and the aligner's `_map_to_template` /
`_brute_force_symmetries` / `canonical_eta_set_representative` / `_orientation_symmetry_graph`. Counts
are exact and contention-robust. `--mode profile` adds cProfile, explicitly labelled in its own output
as "relative attribution only; host may be contended". It prints `loadavg` at start **and** end, so
every artifact records the conditions it was taken under.

**`tools/perf_encode_keys.py`** — how coarse a memo key `AC2BO`'s inner loop tolerates, and how much
it collapses. It states four structural claims read off the source and then measures each: that
`get_UA_pairs` reads `DU` only through `du > 1`; that `get_bonds` is a function of `(AC, tuple(UA))`
alone; that `get_BO` is a function of `(AC, tuple(valences))`; and that the ladder calls `AC2BO`
repeatedly on the same AC so candidates recur across arms. It prints total vs distinct per key — the
**exact upper bound on what memoization can remove** — plus the per-`AC2BO`-call table (AC sha,
charge, `allow_carbenes`, wall, candidate count) that appears in the diagram above. Deterministic;
unaffected by host load.

**`tools/perf_encode_redundancy.py`** — the narrower companion: total vs distinct argument keys for
`AC2BO` and `get_UA_pairs`, plus the count of valence combos each `AC2BO` call actually iterates
(via a counted `get_UA`) and per-call wall. It replicates `get_UA_pairs`' cheap early-out test to
separate calls that genuinely run the matching from those that do not, without running it twice.

**`tools/enc_byte_identity_ab.py`** — the acceptance gate. Prints one
`name<TAB>sha256<TAB>len<TAB>eta` line per molecule, sorted, plus a sha256 over the whole manifest so
two revisions can be compared as a single string. `eta` flags a haptic OIN via the `\{\d+[<>]\}`
detector **so the set is visibly not just cisplatin**. It clears the `AC2BO` memo
(`loc._ac2bo_memo_clear()`) between molecules, and it treats an exception as part of the contract:
an error is recorded as `ERROR:<Type>:<msg>` and must be the **same** error on both arms.

## Dead ends, refutations, and instrument failures

### The measured NEGATIVE table — four "obviously right" optimisations that were wrong

Recorded because each looked correct and was measured wrong. Method: same-process A/B on the real
37-atom ligand AC, **both implementations interleaved in one process so the ratio is immune to load
drift**, and asserting the two agree on every input.

| change | speedup | verdict |
|---|---|---|
| `charge_is_OK` rewrite | **4.45x** | kept (agrees on all 320 cases) |
| `get_UA` + `.tolist()` input | **2.36x** | kept |
| `valences_not_too_large` `.tolist()` | 1.09x | kept (marginal) |
| numpy-vectorising `charge_is_OK` with `np.select` | **~1.0x** | **rejected** |
| numpy-vectorising `get_UA` (`flatnonzero`) | **0.98x** | **rejected** |
| `BO.sum() - AC.sum()` for `BO_is_OK`'s `check_sum` | **0.72x** | **rejected** |

**The lesson is uniform: at ~40-60 atoms a numpy round trip costs more than the Python loop it
replaces.** The first version of this lane's `charge_is_OK` *was* the `np.select` one, and a
whole-encode measurement showed **14.35 s -> 14.17 s — no gain**. Vectorising the small per-atom
loops is a reflex worth resisting here; **unboxing `np.int64` is the change that actually pays.**

Note the deliberate methodological asymmetry, and keep it: **same-process interleaving is the right
design for an *agreement* check and for a *ratio between two implementations of one function*,
because it cancels load drift and lets you assert equality on every input. It is the WRONG design
for a whole-program timing A/B**, where the second arm inherits warm caches — an error made
elsewhere in this release, which reported levers-ON as faster (42.8 s vs 54.4 s) for exactly that
reason.

### Two levers that are real but **not** byte-identical, so they were not touched

* **Lowering `_VALENCE_FALLBACK_TRIES`.** `best_BO` is whichever candidate satisfied
  `BO.sum() >= best_BO.sum() and valences_not_too_large and charge_OK` **among those tried**, so
  cutting the cap changes the perceived bond orders in principle. It cannot be tuned ungated.
* **Replacing `nx.max_weight_matching`.** The graph carries no weight attributes, so it is really
  maximum *cardinality* matching and a faster algorithm exists — but it can return a **different
  equally-large matching**, hence a different `BO`.

Both were named as real options **behind a lever, with a corpus A/B**; neither was attempted here.
The follow-on `valsearch` lane (`docs/agentic-notes/v0.4.5/VALENCE_SEARCH_v0.4.5.md`) then did exactly that, and its
result sharpens this lane's claim rather than contradicting it — see the discrepancies section.

### The `get_BO` memo: a 55.6% hit rate that was correctly declined

Declining an available win is worth recording as explicitly as taking one. 26 669 N x N integer
arrays is ~370 MB to save ~2.3 s. The right call, and the reason the shipped memo is bounded
(6 slots / 200 000 entries) rather than unbounded.

### The load rule, violated once and over-applied once

* **Violated:** a timing probe was run concurrently with a 6-shard sweep at **load 21-33**. Its
  pass/fail outcomes survived; **its `elapsed_s` values did not**, and the timing half of the result
  was lost. This lane's own wall-clock numbers are reported *with* their load precisely because of
  that.
* **Over-applied:** a **byte-identity A/B was deferred as "needs an idle machine"**. Wrong —
  **string equality is deterministic and load-independent.** Contention makes it slow, not wrong.

Carry both forward as one rule: **classify a measurement as load-sensitive or load-independent
before deciding when to run it.** Counters, string equality, marker presence, parse outcomes, and
distinct-key counts are load-independent. Seconds are not.

### Byte-identity probe hygiene

* **Never `git stash`** for an A/B — it is shared across worktrees and a sibling lane collides with
  it. Use `git show HEAD:<path> > <path>` with an `EXIT` trap restoring the working copy. This
  lane's docstring says so on line 7.
* **Python block-buffers stdout into a pipe.** A `timeout` kill therefore discards the buffered
  lines, and a downstream `sort` still exits 0 — producing an **empty file that looks like
  agreement**. Write directly to the output file with `buffering=1` and finish with a `#DONE <n>`
  sentinel, so "incomplete" is instantly distinguishable from "agreed". A later full-encoder
  byte-identity run was in fact killed at 6 of 61 fixtures and the **missing sentinel is what
  revealed it**.
* **Clear the memo between molecules.** A cross-molecule cache hit must never be what makes two
  revisions agree. `enc_byte_identity_ab.py` does this via `loc._ac2bo_memo_clear()`.

### The n=2 problem, and a correction the follow-on lane supplied

This lane deliberately measured **two molecules** and ran no corpus sweeps, so the *size* of the
affected class was left unknown. The `valsearch` lane later found something this lane could not
have known and that materially weakens generalising from its sample: **`QIDKUL_comp_0` and
`QIDKIZ_comp_0` share the same 37-atom ligand** (both at `combo_size` 1 259 712 > 500 000). The
"two molecules" are **one ligand measured twice**. The profile, the 5-call ladder, the redundancy
ratios and the floor are all properties of that one ligand.

That does not undermine the shipped change — byte-identity was verified on 24 molecules, and the
memo's soundness argument is structural, not statistical — but any claim of the form "encodes are
1.62-2.16x faster" is a claim about **one ligand**, and the corpus-wide distribution of over-cap
ligands is measured in the valsearch lane, not this one: **~0.2% of molecules** (4/1992 at seed 7;
0/100 at seed 1), rising to **1.3%** in the slow/failing cohort.

## Where it landed

Code: `src/oinsmiles/utils/xyz2mol_local.py` — the bounded `AC2BO` LRU
(`_AC2BO_MEMO_SLOTS = 6`, `_AC2BO_MEMO_MAX = 200_000`, `_ac2bo_memo_for` / `_ac2bo_memo_anchor` /
`_ac2bo_memo_entries` / `_ac2bo_memo_clear`), the rewritten `charge_is_OK`, the `get_UA` /
`valences_not_too_large` tidy-ups, and `AC2BO`'s in-place `charge_is_OK` evaluation. **Ungated.**
Merged into local `main`.

Tests: `tests/unit/test_ac2bo_memo.py` (8 guards, including the poisoned-result check that pins
"a hit returns a fresh mutable").

Doc: `docs/agentic-notes/v0.4.5/ENCODER_PERF_v0.4.5.md`. Follow-on: `docs/agentic-notes/v0.4.5/VALENCE_SEARCH_v0.4.5.md`.

Below, `DS=/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset` and the
interpreter is the main checkout's pinned `.venv/bin/python` (rdkit `==2025.9.3`; never `uv sync` in
a worktree).

### `tools/perf_encode_profile.py` — where does an encode's time go

```bash
PYTHONPATH=src .venv/bin/python tools/perf_encode_profile.py \
    --dataset "$DS" --molecule QIDKUL_comp_0 --mode both --top 40 \
    --json-out /path/encode_profile_QIDKUL.json
```

Trust the `=== COUNTERS ===` block for attribution; treat the cProfile block as relative only.

### `tools/perf_encode_keys.py` — the exact upper bound on memoization

```bash
PYTHONPATH=src .venv/bin/python tools/perf_encode_keys.py \
    --dataset "$DS" --molecule QIDKUL_comp_0 \
    --json-out /path/encode_keys_QIDKUL.json
```

Reports `total_calls` vs `distinct_by_coarse_key(AC,UA,du>1)` and the per-`AC2BO`-call ladder table.

### `tools/perf_encode_redundancy.py` — total vs distinct, plus candidates iterated

```bash
PYTHONPATH=src .venv/bin/python tools/perf_encode_redundancy.py \
    --dataset "$DS" --molecule QIDKIZ_comp_0 \
    --json-out /path/encode_redundancy_QIDKIZ.json
```

### `tools/enc_byte_identity_ab.py` — the acceptance gate for any encoder change

```bash
F=src/oinsmiles/utils/xyz2mol_local.py

# arm NEW (working tree)
PYTHONPATH=src .venv/bin/python tools/enc_byte_identity_ab.py \
    --dataset "$DS" --label NEW > /path/enc_new.txt 2>/path/enc_new.err

# arm BASE -- two-directional git show. NEVER git stash (shared across worktrees).
cp "$F" /path/mine.py
trap 'cp /path/mine.py "$F"' EXIT
git show HEAD:"$F" > "$F"
PYTHONPATH=src .venv/bin/python tools/enc_byte_identity_ab.py \
    --dataset "$DS" --label BASE > /path/enc_base.txt 2>/path/enc_base.err
cp /path/mine.py "$F"; trap - EXIT

# the only permitted difference is the --label echo on the MANIFEST_SHA256 line
diff /path/enc_base.txt /path/enc_new.txt
grep MANIFEST_SHA256 /path/enc_base.txt /path/enc_new.txt
# expected on both: 9b679638dfeb3cf81abaf563e2b59cced3e3a5a6f3e0005fb983c9f57c3b014b
```

### Discrepancies between the sources, flagged rather than smoothed

* **`docs/agentic-notes/v0.4.5/VALENCE_SEARCH_v0.4.5.md` states "`swimlane/v045-encspeed` is not merged into `main` —
  `git merge-base --is-ancestor 3b15e26c main` is false".** As of today that check is **true**: the
  branch **is** merged. The statement was accurate when written (that lane was based on `main` @
  `ebd6aabc`, pre-merge) and is now stale. The substantive warning it carries is still live and must
  be preserved: **valsearch's speedups were measured against a no-memo baseline where the same path
  is ~35-50 s per arm, so they OVERLAP this lane's win rather than composing with it. Do not add the
  two ratios.**
* **"Lowering `_VALENCE_FALLBACK_TRIES` is not byte-identical" needs a precision.** It is not
  byte-identical *in principle* — `best_BO` depends on which candidates were tried. **Measured**, on
  `QIDKUL_comp_0`, budgets of 20 000 / 5 000 / 1 000 / 200 all produced the **same `best_BO` sha
  (`1ea2ed6bafe7`) and the same OIN sha (`0d428d9dfc56`)**, at 148.10 s / 27.29 s / 3.50 s /
  **0.92 s**; valsearch reports **0 changed OIN strings across 14 molecules where both arms
  completed**. So the principled objection stands (which is why the knob is a lever,
  `OIN_VALENCE_FALLBACK_TRIES`, default unset = 20 000, **not** a new default), but "it changes the
  answer" is not what was observed.
* **This lane's "eta molecules are slow to *generate*, not to encode" is correct and is easily
  misread.** The complementary fact from the generation side is that eta is **23.3% of molecules but
  35.6% of wall clock**, with the penalty holding inside every size band (0-50 atoms 7.9 vs 3.3 s;
  50-80 15.5 vs 5.0 s; 80-120 49.7 vs 8.1 s) — and eta molecules are *smaller* on median (62 vs 73
  atoms). Both statements are true: encoding eta is cheap (0.03 s of emission), generating it is
  expensive (32 pool attempts vs 0).
* **`QIDKIZ_comp_0` is described as "1641.4 s, still unexplained" in `docs/agentic-notes/v0.4.5/PERF_v0.4.5.md`.** This
  lane explains ~57.6 s of that (its bare encode), and a round trip runs at least two encodes —
  which is ~115 s, not 1641 s. The remainder is still unattributed.

## Open questions / for the next agent

1. **`_rescue_unusable_perception` sweeps up to eight more charges through `AC2mol`.** On a ligand
   like these that is 8 x ~15 s. It did **not** trigger for either molecule measured, but it is the
   same cost multiplied and is the best suspect for the `encode_fail` timeout cohort — **48
   molecules, 11.0% of the capstone gap**. Nobody has profiled a molecule that reaches it.
2. **Is a generated conformer's ligand AC byte-identical to the input's?** If yes, **every SL1
   re-encode in a round trip becomes a memo hit** and the cross-encode reuse measured here (34.30 s
   -> 14.79 s -> 15.73 s, `max_weight_matching` 26 668 -> 0) applies to the generation path too. The
   mechanism is proven for a repeated *identical* input; **the generated-side case is untested** and
   it is the highest-leverage open question in either perf lane.
3. **The corpus-wide distribution of over-cap ligands is now partly known and should be finished.**
   valsearch measured ~0.2% of molecules (0/100 seed 1; 4/1992 seed 7) and 1.3% in the slow cohort,
   with strong enrichment — 4 of the 7 molecules `docs/agentic-notes/v0.4.5/ENCODE_FAIL_v0.4.5.md` lists as unresolved
   are over-cap (`KESWUB`, `BENVOG`, `HICLAG`, `HOHKUL`), while the other three (`FAQYUU`,
   `KEMTED`, `NAKLET`) are **not**, so their cost has a different, unidentified cause.
4. **The `nx.max_weight_matching` replacement remains unexplored as a lever.** The graph has no
   weights, so it is maximum *cardinality* matching; a faster algorithm exists but can return a
   different equally-large matching. `OIN_VALENCE_MATCHER` exists as the lever hook and is default
   OFF. Any attempt needs a corpus A/B on emitted strings, not on timing.
5. **The whole shipped result rests on one ligand.** Re-run `perf_encode_profile.py` /
   `perf_encode_keys.py` on an over-cap molecule with a **different** ligand (`KESWUB`, `BENVOG`,
   `HICLAG` or `HOHKUL` from the valsearch scan) before treating "1.62-2.16x" as a general claim.
6. **The product question this lane surfaced is still open at the top level:** for these ligands
   perception *fails* and `best_BO` is a fallback guess. valsearch showed the guess stops changing
   after ~100 candidates. Whether the encoder should instead **decline** such ligands loudly — the
   position taken for boron/carborane cages, where "an encoder that refuses those inputs is behaving
   correctly, not failing" — has not been decided.
