# Encoder performance — v0.4.5 (`swimlane/v045-encspeed`)

Owner: encspeed swimlane, based on `main` @ `b23decb4`. The git-durable record; the
gitignored `spec/handoffs/v0.4.5/PROGRESS-encspeed.md` has the working narrative.

## Why this lane existed

Every performance wave in this project targeted **generation** — v0.4.0's P2–P11, v0.4.4's
`OIN_EARLY_EXIT`, v0.4.5's perf lane (which found a 48–57 s call running twice per rejected
conformer and halved it). **Nobody had ever profiled the encoder.** The perf lane then
measured that a single eta molecule's bare `XYZToSMILES().convert()` takes 46–71 s, and a
round trip runs at least two encodes, so such a molecule's floor was ~90–140 s before the
generator did anything.

## Headline: the cost is one function, and half of it is redundant

`perception_core.AC2BO` is **99.8 %** of a slow encode. Not "dominant" — essentially all of it.

| | `QIDKUL_comp_0` (eta, 59 atoms) | `QIDKIZ_comp_0` (non-eta, 39 atoms) |
|---|---|---|
| encode wall | 49.43 s | 57.58 s |
| `AC2BO` inclusive | **49.31 s (99.8 %)** | **57.47 s (99.8 %)** |
| `get_oin_string` (the entire OIN emission) | 0.03 s | 0.02 s |
| aligner (`_brute_force_symmetries`) | 0.01 s | 0.01 s |
| `ResonanceMolSupplier` | **0 calls** | **0 calls** |

**This is not an eta phenomenon.** The non-eta 39-atom molecule has the *same* profile and
is *slower*. It is also not a size phenomenon — 39 atoms, 57 s. It is a property of the
ligand's valence combinatorics.

### The mechanism

`perception_tmc.py::_select_lig_mol` is a **charge/carbene ladder**: it calls `AC2mol` → `AC2BO`
up to five times on the *same* adjacency matrix, varying only `charge` and `allow_carbenes`.
(`_rescue_unusable_perception` can add eight more arms on an AC of its own.)

For these ligands the per-atom valence product exceeds `_VALENCE_COMBO_CAP` (500 000), so
each call takes the `itertools.islice(product, _VALENCE_FALLBACK_TRIES)` branch and grinds
**exactly 20 000 candidate valence assignments, finds none valid**, and returns `best_BO`.

The five `AC2BO` calls for `QIDKUL_comp_0`:

| # | AC sha | charge | `allow_carbenes` | wall | candidates |
|---|---|---|---|---|---|
| 1 | `df3f6198` | −1 | True | 0.00 s | 1 |
| 2 | `6fb274fd` | −8 | True | 16.81 s | 20 000 |
| 3 | `6fb274fd` | −8 | **False** | 16.99 s | 20 000 |
| 4 | `6fb274fd` | **−6** | True | 16.73 s | 20 000 |
| 5 | `a914ed00` | 0 | True | 0.00 s | 1 |

Calls 2 and 4 differ **only in charge**. Only `BO_is_OK` and `charge_is_OK` read `charge` —
the candidate-generation half (`get_UA` → `get_UA_pairs` → `get_BO`) is charge-independent
and recomputes bit-identical results. Measured redundancy:

| what | total calls | distinct results | removable |
|---|---|---|---|
| `get_UA_pairs` (matching-running) | 60 001 | 26 668 | **55.6 %** |
| `get_bonds` by `(AC, UA)` | 180 005 | 4 007 | **97.8 %** |
| `get_BO` by `(AC, valences)` | 60 002 | 26 669 | 55.6 % |

## What changed

All in `src/oinsmiles/utils/perception_core.py`. Ungated: dead-work removal and memoization
are the two classes this codebase has previously proved byte-identical, and the byte-identity
A/B below covers 24 molecules.

**1. A bounded LRU memo for `get_UA_pairs` and `get_bonds`.**

The soundness argument is the coarse key: `get_UA_pairs`' result is a function of
`(AC, tuple(UA), tuple(du > 1 for du in DU))`, because `DU` is read **only** through the
`du > 1` predicate that decides how many virtual matching nodes get allocated. Two calls
agreeing on that key are handed the identical edge list in the identical insertion order,
and `nx.max_weight_matching` is deterministic on identical input — the memo can only make
the same answer arrive sooner.

Three details that are load-bearing rather than decorative:

* An LRU over adjacency **matrices** (6 slots, 200 000 total entries) rather than one slot,
  because a round trip re-perceives the same ligand connectivity many times — the input, then
  every generated conformer the SL1 accept-first check re-encodes — each a fresh array with
  identical contents.
* It is read only for an array the cache holds a **live reference** to (`slot["ac"] is AC`).
  `id()` is just a fast index; the identity test is what makes a recycled `id` unable to
  alias a stale entry.
* A hit returns a **fresh mutable** object. Callers mutate what they get: `get_UA_pairs`
  appends the virtual-node edges to `get_bonds`' list. A shared object would poison every
  later hit. `tests/unit/test_ac2bo_memo.py` pins this by poisoning a result and checking the
  next hit is clean.

**2. `charge_is_OK` — same logic, minus the interpreter overhead.** It was the second-largest
cost (14.4–16.7 s). Three sources, all removed:

* `list(BO[i, :]).count(1)` allocated a fresh 59-element list of boxed `np.int64` **per
  carbon, per call**. One `(BO == 1).sum(axis=1)` gives every atom's count at once.
* `get_atomic_charge` was one Python call per atom — **4 440 174 calls in a single encode**.
  Inlined as the identical `elif` ladder. Now 132 calls (from `BO2mol` only).
* `list(BO.sum(axis=1))` yields boxed `np.int64`, so every downstream `1 - v`, `v == 2` and
  `Q += q` went through numpy's scalar machinery. `.tolist()` unboxes once.

The carbon corrections stay a Python loop **on purpose**: `if ns == 3 and Q + 1 < charge`
reads the *running* total, so they cannot be reordered or vectorised without changing the
answer. Its dead `q_list` accumulator (built, never read) is dropped.

**3. `get_UA` / `valences_not_too_large`** — `.tolist()` over `list()`, subtraction done once
instead of twice, `append` bound outside the loop.

**4. `AC2BO`'s `elif` evaluates `charge_is_OK` in place** instead of eagerly on the line
above. It sat behind two cheaper short-circuiting predicates and was computed even when
`status` had already returned. Pure dead-work removal; its only side effect is a DEBUG log
line reachable solely on the `allow_carbenes=False` arm, so the emitted OIN is untouched.

### Measured NEGATIVE — do not re-chase

Recorded because each looked obviously right and was measured wrong. Same-process A/B on the
real 37-atom ligand AC, both implementations interleaved in one process so the ratio is
immune to load drift, and asserting the two agree on every input:

| change | speedup | verdict |
|---|---|---|
| `charge_is_OK` rewrite | **4.45x** | kept (agrees on all 320 cases) |
| `get_UA` + `.tolist()` input | **2.36x** | kept |
| `valences_not_too_large` `.tolist()` | 1.09x | kept (marginal) |
| numpy-vectorising `charge_is_OK` with `np.select` | **~1.0x** | **rejected** |
| numpy-vectorising `get_UA` (`flatnonzero`) | **0.98x** | **rejected** |
| `BO.sum() - AC.sum()` for `BO_is_OK`'s `check_sum` | **0.72x** | **rejected** |

The lesson is uniform: at ~40–60 atoms a numpy round trip costs more than the Python loop it
replaces. The first version of this lane's `charge_is_OK` *was* the `np.select` one, and a
whole-encode measurement showed 14.35 s → 14.17 s — no gain. Vectorising the small per-atom
loops is a reflex worth resisting here; unboxing `np.int64` is the change that actually pays.

## Results

### Deterministic counters (exact, contention-robust), `QIDKUL_comp_0`

| counter | before | after | change |
|---|---|---|---|
| `nx.max_weight_matching` calls | 60 001 | **26 668** | −55.6 % |
| `get_atomic_charge` calls | 4 440 174 | **132** | −99.997 % |
| `get_bonds` calls | 180 005 | 146 672 | −18.5 % |
| `get_bonds` inclusive wall | 4.10 s | **0.46 s** | −89 % |
| `charge_is_OK` inclusive wall | 14.35 s | **3.59 s** | −75 % |
| `BO_is_OK` inclusive wall | 8.45 s | **3.17 s** | −62 % |
| `get_UA_pairs` inclusive wall | 28.31 s | **14.86 s** | −48 % |
| peak RSS | — | 157 MB | memo is not a memory problem |

### Wall clock, with conditions stated

The host was **not quiet** — two sibling lanes were running throughout, load 2–22. Absolute
seconds are therefore not clean. Two measurements, both honest about that:

* Single runs: **49.43 s at load 2.4 (before)** vs **24.82 s at load 12.2 (after)**. The
  after-run was faster despite 5x the load, so the improvement is real and understated.
* **Interleaved paired A/B** — the load-robust measurement. Arms alternate
  BASE/NEW/BASE/NEW so host contention averages across both rather than favouring
  whichever ran first. Load 22–38 throughout, which roughly doubles both arms' absolute
  seconds; the *ratio* is what this establishes.

  | molecule | BASE | NEW | speedup |
  |---|---|---|---|
  | `QIDKUL_comp_0` (eta, 59 atoms) | 104.94 s, 107.10 s | 66.02 s, 64.66 s | **1.62x** |
  | `QIDKIZ_comp_0` (non-eta, 39 atoms) | 123.98 s, 128.13 s | 54.09 s, 62.38 s | **2.16x** |

  All eight runs re-confirmed byte-identity: `QIDKUL` sha `0d428d9dfc56` and `QIDKIZ` sha
  `796221298c6c` on every arm and repetition. That is a third independent confirmation,
  after the 24-molecule manifest and the 3-convert warm-memo check.

### Cross-encode reuse

Three consecutive `convert()` calls on `QIDKUL_comp_0` in one process, load 17:

| | wall | `max_weight_matching` calls | OIN sha |
|---|---|---|---|
| 1st | 34.30 s | 26 668 | `0d428d9dfc56e18e` |
| 2nd | **14.79 s** | **0** | `0d428d9dfc56e18e` |
| 3rd | 15.73 s | 0 | `0d428d9dfc56e18e` |

Peak RSS 157 MB across all three. This matters for the round trip, which re-encodes the same
ligand connectivity repeatedly (input, then every conformer SL1 checks).

### Byte-identity

`tools/enc_byte_identity_ab.py`, two-directional `git show HEAD:… > …` A/B (never `git
stash` — it is shared across worktrees), memo cleared between molecules so a cross-molecule
hit cannot be what makes the arms agree.

**24 molecules — 4 goldens (CisPlatin, TransPlatin, Ferrocene, POJJOP) + 8 more fixtures
(Cis-PtCl2(en), fac/mer-Ir(ppy)3, PdCl2-R-BINAP, Zeise's salt, YESKOZ, FeCO5, CuCN2) + 12
dataset molecules — 7 of them eta, 0 errors.**

Manifest sha256 on **both** arms:
`9b679638dfeb3cf81abaf563e2b59cced3e3a5a6f3e0005fb983c9f57c3b014b`

The only differing line between the two runs is the `--label` echo.

Suite: `discover tests/unit` **605 tests OK** (3 skipped, 3 expected failures) — matching the
documented baseline — plus 8 new guards in `tests/unit/test_ac2bo_memo.py`.
`tests/unit/test_regression_stability.py` green. `ruff check` / `format` clean.

## How close to 30 s, honestly

**The encode now fits under 30 s for these molecules; the round trip still does not.**

* `QIDKUL_comp_0` encode: ~25 s at load 12, ~15 s warm. Under 30 s.
* But a round trip runs **at least two** encodes plus generation. Cold + warm ≈ 25 + 15 =
  40 s of encode alone, before the generator starts. **A round trip for this molecule cannot
  hit 30 s**, and this lane does not change that.

### The class for which 30 s is unreachable, and the evidence

**Ligands whose per-atom valence product exceeds `_VALENCE_COMBO_CAP`.** For those, every arm
of the charge/carbene ladder runs `_VALENCE_FALLBACK_TRIES = 20 000` candidates and finds
none valid. What remains after this lane's work is the **charge-dependent** half —
`BO_is_OK` + `charge_is_OK` + `get_BO` over 60 000 candidates — which is ~13–15 s and is
*not* memoizable, because it is exactly the part that differs per arm. That is the floor.

The obvious lever, lowering `_VALENCE_FALLBACK_TRIES`, is **not byte-identical**: `best_BO`
is whichever candidate satisfied `BO.sum() >= best_BO.sum() and valences_not_too_large and
charge_OK` among those *tried*, so cutting the cap changes the perceived bond orders. It
cannot be tuned ungated. Same for replacing `nx.max_weight_matching`: the graph carries no
weight attributes, so it is really maximum *cardinality* matching and a faster algorithm
exists — but it can return a different equally-large matching, hence a different `BO`. Both
are real options **behind a lever**, with a corpus A/B; neither was attempted here.

Note the encoder's own honest framing: for these ligands perception **fails** (no valid Lewis
structure among 20 000 candidates) and `best_BO` is a fallback. The 20 000-candidate search
is buying a *guess*. Whether that guess is worth 15 s per ladder arm is a product question
this lane surfaces rather than answers.

## What was ruled out (with evidence)

- **`ResonanceMolSupplier` is not on the encode hot path** — 0 calls in a full `convert()` of
  either slow molecule. The v0.4.4 SL5 forked, CPU-bounded wrapper never fires here.
- **The aligner is not the encoder's cost.** `_map_to_template` 3 calls / 0.00 s,
  `_brute_force_symmetries` 1 call / 0.01 s. "Factorial in coordination number" is real in
  principle and 0.02 % of this encode in practice.
- **Eta canonicalisation is not the encoder's cost.** All of `get_oin_string` — eta winding,
  ring-rotation canonicalisation, inline emission — is **0.03 s**, 0.06 % of the encode. Eta
  molecules are slow to *generate*, not to encode.
- **Repeated `SanitizeMol` / `MolToSmiles` on an unchanged mol is not the cost** — 19 and 9
  calls, 0.00 s total.
- **`xyz2AC_obabel` is not the cost** — no calls; AC comes from `GetAdjacencyMatrix`.
- **Memoising whole `AC2BO` results wins nothing** — 5 calls, 5 distinct argument tuples.
- **`get_BO` memo declined** despite a 55.6 % hit rate: caching 26 669 N×N integer arrays is
  ~370 MB for ~2.3 s. The cost/benefit is wrong.
- **Encoder slowness is not eta-specific and not size-driven** — the non-eta 39-atom molecule
  is slower than the eta 59-atom one, with an identical profile.

## Open, not chased

- `_rescue_unusable_perception` sweeps up to **eight more charges** through `AC2mol`. On a
  ligand like these that is 8 × ~15 s. It did not trigger for either molecule measured, but
  it is the same cost multiplied and is a good suspect for the `encode_fail` timeout cohort
  (48 molecules, 11.0 % of the capstone gap).
- Whether a generated conformer's ligand AC is byte-identical to the input's, which would
  make every SL1 re-encode in a round trip a memo hit. The cross-encode reuse above proves
  the mechanism works for a repeated *identical* input; the generated-side case is untested.
- Corpus-wide distribution of "over `_VALENCE_COMBO_CAP`" ligands. This lane measured two
  molecules deliberately (no corpus sweeps), so the *size* of the affected class is unknown.
