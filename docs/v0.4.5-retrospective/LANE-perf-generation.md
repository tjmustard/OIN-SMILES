# Lane — generation-side performance

Finding and removing the largest measured cost in the round-trip slow tail: a full re-perception of
the generated geometry that was running twice per rejected conformer.

## ELI5

Turning an OIN string back into 3D coordinates works by guessing a shape, checking it, and guessing
again until one passes. The check is expensive — it re-reads the guessed molecule from scratch and
re-derives its string, which on a hard molecule takes about fifty seconds. This lane found that the
same check was being run **twice on the same guess**: once while the guess was being collected, and
again immediately afterwards while the collected guesses were being ranked. The second run could
only ever recompute an answer already known, so it was removed by remembering the first answer. That
exactly halved the expensive calls on the test case, and the emitted strings are byte-for-byte
identical before and after — which for this project is the only acceptable outcome for a speed
change.

## The work, visually

```
  ROUND TRIP, DEFAULT PATH (OIN_EARLY_EXIT=1, default-ON since v0.4.4)
  =====================================================================

   OIN string
       |
       v
   MetalloGenAdapter.generate()
       |
       +--> generate_3d_structures(..., accept_fn=<closure>)      generator3d/__init__.py
       |        |
       |        |  POOL FILL LOOP   (target_pool: eta 32 / non-eta 10)
       |        |  for attempt i in 0..target_pool-1:
       |        |      embed  ---------------------------- ~1-2 s   <-- CHEAP (measured)
       |        |      _try_accept(conformer)
       |        |          accept_fn(m)
       |        |              _reencode_key_matches(...)
       |        |                  _reencode_oin(m)
       |        |                      XYZToSMILES().convert()  ### 48-57 s ###  [CALL 1]
       |        |          if matched:  return [early_hit]   <-- SHORT-CIRCUITS the whole fill
       |        |                                                generator3d/__init__.py:618
       |        v
       |    mols  (every member here is a PROVEN NON-MATCH -- a match would have
       |           returned early above and never reached this list)
       |        |
       +--> _select_by_geometry(parsed, mols, early_exit=True)     metallogen_adapter.py
                |
                |  its OWN early-exit RE-SCAN over the same pool:
                |      _reencode_key_matches(...)
                |          _reencode_oin(m)
                |              XYZToSMILES().convert()  ### 48-57 s ###  [CALL 2]
                |                                            ^^^^^^^^^^^^^^^^^^^
                |                                            recomputes a KNOWN answer
                v
            selected conformer -> GeneratedStructure(xyz, mol)


  THE FIX  (commit c89f0bf7) -- one dict per generate() call, keyed on id(m)
  ==========================================================================

   MetalloGenAdapter.generate():
       _reencode_cache: dict = {}                       metallogen_adapter.py:1900
              |                       |
     passed into accept_fn      passed into _select_by_geometry(..., reencode_cache=)
       (POPULATES it)                (READS it back)
              |                       |
              v                       v
        [CALL 1] MISS  ~56 s     [CALL 2] HIT  0.04-0.05 s
              \___________  ___________/
                          \/
              _reencode_key_matches(..., cache=None)  <- default: byte-identical to pristine
              cache MISS  ==> recomputes exactly as before (same fn, same inputs)


  WHY OBJECT IDENTITY SURVIVES (and where it deliberately does not)
  =================================================================
   default FF-only path (optimizer=None)  : _try_accept returns the SAME mol object it
                                            appends to successful_mols  --> id(m) matches --> HIT
   optimizer=xtb / g-xTB                  : ASEOptimizer.optimize builds NEW objects --> MISS
   stereo_rejects fallback pool           : never seen by accept_fn at all         --> MISS
   both non-default cases degrade to "recompute", i.e. pre-fix behaviour. Never a wrong answer.


  THE ETA TAIL, MEASURED (successor work, v0.4.6, counter `pool.attempts_spent`)
  ==============================================================================
    molecule            attempts  accepted  target_pool
    Ferrocene   (eta)      32        28         32        <- runs the ENTIRE pool, never
    CisPlatin (non-eta)     0         1         10           short-circuits
                                                          <- accepts on attempt 0

    WHY: accept_fn is handed RAW pool conformers. The round-trip key only matches AFTER
         optimization, and optimization happens AFTER the fill loop. An eta ring's winding
         is not right until relaxation, so at the moment accept_fn sees it, the conformer
         genuinely is not yet acceptable.
    => NO acceptance-predicate change can shorten the eta fill. The widened pool is doing
       real work; the runtime is the price of correctness.

  LEGEND
    ###  ...  ###   the dominant measured cost
    [CALL n]        the two sites the lane collapsed into one
    -->             control flow / consequence
    <--             annotation on the line above
```

## Initial assumptions and hypothesis

The lane inherited a named hypothesis from the release handoff, and it was wrong about the
mechanism.

1. **Goal:** per-molecule round trip **under 30 s** without moving accuracy. Baseline from the
   6,719-molecule capstone: `p50 = 29.1 s`, `p90 = 204.9 s`, `p95 = 308.1 s`, `max = 1803 s`, with
   only **51.1%** of molecules already finishing under 30 s.

2. **Why the goal was believed to be accuracy-aligned:** `hard_fail` is **306 of the 436 gap
   molecules (70.2%)** and **49.3% of those are timeout-shaped (>=290 s)**. Speed converts a timeout
   into either a pass or a fast, informative failure. This framing held up.

3. **The handoff's "prior art" attribution** named the haptic machinery as the big unclaimed prize:
   the `scales_for_haptic` 4x sweep (`generator3d/embed.py:1094`, `:1208`), the widened winding pool
   (`metallogen_adapter.py`), and the C++ batch path refusing haptic (`embed.py:1297`), with an
   estimated upper bound of **~23% off total runtime**. The supporting statistic was real and still
   is: **eta is 23.3% of molecules but 35.6% of wall clock**, and ~3x slower at equal atom count.

4. **The lane's own working assumption** was therefore that it would be tuning haptic embedding. It
   started by tracing per-call costs on a single eta molecule to confirm where the time went, which
   is the step that overturned the premise.

5. **A standing methodology assumption, stated up front and inherited from `tools/v045_state.sh`:**
   **wall-clock timing is meaningless above roughly load 12 on this machine.** Measure deterministic
   counters. `tools/perf_attribute.py`'s module docstring says so in its second paragraph.

## What was actually found

### Confirmed

* **The dominant slow-tail cost is a duplicate expensive re-encode, not the haptic machinery.**
  `metallogen_adapter.py`'s SL1 accept-first early exit (`OIN_EARLY_EXIT`, default-ON since v0.4.4)
  tests every newly embedded conformer against the requested OIN's canonical key via
  `_reencode_key_matches`, whose confirming step `_reencode_oin` is a **full
  `XYZToSMILES().convert()` re-perception of the generated geometry — measured at 48-57 seconds per
  call** on `QIDKUL_comp_0` (eta, 59 atoms). That runs once inside the pool-fill loop's `accept_fn`
  per accepted conformer, and then **again**, redundantly, inside `_select_by_geometry_impl`'s own
  early-exit re-scan over the same returned pool.

* **The second call is provably redundant, by a structural argument over every path between the two
  sites.** A conformer that HAD matched would have short-circuited the whole pool fill and been
  returned alone (`return [early_hit]`, `generator3d/__init__.py`), so anything surviving to
  `_select_by_geometry` is *provably an already-confirmed non-match*. The re-scan recomputes a known
  answer. (The only mols in that list not covered by this argument are ones drawn from the untested
  `stereo_rejects` fallback, which were never tested by `accept_fn` in the first place and are
  handled by a clean cache miss.)

* **The embed itself is cheap.** Direct per-call tracing on `QIDKUL_comp_0` shows the embed and
  scale sweep cost **~1-2 s per pool-fill attempt**. The ~50-210 s per attempt comes almost entirely
  from the acceptance-check re-encode.

* **Measured effect of the fix** (`QIDKUL_comp_0`, bounded `max_attempts=2`):

  | what | before | after | note |
  |---|---|---|---|
  | `reencode_oin_full_calls` | 4 | **2** | exact halving, deterministic counter |
  | cost of the 2 eliminated calls | ~56 s each | **0.04-0.05 s** (cache hit) | same session, same code path |
  | `generate_s` | 212.8 s | 116.0 s | contended host (load ~7-13); same-code-version before/after, **not** cross-host |

  A second independent bounded run reproduced `reencode_oin_full_calls = 2` post-fix.

* **Byte-identity confirmed, two-directional `git show` A/B:** CisPlatin, TransPlatin, Ferrocene,
  POJJOP and `QIDKUL_comp_0` all match SHA256 before and after.

* **Suite:** full `discover tests/unit` (605 tests) green, `tests/integration/test_roundtrip_smoke.py`
  green, `ruff check` / `format` clean.

* **The duplication is not eta-specific in mechanism.** It fires for any molecule on the default
  `early_exit=True` path with a non-trivial rejected-conformer count. Eta molecules generate more
  rejected conformers before a winding-correct one appears, so they pay it most often and most
  visibly.

* **The eta correlation is real and is NOT a size artifact — it is the opposite.** Measured over the
  **634 succeeding molecules** of the 936-molecule re-baseline:

  | atoms | eta median | non-eta median | ratio |
  |---|---|---|---|
  | 0-50 | 7.9 s (n=47) | 3.3 s (n=86) | **2.4x** |
  | 50-80 | 15.5 s (n=91) | 5.0 s (n=182) | **3.1x** |
  | 80-120 | 49.7 s (n=31) | 8.1 s (n=176) | **6.2x** |

  Eta molecules are **smaller** on median (62 vs 73 atoms) while being **3.2x slower overall**, and
  the penalty holds inside every band. Size does not explain it; controlling for size *strengthens*
  the effect.

* **Eta enrichment in the tail, by marker presence** (deliberately a string-presence count, not a
  timing measurement — see the load rule below): `<=30 s` 530 molecules, **19.8%** carry an eta
  winding marker; `>30 s` 104 molecules, **63.5%**; `30-60 s` 51 at 58.8%; `60-120 s` 30 at
  **76.7%**; `>120 s` 23 at 56.5%. A **3.2x enrichment** in the tail.

* **The runtime distribution today**, over the 634 succeeding molecules of the re-baseline:
  **530 (83.6%) complete in <= 30 s**, median **6.9 s**, p90 **50 s**, p95 **101 s**, max **277 s**;
  **104 exceed 30 s**. (Measurement note that costs a whole analysis if missed: `elapsed_s` lives
  **inside** each report's `metrics` dict — `tools/test_dataset_roundtrip.py:766` writes
  `report.setdefault("metrics", {})["elapsed_s"]`. Reading it from the top level of a report
  silently yields 0 for every row.)

* **The eta mechanism, finally measured with a deterministic per-attempt counter**
  (`pool.attempts_spent`, recorded at *every* return of `generate_3d_structures`, not just the
  last):

  | molecule | attempts | accepted | target_pool |
  |---|---|---|---|
  | Ferrocene (eta) | **32** | 28 | **32** |
  | CisPlatin (non-eta) | **0** | 1 | 10 |

  The eta molecule runs the **entire** pool and never short-circuits; the pool is itself widened
  (`ETA_SELECT_POOL = 16` vs `DEFAULT_SELECT_POOL = 5`, doubled for the UFF pre-pool -> 32 vs 10).
  So it pays ~32 attempts where a comparable molecule pays 1.

* **And the root cause of that, which closes the item:** `accept_fn` is handed **RAW** pool
  conformers, and the round-trip key only matches **AFTER** optimization — relaxation happens after
  the fill loop. For CisPlatin the raw conformer already satisfies the key, which is why it accepts
  on attempt 0. For an eta molecule the ring winding is not right until relaxation, so at the moment
  `accept_fn` sees the conformer it genuinely is not yet acceptable. **Therefore no
  acceptance-predicate change can shorten the eta fill.** The widened pool exists so that *after*
  optimization at least one member has the requested face. **The eta runtime cost is the price of
  correctness, not a defect.**

### Refuted

* **REFUTED: the slow tail is the haptic machinery** (4x `scales_for_haptic` sweep, pool widening,
  no C++ batch path for haptic), estimated at ~23% of total runtime. The embed and scale sweep are
  ~1-2 s per attempt. Those costs are real and still unclaimed, but they were not the dominant term.
* **REFUTED: this fix reaches the 30 s target.** It does not, and the lane says so in its own
  document. It removes the single largest identified redundancy — roughly halving the
  reencode-dominated cost — not the whole tail.
* **REFUTED (successor work): the eta tail is closable by cheap means.** Six hypotheses, each killed
  by one measurement: low acceptance rate -> pool widening -> cost-per-attempt (inferred from the
  ratio-vs-size slope) -> attempt-driven (correct) -> selection-side predicate alignment (correct
  but ineffective) -> fill-loop `accept_fn` (unsound, per the raw-conformer finding above). The
  `OIN_ETA_EARLY_EXIT` lever was implemented, **fires**, and moves the attempt count by **zero**
  (32 -> 32) because `generate_3d_structures` fills the entire pool before `_select_by_geometry` is
  ever called; a selection-side early exit can only shorten the selection *scan*, which is cheap
  beside 32 embeds. It stays default-OFF and documented as ineffective.
* **REFUTED: "more compute" is the main accuracy lever.** 24 timeout molecules (seed 42, stratified
  12 `UFF_1` + 12 `g-xTB_1`) re-run on the cheaper `--quick` path: 6 SUCCESS (25%), 6 string
  mismatch, 6 atom-count mismatch, 6 MetalloGen failed, **0 timed out again**. Only ~25% is
  compute-limited; more compute buys ~44 of 936 molecules (~4.7%), not the 174 the timeout count
  implies.

## What was done

Lane branch `swimlane/v045-perf`, based on `e8b603d5`. Git-durable record:
`docs/PERF_v0.4.5.md`. Note that `git log --oneline main..swimlane/v045-perf` returns **nothing** —
the branch is fully merged into local `main`, so its tip is the merge-base. The lane's commits are:

```
git log --oneline e8b603d5..swimlane/v045-perf
1e5399c8  docs(perf): record the completed byte-identity confirmation and honest scope
c89f0bf7  perf: memoize the SL1 accept-first re-encode across pool-fill and re-scan
d7063970  tools(v0.4.5): parameterized, version-controlled sweep runner
572fda1f  docs(v0.4.5): measured encoder instability under input atom renumbering
20044883  tools(v0.4.5): rotation/renumbering canonicality probe (replaces the failed trust gate)
7b85e123  tools(v0.4.5): frozen, reproducible sweep-cohort builder
19d20042  tools(v0.4.5): generator-free canonicality A/B instrument
```

The last five are Wave 0's instrument commits, which the perf lane's branch carries because it was
cut after they landed. The perf change proper is `c89f0bf7` + `1e5399c8`.

### The change, `src/oinsmiles/generation/metallogen_adapter.py`

* `_reencode_key_matches(..., cache=None)` — new optional parameter. When provided, the expensive
  `full = _reencode_oin(m)` step is memoized in `cache`, keyed on `id(m)`.
* `_select_by_geometry_impl` / `_select_by_geometry(..., reencode_cache=None)` — thread the cache
  through to their own `_reencode_key_matches` call.
* `MetalloGenAdapter.generate()` — creates **one** `_reencode_cache = {}` per `generate()` call,
  passes it into the `accept_fn` closure (which populates it) and into
  `_select_by_geometry(..., reencode_cache=_reencode_cache)` immediately after (which reads it back).

### Why it is byte-identical, and why it landed **ungated**

* `_try_accept` returns the *same* mol object it appends to `successful_mols`, and `accept_fn` is
  called on that exact object, so object identity survives from the fill loop into the `mols` list
  for the default FF-only path (`optimizer=None`).
* Default `cache=None` at any call site that does not opt in — behaviour unchanged whenever caching
  is not wired up.
* A cache **miss** recomputes exactly as before, same function, same inputs. The memo can only make
  an answer arrive faster, never make it different.
* It degrades gracefully wherever identity does *not* survive: the `optimizer` (xtb/g-xTB) path
  constructs new mol objects via `ASEOptimizer.optimize` so `id(m)` will not match, and the
  `stereo_rejects` fallback mols were never tested by `accept_fn` at all. Both just recompute.
* Same pattern as two previously validated memos in this codebase: the PuLP topology memo
  (`compute_chg_and_bo_pulp.py`) and `alt_cache` in `generator3d/embed.py`.
* **Ungated** under the release's own rule: *a pure dead-work removal that is provably
  byte-identical may land ungated.* The contrast the lane document draws is worth keeping —
  reordering `scales_for_haptic` would change **which conformer wins** and must be gated if ever
  attempted. Memoizing a pure function is not in the same category as reordering a search.

### The tools, and exactly what each measures

**`tools/perf_attribute.py`** — deterministic-counter attribution for one molecule's generation
step. Monkeypatches the hot call sites to **count, not to change behaviour**; counters are exact
regardless of host contention. It counts `AllChem.EmbedMolecule` (total and rc == -1),
`embed.get_embedding` (== outer pool attempts consumed), `embed.get_embeddings_batch` (should be 0
at `num_threads=1`), CBC/PuLP `actualSolve` (real solver spawns), PuLP memo hits/misses via
`pulp_cache_stats()`, final pool size, and the SL1 re-encode path — `build_contract_mol`,
`_reencode_key_matches`, `_reencode_oin_fast`, `_reencode_oin`. It also prints a per-call stderr
timing line for each of those, so **a single run shows exactly which calls are expensive and which
are cache hits** — that is how the 48-57 s figure and the 4 -> 2 halving were both read off one
invocation. `--max-attempts` is the preferred triage knob because it bounds call **count**, not wall
clock, so it is contention-robust.

**`tools/perf_byte_identity_ab.py`** — runs the full XYZ -> OIN -> XYZ pipeline for one molecule and
prints a sha256 of the generated XYZ, so two revisions can be diffed. Deterministic given a fixed
seed (42, project convention) and a fixed `--max-attempts`. The docstring states why capping
attempts is legitimate *for this specific change*: a memo keyed on object identity either reproduces
a prior result exactly or, on a miss, recomputes it exactly — it never alters which attempt wins or
how many run. That argument does not transfer to a change that alters the search.

## Dead ends, refutations, and instrument failures

### The handoff's eta attribution was the lane's own starting premise, and tracing killed it

Named costs: the 4x `scales_for_haptic` sweep, `ETA_SELECT_POOL` widening, and no C++ batch path for
haptic; upper bound ~23% off total runtime. Direct per-call tracing on `QIDKUL_comp_0` put the
embed/scale sweep at ~1-2 s per attempt. The lane document is explicit that those costs remain real
and unclaimed — they were simply not the dominant term, and the bigger confound had to come out of
the measurement before they could be evaluated at all.

### The stratified 16-molecule sweep was started and deliberately abandoned after one molecule

Each bounded attempt still costs 1-3 minutes post-fix, and the mechanism was already established by
**structural reading of every code path between the two call sites**, not by sampling. This is a
deliberate cost/evidence trade, not an oversight — but it means there is **no per-molecule table
across the 8-eta / 8-non-eta sample**. `run_attrib_final.sh` (referenced in the lane's gitignored
progress file) exists to produce one.

### Violating the load rule: a timing probe run concurrently with a 6-shard sweep

The standing rule is `tools/v045_state.sh`'s: **wall-clock timing is meaningless above roughly load
12 on this machine; measure deterministic counters.** The 24-molecule quick-mode timeout probe was
run alongside the 6-shard 5k sweep at **load 21-33**. Its pass/fail outcomes survived (and were the
point), but **its `elapsed_s` values are unusable** and the timing half of the result was lost. The
same trap then recurred when characterizing the eta tail, which is why that measurement was
deliberately reduced to *marker presence* — identical on an idle and a loaded machine.

### Over-applying the load rule in the opposite direction

A **byte-identity A/B was deferred** on the grounds that it needed an idle machine. That was wrong:
**string equality is deterministic and load-independent.** A contended host makes it *slow*, not
*wrong*. The corrected sequencing, recorded in the successor doc: run the byte-identity A/B **now**
because it is load-independent; let the sweep finish before quoting any absolute timing.

Both errors are worth carrying forward together, because they are the same misjudgement in mirror
image: **classify each measurement as load-sensitive or load-independent before deciding when to
run it.** Counters, string equality, marker presence and parse outcomes are load-independent.
Seconds are not.

### The one-interpreter timing A/B

While attributing `VIMZUN` (the single undetermined re-baseline regression), the first A/B ran both
arms in **one interpreter**. The second arm inherited warm caches and the `AC2BO` memo and reported
levers-ON as **faster: 42.8 s vs 54.4 s**. Meaningless — and note the direction: the artifact
flattered the change under test. A valid separate-process A/B needs ~4 x 48 s on a host too
contended for the number to mean anything, so `VIMZUN` was recorded as **UNDETERMINED**, neither a
regression nor an artifact.

### Byte-identity probe hygiene, learned the hard way

* **Never `git stash`** for an A/B — the stash is shared across worktrees and a sibling lane will
  collide with it. Use `git show HEAD:<path> > <path>` with an `EXIT` trap that restores the working
  copy.
* **Python block-buffers stdout into a pipe.** A `timeout` kill therefore discards the buffered
  lines, and a downstream `sort` still exits 0 — producing an **empty file that looks like
  agreement**. Write directly to the output file with `buffering=1` and end with a `#DONE <n>`
  sentinel, so "incomplete" is immediately distinguishable from "agreed". A later full-encoder
  byte-identity run was in fact killed mid-flight at 6 of 61 fixtures and was caught precisely
  because the sentinel was absent.

### The 30 s target: what this lane can and cannot reach

Stated honestly in the lane document rather than claimed: the encode side alone measured
**46-71 s** for `QIDKUL_comp_0` across repeated runs, and a round trip runs at least two encodes.
**An eta molecule whose bare encode exceeds 30 s cannot hit a 30 s round trip regardless of any
generation-side fix.** That observation is what spawned the encoder-speed lane (see
`docs/v0.4.5-retrospective/LANE-encoder-speed.md`).

### Discrepancies between the sources, flagged rather than smoothed

* **`docs/PERF_v0.4.5.md` cites `generator3d/__init__.py:590` for `return [early_hit]`.** In current
  `main` that statement is at **line 618** (line 441 of the same file carries the explanatory
  comment). The structural claim is unchanged; the line number has drifted with later commits.
  Likewise `V045_STATUS_2026-07-25.md` cites `metallogen_adapter.py:1311` for the widened winding
  pool; the pool-sizing block is now around **:1843**. Grep for the identifier, not the line.
* **"eta is 23% of molecules but 35.6% of CPU".** The repo's own figure is **23.3% of molecules /
  35.6% of wall clock** (`docs/V045_STATUS_2026-07-25.md`). The rounded "23%" also collides
  confusingly with the handoff's separate *"~23% off total runtime"* upper-bound estimate. They are
  two different numbers about the same subject and should not be conflated.
* **`docs/PERF_v0.4.5.md` lists `QIDKIZ_comp_0` as "39 atoms, non-eta, 1641.4 s — still
  unexplained"** and proposes it as a test of whether the same duplication mechanism fires on a
  non-eta molecule. The encoder-speed lane subsequently profiled that exact molecule and found its
  **bare encode** takes 57.58 s with a profile identical to the eta case, which partly explains it —
  but the 1641.4 s figure itself was never re-attributed. Treat it as open.
* **Pool sizes.** The successor doc quotes `target_pool` 32 (eta) vs 10 (non-eta). The constants in
  `metallogen_adapter.py` are `ETA_SELECT_POOL = 16` and `DEFAULT_SELECT_POOL = 5`; the quoted
  figures are those doubled for the UFF pre-pool (`uff_pool_size = max(uff_pool_size, 2 * pool_n)`).
  Not a contradiction, but the two numbers describe different objects.

## Where it landed

Code: `src/oinsmiles/generation/metallogen_adapter.py` — `_reencode_key_matches(..., cache=None)`,
`_select_by_geometry_impl` / `_select_by_geometry(..., reencode_cache=None)`, and the per-`generate()`
`_reencode_cache` dict. Merged into local `main`. **No env lever** — pure dead-work removal.

Docs: `docs/PERF_v0.4.5.md` (this lane), `docs/GENERATION_PIPELINE.md` (pipeline context),
`docs/V046_HFAITHFUL_FINDINGS.md` (the six-hypothesis eta sequence and its resolution).

Counter: `pool.attempts_spent`, emitted at **every** return of `generate_3d_structures`
(`src/oinsmiles/generator3d/__init__.py`). Keep it. It is what made each eta refutation cost one A/B
instead of a cycle of argument.

Lever: `OIN_ETA_EARLY_EXIT`, **default OFF and documented as ineffective** in
`src/oinsmiles/oin/levers.py::_HELD_OFF`. Kept because it is correct and harmless, and because it
marks exactly where the boundary between the two mechanisms lies.

### `tools/perf_attribute.py` — deterministic counters for one molecule's generation

```bash
PYTHONPATH=src .venv/bin/python tools/perf_attribute.py \
    --dataset /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset \
    --molecule QIDKUL_comp_0 --max-attempts 2 --timeout 600 \
    --out /path/perf_attrib.jsonl
```

Prefer `--max-attempts` over `--timeout` for triage: it bounds call count, which is contention-robust.
Read the `[progress] _reencode_oin (FULL XYZToSMILES) call#N took ...s` stderr lines — that is the
attribution.

### `tools/perf_byte_identity_ab.py` — prove a generation-side change byte-identical

```bash
M=QIDKUL_comp_0
DS=/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset
F=src/oinsmiles/generation/metallogen_adapter.py

# arm NEW (working tree)
PYTHONPATH=src .venv/bin/python tools/perf_byte_identity_ab.py \
    --dataset "$DS" --molecule $M --max-attempts 8 > /path/new.txt

# arm BASE -- two-directional git show, NEVER git stash (shared across worktrees)
cp "$F" /path/mine.py
trap 'cp /path/mine.py "$F"' EXIT
git show HEAD:"$F" > "$F"
PYTHONPATH=src .venv/bin/python tools/perf_byte_identity_ab.py \
    --dataset "$DS" --molecule $M --max-attempts 8 > /path/base.txt
cp /path/mine.py "$F"; trap - EXIT

diff /path/base.txt /path/new.txt      # must be empty apart from any label echo
```

## Open questions / for the next agent

1. **`QIDKIZ_comp_0` (39 atoms, non-eta) at 1641.4 s is still not fully attributed.** Its bare
   encode is 57.58 s (encoder lane), which accounts for ~115 s of a round trip, not 1641 s. It is
   also the best available test of whether the accept-first duplication mechanism fires on a non-eta
   molecule with an unusually large rejected-conformer count. Run `tools/perf_attribute.py` on it
   with `--max-attempts` bounded.
2. **The handoff's original eta costs remain unclaimed and are now measurable cleanly**, since the
   bigger confound is out of the measurement: the 4x `scales_for_haptic` sweep, `ETA_SELECT_POOL`
   widening (now understood as *irreducible via acceptance*, but not necessarily irreducible via
   construction), and the C++ batch path refusing haptic. Any change here alters **which conformer
   wins** and must therefore be **gated**, with a corpus A/B — unlike this lane's memo.
3. **The two genuinely expensive eta options, stated so nobody re-derives them:** make the embed
   produce the requested ring face *before* relaxation (construction over selection — three prior
   negative results in this project), or relax fewer candidates by scoring winding
   *pre*-relaxation, which needs a pre-relaxation winding predicate that does not currently exist.
4. **The per-molecule attribution table was never produced.** `run_attrib_final.sh` over the
   stratified 8-eta / 8-non-eta sample would give one, at ~1-3 min per bounded attempt. Everything
   this lane concluded rests on `QIDKUL_comp_0` plus a structural code reading.
5. **Whether a generated conformer's ligand AC is byte-identical to the input's** is untested, and
   it is the highest-leverage open question at the seam between the two perf lanes: if it is, every
   SL1 re-encode in a round trip becomes an `AC2BO` memo hit and the encoder lane's cross-encode
   reuse (34.30 s -> 14.79 s -> 15.73 s) applies to generation too.
6. **`VIMZUN` needs one quiet-host, separate-process, paired run** to be resolved either way.
7. **Do not re-derive the eta acceptance story.** Six hypotheses, six refutations, one line of
   resolution: `accept_fn` sees **raw, unrelaxed** conformers, so it cannot test a property the
   conformer does not have yet. Read
   `docs/V046_HFAITHFUL_FINDINGS.md` § "RESOLVED: the eta pool cost is irreducible via acceptance"
   before proposing anything in that area.
