# Wave 0 — the measurement instrument

How v0.4.5 built the ruler that gated every other lane, why the first ruler was thrown away
after failing its own trust gate, and what the replacement can and cannot see.

## ELI5

v0.4.5's headline claim was "the string the encoder emits is now canonical" — same molecule in,
same string out, no matter how the input file happens to be written. You cannot make that claim
with a ruler you have not checked, because a broken ruler does not look broken: it produces a
number that is plausible, quotable, and wrong. So the first job of the release was not to fix the
encoder, it was to build something that could tell whether a fix helped, and then to *test the
tester* against a result already known. The first instrument failed that test and was demoted; a
second, sharper one replaced it and became the gate every lane had to pass. Five separate
instrument defects were caught this way during the release, and each one would have shipped a
confident wrong number if the checks had not existed.

## The work, visually

```
                            THE QUESTION W0 HAD TO MAKE MEASURABLE
             "Does the encoder emit ONE string for one molecule, however presented?"
  ============================================================================================

  [1] PLANNED INSTRUMENT                    [2] TRUST GATE                  [3] VERDICT
  ----------------------                    --------------                  -----------
  tools/reencode_ab.py                      Arm A runs on UNMODIFIED        ***FAILED***
   * reads a finished sweep's reports        `main`. It must reproduce
   * re-encodes the stored pair              the 2026-07-15 capstone        structural
       input_xyz                            within noise:                    expected 1.04%
       structures/<mol>_generated.xyz          byte_exact 81.19%              measured  19.3%
   * NON-MUTATING (fresh results dir)         key_equal  12.32%               ~19x inflated
   * GENERATOR-FREE -> no timeout                                           key_equal
     confound                                                                 12.04% vs 12.32%
   * 6,404 stored *_generated.xyz                                             -> tracks fine
   * 0.178 s/structure  (~7 min, 6 shards)
   * all 828 key_equal mols (500+315+13)
     have a stored structure
                |
                |  ROOT CAUSE (read off the harness, not guessed):
                |  tools/test_dataset_roundtrip.py:186
                |      oin2_string = get_oin_string(mol_gen_bonded, xyz_coords)
                |                                   ^^^^^^^^^^^^^^ gen_result.mol
                |                                   the GENERATOR'S OWN bond graph
                |  only on exception:  xyz_to_smiles.convert(gen_xyz_path)
                |
                |  A stored .xyz on disk has no gen_result.mol, so re-encoding it must
                |  RE-PERCEIVE connectivity with xyz2AC_obabel. On a slightly distorted
                |  generated geometry that returns a genuinely different graph.
                |  => the tool measured perception drift, not serialization drift.
                v
  [4] REPLACEMENT INSTRUMENT — became the primary
  ---------------------------------------------
  tools/canonicality_probe.py
      take ONE input structure; hold the MOLECULAR GRAPH FIXED; re-present it:
          rotate    random PROPER rotation (det = +1; improper would mirror a chiral mol)
          renumber  permute the order atoms appear in the XYZ file
          both      renumber, then rotate
      expected answer is BYTE-IDENTICAL  <-- known ground truth, not an inferred baseline
      any difference is classified with the sweep's own taxonomy
          slot_renumber / winding_star_drift / fragment_reorder / rdkit_canonical
      and re-tested at comparison-KEY level (canonical_roundtrip_key)
                |
                v
  [5] PROMOTION GATE  (docs/agentic-notes/v0.4.5/PROMOTION_GATE_v0.4.5.md)  300 mols, seed 42, 2 trials/mode
  ----------------------------------------------------------------------------------
      arm                byte-stable          comparison KEY broken
      all 6 levers OFF   173/298  58.1%       60   20.1%
      all 6 levers ON    208/299  69.6%       16    5.4%
      delta              +11.5 pts (+35)      60 -> 16  (-73%)
      drift by transform renumber 176->107, both 165->119, rotate 0 in BOTH arms
                |
                +-- VETOES (each had to pass, or the lever does not ship)
                |     facmer/cis-trans stay distinct RAW and at KEY level ......... PASS
                |     goldens byte-identical on the levers-OFF path ............... PASS
                |     OIN_STABLE_STEREO mirrors still DIFFER (10/10) ............. PASS
                |     tools/geometry_tag_shift.py: [M_XXX] moved 0/298 ........... PASS
                v
      recommend all six default-ON; hold the three that ADD descriptors

  ============================================================================================
  WHAT EACH INSTRUMENT CAN AND CANNOT SEE
  ---------------------------------------
                                    reencode_ab   canonicality_probe   run_sweep.sh (harness)
   encoder canonicality (string)     ~ partial     + DIRECT             ~ diluted
   comparison-key drift              + tracks      + DIRECT             + but it IS the predicate
   generator changes                 - blind       - blind              + only place it shows
   generation timeout confound       + immune      ++ none at all       XX DOMINANT (67.4% of
                                                                          failures are timeouts)
   perception drift on a generated   XX CONTAMINATES   + impossible by  ~ present but folded
     geometry                           the number       construction
   new descriptors (OIN_EMIT_AXIAL)  - blind       - blind              + as loud false negatives
   needs the gitignored dataset      no (reports)  YES                  YES

  LEGEND
    +  measures it directly / correctly        ~  partial or diluted signal
    -  structurally blind to it                XX actively corrupts the number
    ++ the confound cannot exist here          *** a gate outcome
    [n] pipeline stage                         -> transition / becomes
```

## Initial assumptions and hypothesis

1. **The quantity to move was named before the wave started.** The capstone sweep (6,719
   molecules, 2026-07-15) reported **12.32% of molecules round-tripping to the same isomer but a
   different string** — the `key_equal` bucket. `oin/compare.py` was already canonical; the
   emitted string was not. v0.4.5's plan was to promote that machinery upstream into the encoder
   (`docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md`).

2. **The batch harness was assumed unusable as the primary instrument**, and this assumption
   survived. From the Y3 audit: **67.4% of round-trip failures are generation timeouts**, so any
   code change that shifts runtime moves the pass rate for reasons that have nothing to do with
   the notation. A canonicality lane measured through the harness would be measuring compute.

3. **The hypothesis for the instrument** was that a *generator-free* measurement existed for
   free, because the capstone had already stored everything needed: every report records an
   absolute `input_xyz`, 6,404 `*_generated.xyz` structures exist on disk, all 828 `key_equal`
   molecules (500 + 315 + 13) have a stored structure, and re-encoding both sides of a pair was
   measured at **0.178 s/structure** — about 7 minutes for all 6,404 at 6 shards. Freeze the
   geometry, vary only the encoder, and the timeout confound disappears.

4. **Assumed, wrongly, that "re-encode the stored generated structure" is what the sweep did.**
   This is the single load-bearing assumption of the whole wave and it was false. See below.

5. **Assumed the frozen `bucket_report.json` was a usable baseline.** Also false, four
   independent ways.

## What was actually found

### Confirmed

* **The generator-free idea is right; the specific implementation was not.** A geometry-frozen,
  generator-free measurement of the encoder is achievable and cheap. `canonicality_probe.py`
  delivers it at roughly 1 molecule/second-scale throughput with no generator invocation at all.
* **`reencode_ab.py`'s `key_equal` sub-split does track the sweep**: **12.04%** against the
  capstone's **12.32%**. That half of the tool was sound and it is why the tool was kept rather
  than deleted.
* **The encoder is fully orientation-invariant.** Across 225 molecules x 2 trials, `rotate`
  produced **exactly 0** drift; in the 300-molecule promotion sample it was **0 in both arms**.
  `_align_to_pai` does its job for orientation. Every defect W0 found is an *atom-numbering*
  dependence.
* **A soundness defect nobody had instrumented, found by the replacement instrument on its first
  real run** (`docs/agentic-notes/v0.4.5/RENUMBERING_INSTABILITY_v0.4.5.md`, `main` @ `20044883`, all levers OFF, 225
  molecules seed 42, 223 encoded):

  | | count | share |
  |---|---:|---:|
  | byte-stable across all transforms | 125/223 | **56.1%** |
  | drifted | 98/223 | **43.9%** |
  | of which the comparison KEY also changed | 47/223 | **21.1%** |

  Severity classes under pure renumbering: **29 molecules (13.0%) are stereo-only flips** —
  byte-identical once `@`/`@@` are deleted, so at least one of the two encodings is *wrong*; 79
  (35.4%) other drift; 1 (0.4%) changed its `[M_XXX]` geometry classification; 8 (3.6%)
  aromaticity collapse. Worked examples: `FEQFIS_comp_0` (verified a pure permutation — identical
  sorted coordinate multiset), `DUDREA_comp_0` (`[Y_SPY]` -> `[Y_TET]`), `CEBVIR_comp_0`
  (aromaticity collapses entirely).
* **Why no prior instrument could have found it.** The Y1/Y2/Y3 injectivity audit asked "does the
  encoder *separate* two enantiomers?" (mirror-twin probes: two *different* structures). It never
  asked "does the encoder *consistently report* one enantiomer?" (one structure against itself).
  The round-trip sweep cannot see it either: the generator builds coordinates *from* the OIN
  string, so both sides inherit whatever configuration the encoder chose on that input ordering —
  the error is common-mode and cancels.
* **The promotion gate result** (`trial/v045-integration`, `--n 300 --trials 2`, seed 42 so every
  arm samples the same molecules): byte-stability **58.1% -> 69.6%** (+11.5 pts, +35 molecules);
  comparison-key instability **60 -> 16 molecules**, a 73% reduction. Subclass movement:
  `rdkit_canonical` 91 -> 18, `slot_renumber` 42 -> 74, `encode_fail` 1 -> 0.
* **`slot_renumber` rising is reclassification, not regression** — established independently by
  Lane 2 with per-molecule accounting (0 molecules broken in either arm). A molecule drifting in
  *both* body and slots counts under `rdkit_canonical` (first-matching subclass); once the body is
  canonical its residual slot drift reclassifies. Byte-stability up while `slot_renumber` up is
  the expected signature.
* **Every veto passed.** fac/mer and cis/trans stay distinct raw *and* at key level with all six
  levers on; goldens byte-identical on the opt-out path; `OIN_STABLE_STEREO` mirrors differ 10/10
  (nothing collapsed); `geometry_tag_shift` **0/298** `[M_XXX]` changes, 0 coordination-number
  changes, 0 transitions.
* **The dataset needs deduping before any cohort is built.** 26,230 `.xyz` files, only **25,197
  unique basenames** — **1,033 names exist in BOTH `cat/` and `photo/`**. The harness keys reports
  by basename, so the raw tree makes those 1,033 double-match and race each other's report writes.
  (Two independent counts agree: `uniq -d` over basenames, and 26,230 - 25,197.
  `run_regression_sweep.sh`'s "239 names" is stale.)

### Refuted

* **REFUTED: re-encoding a stored generated `.xyz` reproduces what the sweep measured.** It does
  not, and the gap is large — `structural` **1.04% expected vs 19.3% measured**, roughly 19x. The
  harness never encodes the stored file.
* **REFUTED: the frozen `bucket_report.json` is a trustworthy baseline.** Wrong four independent
  ways: **stale** (29 `rmsd_gate` molecules were fixed a week *after* the snapshot was taken),
  **misattributed** (`HOCVAY`/`WEFZAL` are generation-side deaths bucketed as encoder failures),
  **hiding a regression** (`XOSTUW_comp_0` passed then and fails now), and **understating
  `atom_count`**. Consequence: only **per-molecule transitions** are reportable from it, never an
  aggregate delta. This is `tools/rebaseline_report.py`'s entire reason to exist.
* **REFUTED: the 13% stereo-flip class is downstream of order-dependent bond-order perception.**
  This was reasoned out and put in writing to the Lane 8 agent. Measured on the three worked
  examples, 3 trials x 3 transforms: levers all OFF 0/3 stable, 3 key-level defects;
  `OIN_STABLE_METAL_AC` alone 2 defects; `+ OIN_CANONICAL_PERCEPTION` 2; all four levers 2. Only
  the metal-AC fix closes anything. `FEQFIS_comp_0`'s stereo flip and `CEBVIR_comp_0`'s
  aromaticity collapse survive all four. Lane 8's independent work was therefore required, not
  redundant. (Scope caveat stated in the source: n=3, deliberately the hardest hand-picked cases —
  it refutes the specific claim, not Lane 1's general 6-fixed/0-regressed result over 250
  molecules.)
* **REFUTED: `pending_g-xtb` means "no xtb binary".** xtb **6.7.1 is installed** at
  `.venv/bin/xtb`; the status means the tier-2 g-xTB refinement is **queued or running**. The
  wrong reading was inferred from the status *name* and published as a pass-rate range.
* **REFUTED: `OIN_STABLE_METAL_AC` would reclassify geometries corpus-wide.** This veto was set
  *against* the author's own fix and it lost: 0/298.

## What was done

Five commits landed on local `main`, all tooling and findings — no encoder default changed in
Wave 0 itself:

| commit | artifact |
|---|---|
| `19d20042` | `tools/reencode_ab.py` — generator-free re-encode A/B (the instrument that failed its gate) |
| `7b85e123` | `tools/build_sweep_cohort.py` — frozen, reproducible cohort builder |
| `20044883` | `tools/canonicality_probe.py` — the corrected canonicality instrument |
| `572fda1f` | `docs/agentic-notes/v0.4.5/RENUMBERING_INSTABILITY_v0.4.5.md` — the soundness finding |
| `d7063970` | `tools/run_sweep.sh` — parameterized, version-controlled sweep runner |

Later in the release, still Wave-0 infrastructure: `tools/geometry_tag_shift.py` (the perception
veto), `tools/rebaseline_report.py` (per-molecule transitions), `tools/v045_state.sh` (state
computed from the repo).

### What each tool measures, precisely

**`tools/canonicality_probe.py`** — the primary instrument. Reads the corpus, dedups by basename
(first occurrence wins), samples `--n` with a hardcoded `SEED = 42` so arms are comparable, then
for each molecule encodes the original and `--trials` random variants of each of `rotate`,
`renumber`, `both`. `random_rotation()` builds a **proper** rotation via QR then forces `det = +1`
by flipping a *column* — never reflecting the result, because an improper operation mirrors the
structure and legitimately changes a chiral molecule's encoding. Drift is classified by
`_subclass()`, which is deliberately the *same* taxonomy as
`tools/roundtrip_bucket_report.py::_key_equal_subclass`, so a number from the probe is directly
comparable to a `key_equal` sub-split from a real sweep. Each drift is then re-tested with
`canonical_roundtrip_key` to separate string drift from isomer-level drift. Output JSON records
`commit_id` and **every `OIN_*` env var present**, so an arm's configuration is inside its own
artifact rather than in someone's shell history.

**`tools/reencode_ab.py`** — retained, demoted. Non-mutating sibling of
`tools/recalculate_oin_smiles.py` (which mutates source reports in place and applies RMSD/atom-count
geometry gates). This one writes a fresh results dir and deliberately drops the geometry gates,
stamping `status: success` whenever both strings exist so that
`roundtrip_bucket_report.classify()` routes purely on the STRINGS. Exclusions are counted, never
dropped: `skipped_no_input_xyz`, `skipped_no_generated_xyz`. It is a *lower bound* on the drift a
real sweep sees, because rows that timed out have no second geometry.

**`tools/geometry_tag_shift.py`** — the veto on `OIN_STABLE_METAL_AC`. That lever changes
**perception** (capping highest-Z first), and the change is *asymmetric*: capping the metal before
the light atoms bridging to it can only **keep** metal bonds the old atom-order iteration
discarded. More metal bonds -> higher coordination number -> the geometric template fit can land on
a different polyhedron. And `[M_XXX]` is not cosmetic: it selects the vertex table, hence the
rotation group, hence the canonical slot labelling, hence the comparison key's entire vertex
signature. The tool encodes each molecule twice (lever off, lever on), reports the `[M_XXX]`
transition matrix by direction, and reports coordination number as the count of **distinct** `{n}`
slots so a tag change can be attributed to a genuine new donor rather than a template tie flipping.
Its own printed guidance: "a tag change WITH a CN change is the lever keeping a real metal bond; a
tag change WITHOUT one is a template tie flipping — less defensible."

**`tools/build_sweep_cohort.py`** — materializes a cohort as a **deduped symlink directory plus a
committed manifest**, because of the 1,033 colliding basenames above, and because
`test_dataset_roundtrip.py --random` seeds from **system time** so two arms would not see the same
molecules. Seeds explicitly (default 42), samples from the *sorted* name list so the draw depends
only on `(names, seed)`, refuses to overwrite an existing cohort dir, and reports overlap against
prior results dirs — the only slice that is identical-molecule diffable. Produced the frozen
5,000-molecule cohort at `tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k` (+ `_manifest.json`),
overlapping **1,319** capstone and **790** v0.4.4-regression molecules.

**`tools/run_sweep.sh`** — replaces the one-off `run_*_sweep.sh` scripts previous releases dropped
into their gitignored results dirs, where they were lost. Resolves the interpreter in a stated
order (`OIN_SWEEP_PYTHON` -> `$REPO/.venv` -> a sibling checkout's), **pins `PYTHONPATH` to
`$REPO/src` explicitly**, records `run_config.json` (commit, cohort, molecule count, shards,
timeout, every `OIN_*` lever) beside the results, then rebuilds the summary and runs the bucket
report. Its header carries the timeout-budget warning in the first-person.

**`tools/rebaseline_report.py`** — per-molecule transitions against the capstone snapshot, for two
deliberately different populations: the **436 gap molecules** (was failing — how many now pass?)
and a **500-molecule seed-42 passing guard** (did anything REGRESS? — the half a gap-only re-run
structurally cannot see, and where `XOSTUW`-class regressions hide). `pending_g-xtb` gets its own
column: calling it a pass or a fail would both be dishonest.

**`tools/v045_state.sh`** — written after an account monthly spend limit terminated eight parallel
lane agents mid-task. Its lesson is in its own header: *a hand-maintained status file is stale the
moment work resumes*, so it **computes** state from the repo — branch tips, commits ahead of the
v0.4.5 base `e8b603d5`, uncommitted file counts per lane worktree, `WIP` detection from the tip
subject, every `OIN_(CANONICAL|STABLE)_*` lever grepped out of each lane branch's `src/`, active
`systemd --user` units, load, cores — and then prints the next action. `--next` prints only that.
It ends with a "RULES THAT BIT US ALREADY" block.

### The result this instrument fed

`docs/agentic-notes/v0.4.5/PROMOTION_GATE_v0.4.5.md`, reproduced in the diagram above: **58.1% -> 69.6%** byte-stability,
key instability **60 -> 16 molecules**, all vetoes green, recommend all six canonicality levers
default-ON via the `OIN_EARLY_EXIT` template (`metallogen_adapter.py:1636-1653`) and hold the three
levers that *add* descriptors (`OIN_EMIT_AXIAL`, `OIN_EMIT_BOUND_AMINE`, `OIN_RESCUE_STUCK_RING`).
The distinction that makes that coherent: the six promoted levers **repair a renumbered
presentation without rewriting the canonical answer** — which is exactly why the corpus shows no
churn. The three held levers **add information**, which converts a silent false positive into a
loud false negative — right direction, separate product call.

## Dead ends, refutations, and instrument failures

This is the core of the wave. Each is named as a failure *mode*, because the specific bug matters
less than the shape.

### FM-1. The trust gate fired — a wrong instrument caught before any lane read it

**What was built:** `tools/reencode_ab.py`, on the reasoning laid out under "Initial assumptions"
above.

**The gate:** Arm A had to run on unmodified `main` and reproduce the capstone's `byte_exact`
**81.19%** and `key_equal` **12.32%** within noise before *any* Arm B number would be believed.

**What happened:** `structural` came out at **19.3%** against an expected **1.04%** — roughly 19x
inflated. `key_equal` was fine (12.04% vs 12.32%).

**Root cause, read off the harness rather than guessed:** at
`tools/test_dataset_roundtrip.py:186` the harness computes the second OIN as
`get_oin_string(mol_gen_bonded, xyz_coords)` where `mol_gen_bonded = gen_result.mol` — the
**generator's own bond graph**, held in memory. Only on exception does it fall back to
`xyz_to_smiles.convert(gen_xyz_path)`. A stored `.xyz` on disk has no `gen_result.mol`, so
re-encoding it must re-perceive connectivity through `xyz2AC_obabel`; near the covalent-radius +
0.45 Å cutoff a slightly distorted generated geometry yields a genuinely **different graph**. The
excess is *connectivity-perception drift, not serialization drift* — a real problem, but a
different problem, and out of scope for the canonicality lanes.

**Why this is the best thing that happened in the wave:** the gate did exactly its job. It caught
a wrong instrument BEFORE any lane read a number off it. Had the gate been skipped, every
canonicality lane would have been optimising against a 19x-inflated `structural` bucket.

**Generalises to:** *never believe an instrument's Arm B until Arm A has reproduced a known result
on unmodified code.* And: *"my tool re-does what the pipeline does" is a claim about someone else's
code — go read that code.*

### FM-2. The serene `0/0` — a harness that measured nothing and reported a clean pass

`canonicality_probe.py`'s own first failure, and it cost a full A/B run. The dataset directory is
**gitignored**, so it exists in the main checkout and **not** in a fresh git worktree — and lanes
work in worktrees. With no files found, the probe printed `0/0 byte-stable, 0/0 drifted` and
exited 0. That does not read as a misconfiguration; it reads as a clean result. It was caught only
by noticing that an A/B had returned zeros for *both* arms.

Hardened: the probe now `sys.exit()`s with an explanatory message naming the gitignore cause and
the `--dataset` fix, and it also exits non-zero when `--only` matches nothing and when `--shard`
selects nothing. `geometry_tag_shift.py` carries the same guard.

**Generalises to, and state it as a rule:** *a harness that measures nothing almost always exits 0,
and an empty result set reads as agreement.* Every measurement tool needs a floor assertion — "n
must be > 0" — and every A/B needs its denominators eyeballed before its ratios.

### FM-3. The mixed-code-tree sweep — two `src` trees in one process

`tools/run_sweep.sh` had two defects at once, and the second is nastier than the first.

1. **No worktree has a `.venv`.** The project rule is to use the MAIN checkout's pinned venv
   (rdkit `==2025.9.3`) and never `uv sync` in a worktree, so `$REPO/.venv` is normally absent when
   sweeping from a lane tree — the script could not run there at all.
2. **Running main's copy of the script against another tree's `src` MIXED TWO CODE TREES.** The
   harness appends its *own* `../src` to `sys.path`, so a copy executed from a different checkout
   than the code under test silently imports from both. Shard logs showed both
   `OIN-SMILES/src` and `oin-v045-trial/src`.

Caught by checking `oinsmiles loaded from` in the shard logs about **30 s** after launch; the run
was killed and its output deleted. Both fixed: explicit interpreter resolution order with a loud
error, and `export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"` before launching shards.

**Generalises to:** *in a multi-worktree checkout, "which code am I running" is a measurement, not
an assumption.* Print the resolved module path in the first lines of any long job's log, and read
it before walking away.

### FM-4. The one-interpreter timing A/B — the second arm inherits warm caches

While attributing the single undetermined re-baseline regression (`VIMZUN`), the first A/B ran
**both arms in one interpreter**. The second arm inherited warm caches and the `AC2BO` memo, and
duly reported **levers-ON as faster: 42.8 s vs 54.4 s**. Meaningless — and note the direction: the
artifact flattered the change under test.

`VIMZUN` was recorded as **UNDETERMINED**, neither a regression nor an artifact, because a valid
separate-process A/B needs ~4 x 48 s on a machine too contended for the number to mean anything.
Naming it beat guessing either way.

**Generalises to:** *any timing A/B must be separate processes, and it should alternate arms
(BASE/NEW/BASE/NEW) so contention averages across both rather than favouring whichever ran first.*
Same-process is only legitimate for **agreement** checks (does the fast path return what the slow
path returns), never for ratios — except where the interleaving itself is the load defence, as in
the encoder lane's paired A/B.

### FM-5. The timeout-budget artifact — a warning walked into by its own author

`tools/run_sweep.sh`'s header says, in the first person, *do not compare pass rates across sweeps
with different budgets*, and cites v0.4.4's **11 apparent "regressions"** which were all 300 s
timeouts against an arm run under quick mode's cap — a config artifact, zero wrong answers.

The v0.4.5 re-baseline then ran at `--mol-timeout 300` against a capstone baseline that ran at
**1800**. A molecule legitimately needing 300–1800 s **passed then and cannot pass now regardless
of code quality**, and it surfaces as a REGRESSION.

Measured on the first five flagged rows, using the capstone's own `metrics.elapsed_s`: **303.6,
331.5, 344.6, 303.3 and 1117.2 s** — **5/5 already exceeded 300 s while passing; ZERO were
correctness regressions.** In the full re-baseline the same check resolved 8 of 9 flagged rows as
artifacts (`PUMLEP` 1117.2, `YOQMAT` 356.6, `CALPOX` 344.6, `AXIDUH` 331.5, `REXROC` 324.5,
`AHUKOF` 303.6, `LUYYOT` 303.3, `XEDNUQ` 302.5) and left `VIMZUN` (100.0 s) undetermined.

The warning now lives as a footer block in `tools/rebaseline_report.py:158-174`, ending: *"I wrote
that warning and then walked into it. So: for every REGRESSED row, check the capstone `elapsed_s`
before believing it."*

**Generalises to:** *a regression whose old `elapsed_s` exceeds the new budget is an artifact, not
a defect* — and writing a warning down is not the same as being protected by it. Also note the
adjacent trap: `tools/injectivity/missed_success_audit.py` hardcodes `TIMEOUT_S = 300.0`, so a
sweep at any other budget is mis-attributed by it.

### FM-6. Misreading `pending_g-xtb` — a pass rate published from a sweep still executing

`pending_g-xtb` was read as "the tier-2 g-xTB refinement was skipped because there is no `xtb`
binary", and a headline range was published on that basis. Wrong: **xtb 6.7.1 is installed** at
`.venv/bin/xtb`, an `xtb struc.xyz --gxtb --opt` process was live at the time, and the pending
count was **draining** (109 -> 88 while it was being looked at). The status means tier 2 is
**queued or in progress**.

The compounding failure: **936/936 individual reports existed while tier 2 was still running.**
Report count is not completion.

The rule, now written into `tools/rebaseline_report.py`'s `verdict()` comment and printed by the
tool itself: **poll `systemctl --user is-active <unit>` — never the report count.** `pending` gets
its own column and the tool prints "The sweep is NOT finished while this is non-zero."

### FM-7. Three more instrument defects from the same release, recorded for the pattern

Not W0's own tools, but the same failure class and a future agent should recognise them:

* **`oracle.py`** called `ROGYAO_comp_0` "distinct, ENCODER-BLIND (total)" at 2.586 Å. The
  molecule is **achiral** (0.423 Å). Cause: the automorphism cap was starved by four methyls on
  the H-explicit graph. Caught by Lane 7's corrected torsion oracle, built for a different purpose
  — i.e. by a *second, independently built* instrument for the same quantity.
* **`adapter_scan`** replays **frozen** OIN strings, so it is structurally blind to encoder fixes.
  It reported `MEGZIH` unfixed *after* it had been fixed. Caught by the atom_count lane
  cross-checking end-to-end instead of trusting the scan.
* **A suite that reported `Ran 623 tests, OK` and was still wrong.** 623 is one short of the 624
  the loader collects; a test file had been edited **while that suite run was in flight in the same
  worktree**, and `discover` imports at collection time. Green with a slightly wrong denominator
  gets quoted; red gets investigated. Rule: do not edit a test file while a suite run is in flight
  in the same worktree.

### What generalises across all of them

In every case the instrument answered a subtly **different question** than the one asked, and the
wrong answer was never implausible — it was always a number you could put in a release note. The
four defences that actually worked:

1. a **trust gate**: reproduce a known result on unmodified code before believing any delta;
2. **loud failure instead of an empty-set zero**;
3. a **second, independently built instrument** for the same quantity;
4. **end-to-end verification** rather than checking that a symptom disappeared.

The same discipline also caught **four** cases where a real-looking defect did not exist at all
(the 7 "wrong-donor" molecules were related by a proper rotation; the trivalent-P gap was already
handled; `ROGYAO` was achiral; three `atom_count` molecules have a CH2 the crystal never located).

### Discrepancies between the sources, flagged rather than smoothed

* **`structural` inflation: 19.3% or 19.6%?** `docs/agentic-notes/v0.4.5/V045_STATUS_2026-07-25.md:878` records
  `1.04% -> 19.3%`; `tools/canonicality_probe.py`'s docstring records "19.6% where the capstone
  sweep reports 1.0%". Same finding, two readings of the same run (or two runs). The conclusion —
  ~19x inflated, cause identified, tool demoted — is unaffected. **Not resolved here.**
* **"1 in 25" vs "1 in 19".** `docs/agentic-notes/v0.4.5/PROMOTION_GATE_v0.4.5.md` §1 and §4 both say key instability
  went from "1 molecule in 5" to "1 in 25". 16/298 = 5.4%, which is **1 in 18.6**;
  `docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md` says "1-in-5 to roughly 1-in-19". The **1-in-19 figure is the
  arithmetically correct one**; "1 in 25" in the promotion gate appears to be an error and should
  not be quoted.
* **Two different baseline levels for byte-stability: 56.1% and 58.1%.** Not a contradiction —
  `RENUMBERING_INSTABILITY` measured 225 molecules (223 encoded) at `main` @ `20044883`;
  `PROMOTION_GATE` measured a 300-molecule seed-42 sample (298 encoded) on
  `trial/v045-integration`. Different samples, different sizes. Quote the sample with the number.
* **The promotion gate document itself carries a self-correction worth preserving:** its first
  publication quoted 62.0% -> 73.5% and key-broken 39 -> 8 (-79%) from **2 of 3 shards**. The third
  shard held harder molecules, so the absolute levels are lower and the key reduction is -73%.
  **The delta held at exactly +11.5 points across both reads** — which is the quantity the decision
  rests on — but the partial absolutes were replaced rather than left standing.
* **`geometry_tag_shift` is quoted at both 0/200 and 0/298.** `V045_STATUS` §"the geometry veto
  PASSED" reports 0/200 (2 of 3 shards); the promotion gate reports the completed 0/298. The
  larger, completed number is the one to cite.

## Where it landed

All paths relative to the repo root. The main checkout's pinned interpreter is
`/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python` (rdkit `==2025.9.3`); **never
`uv sync` in a worktree**. The dataset lives at
`/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset` and is **gitignored**, so a
worktree must pass `--dataset` explicitly.

### `tools/canonicality_probe.py` — PRIMARY canonicality instrument

Measures: is the emitted OIN byte-identical when the same graph is re-presented under proper
rotation and/or atom renumbering, and if not, in which subclass and does the comparison key move
too.

```bash
PYTHONPATH=src .venv/bin/python tools/canonicality_probe.py --n 300 --trials 2 \
    --out tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-canonicality/baseline-main
```

Single-molecule debugging, and the three worked examples:

```bash
PYTHONPATH=src .venv/bin/python tools/canonicality_probe.py \
    --only FEQFIS_comp_0,DUDREA_comp_0,CEBVIR_comp_0 --trials 3 -v
```

Arm-vs-arm (seed is fixed internally, so both arms see the same molecules; shard for parallelism):

```bash
env -u OIN_CANONICAL_BODY -u OIN_CANONICAL_PERCEPTION -u OIN_CANONICAL_SLOTS \
    -u OIN_CANONICAL_ETA_WINDING -u OIN_STABLE_METAL_AC -u OIN_STABLE_STEREO \
  PYTHONPATH=src .venv/bin/python tools/canonicality_probe.py \
    --dataset /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset \
    --n 300 --trials 2 --shard 1:3 --out /path/armOFF
```

### `tools/geometry_tag_shift.py` — perception veto

Measures: does a perception lever move the `[M_XXX]` coordination-geometry tag, and is any move
accompanied by a real coordination-number change.

```bash
PYTHONPATH=src .venv/bin/python tools/geometry_tag_shift.py \
    --lever OIN_STABLE_METAL_AC --n 300 \
    --dataset /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset \
    --out /path/veto
```

### `tools/build_sweep_cohort.py` — frozen cohort

Measures nothing; **guarantees** that two arms see the same molecules and that no two molecules
share a report filename.

```bash
PYTHONPATH=src .venv/bin/python tools/build_sweep_cohort.py \
    --n 5000 --seed 42 \
    --out tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k \
    --overlap-with tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042
```

### `tools/run_sweep.sh` — reproducible sweep runner

```bash
tools/run_sweep.sh tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k \
    tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-5k 6 300
# override the interpreter when running from a worktree:
OIN_SWEEP_PYTHON=/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python \
  tools/run_sweep.sh <cohort> <out> 6 300
```

Wrap long runs so a harness timeout cannot kill them, then poll:

```bash
systemd-run --user -p OOMPolicy=continue -p MemoryMax=14G --unit v045-sweep \
  tools/run_sweep.sh <cohort> <out> 6 300
systemctl --user is-active v045-sweep      # <- NOT the report count
```

### `tools/rebaseline_report.py` — per-molecule transitions

```bash
PYTHONPATH=src .venv/bin/python tools/rebaseline_report.py \
    --sweep tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-rebaseline \
    --cohort spec/handoffs/v0.4.5/rebaseline_cohort.json \
    --out /path/transitions.json
```

Read its footer block before calling anything a regression. Note `elapsed_s` lives inside the
report's `metrics` dict (`tools/test_dataset_roundtrip.py:766`), not at the top level.

### `tools/v045_state.sh` — state computed, not recorded

```bash
tools/v045_state.sh          # full state: branch tips, WIP, levers, jobs, load, next action
tools/v045_state.sh --next   # just the next action
```

### `tools/reencode_ab.py` — retained, DEMOTED

Use **only** for the `key_equal` sub-split, which tracks the sweep (12.04% vs 12.32%). Its
`structural` and `byte_exact` numbers are **not comparable to a sweep's** — see FM-1.

```bash
CAP=tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042
for i in 1 2 3 4 5 6; do
  PYTHONPATH=src .venv/bin/python tools/reencode_ab.py \
    --results-dir $CAP --out /path/armA --shard $i:6 &
done; wait
PYTHONPATH=src .venv/bin/python tools/roundtrip_bucket_report.py --results-dir /path/armA
```

### Reproducing the lane's git history

`git log --oneline main..swimlane/v045-perf` returns **nothing** — the branch is fully merged, so
its tip *is* the merge-base. The W0 tooling commits are reachable as:

```bash
git log --oneline e8b603d5..swimlane/v045-perf   # W0 tools + the perf lane, base = v0.4.5 base
```

## Open questions / for the next agent

1. **The 19.3% `structural` excess is a real, unowned problem.** It is connectivity-perception
   drift: near the covalent-radius + 0.45 Å cutoff, `xyz2AC_obabel` can perceive a genuinely
   different graph for the input vs the generated geometry. `OIN_STABLE_METAL_AC` root-caused one
   mechanism inside it (`DUDREA_comp_0`). The rest is unattributed and it bounds what any
   canonicalization can achieve. Named as a separate perception-hardening problem, not scheduled.
2. **Two of the three renumbering mechanisms remain open at the string level.**
   `FEQFIS_comp_0`'s stereo flip is closed for 8/10 of the sampled stereo class by
   `OIN_STABLE_STEREO`; `CEBVIR_comp_0`'s aromaticity collapse survives all four levers tested and
   has no owner. Suspects, in the order recorded: `AC2BO`/`get_UA_pairs` order-dependence
   (`utils/xyz2mol_local.py:800`, `:542`), `core/chirality.py`'s ordering assumptions, and
   `_align_to_pai`'s index-dependent pivot and `(i+1)**3` Z-sign (`utils/xyz2mol.py:941`, `:971`).
3. **The residual 16 key-broken molecules are characterized but not closed.** Lane 2 showed
   **32/32 of its residual pairs are `same_vcolor_identical`**, so no relabeling at that seam can
   close them (`docs/agentic-notes/v0.4.5/CANONICAL_SLOTS_v0.4.5.md` §7a), plus the 7 wrong-donor molecules Lane 9 owns.
   Do not re-attempt slot relabeling for these without reading §7a first.
4. **`VIMZUN` needs one quiet-host, separate-process, paired run.** It is the single undetermined
   re-baseline regression: 100.0 s at capstone, >300 s now, encode alone ~48 s, sweep ran at load
   30–40 on 12 cores. Everything else flagged as a regression in v0.4.5 is a proven budget
   artifact.
5. **The 5,000-molecule cohort is built and unused.** `cohort-v0.4.5-5k` + manifest, overlapping
   1,319 capstone and 790 v0.4.4-regression molecules — two identical-molecule diffable slices.
   The clean absolute v0.4.5 number that the boron promotion (and anything after it) needs to diff
   against has to come from a run on this cohort at a **stated, matched** timeout budget.
6. **`OIN_EMIT_AXIAL`'s evidence predates canonical perception.** The Y2 cohort A/B (single-axis
   22/22 vs baseline 8/22) was measured before the six levers were promoted. Re-measure under
   canonical perception before acting on it.
7. **Consider generalising `canonicality_probe.py` to a fourth transform.** It currently varies
   presentation only. A conformer transform would overlap the parallel conformer-invariance work,
   and the two instruments should be reconciled rather than duplicated.
8. **`tools/rebaseline_report.py` reads a hardcoded capstone path.** `CAPSTONE` is fixed to
   `tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042/bucket_report.json` — the snapshot known wrong
   four ways. When a clean v0.4.5 baseline exists, that constant should become a flag, and the
   snapshot's four defects should be re-checked against the replacement rather than assumed fixed.
