# Boron cage generation ceiling — v0.4.7 (L4-boronfast)

`OIN_BORON_CAGE` (promoted default-ON in v0.4.6, `docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md`) fixed a genuine
**encoder** ceiling: 34 boron-cluster molecules that used to fail to encode at all now emit a
correct, coordinated, single-bonded deltahedral cage. §10 of that document (added 2026-07-26,
this release's base) measured a 10-molecule sample of the arm nobody had checked -- does the
string that now encodes also *generate* a 3D structure -- and found 0/10. This document extends
that to the full 34-molecule class plus a 14-molecule control group, attributes the cost to a
specific call chain with deterministic counters (not clock-only), ships a fail-fast lever with a
measured safety proof including one confirmed counter-example, and pins the unbounded-budget
mechanism precisely.

**Verdict up front:** the class is real and it is exactly as big as the encoder fix (34 molecules)
plus the population it silently corrected alongside it (14 more, `docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md` §5a)
-- **48 molecules, 1 of which currently generates.** The other 47 either already fail instantly
(7) or burn the whole embed budget finding nothing (40). `OIN_BORON_GEN_FASTFAIL` (default OFF)
converts those 40 from "burn the budget" to "fail in the time it takes to parse the OIN string,"
with zero effect on the 8 that are already fine.

---

## 1. Class size and failure-shape breakdown

Sample: all 34 `OIN_BORON_CAGE` `encode_fail`-class molecules (`docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md` §2)
plus all 14 molecules from that document's §5a ("silently wrong before, correctly cage-encoding
now" -- a *different* population that carries the same motif and is affected the same way, never
measured at the generation level before this). `optimizer=None`, `ensemble_size=1`, 1 worker,
`GEN_CAP=30` (soft, internal) wrapped in a **hard external `timeout -k 10 120`** (120s, 4x the
ask) -- because, per §4 below, the internal budget is not a hard bound. Tools:
`tools/boron_gen_sweep34.sh`, `tools/boron_gen_sweep_14silent.sh`, both new this lane, both
extending the existing `tools/boron_gen_time.py` rather than replacing it.

| outcome | n / 48 | detail |
|---|---:|---|
| **produced a 3D structure** | **1** | `RAWJEG_comp_0` ([Hg_LIN]), 2.5s |
| instant loud failure, geometry code unsupported | 4 | `AVOFIB`, `BEKLUA`, `BEKMIP`, `PEKQUU` -- 0.00-0.01s, pre-existing, unrelated to this lever |
| instant loud failure, cage is an uncoordinated fragment | 3 | `MODZUA`, `RANCIU`, `RULBUV` -- `UncoordinatedFragmentError` in ~0.01s, pre-existing |
| **burned the soft cap (30s), produced nothing** | 36 | 30.4s - 111.1s actual (1.0x - 3.7x the ask) |
| **burned the 120s HARD cap, produced nothing** (never returned inside 4x the ask) | 4 | `XUKRIF`, `HAXJOG`, `COZCEZ`, `OZAREO` |

Class outcomes (produced/not, exception type) are load-independent and hold at any host load;
the `gen_s` figures above are ADVISORY, measured on a box shared with a 5k sweep and sibling
lanes (load 34-58 across this session).

Read plainly: **97.9% (47/48) do not produce a structure. 83.3% (40/48) get there only after
burning a large, uncapped amount of compute.** The 7 already-instant cases were never this
lever's problem; they are listed here to be precise about the denominator.

**Reproduce:**
```bash
V=.venv/bin/python; export PYTHONPATH=$PWD/src
GEN_CAP=30 HARD_CAP=120 tools/boron_gen_sweep34.sh > scratchpad/boron_gen_sweep34.jsonl
GEN_CAP=30 HARD_CAP=120 tools/boron_gen_sweep_14silent.sh > scratchpad/boron_gen_sweep_14silent.jsonl
```

## 2. Where the time actually goes (counters, not clock)

`tools/perf_attribute.py` (pre-existing, v0.4.5 perf lane) monkeypatches the hot call sites to
*count*, not to change behaviour, so the numbers below are exact regardless of host contention.
Run on `XIQKOY_comp_0` (one of the §5a 14, `[Zr_PBP]`, a haptic cage arm) with tight, deterministic
caps:

```bash
V=.venv/bin/python; export PYTHONPATH=$PWD/src
$V tools/perf_attribute.py --dataset $DATASET --molecule XIQKOY_comp_0 --timeout 100 --max-attempts 8
```

Four calls to `embed.get_embedding` (the outer pool-fill loop's per-attempt unit of work) were
captured before the 100s deadline stopped a fifth attempt from starting:

| attempt | wall time | `AllChem.EmbedMolecule` calls (cumulative) | notes |
|---|---:|---:|---|
| #1 (option 0, first touch) | 27.61s | 4 | alt-cache MISS: primes `option 0`'s dummy-atom candidates via a fresh PuLP/CBC solve |
| #2 (option 1, first touch) | 54.91s | 12 (+8) | alt-cache MISS: primes `option 1` |
| #3 (option 2, first touch) | 13.54s | 16 (+4) | alt-cache MISS: primes `option 2` |
| #4 (option 0, second touch) | **91.45s** | 32 (+16, all `rc=-1`) | alt-cache **HIT** -- no new PuLP solve, yet this is the SLOWEST attempt |

Total: `generate_s=188.289` against a 100s ask (88% over), `actual_solve_calls=56`,
`get_valid_molecule_calls=19`, `pulp_cache_hits=7` / `misses=12`, encode itself `0.944s`.

Two independent mechanisms are visible, and **both are unbounded, and neither is capped by
`embed_time_budget`**:

1. **The PuLP/CBC bond-order/charge solve.** `embed.get_alternative_molecule` builds every
   dummy-atom candidate combination for the complex (`itertools.product` over
   `get_dummy_atom_list(cn)` per haptic/multidentate binding site -- `cn=5` returns **two**
   candidates, `[Fe, S]`) and calls `chem.Molecule.get_valid_molecule(method="pulp")` -- a
   full-COMPLEX (metal + every ligand + dummy atoms, not just the cage) bond-order/charge ILP --
   **once per candidate combination**. `_alt_mol_cached` (`embed.py`) memoizes this per `option`
   (0/1/2) for the lifetime of one `generate()` call, so it only pays this ONCE per option --
   but that one payment (attempts #1-#3 above) is a genuinely expensive, unbounded solve: a cage
   vertex has no 2c-2e Lewis structure, so `compute_chg_and_bo`'s multi-objective sequential
   solve (`moSolve`, `generator3d/utils/compute_chg_and_bo_pulp.py`) is fighting a constraint set
   with no clean optimum, which is presumably why it is slow rather than terminating quickly.
   **`pl.LpSolverDefault.actualSolve` is called with no `timeLimit`, no thread count, nothing**
   (`compute_chg_and_bo_pulp.py:254-260`, `moSolve`) -- the CBC solve itself has zero wall-clock
   bound.
2. **The nested embed sweep is NOT memoized and is NOT itself bounded.** Attempt #4 above is an
   alt-cache HIT (no new PuLP work) and still cost 91.45s, entirely from 16 fresh
   `AllChem.EmbedMolecule` calls (~5.7s/call average) that all returned `rc=-1`. This loop --
   `for alternative_ace_mol in alternative_ace_mol_list: for haptic_scale in scales_for_haptic:`
   inside `get_embedding` -- reruns in full on **every** attempt of the outer pool-fill loop
   (`generate_3d_structures`'s `for i in range(max_attempts)`, default `max_attempts =
   max(pool*5, 250)`), cache hit or miss, with no internal time check at all.

So "the cost" is not one thing: it is (a) an uncapped ILP solve, paid up to 3 times per
`generate()` call, plus (b) an uncapped nested RDKit-embed sweep, repeated on **every** attempt
up to `max_attempts` or the between-attempts deadline -- whichever mechanism a given attempt
happens to hit is what makes that attempt slow, and there is no per-attempt cap on either.

## 3. The fail-fast predicate

**Goal:** detect that a cage cannot be assembled before spending the budget on it, converting the
40-molecule "burn the budget" class into a "fail in ~0.01s" class with the same final outcome
(no structure produced) and the same, already-correct message substrings the harness classifies
on ("MetalloGen failed", "failed to generate any conformers").

### 3.1 The naive version, and why it is wrong

The obvious predicate is "does any ligand fragment carry the `OIN_BORON_CAGE` motif" (the same
B-B-B triangle test, `_has_boron_cage` / `boron_cage_vertices`, `utils/xyz2mol.py`). **This is
measured WRONG.** `RAWJEG_comp_0` -- `[Hg_LIN].[BH]12[BH]345C{0}6...` -- carries a genuine cage
(confirmed: `_has_boron_cage` fires) and **produces a structure in 2.5s.** A predicate keyed only
on the motif would have fast-failed a molecule that works today -- exactly the class of mistake
this document exists to avoid making.

Root cause, chased with the same counter/graph tools as §2: `RAWJEG` is `[Hg_LIN]` -- the **only**
CN=2 (linear) molecule anywhere in the 48-molecule sample -- with a **monodentate** cage (one
carbon donor, one binding vector) and no haptic ligand anywhere in the complex. Two more
candidate discriminators were tried and **also refuted** by a counter-example in the same sample,
which is worth recording precisely because each one looked clean until the next check:

- **"Monodentate cages are safe"** -- refuted by `HAXJOG` (`[Rh_TET]`): its cage is ALSO
  monodentate (1 vector, similar atom count to RAWJEG's cage: 13 vs 12) and it still burns the
  hard cap. (`HAXJOG` additionally carries a haptic Cp* ring and a second cage, so its complex is
  much bigger overall -- size/hapticity elsewhere in the complex plausibly explains the
  difference, but "the cage itself is monodentate" alone does not.)
- **"Small total atom count is safe"** -- refuted by `PEKQII` (37 atoms), `SEMTOV` (41), `VEJXOZ`
  (43): all smaller than every failing molecule in the 34-set (smallest there: `RIWKAK`/`RIWKEO`
  at 55) and all three still fail. `RAWJEG` is 25 atoms, well below all of these, so atom count
  alone does not separate the one success from the many failures either.

Geometry code is the one discriminator that survived every counter-example check run in this
sample.

### 3.2 The shipped predicate

`metallogen_adapter._parsed_oin_has_boron_cage(parsed)` fires only when **all** of:

1. A non-metal fragment carries the cage motif (`_has_boron_cage`, reused as-is).
2. That fragment is **coordinated** -- has `>=1` entry in `parsed.vectors`. `MODZUA`/`RANCIU`
   carry the motif on an outer-sphere counterion (0 vectors) and already hit the existing,
   differently-labelled `UncoordinatedFragmentError` just as fast; requiring a vector here means
   this predicate never relabels that already-correct fast failure under a less precise
   exception.
3. `parsed.geo_code` is not in `_BORON_GEN_FASTFAIL_SAFE_GEOMETRIES = {"LIN"}` -- the one
   confirmed-safe case from §3.1.

**Read only from the parsed OIN string's own data** (`parsed.fragments`, re-parsed fresh with
`Chem.MolFromSmiles(..., sanitize=False)`; `parsed.vectors`; `parsed.geo_code`) -- never from
anything the generator constructs (no `metal_complex`, no embedded conformer, no dummy-metal
graph). This is deliberate: a predicate that reads connectivity the *generator* built (e.g. a
metal object's own `GetBonds()`) can certify exactly the defect it exists to catch, because a
mis-attached or detached fragment still owns a bond object in that graph -- a failure mode a
sibling v0.4.7 lane hit this same release with a bond-derived coordination check. Reading the
encoder's own text instead means the predicate answers "did the encoder emit a cage," which is
knowable in about the time it takes to parse a string.

### 3.3 Safety proof

Full 48-molecule validation (`_parsed_oin_has_boron_cage` vs. the measured outcome of §1, every
molecule checked, not a sample of the easy ones):

| | flagged | not flagged |
|---|---:|---:|
| burns the budget (soft or hard cap) | **40 / 40** | 0 |
| produces a structure (`RAWJEG`) | 0 | **1 / 1** |
| already-instant (`NON` geometry / uncoordinated) | 0 | **7 / 7** |

**48/48, zero mismatches.** The predicate catches every currently-slow molecule in the sample and
none of the 8 that are already fine.

End-to-end proof that the lever changes nothing for a working molecule and genuinely fast-fails a
slow one (`OIN_BORON_GEN_FASTFAIL` toggled, everything else identical):

| molecule | lever | `got_mol` | wall time |
|---|---|---|---:|
| `RAWJEG` | OFF | True | 2.56s |
| `RAWJEG` | **ON** | **True** (unchanged) | 3.58s |
| `CAKBEW` | OFF | False | 28.35s |
| `CAKBEW` | **ON** | False (same outcome) | **0.00s** |

**`got_mol=True` alone is not proof of a correct structure** -- a harness that re-encodes through
the generator's own bond graph can score a wrong-graph structure as a pass, which is exactly the
failure mode `docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md` §5a documents for 14 OTHER molecules at the encoder
layer. So `RAWJEG`'s xyz was written to a fresh file and independently re-perceived with a brand
new `XYZToSMILES().convert()` call -- coordinates in, nothing from the generator's own mol object
reused -- and compared against the original OIN:

```
ORIGINAL OIN:   [Hg_LIN].[BH]12[BH]345C{0}6[BH]789C%10[BH]1%11%12[BH]231[BH]423[BH]567[BH]824[BH]9%10%11[BH]%12134.[Cl]{1}
RE-ENCODED OIN: [Hg_LIN].[BH]12[BH]345C{0}6[BH]789C%10[BH]1%11%12[BH]231[BH]423[BH]567[BH]824[BH]9%10%11[BH]%12134.[Cl]{1}
BYTE-IDENTICAL: True
```

Byte-identical. `RAWJEG`'s success is a genuine round trip, not a wrong-graph pass -- the `LIN`
exclusion in §3.2 is not built on a false positive.

Pinned in `tests/unit/test_boron_gen_fastfail.py`: the predicate against every confirmed-failing
geometry class checked (`KIXXOF` TET, `VEJXOZ` TPL, `OZAREO` TPY) and both specificity controls
(`ASUVIV`, `AROTAE` -- boron-rich, no cage motif, must never fire regardless of the lever); an
end-to-end test that the lever raises before embedding a confirmed-slow cage (bounding its own
runtime, since a real run would burn most of a 300s budget if the short-circuit ever regressed);
a monkeypatch proving the predicate is never even *called* with the lever off; and that the
predicate does not fire on a non-cage boron molecule with the lever on.

### 3.4 Ship state

`OIN_BORON_GEN_FASTFAIL`, **default OFF**, registered in `oin/levers.py::_HELD_OFF` with the full
justification (including the refuted-hypothesis account above) inline. Held off because:

- The false-positive check so far covers this 48-molecule sample, not a corpus-wide scan for some
  *other* boron motif or geometry that might trip the same test and currently succeed.
- The `LIN` exclusion rests on `n=1` (`RAWJEG`). `RULBUV` (also `[Hg_LIN]`, two monodentate
  coordinated cages) was checked as a second data point but is **uninformative**, not
  confirmatory: it fails today for an entirely unrelated reason (`UncoordinatedFragmentError` on
  a separate, unattached pyridyl fragment) before the cage question is ever reached -- so it
  neither confirms nor refutes the LIN hypothesis, and the predicate correctly does not touch it
  either way (excluded by geometry, independent of the actual failure cause).

**Promotion gate:** a corpus-wide scan of every molecule carrying >=3 boron, confirming zero of
them produce a structure once flagged (or an explicit accepted-risk note if that scan is not
run), plus at least one more confirmed CN=2 success to move the `LIN` exclusion past `n=1`.

## 4. The unbounded-budget bug

`generate_3d_structures`'s `embed_time_budget` (`generator3d/__init__.py`) is a `deadline =
time.monotonic() + embed_time_budget`, checked **only** at the top of the `for i in
range(max_attempts)` loop (lines ~456 and ~528 for the serial/batched variants) -- i.e. *between*
whole attempts. It cannot interrupt an attempt already in flight, and per §2, a single attempt's
cost is unbounded on two independent axes:

- The PuLP/CBC solve inside a fresh `option`'s `get_alternative_molecule` call has no `timeLimit`
  (`compute_chg_and_bo_pulp.py:254-260`) -- confirmed by reading the solver construction, not
  inferred from timing.
- The nested `alternative_ace_mol_list x scales_for_haptic` `EmbedMolecule` sweep inside
  `get_embedding` has no internal time check either, and reruns in full on every attempt
  regardless of alt-cache status (measured: a cache-HIT attempt cost 91.45s, §2 attempt #4).

Measured consequence: a 100s ask produced a 188.3s run (88% over) where the overrunning attempt
*alone* (91.45s) exceeded the entire remaining budget at the moment it started (t=96.8s, 3.2s of
budget left). The only thing that actually enforces a budget on this class today is an external
SIGKILL watchdog -- confirmed by needing a **hard external `timeout -k 10 120`** wrapper around
an internal `GEN_CAP=30` ask for four molecules (`XUKRIF`, `HAXJOG`, `COZCEZ`, `OZAREO`) that
still had not returned at 4x the requested budget.

**This matters beyond boron:** any per-molecule timing collected without an external SIGKILL
watchdog understates the tail for whatever else can make a single attempt expensive, not just a
boron cage.

### Can a bound be added safely? Not in this lane, and here is why

Two candidate fixes exist, both identified precisely, neither implemented here:

1. **Add a `timeLimit` to the CBC solve** (`compute_chg_and_bo_pulp.py`'s `pl.LpSolverDefault`
   construction). This is the more surgical of the two -- `get_alternative_molecule` already has
   a graceful degradation chain when a solve fails (`pulp` -> `xyz2mol` -> raw single-bonded
   adjacency with charge 0, `embed.py` ~285-296), so a solver that gives up on a time limit would
   very likely be absorbed by existing fallback code, not crash outright. **Not done**: this
   solver call is a shared hot path for every bond-order/charge perception in the generator, not
   just boron. Validating that a time limit does not silently truncate some OTHER molecule
   class's legitimately-slow-but-correct solve requires a full non-boron corpus A/B, which is
   outside this lane's scope and box-time budget.
2. **Thread `deadline` into `get_embedding` and check it inside the nested candidate/scale
   loops.** More invasive (touches three files' function signatures on a byte-identical-guaranteed
   hot path) and has the same validation requirement.

Given the scope of this lane, the fail-fast predicate in §3 is the complete, safe, targeted fix
for the boron-cage class specifically -- with the lever on, a cage molecule never reaches
`get_embedding`/PuLP at all, so the unbounded region is simply never entered for this class. The
general mechanism remains open for whoever next owns generator perf; the exact patch locations
above are recorded so that work does not start from a re-derivation of this section.

## 5. What this does not claim

- **Not a fix for the 40 molecules' accuracy.** Nothing here makes a boron cage assemblable; it
  makes the discovery that it cannot be assembled cost ~0.01s instead of up to the full budget.
  The 1 success (`RAWJEG`) was already working before this lane; nothing here changes it.
- **Not a general perf fix.** §4's unbounded-budget mechanism is a property of the shared
  PuLP/CBC solve and the nested embed loop; both are used far beyond boron chemistry, and fixing
  either safely needs a corpus-wide A/B this lane did not run.
- **Not proven past this sample.** The predicate's specificity is proven against 48 boron-cage
  molecules plus 2 non-cage boron controls (`ASUVIV`, `AROTAE`) -- 50 checks, 0 mismatches -- not
  against every boron-containing molecule in the corpus. See the promotion gate in §3.4.

## 6. Reproduce

```bash
V=.venv/bin/python; export PYTHONPATH=$PWD/src
D=/home/tjmustard/Documents/GitHub/tmCat-tmPhoto/tmCAT-tmPHOTO_xyz_dataset

# SS1: class-outcome sweep, 1 worker, hard external timeout (34-molecule class + 14-control)
GEN_CAP=30 HARD_CAP=120 tools/boron_gen_sweep34.sh
GEN_CAP=30 HARD_CAP=120 tools/boron_gen_sweep_14silent.sh

# SS2: deterministic counter attribution on one haptic-cage molecule
$V tools/perf_attribute.py --dataset $D --molecule XIQKOY_comp_0 --timeout 100 --max-attempts 8

# SS3.3: the lever's own tests (predicate correctness + end-to-end safety proof)
$V -m unittest tests.unit.test_boron_gen_fastfail tests.unit.test_boron_cage tests.unit.test_levers
```

## 7. Commits

| SHA | subject |
|---|---|
| `8e9d8ffd` | feat(gen): fail fast on a boron cage the generator cannot embed (opt-in) |
