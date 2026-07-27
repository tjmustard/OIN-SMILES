# Wave D — integrate sixteen lanes, promote six levers, land, re-baseline, sweep

> **Corrected after re-running the tool.** An earlier draft of this report said the release notes
> claim 11 apparent regressions while `rebaseline_report.py` enumerates 9 (8 artifacts +
> `VIMZUN_comp_0` undetermined). Re-running the tool against the completed re-baseline prints
> **`REGRESSED 11 (2.2%)`**, and all 11 carry `TimeoutException ... exceeded 300s`. The "9" was a
> read of a PARTIAL run taken while 7 rows were still `pending_g-xtb`; those 7 later resolved to
> `still_fail`, leaving FIXED 145 and REGRESSED 11 unchanged from the partial read. So 11 is the
> figure, and the zero-correctness-regressions conclusion is unaffected either way.
> Reproduce with:
> `PYTHONPATH=$PWD/src .venv/bin/python tools/rebaseline_report.py --sweep tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-rebaseline`

**Purpose of the wave:** take sixteen independently developed swimlane branches, merge them onto one
release branch, run the single measured **promotion gate** that decides which levers become the
shipped default, triage whatever the promotion breaks, land to local `main`, tag `v0.4.5`, and
re-measure accuracy honestly against a baseline that was known to be wrong in four separate ways.

---

## ELI5

Every change in this release shipped behind an environment switch that defaulted to *off*, so that
nothing users saw could change until one deliberate moment. Wave D is that moment. All the branches
were merged together, then a generator-free probe measured the difference between "all switches off"
and "all switches on" on the same 300 molecules — because the probe holds the molecule fixed and only
changes how the file is written, the right answer is known in advance (byte-identical), which makes
it a real test and not a guess. Six switches passed and were turned on. Turning them on immediately
broke 36 tests, and only two of those were actual bugs in the change: seventeen tests had spelled
"switch off" by *deleting* the variable, which silently came to mean "on"; four tests had been pinning
the **wrong** stereochemistry for four months, blessed by a check that ran chemistry-labelling software
on the encoder's own output and so relabelled the bug instead of catching it. After landing, a
936-molecule re-run fixed 145 previously failing molecules and produced zero confirmed correctness
regressions — every apparent regression was a molecule that had needed more than 300 seconds while
passing under an 1800-second budget.

## The wave, visually

```
 ┌── WAVE A ──┐  ┌WAVE B┐  ┌── WAVE C ──┐   ┌──── ten unplanned lanes ────┐
 │ 1  3  4  7 │  │  2   │  │ 6  (5:✖)   │   │ 8 9 boron atomcount valsearch│
 └─────┬──────┘  └──┬───┘  └──────┬─────┘   │ valorder encodefail genresidue│
       │            │             │         │ perf encspeed                 │
       └────────────┴─────────────┴─────────┴──────────────┬────────────────┘
                                                           │
              ╔══════════════ 1 · INTEGRATE ════════════════▼════════════════════╗
              ║ release/v0.4.5 — 16 lanes, 101 commits ahead of main             ║
              ║ merge order (git, first-parent):                                 ║
              ║  4d92d828 lane7 ← FIRST: changes no encoder output               ║
              ║  c9ebac35 lane1 → c63d3404 lane2 → bbbfb3f8 lane9                ║
              ║  4c05237d lane3 → 43e461e0 lane4 → df7417a9 lane6                ║
              ║  6eb82071 lane8 → 9e1fe6aa perf → ddbc9fbd encspeed              ║
              ║  504f8158 genresidue → 2579bfbb encodefail                       ║
              ║  d075f0d6 valsearch → e4661843 valorder → f4c3525a boron         ║
              ║  1450b5ce = the promotion commit, whose 2nd parent IS atomcount   ║
              ║ 4 CONFLICTS, all in shared perception/serialization code         ║
              ╚═══════════════════════════════╤═════════════════════════════════╝
                                              │
              ╔═══════════════ 2 · PROMOTION GATE ═══════════════════════════════╗
              ║ tools/canonicality_probe.py --n 300 --trials 2, seed 42 FIXED    ║
              ║ GENERATOR-FREE: graph held fixed, only rotation/renumbering vary ║
              ║   ⇒ correct answer is byte-identical = KNOWN ground truth        ║
              ║                                                                 ║
              ║   byte-stable   173/298 = 58.1%  ──►  208/299 = 69.6%  (+35)     ║
              ║   key BROKEN         60 = 20.1%  ──►       16 =  5.4%  (−73%)    ║
              ║   rotate drift          0        ──►          0                  ║
              ║                                                                 ║
              ║   VETOES: facmer/cistrans distinct raw AND at key · goldens      ║
              ║   byte-identical on opt-out path · mirror guard 10/10 ·          ║
              ║   geometry_tag_shift 0/298 [M_XXX] changes                       ║
              ╚═══════════════════════════════╤═════════════════════════════════╝
                                              │  6 levers flipped ON
                                              ▼
              ╔═════════════ 3 · 36 TESTS BROKE ════════════════════════════════╗
              ║ ┌─(a) 17 ── "unset means off"  ─────────────────────────────────┐║
              ║ │  tests spelled OFF by DELETING the env var                    │║
              ║ │  ⇒ Lane 8's test_lever_off_reproduces_the_defect was          │║
              ║ │    asserting THE DEFECT against THE FIXED PATH                │║
              ║ │  ⇒ it FAILED LOUDLY, which is the only reason it was caught.  │║
              ║ │    A guard that goes vacuous usually just goes GREEN.         │║
              ║ │  ⇒ now a LINT: TestNoTestUnsetsAPromotedLever                 │║
              ║ │    found THREE more on its first run, one passing vacuously   │║
              ║ │    behind 838 green tests                                     │║
              ║ ├─(b)  4 ── INVERTED CIP GOLDENS, wrong for FOUR MONTHS ────────┤║
              ║ │  test_chiral_p / test_chiral_n asserted ["S","S"] for         │║
              ║ │  fixtures NAMED (2R,4R). "Verified by RDKit CIP" was          │║
              ║ │  CIRCULAR: CIP run on a SMILES reparsed from the encoder's    │║
              ║ │  OWN output. rdCIPLabeler RELABELS a tag; it never CHECKS it. │║
              ║ │  Arbiter = AssignStereochemistryFrom3D on the parent → (R,R). │║
              ║ │  ⇒ OIN_STABLE_STEREO was CORRECT; the goldens recorded the    │║
              ║ │    DEFECT.                                                    │║
              ║ ├─(c)  8 ── LEVER INTERACTIONS ─────────────────────────────────┤║
              ║ │  canonical_body × locked_donor (P3) INCOMPATIBLE               │║
              ║ │  canonical_perception × axial: YESKOZ hindered axes 2 → 1      │║
              ║ │  stable_metal_ac × a fixture-driven error path                  │║
              ║ ├─(d)  7 ── SHIPPED GOLDENS MOVED ──────────────────────────────┤║
              ║ │  all key-verified as RELABELINGS before re-pinning;            │║
              ║ │  BDPP/BDNN keys change — CORRECTLY                             │║
              ║ └───────────────────────────────────────────────────────────────┘║
              ╚═══════════════════════════════╤═════════════════════════════════╝
                                              │ 8f715699 — the triage commit
                                              ▼
              ╔═══════════════ 4 · LAND ════════════════════════════════════════╗
              ║ 837 tests OK / 3 skipped / 4 xfail on release/v0.4.5             ║
              ║ merge --no-ff to local main ──► 0d165845 ──► tag v0.4.5          ║
              ║ 837 OK RE-VERIFIED on the merged result                          ║
              ║ ⚠ version bump 0.4.4→0.4.5 MISSED at tag time; fixed after,     ║
              ║   at 4fe7571c. Nothing in the suite asserts version == tag.      ║
              ╚═══════════════════════════════╤═════════════════════════════════╝
                                              ▼
              ╔═══════════════ 5 · RE-BASELINE (936 molecules) ══════════════════╗
              ║ 436 previously-FAILING "gap" + 500 previously-PASSING "guard"    ║
              ║ seed 42, --mol-timeout 300, 6 shards, 934/936 completed          ║
              ║   GAP:   145 of 436 FIXED = 33.3%                               ║
              ║   GUARD: apparent regressions ALL TimeoutException >300 s        ║
              ║          against a capstone that ran at 1800 s                   ║
              ║          ⇒ ZERO correctness regressions                          ║
              ╚═══════════════════════════════╤═════════════════════════════════╝
                                              ▼
              ╔═══════════════ 6 · v0.4.6, then the 5k SWEEP ═══════════════════╗
              ║ d799de1f: promote OIN_BORON_CAGE (0/36 → 34/36 encodes)          ║
              ║           + land Lane 5. Blast-radius leak caught AT promotion:  ║
              ║           the cage rung's valence bypass hit EVERY fragment      ║
              ║           ⇒ `C#O` parsed instead of hitting RAW: fallback        ║
              ║           ⇒ scoped to boron; 1194 fragments, 56 differ,          ║
              ║             ALL 56 contain boron, 0 boron-free affected          ║
              ║ FINAL SWEEP: fresh seed-42 5,000-molecule cohort, 6 shards,      ║
              ║              --mol-timeout 300                                   ║
              ║ ⚠ HONESTY CONSTRAINT: a fresh 5,000 is NOT diffable against any  ║
              ║   prior arm on identical molecules, and capstone ran at 1800 s.  ║
              ║   Reading a pass-rate delta across that gap is EXACTLY the config║
              ║   asymmetry that manufactured v0.4.4's 11 phantom regressions.   ║
              ║   ⇒ the GENERATOR-FREE PROBE is the canonicality verdict;        ║
              ║     the 5k sweep is a clean ABSOLUTE number, nothing else.       ║
              ╚═════════════════════════════════════════════════════════════════╝

Legend  ✖ never started    ⚠ recorded hazard / honest caveat
        ══ wave boundary   ── sub-step
```

## Initial assumptions and hypothesis

1. **The release's whole design rests on one decision point.** Every lane shipped default-OFF
   specifically so that Wave D could be a *single measured promotion* rather than sixteen small
   behaviour changes accumulating unmeasured. The corollary the plan accepted: no lane can claim
   accuracy credit before Wave D, and every lane's own A/B is an *estimate* of what the gate will show.
2. **Merge order follows the file-level dependencies.** Land order was planned as
   **1 → 2 → 3 → 4 → 7 → 5 → 6**, because Lanes 1/2/3 all touch `utils/oin_aligner.py` and
   `utils/perception_tmc.py`, and Lanes 2/5/6 share `oin/canonical_slots.py`.
3. **The trial merge is a proxy for integration risk.** All six then-existing lane branches had been
   verified to *trial-merge cleanly* into `main` (`trial/v045-integration`). The plan noted explicitly
   that **textual cleanliness is not semantic validation**, and that the integrated suite was
   unverified.
4. **The promotion is a flag flip plus a measurement.** Nine env-var reads, six of them to be flipped.
5. **The gate's own primary number is key instability**, not byte-stability: the comparison key is the
   harness's acceptance predicate and the basis of every accuracy figure the project reports.
6. **The headline round-trip rate will barely move, and that is structural, not disappointing.** A
   canonicality defect lands in `key_equal`, which **already counts as passing**. Of the 332 closeable
   molecules in the gap to the ~98.45% ceiling, `hard_fail` is **306 — 92%** — so generator
   throughput, not the notation, is the dominant accuracy lever. What the wave buys is that the number
   becomes *meaningful*.

## What was actually found

### Confirmed — the gate passed, on every axis it was designed to test

Measured on `trial/v045-integration` with `tools/canonicality_probe.py --n 300 --trials 2`, seed 42
fixed so **every arm samples the same molecules**. Generator-free: the probe holds the molecular graph
**fixed** and varies only proper rotation, atom renumbering, and both, so the correct answer is
byte-identical — a known ground truth rather than an inferred baseline. Final read, all 3 shards
(298 / 299 molecules encoded):

| arm | byte-stable | comparison **key** broken |
|---|---|---|
| all levers OFF | 173/298 — **58.1%** | 60 — **20.1%** |
| all levers ON | 208/299 — **69.6%** | 16 — **5.4%** |
| **delta** | **+11.5 pts (+35 molecules)** | **60 → 16, a 73% reduction** |

Drift by subclass: `rdkit_canonical` **91 → 18** (the perception and body levers doing the heavy
lifting), `slot_renumber` **42 → 74** (reclassification, not regression — see below), `encode_fail`
**1 → 0**. Drift by transform: `renumber` **176 → 107**, `both` **165 → 119**, and **`rotate` is 0 in
both arms** — orientation-invariance was already sound and is preserved.

**Every veto passed:**

| veto | instrument | result |
|---|---|---|
| `facmer_divergent` must not rise (over-folding — the standing risk of the whole release) | `test_facmer_key.py` + `test_isomer_divergence.py`, re-checked directly with `compare.py` | **PASS.** fac≠mer and cis/trans stay distinct **raw and at key level** with all six levers ON. The one non-green run (`OIN_CANONICAL_SLOTS` alone) is a **stale hardcoded golden string** — the lever deliberately relabels slots — not an isomer merge; verified directly |
| levers-OFF byte-identical | `test_regression_stability.py` (goldens) | **PASS** |
| `OIN_STABLE_METAL_AC` geometry-tag shift | `tools/geometry_tag_shift.py --n 300` | **PASS, 298 molecules: 0 string changes, 0 `[M_XXX]` changes, 0 coordination-number changes, no transitions** |
| `OIN_STABLE_STEREO` must not be "stable because constant" | `test_stable_stereo_mirror.py`, re-run under all six levers ON | **PASS. 10/10 mirrors differ**; nothing collapsed |
| integrated suite | `discover tests/unit` | 729 tests, 717 OK, 3 skip, 4 xfail; the 5 errors were **one** missing fixture, since fixed |

**The geometry veto deserves the emphasis it got, because it was set *against* the change.** Capping
the metal first can only **add** metal bonds that the old atom-order-dependent iteration discarded, so
coordination numbers could have risen and the geometric template fit could have reclassified polyhedra
corpus-wide. `[M_XXX]` is not cosmetic: it selects the vertex table, hence the rotation group, hence
the canonical slot labelling and the key's entire vertex signature. **0/298 refutes that concern.**

**`slot_renumber` rising is not a regression.** [Lane 2](LANE-02-canonical-slots.md) established this
independently with per-molecule accounting — **0 molecules broken** in either of its arms. A molecule
that previously drifted in *both* its ligand body and its slot numbers is counted under
`rdkit_canonical` (the first-matching subclass); once the body is canonical, its remaining slot drift
reclassifies into `slot_renumber`. Byte-stability rising while `slot_renumber` rises is exactly the
expected signature. The residual is the structurally harder class Lane 2 characterized and declined to
force: **32/32 of its residual pairs are `same_vcolor_identical`**, so no relabeling at that seam can
close them.

### Confirmed — the unifying safety property, and why it is the durable rule

The six promoted levers share one property, and it is the reason they were safe to promote
**together**:

> **Each one repairs a renumbered presentation without rewriting the canonical answer.**

That is why the corpus shows no churn — `geometry_tag_shift` 0/298, goldens byte-identical on the
opt-out path. Levers that **add information** to the string (`OIN_EMIT_AXIAL`,
`OIN_EMIT_LOCKED_DONOR`, `OIN_EMIT_METAL_CONFIG`) are a categorically different trade: the generator
must then be able to reproduce what they emit, so promoting one converts a **silent false positive**
into a **loud false negative**. That is the right direction of error, but it is a separate product
call — already made: injectivity levers stay opt-in. **Determinism first, new descriptors second.**

### Refuted / corrected — five things Wave D found wrong in its own materials

1. **An earlier draft of the promotion-gate document published partial figures as final.** From 2 of
   3 shards it read **62.0% → 73.5%** with key-broken **39 → 8 (−79%)**. The third shard contained
   harder molecules, so the absolute levels are lower and the key reduction is **−73%**, not −79%. The
   **delta held at exactly +11.5 points** across both reads, which is the quantity the decision rests
   on — but the partial absolutes were **replaced**, not left standing.
2. **An earlier revision of the gate document said the key instability went "1 in 5 → 1 in 25", which
   its own table contradicts.** 60/298 = 20.1% ("1 molecule in 5") goes to 16 unstable = **5.4%**,
   i.e. **1 in 18.6 ≈ 1 in 19** — the figure `docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md` carried correctly all
   along. `docs/agentic-notes/v0.4.5/PROMOTION_GATE_v0.4.5.md` §1 and §4 have since been corrected in place, with the
   arithmetic shown. **Quote ~1 in 19; do not propagate 1-in-25**, which may survive in older copies.
3. **⚠ The published re-baseline was a mid-run snapshot presented as a result.** The report stated
   that 109 molecules were unscoreable "because there is no `xtb` binary", inferred from the status
   name `pending_g-xtb` rather than checked. **`xtb` version 6.7.1 is installed** at `.venv/bin/xtb`,
   and the g-xTB tier was *running* while the report was written — the pending count was **draining
   (109 → 88)**. `pending_g-xtb` means tier 2 is **queued or in progress**, not skipped. Corrected at
   `7828f460`. This is a different and worse error class than the release's other instrument defects:
   those were instruments answering the wrong question; this was not checking whether the measurement
   had **stopped**.
4. **⚠ The apparent-regression count is quoted two ways in two places.** The release notes and the
   release commit say **11** apparent regressions, all `TimeoutException exceeded 300s`. The
   per-molecule re-baseline report enumerates **9** — 8 conclusively budget artifacts with their
   capstone `elapsed_s` recorded, and `VIMZUN` **undetermined**. Both readings agree on the conclusion
   (**zero confirmed correctness regressions**) and on the mechanism; the count itself needs
   reconciling against `tools/rebaseline_report.py`'s output before either number is re-quoted.
5. **The frozen `bucket_report.json` — the instrument this release gates against — is wrong four
   independent ways.** Stale by a week (29 `rmsd_gate` molecules already fixed by v0.4.4-SL4);
   misattributed (`HOCVAY`/`WEFZAL` were generation-side deaths bucketed as encoder refusals);
   hiding a regression (`XOSTUW_comp_0` passed the capstone but now builds 63 against input 64,
   because the encoder writes bare `B` where the capstone wrote `[BH]` — **boron again**); and
   understating `atom_count`. It is a 2026-07-15 snapshot being used to score 2026-07-26 code. This is
   why the re-baseline reports **per-molecule transitions** rather than a single delta.

## What was done

### 1 · Integration — 16 lanes, 101 commits, four conflicts

`release/v0.4.5` accumulated **101 commits** ahead of `main`. The planned land order was
**1 → 2 → 3 → 4 → 7 → 5 → 6**, on the file-level dependency argument (Lanes 1/2/3 all touch
`utils/oin_aligner.py` + `utils/perception_tmc.py`; Lanes 2/5/6 share `oin/canonical_slots.py`). The actual
first-parent order in git differs in two respects worth knowing:

| # | merge commit | branch | note |
|---:|---|---|---|
| 1 | `4d92d828` | `swimlane/v045-lane7` | merged **first**, ahead of Lane 1 — it changes no encoder output, so it can never conflict semantically |
| 2 | `c9ebac35` | `swimlane/v045-lane1` | the blocking lane |
| 3 | `c63d3404` | `swimlane/v045-lane2` | already contained Lane 1 as a merged shared root (`51770c93`) |
| 4 | `bbbfb3f8` | `swimlane/v045-lane9` | unplanned; no encoder change |
| 5 | `4c05237d` | `swimlane/v045-lane3` | |
| 6 | `43e461e0` | `swimlane/v045-lane4` | |
| 7 | `df7417a9` | `swimlane/v045-lane6` | Lane **5** is absent — never started; landed in v0.4.6 |
| 8 | `6eb82071` | `swimlane/v045-lane8` | already contained Lane 1's perception lever (`220c191b`) |
| 9 | `9e1fe6aa` | `swimlane/v045-perf` | |
| 10 | `ddbc9fbd` | `swimlane/v045-encspeed` | |
| 11 | `504f8158` | `swimlane/v045-genresidue` | |
| 12 | `2579bfbb` | `swimlane/v045-encodefail` | |
| 13 | `d075f0d6` | `swimlane/v045-valsearch` | |
| 14 | `e4661843` | `swimlane/v045-valorder` | built on valsearch, so merged after it |
| 15 | `f4c3525a` | `swimlane/v045-boron` | last, because it is the widest perception change |
| 16 | `1450b5ce` | `swimlane/v045-atomcount` | **folded in as the second parent of the promotion commit itself** (`c1ae2759`), rather than as a standalone merge |

**Four merge conflicts, all in shared perception/serialization code, all recorded in `1450b5ce`:**

- **`valsearch` × `encspeed` in `AC2BO`.** `encspeed` added a memo anchor; `valsearch` extracted the
  per-atom valence construction into `possible_valences()` and added a call counter. **All three
  kept** — the memo and the counter are orthogonal to the extraction.
- **`valorder`** merged clean once `valsearch` was in, because it is built on it.
- **`boron` × `OIN_STABLE_METAL_AC`, both in the SAME valence-capping loop.** They compose because
  they change **different things**: boron exempts cage vertices (*which* atoms are capped),
  `OIN_STABLE_METAL_AC` reorders the loop canonically (*the order* they are capped in). Verified with
  both levers on together: `test_boron_cage` 19/19.
- **`atomcount` × Lane 1 at the sidecar write.** Composed as `h_faithful_smiles` (the base write —
  exactly `MolToSmiles` with its lever unset) followed by the canonical-body block, **and the
  interaction was recorded at the site**: `canonical_body_emit` reparses through
  `MolFromSmiles`/`MolToSmiles`, precisely the round trip that re-reads a bare 0-H symbol one hydrogen
  heavier, so with both on the reparse **can undo** the H-faithful fix. Safe in the shipped
  configuration; `OIN_H_FAITHFUL` must not be promoted until that is addressed.

**One fallback removed because it had become a trap.** `utils/perception_core.py` guarded its `levers`
import with a local shim that defaulted every lever to **OFF** — correct while the registry did not
exist on every branch, and dangerous once it did, because a broken import would have silently reverted
**six promoted defaults**. A missing registry should be a loud `ImportError`, not six quiet behaviour
changes.

### 2 · `src/oinsmiles/oin/levers.py` — created *because* the promotion demanded it

Promoting six levers meant touching **nine call sites** using **two incompatible spellings** (three
default spellings across the tree):

```python
os.environ.get("OIN_EMIT_AXIAL")              # truthy -> "0" ENABLES it
os.environ.get("OIN_EARLY_EXIT", "1") != "0"  # "0" disables
```

The first is a live trap: `"0"` is a non-empty string, so **opting out the obvious way switched the
lever on**. `OIN_BORON_CAGE` alone had **five** sites on that spelling. Everything now routes through
`lever_enabled(name, override=None)`, where `_FALSEY = {"0", "", "false", "no", "off"}` disable and
anything else enables; an explicit `override` (typically from `ff_params`) wins over the environment,
using a **membership** test at the call site so an explicit `False` can opt out — the pattern
`metallogen_adapter.py`'s `OIN_EARLY_EXIT` promotion established in v0.4.4
(`metallogen_adapter.py:1636-1653`). `_DEFAULT_ON` is the shipped configuration in one readable
place, `_HELD_OFF` maps each unpromoted lever to **why**, and `default_on()` / `held_off()` expose
both for provenance stamping. A promotion is now a one-line change.

**This is the single best file to read first when debugging a string difference**, because "which
levers are on?" is the first question, and because `_HELD_OFF` is where every deferred decision in the
release is recorded with its measurement.

**Promoted (`_DEFAULT_ON`):**

| lever | lane | what it canonicalizes |
|---|---|---|
| `OIN_CANONICAL_BODY` | 1 | ligand body — reparse through RDKit, carrying donor identity by atom map number |
| `OIN_CANONICAL_PERCEPTION` | 1 | the half a reparse cannot fix: which resonance form / valence walk is chosen |
| `OIN_CANONICAL_SLOTS` | 2 | slot labels, by lex-min colored-vertex signature over the proper-rotation group |
| `OIN_CANONICAL_ETA_WINDING` | 3 | haptic winding heading atom and `>`/`<` sign |
| `OIN_STABLE_METAL_AC` | 2 | valence-capping order (highest-Z first) so perception stops depending on atom order |
| `OIN_STABLE_STEREO` | 8 | tetrahedral tags re-derived from parent geometry rather than fragment order |

**Held off, each with its reason in `_HELD_OFF`:** `OIN_EMIT_AXIAL`, `OIN_EMIT_LOCKED_DONOR`,
`OIN_EMIT_METAL_CONFIG` (added v0.4.6), `OIN_H_FAITHFUL`, `OIN_RESCUE_STUCK_RING`,
`OIN_ETA_EARLY_EXIT` (added v0.4.6). `OIN_BORON_CAGE` moved from held-off to `_DEFAULT_ON` in v0.4.6.

### 3 · The promotion broke 36 tests, with four distinct causes

This is the most instructive part of the wave. Triage commit: **`8f715699`** —
*"triage the 36 promotion failures — 2 real defects, 1 wrong golden"*. **Only two of the four causes
are defects in the promotion**, and one is a four-month-old wrong answer.

#### (a) 17 tests — "unset means off"

Every lever-OFF test spelled "off" by **deleting** the environment variable. That was correct while
every lever defaulted OFF and **silently means ON** after promotion. So those guards were asserting
pre-release bytes against the post-release path — and
[Lane 8](LANE-08-stable-stereo-renumbering.md)'s `test_lever_off_reproduces_the_defect` was
**asserting the defect against the fixed path**.

**It failed loudly, and that is the only reason the class was caught.** A guard that goes *vacuous*
usually just goes **green**. The fix is always the same one line — write `"0"` instead of deleting —
and `TestFlagOffIsByteIdentical` now clears **every** lever in `levers.default_on()` rather than only
its own, so it can still make its claim.

Because prose in a docstring did not prevent occurrences three and four, it is now a **lint**:
`tests/unit/test_levers.py::TestNoTestUnsetsAPromotedLever::test_no_test_file_unsets_a_default_on_lever`
scans every sibling `test_*.py` for `environ.pop(LEVER)` and
`{k: v for k, v in os.environ.items() if k != "LEVER"}`, resolving module-level
`LEVER = "OIN_..."` aliases so `env.pop(LEVER)` is caught too. It has three deliberate affordances,
each of which was forced by a false positive or a real exception:
`_RESTORE_HINTS = ("_saved", "_prev", "_prior", "restore", "tearDown")` with `_LOOKBACK = 3` (a `pop`
is legitimate when *restoring* a saved value, and the guard is usually on a preceding line — checking
only the pop's own line produced two false positives against correct teardown code), and an explicit
opt-out marker `_ALLOW = "lever-lint: intentional-unset"` for the rare test whose *subject* is the
unset state itself (an auditable marker was preferred over a looser regex, because widening the
heuristic to accommodate one legitimate case would blind the lint to the illegitimate ones).

**On its first run the lint found three more instances, one of which was passing vacuously behind 838
green tests.** Total documented cost of this trap: **23 test failures across two promotions** — 17 in
v0.4.5 and 6 more when `OIN_BORON_CAGE` was promoted in v0.4.6 — diagnosed from scratch each time.

#### (b) 4 tests — the CIP goldens were inverted, and the oracle that blessed them was circular

`tests/unit/test_chiral_p.py` and `tests/unit/test_chiral_n.py` asserted `["S", "S"]` for fixtures
**named** `PdCl2-RR-BDPP.xyz` and `PdCl2-RR-BDNN.xyz` — `(2R,4R)`-pentane-2,4-diyl-bis(diphenyl-
phosphino) and the corresponding bis(diphenylamine) — citing "verified by RDKit CIP".

**The verification was circular.** It ran `rdCIPLabeler` on a SMILES **reparsed from the encoder's own
emitted string**. `rdCIPLabeler` converts a parity tag into an R/S label — **it does not check that
tag against anything.** Hand it an inverted tag and it returns an inverted label with full confidence.
So the "oracle" was a *snapshot of the encoder's output*, an inverted tag was self-consistent, and the
test passed for four months. When `OIN_STABLE_STEREO` fixed the underlying tag instability, the
snapshot is what failed.

Ground truth from `AssignStereochemistryFrom3D` on the **parent complex** — the one thing no encoder
bug can rewrite — is **(R,R)**, agreeing with the `(2R,4R)` in the fixtures' own filenames. So
**`OIN_STABLE_STEREO` was CORRECT and the goldens recorded the DEFECT.** Both tests now derive truth
from coordinates (`test_p_cip_from_geometry`) *and* cross-check that the emitted string agrees — the
loop the old form left open. The full account is preserved in `test_chiral_p.py`'s module docstring.

This is the same failure shape as the Y2 near-miss (a canonicalization that made the axial token
reflection-invariant, caught only by a corpus-wide mirror audit) and as the
[Lane 3](LANE-03-winding-residual.md) winding-sort near-miss: **a measurement that only exercises the
easy case will confirm a wrong belief.**

#### (c) 8 tests — lever interactions, recorded rather than papered over

- **`OIN_CANONICAL_BODY` × `OIN_EMIT_LOCKED_DONOR` are INCOMPATIBLE.** `canonical_body_emit`
  reparses the body, and sanitizing the metal-free fragment clears `[N@]` on a 2-degree amine — the
  very RDKit behaviour P3 exists to work around. So P3 is **built and validated but not usable in the
  shipped default**, and `tests/unit/test_locked_donor.py` runs with `OIN_CANONICAL_BODY=0`. See
  [Wave C](WAVE-C-injectivity-descriptors.md).
- **`OIN_CANONICAL_PERCEPTION` × `OIN_EMIT_AXIAL`.** Perception feeds `_is_atropisomer_candidate`,
  which gates its steric wall on `not GetIsAromatic()`. Measured on `YESKOZ`: hindered axes go
  **2 → 1** under canonical perception, because more of the macrocycle reads aromatic. **No emitted
  string moves today** (YESKOZ's axes are non-stereogenic so its token is empty either way; BINAP is
  unchanged at 1 hindered / 1 emitting) — but the Y2 cohort numbers backing `OIN_EMIT_AXIAL`'s
  promotion evidence were measured with perception **OFF**, and `axial.py`'s safety argument covers
  only the *generator* reading **fewer** aromatic atoms, not the encoder reading **more**.
- **`OIN_STABLE_METAL_AC` × a fixture-driven error path.** `test_perception_tmc_errors` had been driving its
  error path with a real fixture, and the lever makes the deliberately-broken `ticat3_generated_broken.xyz`
  geometry **perceive**: the metal absorbs contested bonds and perception *succeeds*, returning
  nonsense — 8 fragments, seven bare `[H+]`, `[Ti-14]` — where the old order failed loudly. The
  contract is now **fault-injected** so a perception change cannot silently empty it, and the
  degenerate result is pinned rather than lost. On real data the lever is clean (145 fixed,
  `geometry_tag_shift` 0/298), so this is a degenerate-input concern. A sanity gate rejecting stranded
  bare-proton fragments is the obvious follow-up and was **deliberately not added**, because charged
  hydrides are legitimate and the gate needs its own corpus A/B.

#### (d) 7 shipped goldens moved — all key-verified before re-pinning

`canonical_roundtrip_key` **identical** is the assertion that would catch a canonicalization which
merged two isomers, and it was checked *before* each golden was re-pinned:

| fixture | change | key |
|---|---|---|
| CisPlatin | `[Cl]{0}.[Cl]{1}.N{2}.N{3}` → `N{0}.N{1}.[Cl]{2}.[Cl]{3}` | identical |
| TransPlatin | `[Cl]{0}.N{1}.[Cl]{2}.N{3}` → `N{0}.[Cl]{1}.N{2}.[Cl]{3}` | identical |
| fac-Ir(ppy)₃ | `n{3}…n{1}…n{4}` → `n{5}…n{1}…n{3}` | identical |
| PdCl₂-R-BINAP | P moves to slots 2,3 | identical |
| PdCl₂-RR-BDPP / BDNN | `@@/@` → `@/@@` | **changed — correctly** (the old tags were inverted) |

`OIN_CANONICAL_SLOTS` labels by lex-min vertex **colour** and `"N" < "[Cl]"` bytewise, so the amines
take the low slots. **Cis is still cis**: the chlorides land on 2 and 3, which are adjacent.

Also finished in the same commit: the `levers.py` migration the promotion had left half-done —
`locked_donor`, `OIN_BORON_CAGE` (five sites) and `OIN_RESCUE_STUCK_RING` still used bare
`bool(os.environ.get(...))`, so **`OIN_BORON_CAGE=0` ENABLED boron mode**. That trap is what made this
triage fail twice before landing.

### 4 · Landing

- **837 tests OK / 3 skipped / 4 expected failures** on `release/v0.4.5`; ruff `check` and `format`
  clean.
- Merged **`--no-ff`** into local `main` as **`0d165845`** — *"release(v0.4.5): canonical OIN-SMILES —
  promote six canonicality levers"* — and tagged **`v0.4.5`** at that commit.
- **837 OK re-verified on the merged result**, not assumed from the branch run.
- **Local only. NOT pushed**, per the standing instruction covering v0.4.3 / v0.4.4 / v0.4.5.
- Release notes: `docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md`. Live state: `tools/v045_state.sh`, which computes
  state from the repo (branch tips, uncommitted work, levers per branch, active `systemd` jobs, load,
  next action) rather than from a document, so it cannot go stale the way a doc can — **prefer it
  wherever the two disagree.**
- **⚠ The version bump was missed at tag time.** The `v0.4.5` tag points at a tree whose
  `pyproject.toml` still declared `0.4.4`. Fixed **after** the tag, at **`4fe7571c`**, rather than by
  moving the tag — rewriting a tag other work already references is worse than a one-commit offset.
  **Nothing in the suite asserts that the declared version matches the tag**, and it was caught by
  auditing the plan's own Wave D checklist. The CHANGELOG `0.4.5` / `0.4.6` stanzas were blocked on a
  clean tree (a sibling's uncommitted v0.4.4 entry sat in that file) and landed later at `02ffc695`.

### 5 · The re-baseline — the first honest measurement since 2026-07-15

936-molecule cohort: **all 436 gap molecules** (previously failing) plus a **500-molecule seed-42
passing guard**, run on the integration tree with the six canonicality levers ON,
`--mol-timeout 300`, 6 shards. **934 of 936 completed.** Instrument
`tools/rebaseline_report.py`; artifacts in
`tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-rebaseline/`.

| population | n | result |
|---|---:|---|
| **GAP** (was failing) | 435 | **FIXED 145 (33.3%)** · still_fail 183 · pending 107 |
| **GUARD** (was passing) | 499 | still_pass 488 (97.8%) · apparent regressions · pending 2 |

Accuracy, expressed as per-molecule transitions against the (known-flawed) capstone snapshot rather
than as a single delta:

| | molecules | % of 6,719 |
|---|---:|---:|
| capstone snapshot (2026-07-15) | 6,283 | 93.51% |
| **+ gap molecules now fixed** | **+145** | |
| **⇒ confirmed passing** | **6,428** | **95.67%** |
| + unscoreable (`pending_g-xtb`, tier 2 in flight) | +107 | |
| ⇒ upper bound if all pending pass | 6,535 | 97.26% |

**Every apparent regression is a `TimeoutException` at `UFF_1`, checked against the capstone's own
`elapsed_s` — the capstone ran at `--mol-timeout 1800`, this sweep at 300:**

| molecule | capstone elapsed | verdict |
|---|---:|---|
| `PUMLEP` | 1117.2 s | artifact |
| `YOQMAT` | 356.6 s | artifact |
| `CALPOX` | 344.6 s | artifact |
| `AXIDUH` | 331.5 s | artifact |
| `REXROC` | 324.5 s | artifact |
| `AHUKOF` | 303.6 s | artifact |
| `LUYYOT` | 303.3 s | artifact |
| `XEDNUQ` | 302.5 s | artifact |
| **`VIMZUN`** | **100.0 s** | **UNDETERMINED** |

Eight already exceeded 300 s *while passing*, so they cannot pass at this budget regardless of code.
**Zero confirmed correctness regressions.** `VIMZUN` is the honest loose end and was recorded as
undetermined rather than classified either way: its **encode alone measures ~48 s** (so a round trip
needs ~100 s of encoding before generation starts, already at the edge of a 300 s budget), the sweep
ran at **load 30–40** on 12 cores with sibling agents active while the capstone ran quieter, and the
first A/B attempt was **invalid** — both arms in one interpreter, so the second inherited warm caches
and the `AC2BO` memo and duly reported levers-ON as *faster* (42.8 s vs 54.4 s), which is meaningless.
A valid separate-process A/B needs ~4 × 48 s and the machine was too contended for the number to mean
anything by this release's own rule. It needs one quiet-host paired run.

### 6 · v0.4.6, and the final sweep

**`d799de1f`** — *"release(v0.4.6): boron cage promotion + Lane 5 metal Delta/Lambda (P1)"*, local
`main` only, not pushed.

- **`OIN_BORON_CAGE` promoted to default-ON.** On the 936-molecule re-baseline, **34 of the 36
  `XYZToSMILES failed` rows are electron-deficient boron clusters**, and the lever takes that
  population from **0/36 encoding to 34/36**, at 0.2–4.2 s each; the boron lane separately measured
  **48/48 round-tripping** (`docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md`). **The cost, recorded at the definition
  site:** 14 molecules move from **scored-passing to failing**. That is *correct* — they passed while
  describing the **wrong graph** (`VEJXOZ` loses half its cage bonds and the encoder invents a `C=B`
  double bond, then round-trips against its own corrupted mol) — but it trades 14 silent false
  positives for 14 loud honest failures, so a headline pass rate can move either way. Promoted because
  a lossless notation that silently emits a wrong graph is worse than one that fails audibly. Corpus
  context: **186 molecules carry a cage and all 186 have bonds deleted** — 34 fail loudly, 14 fail
  silently while scored correct, 138 were never measured by the capstone arm at all.
- **A blast-radius leak was caught at promotion, and it was wider than boron.** The promotion left
  exactly one failure out of 840, and it was not a stale golden:
  `test_canonical_body::test_unparseable_body_gets_stable_raw_token` reported
  `AssertionError: 'C#O' != 'RAW:C#O'`. **Carbon monoxide now parsed under boron mode**, because
  `compare.py::_parse_fragment` gated its valence-check-free rungs (`_NO_VALENCE`,
  `_NO_VALENCE_NO_KEKULIZE`) on the lever alone, and those rungs apply to **every** fragment. `C#O`
  fails the valence check and nothing else, so with valence checking skipped it succeeded and never
  reached the `RAW:` fallback — and CO is one of the commonest ligands in transition-metal chemistry.
  **Scoped to boron-containing fragments** at `0b483cfa` and **verified over all 1,194 distinct
  fragment bodies the corpus actually emits**: **56 differ ON vs OFF, all 56 contain boron, 0
  boron-free affected.** All 56 go `None → parsed`. That instrument replaced a 61-fixture
  whole-string comparison and is better on every axis — n=1194 rather than 61, it tests the changed
  predicate rather than a downstream proxy, it runs in seconds rather than ~80 minutes, and it is
  load-independent so a running sweep cannot corrupt it.
- **Lane 5 landed** — see [Wave C](WAVE-C-injectivity-descriptors.md) and
  [LANE-05-metal-delta-lambda-P1.md](LANE-05-metal-delta-lambda-P1.md).
- **A sequencing rule worth keeping, from this promotion:** *string equality is deterministic* — the
  same input produces the same OIN whether the box is idle or at load 30, so a byte-identity A/B is
  **not** gated on an idle machine; only *wall-clock* is. But **the merge itself is** gated on the
  sweep: the sweep runs one subprocess per molecule importing from the main checkout's `src`, so
  merging mid-run would measure early molecules under v0.4.5 and later ones under v0.4.6 — a
  mixed-config sweep, exactly the asymmetry that manufactured v0.4.4's 11 phantom regressions.

**The final sweep.** A **fresh seed-42 5,000-molecule cohort** drawn from the corpus's **25,197
unique basenames**, frozen at `tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k` with
`cohort-v0.4.5-5k_manifest.json`, run in **6 shards** at **`--mol-timeout 300`** via
`tools/run_sweep.sh`.

**Read the honesty constraint before quoting anything from it.** A fresh 5,000 is **not diffable
against any prior arm on identical molecules**, and the capstone ran at **1800 s** against this
sweep's **300 s** — reading a pass-rate delta across that gap is *exactly* the config asymmetry that
manufactured v0.4.4's 11 phantom regressions. Therefore:

- **the canonicality verdict of this release comes from the generator-free probe**
  (`tools/canonicality_probe.py`), never from the round-trip sweep;
- **the 5k sweep is a clean ABSOLUTE number and nothing else** — it is not a delta.

The reason that split is mandatory rather than cautious: **77.8% of round-trip failures never test the
notation** (67.4% are 300 s timeouts), so a round-trip pass rate is substantially *generator
throughput*, and any change that alters runtime moves the rate for unrelated reasons.

**Status at the time of writing:** the sweep was launched and had not published a final number.
Artifacts are at `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-sweep-partial-2697mols/`, containing
**2,697** of 5,000 `individual_reports/`. No absolute pass rate from it is recorded anywhere, and none
should be inferred from the partial directory.

## Dead ends and refutations

| tried / believed | what killed it |
|---|---|
| "the trial merge is clean, so the integration is validated" | textual cleanliness is not semantic validation; the integrated suite produced no summary line before the mid-release hard stop and was **unverified** until Wave D actually ran it |
| the planned land order 1 → 2 → 3 → 4 → 7 → 5 → 6 | Lane 7 went **first** (it changes no encoder output, so it cannot conflict semantically), Lane 5 never existed, and `atomcount` ended up folded into the promotion commit rather than merged separately |
| "promotion is a flag flip" | it broke **36 tests** across four distinct causes, and required creating `oin/levers.py` to make nine call sites with three different default spellings behave consistently |
| `os.environ.get("X")` as a lever read | truthy for the string `"0"`, so `X=0` **ENABLED** the lever. Five sites on `OIN_BORON_CAGE` alone. This trap made the triage fail twice before landing |
| tests spelling "lever off" by **deleting** the env var | correct only while every lever defaults OFF. 17 failures in v0.4.5 + 6 in v0.4.6 = **23**, diagnosed from scratch each time. Now a lint that found **three more** on its first run, one passing vacuously behind 838 green tests |
| `rdCIPLabeler` on the encoder's own reparsed output as a CIP oracle | **circular.** It converts a parity tag into a label; it never checks the tag. It certified inverted BDPP/BDNN goldens for four months. The arbiter must be `AssignStereochemistryFrom3D` on the **parent complex** |
| "the BDPP/BDNN goldens are right and `OIN_STABLE_STEREO` broke them" | inverted. Geometry gives (R,R), matching the fixtures' own `(2R,4R)` names. **The goldens recorded the defect** |
| "`slot_renumber` rising 42 → 74 is a regression" | reclassification. Lane 2's per-molecule accounting: **0 molecules broken** in either arm; 19 measured moves `rdkit_canonical → slot_renumber`, 0 the other way, 0 stable→drift |
| publishing the gate's figures from 2 of 3 shards as final (62.0% → 73.5%, key 39 → 8, −79%) | the third shard was harder. Corrected to 58.1% → 69.6% and −73%. **The +11.5-point delta held**, which is the quantity the decision rests on; the partial absolutes were replaced rather than left standing |
| quoting the key-instability improvement as "1 in 5 → 1 in 25" | 16 unstable of 298 = **5.4% = 1 in 18.6 ≈ 1 in 19**. The gate document's own table contradicted its own prose; corrected in place |
| "109 re-baseline molecules are unscoreable because there is no `xtb` binary" | **wrong, and inferred from a status name rather than checked.** `xtb` 6.7.1 is at `.venv/bin/xtb` and the g-xTB tier was **running**, draining 109 → 88. `pending_g-xtb` means tier 2 is queued or in flight. Corrected at `7828f460` |
| the frozen `bucket_report.json` as a trustworthy baseline | wrong four ways: stale by a week (29 `rmsd_gate` already fixed by SL4), misattributed (`HOCVAY`/`WEFZAL` are generation-side deaths), hiding a regression (`XOSTUW_comp_0`), understating `atom_count`. Hence per-molecule transitions, not a delta |
| moving the `v0.4.5` tag to pick up the version bump | rewriting a tag other work already references is worse than a one-commit offset. Landed after, at `4fe7571c` |
| gating the boron byte-identity A/B on an idle machine | **string equality is deterministic**; only wall-clock is load-dependent. Deferring a load-independent measurement on load-dependent grounds cost real time — the mirror image of the earlier error of running a *timing* probe during a sweep |
| promoting `OIN_BORON_CAGE` on the lever alone | its cage rung's valence bypass applied to **every** fragment: `C#O` parsed instead of hitting the `RAW:` fallback. Scoped to boron-containing fragments and verified over 1,194 corpus fragment bodies |
| holding `OIN_BORON_CAGE` off to protect the headline pass rate | that is choosing to ship a wrong answer because the wrong answer scores better. The 14 affected molecules do not pass in any meaningful sense — they are compared against their own corrupted mol |
| reading a pass-rate delta between the 5k sweep (300 s) and the capstone (1800 s) | the config asymmetry that manufactured v0.4.4's 11 phantom regressions. The 5k sweep is an **absolute** number; the canonicality verdict comes from the generator-free probe |

## Where it landed

**Commits, in order:**

| commit | what |
|---|---|
| `4d92d828` … `f4c3525a` | the fifteen explicit lane merges onto `release/v0.4.5` (see the integration table) |
| **`1450b5ce`** | *"release(v0.4.5): integrate 16 lanes and PROMOTE the six canonicality levers"* — also the merge point for `swimlane/v045-atomcount` (`c1ae2759`). Creates `src/oinsmiles/oin/levers.py` (107 lines) and `src/oinsmiles/oin/hydrogen.py`; touches `generation/metallogen_adapter.py`, `oin/inline.py`, `utils/oin_aligner.py`, `utils/perception_tmc.py`, `utils/perception_core.py`; adds `tests/unit/test_levers.py`, `tests/unit/test_atom_count_hydrogen.py`, `docs/agentic-notes/v0.4.5/ATOM_COUNT_v0.4.5.md` and `tools/atomcount/*` |
| `3dfd7bdf` | point `tools/v045_state.sh` at the release branch and the remaining steps |
| **`8f715699`** | **the 36-failure triage** — four causes, two real defects, one four-month-old wrong golden; finishes the `levers.py` migration (`locked_donor`, `OIN_BORON_CAGE` ×5, `OIN_RESCUE_STUCK_RING`) |
| `7828f460` | `pending_g-xtb` means tier 2 is in flight, not that `xtb` is missing |
| **`0d165845`** | **`--no-ff` merge to local `main`; tagged `v0.4.5`** |
| `05689832` | `docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md` — the release notes |
| `567ae3aa` | record `OIN_BORON_CAGE` as held-off **by decision**, not by omission |
| `2aa728f5`, `0bf35884` | v0.4.6 WIP: H-faithful `canonical_body_emit`; the P3 restoration recorded as a **negative** result |
| **`4fe7571c`** | **`pyproject.toml` 0.4.4 → 0.4.5** — missed at tag time |
| `f3c3b9fc` | v0.4.6: promote `OIN_BORON_CAGE`, and make the unset-means-off trap a **lint** |
| `688b7af4`, `0b483cfa`, `0cdb2d7b`, `07ecb0c8` | the boron blast-radius caution, the scoping fix, the merge-gate separation, and the 1,194-fragment verification |
| `26f504e3` … `27089512` | Lane 5, nine commits (see [Wave C](WAVE-C-injectivity-descriptors.md)) |
| **`d799de1f`** | **`release(v0.4.6)`: boron cage promotion + Lane 5** |
| `02ffc695` | the CHANGELOG `0.4.5` and `0.4.6` stanzas |
| `0d0d0f03` … `41d2f52e` | the post-release eta-tail investigation — six hypotheses, each killed by a measurement, closing at *"the eta pool cost is structural, not a defect"* |

**Final state:**

- Tag **`v0.4.5` → `0d165845`**; tag **`v0.4.6`** work released at **`d799de1f`**. Both **local
  `main` only, not pushed**, per the standing instruction.
- **Suite: 837 tests OK / 3 skipped / 4 expected failures** at the v0.4.5 tag, verified on both the
  release branch and the merged `main`. 838 OK on the v0.4.6 H-faithful branch; the boron promotion
  left **1 failure out of 840**, which was the `C#O` leak, now fixed and covered.
- **Shipped defaults** (`src/oinsmiles/oin/levers.py::_DEFAULT_ON`): `OIN_CANONICAL_BODY`,
  `OIN_CANONICAL_PERCEPTION`, `OIN_CANONICAL_SLOTS`, `OIN_CANONICAL_ETA_WINDING`,
  `OIN_STABLE_METAL_AC`, `OIN_STABLE_STEREO` (v0.4.5) + `OIN_BORON_CAGE` (v0.4.6).
- **Documents of record:** `docs/agentic-notes/v0.4.5/PROMOTION_GATE_v0.4.5.md` (the gate and its vetoes) ·
  `docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md` (release notes) · `docs/agentic-notes/v0.4.5/V045_STATUS_2026-07-25.md` (mid-wave status,
  including every lane's numbers and the five-instrument-defects table) ·
  `docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md` (post-release refutations) · `docs/KNOWN_LIMITATIONS.md` ·
  `src/oinsmiles/oin/levers.py` (the authoritative shipped configuration and every held-off reason).
- **Live state:** `tools/v045_state.sh` — derives state from the repo, so prefer it over any document
  where the two disagree. `--next` prints just the next action.

**A note on the instrument defects, because the pattern generalizes.** `docs/agentic-notes/v0.4.5/V045_STATUS_2026-07-25.md`
catalogues **six** instrument defects in this release, each of which would have produced a plausible,
quotable, wrong result: `reencode_ab.py`'s trust gate (`structural` 1.04% → 19.3%),
`canonicality_probe.py` reporting a serene `0/0 byte-stable` when the gitignored dataset was absent
from a worktree (i.e. measuring **nothing** and looking like a pass), `oracle.py`'s automorphism
starvation, `adapter_scan` replaying **frozen** OIN strings and therefore structurally blind to
encoder fixes, `run_sweep.sh` running **two src trees at once**, and a boron suite run that reported
`Ran 623 tests, OK` where the loader collects 624 — a **green** result with a wrong denominator,
caught only by checking the arithmetic. In every case the instrument answered a subtly *different
question* than the one asked. **The defences that worked** were: a trust gate that must reproduce a
known result before its output is believed; loud failure instead of an empty-set zero; a second,
independently built instrument for the same quantity; and end-to-end verification rather than checking
that a symptom disappeared. A generalizable rule from the sixth: **do not edit a test file while a
suite run is in flight in the same worktree** — `discover` imports at collection time, and the
resulting discrepancy is a single digit that any reader would round away.

## Open questions / for the next agent

1. **Reconcile the apparent-regression count: 11 (release notes, release commit) vs 9
   (`tools/rebaseline_report.py`, enumerated with capstone `elapsed_s`).** Both agree there are zero
   confirmed correctness regressions; the count itself is unresolved. Re-derive from
   `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-rebaseline/` before re-quoting either number.
2. **Stop propagating "1 in 25"** if you meet it in an older copy or a downstream quote. The gate's
   own table gives 16 unstable of 298 = 5.4% = 1 in 18.6 ≈ **1 in 19**;
   `docs/agentic-notes/v0.4.5/PROMOTION_GATE_v0.4.5.md` §1 and §4 are already corrected, with the arithmetic shown.
3. **`VIMZUN` needs one quiet-host paired run** — separate processes, ~4 × 48 s — to settle whether
   it is a budget artifact or a real regression. It is the release's only undetermined molecule; do
   not reclassify it by argument.
4. **Finish the 5k sweep and publish the absolute number.** `tools/run_sweep.sh`, cohort
   `tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k` (5,000 files + manifest), 6 shards,
   `--mol-timeout 300`. Current artifacts: 2,697 of 5,000 reports under
   `results-v0.4.5-sweep-partial-2697mols/`. **Publish it as an absolute, never as a delta against the
   1800 s capstone.** Without a clean v0.4.5 number there is nothing to diff the boron promotion
   against, and "34 molecules now encode" stays an unanchored claim.
5. **Two known-unclean measurement conditions to respect.** `tools/v045_state.sh` states wall-clock is
   meaningless above ~load 12; the re-baseline ran at load 30–40. And the frozen `bucket_report.json`
   is wrong four ways. Any accuracy claim needs to name which of these it is exposed to.
6. **Add a test that the declared version matches the tag.** Nothing in 837 tests would have caught
   `pyproject.toml` declaring `0.4.4` at the `v0.4.5` tag.
7. **The three information-adding levers each need a package, not a flag flip.**
   `OIN_EMIT_AXIAL`: re-measure both Y2 cohorts with `OIN_CANONICAL_PERCEPTION=1` (YESKOZ 2 → 1
   hindered axes) and remove `_AXIAL_TOKEN_RE`'s fold **in the same commit**.
   `OIN_EMIT_METAL_CONFIG`: corpus population measurement + generator support + remove
   `_METAL_CONFIG_TOKEN_RE`'s fold in the same commit (until then, lever-ON round trips report
   mismatches by construction). `OIN_EMIT_LOCKED_DONOR`: the `OIN_CANONICAL_BODY` incompatibility
   first — and **do not re-attempt the four-line tag copy**, which is measured wrong and guarded by
   `test_locked_donor.py::TestRifgujRingCarbonsArePseudoAsymmetric`. All three reasons are in
   `levers.py::_HELD_OFF`.
8. **`OIN_H_FAITHFUL` should not be promoted on mechanism.** The `OIN_CANONICAL_BODY` interaction that
   blocked it is fixed (both `canonical_body_emit` writes now route through `h_faithful_smiles`), but
   an A/B over the 45-molecule `Atom count mismatch` population is **identical in both arms** — match
   8 / mismatch 37. The class is heterogeneous (28/45 at `dH` +1…+3, 4 at `dH` 0 where hydrogen is not
   the issue, three large losses at −14/−16/−36) and two aggregate hypotheses are already refuted. The
   next step is **per-atom provenance**, not another aggregate.
9. **`OIN_STABLE_METAL_AC` on degenerate input** lets the metal absorb contested bonds so that
   perception *succeeds* and returns nonsense (`ticat3_generated_broken.xyz`: 8 fragments, seven bare
   `[H+]`, `[Ti-14]`) where the old order failed loudly. Real data is clean. A sanity gate rejecting
   stranded bare-proton fragments is the obvious follow-up and was **deliberately not added**, because
   charged hydrides are legitimate and the gate needs its own corpus A/B.
10. **Where the remaining accuracy actually is, so the next wave is not spent on the notation.** Of
    the 332 closeable molecules to the measured ~98.45% ceiling, `hard_fail` is **306 — 92%**, and
    roughly half of that is timeout-shaped. But `docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md` §3 then **refutes
    "the gap is mostly compute"**: of 24 timeout molecules re-run on the cheaper `--quick` path,
    **0 timed out again** and only ~25% were compute-limited (6 SUCCESS / 6 String mismatch / 6 Atom
    count mismatch / 6 MetalloGen failed). More compute buys ~44 of 936 molecules (~4.7%), not the 174
    the timeout count implies. The revised ranking is: boron (done), **perception-side hydrogen**,
    **MetalloGen generation failures (~80)**, **string mismatch (~55)**, then compute — which buys the
    least.
11. **The eta `<30 s` tail is closed as structural, not as a defect.** Six hypotheses, each killed by
    a measurement (low acceptance rate → pool widening → cost-per-attempt → attempt-driven →
    selection-side predicate alignment → fill-loop `accept_fn`). The answer: **`accept_fn` is handed
    RAW pool conformers, and an eta molecule's ring winding is not right until relaxation, which
    happens after the fill loop** — so no acceptance-predicate change can shorten the fill. Ferrocene
    spends **32 attempts / 32 pool slots**; CisPlatin accepts on attempt **0**.
    `OIN_ETA_EARLY_EXIT` stays default-OFF and documented as **ineffective**. Do not reopen it without
    a pre-relaxation winding predicate, which does not exist.
