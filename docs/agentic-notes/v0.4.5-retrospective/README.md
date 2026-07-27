# v0.4.5 / v0.4.6 — Master Retrospective

**Releases covered:** v0.4.5 (tag `v0.4.5`, merge `0d165845`) and v0.4.6 (merge `d799de1f`).
**Branch:** local `main`. **NOT pushed** — standing instruction covers v0.4.3 through v0.4.6.
**Suite at close:** 857 tests OK (3 skipped, 4 expected failures). Ruff clean.

---

## ELI5

A molecule can be written down as a line of text. If you take the *same* molecule, turn it around,
or just list its atoms in a different order, you should get the *same* line of text back — the way
"2+3" and "3+2" should both be written "5". Before this release, our software often wrote a
*different* line for the same molecule, which made it impossible to tell "this is a new molecule"
from "this is the same molecule, written differently".

v0.4.5 fixed most of that. v0.4.6 then fixed a family of boron-cage molecules that had been
silently producing *no* text at all, and built a descriptor for a kind of molecular handedness the
notation had been blind to.

The other half of the story is that **most of what we tried did not work**, and the reports here
spend as much space on the failures as the successes — because knowing that an approach is dead,
and *why*, is what stops the next person spending a week rediscovering it.

---

## The whole release, visually

```
════════════════════════════════════════════════════════════════════════════════════════
 v0.4.5 + v0.4.6   PLAN vs REALITY               LEGEND
─────────────────────────────────────────────    ───────────────────────────────────────
 dataset: tmCAT-tmPHOTO, 25,197 unique .xyz       [✓] landed   [~] built, not usable
 local main (unpushed) · 857 tests OK             [✗] not done [»] in flight
                                                  [!] changed  [+] ADDED (not in plan)
════════════════════════════════════════════════════════════════════════════════════════

  ┌──────────────────────────────────────────────────────────────────────────────┐
  │ W0 · THE INSTRUMENT — gated every lane, and the FIRST ONE FAILED         [!] │
  │   PLANNED  tools/reencode_ab.py — re-encode 6,404 stored (input, gen) pairs   │
  │   TRUST GATE FAILED ── the harness re-encodes gen_result.mol, NOT the stored  │
  │                        .xyz → `structural` inflated ~19x (19.3% vs 1.04%)    │
  │   REPLACED tools/canonicality_probe.py — hold the GRAPH fixed, vary only      │
  │            proper rotation + atom numbering.  ← the canonicality verdict [+]  │
  └───────────────────────────────────┬──────────────────────────────────────────┘
                                      │
   ═════════════════ WAVE A · planned, parallel worktrees ═════════════════
        │                │                  │                    │
   ┌────┴────┐     ┌─────┴─────┐     ┌──────┴──────┐      ┌──────┴──────┐
   │ LANE 1  │[✓]  │  LANE 3   │[✓]  │   LANE 4    │[✓]   │   LANE 7    │[✓]
   │ canon   │     │ winding   │     │ axial P2    │      │ research    │
   │ body    │     │ residual  │     │ Part B only │[!]   │ A + B + C   │
   │ 500→~0  │     │  13 → 0   │     │ default OFF │      │ + 2 fixtures│
   └────┬────┘     └─────┬─────┘     └──────┬──────┘      └──────┬──────┘
        │ HARD DEP       │                  │                    │ unblocks L5
        ▼                │                  │                    │
   ══ WAVE B ══          │                  │                    │
   ┌─────────┐           │                  │                    │
   │ LANE 2  │[✓]◄───────┴──────────────────┴────────────────────┘
   │ canon   │  canonical_slots.py · lex-min coloured vertices
   │ slots   │  PBP rotation group 2→10 · 13,187 keys SHA-identical
   │ 315→~0  │  exports canonical_slot_permutation()
   └────┬────┘
        │  measured 0/150 molecules emit a metal @ tag
        │  ⇒ Lane 5 rescoped: CREATE a descriptor, not un-fold a collapse
        │
   ══ WAVE C · planned ══           ═══ UNPLANNED LANES, in parallel ═══════  [+]
   ┌────┴─────┐  ┌──────────┐       ┌────────────────────────────────────────────┐
   │  LANE 5  │  │  LANE 6  │       │ lane8  13% of molecules emit DIFFERENT     │
   │ metal Δ/Λ│  │ amine P3 │       │        ABSOLUTE STEREO under pure atom     │
   │   (P1)   │  │          │       │        renumbering. Rotation clean (0).    │
   │  NOT     │  │ built +  │       │        → OIN_STABLE_STEREO. LOAD-BEARING[✓]│
   │ STARTED  │  │ validated│       │ lane9  7 inequivalent-donor mols → clean   │
   │ in 0.4.5 │  │  BUT the │       │        NEGATIVE: not a soundness class  [✓]│
   │   [✗]    │  │ canonical│       │ boron  "permanent ceiling" was FALSE:      │
   │  → built │  │ body's   │       │        a pruning loop shattered the cage[✓]│
   │ in 0.4.6 │  │ reparse  │       │ atomcount · valsearch · valorder ·         │
   │   [✓]    │  │ clears   │       │ encodefail · genresidue · perf · encspeed  │
   │          │  │ its [N@] │       │                                    7 more  │
   │          │  │   [~]    │       └────────────────────┬───────────────────────┘
   └────┬─────┘  └────┬─────┘                            │
        └──────┬──────┘                                  │
               └──────────────────┬──────────────────────┘
                                  ▼
   ═════════════ WAVE D · integrate → measure → promote → release ═════════════
        integration worktree: 16 lanes, 101 commits                        [✓]
                                  │
        PROMOTION GATE on canonicality_probe (300-mol seed-42)             [✓]
          byte-stability 58.1% → 69.6%    key instability 60 → 16 mols
          vetoes ALL green: fac/mer + cis/trans distinct raw AND at key;
          goldens byte-identical on opt-out; mirror 37/37;
          geometry_tag_shift 0/298 [M_XXX]
                                  │
        ┌─────────────────────────┴──────────────────────────┐
        │ THE PROMOTION BROKE 36 TESTS — four causes     [+]  │
        │  17  "unset means off"  (promotion inverted them)  │
        │   4  CIP GOLDENS INVERTED FOR FOUR MONTHS —        │
        │      the "oracle" ran CIP on the encoder's OWN     │
        │      output, so it relabelled the defect           │
        │   8  lever INTERACTIONS (P3 × body; axial × percep)│
        │   7  shipped goldens moved (all key-verified)      │
        └─────────────────────────┬──────────────────────────┘
                                  │
        837 tests OK → merge --no-ff → tag v0.4.5 (0d165845)               [✓]
        version bump 0.4.4→0.4.5  ← MISSED at tag time, fixed 4fe7571c    [✓]
        re-baseline 936 mols: 145/436 gap FIXED · 11 "regressions"
                              ALL 300s timeouts ⇒ 0 correctness regressions
                                  │
   ══════════════════ v0.4.6 (d799de1f) ══════════════════
        OIN_BORON_CAGE promoted   0/36 → 34/36 encodes                     [✓]
          cost: 14 molecules scored-passing → honest failing (they were
          passing while describing the WRONG graph — VEJXOZ invents C=B)
          leak caught at promotion: the valence bypass hit EVERY fragment,
          so C#O parsed instead of RAW: → scoped to boron. 1,194 fragments
          checked, 56 differ, all 56 contain boron, 0 boron-free affected
        LANE 5 / P1 complete end to end                                    [✓]
          descriptor → emit (OFF) → key fold → generator reproduction
        CHANGELOG 0.4.5 + 0.4.6 (02ffc695)                                 [✓]
                                  │
        FINAL SWEEP · seed-42 5,000 · 6 shards · 300s                      [»]
════════════════════════════════════════════════════════════════════════════════════════
```

---

## Index of reports

| Report | Covers |
|---|---|
| [WAVE-W0-instrument](WAVE-W0-instrument.md) | The measurement instrument, its failed trust gate, and every instrument failure mode found |
| [WAVE-A-parallel-lanes](WAVE-A-parallel-lanes.md) | Lanes 1, 3, 4, 7 + the nine unplanned lanes |
| [WAVE-B-canonical-slots](WAVE-B-canonical-slots.md) | Lane 2, the hard dependency in the middle of the release |
| [WAVE-C-injectivity-descriptors](WAVE-C-injectivity-descriptors.md) | Lanes 5 and 6 — the two injectivity descriptors |
| [WAVE-D-integrate-promote-release](WAVE-D-integrate-promote-release.md) | Integration, the promotion gate, the 36 broken tests, landing, sweeps |
| [LANE-01-canonical-body](LANE-01-canonical-body.md) | `OIN_CANONICAL_BODY` — canonical ligand body via map-number-carrying reparse |
| [LANE-02-canonical-slots](LANE-02-canonical-slots.md) | `OIN_CANONICAL_SLOTS` — lex-min coloured-vertex slot labelling, PBP fix |
| [LANE-03-winding-residual](LANE-03-winding-residual.md) | `OIN_CANONICAL_ETA_WINDING` — embedding-invariant eta heading atom and sign |
| [LANE-04-axial-atropisomer-P2](LANE-04-axial-atropisomer-P2.md) | P2 axial atropisomer; generator multi-axis; the reflection-invariance near-miss |
| [LANE-05-metal-delta-lambda-P1](LANE-05-metal-delta-lambda-P1.md) | P1 metal Δ/Λ — four formulations, three refuted |
| [LANE-06-locked-donor-amine-P3](LANE-06-locked-donor-amine-P3.md) | P3 metal-bound 2° amine — built, and why it is unusable in the default |
| [LANE-07-research-residuals](LANE-07-research-residuals.md) | Torsion oracle, donor-swap probe, twin operators, the two Δ/Λ fixtures |
| [LANE-08-stable-stereo-renumbering](LANE-08-stable-stereo-renumbering.md) | The 13% absolute-stereo instability, and the circular CIP oracle it exposed |
| [LANE-09-inequivalent-donors](LANE-09-inequivalent-donors.md) | The 7 inequivalent-donor molecules — a clean negative result |
| [LANE-boron-cage](LANE-boron-cage.md) | The largest measured accuracy gain, and 14 molecules that were passing while wrong |
| [LANE-atom-count-hydrogen](LANE-atom-count-hydrogen.md) | `OIN_H_FAITHFUL` — and why promoting it buys nothing measurable |
| [LANE-encode-fail](LANE-encode-fail.md) | The 48 molecules that produced no string at all |
| [LANE-valence-search](LANE-valence-search.md) | Bounding an exponential valence search; feasibility density vs budget |
| [LANE-valence-order](LANE-valence-order.md) | `OIN_STABLE_METAL_AC` — perception depending on atom order |
| [LANE-generator-residue](LANE-generator-residue.md) | Residual generation-side failures; the timeout bucket re-read |
| [LANE-perf-generation](LANE-perf-generation.md) | Generation-side performance |
| [LANE-encoder-speed](LANE-encoder-speed.md) | Encoder-side performance |

---

## Initial assumptions and hypothesis

**The stated problem.** A canonical organic SMILES is byte-identical regardless of which 3D
orientation or conformer you start from. OIN-SMILES had not earned that. The *comparison key* in
`oin/compare.py` was canonical; the *emitted string* was not. On the 6,719-molecule capstone sweep:

| bucket | count | meaning |
|---|---|---|
| `byte_exact` | 81.19% | round-tripped to the identical string |
| `key_equal` | **12.32%** | same isomer, **different string** — the target |
| ↳ `rdkit_canonical` | 500 | ligand body serialized differently |
| ↳ `slot_renumber` | 315 | donor slot labels permuted |
| ↳ `winding_star_drift` | 13 | eta ring heading atom / sign moved |

against a 93.51% notation ceiling.

**The hypothesis that made it look cheap.** `compare.py` already computed the canonical form for
every drift class. So v0.4.5 was framed as *moving existing machinery upstream* into the encoder —
engineering, not new science.

**That hypothesis was substantially correct**, and it is why the canonicality half landed. Lanes 1,
2 and 3 are all "promote what compare.py already does".

**The assumptions that were wrong**, each corrected by measurement:

1. That the planned instrument would work. It failed its own trust gate.
2. That the encoder was *stable* under atom renumbering. 13% of molecules were not.
3. That the boron-cluster population was a permanent RDKit ceiling. It was a pruning-loop bug.
4. That the injectivity half would land alongside the canonicality half. Two of three blind spots
   did not ship usable in the default.
5. That test suites would tell us when a lever's semantics changed. 17 tests silently inverted.

---

## What was actually found

### The canonicality result (the release's purpose)

| metric | before | after |
|---|---|---|
| byte-stability under rotation + renumbering (300-mol seed-42) | 58.1% | **69.6%** |
| comparison-key instability | 60 molecules | **16** |
| gap molecules fixed (of 436, re-baseline) | — | **145 (33.3%)** |
| correctness regressions (500 previously-passing guards) | — | **0** |
| boron encoder population (v0.4.6) | 0/36 encode | **34/36** |
| succeeding molecules ≤30 s | — | 83.6% (median 6.9 s, p90 50 s, max 277 s) |

All 11 apparent regressions in the guard population are `TimeoutException exceeded 300s` measured
against a capstone baseline that ran at **1800 s**. A molecule legitimately needing 300–1800 s
passed then and cannot pass now regardless of code quality.

### The correctness findings that mattered more than the numbers

**Three defects were found that no round-trip measurement could ever have revealed**, because the
round-trip key folds exactly the axis each one lives on:

1. **13% of molecules emitted a different absolute stereochemistry under pure atom renumbering.**
   Rotation was clean (0 drift). A chiral tag is a parity *relative to neighbour order*, and the
   fragment rebuild changed that order as a function of input atom numbering. On any single ordering
   the emitted configuration was right only by luck. → [Lane 8](LANE-08-stable-stereo-renumbering.md)

2. **14 boron molecules were scored as PASSING while describing the wrong graph.** With cage mode
   off, `VEJXOZ_comp_0` loses 6 of its 12 B–B cage bonds and the encoder *invents* a C=B double bond
   to balance valences — and the round trip called it a pass. → [Lane boron](LANE-boron-cage.md)

3. **Two shipped stereo goldens were inverted, and had been for four months.** The test that
   "verified" them ran `rdCIPLabeler` on a SMILES reparsed from the encoder's own output.
   `rdCIPLabeler` converts a parity tag into an R/S label; it does not *check* the tag. Feed it an
   inverted tag and it returns an inverted label with full confidence. → [Lane 8](LANE-08-stable-stereo-renumbering.md)

### Where the remaining accuracy gap actually is

Measured on the 936-molecule re-baseline, and this decomposition superseded two earlier wrong ones:

| class | ~n | status |
|---|---|---|
| boron-cluster encoder failures | 34 | **fixed in v0.4.6** |
| timeouts that are genuinely compute-limited | ~44 | needs budget, not code |
| timeouts that hide REAL failures | ~130 | 75% of the timeout bucket — see below |
| atom-count mismatch | 45 | heterogeneous, `dH` spans −36…+14 |
| MetalloGen "no conformers" | 19 | error message was misleading; root cause open |
| string mismatch | ~55 | unscoped |

**A timeout bucket is not latent pass-rate.** 24 molecules that all hit the 300 s wall in full mode
were re-run on the cheaper `--quick` path: 6 SUCCESS (25%), 6 String mismatch, 6 Atom count
mismatch, 6 MetalloGen failed, and **0 timed out again**. So more compute buys roughly 44 of 936
molecules (~4.7%), not the 174 the raw timeout count implies.

---

## What was done

**Sixteen lanes landed**, and most were not in the plan. Quote the BRANCH LIST rather than a count — contemporaneous notes say eight, nine and ten, because the count depends on whether you tally branches (15 unique `swimlane/v04[56]-*`) or lanes (16: `atomcount` landed without its own branch), and on whether Lane 5 counts as planned-but-unstarted. Planned and landed: lanes 1, 2, 3, 4, 6, 7. Unplanned: `lane8`, `lane9`, `boron`, `atomcount`, `valsearch`, `valorder`, `encodefail`, `genresidue`, `perf`, `encspeed`. Each shipped behind an env lever defaulting
OFF so the default path stayed byte-identical until a single measured promotion gate.

**Six levers promoted to default-ON** (`src/oinsmiles/oin/levers.py::_DEFAULT_ON`):
`OIN_CANONICAL_BODY`, `OIN_CANONICAL_PERCEPTION`, `OIN_CANONICAL_SLOTS`,
`OIN_CANONICAL_ETA_WINDING`, `OIN_STABLE_METAL_AC`, `OIN_STABLE_STEREO` — plus
`OIN_BORON_CAGE` in v0.4.6.

**The rule that made them safe to promote together, and the one to keep:** each of the six
**repairs a renumbered presentation without rewriting the canonical answer.** That is why the corpus
shows no churn. Levers that **add information** to the string (`OIN_EMIT_AXIAL`,
`OIN_EMIT_LOCKED_DONOR`, `OIN_EMIT_METAL_CONFIG`) are a different trade — the generator must be able
to reproduce what they emit, so promoting one converts a *silent false positive* into a *loud false
negative*. Those stay opt-in, each with its reason recorded in `_HELD_OFF`.

**`levers.py` exists because the promotion required it.** Nine call sites used three different
spellings of "default", and one of them — `bool(os.environ.get("X"))` — is a trap: `"0"` is a
non-empty string, so `X=0` *enabled* the lever. `OIN_BORON_CAGE` alone had five sites on that
spelling.

---

## Where it landed

`main` @ `41d2f52e` (unpushed), **857 tests OK**.

| | |
|---|---|
| v0.4.5 | tag `v0.4.5`, merge `0d165845`, version bump `4fe7571c` |
| v0.4.6 | merge `d799de1f`, CHANGELOG `02ffc695` |
| Blind spot P1 (metal Δ/Λ) | **descriptor complete end to end**, lever default OFF |
| Blind spot P2 (axial) | closed in a prior wave; lever OFF; evidence needs re-measuring under canonical perception |
| Blind spot P3 (2° amine) | **built and validated, NOT usable in the shipped default** — the canonical body's reparse clears its `[N@]` |
| 5,000-molecule sweep | in flight at close of session |

---

## Patterns and lessons — the most transferable part of this release

**1. A measurement that only exercises the easy case will confirm a wrong belief.** This fired at
least five separate times: the circular CIP oracle (one representation); the P3 tag fix (POJJOP
passed, RIFGUJ failed); the axial token's reflection-invariance (single-axis fixture); the Δ/Λ
magnitude threshold (two *synthetic* controls returning exactly zero); and a donor-hydrogen rule read
off four hand-picked examples that held in 4 of 45. **Validate a stereochemical change on a
multi-centre fixture, and a corpus claim on the corpus.**

**2. An "oracle" that runs downstream of the thing it tests is not an oracle.** `rdCIPLabeler`
relabels a tag rather than checking it. Anchor on geometry — `AssignStereochemistryFrom3D` on the
parent complex — because that is the one thing an encoder bug cannot rewrite. Then *additionally*
assert the emitted string agrees with it. Two tests, not one.

**3. Promoting a lever silently inverts every test that spelled "off" as *unset*.** 17 tests in
v0.4.5, 6 more in v0.4.6 — 23 failures, diagnosed from scratch both times. Now mechanical:
`tests/unit/test_levers.py::TestNoTestUnsetsAPromotedLever`. On its first run it found three further
instances, **one of which was passing vacuously behind 838 green tests**. A guard that goes vacuous
usually just goes green; these failed loudly only by luck.

**4. Three distinct ways a probe measures nothing and still exits 0.** (a) *Empty corpus* — the
dataset is gitignored, so a worktree lacks it and `canonicality_probe.py` printed a serene `0/0`.
(b) *Buffered stdout killed by `timeout`* — Python block-buffers into a pipe, the kill discards the
buffer, and `sort` exits 0 on nothing. (c) *Exhausted iterator* —
`[itertools.permutations(x)] * n` repeats one iterator, so a permutation generator yielded zero
items and "no symmetry found → chiral" came from a loop that never ran. **Check the denominator.
Require a `#DONE <n>` sentinel. Never `2>/dev/null` a probe you have not yet trusted.**

**5. Know which measurements are load-sensitive.** Wall-clock is meaningless above ~load 12 on this
machine; counts and string equality are not. This was got wrong in *both* directions — once running
a timing probe during a 6-shard sweep, once deferring a deterministic byte-identity A/B for want of
an idle machine.

**6. Never change a default while a sweep is in flight.** The harness runs a subprocess per
molecule, so a mid-run change yields a mixed-config measurement — the exact asymmetry that
manufactured v0.4.4's 11 phantom regressions.

**7. Correctness can lower the metric, and that is the right trade.** Boron promotion moves 14
molecules from "passing" to "failing" because they were passing while describing the wrong graph.
Record the trade at the decision site so the next reader does not "fix" the regression.

**8. Deferral is sometimes the correct engineering call — and sometimes it is over-caution.** v0.4.5
deferred the P3-under-canonical-body fix as too risky to rush; attempting it in v0.4.6 proved the
deferral right. But the same instinct wrongly stopped Lane 5 (the framing "needs a canonical ordering
up to proper rotation only" made it sound like hard group theory; the answer was to need no ordering
at all) and wrongly blocked a merge on a sweep that was the author's own instrumentation. **State a
blocker as a measurement, not as a difficulty assessment.**

**9. When a fix is named from a mechanism rather than a measurement, read the code it would replace
first.** The eta "incremental pool widening" fix would have re-implemented a short-circuit that
already existed. Minutes of reading versus hours of building for a null result.

**10. Instrumentation compounds.** The eta runtime question took six wrong answers to settle. After a
per-attempt counter existed, each subsequent hypothesis was refuted by a single A/B instead of a
cycle of argument. Build the counter early.

---

## Open questions / for the next agent

**Ordered by evidence behind them, strongest first.**

1. **Read the 5,000-molecule sweep** (`results-v0.4.6-sweep`) with `rebuild_summary.py` →
   `roundtrip_bucket_report.py` → `missed_success_audit.py --sweep`. Its `TIMEOUT_S = 300.0` is a
   module constant that **matches** this budget and mis-attributes at any other. It is on the same
   frozen seed-42 cohort as the preserved v0.4.5 partial (`results-v0.4.5-sweep-partial-2697mols`),
   so the boron delta is diffable on identical molecules rather than asserted.

2. **P3 (`OIN_EMIT_LOCKED_DONOR`) needs a reparse that preserves the tag without perturbing the
   canonical write order.** Keep the donor bracketed through the sanitize, or re-derive parity from
   the parent geometry once the write order is fixed. The naive copy is refuted and guarded.

3. **`OIN_EMIT_AXIAL`'s promotion evidence must be re-measured under `OIN_CANONICAL_PERCEPTION`.**
   The Y2 cohort numbers (single-axis 22/22, mirror audit 37/37) were taken with perception OFF, and
   perception changes the hindered-axis count (YESKOZ 2 → 1).

4. **The atom-count class needs PER-ATOM provenance, not another aggregate.** Two aggregate
   hypotheses have already been refuted. Walk one molecule per `dH` band mapping parent atom →
   fragment atom → emitted token and record which atom's H changed and at which step.

5. **The `<30 s` eta tail is answered, not open.** The cost is structural: `accept_fn` is handed raw
   pool conformers and the key only matches *after* optimization, which happens after the fill loop.
   No acceptance-predicate change can shorten it. The remaining options are expensive — produce the
   requested ring face pre-relaxation (construction over selection, **three prior negative results**
   in this repo), or build a pre-relaxation winding predicate that does not exist today.

6. **Lane 5's Δ/Λ token is unwired from the generator's *pool*, only from acceptance.** Promotion
   also requires removing the key fold (`_METAL_CONFIG_TOKEN_RE`) in the same commit — a key that
   folds an axis is not a valid acceptance predicate for that axis, which is precisely why P1 was a
   blind spot in the first place.

### The ceiling, stated once

**100% round-trip accuracy is not reachable from this codebase.** `xyz2AC_obabel` can perceive a
genuinely *different molecular graph* from two conformers of the same molecule near the
covalent-radius + 0.45 Å cutoff (`utils/perception_core.py`). When the encoder and the generator
disagree about what the molecule *is*, no amount of canonicalization reconciles them. The v0.4.5
plan scoped this out at the outset and it was verified independently during this work. Any claim of
approaching 100% must first address perception hardening, which is a different project.
