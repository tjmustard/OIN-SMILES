# Census C2 — the encoder, judged against itself

**Date:** 2026-09-17 · **Branch:** `research/census` · **Plan:** [`CENSUS_PLAN.md`](CENSUS_PLAN.md) ·
**Before this:** [`CENSUS_C1_generator_vs_input.md`](CENSUS_C1_generator_vs_input.md)
**Instrument:** `tools/census/e_selfconsistency.py` (I2 runner, 11 encodes per molecule) and
`tools/census/veto_probe.py` (the mechanism, with controls).
**Population:** all 5,000 cohort inputs. No generator anywhere. 55,000 encodes, 17.76 CPU-h,
8,017 s wall on 10 workers (the plan forecast 15 CPU-h; the wall is the tail — 1% of base
encodes take > 26 s, the slowest 261 s).
**Outputs (gitignored, freeze at C4):** `results-census/e_selfconsistency.jsonl` (one row per
molecule, every differing string kept), `veto_probe.jsonl`, `veto_probe_chiral.jsonl`,
`veto_probe_renumber.jsonl`, and the `.log` beside each.

## Headline

1. **`E` is deterministic, orientation-invariant and unchanged since v0.4.14** — 4,987/4,987 on
   each of: encode twice, re-written file, random proper rotation, and `base == smiles_1` of the
   v0.4.14 sweep (a different process, a different hash seed, a tree three releases older). Every
   defect below is therefore a property of the encoder, not of the run.
2. **`E` is not canonical for 20.8% of achiral molecules and not injective for 8.3% of chiral
   ones.** The 2×2 the plan asked for, over all 4,987 encodable inputs, joined to the I1 ruler's
   chirality verdict (which never imported `oinsmiles`):

   | ruler says | `E(mirror x) == E(x)` | `E(mirror x) != E(x)` |
   |---|---|---|
   | **achiral** (3,218) | 2,548 | **670 = NON-CANONICAL** |
   | **chiral** (1,769) | **146 = NON-INJECTIVE** | 1,623 (the column fires) |

   C1's sample estimated ~470 and 10.5%; the census says 670 and 8.3%.
3. **The 670 have one named mechanism, and it is the parity veto.** `oin/fold_parity.py::resolve`
   vetoes the donor fold when `S_rot(x) ≠ S_rot(mirror x)` and `S_fold(x) == S_fold(mirror x)`,
   and reads the first inequality as "these are two enantiomers". For an achiral molecule that
   premise is false: the rotation-only labeling differing between the two hands *is* the defect
   and the fold *is* the repair. Worse, the veto fires on **one hand only** — the other hand's raw
   labeling already equals the folded form, so there the fold is inactive and nothing is vetoed
   (`{fold_inactive, vetoed_collapse}` is the modal outcome pair). One molecule, two strings.
   **568 of the 670 are vetoed on a hand, and with the veto off all 568 encode identically**;
   the same asymmetry explains 464 of the 494 slot-level renumber drifters. Against that, the
   veto fires on **36** of the 1,623 ruler-chiral molecules whose hands `E` separates — the
   pairs it actually protects. Sixteen achiral molecules split for every chiral pair kept.
   **The veto has no achirality test. `E` has the 3D to make one; the I1 ruler already does.**
   The 102 residual are the same root cause one level up:
   the slot canonicalization minimizes over *proper* rotations only, which is right for a
   chiral complex and wrong for an achiral one.
4. **18.8% of molecules (937) change their string under a pure atom renumbering or 0.02 Å of
   noise; for 5.0% (250) the renumbering changes the comparison KEY — perception, not labels.**
   The generator returns atoms in its own order and its own geometry, so every one of these is a
   round trip that can fail with a correct structure. They are 12× enriched among the failures
   that produce nothing: **51% of `hard_fail` is perception-fragile vs 4% of passes.**
5. **59% of the encodable failures (666/1,130) carry an encoder flag** — unstable under a proper
   transform, or non-canonical / non-injective under reflection — before any generator evidence
   is consulted. For `key_equal/slot_renumber` it is 246/252.

## What the instrument does, and what a broken one would print

Each input is encoded eleven ways; every way has a known right answer:

| tag | transform | must | result (n = 4,987 encodable) |
|---|---|---|---|
| `again` | same file, second `XYZToSMILES()` | = base | **4,987** |
| `rewrite` | file re-written by this script, untouched | = base | **4,987** (the writer is not a variable) |
| `rot0` | random proper rotation, det = +1 | = base | **4,987** |
| `renum0-2` | random non-identity atom permutation ×3 | = base | 4,238 byte / 4,732 key |
| `noise0-2` | + N(0, 0.02 Å) ×3 | should = base | 4,753 byte / 4,756 key |
| `mirror` | z → −z | = base iff ruler-achiral | 2,694 byte / 3,862 key |

A probe whose variants were never really transformed prints the same "100% stable" a perfect
encoder prints, so every variant is **asserted to be its transform** before it is encoded (0
assertion failures), and three independent checks show the columns fire: (a) the chiral half of
the 2×2 — 1,623/1,769 mirror strings differ; (b) agreement with C1's separately written
`mirror_probe.py` on the 1,044 molecules both encoded — **1,044/1,044**; (c) a must-drift control
for `renumber`: `FEQFIS_comp_0`, a v0.4.5 drifter, drifts 2/3 (`structural`) with the four
v0.4.5 canonicality levers held OFF and 0/3 with them ON.

13 base encodes hit the 300 s kill (= the harness's own encode budget). 12 are the sweep's 12
`encode_fail`; the 13th, `UDEJUI_comp_0`, passed the sweep and timed out here under 10-way
load — the budget is load-sensitive at the margin, as v0.4.14 already recorded. Each *encode*
is under its own alarm, so a hung variant costs one tag, not the molecule.

## The mirror column — non-canonical and non-injective

**670 achiral & differs.** By relation between the two strings: 611 `key_equal/slot_renumber`
(slot labels swapped between symmetry-equivalent donors — exactly C1's description), 46
`facmer_divergent`, 13 `structural`. 505 of the 670 also drift under renumbering; **142 are stable
under every proper transform and differ only under reflection** — a pure parity defect. By the
v0.4.14 honest bucket: 280 passes, 239 `slot_renumber`, 81 `structural`, 52 `hard_fail`, 13
`rdkit_canonical`, 5 `facmer_divergent`. The 280 passes are the at-risk population C1 described:
they passed because the generator happened to return the hand whose labeling `E` was going to
emit anyway.

**146 chiral & same.** 116 are passes — and the round trip *cannot* have verified them, because
the string does not carry the hand. The I1 ruler's real verdict on their generated structures:
23 `ISO + MIRROR` (the generator built the enantiomer and the metric said PASS — the plan's row 7,
*metric false pass*), 19 `ISO + SAME`, 18 `ISO + NO_STEREO` (ligand-level chirality the ruler
does not read), 11 `ISO_CLASH`, and **17 with a different coordination sphere or ligand set
(`SPHERE_DIFF` 9, `LIGAND_DIFF` 8) that still scored `byte_exact`.** The other 30 are failures
(16 `structural`, 10 `hard_fail`, 3 `rdkit_canonical`, 1 `facmer`) whose fault table row must be
5b, not 5: `E` never told `G` which hand to build.

## The mechanism — `tools/census/veto_probe.py`

Three arms, each encoding x and its mirror and recording `fold_parity.last_outcome()`:
`default` (fold ON, veto ON — shipped); `veto_off` (fold alone); `fold_off` (rotation-only).

| population | n | `E(mirror)==E(x)`: default | veto off | fold off | default outcomes {x, mirror} |
|---|---|---|---|---|---|
| **achiral & differs** (the defect) | 670 | 0 | **568** | 0 | 486 `{fold_inactive, vetoed_collapse}` · 82 `{vetoed, vetoed}` · 81 `{inactive, inactive}` · 21 `allowed_separation_survives` |
| achiral & same (control) | 150 | 150 | 150 | 150 | 142 `{inactive, inactive}` · 8 `allowed_preexisting_fold` |
| chiral & differs (control, sample) | 150 | 0 | 0 | 0 | 131 `{inactive, inactive}` · 19 `allowed_separation_survives` |
| **chiral & differs (ALL)** | 1,623 | 0 | **36** | 0 | 1,393 `{inactive, inactive}` · 194 `allowed_separation_survives` · **36 `vetoed_collapse` on a hand** |

Read the controls first. `achiral & same` does not move under any arm: the arms change nothing
that was not already broken. `chiral & differs` is where a broken instrument would show: the
lever is demonstrably reaching the encoder (568 strings flip in the defect row), and on
ruler-chiral molecules the fold does **not** collapse the pair — `allowed_separation_survives`
is the outcome, as designed. The pairs the veto protects are the `vetoed_collapse` outcomes on
chiral molecules: **36 of 1,623** over the whole cohort, and exactly those 36 collapse with it
off. Set that against the 568 achiral molecules it splits.

Two caveats on the 568, both in the ruler's favour but stated. The ruler's blind spots (E/Z,
atropisomers, planar chirality of a substituted η ring, ring cis/trans with symmetric paths)
can only produce a false *achiral*, so 568 is an upper bound on the wrongly-split count; the
fold's own swaps are between donors of one symmetric ligand, which none of those blind spots
makes chiral on its own, so the bound is tight in the cases that matter. And v0.4.11's finding
that the fold "collapses enantiomers in 221 of 393 gains" was measured by a string-level mirror
audit that shares the veto's premise — a string that changes under reflection was read as a
chiral molecule. By the ruler, the fold collapses 36 chiral pairs in this cohort; the rest of
the 221 were achiral molecules being correctly unified. The veto was built to a number that was
mostly the instrument.

Then the defect row: 568 of the 670 have `vetoed_collapse` on at least one hand, and turning
the veto off makes **568/568** encode the same string for both hands. The residual 102 are
`{fold_inactive, fold_inactive}` (81) or `allowed_separation_survives` (21) on both hands and
differ under all three arms — the fold never reaches them. By relation they are 46
`facmer_divergent`, 43 `slot_renumber`, 13 `structural`, and by metal geometry they are
disproportionately low-symmetry polyhedra (`SPY`, `TPY`, `TPL`). Read examples: `AHEVAM_comp_0`
`[Mo_OCT]` swaps `O{4}`/`[OH2]{5}` — two *different* monodentate ligands on enantiotopic
positions, which no within-ligand fold can relate; `BUCDEJ_comp_0` `[Fe_SPY]` swaps the imine N
and the benzimidazole N of one pincer — not the same symmetry class, so again outside the fold;
`CASCIK_comp_0` `[Ni_TPL]` flips an η-winding `B{0<}`/`B{0>}` that is not a stereo element of an
achiral complex, and the key does not fold it either (hence `facmer_divergent`). All three are
relabelings realised by a *reflection* of the whole complex. The canonical slot map
(`oin/canonical_slots.py`) minimizes over `derive_rotation_group`, `det > 0` only; for a
chiral complex that is exactly right, and for an achiral complex the improper operations are
symmetries too and must be in the group. The fold is the special case of this that v0.4.12
built, and the veto is the special case it then blocked on one hand.

**Why one hand only.** `resolve` computes `S_rot = canonicalize(inline, fold OFF)` and
`S_fold = canonicalize(inline, fold ON)`. When the incoming presentation's rotation-only
lex-min already coincides with the fold's lex-min, `S_rot == S_fold`, the function returns early
as `fold_inactive` and emits the folded form. When they differ it builds the mirror, finds
`S_rot ≠ S_rot_m` (true for an achiral molecule too — a reflection is not in the proper rotation
group the rotation-only canonicalization minimizes over) and `S_fold == S_fold_m` (the fold
correctly unifies them), concludes "enantiomer pair", and emits `S_rot`. Which hand takes which
path depends on the input atom order, which is also why 505 of the 670 drift under renumbering.

**The same asymmetry explains the slot-level renumber drift** (`--mode renumber`, same three
arms, the same three permutations as I2):

| population | n | all 3 renumberings == base: default | veto off | fold off | vetoed on some presentation |
|---|---|---|---|---|---|
| **renumber drifts at slot level** | 494 | 0 | **464** | 0 | 471 |
| renumber changes the key (control) | 149 | 0 | 0 | 0 | 25 |
| renumber stable (control) | 150 | 150 | 150 | 120 | 2 |

The key-changing control does not move: no lever setting touches a perception difference.
The stable control moves only under `fold off` (30/150 start drifting), which is the fold doing
what it was built for. The question row: 471/494 are vetoed on some presentation and the veto
alone accounts for 464 of the 494 (94%). Which presentation gets vetoed is decided by the
incoming atom order — the very thing a canonical encoder is supposed to be blind to.

**What this is not.** It is not a claim that the fold or the veto is wrong to exist: 36 real
pairs collapse without the veto, and the stable control shows the fold repairing 30/150 with
it off. It is a claim that the predicate lacks the one input that distinguishes the two cases —
whether x is superimposable on its mirror image under a proper rotation plus a graph
automorphism — and that this is computable from the 3D `E` already holds (the ruler does it in
~1 ms per molecule with no perception at all). With that bit in hand the fold needs no veto on
an achiral molecule, the slot group can include the improper operations for the 102, and the
veto keeps its 36. That is a v0.4.17 lane, not a census chunk; it is not built here.

## The proper-transform columns — perception fragility

| honest bucket | n | renumber changes key | noise changes key | either | share |
|---|---|---|---|---|---|
| byte_exact | 3,857 | 48 | 112 | 158 | 4.1% |
| structural | 484 | 110 | 36 | 138 | 28.5% |
| hard_fail | 266 | 76 | 67 | 136 | **51.1%** |
| key_equal/slot_renumber | 252 | 11 | 8 | 19 | 7.5% |
| key_equal/rdkit_canonical | 113 | 10 | 5 | 15 | 13.3% |
| facmer_divergent | 15 | 0 | 3 | 3 | 20.0% |

"Changes key" means the two presentations of one crystal structure have different comparison
keys — different bond orders, different donor set, different geometry descriptor. Three read
examples: `AFIROW_comp_0` (`hard_fail`) loses aromaticity of a thioester-phenyl under one
permutation (`c1ccc(SC(C)=O)cc1` → `C1[CH][CH]C(=[SH]C(C)=O)[CH][CH]1`); `ALITEU_comp_0`
(`hard_fail`) is perceived with a *different ligand set* (`CC#N`, `FP(F)F` appear); `AHAZOZ_comp_0`
(a pass) moves the η-ring heading marker `{0<}` to a different ring atom. None of these is a
labeling problem; they are `perception_tmc` reading the same coordinates differently depending
on the order the atoms arrive in — which is the order the *generator* chooses on the way back.

The enrichment is the finding: a molecule whose string depends on atom order is 12× more likely
to be one the generator produces nothing for. The plan's row 3 (`G-nothing`) will need a P
sub-split: the string handed to `G` for half of these was one of several `E` can emit.

Noise at 0.02 Å — a tenth of the bond-length error C1 measured in the generator's own output —
changes the key of 234 molecules (4.7%). That is the floor for "the generator was right and the
re-encode read it differently" (C1's 142 generator-correct `structural`), measured from the
input side.

## What this changes in the plan

- **Row 5/5b of the fault table is now decidable per molecule** from `rel.mirror` joined to the
  ruler: 146 chiral molecules are 5b by construction.
- **Row 6 gets its sub-split**: 937 molecules carry an I2 instability flag; for a failure with a
  correct generated structure (C1's 439), whether the input itself is presentation-unstable is
  now a column, not a guess.
- **Row 3 needs a P sub-split** (above).
- **v0.4.17 L1 is re-aimed a second time.** C1 moved it from "build the enantiomer" to "the
  encoder labels achiral molecules by hand"; C2 names the line — `fold_parity.resolve`'s left
  conjunct — and the missing input, a geometric achirality test at encode time. Expected reach
  if fixed with no regression: the 237 vetoed `slot_renumber` failures (4.74 pts) plus a share
  of the 63 `structural` and 36 `hard_fail` in the vetoed 568 — bounded above by C1's 8.78 pts
  of failures with a correct structure. The 219 vetoed passes stop being at risk. Regression
  budget: the 36 chiral pairs, each of which the ruler can name in advance.
- **The 116 blind passes are a metric problem**, not an encoder one, and belong with C1's 380.

## Traps this chunk added or confirmed

- `os._exit(0)` after a piped `print` **loses the whole summary**: the interpreter's stdout
  buffer is never flushed. The first veto probe printed one line; the JSONL had everything.
  Flush before `_exit`, or write the summary to a file first.
- A `ProcessPoolExecutor` parent that has already written its `#DONE` can still be alive in
  `pipe_read` two hours later. Kill by PID from `ps`, never `pkill -f`.
- `fold_parity.last_outcome()` is sticky per thread: a veto-off encode never writes it, so a
  probe that reads it after such an encode reports the *previous* molecule's verdict. Reset
  `_state.outcome` before every encode.
- The encode budget is load-sensitive at the margin: one pass of the sweep timed out under
  10-way load. Do not read a single `TIMEOUT` as a regression without the load recorded.
- `grep -c` exits 1 on zero matches; a waiter that ends with it reports "failed" on success.

## Next (C3)

Parse-back (`smiles_1` → `generation/oin_parser` → contract mol → neutral graph vs input),
collision scan (same string / different molecule; same molecule / different string across
`_comp_n` siblings and the cat/photo twins), perception plausibility flags. Offline, minutes.
Then C4 joins C1–C3 into the attribution table and freezes `results-census/`.
