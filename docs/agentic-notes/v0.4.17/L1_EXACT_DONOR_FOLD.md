# v0.4.17 L1 — the donor fold swaps look-alikes, not symmetries

**Branch** `research/v0417-encoder-canonicality` (worktree `../oin-v0417`, off `research/census`).
**Lever** `OIN_EXACT_DONOR_FOLD`, default **OFF**. **Status:** built, unit-tested, gated offline
(exact) and live (full cohort). Not promoted; no sweep owed yet.

## Headline

The census (C2) located the `E1_NONCANONICAL` defect — 266 failures, 5.32 pts — in
`fold_parity.resolve`'s left conjunct and prescribed a geometric achirality test. The defect is
one level further down, and fixing it there needs **no achirality test, no veto and no mirror
encode**:

> `_donor_swap_permutations` permutes every symmetry class of a fragment **independently**. A
> ligand automorphism moves all of its classes **at once**. A one-class swap is therefore not a
> symmetry of the ligand; the labeling it produces describes a *different arrangement*.

For a linear tetradentate `t1–i1–i2–t2` wrapped cis-α, exchanging the terminal pair alone yields
the mirror image. That single over-fold is why the fold collapsed enantiomers in v0.4.11, why
v0.4.12 needed a veto, and why that veto — asked to tell an enantiomer pair from an achiral
molecule using strings alone — fires on one hand of 568 achiral molecules.

Fold over **true fragment automorphisms** instead
(`canonical_slots._donor_automorphism_permutations`) and every candidate labeling describes the
molecule in hand. Then a chiral complex and its mirror have disjoint candidate sets (nothing to
veto), and an achiral complex's mirror labeling is `g·L·a` for a rotation `g` and an automorphism
`a` — already inside its own candidate set (one string, no test).

| | shipped (bucket fold + veto) | exact fold, no veto |
|---|---|---|
| the 568 achiral pairs a hand is vetoed on: `E(mirror x) == E(x)` | 0 | **550** |
| `E1_NONCANONICAL` failures whose pair is unified | 0 / 266 | **252 / 266** (5.04 pts upper bound) |
| the 36 ruler-chiral pairs the veto protects: still split | 36 | **19 by construction**, 17 unified — see § "The 17" |
| extra encodes on a fold-active molecule | 2 (self-check + mirror) | **0** |

## The instrument, and what a broken one would print

`tools/v0417/autofold_audit.py`. The slot post-pass is a pure function of the inline string and
its candidate set is an orbit, so applying a fold to the census's stored rotation-only strings
(`veto_probe*.jsonl`, arm `fold_off`) gives byte-for-byte what the encoder emits. It calls the
`src/` implementation with the lever forced — the first draft carried its own copy of the
automorphism search and was rewired before any number here was taken.

**Positive control, read first:** the SHIPPED bucket fold applied to `fold_off` must reproduce the
stored `veto_off` string. **5,792 / 5,792** (x and mirror, all four populations). A path that did
not reach the encoder's code, or a stale string, would print less than that.

**The controls that must not move:** `achiral & same` 150/150 unified under both folds;
`chiral & differs` sample 0/150 unified under both. A fold that unified everything, or nothing,
would fail one of them.

## Results (offline, exact)

`E(mirror x) == E(x)`, bucket fold × exact fold:

| population (C2's, by the ruler) | n | both | **bucket only** | exact only | neither |
|---|---:|---:|---:|---:|---:|
| achiral & differs (the defect) | 670 | **550** | **18** | 0 | 102 |
| achiral & same (control) | 150 | 150 | 0 | 0 | 0 |
| chiral & differs (control sample) | 150 | 0 | 0 | 0 | 150 |
| chiral & differs (ALL) | 1,623 | 17 | 19 | 0 | 1,587 |
| — of which the 36 protected pairs | 36 | 17 | 19 | 0 | 0 |

"exact only" is 0 everywhere: the exact candidate set is a subset of the bucket one, as it must be.

By census fault, the 550: `E1_NONCANONICAL` **252**, `NONE` (passes at risk) 195, `G_CONSTRUCTION`
55, `G_NOTHING` 15, `E1_GRAPH` 13, `E1_HCOUNT` 12, other 8. The 266 `E1_NONCANONICAL` failures
split **252 unified · 5 kept split · 9 in the 102**.

## The 18 — a correction to the census, not a miss

Eighteen pairs the ruler calls achiral, the bucket fold unifies, and the exact fold keeps apart.
Read one by one they are **chiral by how a ligand wraps**:

- six Hf/Zr/Y bis(phenolate)–bis(ether) `[OOOO]` catalysts (`VOLGOV`, `VOLHAI`, `VOLHIQ`,
  `VOLPUK`, `YUXREP`, `ZAYWOO`): phenolates trans, ethers cis — the textbook C2-symmetric cis-α
  tetradentate, chiral and used as such. x and mirror differ by the **ether pair alone**;
- helical pentadentates (`GEFGUW` PNNNP, `VIMQAJ` SNNNS, `ULIDIE`), a tetrahedral bis(N,O)
  (`AMAROX`), three boron-capped clathrochelates labelled `OCT` (`FOSNEI`, `HUTCIK`, `HUTCOQ`).

The ruler's sphere test (`neutral_graph.sphere_fit`) matches donors by **symmetry class and
ligand id**, not by automorphism, so inside one ligand it allows exactly the one-class swap the
bucket fold allows. **The ruler and the bucket fold share a blind spot**, which is why C2's
"568/568 agree with the veto off" could not see these. It is not in the ruler's stated list and
should be added: *chirality carried by chelate connectivity within one ligand reads as achiral.*

Consequences for the census table: five `E1_NONCANONICAL` rows (`FEDYUI`, `JOWBOP`, `TIZLAQ`,
`ULIDIE`, `VIMQAJ`) are chiral molecules whose failure is something else (most plausibly
handedness), and C2's 670 "non-canonical" is at most 652. The table is frozen; the correction
is recorded here, not applied.

## The 17 — the open question, stated rather than buried

The owner's gate was "the 36 protected pairs must stay split". Nineteen do, by construction. The
other seventeen are unified, and I am not going to claim that is free. What they are:

- **the two strings name the same molecule by any reading.** `PEBVEZ`:
  `CN{0}(C)CCN{1}(C)C` vs `CN{1}(C)CCN{0}(C)C` — tmeda written forwards and backwards. `JOCGOZ`:
  the two pyrazolyl arms of one tripod. `CUKPUV`: the two amido N of a dimethyl-cyclen. The only
  difference is which of two **automorphic** atoms carries which integer, and that names no
  property of a molecule;
- **the chirality is real but lives where the string does not look**: a metal-bound N–H or C–H
  stereocentre whose tag is cleared (`JOCGOZ`, `LISNUX`, `GIMKIY`, `PEBVEZ`, `SIDNOK/SIDNUQ` — the
  Y1 audit's P3 blind spot), an η² alkene face (`FOJZOU`), or a conformational twist of a sphere
  that is achiral as a configuration (`ODIHOY` bis-nacnac D2d, `YEGDIA` neocuproine, `CUKPUV`);
- **the split was not stable.** With the SHIPPED encoder **13 of the 17 already change string
  under three random renumberings** (C2's own `e_selfconsistency.jsonl`). A label that atom order
  flips was not protecting anything.

So these move census "chiral & same" from 146 to 163. The honest reading is that they were always
in that class and the veto was hiding it behind presentation noise; the fix is to encode the
centre (Y1 P3), not to keep the label. **This is an owner decision** and is listed as one below.

## Side finding — the shipped string can describe a different arrangement

If the bucket fold's lex-min lands on a non-automorphic candidate and the mirror test does not
veto it, the encoder emits a labeling **outside the input's own orbit**: a string for an isomer
the input is not (often a geometrically impossible one). Among the 2,443 audited molecules the
shipped `E(x)` is outside x's orbit for **33** — and they are `hard_fail` **13 times (39% vs a
6.5% base rate)**, `G_NOTHING` 12, `G_CONSTRUCTION` 7. The generator was handed a complex that
cannot be built. The exact fold cannot emit such a string. The census put these under **G**; the
fault is **E1**. The un-audited remainder (~2,550, mostly achiral & same) is covered by the live
run below.

## Live, full cohort

`tools/census/e_selfconsistency.py` re-run unchanged with `OIN_EXACT_DONOR_FOLD=1` — same seeds,
same eleven presentations per molecule as C2, so every column is comparable line by line.
Output `results-v0.4.17-exactfold/e_selfconsistency_exact.jsonl`.

`#DONE 5000`, 6,909 s wall, stderr empty, `INSTRUMENT` 0, `MISSING` 0. Report:
`tools/v0417/exact_fold_live_report.py` (reads both files, encodes nothing).

**Read this row first — it is what licenses the offline audit.** For all 2,443 audited
molecules the string the live encoder emitted equals the audit's prediction, for x **and** for
its mirror: **2,443/2,443 and 2,443/2,443**. Molecules where the live string equals the SHIPPED
one although the audit predicted a change — what a lever that never reached the encoder would
print: **0**.

| | shipped | exact fold | |
|---|---:|---:|---|
| base encoded (of 5,000) | 4,987 | 4,989 | two fewer 300 s timeouts — no veto encodes |
| again / rewrite / rotate not byte-stable | 0 | 0 | determinism floor unchanged |
| **achiral & `E(mirror) ≠ E(x)`** — NON-CANONICAL | **670** | **120** | −550, the audit's number exactly |
| chiral & `E(mirror) == E(x)` — NON-INJECTIVE | 146 | 163 | +17, the pairs of § "The 17" |
| **renumber drift at slot level** (3 renumberings) | **494** | **33** | −461; C2 attributed 464 to the veto |
| renumber drift at KEY level (perception) | 255 | 255 | membership changed on **0** molecules |
| `E1_NONCANONICAL` failures now mirror-stable | 0 / 266 | **252 / 266** | |
| CPU-h for the 55,000 encodes | 17.76 | 13.29 | ⚠ not like-for-like load (10 vs 9 workers + a unit suite); direction only |

The key-level row is the control that matters: a lever that touched perception would move it,
and it moved on no molecule. The 120 that remain are the 102 (§ "Not done") plus the 18.

**Movers — the generator A/B population, derived from coordinates.** `E(x)` changes for
**392 / 4,987**: `key_equal` 166, `byte_exact` 132, `structural` 54, `hard_fail` 40; by census
fault `E1_NONCANONICAL` 158, `NONE` 117, `G_CONSTRUCTION` 48, `G_NOTHING` 25, other 44. List:
`results-v0.4.17-exactfold/exact_fold_movers.txt`. This is the INPUT side only. A molecule whose
input string does not move can still change verdict because the re-encode of its generated
structure does — 252 − 158 = at least 94 `E1_NONCANONICAL` failures are of that kind — so the A/B
population is this list ∪ the generated-side movers (`tools/lever_string_movers.py`).

**What this run does not say.** It is an encoder-only run. It cannot show a generator-side loss
among the 132 passing movers, which now hand the generator a differently-labelled string
(v0.4.14 measured 7 such losses in 182). That is L1b.

## What was built

| file | change |
|---|---|
| `src/oinsmiles/oin/canonical_slots.py` | `_automorphism_query`, `_donor_automorphism_permutations`; `canonical_slot_relabeling` selects it under the lever |
| `src/oinsmiles/oin/fold_parity.py` | `resolve` returns the exact fold with outcome `exact_fold`, no veto, no mirror |
| `src/oinsmiles/oin/levers.py` | `_HELD_OFF["OIN_EXACT_DONOR_FOLD"]` with the evidence and the open question |
| `tests/unit/test_exact_donor_fold.py` | 9 tests on six real cohort pairs, incl. the defect pinned (bucket collapses `VOLGOV`, exact does not) |
| `tools/v0417/autofold_audit.py` | the offline gate |

Two implementation facts worth keeping:

- **Radicals.** `CanonicalRankAtoms` ignores radical electrons; the substructure matcher does
  not. `LAMTAX`'s corrole has skeleton classes A~B, C~D and *zero* non-trivial self-matches until
  the `[CH][CH]` radicals are zeroed. Five pairs moved on that one line.
- **Pendant explosion.** Eight t-Bu groups are 6⁸ automorphisms acting identically on the donors.
  Non-donor leaves are stripped to a fixpoint and each surviving atom's full-graph symmetry class
  is pinned as its isotope, so the stripped graph cannot gain a symmetry. Over the cap the
  fragment contributes the identity — a missed fold, the safe direction.

## Not done, on purpose

- **The 102.** Not one mechanism: 63 differ in slot digits only (42 are `SPY` — a cis-dioxo /
  planar-tridentate sphere that is Cs in reality and formally chiral once one oxo is called
  "apical"; a template artifact), 15 fragment order, 11 η winding char, 13 carry inverted `@`
  tags (chiral by the string's own account). The first three groups are where C2's achirality
  bit is genuinely needed — and it must be **automorphism-based**, or it will collapse the 18.
- **Generator A/B.** The string is the generator's input (v0.4.14's hole); an offline audit
  cannot express a loss. The mover set is derivable from the live run (`base` moved).
- **Promotion.** Owes the A/B, both numbers (verified 69.24 / self-consistent 77.16), and the
  owner's call on the 17.

## Traps this lane added

- A ruler and the thing it measures can share a blind spot **without sharing code** — here, by
  sharing a *method* (match by class). Reading the 18 disagreements is what found it; the counts
  alone said "568/568".
- "Must stay split" is only a gate if the split is a property of the molecule. Check it under
  renumbering before defending it.
- `GetSubstructMatches` on a `Mol` query reads radicals, charges and isotopes;
  `CanonicalRankAtoms` does not read radicals. Two "symmetry" notions in one function will
  disagree on exactly the molecules perception mangled.
