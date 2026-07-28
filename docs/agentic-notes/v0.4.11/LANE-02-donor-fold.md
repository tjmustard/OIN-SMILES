# v0.4.11 Lane 2 — the within-fragment donor fold: **built, measured, REFUTED**

> 🔴 **The fix v0.4.5 specified is not safe as specified. It collapses enantiomers.**
> `OIN_CANONICAL_DONOR_FOLD` ships **default OFF** and must not be promoted.

The fold works exactly as designed and buys **+7.86 `byte_exact` points** (393 molecules,
movement in one direction only, comparison key untouched). Then a corpus mirror audit found it
merges a structure with its mirror image in **19 of 250 molecules (7.6%)** on a uniform draw —
and, run directly on the 393 molecules the fold claims as wins, in **221 of 393 (56.2%)**.

> **More than half the gain IS the damage.** An independent geometric oracle confirms 18/19 of
> the uniform-draw collapses and 26/30 of a sample of the 221 are genuinely chiral. At most
> **~172 of the 393 (43.8%, ≈3.44 points)** are safe, and even those are unproven without a
> parity filter.

This is the fourth release in this project to end by refuting its own plan.

---

## 1. What was built

v0.4.5 Lane 2 specified the fix in writing and declined to implement it:

> fold, in addition to the rotation group, permutations exchanging donors **within one
> fragment** that are (a) in the same `breakTies=False` symmetry class and (b) the same colour.

`canonical_slots._donor_swap_permutations` builds that swap set;
`canonical_slot_relabeling` composes it with the rotation group (`mapping[s] = p(d(s))`), so
the donor fold runs *inside* the freedom the rotation fold leaves. Identity is always in the
set, so the candidate set is a superset of the rotation-only one — which is why lever-OFF
output is byte-identical by construction, not by luck.

## 2. What it buys — measured, and real

Offline simulation over the 5000 stored string pairs of `results-v0.4.8-honest` (the v0.4.8
precedent: same stored strings through both arms, so nothing generator-side moves underneath).

| before | after | n |
|---|---|---:|
| **`key_equal/slot_renumber`** | **`byte_exact`** | **393** |
| `key_equal/slot_renumber` | unchanged | 103 |
| every other bucket | unchanged | — |
| **moved in any other direction** | | **0** |

`byte_exact` **3623 → 4016 = 72.46% → 80.32%, +7.86 points.** `facmer_divergent` **16 → 16**.

Predicted 377, delivered 393. Cross-tabulated against Lane 1's verdict: **377/377 of the
`automorphism` class fixed** — the prediction is exact where it applies — plus **16** of the
`distinct_donors_LOCAL` class, which are multi-fragment molecules where a *global rotation*
resolves the locally-distinct fragment while donor swaps resolve the others. That is the
mechanism `atom_verdict`'s own docstring names (RUBTIS: *"a rotation acts on all fragments at
once"*), so a per-fragment verdict is a **lower bound** on reachability and behaved as one.

Supporting invariants, all measured over the 992 strings of the population:

- comparison key changed: **0**; failed to re-key: **0**
- combinatorial `cap=4096` tripped: **0** (largest swap set **16** — a 256× margin, so no
  silent degradation)
- ARM 1 **#DONE 62, PASS byte-identical**; ARM 2 **#DONE 90, PASS 90 gated** (lever OFF)

## 3. 🔴 What it costs — the refutation

```bash
PYTHONPATH=$PWD/src .venv/bin/python tools/mirror_audit_donor_fold.py \
    --dataset tmCAT-tmPHOTO_xyz_dataset/cat --n 250 --seed 7
```

| verdict | n | share |
|---|---:|---:|
| `achiral_or_preexisting_fold` | 157 | 62.8% |
| `distinct_both_arms` | 73 | 29.2% |
| 🔴 **`REGRESSION_raw_collapsed`** | **19** | **7.6%** |
| `encode_failed` | 1 | 0.4% |

The regression criterion is an implication, not a pass count — `OFF_distinct and not
ON_distinct` — so a pair the shipped encoder *already* folds (metal Δ/Λ, whose descriptor is
held off) is correctly **not** counted against this lever.

### The 19 are not an artifact — an independent oracle says 18 are genuinely chiral

`tools/injectivity/oracle.py::is_distinct_enantiomer` decides chirality from **geometry and
topology only** (minimum proper-rotation mirror RMSD over graph automorphisms), so it shares
no machinery with the encoder under test.

| | n |
|---|---:|
| mirror is a genuinely distinct isomer ⇒ **real regression** | **18** |
| achiral ⇒ the fold was *right* to collapse it | 1 (`IPEQOJ_comp_0`) |

15 of the 18 hit the oracle's 4000-automorphism cap, which it self-flags as *"verdict may be
unreliable"*. **Three do not**, and those alone settle it:

| molecule | mirror RMSD | automorphisms | cap hit |
|---|---:|---:|---|
| `BIWDIV_comp_0` | 1.52 Å | 8 | no |
| `CIHVAT_comp_0` | 1.22 Å | 16 | no |
| `OJEKET_comp_0` | 2.52 Å | 864 | no |

On the most conservative reading available, the fold provably collapses real enantiomers.

### Run on the 393 claimed gains, the damage rate is 56.2%

The uniform draw asks *"how often does the fold hurt a random molecule?"*. The sharper question
is *"how much of the gain is the damage?"* — so the same audit was run on exactly the 393
molecules the fold moves into `byte_exact`:

| verdict over the 393 gains | n | share |
|---|---:|---:|
| 🔴 **`REGRESSION_raw_collapsed`** | **221** | **56.2%** |
| `distinct_both_arms` (safe gain) | 95 | 24.2% |
| `achiral_or_preexisting_fold` (safe gain) | 77 | 19.6% |

Oracle-checked on a random sample of 30 of the 221: **26 genuinely chiral, 4 achiral** — so
roughly **192 of the 221 are real chirality losses**, not classification noise.

**Net: at most ~172 of 393 gains (≈3.44 points) are safe; the remaining ≈4.42 points of the
headline are paid for in destroyed stereochemistry.** A canonicalization that buys byte-identity
by merging enantiomers is precisely what §4's Rule forbids:

> *This release may impose a choice where the encoder currently has none. It may not merge two
> things that differ.*

`BIWDIV_comp_0`, a Co(III) bis(tridentate) — the mirror differs *only* by the `{4}`/`{5}`
exchange on the second ligand, and the fold canonicalizes exactly that away:

```
off self  : [Co_OCT].…N{0}…N{1}…n{2}1.…N{4}…N{5}…n{3}1
off mirror: [Co_OCT].…N{0}…N{1}…n{2}1.…N{5}…N{4}…n{3}1     <- distinct
on  self  : [Co_OCT].…N{0}…N{1}…n{2}1.…N{4}…N{5}…n{3}1
on  mirror: [Co_OCT].…N{0}…N{1}…n{2}1.…N{4}…N{5}…n{3}1     <- COLLAPSED
```

## 4. Why the safety argument was wrong — the transferable lesson

The argument in `_donor_swap_permutations`' first draft was:

> two donors in the same `breakTies=False` symmetry class of their fragment ⇒ exchanging their
> slots denotes the same molecule.

**False.** `CanonicalRankAtoms` computes the symmetry of the **isolated ligand graph**. Two
donors automorphic in the free ligand occupy two *distinct vertices* whose relationship to the
other ligands is chirality-bearing, so the vertex permutation the exchange induces can be an
**improper** operation on the coordination polyhedron. Folding by it merges a structure with
its mirror image.

> **A fragment's automorphism says nothing about the PARITY of the vertex permutation it
> induces.** v0.4.5's restriction to proper rotations was not conservatism — it was the
> load-bearing correctness condition, and narrowing the scope *within a fragment* is not a
> substitute for it.

The three scope conditions are all necessary and jointly **insufficient**. What is missing is a
**reflection-parity** filter: admit `d` only when the labeling it produces is related to the
original by a *proper* operation on the whole coordination sphere. That is the concrete next
step this lane hands forward, and it is testable with the tooling built here.

### Why every other signal said "clean"

| check | verdict | why it could not see this |
|---|---|---|
| transition matrix, 5000 molecules | 393 gains, 0 elsewhere | `byte_exact` cannot see chirality |
| comparison key, 992 strings | 0 changed | the key **deliberately** folds this axis (`_parse_vertex_colors` is colour-blind to which donor holds which slot) |
| ARM 1 / ARM 2 | 62/62, 90/90 PASS | lever OFF — correct, and silent about ON |
| `ZUMNEC` + `fac-Ir(ppy)₃` fixtures | pass | the wrong fixtures; neither has the vulnerable motif |
| stratified audit (slow band, 300) | 0 regressions in the first 200 | **the wrong sample** |
| **uniform draw (250)** | **19 collapses** | the only instrument that could |

> 🔴 **The +7.86 points looked free precisely because the metric and the key are both blind to
> the axis being destroyed.** For a roadmap whose target *is* `byte_exact` 100%, that is the
> finding that outranks the points: **this metric can be raised by deleting information.** Any
> future canonicality lever must be mirror-audited on a uniform draw before its points are
> quoted.

**Sampling matters as much as the check.** The stratified runtime cohort is selected for being
*slow*, not for carrying stereochemistry; it returned a clean bill for 200 molecules while a
uniform draw of 250 found 19 collapses. A corpus audit on the wrong stratum is not a corpus
audit.

## 5. Status and disposition

- `OIN_CANONICAL_DONOR_FOLD` stays in `_HELD_OFF`, with the refutation recorded in its entry.
- Lever OFF is byte-identical on both gate arms, so the shipped encoder is **unchanged**.
- `tools/mirror_audit_donor_fold.py` is the durable deliverable: it is the instrument that
  caught this, and it should gate every future canonicality lever.
- `TestResidualClassIsOutOfReachByDesign` remains **inverted** — the class *is* reachable, which
  is true; what is not true is that reaching it this way is safe. Its docstring says so.
- **Do not promote without a reflection-parity filter**, and re-run the uniform mirror audit
  after building one.

## 6. Reproducing

```bash
cd <repo>; V=$PWD/.venv/bin/python; export PYTHONPATH=$PWD/src
D=$PWD/tmCAT-tmPHOTO_xyz_dataset

# the refutation                                                  (~35 min)
$V tools/mirror_audit_donor_fold.py --dataset $D/cat --n 250 --seed 7
#   expect: 19 REGRESSION_raw_collapsed of 250; exit status 1

# the independent chirality verdict on those 19
$V -c "from tools.injectivity.oracle import is_distinct_enantiomer as f; \
       print(f('$D/cat/.../BIWDIV_comp_0.xyz'))"

# the gain, for the record                                        (~6 min, offline)
#   see the transition matrix in section 2; regenerate with the snippet in
#   docs/agentic-notes/v0.4.11/CANONICALITY_v0.4.11.md
```
