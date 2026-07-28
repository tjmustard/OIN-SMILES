# v0.4.11 — `slot_renumber`: the fold that works, and must not ship

**Headline: the fix v0.4.5 specified for the largest block in the gap is not safe as specified.
It buys +7.86 `byte_exact` points by collapsing enantiomers, and neither `byte_exact` nor the
comparison key can see the damage.**

`OIN_CANONICAL_DONOR_FOLD` ships **default OFF**. The shipped encoder is byte-for-byte
unchanged (ARM 1 62/62, ARM 2 90/90).

---

## 1. What the release set out to do

The v0.4.8 honest baseline puts **496 molecules / 9.92 points** — the largest single block of
the 27.54-point gap — in `key_equal → slot_renumber`: pairs that are byte-identical *once the
`{n}` slot integers are blanked*. Same fragments, same order, same winding, same comparison
key. Not a chemistry problem; a determinism gap.

v0.4.5 Lane 2 measured its shape, **specified the fix in writing, and declined to implement
it**, giving three reasons — one of which v0.4.5 Lane 9 later refuted. v0.4.11 was re-pointed
at this block by the owner-accepted `LADDER DECISION 2026-07-27`, taking it from a 0.30-point
release to a 9.92-point one.

## 2. Lane 1 — the 496, classified for the first time

The mechanism taxonomy had only ever run on *re-presentation* pairs (renumber/rotate), never on
*round-trip* pairs whose second string comes from a generated conformer. The charter's stated
risk: if `diff_occupancy` dominates, the 9.92-point canonicality label is wrong.

**It does not.** All 496 are `same_vcolor_identical`; `diff_occupancy`, `diff_geometry`,
`diff_colors` and `postpass_BUG_diverges` are all **0**.

| atom-level verdict | n | pts | belongs to |
|---|---:|---:|---|
| `automorphism` | 377 | 7.54 | v0.4.11 (this fold) |
| `distinct_donors_LOCAL` | 118 | 2.36 | split below |
| `unparsable` (`ZIBJOM`, borane) | 1 | 0.02 | — |
| **sum** | **496** ✓ | **9.92** ✓ | |

**New result — what the 118 actually are.** Re-tested against a resonance-*insensitive* copy of
the fragment (all bonds single, charges zeroed), **90 of 118 become equivalent**. The archetype
is acetylacetonate: acac binds through two chemically equivalent oxygens, but perception freezes
one localized resonance form — one O a ketone, one an enol — so `CanonicalRankAtoms` separates
them. 50 of the 118 carry the delocalizable `O=C-C=C-O` motif outright and 114/118 exchange
donors of the *same element*.

> The 118 are a **ligand-body** canonicalization gap, not a slot-fold problem — `rdkit_canonical`'s
> mechanism, i.e. **v0.4.14**, which is therefore worth its own 2.28 points **plus** these 1.80.

Full detail: `LANE-01-classify-496.md`.

## 3. Lane 2 — the fold works

`_donor_swap_permutations` exchanges donors within one fragment that share a
`breakTies=False` symmetry class and a vertex colour, composed with (not replacing) the
rotation group. Offline simulation over the 5000 stored string pairs:

| | before | after |
|---|---:|---:|
| `byte_exact` | 3623 (72.46%) | **4016 (80.32%)** |
| `key_equal/slot_renumber` | 496 | 103 |
| every other bucket | — | **unchanged** |
| moved in any other direction | | **0** |

`facmer_divergent` **16 → 16**. Comparison key changed on **0 of 992** strings. Combinatorial
cap tripped **0** times. Predicted 377, delivered 393 — **377/377** of the predicted class plus
16 multi-fragment cases where a *global rotation* resolves the locally-distinct fragment, the
mechanism `atom_verdict` documents (RUBTIS). A per-fragment verdict is a lower bound and
behaved as one.

## 4. 🔴 And it must not ship

```bash
PYTHONPATH=$PWD/src .venv/bin/python tools/mirror_audit_donor_fold.py --dataset <dir> --n N
```

| population | collapses | rate |
|---|---:|---:|
| uniform draw, 250 molecules | **19** | 7.6% |
| stratified runtime cohort, 300 molecules | **31** | 10.3% |
| **the 393 molecules the fold claims as wins** | **221** | **56.2%** |

**Both independent cohorts find it**, at comparable rates. The damage is not an artifact of one
sample.

Independent verdict from `tools/injectivity/oracle.py` (chirality from geometry + topology
alone, no shared machinery with the encoder): **18 of 19** uniform-draw collapses and **26 of a
30-sample** of the 221 are genuinely chiral. Three — `BIWDIV`, `CIHVAT`, `OJEKET` — are
cap-free, so the verdict does not rest on the ones the oracle self-flags as unreliable.

**At most ~172 of 393 gains (~3.44 of the 7.86 points) are safe.** The rest is paid for in
destroyed stereochemistry — exactly what the release's own Rule forbids: *may impose a choice
where the encoder has none; may not merge two things that differ.*

### Why the safety argument failed

> two donors in the same `breakTies=False` symmetry class of their fragment ⇒ exchanging their
> slots denotes the same molecule

**False.** `CanonicalRankAtoms` computes the symmetry of the **isolated ligand graph**. The two
donors occupy *distinct vertices* whose relation to the other ligands is chirality-bearing, so
the vertex permutation the exchange induces can be **improper**.

> **A fragment's automorphism says nothing about the parity of the vertex permutation it
> induces.** v0.4.5's restriction to proper rotations was not conservatism — it was the
> load-bearing correctness condition. The three scope conditions are each necessary and jointly
> **insufficient**.

## 5. The finding that outranks the points

| check | said | why it was blind |
|---|---|---|
| transition matrix, 5000 molecules | 393 gains, 0 elsewhere | `byte_exact` cannot see chirality |
| comparison key, 992 strings | 0 changed | the key folds this axis **by design** |
| ARM 1 / ARM 2 gates | 62/62, 90/90 PASS | lever OFF — correct, and silent about ON |
| `ZUMNEC`, `fac-Ir(ppy)₃` fixtures | pass | wrong fixtures; neither carries the vulnerable motif |

> 🔴 **`byte_exact` can be raised by deleting information.** The +7.86 points looked free
> because the metric and the key are blind to the axis being destroyed. For a roadmap whose
> target *is* `byte_exact` 100%, this is a standing hazard, not a one-off:
>
> **Mirror-audit every future canonicality lever before quoting its points.**

### ⚠ A methodology error made during this release, recorded because it nearly became a finding

Mid-run, the stratified audit was read as **"0 regressions in the first 200"** and written up as
evidence that the runtime-stratified cohort was *the wrong sample* — too slow-selected to carry
stereochemistry. **That was wrong.** `mirror_audit_donor_fold.py` prints a per-molecule verdict
only every 50th molecule, so four clean progress lines were mistaken for 200 clean molecules.
Run to completion the same cohort reports **31 collapses (10.3%)** — a *higher* rate than the
uniform draw.

> **A partial run is not a result, and a progress line is not a tally.** The conclusion drawn
> from it — a "wrong stratum" lesson — was an artifact of reading a sampled log as a census. Both
> cohorts detect the defect; nothing here supports a claim about stratum choice.

## 6. Predicted vs actual

| | predicted | actual |
|---|---|---|
| `byte_exact` | **UP**, size conditional on Lane 1's classification | **FLAT — 72.46%.** The lever that would raise it is unsafe and ships OFF |
| `key_equal` down by the same count | yes | flat, for the same reason |
| `facmer_divergent` 16 → 16 | yes | **16 → 16** ✓ (measured with the lever on) |
| Lane 1 reachable fraction | unknown | **377 / 496 = 7.54 pts**, and 377/377 of it was reached |
| suite | ≥ 946 | **954 OK** |

The charter predicted the release would either raise `byte_exact` or refute the 9.92-point
label. It did **neither**: the label is right, the block is reachable, and reaching it this way
costs more than it pays.

## 7. What v0.4.12+ inherits

1. **A concrete next step:** filter the swap set by **reflection parity** — admit a swap only
   when the labeling it produces is related to the original by a *proper* operation on the whole
   coordination sphere. Then re-run the uniform mirror audit **before** quoting points.
   `TestDonorFoldCollapsesEnantiomers` pins the current defect and must be **inverted, not
   deleted**, when that lands.
2. **A re-sized v0.4.14:** 2.28 points of `rdkit_canonical` **+ 1.80** points of frozen-resonance
   `slot_renumber` = ~4.08.
3. **A standing gate:** `tools/mirror_audit_donor_fold.py`, on a uniform draw, for any lever
   that touches canonicalization.
4. **The carry-forward licence is INTACT.** No lever was promoted and no default answer moved,
   so `BASELINE.md` §1's licence survives and v0.4.12 does **not** owe a 55 CPU-h sweep.
