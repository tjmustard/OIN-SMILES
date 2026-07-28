# v0.4.11 Lane 1 — the 496 `slot_renumber` molecules, classified by mechanism

**The charter's headline risk is refuted, in the release's favour.** It said: *"the taxonomy
has only ever run on re-presentation pairs, never on round-trip pairs — if `diff_occupancy`
dominates the 496, the roadmap's 9.92-point canonicality label is wrong, and saying so is this
release's result."*

`diff_occupancy` is **zero**. All 496 are `same_vcolor_identical` — the exact class the
v0.4.11 fold targets. Nothing in this block is upstream, and nothing in it is a defect.

Every number below is reproducible with the command that produced it.

---

## 1. What was built

`tools/slot_drift_mechanism.py` gained a `--roundtrip <results-dir>` input mode. The taxonomy
(`mechanism()`) and the atom-level verdict (`atom_verdict()`) are untouched — extending the
shipped classifier rather than writing a second one, because two classifiers that must agree
are two classifiers that will drift.

Selection delegates to the shipped predicate, `roundtrip_bucket_report.classify(rep,
score="honest")`, so this tool and the bucket report cannot disagree about what is in the
bucket. `--expect N` aborts unless exactly `N` pairs are selected, which pins the population.

```bash
PYTHONPATH=$PWD/src .venv/bin/python tools/slot_drift_mechanism.py \
    --roundtrip tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest \
    --subclass slot_renumber --expect 496 --explain-distinct
```

Determinism: two runs are byte-identical (`diff` empty).

## 2. The mechanism histogram — all 496, no sampling

| class | n | share | reachable by a within-fragment fold? |
|---|---:|---:|---|
| `same_vcolor_identical` | **496** | **100.0%** | **yes — this is the fold's target class** |
| `diff_geometry` | 0 | — | no (upstream, the 3D fit) |
| `diff_occupancy` | 0 | — | no (upstream, and legitimate) |
| `diff_colors` | 0 | — | no (downstream of the canonical body) |
| `postpass_BUG_diverges` | **0** | — | would be a defect |
| **sum** | **496** ✓ | | |

`postpass_BUG_diverges` is **0**, so the post-pass itself is not at fault anywhere in this
population.

## 3. The atom-level verdict — and the decomposition that sizes the ladder

| verdict | n | points | where it belongs |
|---|---:|---:|---|
| `automorphism` | **377** | **7.54** | **v0.4.11** — the within-fragment donor fold |
| `distinct_donors_LOCAL` | 118 | 2.36 | split below |
| `unparsable` | 1 | 0.02 | `ZIBJOM_comp_0`, a borane cluster |
| **sum** | **496** ✓ | **9.92** ✓ | |

> **N of 496 reachable: 377 of 496 (76.0%, 7.54 points) are `same_vcolor_identical` with an
> atom-level `automorphism` verdict and are therefore reachable by a within-fragment fold.**

That is the sentence Lane 2 quotes. It is the release's predicted size.

### The 118, re-tested resonance-insensitively — a new result

`distinct_donors_LOCAL` means the per-fragment, string-only test could not show the two donors
interchangeable. It has never been asked *why*. Re-running the same verdict against a
**resonance-insensitive** copy of the fragment (every bond single, formal charges zeroed —
`flattened_ranks()`) splits it cleanly:

| | n | share of the 118 | points |
|---|---:|---:|---:|
| equivalent once bond orders are flattened | **90** | **76.3%** | 1.80 |
| still distinct after flattening | 28 | 23.7% | 0.56 |

The archetype is `ALAJON_comp_0`, an acetylacetonate:

```
smiles_1       = [Pt_SPL].CC(=O{0})C=C(C)O{1}.Cc1cc(-c2c{2}ccc3ccccc23)n{3}c2ccccc12
smiles_2_indep = [Pt_SPL].CC(=O{1})C=C(C)O{0}.Cc1cc(-c2c{2}ccc3ccccc23)n{3}c2ccccc12
```

acac binds through two chemically equivalent oxygens, but perception freezes one localized
resonance form — writing one O as a ketone and one as an enol — so `CanonicalRankAtoms` puts
them in different symmetry classes. **50 of the 118 carry the delocalizable `O=C-C=C-O`
motif outright**, and 114 of the 118 exchange donors of the *same element*.

**So the 118 are not a soundness class and not a slot-fold problem — they are a ligand-BODY
canonicalization gap**, which is `rdkit_canonical`'s mechanism and therefore **v0.4.14**. This
confirms the ladder's placement and sizes it: v0.4.14 is worth its own 2.28 points **plus**
these 1.80, i.e. ~4.08.

⚠ **Flattening is a diagnostic, not a folding criterion.** It also erases distinctions that
are real (an amide N against an amine N), so a fold keyed on it would merge donors that
genuinely differ — exactly what the v0.4.11 Rule forbids. It is used here only to say *where*
the residual belongs.

## 4. Probe vs round-trip, side by side — the shapes agree

The charter's concern was that a re-presentation population and a round-trip population would
not look alike. They do.

```bash
.venv/bin/python tools/canonicality_probe.py --n 150 --trials 2 \
    --dataset <main-checkout>/tmCAT-tmPHOTO_xyz_dataset --out probe150
.venv/bin/python tools/slot_drift_mechanism.py probe150 --subclass slot_renumber --explain-distinct
```

| | probe (re-presentation) | round-trip (this release) |
|---|---:|---:|
| pairs | 35 | 496 |
| `same_vcolor_identical` | 100.0% | 100.0% |
| `automorphism` | 80.0% | 76.0% |
| `distinct_donors_LOCAL` | 20.0% | 23.8% |
| `unparsable` | 0 | 0.2% |

**Mechanism for the small divergence:** a generated conformer *can* place a donor on a
different vertex, which is why the round-trip share of `distinct_donors_LOCAL` is a few points
higher. It cannot do so often — `diff_occupancy` is 0 — because the round-trip pairs were
pre-filtered to `key_equal`, and a moved arrangement generally breaks the key.

### ⚠ Two of the charter's prior numbers are stale — corrected here

`BASELINE.md` §4 quotes the v0.4.5 probe as **32 pairs / 23 `automorphism` / 7
`distinct_donors_LOCAL` / 2 unparsable**. Re-run today it is **35 / 28 / 7 / 0**.

- The **7 `distinct_donors_LOCAL` reproduce exactly**, and 6 of the 7 are the same molecules
  the charter names (`IMACOO`, `RUBTIS`, `VEXHIR`, `HAVGIW`, `KISQAG`, `ZOSNUS`; `PAGHOV`
  replaces `ZACFER`).
- The **2 unparsable are gone because `OIN_BORON_CAGE` was promoted to default-ON in v0.4.6**,
  after the prior was measured. The v0.4.5 figure was recorded with it off.

Neither correction changes any conclusion; both are recorded so the next reader does not
re-derive them.

## 5. What this lane hands to Lane 2

- **Predicted size: 377 molecules, 7.54 points.** Not 9.92 — the roadmap's block is real but
  only 76% of it is reachable at this seam.
- **Fixture archetype:** `AGUKOD_comp_0` (Rh, SPL, COD across two cis vertices; bodies
  byte-identical, `vcolor` identical, only which alkene arm carries which integer moves).
- **Negative controls** — the fold must **not** move these: the 118 `distinct_donors_LOCAL`
  (they are in different symmetry classes by construction, so a correctly-scoped fold cannot
  reach them) and `ZIBJOM_comp_0` (unparsable).
- **`postpass_BUG_diverges` is 0**, so there is no defect outranking the fix.

## 6. Reproducing everything above

```bash
cd <repo>; V=$PWD/.venv/bin/python; export PYTHONPATH=$PWD/src
D=$PWD/tmCAT-tmPHOTO_xyz_dataset

# the classification, all 496                                     (~4 min, offline)
$V tools/slot_drift_mechanism.py --roundtrip $D/results-v0.4.8-honest \
    --subclass slot_renumber --expect 496 --explain-distinct --out-json l1_496.json

# determinism
$V tools/slot_drift_mechanism.py --roundtrip $D/results-v0.4.8-honest > /tmp/a
$V tools/slot_drift_mechanism.py --roundtrip $D/results-v0.4.8-honest > /tmp/b
diff /tmp/a /tmp/b        # empty

# the probe side-by-side                                          (~9 min)
$V tools/canonicality_probe.py --n 150 --trials 2 --dataset $D --out probe150
$V tools/slot_drift_mechanism.py probe150 --subclass slot_renumber --explain-distinct
```

⚠ The dataset is gitignored, so a **git worktree does not have it** — `canonicality_probe.py`
needs `--dataset <main-checkout>/tmCAT-tmPHOTO_xyz_dataset` explicitly or it exits with
"no .xyz found".
