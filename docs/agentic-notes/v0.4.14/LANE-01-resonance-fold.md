# v0.4.14 Lane 1 — `OIN_RESONANCE_DONOR_FOLD`, built, measured, promoted

**Result: `byte_exact` 75.88% → 77.44%, +1.56 points, 78 molecules, 0 in a bad direction.**

The lane's subject was re-derived before anything was built, and the charter's attribution turned
out to be wrong in both halves — see [`RESIDUE_RESCOPED_v0.4.14.md`](RESIDUE_RESCOPED_v0.4.14.md)
(what the buckets actually contain) and
[`VETO_RESIDUE_OWNERSHIP_v0.4.14.md`](VETO_RESIDUE_OWNERSHIP_v0.4.14.md) (who owns the residue).
This file is the lane's own record: what was built, what gates it cleared, and what it costs.

```bash
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
V=$MAIN/.venv/bin/python
SW=$MAIN/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest
R=$MAIN/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1
export PYTHONPATH=$PWD/src
find $MAIN/tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k -xtype l | wc -l   # 0 ✓
```

---

## 1. The change

`canonical_slots._skeleton_ranks` ranks a fragment's **constitutional skeleton** — bond orders,
aromatic flags, formal charges and hydrogen counts erased; **connectivity, element and chiral tag
kept** — and `_merge_classes` unions that partition with the strict one by union-find, so the
grouping can only get **coarser**. A composite key would have been the obvious spelling and is
wrong: it could move a slot into a different bucket and *lose* a labeling the shipped encoder
already reaches. A widening must be a superset or it is not a widening.

Safety is **inherited, not re-argued**. The lever only widens a candidate set that
`OIN_FOLD_PARITY_VETO` already polices per molecule, so it carries the same coupling invariant and
cannot be promoted ahead of the veto —
`test_resonance_donor_fold::TestResonanceFoldInheritsTheVetoCoupling`.

## 2. Points

```bash
$V tools/resonance_transition_sim.py --sweep $SW \
    --dataset $MAIN/tmCAT-tmPHOTO_xyz_dataset/cat \
    --dataset $MAIN/tmCAT-tmPHOTO_xyz_dataset/photo \
    --baseline-byte-exact 3794 --out-json $R/resonance_transition.json
```
```
movers measured : 168/179  (unavailable 11, drift 0)
bucket changed  : 78
  key_equal/slot_renumber -> byte_exact    78
byte_exact 3794 -> 3872  (75.88% -> 77.44%, +1.56 points)
moved in a BAD direction: 0
```

**Both arms re-encode from coordinates**, because the state being compared against is v0.4.13's
shipped default path, which runs the parity veto and therefore needs a conformer. That is why this
is a new tool rather than `fold_transition_sim.py` with a flag.

### The reconciliation, which is what makes the number trustworthy

| the 78 gains land in… | n |
|---|---:|
| class A — the 103 the v0.4.13 fold cannot reach | **78** |
| class B — the 222 the veto reverts | 0 |
| `rdkit_canonical` | 0 |
| anywhere else | 0 |

The 11 unavailable molecules are in **none** of A, B or `rdkit_canonical`, so they cannot be
understating the gain. Exactly **1** class-B molecule moved a string and it did **not** change
bucket — the veto held against the widening, which was the specific regression mechanism this
lane was most exposed to.

### Residue characterized, not merely reduced

103 − 78 = **25** remain: **20** are still a slot swap the skeleton deliberately keeps apart, **5**
never fire the widening at all.

## 3. No sweep was owed, and it was a test rather than a judgement call

```bash
$V tools/fold_key_invariance.py --sweep $SW \
    --lever OIN_RESONANCE_DONOR_FOLD --holding OIN_CANONICAL_DONOR_FOLD
```
```
strings compared : 9669 · LEVER MOVED : 228 · KEY CHANGED : 0 · skipped : 0
✅ GENERATOR-NEUTRAL
```
The `moved` count is load-bearing: a lever that never fired would also print 0 key changes.

## 4. 🔴 Mirror audit — and its COVERAGE, stated before its verdict

Two draws, and only one of them can see this change.

| draw | coverage of the 179 movers | control (reso OFF) | resonance (reso ON) |
|---|---:|---|---|
| uniform `cat/`, n=250, seed 7 | **1 / 250 = 0.4%** | 157 / 92 / 1 · CLEAN | 157 / 92 / 1 · CLEAN |
| **mover-enriched cohort, n=179** | **179 / 179 = 100%** | 108 / 71 · **0 REGRESSIONS** | 108 / 71 · **0 REGRESSIONS** |

*(`distinct_both_arms` / `achiral_or_preexisting_fold` / `encode_failed`)*

**The uniform draw's identical tallies are not safety evidence.** At 0.4% coverage they say only
that the lever does not damage molecules it never touches — this is v0.4.13's ARM 1 failure
(0 of 62 fixtures were fold-movers) reproduced exactly, and it would have been quoted as a clean
result if the coverage had not been computed. It is still worth running: it reproduces v0.4.13's
frozen `mirror_cat_promoted.json` **exactly** (157/92/1), which is what proves the instrument and
the draw are alive rather than silently broken.

On the draw that *can* see the change:

- **per-molecule verdicts differing between arms: 0.** A tally can stay flat while individual
  verdicts churn underneath it, so the comparison is per molecule, not per bucket.
- **`achiral_or_preexisting_fold` unmoved at 71 in both arms** — the check that the zero is not
  bought by declining everything.
- **33 of the 78 gains are `distinct_both_arms`** — molecules the encoder resolves as chiral, and
  still resolves after the widening. The other 45 are achiral, where folding cannot collapse
  anything. This is the number that separates "safe" from "never fired on anything chiral", and
  without it the 0 regressions would be compatible with a lever that only ever touches achiral
  molecules.

## 5. Gates

| | result | coverage of the moved population |
|---|---|---|
| ARM 1 (`gate_v047.sh arm1`) | **PASS — byte-identical**, `#DONE 62` | **0 of 62** |
| `gate_v047_arm2_golden.tsv` | 2 rows re-frozen | 2 of 100 |
| `gate_v049_arm2_golden.tsv` | 12 rows re-frozen | 12 of 325 |
| `tests/unit` | **1006 → OK** (skipped 3, xfail 5) | — |
| ruff 0.15.20 | format + check clean, 309 files | — |

⚠ The arm2 goldens carry a `MANIFEST_SHA256` line that **arm2 does not verify**, so a stale one
goes unseen. It is recomputed with the re-frozen rows.

## 6. What this lane does NOT claim

- **Chiral-tag retention is not a guarantee against stereochemical collapse.** It means the
  widening does not *discard* what the strict ranking used. Measured on a diol pair, the
  C2-symmetric `(R,R)` arms do **not** merge (over-conservative — a missed fold, the safe
  direction) and the meso `(R,S)` arms **do**, because they are enantiotopic; that case is the
  veto's job. The test docstring states this limit rather than leaving it to be rediscovered.
- **The fixture tests prove three ligand motifs, not corpus safety.** Their coverage is acac,
  carboxylate/sulfonate, and three negative controls. The corpus instrument is the mirror audit.
- **`rdkit_canonical`'s 114 molecules were not addressed** and no lever was built for them. They
  are 80.7% η-set perception drift, which no string canonicalization can reach.

## 7. Predicted vs actual

| | charter | actual |
|---|---|---|
| `byte_exact` | **UP ~3.5–4.1 pts** | **+1.56** |

**A miss, and the reason is the re-scoping rather than the lever underperforming.** The charter's
~4.08 assumed both blocks were reachable encoder-side: `rdkit_canonical` (114 / 2.28) plus a
~90-molecule resonance class (1.80). Measured, `rdkit_canonical` contributes **0** — it is haptic
denticity drift, not canonicalization — and the resonance class is 103 rather than 90, of which
**78 converted**. So the lane took 78 of a genuinely-reachable 103, and the shortfall against the
prediction is 114 molecules that were never reachable by this kind of change.
