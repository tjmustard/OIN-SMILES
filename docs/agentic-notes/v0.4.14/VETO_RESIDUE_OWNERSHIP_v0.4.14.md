# v0.4.14 — 183 molecules of the gap are a GENERATOR chirality error filed under a canonicalization bucket

**The `key_equal` bucket is documented as "benign canonicalization — the win reclaimed". For 183
of its 439 members (41.7%) that is false.** The round trip built the **enantiomer**, and the bucket
cannot see it because the comparison key folds reflection deliberately.

---

## 1. The claim this measurement exists to stop being made without evidence

`tools/veto_outcome_audit.py` established that all **222** molecules the parity veto reverts are
`vetoed_collapse` rather than `declined_*` (§1 of [`RESIDUE_RESCOPED_v0.4.14.md`](RESIDUE_RESCOPED_v0.4.14.md)).

The tempting next sentence — *"so those 4.44 points are a generator problem"* — **does not follow**.
`vetoed_collapse` says only: *this structure's mirror encodes differently today, and folding would
make it encode identically.* That is a property of **one structure**. It says nothing about the
relationship between the **input** and its **round trip**, which is the thing being re-filed on the
roadmap.

Two possibilities, with completely different owners:

| | meaning | owner |
|---|---|---|
| **MIRROR_MATCH** | the generated structure really is the input's mirror image | the **generator** |
| **NOT_A_MIRROR** | same enantiomer, two labelings this fold cannot unify | a future **encoder** lane |

## 2. The test

`tools/veto_residue_chirality.py` — new this release. All three strings are re-derived from
coordinates with the fold OFF (the rotation-only labeling):

```
s1   = encode(input.xyz)
s2   = encode(generated.xyz)          the sweep's stored structure
s1_m = encode(mirror(input.xyz))      z-negation, matching fold_parity._mirror_coords

s1_m == s2  ->  MIRROR_MATCH
s1_m != s2  ->  NOT_A_MIRROR
s1  == s2   ->  UNEXPECTED_IDENTICAL  (reported, never dropped -- dropping it would
                                       flatter whichever verdict remained)
```

```bash
$V tools/veto_residue_chirality.py --outcomes $R/veto_outcomes.json --sweep $SW \
    --dataset $MAIN/tmCAT-tmPHOTO_xyz_dataset/cat \
    --dataset $MAIN/tmCAT-tmPHOTO_xyz_dataset/photo --out-json veto_residue_chirality.json
```

**The control:** the fold-OFF re-encode must reproduce the frozen `smiles_1` / `smiles_2_indep`, or
the verdict is describing encoder drift since v0.4.8 rather than the round trip.

## 3. Result

```
=== VETO-REVERTED RESIDUE, n=222 of 222 ===
  excluded: unavailable 0, drift 0
  MIRROR_MATCH             183  (82.4%)
  NOT_A_MIRROR              39  (17.6%)
```

**Full coverage. No exclusions.**

| | n | pts | owner |
|---|---:|---:|---|
| `MIRROR_MATCH` — the round trip built the enantiomer | **183** | **3.66** | generator |
| `NOT_A_MIRROR` — same enantiomer, unfoldable labeling | **39** | **0.78** | encoder (parity-aware canonicalization) |
| | 222 | 4.44 | |

### Why `MIRROR_MATCH` cannot be a false positive here

The verdict is only reachable for a molecule the encoder already treats as **chiral**, and that is
forced twice over rather than assumed:

- `UNEXPECTED_IDENTICAL` is checked first, so `s1 != s2` holds on every classified molecule. With
  `s1_m == s2` that gives `s1_m != s1` — the mirror encodes differently from the original, which is
  the encoder's own definition of a resolved enantiomer pair.
- Independently, `vetoed_collapse` **requires** `s_rot != s_rot_m` (`fold_parity.resolve`'s left
  conjunct). Every molecule in this population had already satisfied it, through a different code
  path, before this tool ran.

So `MIRROR_MATCH` is: *a chiral complex whose round trip encodes byte-for-byte as its own mirror.*

## 4. What this corrects

`roundtrip_bucket_report`'s docstring calls `key_equal` *"benign canonicalization — the win
reclaimed"*, and the roadmap's gap decomposition carries it as a **string** problem. For 183
molecules that description is wrong in the most consequential possible direction: **`byte_exact`
failing on them is CORRECT**, and any change that "recovers" them would be raising the headline by
deleting the fact that the generator produced the wrong enantiomer.

The mechanism is the same blindness v0.4.11 documented, one level up:
`compare._parse_vertex_colors` folds reflection deliberately, so an enantiomer pair has an **equal
comparison key** and lands in `key_equal` by construction. `accept_fn` decides by that key, which
means **the generator accepted a mirror-image structure and the harness recorded it as a
same-isomer string difference.**

This is the third time this project has found its headline resting on a metric that folds the axis
under test — after v0.4.8's scored-vs-honest re-baseline and v0.4.11's fold. The transferable form:

> **A lossy key must never be reused as an acceptance predicate for an axis it folds** — and when a
> bucket's *name* asserts a cause ("benign canonicalization", "rdkit_canonical"), that name is a
> hypothesis, not a measurement.

## 5. Roadmap consequence

The 24.12-point gap needs re-filing. `slot_renumber` was carried as one 325-molecule
canonicalization block; it is three things with three owners:

| | n | pts | owner |
|---|---:|---:|---|
| frozen resonance forms | 103 | 2.06 | **encoder — closed this release** (78 taken, 25 residue) |
| generator built the enantiomer | 183 | 3.66 | **generator** — new lane, does not belong to the encoder ladder |
| unfoldable same-enantiomer labeling | 39 | 0.78 | encoder — needs a parity-aware canonicalization |

Add the separate re-scoping of `rdkit_canonical` (114 / 2.28 pts → **η-set perception drift**, not
canonicalization at all) and **5.94 of the gap's 24.12 points move off the encoder ladder
entirely.**

⚠ **None of this changes `byte_exact`.** It changes who owns the remaining distance to 100%, which
is the number the ladder is scheduled against.
