# v0.4.14 — the `key_equal` residue, re-scoped

**The charter's arithmetic was right and its attribution was wrong.** v0.4.14 was chartered on
`rdkit_canonical` (114) plus "~90 frozen resonance forms", **~4.08 points combined**. Both blocks
are re-derived here from `results-v0.4.8-honest`, and neither is what its name says.

Every number below states the command that produced it. Nothing here is a sweep: the change under
test is encoder-side and **generator-neutral**, verified before anything was measured.

```bash
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
V=$MAIN/.venv/bin/python
SW=$MAIN/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest
export PYTHONPATH=$PWD/src
find $MAIN/tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k -xtype l | wc -l   # 0 ✓
```

---

## 1. The 325 residual `slot_renumber` is two populations, not one

`tools/veto_outcome_audit.py` — new this release. It answers the question **no bucket report can**,
because every reverted molecule lands in the same bucket regardless of why it reverted.

```bash
$V tools/veto_outcome_audit.py --sweep $SW \
    --dataset $MAIN/tmCAT-tmPHOTO_xyz_dataset/cat \
    --dataset $MAIN/tmCAT-tmPHOTO_xyz_dataset/photo --out-json veto_outcomes.json
```

| | n | pts | what it is |
|---|---:|---:|---|
| the fold cannot reach it | **103** | 2.06 | 101 are `same_colour_DIFFERENT_rank` — the frozen resonance form |
| the fold reaches it, the **veto reverts** it | **222** | 4.44 | `vetoed_collapse`, 222 of 222 |
| | **325** | **6.50** | = 496 base − 171 v0.4.13 gains |

**The accounting closes exactly against v0.4.13**, which is what makes it trustworthy rather than
merely plausible: `496 − 393` fold-reached `= 103`, and `393 − 171` veto-kept `= 222`.

### Why the split matters more than the sizes

`fold_parity.resolve` has **five** outcomes and only one of them means "the fold is unsafe here".
Three mean *"the veto never got a reading"* — and a `declined_*` molecule emits the **byte-identical
conservative string** that a `vetoed_collapse` molecule does, so no bucket, no golden and no mirror
audit can tell them apart. `fold_parity`'s own docstring records a build of that module which
declined on 18 of 18 movers while three fixture tests passed.

`levers.py` asserts that *"a corpus reading of `declined_* = 0` is what proves it alive"*. **That is
now measured rather than asserted:**

```
=== VETO OUTCOMES, n=393 of 393 movers ===     (0 unavailable, 0 drift)
  kept (fold survived):  171
  reverted by the veto:  222

  decided_against      222  (100.0%)
  no_evidence            0  (0.0%)
```

Full coverage, no exclusions, and it independently reproduces v0.4.13's 171/222 split to the
molecule — an instrument that was broken could not have landed on that number.

---

## 2. `rdkit_canonical` (114) is not RDKit's write order. It is haptic denticity drift.

The bucket is the `key_equal` sub-split's **residual** — "not `fragment_reorder`, not
`slot_renumber`, not `winding_star_drift`" (`roundtrip_bucket_report._key_equal_subclass`). The
name was a guess about what falls through, and it is wrong.

Counting how many atoms carry each `{n}` marker in the two strings:

| | n | share |
|---|---:|---:|
| **denticity differs** — the two strings mark a different number of atoms for a slot | **92** | 80.7% |
| …and every one of those 92 involves a **haptic** slot (a slot on >1 atom) | 92 | 100% |
| same multiplicities, difference is elsewhere | 22 | 19.3% |

```
ABOZEW_comp_0
  1: ...[cH]{2>}1[cH]{2}[cH]{2}c{2}(...)[cH]{2}1     <- 5 ring atoms in the eta set
  2: ...c1c[cH]{2}c{2}(...)[cH]{2>}1                 <- 3
BARRES_comp_0
  1: c1ccc{3}2[cH]{3}[cH]{3>}[cH]{3}c{3}2c1          <- 6
  2: c1ccc{3}2cccc{3>}2c1                            <- 2
```

This is a **perception** difference, not a serialization one: the generated geometry has the π-ring
slipped or tilted, so fewer ring carbons fall inside the metal-bonding cutoff. It is the charter's
own carried trap — *"the deepest risk is connectivity, not serialization … 210/572 generated
structures hold a ligand within 0.1 Å of that cutoff"* — showing up as a named bucket.

**Consequence for the ladder: these 114 are not an encoder-canonicalization lane.** Re-canonicalizing
the string cannot change how many atoms the perceiver bonded to the metal. Two honest routes exist,
and both are somebody else's release: make the η-set determination robust to slip (perception), or
make the generator reproduce the input's ring geometry (generator). **No lever was built for them
here**, and the charter's 2.28 points for this block should not be counted against v0.4.14.

---

## 3. What was built: `OIN_RESONANCE_DONOR_FOLD` (**promoted default-ON**)

> Built default-OFF, measured, and promoted on the evidence in
> [`LANE-01-resonance-fold.md`](LANE-01-resonance-fold.md): **+1.56 points, 78 molecules, 0 in a
> bad direction**, with a mirror audit at **100% coverage of the moved population**.

The 103 the fold cannot reach fail on one specific condition — 101 of 103 — and it is a
bookkeeping artifact rather than chemistry:

```
ALAJON_comp_0   CC(=O{0})C=C(C)O{1}          acac      strict ranks 2 / 3
AROKUP_comp_0   O{0}S(=O)(=O{2})c1ccccc1...  sulfonate strict ranks 2 / 0
```

`_donor_swap_permutations` requires both donors to sit in one
`CanonicalRankAtoms(breakTies=False)` class. Which oxygen the perceiver wrote as the ketone is a
property of the **Kekulé structure**, not of the ligand — the real ligand is delocalized and its two
donors are the same atom.

`canonical_slots._skeleton_ranks` ranks the fragment's **constitutional skeleton**: bond orders,
aromatic flags, formal charges and hydrogen counts erased; **connectivity, element and chiral tag
kept**. So it merges what resonance makes equivalent and still refuses what constitution makes
different:

| | strict | skeleton | |
|---|---|---|---|
| acac O / O | separate | **merge** | resonance pair |
| carboxylate O / O | separate | **merge** | resonance pair |
| sulfonate O / O | separate | **merge** | resonance pair |
| ester `-O-` / `=O` | separate | separate | 2-connected vs terminal |
| ether O / ketone O | separate | separate | different neighbourhoods |
| amide N / O | separate | separate | different elements |

Grouping is **union-find** (`_merge_classes`), not a composite key, so the partition can only get
*coarser*. A composite key could move a slot into a different bucket and **lose** a labeling the
shipped encoder already reaches — the widening must be a superset or it is not a widening.

### What the chiral-tag retention does and does not buy

Stated because the attractive reading is false, and this project has shipped that mistake before.
Retention means the widening does not **discard** stereochemical information the strict ranking
used. It does **not** mean stereochemically distinct donors can never merge. Measured on a diol
pair: the C2-symmetric `(R,R)` arms do **not** merge (over-conservative — a missed fold, the safe
direction), while the meso `(R,S)` arms **do**, because they are enantiotopic. Folding an
enantiotopic pair is a reflection, and the guard against that is `fold_parity`'s per-molecule veto.
`test_resonance_donor_fold::test_the_skeleton_ranking_consumes_those_tags` states the limit in its
own docstring rather than leaving it to be rediscovered.

Safety is **inherited, not re-argued**: the lever only widens a candidate set the veto already
polices, so it carries the same coupling invariant and cannot be promoted ahead of the veto.

---

## 4. No sweep was owed, and there is a test rather than a judgement call

```bash
$V tools/fold_key_invariance.py --sweep $SW \
    --lever OIN_RESONANCE_DONOR_FOLD --holding OIN_CANONICAL_DONOR_FOLD
```
```
  strings compared          : 9669
  strings the LEVER MOVED   : 228   <-- must be > 0 or the lever is not wired
  strings whose KEY CHANGED : 0
  ✅ GENERATOR-NEUTRAL
```

`--lever` / `--holding` are new. **`--holding` is not a convenience**: `OIN_RESONANCE_DONOR_FOLD`
only widens a candidate set `OIN_CANONICAL_DONOR_FOLD` creates, so measuring it against a fold-OFF
baseline would report the *fold's* movement as the widening's — a larger, more attractive, wrong
number.

---

## 5. Gate coverage of the moved population — the number, not the verdict

179 molecules move a string under the widening.

| golden | movers covered | |
|---|---:|---|
| `gate_v047_arm1_golden.tsv` | **0 of 62** (0.0%) | a PASS here is "no regression", **not** "it works" |
| `gate_v047_arm2_golden.tsv` | 2 of 100 (2.0%) | |
| `gate_v049_arm2_golden.tsv` | **12 of 325** (3.7%) | these rows need individual re-freezing on promotion |

ARM 1 is blind to this change for the same reason it was blind to v0.4.13's, and stating that up
front is cheaper than discovering it after quoting the gate.

---

## 6. An instrument caught printing plausible nothing — the fifth in two releases

The gate-coverage table above was first computed from the **main checkout**, whose `src/` does not
carry the lever. It printed:

```
string-level movers: 0
  gate_v047_arm1_golden.tsv   0 of 62 fixtures are movers (0.0%)
```

A clean, well-formatted, entirely false table — and "0 movers" is *also* what a genuinely inert
lever prints. The re-run from the worktree reads **179**. The guard now used in every ad-hoc probe
in this release, and worth more than the finding it protected:

```python
import oinsmiles.oin.canonical_slots as cs
assert hasattr(cs, "_skeleton_ranks"), "WRONG CHECKOUT: this src/ has no _skeleton_ranks"
print("src under test:", cs.__file__)
```

> **Before quoting an instrument, ask what a BROKEN version would print. If that is the same
> thing — or worse, something more attractive — you have not measured anything yet.**

---

## 7. Roadmap consequences

1. **`rdkit_canonical` 114 / 2.28 pts is misfiled.** It is η-set perception drift, not
   canonicalization. It cannot be closed by an encoder-string change.
2. **222 molecules / 4.44 pts of `slot_renumber` are `vetoed_collapse`.** They are not reachable by
   *this* fold, and any change that "recovers" them re-does v0.4.11's refuted work. Whether they are
   a generator chirality error or an encoder gap is settled in
   [`VETO_RESIDUE_OWNERSHIP_v0.4.14.md`](VETO_RESIDUE_OWNERSHIP_v0.4.14.md) — it is **not** implied by
   `vetoed_collapse`, which is a statement about one structure and its mirror, not about the
   relationship between the input and its round trip.
3. **The genuinely encoder-side residue of `slot_renumber` is 103 molecules / 2.06 points**, and
   that is what this release's lever addresses.
