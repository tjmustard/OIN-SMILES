# v0.4.16 Lane 2 — the 3.74 "unexplained" points, classified · **82% is PERCEPTION, not construction**

**Verdict: all 187 molecules land in a named class, 0 unaccounted, every class read.** The block
the roadmap has carried for three releases as *"attachment fine, independent re-perception still
disagrees — nobody knows why"* is **70.9% bond-order/aromaticity perception on a geometry whose
heavy-atom graph is already correct.**

## What was measured

Offline over the frozen `results-v0.4.14-sweep`, for the same reason the attachment split is: the
classification is a property of the stored strings, and re-running the generator would not make it
more true — it would make it a different corpus.

```bash
V=$PWD/.venv/bin/python
$V tools/attach_class_audit.py --results-dir tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep \
    --characterise --out measurements/v0.4.16/audit_characterise_v0416.json
```

The classifier asks **which component of the round-trip key disagrees** — metal element, ligand
body multiset, coordination-geometry code, or the arrangement signature — and for a body
difference, whether the **heavy-atom graph** survives it.

## The result

| class | INTACT (MEDZUR) | BOUNDARY | combined | pts | mechanism that could reach it |
|---|---:|---:|---:|---:|---|
| `PERCEPTION` | 79 | 43 | **122 (70.9%)** | 2.44 | **perception** — heavy-atom graph agrees; bond orders / aromaticity / charge / H differ |
| `LIGAND_SPLIT` | 9 | 11 | 20 (11.6%) | 0.40 | **construction** — a bond broke, one ligand became several |
| `GEOM_CODE` | 15 | 4 | 19 (11.0%) | 0.38 | **perception / selection** — the polyhedron was read differently |
| `STEREO_INVERSION` | 4 | 2 | 6 (3.5%) | 0.12 | **construction / selection** — same connectivity, inverted stereocentre |
| `SKELETON` | 3 | 1 | 4 (2.3%) | 0.08 | **construction** — graph differs even ignoring stereo |
| `LIGAND_MERGE` | 0 | 1 | 1 (0.6%) | 0.02 | **construction** — fragments fused |
| **total** | **110** | **62** | **172** | **3.44** | |

**By family: PERCEPTION 141 (82.0%) · CONSTRUCTION 25 (14.5%) · STEREO 6 (3.5%).**

`facmer_divergent` (15, 0.30 pts) is **100% `ARRANGEMENT_ONLY`** — metal and ligand bodies agree
and only the arrangement signature differs, which is exactly what the bucket name means. It needs
no new class and no new lane; it is the fac/mer axis and nothing else.

### Read examples — one per class, actually read

| class | molecule | what the strings show |
|---|---|---|
| `PERCEPTION` | `BOVCUM_comp_0` | `CC(=O)c1sc2ccccc2c1O` → `CC(O)=C1Sc2ccccc2C1=O` — a **ketone/enol tautomer shift**. Same atoms, same bonds, different bond *orders*. |
| `PERCEPTION` | `AHEBEV_comp_0` | `C#O` → `[C-]#[O+]` — carbon monoxide, neutral vs charge-separated. |
| `GEOM_CODE` | `DUNHEB_comp_0` | `Ni SPL` → `Ni OCT`, bodies **identical**. The coordination number was read as 6 where the input says 4. |
| `LIGAND_SPLIT` | `DAGNIL_comp_0` | 3 → 5 fragments: two benzyls each shed an H (`[CH2]c1ccccc1` → `C=C1C=C[CH]C=C1` + `[H]`). |
| `STEREO_INVERSION` | `DIVZOY_comp_0` | `C[P@@](…)` → `C[P@](…)`. **One inverted phosphorus.** Nothing else differs. |
| `SKELETON` | `MUYTEH_comp_0` | `C[C@]1(O)[CH]C(=O)CCC1` → `[CH]C(=O)CCCC(C)=O` — **the ring opened.** |
| `SKELETON` | `DIDROZ_comp_0` | a pyridine became a **fused bicyclic**. |
| `LIGAND_MERGE` | `XOYCOE_comp_0` | a chloride **fused into a ligand**. |
| `ARRANGEMENT_ONLY` | `ADANEB_comp_0` | metal same, bodies same, signature differs. |

## 🔴 Three corrections this lane makes

### 1. The charter's own framing of `BOUNDARY` is refuted by its control

The charter files the 62 `structural`→`BOUNDARY` as *"the attachment call itself is inside the
tolerance band"* — i.e. treats BOUNDARY as a **cause**. Over the 3858 molecules that round-trip
**perfectly**:

```
byte_exact BOUNDARY   1367 / 3858  (35.4%)
byte_exact INTACT     2431 / 3858  (63.0%)
```

**BOUNDARY is the modal state of a passing molecule.** Neither class discriminates, so neither
explains a failure. And the classification bears that out: BOUNDARY and INTACT decompose almost
identically (PERCEPTION 69.4% vs 71.8%). *A bucket name that asserts a cause is a hypothesis, not
a measurement* — the rule this release was chartered on, applied to the charter itself.

### 2. The first version of this instrument was wrong 57 times in 109, and said so confidently

The heavy-atom comparison was first written as string normalization — strip brackets, drop bond
symbols, uppercase. Checked against an RDKit canonical heavy-graph it **disagrees on 57 of 109**
molecules, a coin flip, because SMILES ring-closure digits and atom ordering are arbitrary labels
and an identical graph written two ways reads as different.

It reported `BEDLII_comp_0` as CONSTRUCTION. Its heavy graph is unchanged; the difference is
aromatic perception. **The broken version printed SKELETON 74/172 (43%) and the corrected one
prints 4/172 (2.3%)** — the two differ by a factor of 18 and *both* look like a finished
measurement. Caught only by validating the normalizer against a canonical comparison instead of
trusting it.

Two smaller instances of the same shape, both caught the same way:
- 23 bodies carry a `RAW:` **sentinel** the key builder prefixes when canonicalization failed.
  Parsed verbatim they read as "RDKit cannot parse this molecule". They were 13.4% of the
  population, reported as UNCLASSIFIED, and the truth was a missing four-character strip.
- Chiral tags survive `MolToSmiles`, so an **inverted stereocentre read as a broken graph**.
  `DIVZOY_comp_0` differs only as `[P@@]` vs `[P@]`. Comparing with and without stereo separates
  them, and STEREO_INVERSION is a distinct mechanism — the same family v0.4.17's enantiomer lane
  owns.

### 3. The reachability map's mechanism column needs amending for these 172

v0.4.15's map marks `structural`→INTACT and →BOUNDARY as **UNMEASURED** and the roadmap sizes
v0.4.17 as construction against the whole `structural` block. For these 172, that is now measured
and it is **mostly wrong**: 82% is perception, 14.5% construction.

⚠ **This does not move `byte_exact` by itself, and it is not a lane.** Perception being the
mechanism does not mean a fix is cheap — a tautomer shift on a generated geometry may be the
geometry's fault rather than the perceiver's, and this lane did not measure that. What it
establishes is **which mechanism to point the next lane at**, which is precisely what the block
lacked. The distinction between "reachable" and "reachable *by this mechanism*" is the one this
project has now got wrong in both directions in consecutive releases.

## Caveats, stated rather than buried

- **The classes are priority-ordered and a molecule can carry more than one.** Fragment count is
  tested before geometry code, so `DAGNIL_comp_0` (which is both `LIGAND_SPLIT` and `TBP`→`TPY`)
  is filed as `LIGAND_SPLIT`. The counts are therefore "the first mechanism that applies", not
  "the only one".
- **`PERCEPTION` is a mechanism, not a diagnosis.** It says the heavy-atom graph survived and the
  decoration did not. Whether the cause is the perceiver or a distorted geometry is unmeasured.
- **`METAL_ELEMENT` is 0 of 172** — the metal is never re-perceived as a different element. Kept
  as a class because a measured zero and a silent absence are different claims.
