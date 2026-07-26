# The `atom_count` hard-fail class — v0.4.5

74 molecules (1.10% of the 6,719-molecule capstone corpus) fail the round trip with
`Atom count mismatch at <tier>. Input N != Gen M`. This is the largest non-timeout
`hard_fail` class and had not been looked at this release.

Everything below is measured. The reproduction commands are in §6.

---

## 1. Three claims in the charter that the measurement overturns

**"Atoms went missing."** No. **60 of the 74 GAIN atoms; only 14 lose them.**

| delta (gen − input) | −2 | −1 | +1 | +2 | +3 | +4 | +5 | +6 | +7 | +8 | +10 | +12 | +13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| molecules | 4 | 10 | **30** | 13 | 6 | 2 | 2 | 1 | 1 | 2 | 1 | 1 | 1 |

**"The generator built a smaller molecule than the one it was asked for."** No — the
generator is **exonerated completely**. Running the recorded OIN string through the
generator's own front half (`OINParser` → `metallogen_adapter._prepare_ligand_fragments`)
and counting the atoms MetalloGen is *asked* to build reproduces the harness's reported
generated atom count in **59 of 59** auditable molecules, exactly. The wrong atom count is
already fixed before any 3D work happens. Nothing is lost in embedding, assembly, or
optimization.

**"The harness correctly refused it."** Only in the narrowest sense. The atom-count check is
the **last** gate in `tools/test_dataset_roundtrip.py:255-263` — reached only *after*
`canonical_roundtrip_key` has already declared the round trip **equal**. That key
deliberately folds implicit-vs-explicit hydrogen (`oin/compare.py:441`). So this class is
precisely the residue the lossy key folds: the two instruments disagree, by construction,
about hydrogen — and the atom count is the honest one.

## 2. It is hydrogen. Only hydrogen.

Of the 74, 27 have a stored generated structure (exactly the rows that failed at
`FF_reroll_5`; `UFF_1` rows save nothing). Element histogram, input vs generated:

- **27/27 differ in H and nothing else. `heavy_delta == 0` in every single one.**
- Net over the class: **H +28**, every other element **0**.

No dropped counter-ion, no truncated ligand, no missing solvent fragment. The charter's
"uncoordinated fragment being dropped" and "a ligand truncated during m-SMILES assembly"
hypotheses are both **refuted** by this: either would show a heavy-atom deficit.

## 3. Root cause A — the SMILES writer's bare-symbol asymmetry (the GAIN rows)

RDKit's SMILES writer emits a **bare** organic-subset symbol whenever it judges brackets
unnecessary. A bare symbol is read back as *"fill to the next allowed valence with
hydrogen."* `SetNoImplicit(True)` does **not** force a bracket. So an atom the encoder
believes carries 0 H serializes bare and re-parses with 1 (or 2) H — a phantom hydrogen,
baked into the OIN string, which the generator then faithfully builds.

Measured, `CIDDAU_comp_0` (18 → 19):

```
encoder fragment SMILES : N#C[C@H]1C=CC=S1        intent H = 4   (matches the input exactly)
re-parsed               : N#C[C@H]1C=CC=[SH]1     read-back H = 5
OIN emitted             : [Mn_SPY].C{0}#O.N#C[C@H]1[CH]{1}=[CH]{1>}C=[SH]{4}1.C{2}#O.C{3}#O
```

That `[SH]{4}` is a thiophene sulfur bonded to Mn and two ring carbons and **no hydrogen at
all** in the input XYZ. Its perceived valence is 3, which sits strictly *between* sulfur's
allowed valences 2 and 4, so RDKit adds one H to climb to 4.

The same writer asymmetry bites a **second, independent site**, `oin/inline.py`'s
`MolToSmiles(mol, canonical=False)`, which silently **de-brackets** a correctly-written
0-H atom. Measured, `INENOF_comp_0` (58 → 60):

```
encoder fragment  : CCCCN1[C]N([C]c2c3ccccc3cc3ccccc23)C=C1   intent H = 20
after inline.py   : CCCCN1CN(Cc2c3ccccc3cc3ccccc23)C=C1       read-back H = 24
```

Both `[C]` brackets are destroyed. The adapter's bare-donor strip heuristics
(`metallogen_adapter.py:190-250`) recover 2 of the 4 — the ones carrying a `{slot}` marker —
and the non-binding benzylic `[C]` has no heuristic to save it. Net +2, exactly the reported
delta. The encoder's own intent, 58 atoms, is **exactly the input**; the information is lost
between the aligner and the string.

Note what this means: those adapter strip heuristics, and `oin_aligner.py:311-350`
(step 1b, the earlier COLWIK/ACOXEX fix), are both **compensations for this asymmetry**,
applied per-motif. This class is the residue of the motifs not yet enumerated.

## 4. Root cause B — the kekulization rescue destroys an innocent C–H (the LOSS rows)

`metallogen_adapter._prepare_ligand_fragments` rescues a fragment that will not kekulize
(the Cp-anion case) like this:

```python
except Exception:
    for ring in ring_info.AtomRings():
        if len(ring) == 5 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            mol.GetAtomWithIdx(ring[0]).SetFormalCharge(-1)
```

It puts −1 on `ring[0]` of **every** 5-membered all-aromatic ring in the fragment — including
rings that kekulize perfectly well and are not the reason the fragment failed. And a −1 on a
**bare** aromatic carbon flips its implicit H count from 1 to 0. The hydrogen is destroyed
silently.

Measured, `QOBFOF_comp_0` (31 → 30). Its ligand carries a thiophene that is perfectly
kekulizable; the fragment fails only because of the odd-valence carbene donor:

```
encoder intent    : C[C@]1(O)[C]=C(c2ccsc2)n2ccccc21          H = 11  (total 31 = input)
adapter prepared  : C[C@@]1(O)c2ccccn2C(c2[c-]csc2)=[C:1]1    H = 10
                                          ^^^^ thiophene C-H charged, H gone
```

**This explains 13 of the 13 auditable LOSS rows — 100%.** A Cp ring is unaffected because
its carbons are already `[cH]` brackets (`NoImplicit` set), so the charge cannot change their
H count; only a bare aromatic carbon is vulnerable.

## 4b. Is the encoder's 0-H reading correct? Geometry says yes, 763 times out of 766

Cause A rests on the encoder recording carbons and heteroatoms with **0 H**, and the whole
class turns on whether that reading is right. If those inputs were simply missing their
hydrogens -- ordinary in CSD-derived data -- then the generator adding them would be
chemically *correct*, the notation would be fine, and there would be no defect to fix. That
had to be settled before touching the encoder, and the input geometry settles it without any
chemistry model (`tools/atomcount/hybridization_probe.py`):

* a 3-coordinate 0-H carbon that is **planar** is sp2 -- a double bond perception missed --
  and legitimately carries no hydrogen;
* the same carbon **pyramidal** is sp3, so the hydrogen was never located and the input is
  incomplete;
* a 2-coordinate 0-H carbon at ~180 deg is sp/nitrile (0 H correct), at ~109 deg is sp3
  (two hydrogens absent).

Over all 74 molecules, every ambiguous 0-H carbon:

| verdict | sites | reading |
|---|---|---|
| 3-coordinate **sp2**, out-of-plane <= 0.07 A, angle sums 359.3-360.0 deg | **727** | 0 H is CORRECT |
| 2-coordinate **sp** (alkyne / nitrile) | **36** | 0 H is CORRECT |
| 2-coordinate **sp3**, angle ~112 deg | **3** | 2 H genuinely **missing from the input** |

**763 of 766.** Not one pyramidal 3-coordinate carbon in the entire class. The encoder's 0-H
reading is right and the *string* is what is wrong, which is what makes cause A a genuine
losslessness defect rather than a disagreement about chemistry.

### The three that are not defects

`GOFTUQ_comp_0` (C7, neighbours Rh and O), `INENOF_comp_0` (C18, neighbours C and N) and
`TESFIH_comp_0` (C14, neighbours N and C) each have one sp3 carbon at ~112 deg with **no
hydrogen in the input file** -- a CH2 whose hydrogens the crystal structure never located.
For these the generator's extra hydrogens are *chemically correct* and the input XYZ is the
incomplete artifact. Their deltas (+2, +2, +1) are consistent with exactly that.

They can still be made to *pass*, because encoding the input faithfully means encoding a 0-H
carbon, and then the counts agree -- but the molecule the notation then describes is a carbene
rather than a methylene. That is lossless with respect to the input and wrong with respect to
chemistry, and it is worth saying so plainly rather than counting them as wins.

## 5. Root-cause histogram over the 74

All **74** are now audited (the first version of the probe reported 15 as "unparseable"
purely because a carbonyl written `[C]#O` is a carbon radical no sanitize accepts — a probe
limitation, since fixed, that would have read as 15 unexplained molecules).

| cause | molecules | evidence |
|---|---|---|
| B — kekulize rescue destroys a bare aromatic C–H | **14** — every LOSS row | a new aromatic `[c-]` in the adapter fragment that is absent from the OIN fragment, in 14/14; and with the fix the adapter-implied count equals the input in exactly those 14 |
| A — bare-symbol write/read asymmetry | **60** — every GAIN row | adapter-implied count == reported generated count in 60/60; per-site traces close each delta exactly |

The split is exactly the GAIN/LOSS split measured in §1 — 60 and 14 — which is the strongest
evidence that the two causes are disjoint and that nothing in the class is unaccounted for.
The GAIN rows where the rescue also fires are still cause A: the rescue only ever *removes*
hydrogen, so it can never be the source of a net gain.

**Neither cause is a generator defect.** Cause B lives in the generator's front-end SMILES
preparation but is pure bookkeeping — it never reaches the 3D machinery. Cause A is an
encoder losslessness defect: the OIN string does not record the input's atom count.

## 5b. The two fixes, and what each is scoped to

### Fix B — default ON, no lever

`generation/metallogen_adapter.py`: the kekulization rescue now charges only **all-carbon**
5-membered aromatic rings. A ring of five carbons cannot kekulize neutral and the `-1` is
what makes it a legal aromatic anion; a thiophene, pyrrole, furan or pyrazole kekulizes as it
stands and never needed charging. Rings that need the charge keep **exactly** the old
behaviour, so no Cp/indenyl/fluorenyl path moves.

It is default-on because it cannot change any OIN string — it is entirely on the generator's
SMILES-preparation side — and because the only behaviour it removes is charging a ring that
had no reason to be charged.

Measured: the adapter-implied atom count now equals the input for **all 14 LOSS molecules**,
and **13/13 of the ones run end-to-end come back `status: success`** at tier `UFF_1` — the
canonical key, the RMSD mapping and the atom count all passing, not merely the count.

### Fix A — behind `OIN_H_FAITHFUL`, default OFF

`oin/hydrogen.py::h_faithful_smiles` writes a fragment, reads it back, compares hydrogen
counts atom by atom against RDKit's own `_smilesAtomOutputOrder`, and forces a bracket on any
atom that came back different. Applied at the three **encoder-side** serializations
(`oin_aligner`, `xyz2mol`, `oin/inline.py`), plus a narrower bracket-preserving step in the
adapter.

Default off because it changes OIN strings. Verified byte-identical with the lever unset, and
the four golden fixtures pass with it **both off and on**.

Measured, and this is the number that matters: with the lever the only difference in one
session, the same molecules go from **4/4 failing on the exact atom-count mismatch** to
**4/4 `status: success`** (CIDDAU, FABPEG, HIKDIR, INENOF).

### The measurement that saved the fix from being cosmetic

Patching the three encoder-side serializations alone moved **45 of the 74 OIN strings and
changed the built atom count for exactly 0 of them.** The adapter's own final `MolToSmiles`
re-de-brackets everything the string had preserved. Without that A/B this would have shipped
as "45 strings fixed" with zero molecules passing.

And one thing deliberately *not* done: `h_faithful_smiles` is **not** used at the adapter's
final serialization. It makes a string faithful to the molecule it is handed, and that
molecule's H counts are not ground truth — a bare `c{n}` haptic ipso carbon arrives already
carrying a phantom implicit H, so "faithful" preserves the error being removed. It broke the
eta indenyl and every LOSS fixture. Ground truth lives upstream, in the encoder, where counts
come from the input geometry.

## 6. Reproduction

```bash
cd /home/tjmustard/Documents/GitHub/oin-v045-atomcount
export PYTHONPATH=$PWD/src
V=/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python
W=/home/tjmustard/Documents/GitHub/OIN-SMILES/spec/handoffs/v0.4.5/hard_fail_worklists.json

$V tools/atomcount/classify.py                       # element histograms, direction, deltas
$V tools/atomcount/adapter_scan.py --worklist $W     # authoritative counts + cause B flag
$V tools/atomcount/hfaithful_scan.py --worklist $W   # cause A, per-atom write/read drift
$V tools/atomcount/stage_trace.py QOBFOF_comp_0      # one molecule through all five stages
```
