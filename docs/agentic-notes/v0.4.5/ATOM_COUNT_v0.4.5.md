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

### Where that leaves the class

| | molecules | status |
|---|---|---|
| LOSS, cause B | 14 | **13/13 run end-to-end pass.** Default-on. |
| GAIN, cause A | 60 | **11/12 sampled pass end-to-end** behind the lever. Not extrapolated to the other 48. |
| of which not really defects | 3 | `GOFTUQ`, `INENOF`, `TESFIH` — the input is missing hydrogens (Sec 4b) |
| residual, open | 1 of the sample | `LOCGAL_comp_0` (+2) — diagnosed, left open as a notation question |

### Fix A — behind `OIN_H_FAITHFUL`, default OFF

`oin/hydrogen.py::h_faithful_smiles` writes a fragment, reads it back, compares hydrogen
counts atom by atom against RDKit's own `_smilesAtomOutputOrder`, and forces a bracket on any
atom that came back different. Applied at the three **encoder-side** serializations
(`oin_aligner`, `xyz2mol`, `oin/inline.py`), plus a narrower bracket-preserving step in the
adapter.

Default off because it changes OIN strings. Verified byte-identical with the lever unset, and
the four golden fixtures pass with it **both off and on**.

Measured, and this is the number that matters. Twelve GAIN molecules run end-to-end through
the real harness in one session, with the lever as the **only** difference:

| arm | result |
|---|---|
| `OIN_H_FAITHFUL` unset | **0 of 12 pass** — every one fails, 11 on the exact atom-count mismatch and `LOCGAL` on a 300 s timeout |
| `OIN_H_FAITHFUL=1` | **10 of 12 `status: success`** |
| `OIN_H_FAITHFUL=1` + the bare-P strip branch | **11 of 12** — `MEGZIH_comp_0` also passes, tier `UFF_1`, `oin1 == oin2` byte-identical |

`status: success` is the full contract, not an atom count: the canonical round-trip key
matched, the RMSD mapping succeeded, *and* the atom count agreed. Passing: CIDDAU, FABPEG,
HIKDIR, INENOF, JAPCOT, JOTJEK, KIKROO, KILQAZ, UQUXAG, XAJBIW.

### The two residuals, diagnosed

Both are genuine — the adapter-implied count is wrong for each in *both* arms (+2, +1), so
neither is an artifact of `LOCGAL`'s 300 s timeout masking its gate.

**`MEGZIH_comp_0` (+1) — a bare phosphorus donor had no strip branch. Now fixed.**
The adapter's bare-donor reconciliation had branches for C, O/S and N and **none for P**, so a
bare `P{n}` kept a phantom implicit H. Added one, on exactly the N branch's argument:
`replace_map` de-brackets a binding atom only when the bracket content is a bare
organic-subset symbol, so `[PH]`/`[PH2]` keeps its bracket and takes the explicit branch —
therefore a bare `P{n}` always means 0 H. Exact, not a heuristic.

The phantom only exists when perception gave the phosphorus valence 4 (a phosphaalkene `C=P`,
an ylide) and RDKit climbed to P's next allowed valence, 5, with one hydrogen. A tertiary
phosphine sits at valence 3 with 0 implicit H already, so PPh3/PMe3/dppe are untouched and the
many phosphine complexes that already round-trip cannot move. Measured: 0/100 regressions.

This one needs the lever **and** the branch together, and the reason is worth recording. With
the lever off the encoder writes `[PH]{1}` — the phantom has been **promoted into an explicit
bracket**, which the adapter is right to treat as authoritative, so nothing downstream can
recover it. With the lever on the encoder keeps `[P]`, `replace_map` renders it `P{1}`, and the
new branch strips it: **73 = input**.

**`LOCGAL_comp_0` (+2) — still open.** A Ru allenylidene, `C{2}=C=C(c1ccccc1)c1ccccc1`. The
donor carbon is bare with a double bond and **one** heavy neighbour, and the adapter's carbon
branch only strips a non-aromatic sigma-carbon at `heavy >= 2`, so it keeps two phantom
hydrogens. Not fixed here on purpose: a bare `C{n}=` with one heavy neighbour is genuinely
ambiguous between a 0-H carbene/vinylidene and a `=CH2` methylene ligand, so resolving it is
an OIN notation question like the nitride/ammine case, not a heuristic to guess at.

### A methodology note that cost real time

`tools/atomcount/adapter_scan.py` replays the OIN strings the **frozen sweep** recorded, so it
is blind to encoder-side improvements — it kept reporting `MEGZIH_comp_0` at +1 after the fix,
because the frozen string still carries the old encoder's `[PH]` phantom. Anything that
changes the encoder has to be measured with `encode_adapter_scan.py`, which re-encodes.

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

## 5c. Blast radius, measured

A change that can only affect molecules which currently **fail** cannot regress a passing
one -- but that has to be measured, not asserted. A passing round trip means the generated
structure had exactly as many atoms as the input, so for any passing molecule the
adapter-implied count computed with the current code must still equal its input count. The
frozen sweep's own `status: success` verdict is therefore the "before" arm, with no second
checkout needed.

100 currently-passing capstone molecules, seed 42 — **60 carrying the exact population Fix B
changes** (a 5-membered aromatic ring with a heteroatom: thiophene, pyrrole, furan, pyrazole,
imidazole — the rings it stops charging) and 40 controls with no such ring:

| arm | regressions |
|---|---|
| lever OFF (Fix B only) | **0 / 100** |
| lever ON (Fix B + Fix A's adapter half) | **0 / 100** |

Fix A's *encoder* half needs a re-encode to exercise, since the scan above reuses the stored
OIN strings; `tools/atomcount/passing_reencode_ab.py` does that in both arms and asserts the
built count still equals the input.

### Guards

- `tests/unit/test_regression_stability.py` (the four goldens): **green with the lever both
  off and on**.
- `tests/unit/test_haptic_carbon_hcount.py`, `tests/unit/test_bare_donor_hydrogens.py`,
  `tests/unit/test_atom_count_hydrogen.py`: 28 tests, **green in both arms**.
- Unit-suite baseline at `b23decb4`: `Ran 605 tests, OK (skipped=3, expected failures=3)`.
  Final code including the bare-P strip branch: **`Ran 611 tests, OK (skipped=3, expected
  failures=3)`, 0 failures** — the baseline plus the 6 new fixtures, nothing regressed.
- Key guards re-run in **both** lever arms after the P branch (four goldens, haptic carbon
  H-count, bare-donor hydrogens, atom-count fixtures, trivalent P): **31/31 OK in each arm**.
- Ruff `check` and `format`: clean.

## 5d. Two things the re-encode A/B turned up that are worth keeping

Re-encoding 25 currently-passing molecules in both lever arms produced **25/25
byte-identical strings** and **24/25 built counts still equal to the input**. The one
exception, `XOSTUW_comp_0`, is not a lever regression — both arms give 63 against an input
of 64 — but it is informative twice over.

### (a) A pre-existing boron regression, not from this lane

`XOSTUW_comp_0` passed the capstone at tier `UFF_1`. Its stored OIN wrote the boron as
`[BH]`; the current encoder writes bare `B`:

```
capstone : ...c1cn{2}n([BH](n2cccn{4}2)n2cccn{0}2)c1...   -> adapter builds 64 = input
current  : ...c1cn{2}n(B(n2cccn{4}2)n2cccn{0}2)c1...      -> adapter builds 63
```

Feeding the *stored* string through the **current** adapter still gives 64, so the adapter
and Fix B are not involved. The change is encoder-side and predates this lane: my only
string-affecting change is lever-gated, and both arms here are byte-identical. It therefore
sits somewhere between the capstone commit `58bba7ad` and current `main`.

The consequence matters for planning: **the frozen 74-molecule worklist understates this
class.** A fresh sweep would surface `atom_count` failures that are not in it.

### (b) Cause A also runs in the LOSS direction, and Fix A deliberately declines boron

This is the same write/read asymmetry, losing a hydrogen instead of gaining one. The
encoder's own molecule has the boron right — `exp=1, tot=1, val=4` — and the writer emits
bare `B`, which re-reads as 0 H.

`h_faithful_smiles` *detects* this correctly and then **refuses to repair it**, because the
repaired fragment cannot be verified: a neutral tetravalent boron fails RDKit's valence
check, so `[BH]`-bearing fragments do not sanitize at all.

```
c1cnn(B(n2cccn2)n2cccn2)c1      full sanitize = True    (B reads 0 H)
c1cnn([BH](n2cccn2)n2cccn2)c1   full sanitize = False   <-- cannot be verified
```

That is the verification gate working as designed rather than a bug — emitting a string that
does not sanitize would trade a wrong count for a broken fragment. But it means **boron is
outside Fix A's reach**, consistent with the boron ceiling already recorded for v0.4.4 SL5.
Closing it needs a decision about how OIN should represent tetravalent boron, which is a
notation question, not a serialization one.

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
