# Lane — Atom count / hydrogen (`OIN_H_FAITHFUL`)

**The failure population:** the `atom_count` hard-fail class — round trips that fail with
`Atom count mismatch at <tier>. Input N != Gen M`. **74 molecules (1.10%)** of the 6,719-molecule
v0.4.5 capstone corpus, the largest non-timeout `hard_fail` class and untouched before this
release. Re-measured in v0.4.6 as a **45-molecule** `Atom count mismatch` population on the
936-molecule re-baseline. It is **100% hydrogen** — and it is exactly the residue the round-trip
key is built to fold.

---

## ELI5

SMILES — the text format OIN is built on — has two ways to write an atom. Written bare, `C` means
*"carbon, and fill up the rest of its bonding capacity with hydrogen"*. Written in brackets,
`[C]` means *"carbon with exactly the hydrogens I said"*. RDKit's text writer drops the brackets
whenever it thinks they are unnecessary, so an atom the encoder knew had **zero** hydrogens gets
written bare and read back carrying one. That phantom hydrogen is then baked into the OIN string,
the 3D builder faithfully builds a molecule one atom too big, and the round trip fails its very
last check. The confusing part is that the *main* comparison — the round-trip key — deliberately
ignores hydrogen bookkeeping, so it had already declared the round trip a success. This class is
precisely where the two instruments disagree, by design, and the atom count is the honest one.

## The work, visually

```
   ┌───────────┐   the input crystal structure. GROUND TRUTH for this lane.
   │  input_H  │   read straight off the XYZ file (tools/atomcount/classify.py)
   └─────┬─────┘
         │
         │   perception:  xyz2AC_obabel → AC2BO → get_lig_mol → get_tmc_mol
         │   ┌──────────────────────────────────────────────────────────────┐
         │   │ HYPOTHESIS 1 said the divergence was HERE.   ✘ REFUTED       │
         │   │ perceived_H == input_H in 36 of 45                           │
         │   └──────────────────────────────────────────────────────────────┘
         ▼
   ┌─────────────┐  the encoder's own molecule — "intent". 0-H atoms are 0-H
   │ perceived_H │  ON PURPOSE, and the geometry agrees: 763 of 766 ambiguous
   └──────┬──────┘  sites are genuinely 0-H (727 planar sp2 + 36 sp)
          │
          │   ★★★  THE DIVERGENCE LIVES IN THIS SPAN  ★★★
          │        oin_H != input_H in 41 of 45.  dH spans −36 … +14.
          │
          │   serialization, three encoder-side writes + one emit:
          │   ┌─ utils/oin_aligner.py:378        h_faithful_smiles(...)   ── site 1
          │   ├─ utils/xyz2mol.py:1725           h_faithful_smiles(...)   ── site 2
          │   ├─ oin/inline.py:268               h_faithful_smiles(...)   ── site 3
          │   │     MolToSmiles(canonical=False) silently DE-BRACKETS `[C]` → `C`
          │   └─ oin/canonical_body.py:193, :278 h_faithful_smiles(...)   ── site 4
          │         reparse round trip; plain MolToSmiles here UNDID sites 1-3
          │         (fixed in v0.4.6 — both writes now route through h_faithful)
          ▼
   ┌────────┐   the emitted OIN string. What a consumer actually receives.
   │ oin_H  │   [Mn_SPY].C{0}#O.N#C[C@H]1[CH]{1}=[CH]{1>}C=[SH]{4}1.C{2}#O.C{3}#O
   └───┬────┘                                             ↑
       │                                    a thiophene S with NO hydrogen in the
       │                                    input; perceived valence 3 sits strictly
       │                                    between S's allowed 2 and 4, so RDKit
       │                                    climbs to 4 by adding one H
       │
       │   generation/metallogen_adapter.py::_prepare_ligand_fragments
       │   ├─ CAUSE B lived here: the kekulize rescue charged ring[0] of EVERY
       │   │  5-membered all-aromatic ring; a −1 on a BARE aromatic carbon flips
       │   │  its implicit H from 1 → 0 and destroys an innocent C–H
       │   └─ bare-donor strip heuristics recover SOME phantoms, per motif (C, O/S,
       │      N, and — added this lane — P)
       ▼
   ┌────────────┐  what MetalloGen is ASKED to build. Reproduces the harness's
   │ adapter_H  │  reported generated count in 59 of 59 auditable molecules.
   └─────┬──────┘  ⇒ the wrong count is fixed BEFORE any 3D work happens.
         │
         │   embed / assemble / optimize      ── GENERATOR FULLY EXONERATED
         ▼
   ┌───────────────┐
   │  generated_H  │──► harness gate order:
   └───────────────┘      1. canonical_roundtrip_key  ── says EQUAL (folds implicit-vs-explicit H)
                          2. RMSD mapping             ── passes
                          3. atom count  ◄── LAST GATE, tools/test_dataset_roundtrip.py:256-263
                                              the only instrument that can see this class
```

Legend — `★` = where the defect is. `✘` = a hypothesis killed by measurement. `dH` = `oin_H −
input_H`. "site N" = one of the four serializations `h_faithful_smiles` was wired into. The
`36/45`, `41/45`, `−36…+14` figures are the v0.4.6 re-baseline population; the `59/59`, `763/766`
figures are the v0.4.5 74-molecule cohort. **The two cohorts are different sets** — do not mix a
denominator from one with a numerator from the other.

## Initial assumptions and hypothesis

**Three claims in the lane charter, all overturned by the first measurement:**

1. **"Atoms went missing."** Assumed a deficit. **No** — 60 of the 74 *gain* atoms; only 14 lose
   them.
2. **"The generator built a smaller molecule than the one it was asked for."** Assumed a
   generation defect.
3. **"The harness correctly refused it."** Assumed the gate was doing straightforward work,
   rather than being the *only* instrument that could see a class the primary instrument folds by
   design.

**Two sub-hypotheses from the charter, both refuted the same way:** an uncoordinated fragment
being dropped, and a ligand truncated during m-SMILES assembly. Either would show a heavy-atom
deficit.

**One assumption that had to be settled before touching the encoder:** that the encoder's 0-H
reading is *correct*. If those inputs were simply missing their hydrogens — utterly ordinary in
CSD-derived data — then the generator adding them would be chemically **right**, the notation
would be fine, and there would be no defect to fix.

**Then, in v0.4.6, two further hypotheses of my own — both plausible, both wrong:**

- **H1:** "the divergence lives in perception (`get_lig_mol` / `AC2BO`)". Inferred from a flat A/B
  without checking *where* the hydrogen entered.
- **H2:** "one implicit H per severed metal–donor bond — a bare `N{1}` gains the hydrogen the metal
  bond was spending". Reads convincingly off examples.

## What was actually found

### Confirmed — the direction is GAIN, not loss

| dH (gen − input) | −2 | −1 | +1 | +2 | +3 | +4 | +5 | +6 | +7 | +8 | +10 | +12 | +13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| molecules | 4 | 10 | **30** | 13 | 6 | 2 | 2 | 1 | 1 | 2 | 1 | 1 | 1 |

**60 gain / 14 loss.**

### Confirmed — it is hydrogen, and only hydrogen

Of the 74, **27** have a stored generated structure (exactly the rows that failed at `FF_reroll_5`;
`UFF_1` rows save nothing). Element histogram, input vs generated:

- **27/27 differ in H and nothing else. `heavy_delta == 0` in every single one.**
- Net over the class: **H +28**, every other element **0**.

That single table refutes both charter sub-hypotheses: no dropped counter-ion, no truncated ligand,
no missing solvent fragment — either would show a heavy-atom deficit.

### Confirmed — the generator is exonerated completely

Running the recorded OIN string through the generator's own **front half** (`OINParser` →
`metallogen_adapter._prepare_ligand_fragments`) and counting the atoms MetalloGen is *asked* to
build reproduces the harness's reported generated atom count in **59 of 59** auditable molecules,
exactly. The wrong atom count is already fixed before any 3D work happens. Nothing is lost in
embedding, assembly or optimization.

### Confirmed — the encoder's 0-H reading is right, 763 times out of 766

`tools/atomcount/hybridization_probe.py` settles this from the input geometry alone, with no
chemistry model: a 3-coordinate 0-H carbon that is **planar** is sp2 (a missed double-bond
perception) and legitimately carries no hydrogen; the same carbon **pyramidal** is sp3, so the
hydrogen was never located and the input is incomplete; a 2-coordinate 0-H carbon at ~180° is
sp/nitrile (0 H correct), at ~109° is sp3 (two hydrogens absent).

Over all 74 molecules, every ambiguous 0-H carbon:

| verdict | sites | reading |
|---|---:|---|
| 3-coordinate **sp2**, out-of-plane ≤ 0.07 Å, angle sums 359.3–360.0° | **727** | 0 H is CORRECT |
| 2-coordinate **sp** (alkyne / nitrile) | **36** | 0 H is CORRECT |
| 2-coordinate **sp3**, angle ~112° | **3** | 2 H genuinely **missing from the input** |

**763 of 766.** Not one pyramidal 3-coordinate carbon in the entire class. So the *string* is what
is wrong, which is what makes this a genuine losslessness defect rather than a disagreement about
chemistry.

The three exceptions are worth naming because they are **not defects**: `GOFTUQ_comp_0` (C7,
neighbours Rh and O), `INENOF_comp_0` (C18, neighbours C and N) and `TESFIH_comp_0` (C14,
neighbours N and C) each have one sp3 carbon at ~112° with no hydrogen in the input file — a CH₂
whose hydrogens the crystal structure never located. For these the generator's extra hydrogens are
*chemically correct* and the input XYZ is the incomplete artifact (deltas +2, +2, +1, consistent
with exactly that). They can be made to **pass**, because encoding the input faithfully means
encoding a 0-H carbon and then the counts agree — but the molecule the notation then describes is a
carbene rather than a methylene. That is lossless with respect to the input and wrong with respect
to chemistry, and it should not be counted as a win.

### Confirmed — two disjoint causes, splitting exactly 60 / 14

All **74** audited. (The probe's first version reported 15 as "unparseable" purely because a
carbonyl written `[C]#O` is a carbon radical no sanitize accepts — a probe limitation, since fixed,
that would have read as 15 unexplained molecules.)

| cause | molecules | evidence |
|---|---:|---|
| **A** — bare-symbol write/read asymmetry | **60** — every GAIN row | adapter-implied count == reported generated count in 60/60; per-site traces close each delta exactly |
| **B** — kekulize rescue destroys a bare aromatic C–H | **14** — every LOSS row | a new aromatic `[c-]` in the adapter fragment that is absent from the OIN fragment, 14/14; with the fix the adapter-implied count equals the input in exactly those 14 |

The split is *exactly* the GAIN/LOSS split, which is the strongest evidence that the causes are
disjoint and nothing in the class is unaccounted for. GAIN rows where the rescue also fires are
still cause A: the rescue only ever *removes* hydrogen, so it can never source a net gain.

**Cause A, measured on `CIDDAU_comp_0` (18 → 19):**

```
encoder fragment SMILES : N#C[C@H]1C=CC=S1        intent H = 4   (matches the input exactly)
re-parsed               : N#C[C@H]1C=CC=[SH]1     read-back H = 5
OIN emitted             : [Mn_SPY].C{0}#O.N#C[C@H]1[CH]{1}=[CH]{1>}C=[SH]{4}1.C{2}#O.C{3}#O
```

That `[SH]{4}` is a thiophene sulfur bonded to Mn and two ring carbons and **no hydrogen at all**
in the input XYZ. Its perceived valence is 3, strictly between sulfur's allowed 2 and 4, so RDKit
adds one H to climb to 4. `SetNoImplicit(True)` does **not** force a bracket.

**Cause A bites a second, independent site** — `oin/inline.py`'s `MolToSmiles(mol, canonical=False)`
silently **de-brackets** a correctly-written 0-H atom. Measured on `INENOF_comp_0` (58 → 60):

```
encoder fragment  : CCCCN1[C]N([C]c2c3ccccc3cc3ccccc23)C=C1   intent H = 20
after inline.py   : CCCCN1CN(Cc2c3ccccc3cc3ccccc23)C=C1       read-back H = 24
```

Both `[C]` brackets are destroyed. The adapter's bare-donor strip heuristics recover 2 of the 4 —
the ones carrying a `{slot}` marker — and the non-binding benzylic `[C]` has no heuristic to save
it. Net +2, exactly the reported delta. The encoder's own intent, 58 atoms, is **exactly the
input**; the information is lost between the aligner and the string.

Note what this implies: the adapter's strip heuristics, and `oin_aligner.py` step 1b (the earlier
COLWIK / ACOXEX fix), are both **per-motif compensations for this same asymmetry**. This class is
the residue of the motifs nobody had enumerated yet.

**Cause B, measured on `QOBFOF_comp_0` (31 → 30).** `_prepare_ligand_fragments` rescued a fragment
that will not kekulize (the Cp-anion case) by putting `-1` on `ring[0]` of **every** 5-membered
all-aromatic ring in the fragment — including rings that kekulize perfectly well and are not the
reason the fragment failed. A `-1` on a **bare** aromatic carbon flips its implicit H from 1 to 0.

```
encoder intent    : C[C@]1(O)[C]=C(c2ccsc2)n2ccccc21          H = 11  (total 31 = input)
adapter prepared  : C[C@@]1(O)c2ccccn2C(c2[c-]csc2)=[C:1]1    H = 10
                                          ^^^^ thiophene C-H charged, H gone
```

A Cp ring is unaffected because its carbons are already `[cH]` brackets with `NoImplicit` set, so
the charge cannot change their H count; **only a bare aromatic carbon is vulnerable.**

### Confirmed — Fix B works end-to-end; Fix A works on the v0.4.5 sample

| | molecules | status |
|---|---|---|
| LOSS, cause B | 14 | adapter-implied count now equals the input for **all 14**; **13/13 run end-to-end come back `status: success`** at tier `UFF_1` |
| GAIN, cause A | 60 | **11 of 12 sampled pass end-to-end** behind the lever. **Not extrapolated to the other 48.** |
| of which not really defects | 3 | `GOFTUQ`, `INENOF`, `TESFIH` — the input is missing hydrogens |
| residual, open | 1 of the sample | `LOCGAL_comp_0` (+2) — diagnosed, left open as a notation question |

The 12-molecule GAIN A/B, run through the **real harness** in one session with the lever as the
**only** difference:

| arm | result |
|---|---|
| `OIN_H_FAITHFUL` unset | **0 of 12 pass** — every one fails, 11 on the exact atom-count mismatch and `LOCGAL` on a 300 s timeout |
| `OIN_H_FAITHFUL=1` | **10 of 12 `status: success`** |
| `OIN_H_FAITHFUL=1` + the bare-P strip branch | **11 of 12** — `MEGZIH_comp_0` also passes, tier `UFF_1`, `oin1 == oin2` byte-identical |

`status: success` is the full contract, not an atom count: the canonical round-trip key matched,
the RMSD mapping succeeded, *and* the atom count agreed. Passing: CIDDAU, FABPEG, HIKDIR, INENOF,
JAPCOT, JOTJEK, KIKROO, KILQAZ, UQUXAG, XAJBIW (+ MEGZIH with the P branch).

### ★ Confirmed in v0.4.6 — promoting the lever buys NOTHING measurable

The v0.4.5 blocker was real and it was fixed: `canonical_body_emit` had two plain
`Chem.MolToSmiles` writes — the intermediate that feeds the reparse, and the **final emit whose
output becomes the body** — so `xyz2mol.py:1725` computed an H-faithful string and the canonical
body then overwrote it. That is exactly how `OIN_CANONICAL_BODY` (default-**ON**) "undid"
`OIN_H_FAITHFUL`. Both writes now route through `h_faithful_smiles`
(`oin/canonical_body.py:193` and `:278`), verified **byte-identical on all 61 fixtures**.

**And it changed nothing.** A/B over the 45-molecule `Atom count mismatch` population of the
936-molecule re-baseline, comparing the OIN-implied atom count against the input XYZ
(generator-free):

| arm | match | mismatch |
|---|---|---|
| `OIN_H_FAITHFUL=0`, canonical body ON | 8 | 37 |
| `OIN_H_FAITHFUL=1`, canonical body ON | **8** | **37** |

Identical. `h_faithful_smiles` guarantees only that a string re-reads with the hydrogens it was
*written* with, so whatever moves the count on this population is not that. **The reason the lever
stays off therefore changed** — from "a default-ON lever undoes it" to "it does not move a real
population" — and that is the better reason, because the first one was a bug and the second one is
a fact about the class.

### Confirmed — where the count actually diverges

| measurement | result | reading |
|---|---|---|
| `perceived_H` (parent mol) vs `input_H` (XYZ) | **agrees in 36 / 45** | perception is right |
| `oin_H` (implied by the emitted string) vs `input_H` | **diverges in 41 / 45** | the string is wrong |
| `dH` vs count of **bare donor atoms** in the string | **matches in only 4 / 45** | no donor-cut rule |
| `dH` range | **−36 … +14**, bidirectional | not one mechanism |

So the divergence sits **between the perceived parent and the emitted string** — not in perception,
and not in write/read fidelity.

**What the distribution says instead: the class is HETEROGENEOUS, with at least two mechanisms.**

- **28 of 45** sit at `dH` = **+1…+3** — consistent with something small and donor-adjacent;
- **4** are `dH` = **0** — the mismatch is **not hydrogen at all** for these;
- **3** are large **losses** (**−14, −16, −36**) that no donor-cut story explains — most likely
  haptic/eta bodies or the `[CH]`-radical writing, but that is a guess, not a measurement.

## What was done

### Fix B — default ON, no lever

`src/oinsmiles/generation/metallogen_adapter.py::_prepare_ligand_fragments`. The kekulization
rescue now charges only **all-carbon** 5-membered aromatic rings:

```python
if len(ring) != 5:                                   continue
if not all(a.GetIsAromatic() for a in atoms):        continue
if not all(a.GetAtomicNum() == 6 for a in atoms):    continue   # ← the fix
mol.GetAtomWithIdx(ring[0]).SetFormalCharge(-1)
```

A ring of five carbons cannot kekulize neutral and the `-1` is what makes it a legal aromatic
anion; a thiophene, pyrrole, furan or pyrazole kekulizes as it stands and never needed charging.
**Rings that need the charge keep exactly the old behaviour**, so no Cp / indenyl / fluorenyl path
moves.

Default-on because it cannot change any OIN string — it is entirely on the generator's
SMILES-preparation side — and because the only behaviour it removes is charging a ring that had no
reason to be charged.

**Deliberately NOT done:** preserving the H count of the atom that *does* get charged. On an
all-carbon haptic ring that charge is load-bearing in a second, previously undocumented way — it is
what strips the phantom implicit H off a bare 0-H eta ipso/fusion carbon written `c{n}`. An eta
indenyl depends on it (`tests/unit/test_haptic_carbon_hcount.py`).

### Fix A — `OIN_H_FAITHFUL`, default OFF

`src/oinsmiles/oin/hydrogen.py::h_faithful_smiles(mol, **kwargs) -> str`. It is
`Chem.MolToSmiles(mol, **kwargs)` with a verification loop bolted on:

1. write the SMILES;
2. read it back (`_reparse`: full sanitize, then retry `SANITIZE_ALL ^ SANITIZE_KEKULIZE`, which is
   the one step a metal-stripped aromatic donor ring reliably fails);
3. compare hydrogen counts **atom by atom** against RDKit's own `_smilesAtomOutputOrder`
   (`_divergent_atoms` — `order[canonical_position] = original atom index`, which makes the
   comparison exact rather than a heuristic match);
4. every atom that came back different gets `SetNumExplicitHs(intended)`, `SetNoImplicit(True)`,
   and **an unpaired electron** — purely to force RDKit to BRACKET it. SMILES has no radical
   syntax, so **nothing about the radical reaches the string**;
5. re-serialize and re-check.

It returns the plain `MolToSmiles` output **unchanged** when the string was already faithful, when
the lever is off, when the fragment cannot be re-parsed, or when the repair fails to verify — so it
can only ever change a string that was **measurably wrong**. That bounds the blast radius to
molecules whose atom count the notation was already getting wrong. `mol` is never mutated
(`_cached_copy` works on a copy, guarded by `tests/…::test_never_mutates_the_caller_mol`).

Two guards inside `_repair` that are load-bearing and easy to remove by accident:

- **The atom ORDER must not move.** `if _output_order(candidate_mol) != _output_order(mol): return
  smiles`. Callers read `_smilesAtomOutputOrder` off the molecule they passed in —
  `xyz2mol.get_oin_string` uses it to decide which character position each `{slot}` marker attaches
  to — and that property records the FIRST serialization. An unpaired electron can in principle
  change the canonical ranking, and if it did, **every slot marker would land on the wrong atom**.
  Coverage lost by declining is visible as a still-failing molecule; a mis-slotted string would not
  be visible at all.
- **Never raise.** `h_faithful_smiles` wraps `_repair` in `except Exception: return smiles`, because
  `oin/inline.py` wraps *its* caller in a bare `except Exception` that falls back to a **different
  slot-tagging strategy** — an exception escaping here would not surface as an error, it would
  silently reroute the encoder.

**Wired into four writes** — three encoder-side plus the emit: `utils/oin_aligner.py:378`,
`utils/xyz2mol.py:1725`, `oin/inline.py:268`, and (v0.4.6) `oin/canonical_body.py:193` + `:278`.
Plus a narrower bracket-preserving step in the adapter.

**Default off because it changes OIN strings.** Verified byte-identical with the lever unset, and
the four golden fixtures pass with it **both off and on**.

### The bare-phosphorus strip branch — default ON, in the adapter

`metallogen_adapter.py:279`. The adapter's bare-donor reconciliation had branches for C, O/S and N
and **none for P**, so a bare `P{n}` kept a phantom implicit H. Added one, on exactly the N
branch's argument: `replace_map` in `oin/inline.py` de-brackets a binding atom only when the
bracket content is a bare organic-subset symbol, so `[PH]` / `[PH2]` keeps its bracket and takes
the explicit branch — therefore **a bare `P{n}` always means 0 H**. Exact, not a heuristic.

The phantom only exists when perception gave the phosphorus valence 4 (a phosphaalkene `C=P`, an
ylide) and RDKit climbed to P's next allowed valence, 5, with one hydrogen. A tertiary phosphine
sits at valence 3 with 0 implicit H already, so PPh₃ / PMe₃ / dppe are untouched and the many
phosphine complexes that already round-trip cannot move. Measured: **0/100 regressions.**

This one needs the lever **and** the branch together, and the reason is worth recording: with the
lever off the encoder writes `[PH]{1}` — the phantom has been **promoted into an explicit bracket**,
which the adapter is right to treat as authoritative, so nothing downstream can recover it. With
the lever on the encoder keeps `[P]`, `replace_map` renders it `P{1}`, and the new branch strips it:
**73 = input**.

### What was deliberately NOT done

**`h_faithful_smiles` is not used at the adapter's final serialization.** It makes a string faithful
to the molecule it is handed, and *that* molecule's H counts are not ground truth — a bare `c{n}`
haptic ipso carbon arrives already carrying a phantom implicit H, so "faithful" would preserve the
error being removed. It broke the eta indenyl and every LOSS fixture when tried. Ground truth lives
**upstream, in the encoder**, where counts come from the input geometry.

## Dead ends and refutations

### "Atoms went missing" — REFUTED

**Killed by the delta histogram:** 60 of 74 GAIN, only 14 lose. Every subsequent piece of reasoning
that assumed a deficit was reasoning about 19% of the class.

### "An uncoordinated fragment is being dropped" / "a ligand is truncated during m-SMILES assembly" — BOTH REFUTED

**Killed by the element histogram over the 27 molecules with a stored generated structure:** 27/27
differ in H and nothing else, `heavy_delta == 0` in every one, net H +28 and every other element 0.
Either mechanism would show a heavy-atom deficit.

### "The generator built a smaller molecule than it was asked for" — REFUTED

**Killed by `tools/atomcount/adapter_scan.py`:** the atoms MetalloGen is *asked* to build reproduce
the harness's reported generated count in **59 of 59** auditable molecules, exactly. The wrong count
is already in the string.

### "The inputs are just missing crystallographic hydrogens, so the generator is right" — REFUTED, 763/766

The hypothesis that would have dissolved the whole lane. **Killed by
`tools/atomcount/hybridization_probe.py`:** 727 planar sp2 + 36 sp = 763 of 766 ambiguous 0-H
carbons are genuinely 0-H, and **not one pyramidal 3-coordinate carbon appears in the entire
class**. The 3 exceptions are named and excluded from the win count.

### ⚠ H1: "the divergence lives in perception (`get_lig_mol` / `AC2BO`)" — REFUTED

Inferred in v0.4.6 from the flat A/B without checking *where* the hydrogen entered — the classic
shape of an over-read. **Killed by `perceived_H == input_H` in 36 of 45.** Perception is right.

⚠ **This refutation has not propagated:** `docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md`'s own "Revised ranking"
still lists item 2 as *"Perception-side hydrogen — the atom-count class, **now correctly located in
`get_lig_mol` / `AC2BO`** rather than in serialization"*, which is the hypothesis §1 of the same
document refutes. The ranking list is the earlier, stale version. Treat §1 as authoritative.

### ⚠ H2: "one implicit H per severed metal–donor bond" — REFUTED

The story: a bare `N{1}` gains the hydrogen the metal bond was spending, because stripping the slot
marker turns a coordinated donor into a free ligand. It is a **real effect** and it reads
convincingly off examples like `CUDBOU` (`N{1}c1ccccc1N{2}`). It does not survive the corpus.
**Killed twice over:** `dH` equals the bare-donor count in only **4 of 45**, and `dH` is
**bidirectional, spanning −36 to +14** — no donor-cut rule can produce a negative `dH`.

### "Patching the encoder-side serializations is enough" — REFUTED, and this one nearly shipped

Patching the three encoder-side writes alone moved **45 of the 74 OIN strings and changed the built
atom count for exactly 0 of them** — the adapter's own final `MolToSmiles` re-de-brackets
everything the string had preserved. Without that A/B this would have shipped as "45 strings fixed"
with **zero molecules passing**.

### A methodology dead end that cost real time

`tools/atomcount/adapter_scan.py` replays the OIN strings the **frozen sweep** recorded, so it is
blind to encoder-side improvements — it kept reporting `MEGZIH_comp_0` at +1 *after* the fix,
because the frozen string still carries the old encoder's `[PH]` phantom. **Anything that changes
the encoder must be measured with `tools/atomcount/encode_adapter_scan.py`, which re-encodes.**

### Boron is outside Fix A's reach — a declined repair, not a bug

Cause A also runs in the LOSS direction. On `XOSTUW_comp_0` the encoder's own molecule has the boron
right (`exp=1, tot=1, val=4`) and the writer emits bare `B`, which re-reads as 0 H.
`h_faithful_smiles` **detects this correctly and then refuses to repair it**, because the repaired
fragment cannot be verified:

```
c1cnn(B(n2cccn2)n2cccn2)c1      full sanitize = True    (B reads 0 H)
c1cnn([BH](n2cccn2)n2cccn2)c1   full sanitize = False   ← cannot be verified
```

A neutral tetravalent boron fails RDKit's valence check, so `[BH]`-bearing fragments do not sanitize
at all. That is the verification gate working as designed — emitting a string that does not sanitize
would trade a wrong count for a broken fragment. Closing it needs a decision about how OIN should
represent tetravalent boron: a notation question, not a serialization one. (Note the cross-lane
tension: `OIN_BORON_CAGE` is now default-ON and relaxes precisely this valence rule for *cage*
fragments. Whether `h_faithful_smiles`'s verification should route through
`aromaticity.py::sanitize_allowing_boron_cage` has **not** been measured.)

### `LOCGAL_comp_0` (+2) — deliberately left open

A Ru allenylidene, `C{2}=C=C(c1ccccc1)c1ccccc1`. The donor carbon is bare with a double bond and
**one** heavy neighbour, and the adapter's carbon branch only strips a non-aromatic sigma-carbon at
`heavy >= 2`, so it keeps two phantom hydrogens. Not fixed on purpose: a bare `C{n}=` with one heavy
neighbour is genuinely ambiguous between a 0-H carbene/vinylidene and a `=CH2` methylene ligand, so
resolving it is an OIN **notation** question like the nitride/ammine case, not a heuristic to guess
at.

## Where it landed

| change | lever | default | file |
|---|---|---|---|
| Fix B — kekulize rescue narrowed to all-carbon rings | none | **ON** | `generation/metallogen_adapter.py` |
| Bare-phosphorus strip branch | none | **ON** | `generation/metallogen_adapter.py:279` |
| Fix A — H-faithful serialization | `OIN_H_FAITHFUL` | **OFF** (held) | `oin/hydrogen.py` + 4 write sites |
| H-faithful canonical body (v0.4.6) | inherits `OIN_H_FAITHFUL` | inert while it is off | `oin/canonical_body.py:193`, `:278` |

**`OIN_H_FAITHFUL` is still held off, and the recorded reason CHANGED in v0.4.6.** The current text
lives in `src/oinsmiles/oin/levers.py::_HELD_OFF["OIN_H_FAITHFUL"]`: the `OIN_CANONICAL_BODY`
interaction is **fixed**; it stays off because promoting it **buys nothing measurable** (8/37 in
both arms over the 45-molecule population), because the class is heterogeneous, and because two
aggregate hypotheses have already been refuted. *"Promote only with evidence that it moves a real
population."*

**How it landed.** Branch `swimlane/v045-atomcount`, tip **`c1ae2759`**, forked from `b23decb4`.
Unlike most lanes it has **no individual `Merge branch 'swimlane/v045-atomcount'` commit** — it
entered `main` directly through the release integration merge **`1450b5ce`**
(`release(v0.4.5): integrate 16 lanes and PROMOTE the six canonicality levers`), which is also the
first merge on `main` that contains it. The v0.4.5 release commit is `0d165845`. Note that
`swimlane/v045-boron` forked from the *same* base `b23decb4` and is a **sibling**, not a
descendant: `c1ae2759` is not an ancestor of the boron merge `f4c3525a`.

**Lane commits on `swimlane/v045-atomcount`, oldest first:**

| commit | what |
|---|---|
| `ddf731ea` | `docs(atomcount): the atom_count class is 100% hydrogen, and the generator is innocent` |
| `bc0c04fc` | `fix(atomcount): stop the kekulization rescue from destroying an innocent aromatic C-H` |
| `125259e2` | `fix(atomcount): OIN_H_FAITHFUL -- make a serialized fragment re-read with the H it was written with` |
| `7cc89f80` | `docs(atomcount): geometry confirms the encoder's 0-H reading, 763 sites out of 766` |
| `2db4371d` | `fix(atomcount): narrow the kekulize rescue to all-carbon rings; carry brackets through the adapter` |
| `63eb85a3` | `docs(atomcount): audit all 74; the two causes split exactly 60/14` |
| `6d9a0b8d` | `docs(atomcount): record what each fix is scoped to, and the A/B that kept Fix A honest` |
| `aee77e33` | `test(atomcount): prove neither fix regresses a currently-passing molecule` |
| `899012a9` | `test(atomcount): re-encode A/B for Fix A's encoder-side blast radius` |
| `d18c4cb5` | `docs(atomcount): record the measured blast radius and the guard status` |
| `3c6e2cce` | `docs(atomcount): a pre-existing boron regression, and why Fix A declines boron` |
| `59ab8a9e` | `docs(atomcount): the end-to-end A/B -- 0/12 pass with the lever off, 10/12 with it on` |
| `f9ae6595` | `fix(atomcount): add the missing bare-phosphorus strip branch; diagnose both residuals` |
| `2be1b8a5` | `docs(atomcount): MEGZIH confirmed end-to-end; sampled GAIN result is 11/12` |
| `c1ae2759` | `docs(atomcount): final guard status -- 611 tests OK with the P branch, 0 failures` |
| `d799de1f` (v0.4.6) | H-faithful canonical body + the revised `_HELD_OFF` reason |

**Blast radius, measured — not asserted.** A change that can only affect molecules which currently
*fail* cannot regress a passing one, but that has to be measured. A passing round trip means the
generated structure had exactly as many atoms as the input, so for any passing molecule the
adapter-implied count computed with the current code must still equal its input count — which makes
the frozen sweep's own `status: success` verdict the "before" arm, with no second checkout needed.

100 currently-passing capstone molecules, seed 42 — **60 carrying the exact population Fix B
changes** (a 5-membered aromatic ring with a heteroatom: thiophene, pyrrole, furan, pyrazole,
imidazole — the rings it stops charging) and 40 controls with no such ring:

| arm | regressions |
|---|---|
| lever OFF (Fix B only) | **0 / 100** |
| lever ON (Fix B + Fix A's adapter half) | **0 / 100** |

Fix A's *encoder* half needs a re-encode to exercise, since that scan reuses stored OIN strings;
`tools/atomcount/passing_reencode_ab.py` does it in both arms. Re-encoding 25 currently-passing
molecules gave **25/25 byte-identical strings** and **24/25 built counts still equal to the input**
— the one exception, `XOSTUW_comp_0`, is 63-vs-64 in **both** arms and is therefore not a lever
regression (see the boron dead end, and the open question below).

**Guard tests:**

| test | pins |
|---|---|
| `tests/unit/test_atom_count_hydrogen.py::TestKekulizeRescueKeepsHydrogen::test_adapter_preserves_input_atom_count` | Fix B on 5 real capstone rows: `QOBFOF_comp_0` (31), `AJODEI_comp_0` (97), `MUXKOH_comp_0` (84), `DUJPIJ_comp_0` (63), `EGUDAL_comp_0` (81) |
| `…::TestKekulizeRescueKeepsHydrogen::test_cp_ring_still_gets_its_charge_and_keeps_its_hydrogens` | the rescue still works for the case it exists for — `UQUXAG_comp_0`'s unsubstituted Cp keeps its `-1` **and** all five C–H (`_atom_total == 10`) |
| `…::TestHFaithfulSmiles::test_lever_off_is_byte_identical_to_moltosmiles` | lever off ⇒ exactly `MolToSmiles` |
| `…::TestHFaithfulSmiles::test_lever_on_makes_the_string_reread_with_the_same_hydrogen` | **asserts the premise too** — `assertNotEqual` on the plain string's H count, so if RDKit ever stops drifting, the suite says so instead of the lever silently becoming pointless |
| `…::TestHFaithfulSmiles::test_never_mutates_the_caller_mol` | radical/explicit-H/no-implicit state of the caller's mol unchanged |
| `…::TestHFaithfulSmiles::test_never_raises_on_a_mol_with_no_property_cache` | `oin/inline.py` passes a `sanitize=False` parse; `GetTotalNumHs()` there *raises* rather than returning a wrong answer |
| `tests/unit/test_haptic_carbon_hcount.py`, `tests/unit/test_bare_donor_hydrogens.py` | the eta indenyl and the per-motif strip branches; green in **both** lever arms |
| `tests/unit/test_regression_stability.py` (4 goldens) | green with the lever both off and on |

**Suite:** baseline at `b23decb4` `Ran 605 tests, OK (skipped=3, expected failures=3)`; final code
including the bare-P branch **`Ran 611 tests, OK (skipped=3, expected failures=3)`, 0 failures** —
the baseline plus the 6 new fixtures, nothing regressed. Key guards re-run in **both** arms after
the P branch (four goldens, haptic carbon H-count, bare-donor hydrogens, atom-count fixtures,
trivalent P): **31/31 OK in each arm.** Ruff `check` + `format`: clean.

## Open questions / for the next agent

### 1. The next step is PER-ATOM ATTRIBUTION, and precisely this

Aggregates have now produced **two plausible-and-wrong answers** (H1 and H2 above). Do not run a
third aggregate. Instead:

> **Walk one molecule from each `dH` band, and for each hydrogen that changes, record the full
> provenance chain: parent atom index → fragment atom index → emitted token, plus WHICH STEP the
> count changed at.**

Concretely — the bands, and a candidate from each:

| band | n | pick | why this band is its own question |
|---|---:|---|---|
| `dH` = 0 | **4** | any | the mismatch is **not hydrogen**. Whatever it is has never been named, and it cannot share a mechanism with the rest. |
| `dH` = +1…+3 | **28** | `CUDBOU` (the H2 exemplar) | the bulk. "Small and donor-adjacent" is the *shape* of the effect, not the mechanism. |
| `dH` = −14, −16, −36 | **3** | the −36 | no donor-cut story explains a large loss. Suspected haptic/eta bodies or the `[CH]`-radical writing — **suspected, not measured.** |

The instrument already exists and is the right one: `tools/atomcount/stage_trace.py <MOLECULE>`
walks a single molecule through all five stages (input XYZ → encoder intent `kmol` → emitted
fragment SMILES re-parsed → the OIN string's fragments after `oin/inline.py` de-bracketing → the
adapter's prepared fragment). `tools/atomcount/frag_h_probe.py` splits explicit-vs-implicit H at
the moment the encoder freezes it, which is how you tell an input-derived H from a phantom one.
Use `encode_adapter_scan.py`, never `adapter_scan.py`, for anything encoder-side.

### 2. `hydrogen_faithfulness_enabled` was never migrated to the lever registry

`src/oinsmiles/oin/hydrogen.py:49` still reads `bool(os.environ.get("OIN_H_FAITHFUL"))` — the
**bare-truthiness spelling**, which means **`OIN_H_FAITHFUL=0` ENABLES the lever**. This is the
exact trap `oin/levers.py` was created to close (its module docstring names it), and
`OIN_BORON_CAGE` had five sites on the same spelling before migration. `_HELD_OFF` documents
`OIN_H_FAITHFUL`, but `lever_enabled` is never called for it. Low risk while the lever is off by
default and nobody sets it to `0` — and a live foot-gun for anyone trying to *disable* it
explicitly, e.g. in a test that wants both arms.

### 3. A stale comment says the v0.4.6 fix has not happened

`src/oinsmiles/utils/xyz2mol.py` (~line 1741) still reads *"Do NOT promote `OIN_H_FAITHFUL` until
`canonical_body_emit` is H-faithful too"*. It **is** H-faithful, as of `d799de1f`. The comment's
advice happens to remain correct for a different reason, which is exactly how a stale comment
survives review.

### 4. The frozen 74-molecule worklist UNDERSTATES this class

`XOSTUW_comp_0` passed the capstone at tier `UFF_1`; its stored OIN wrote the boron as `[BH]` and
the current encoder writes bare `B`, taking the adapter-built count from 64 (= input) to 63.
Feeding the *stored* string through the *current* adapter still gives 64, so neither the adapter nor
Fix B is involved — the change is encoder-side and predates this lane, somewhere between the
capstone commit `58bba7ad` and current `main`. **A fresh sweep will surface `atom_count` failures
that are not on the 74-molecule list.** The 45-molecule v0.4.6 population is the newer measurement;
the two are different cohorts and their numbers must not be mixed.

### 5. Two notation decisions are blocking, not deferrable-forever

- **`LOCGAL_comp_0`**: does a bare `C{n}=` with one heavy neighbour mean a 0-H carbene/vinylidene or
  a `=CH2` methylene? OIN currently cannot say.
- **Tetravalent boron**: `[BH]`-bearing fragments do not sanitize, so `h_faithful_smiles` correctly
  declines to repair them. Now that `OIN_BORON_CAGE` is default-ON and relaxes exactly this valence
  rule for cage fragments, it is worth measuring whether the H-faithful verification should route
  through `aromaticity.py::sanitize_allowing_boron_cage`. **Not measured.**

### 6. The three non-defects should stay excluded from any win count

`GOFTUQ`, `INENOF`, `TESFIH` can be made to pass, but the notation then describes a carbene where
chemistry has a methylene. Lossless with respect to the input, wrong with respect to chemistry.
Counting them is how a lane manufactures 3 free wins.
