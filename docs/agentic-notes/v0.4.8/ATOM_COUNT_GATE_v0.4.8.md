# v0.4.8 · Lane 2 — the `Atom count mismatch` gate

**Verdict: (a). The gate is right, it is load-bearing, and it is not a third error direction.**

For 8 of the 18 molecules it fires on, the OIN string of the *generated* structure is
**byte-identical** to the OIN of the input — while the generated structure carries one or two
extra hydrogens. No string comparison can separate them: not the scored one, not the honest
one introduced by Lane 1, not the canonical key. The atom-count gate is the only instrument in
the harness that can see it, and what it sees is a genuine violation of the losslessness
contract this project exists to satisfy.

**Lane 1 needs no correction from this lane.** All 18 are `hard_fail` in both scoring arms
(`status != success` with an error that is not `String mismatch`), so they never entered the
`byte_exact` / `key_equal` counts the honest baseline reports.

---

## 1. What was asked, and what the charter got wrong

The v0.4.8 charter framed a binary: either **(a)** the gate is right and the key folds
something real, or **(b)** the gate is a third error direction failing 27 correct structures.
It also stated the gate runs *before* the string comparison.

Two corrections, both from the code and the corpus rather than from the 936-molecule cohort
the charter was written against:

| charter | measured |
|---|---|
| the gate runs **before** the string comparison | it runs **after** — key comparison at `tools/test_dataset_roundtrip.py:264`, RMSD, then atom count at `:319`. Every `Atom count mismatch` row is a molecule whose notation **already passed**. |
| population is **27** | at N = 5000 it is **18** (with 10 `String mismatch`; 18 + 10 = exactly the 28 failures that produced a stored structure). The 27 was the 936-molecule rebaseline cohort. |
| the sharpest probe is the **ΔH 0** cases | **there are none at corpus scale.** All 18 are ΔH `+1` or `+2`. The 4 ΔH-0 cases were an artifact of the smaller cohort. |

Reproduce:

```bash
V=$PWD/.venv/bin/python; export PYTHONPATH=$PWD/src
$V tools/atom_count_provenance.py \
    --results-dir tmCAT-tmPHOTO_xyz_dataset/results-v0.4.6-sweep \
    --json-out /path/to/atom_count_provenance.json
```

## 2. The population, per atom

18 molecules. **The element delta is hydrogen only in 18/18** — always gained, never lost, and
no non-hydrogen atom is involved in any of them.

| n | ΔH |
|---:|---|
| 10 | +1 |
| 8 | +2 |

Per-atom provenance is computed from geometry alone: for every heavy atom, the signature
`(element, #H attached, sorted heavy-neighbour elements)`, diffed as a multiset between the
two files. Multiset diffs need no atom correspondence, which matters because the generator does
not preserve heavy-atom order. Each hydrogen is assigned to its **nearest** heavy atom rather
than to everything within a cutoff — a bridging or agostic H would otherwise be double-counted
and manufacture a difference that is not there.

**9 of 18 protonate an atom in place**: the same element with the same heavy neighbours
reappears carrying one more hydrogen. The clean minimal cases are `BUWHAD_comp_1`
(`C(C,C,Cu)` 0H → 1H — a metal-bound carbanion picking up a proton), `ODUJEC_comp_0`
(`C(C,C,Pd)` 0H → 1H), `TOXDIV_comp_0` (`N(N)` 0H → 1H) and `XAJBIW_comp_0` / `XENNIO_comp_0`
(`C(C,C,X)` 0H → 1H).

The other 9 move the heavy-atom skeleton as well, and 3 of those are outright ligand
detachment — `UQUXAG_comp_0` takes Fe from 8 contacts to **0**, a whole ferrocene sandwich
coming apart, and the freed carbons pick up hydrogens where the metal bond used to be.

## 3. The decisive measurement

For each molecule, re-encode the *generated* structure independently — a full
`XYZToSMILES().convert()` of its coordinates — and compare to the input's OIN.

⚠ This comparison must use the independent re-encode, **never** the sweep's `smiles_2`.
`smiles_2` comes from `get_oin_string(gen_result.mol, coords)`, the generator's own bond graph,
so "the strings match" there is circular and proves nothing about the geometry. That
circularity is precisely what Lane 1 exists to correct; reusing it here would have produced a
confident wrong answer.

| | n |
|---|---:|
| independent re-encode **still equals the input OIN** — invisible to every string comparison | **8** |
| independent re-encode differs — the honest arm catches it too | 10 |
| `coordination.py` flags a lost metal contact | 3 |

The 8 that are string-invisible:
`ALEMOT_comp_0`, `DOCPAO_comp_0`, `MIBFEL_comp_0`, `NEFNER_comp_0`, `NOYTUS_comp_0`,
`ULOQIX_comp_0`, `UPABUK_comp_0`, `XAKCAP_comp_0`.

### Why this settles (a) vs (b)

**(b) required that the 18 be correct round trips.** They are not: a structure with two extra
hydrogens is not the structure that was encoded. The notation is simply not injective over
hydrogen count for this class, so its silence is not evidence of agreement.

**The trap that would have produced the opposite answer.** The charter's own reading of (b)
rested on "27 molecules key-match under independent re-perception". They do — and it means
nothing, because the key folds hydrogen count. *A lossy key must never be reused as an
acceptance predicate for an axis it folds.* Checking what the key actually folds, before
concluding anything from `indep_key_match`, is what turns this from a plausible (b) into a
measured (a).

## 4. `XAKCAP_comp_0` — the case that defeats three instruments

```
[Mo_OCT].CN(C)CS{0}.S{1}=c1ccccn{3}1.C{2}#O.C{4}#O.c1ccc(P{5}(c2ccccc2)c2ccccc2)cc1
```

61 atoms in, 63 atoms out, `+2 H`, and:

| instrument | verdict |
|---|---|
| scored string (`gen_result.mol`) | **byte-exact pass** |
| honest string (independent re-perception) | **byte-exact pass** — the same string |
| canonical key | **equal** |
| `coordination.py` | **intact** — this is not a metal-sphere defect |
| **atom count** | **FAIL, 61 != 63** |

Both added hydrogens sit on a heavy atom (nearest-neighbour distance < 1.3 Å), so this is
protonation, not stray unbonded atoms.

Pinned by `tests/unit/test_atom_count_gate.py`, with `KIKSAB_comp_0` pinned alongside it as a
case the honest arm *does* catch — so the lane's conclusion cannot be misread as either "the
gate catches everything" or "the honest metric makes the gate redundant". Neither is true.

## 5. Where the hydrogen comes from, and what is left open

`OIN_H_FAITHFUL` was built for this class and buys **nothing**: A/B over the 45-molecule
population reads match 8 / mismatch 37 with the lever off and the identical 8 / 37 with it on.
That A/B did establish the divergence sits between the *perceived parent* and the *emitted
string* — `perceived_H == input_H` in 36/45 — so perception is not the suspect
(`docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md`).

This lane adds the per-atom half: the atoms that gain hydrogen are **donor atoms whose
metal-bond the SMILES valence model does not represent as a bond**. A carbanion carbon bound to
Pd, an amide nitrogen bound to Ru, a haptic ring carbon — each is written without an explicit
bond to the metal, so RDKit fills the free valence with hydrogen on the way back out. That is
consistent with all 18 being `+1`/`+2` and with the protonated atoms being exactly the
metal-adjacent ones.

**Not fixed here, deliberately.** The fix is an encoder change — making the emitted string
pin hydrogen count on metal-adjacent donors, in the spirit of the existing donor-bracket
convention. Landing an encoder change in the release that re-baselines the accuracy number
would make the before/after unreadable, which is the entire cost v0.4.8 is paying to avoid.
The 18 are a worklist, not this lane's deliverable.

**No code changed in this lane.** The gate is correct as written; the deliverable is the
measurement, the tooling that reproduces it, and the tests that pin it.

## 6. Artifacts

| what | where |
|---|---|
| the tool | `tools/atom_count_provenance.py` |
| the per-molecule data | `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/atom_count_provenance.json` |
| the tests | `tests/unit/test_atom_count_gate.py` |
| the fixtures | `tests/fixtures/atom_count_gate/{XAKCAP,KIKSAB}_comp_0_{input,generated}.xyz` |
