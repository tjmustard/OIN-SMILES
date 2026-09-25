# v0.4.19 E2 lane — perception is order-dependent in three places, found by bisecting ONE molecule

> Branch `research/v0419-e2` (worktree `../oin-v0419-e2`, off `research/v0419-serializer` @ `1f61681c`).
> **Nothing promoted.** Three levers, all `_HELD_OFF`: `OIN_N_VALENCE_2`, `OIN_CANONICAL_CHARGES`,
> `OIN_CANONICAL_RESONANCE`. Test: `tests/unit/test_e2_canonical_perception.py` (fixture
> `tests/fixtures/XIVMEX_comp_0.xyz`). Results: `results-v0.4.19-e2/` (dataset tree, gitignored).

## 0. The question

The census of the v0.4.18 release sweep attributes **105 molecules (2.10 pts) to `E2_P_FRAGILE`**: the
input string changes under renumbering or 0.02 Å noise, so the round trip compares a string that is
not a function of the molecule. The scoping pass (`spec/handoffs/v0.4.19/NEXT.md`) had refuted three
hypotheses read-only — the valence-search truncation (1/105), the canonical ranking (the canonically
permuted AC is identical under random renumbering 105/105) and the `suppress_canonical_perception`
retry (1/60) — and concluded the order-dependence enters DOWNSTREAM of `AC2BO`. This lane bisected
`get_tmc_mol` stage by stage on one molecule with two numberings, then a third and a fourth.

## 1. Instrument: one numbering-invariant summary per stage

`e2_bisect.py` (session scratch; the logic is in the test) runs the base file and the audit's own
`renum0/1/2` variants (same seeds as `tools/census/e_selfconsistency.py`) through twelve stages and
prints, for each, a value that cannot depend on the numbering unless the stage does:

| stage | value compared |
|---|---|
| S1 `get_basic_mol` | canonical SMILES of the all-single-bond mol |
| S2 donors | sorted Z of the metal's AC neighbours |
| S3 `MetalDisconnector` → fragments | canonical SMILES per ligand, sorted |
| S4 `get_proposed_ligand_charge` | the Hückel charge |
| S5 `AC2BO` at q, q±2 | bond-order multiset **and** canonical SMILES of the bare BO graph |
| S6 `AC2mol` | canonical SMILES (bond orders + formal charges) |
| S7 `lig_checks` | the sorted `_resonance_candidate_key` list |
| S8/S9 `_select_lig_mol` / `get_lig_mol` | chosen form + final charge |
| S10/S11 `_get_tmc_mol_impl` / `get_tmc_mol` | canonical SMILES of the TMC mol |
| S12 | the OIN string |

The first stage whose value differs is where the order enters. Each numbering ran in its own
process (the `AC2BO` memo is process-global; and `hash(bytes)` is salted per process — a first
version compared those and printed a false divergence).

## 2. XIVMEX_comp_0 — Zn 5,15-ditolyl-porphyrin, census `E2_P_FRAGILE`, 63 atoms

Shipped (all three levers off): base string is the correct aromatic dianion; `renum0/1/2` give
three DIFFERENT strings, all with `[CH][CH]` radicals; `noise1` gives a quinoid.

**Bisection, shipped:** S1–S5 identical — the bond-order GRAPH `AC2BO` returns is the same under both
numberings (canonical perception works) — **S6 differs**: `AC2mol` charges the same BO differently.

### 2a. Why the BO is a zwitterion at all: nitrogen has no valence 2

```
valence options per (Z, AC valence):  H (1,)  C(3) (4,)  C(4) (4,)  N(2) (3, 4)
combo size 16   q=-2: doubles=0  (16 candidates, none valid)   q=0: doubles=18, valid at candidate 1
```

`atomic_valence[7] = [3, 4]`. A two-coordinate nitrogen can be an imine (3) or an ammonium (4) but
never an **N⁻** — pyrrolide, amide, the two bare pyrrole nitrogens of a porphyrin dianion. At
charge −2 every one of the 16 assignments fails `charge_is_OK`, and `AC2BO` returns `best_BO`,
which is initialised to **the AC itself: zero double bonds**. `set_atomic_charges` then walks that
graph: every ring carbon has BO-valence 3 → −1, flipped to +1 by the running-total rule
(`if n_single == 3 and q + 1 < mol_charge: q += 2; charge = 1`), the four N are −1 → a 33-carbon
zwitterion whose `[c+]/[c-]` pattern is whatever the walk order made it. RDKit then perceives the
charged single-bond rings as aromatic, and `lig_checks`' `ResonanceMolSupplier` enumerates
resonance forms of THAT: from the file's numbering the enumeration reaches the neutral-carbon
aromatic dianion (41 candidates, best key `(-34, 0, 0)`); from `renum0` the supplier returns
**one** form, the zwitterion itself (`(-32, 32, 0)`), and the string ships with radicals.

So the shipped "correct" base string was a lottery ticket. **The census class is Ni 24 / Zn 22
N-macrocycles because those are the anionic-N ligands.**

### 2b. Leak 1 — `OIN_CANONICAL_CHARGES`: the charge walk is in input order

`charge_is_OK` (inside `AC2BO`, run in the canonical labelling under `OIN_CANONICAL_PERCEPTION`)
and `set_atomic_charges` (inside `BO2mol`, run on the un-permuted BO) are "the same loop, same
branch ladder, same order-dependent carbon corrections" — the first decides in canonical order,
the second places in **input** order. Same BO, two zwitterions. The lever passes `AC2BO`'s
permutation (`_canonical_perception_perm`, factored out so the two walks cannot disagree on
whether the canonical path was taken) to `set_atomic_charges` as its walk order.

Alone it makes XIVMEX **stable and wrong**: all four numberings now agree on a quinoid with
`[CH][CH]` — the canonical walk's zwitterion is one the enumeration cannot neutralise. A
reproducible wrong answer; the lever is a stability lever, not a chemistry lever.

### 2c. Leak 0 — `OIN_N_VALENCE_2`: let the search write an N⁻

`possible_valences` (and its mirror in `_valence_search_is_truncated`) append 2 for a nitrogen of
AC-valence ≤ 2, **last**, so a ligand the search already solves keeps its first valid candidate.
XIVMEX at −2 is then solved at the first assignment with two N at 2 (perfect matching, 17 double
bonds): S5/S6 identical across numberings, no radicals from any numbering.

Divergence moves to **S7**: the sorted candidate lists differ.

### 2d. Leak 2 — `OIN_CANONICAL_RESONANCE`: the supplier returns a numbering-dependent subset

Same `AC2mol` output (canonical SMILES identical), `ResonanceMolSupplier`: **167 forms from the
file's order, 131 from `renum0`; 82 vs 94 distinct canonical forms; 16 only in one set, 28 only
in the other**; `maxStructs` is not the cause (uncapped: still 167). `lig_checks` sorts what it
was given and `_select_lig_mol` keeps the first, so the Kekulé form still moves.

First fix tried: `Chem.RenumberAtoms(lig, canonical order)` — **not enough** (`renum1/2` still
differ): the supplier also walks the BONDS by index and `RenumberAtoms` keeps them in insertion
order. Fix: re-parse the ligand's canonical SMILES (`removeHs=False`), which creates atoms AND
bonds in string order, enumerate on that, renumber each form back with the inverse of
`_smilesAtomOutputOrder` before anything reads an index. Forms carry no props; the caller already
restores `__origIdx` by index.

### 2e. Result on XIVMEX

| configuration | base | renum0 | renum1 | renum2 |
|---|---|---|---|---|
| shipped | dianion ✓ | radicals | radicals | quinoid+radicals |
| CC only | quinoid | = base | = base | = base |
| N2 only | dianion ✓ | dianion, other Kekulé | = renum0 | = renum0 |
| N2 + RES | dianion ✓ | = base | = base | = base |
| N2 + RES + CC | dianion ✓ | = base | = base | = base |

`noise1` still moves under every configuration: the extended-Hückel charge proposal is a
coordinate-sensitive threshold (an orbital crossing −10 eV), and with noise it proposes 0, which
the search solves as the quinoid at candidate 1. Not addressed here (the census says 100/105
break under RENUMBERING; noise is the minority axis).

## 3. Measurement

### 3a. The 105 fragile rows × 11 columns under N2 + RES, first pass — an instrument bug

`results-v0.4.19-e2/esc105_e2b_v1_noconformer.jsonl` (unit `oin-v0419-e2-esc-e2b`, `--cpu 2` so the
concurrent candidate sweep's load is barely touched; `tools/v0419/e2_fragility_report.py`):

```
molecules 105  arm encoded 84  shipped encoded 105
fragile (shipped, arm): fixed 68   broken 0   still 16   (21 not encoded)
  renum  fixed 72  broken 0  still 7        noise  fixed 1  broken 0  still 11
base string moved vs sweep smiles_1: 56 of 84
arm encode errors: 'Cannot normalize a zero length vector' 20, 'bgnIdx not connected' 1
```

The 21 failures were the lever's own bug: the re-parsed canonical ligand has **no conformer**, the
supplier's forms are copies of what it is given, and `CIPAssigner.assign_all` reads 3D from the
ligand — every atom at the origin. Fixed by carrying the conformer across in the new order
(`bd4c73a3`, pinned by `test_frame_carries_the_conformer`); ALITEU encodes. Of the 84 that did
encode: **68 stable, 0 newly fragile, renumbering fixed on 72 of 79** — the mechanism holds on the
class, not just on XIVMEX. The 11 still noise-fragile are the Hückel-threshold axis (§2e).

### 3b. Second pass (in progress)

`esc105_e2b.jsonl` (unit `oin-v0419-e2-esc-e2b2`, same arm, conformer fix).
- Then, after the candidate sweep releases the machine: the whole-cohort audit (5,000 × 11) under
  the chosen arm, the changed-string set, the offline re-score of the unchanged rows, the live
  A/B over the changed set with parse-back — the serializer lane's pipeline
  (`tools/v0419/run_serializer_ab.sh` gets a new arm).

## 4. Process note — the sweep I contaminated

The first launch of the v0.4.19 candidate sweep ran from `../oin-v0419` while I edited `src/` in
that worktree for this bisection. The harness spawns a fresh interpreter per molecule, so every
molecule after 19:00 imported the dirty tree. Aborted at 528/5000, parked as
`results-v0.4.19-candidate-sweep.ABORTED-dirty-src/`, relaunched clean at 19:08 from `1f61681c`;
this lane moved to its own worktree. A worktree with a detached job is read-only.
