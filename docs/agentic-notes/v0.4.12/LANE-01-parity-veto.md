# Lane 1 — `OIN_FOLD_PARITY_VETO`: the filter v0.4.11 specified

**Status: BUILT, default-OFF.** Measurement tables below carry the command that produced them.

v0.4.11 built the within-fragment donor fold, measured **+7.86 `byte_exact` points across 393
molecules in one direction with the comparison key untouched**, and then found it **collapses
enantiomers in 221 of those same 393 gains**. Its close-out named one next step: *filter the
swap set by reflection parity, then re-run the uniform mirror audit before quoting points.*
This lane is that filter.

---

## 1. Both of v0.4.11's anchors reproduce exactly on today's tree

Run **before** any code was written, because a lane whose premise has moved is not worth
building.

```bash
PYTHONPATH=$PWD/src .venv/bin/python tools/fold_transition_sim.py \
    --sweep tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest --arm fold
PYTHONPATH=$PWD/src .venv/bin/python tools/mirror_audit_donor_fold.py \
    --dataset tmCAT-tmPHOTO_xyz_dataset/cat --n 250 --seed 7
```

| anchor | v0.4.11 | today | |
|---|---|---|---|
| `key_equal/slot_renumber → byte_exact` | 393 | **393** | ✓ |
| `byte_exact` | 3623 → 4016, +7.86 pts | **3623 → 4016, +7.86 pts** | ✓ |
| `facmer_divergent` | 16 → 16 | **16 → 16** | ✓ |
| moved in any other direction | 0 | **0** | ✓ |
| mirror audit `REGRESSION_raw_collapsed` | 19 / 250 (7.6%) | **19 / 250** | ✓ |
| `achiral_or_preexisting_fold` / `distinct_both_arms` / `encode_failed` | 157 / 73 / 1 | **157 / 73 / 1** | ✓ |

v0.4.11 measured the +393 with an **uncommitted ad-hoc script**, so the release's headline
number could not be re-derived without rewriting it. It is now `tools/fold_transition_sim.py`.

---

## 2. The obvious fix is wrong, and it is wrong *before* it is built

The tempting repair — require each swap to be a proper rotation of the polyhedron, reusing
`canonical_slots.derive_rotation_group`'s `det > 0` test — fails on inspection:

> A donor swap is a **transposition fixing every other vertex**. For a rank-3 polyhedron such
> a permutation does not preserve the Gram matrix at all, so that test rejects **every** swap
> and the fold degenerates to the identity.

The fold is not justified by a symmetry of the polyhedron. It is justified by a symmetry of the
**ligand**, realized as a proper rotation of the whole complex — which is not a property of the
coordination graph and cannot be read off it. Hence:

> ### 🔴 Reflection parity is not a property of the emitted string.
> **The fix cannot live where the defect is visible.**

`tests/unit/test_canonical_slots.py::TestDonorFoldCollapsesEnantiomers` was expected by
v0.4.11's close-out to be *"inverted, not deleted"*. It is **not inverted**: `canonicalize_oin_slots`
still collapses its two hand-written strings and always will. Its docstring now records why,
and the coordinate-level proof lives in `tests/unit/test_fold_parity.py`.

---

## 3. What was built

The veto sits in `get_oin_string`, where the **pristine** conformer is still in hand, and lifts
`tools/mirror_audit_donor_fold.py`'s own implication into the encoder:

```
S_rot   = canonicalize(inline,        fold OFF)    S_rot_m  = canonicalize(mirror, fold OFF)
S_fold  = canonicalize(inline,        fold ON )    S_fold_m = canonicalize(mirror, fold ON )

veto  <=>  (S_rot != S_rot_m)  and  (S_fold == S_fold_m)
emit   =   S_rot if veto else S_fold
```

**The left conjunct is load-bearing.** Without it every *achiral* molecule is vetoed — its
mirror encodes identically with the fold on **or** off — and so is every metal-centred Δ/Λ pair
whose descriptor the shipped encoder already folds (`OIN_EMIT_METAL_CONFIG` is held off). That
is v0.4.16's gap, and this lever must neither be blamed for it nor allowed to hide behind it.

**Cost: one extra encode**, and only where the fold actually fires. `S_rot == S_fold`
short-circuits before the mirror is ever built; the other three strings are pure string
operations. It reads the conformer captured **before** `_align_to_pai`, which can itself
reflect — mirroring an already-reflected conformer composes two improper maps into a proper
one, and the veto would never fire while looking perfectly healthy.

---

## 4. 🔴 Two silent defects in this lever's own construction

Both are recorded because the way they were caught is the transferable part.

### 4.1 `tmc_mol`'s atom order is not the coordinate order

Zipping `GetAtoms()` against `xyz_coords` positionally encoded `BIWDIV_comp_0` as **`[Co_TBP]`
with invented bond orders** instead of `[Co_OCT]` — a chemically different molecule that
encodes perfectly cleanly and whose mirror comparison is therefore fiction. Perception stamps
`__origIdx` for exactly this reason and `get_oin_string` itself reads `xyz_coords[orig_i]`.

Fixed via `__origIdx`, **plus a runtime self-check**: re-encode the reconstruction *unmirrored*
and require it to reproduce the labeling being decided about, else decline to fold.

### 4.2 The mirror was encoded with the fold inherited, and every test still passed

The mirror encode ran under the ambient (ON) fold, so `mirror_oin` was the **folded** mirror
string. That made `s_rot_m` identical to `s_fold_m`, which reduced the left conjunct to *"did
the fold fire?"* — true by construction at that point — and **disarmed the achiral guard**.
The self-check then tripped on **18 of 18** movers.

> **All three fixture tests passed anyway.** Declining to fold *also* separates a mirror pair,
> so a dead veto and a working one are indistinguishable from the output string alone.

This is the project's recurring shape — v0.4.7's attachment check was a silent no-op whose
complete 21-molecule A/B "reported what a genuine null result looks like". The fix is
structural, not a one-off correction:

`fold_parity.resolve` now records **why** it decided, and the five outcomes are kept distinct:

| outcome | meaning |
|---|---|
| `fold_inactive` | the fold never fired; nothing to police, nothing paid |
| `declined_no_conformer` / `declined_no_pairs` / `declined_no_self_encode` / `declined_reconstruction_drift` / `declined_no_mirror` | **no instrument** — emits the rotation-only labeling for lack of evidence |
| `vetoed_collapse` | the evidence says the fold destroys a mirror pair |
| `allowed_preexisting_fold` | the shipped encoder already folds this pair — not this lever's doing |
| `allowed_separation_survives` | the fold fired and the enantiomers stay distinct |

Three of these emit the same string. Only the outcome tells them apart, and
`test_it_separates_them_BY_VETOING_and_not_by_giving_up` now asserts on it.

### 4.3 The same trap, one layer out, in the measuring tool

`tools/fold_transition_sim.py`'s first draft implemented its `veto` arm by setting
`OIN_FOLD_PARITY_VETO=1` and re-canonicalizing **strings**. Since the veto lives one level up,
that arm was byte-identical to the `fold` arm and would have reported a clean, plausible and
completely false *"the veto costs nothing"*. The shipped arm re-encodes from **coordinates**
and carries a drift control.

---

## 5. Outcome distribution — the veto is doing the work, not declining

40 of the fold's own movers, input side (`tools/fold_transition_sim.py` movers list):

| outcome | n |
|---|---:|
| `fold_inactive` | 12 |
| `allowed_separation_survives` | 12 |
| `vetoed_collapse` | 11 |
| `allowed_preexisting_fold` | 5 |
| **`declined_*`** | **0** |

Zero declines is the number that matters: the self-check passes corpus-wide, so the veto proper
— not the fallback — is deciding every case.

**A side effect worth its own line:** the drift control re-encodes the frozen corpus's inputs
and its stored generated structures with today's encoder and requires them to reproduce the
v0.4.8 strings. They do. So v0.4.9 / v0.4.10 / v0.4.11's *"default path byte-identical"* claims
are now **confirmed rather than inherited**, and the carry-forward licence rests on measurement.

---

## 6. 🔴 The gate: 19 → 0 collapses, and nothing else moved

```bash
OIN_FOLD_PARITY_VETO=1 PYTHONPATH=$PWD/src .venv/bin/python \
    tools/mirror_audit_donor_fold.py --dataset tmCAT-tmPHOTO_xyz_dataset/cat --n 250 --seed 7
```

| verdict | baseline | **veto** | |
|---|---:|---:|---|
| 🔴 `REGRESSION_raw_collapsed` | 19 | **0** | the gate |
| `distinct_both_arms` | 73 | **92** | ← the same 19, now correctly separated |
| `achiral_or_preexisting_fold` | 157 | **157** | the achiral guard did **not** over-fire |
| `encode_failed` | 1 | **1** | unchanged |

The movement is exactly one-directional and exactly the intended 19 molecules. **That
`achiral_or_preexisting_fold` is unmoved at 157 is the second-most-important number here**: it
is the direct evidence that the left conjunct works. A veto without it would have swallowed all
157 of those molecules as well, still reported `REGRESSION_raw_collapsed = 0`, and looked like
a total success while destroying the fold's entire benefit.

## 7. What the veto costs — 171 of 393 gains survive, +3.42 points

```bash
PYTHONPATH=$PWD/src .venv/bin/python tools/fold_transition_sim.py \
    --sweep tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest --arm veto
```

| | fold alone | **fold + veto** |
|---|---:|---:|
| `key_equal/slot_renumber → byte_exact` | 393 | **171** |
| `byte_exact` | 3623 → 4016 | **3623 → 3794** |
| points | **+7.86** | **+3.42** |
| `facmer_divergent` | 16 → 16 | **16 → 16** |
| moved in a bad direction | 0 | **0** |
| excluded as `drift` / unavailable | — | **0 / 0** |

**The veto keeps 43.5% of the fold's gain and gives back the rest.** v0.4.11 bounded the safe
set at *"at most ~172 of 393 (~3.44 of the 7.86 points)"* from its collapse count; the filter
independently lands on **171 / +3.42**. Two different measurements — one counting collapses in
a mirror audit, one counting survivors through the shipped predicate — agreeing to a single
molecule is the strongest evidence in this lane that the veto is separating the right set.

**0 drift across all 393** is a second result carried by the same run: today's encoder
reproduces the v0.4.8 strings for every mover's input *and* generated structure, so
v0.4.9 / v0.4.10 / v0.4.11's byte-identity claims — and the carry-forward licence that rests on
them — are **measured, not inherited**.

> The gain is bankable only because §6 read 0. A surviving-gain count says what the veto
> *kept*; only the audit says whether what it kept is safe. Both were run before the number
> above was quoted anywhere.

### 7.1 The default path is provably untouched

```bash
bash tools/gate_v047.sh arm1
#   [gate] python=.../.venv/bin/python (rdkit 2025.09.3)   <- the right interpreter
#   [gate/arm1] sentinel OK: #DONE 62
#   [gate/arm1] PASS -- byte-identical to golden
```

62/62 fixtures byte-identical with both levers off, and the gate resolved the **pinned**
rdkit 2025.09.3 rather than a sibling project's venv — the v0.4.9 trap where a byte-identity
gate silently ran on rdkit 2025.09.2 and reported MISMATCHes that read as code regressions.

Full suite: **988 OK** (skipped 3, expected failures 5).

## 8. Files

`src/oinsmiles/oin/fold_parity.py` (new) · `src/oinsmiles/utils/perception_tmc.py`
(pristine-conformer capture + post-pass call site) · `src/oinsmiles/oin/levers.py` ·
`tests/unit/test_fold_parity.py` (15 tests) · `tests/fixtures/fold_parity/` (3 vendored,
oracle-confirmed chiral) · `tools/fold_transition_sim.py` (new).
