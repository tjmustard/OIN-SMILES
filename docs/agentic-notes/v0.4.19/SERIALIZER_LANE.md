# v0.4.19 — the serializer lane: two string-side levers, measured with the parse-back ruler

> Branch `research/v0419-serializer` (worktree `../oin-v0419`, off `v0.4.18`). Started 2026-09-23
> under the owner's week-long delegation, which covers v0.4.18; **a v0.4.19 lever is MEASURED and
> left default-OFF for the owner.** Read-only map that chose the targets:
> `spec/handoffs/v0.4.19/SERIALIZER_LANE_MAP.md` (gitignored). Results
> `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-serializer/` (gitignored; frozen copies below).

## 0. The ruler, and why the old A/B was blind

The census assigns `E1_HCOUNT` / `E1_GRAPH` by **parse-back**: the string is read back to a graph
with no 3D — at parser level and at the generator adapter's level
(`tools/census/string_sufficiency.py::parseback_one`) — and compared to the input xyz. For an
ENCODER change that ruler is exact: no generator, no seed, no budget. `OIN_H_FAITHFUL` was held
off in v0.4.6 because a **round-trip** A/B over the 45-molecule `Atom count mismatch` population
"moved nothing" (`levers.py`): 55 of the 72 `E1_HCOUNT` rows of that era were `hard_fail`
timeouts, so a repaired string still went nowhere through the generator. The lane map's rule 6.1,
"do not measure this lane with the round trip", is the whole reason the result below exists.

Instrument: `tools/v0419/lever_parseback.py` — fresh process per molecule and arm, levers WRITTEN
("1"/"0"), `sys.path.insert` of the worktree's `src` (asserted on every worker's
`oinsmiles.__file__`), then `parseback_one` on the live string. Populations from the release census
(`results-v0.4.18-release-census/attribution_table.tsv`): every row whose shipped string fails to
describe its input at the adapter level, plus 80 verified passes (40 η / 40 non-η, seed 419) as the
control a canonicality lever must not move.

**Control first:** the live "shipped" arm reproduced the release sweep's `smiles_1` on
**296 / 296** rows. The arm that is supposed to be the baseline IS the baseline.

## 1. Result — 296 molecules × arms (`pop1.jsonl`, `pop1b.jsonl`)

| population | n | shipped | `OIN_H_FAITHFUL=1` | RC1 fail-safe forced | **`OIN_RC1_PROPAGATE=1`** | **both levers** |
|---|---:|---:|---:|---:|---:|---:|
| CONTROL (verified passes) | 80 | 80 ISO&SAME, 0 changed | 80, **0 changed** | 80, **0 changed** | 80, **0 changed** | 80, **0 changed** |
| `E1_HCOUNT` (adapter H ≠ input) | 84 | 0 | **56** (57 strings changed) | 0 | 0 | **56** |
| `E1_GRAPH` / `SPHERE_DIFF` | 76 | 0 | 0 (6 changed) | **14** (14 changed) | **14** (14 changed) | **14** (20 changed) |
| `E1_GRAPH` / `LIGAND_DIFF` | 56 | 0 | 0 (1 changed) | 0 (1 changed) | 0 (1 changed) | 0 (2 changed) |

ISO&SAME = adapter graph ∈ {ISO, ISO_MARGINAL, ISO_CLASH} **and** adapter H decoration SAME.
Zero errors, zero timeouts in every cell. The two repairs are **additive** (70 = 56 + 14) and the
fix and its diagnostic (RC1's own fail-safe, forced by monkeypatch) produce **byte-identical
strings on all 296 rows** — the permutation, propagated, lands exactly where "no permutation"
lands, which is what step 7 of `get_oin_string` (re-sort by slot) predicts.

### 1a. `E1_HCOUNT`: `OIN_H_FAITHFUL` repairs 56 of 84 (H1 = H6 of the map)

The shipped writer neutralises a hypovalent 0-H anion and RDKit's bare-write rule then emits `C`
for a carbon the reader gives one H. With the lever the 0-H carbon is bracketed (`[C]`); the
adapter reads it back with the input's H count. The 56 are 50 non-η / 6 η. The **28 it does not
repair** are reader-side or donor-convention, not the writer: `c{0}` aryl-C donors inside 3+-donor
macrocycles that the adapter's `is_haptic` treats as haptic and gives a phantom H (BUWHAD, JESFIX,
KAXPAA — H3 of the map), benzyl/aryl carbanion donors `C{0}c1ccccc1`, hydrides `[H]{0}`,
`[BH3]{0}`. One string (EBUBAH) changes under the lever and stays DIFF.

### 1b. `E1_GRAPH/SPHERE_DIFF`: the RC1 rank swap never reached the writer (H2 of the map)

`OINDiscreteAligner` step 3b (RC1) re-ranks "same-mass" η fragments by canonical content — and
same-mass means the same **first binding atom** (`chem_id[0]`), so a Cp and an allyl, an indenyl
and a butyne, a Cp* and an ethylene all qualify. Each w-tag entry keeps the local indices of the
fragment it was computed on; `get_oin_string` reads the entry's rank against the **un-permuted**
fragment list. TULTAX (Ru, Cp + allyl + picolinate):

```
shipped  [Ru_TET].O{1}C(=O)c1cc(Cl)ccn{0}1.[CH2]{2>}[CH]{2}=[CH2]{2}.[cH]{3>}1[cH]{3}[cH]{3}cc1
lever    [Ru_TET].O{1}C(=O)c1cc(Cl)ccn{0}1.[cH]{2>}1[cH]{2}[cH]{2}[cH]{2}[cH]{2}1.[CH2]{3>}[CH]{3}=[CH2]{3}
```

The Cp received the allyl's three markers (two ring carbons written bare, `cc1`), the allyl
three of the Cp's five (the surplus dropped at `inline.py`). ULORIY (W): ethylene written
`C=[CH2]{2}` and a Cp* **methyl** carbon marked η-bound, `[CH3]{4>}c{4}1c(C)c(C)c(C)c1C`. The
string names a different coordination sphere from the input; the generator builds it as written.

Fix (`OIN_RC1_PROPAGATE`, `_HELD_OFF`): the aligner records its non-identity permutation
(`OINDiscreteAligner.rc1_rank_map`, `{new_rank: original_rank}`, reset per call) and
`get_oin_string` permutes `fragments_data` to match before reading the w-tag — the map's
recommended fix, not the re-key it warns against. Pinned both ways on TULTAX in
`tests/unit/test_rc1_propagate.py` (+ ferrocene as the no-op control). The 14 repaired are all η
(the lever can touch nothing else by construction); the 62 residual SPHERE_DIFF are 32 non-η /
30 η with `cn_str` unreadable (perception trims, H5 of the map — not this lane).

### 1c. `E1_GRAPH/LIGAND_DIFF`: neither lever (H4 — perception's valence cap)

As predicted by the map: 0 of 56 move. A perception lane, not a serializer one.

## 2. What is still owed before either lever is the owner's to promote

Both are **canonicality levers** on the encoder, so the v0.4.11 rule applies — mirror-audit and
renumbering-audit the whole cohort, never trust a subset — and because `smiles_1` moves, the
sweep of record cannot be re-scored offline for the moved rows; they must be generated live.

1. ✅ running — `tools/census/e_selfconsistency.py` over all 5,000 under both levers
   (`e_selfconsistency_fix2.jsonl`): base string ⇒ the COMPLETE changed-string set; determinism,
   rotation, 3 renumberings, 3 noise, mirror ⇒ the canonicality audit against the shipped run
   (`results-v0.4.17-exactfold/e_selfconsistency_exact.jsonl`, valid for v0.4.18 because
   `smiles_1` moved on 0 rows there).
2. Attribute each changed string to one lever (base-only encode under each single lever, changed
   set only).
3. Offline re-score of every UNCHANGED-string row: re-encode the sweep's generated xyz under the
   levers; a pass whose re-encode now differs is a loss the levers would cause with no generation
   (exact — generation is byte-deterministic at this load, v0.4.18).
4. Live harness A/B on the changed set (shipped vs both levers), scored self-consistent AND
   VERIFIED (`g_vs_input` + parse-back on the generated structure).
5. Freeze; handoff; the owner decides. Sweep only if promoted.

## 3. Not done / not this lane

- H3 (adapter `is_haptic` on macrocycles), H4 (valence cap), H5 (η trims in perception): reader
  and perception defects; each moves the census's verdict and must be re-attributed, not assumed.
- No lever was promoted. No sweep was run.
