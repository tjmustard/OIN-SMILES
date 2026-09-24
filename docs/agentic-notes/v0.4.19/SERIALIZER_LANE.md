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

## 2. The whole cohort: changed set, canonicality, offline re-score (`tools/v0419/serializer_ab_report.py`)

Both are **canonicality levers** on the encoder, so the v0.4.11 rule applies — mirror-audit and
renumbering-audit the whole cohort, never a subset — and because `smiles_1` moves, the moved rows
must be generated live while the rest can be re-scored offline (generation is byte-deterministic at
this load, v0.4.18).

**2a. Changed-string set** — `tools/census/e_selfconsistency.py`, all 5,000, under both levers
(`e_selfconsistency_fix2.jsonl`, 55,000 encodes, 13.85 CPU-h, 11 encode timeouts = the sweep's 11
`P_E1_COVERAGE`): base differs from the release sweep's `smiles_1` on **80 molecules**:
`E1_HCOUNT` FAIL 57 · `E1_GRAPH` FAIL 21 · `E1_GRAPH` false PASS 1 · `P_E1_COVERAGE` 1 (SUGQOA — a
300 s encode timeout in the sweep that finished this time; not the levers). **4,909 unchanged**,
which includes every verified pass. The levers touch the E1 population and nothing else; the two
runs' base strings differ on exactly those 79 (the instrument fired where it should and nowhere
else). One shipped string timed out in the lever run (UVIWIG, load; see §4).

**2b. Canonicality audit** — lever run vs shipped run (`results-v0.4.17-exactfold/
e_selfconsistency_exact.jsonl`, valid for v0.4.18: `smiles_1` moved on 0 rows), per transform:

| transform | shipped byte-stable | lever byte-stable | became unstable | became stable |
|---|---:|---:|---:|---:|
| again (determinism) | 4,989 | 4,988 | 1 (UVIWIG timeout) | 0 |
| rewrite / rotate | 4,989 | 4,989 | 0 | 0 |
| renumber ×3 | 4,777 / 4,770 / 4,777 | same | 0 | 0 |
| noise 0.02 Å ×3 | 4,868 / 4,870 / 4,870 | same | 0 | 0 |
| mirror | 3,262 | 3,262 | 0 | 0 |

Mirror 2×2 against the ruler unchanged: NON-CANONICAL 120 → 120, NON-INJECTIVE 163 → 163;
E2-fragile 502 → 502 (0 became fragile, 0 became stable). **A null on every axis a canonicality
lever could break.** The 44 changed-string molecules that are renumbering-fragile were fragile
before (perception, `E2_P_FRAGILE`).

**2c. Offline re-score of the 4,909 unchanged rows** — `tools/honest_rescore.py` under the levers
(4,758 stored structures re-encoded; 147 no structure): self-consistent **pass→fail 0**,
fail→pass 5; the census fault recomputed per row through `attribution_table.attribute()` (control:
shipped inputs reproduce the table 4,920 / 4,920): **VERIFIED 3,920 → 3,920, 0 losses, 0
gains.** The 5 self-consistent gains (DOCPAO, IFAPUD, MIBFEL, NEFNER, NOYTUS) are strings the
levers do NOT repair whose generated structure now re-encodes to the same wrong string — false
passes, the class the two-number rule exists for. 4 re-encode errors under the levers (AFIROW,
HOBBUY, MUZZUC, RIPDEA) are `hard_fail` rows whose shipped honest re-encode also failed.

## 3. The live A/B on the 80 (`tools/v0419/run_serializer_ab.sh` → `post_serializer_ab.sh`)

ON arm = the harness, 6 shards, `--mol-timeout 300`, under both levers; OFF = the release sweep's
own rows. Scored self-consistent AND VERIFIED, the latter by the census's own `attribute()` with
each arm's bucket, ruler verdict, parse-back and self-consistency records (the string clause is
recomputed per arm — the ON arm's string is the lever's). Controls: the recomputed OFF fault must
equal the census table row for row; `smiles_1` must DIFFER on every row (the encoder inversion of
the dead-lever check).

### 3a. Result — 80 molecules, ON arm commit `d3adba15` (`serializer_ab_report.txt`)

Controls: recomputed OFF fault == census table **80 / 80**; `smiles_1` identical across arms on
**0** rows. Per-lever attribution of the 80 (single-lever encodes): **H_FAITHFUL 64 ·
RC1_PROPAGATE 15 · neither 1** (SUGQOA, the timeout row) — the two levers are disjoint.

| | OFF (sweep rows) | ON | transitions | net |
|---|---:|---:|---|---:|
| self-consistent (`byte_exact`, honest) | 1 | 29 | fail→pass 29, pass→fail 1 | **+28** |
| VERIFIED (census fault `NONE`, each arm's own reader) | 0 | 24 | fail→pass 24 | **+24** |

By lever: VERIFIED gains **H_FAITHFUL 20, RC1_PROPAGATE 4**; self-consistent gains 24 / 5, the
one loss RC1's. That loss is **MOSLEL** — a byte-exact FALSE pass in the sweep (`E1_GRAPH`, the
corrupted string round-tripped through the generator unchanged); with the corrected string the
generator builds the η ring one carbon short in 7 s and the honest re-encode says so
(`G_DIASTEREOMER`). A false pass became an honest fail; VERIFIED did not move.

ON-arm fault of the 56 still failing VERIFIED: `E2_P_FRAGILE` 18 (the string now describes the
input and is renumbering-fragile — it always was; the census's rule order surfaces it once E1
clears), `G_NOTHING` 11, `G_CONSTRUCTION/DETACHED` 9, `E1_GRAPH` 8 (7 H rows whose sphere is a
perception trim, 1 RC1), `E2_P_OTHER` 3, `E1_NONCANONICAL` 2, `SPHERE_DIFF` 2, `G_DIASTEREOMER`
1, `LIGAND_DIFF` 1, `E1_HCOUNT` 1.

### 3b. Runtime: a wrong string costs the generator its whole budget

| lever rows | Σ `elapsed_s` OFF | ON | rows at the 300 s budget OFF → ON |
|---|---:|---:|---:|
| H_FAITHFUL (64) | 19,891 s | 5,417 s | **49 → 8** |
| RC1_PROPAGATE (15) | 4,475 s | 1,538 s | **11 → 2** |

The lane map's last "do not": the "compute floor" for these rows was a wrong molecule, not
throughput. 24,666 → 7,254 s over the 80.

### 3c. 🔴 `OIN_H_FAITHFUL` is a WRITER + READER lever

`generation/metallogen_adapter.py` (the "authoritative" block, gated on
`hydrogen_faithfulness_enabled()`) preserves a bracketed `[C]` only while the lever is on. The
census's parse-back ruler IS that reader. The same 80 ON strings parse back ISO&SAME on **71 rows
under the ON reader and 15 under the shipped reader**; the H verdict depends on the reader's lever
on **56** rows. So the first post-processing pass (shipped reader, by accident) reported VERIFIED
+4, and the arm-consistent pass +24 — both are in `serializer_ab_report.txt`
(`READER COUPLING`), and `parseback_shippedreader.jsonl` is frozen beside `parseback.jsonl`.
For a promotion decision the ON reader is the one that ships; for the notation it means a string
written with the lever and read by a build without it is a different molecule. Promote it
everywhere or nowhere.

## 4. Projection and what the owner decides (the two serializer levers)

| | v0.4.18 release sweep | with both levers (live 80 + offline 4,909, exact by determinism) |
|---|---:|---:|
| self-consistent | 4,347 = **86.94%** | 4,380 = **87.60%** (+28 live, +5 offline false passes) |
| VERIFIED | 3,920 = **78.40%** | 3,944 = **78.88%** (+24, 0 losses anywhere) |

Frozen `measurements/v0.4.19-serializer/` (16 files); `tools/v0419/freeze_stage.py --verify`
re-derives every line above from the public tree alone. **Nothing promoted; no sweep run.** If
the owner promotes either: full sweep, `/refreeze-goldens`, census `EXPECTED`, and — for
`OIN_H_FAITHFUL` — the reader coupling in §3c stated in the CHANGELOG. `OIN_RC1_PROPAGATE` alone is
a pure serializer fix (+4 VERIFIED, one false pass turned honest, 11 → 2 budget rows, nothing
else moves) and has no coupling.


## 6. H4 — `OIN_CAP_IGNORES_METAL`: perception's valence cap (started 2026-09-24, in flight)

`perception_core`'s per-atom valence cap counts a ligand atom's metal contact as a bond; when the
atom is over its maximum valence, `remove_weakest_bond` deletes the neighbour with the largest
excess `d − r_i − r_j`. Whenever the metal contact is the *shorter* by the radii (Pd–Se 2.38 Å
= −0.21, Se–C 1.94 Å = −0.02) the ligand bond goes: KICSUM's PhSe–CH₂ is split into two
fragments, OBILAM's Si–C ring is opened (`E1_GRAPH/LIGAND_DIFF`; the map's H4).

**Rule v1 — "a heavy non-metal atom's cap ignores its metal contact" — REFUTED on the cohort.**
Probe (296 rows): `LIGAND_DIFF` 14/56, `SPHERE_DIFF` 7/76, `HCOUNT` 3/84 repaired, controls 80/80
— after hydrogen was exempted (8 rows died on "Explicit valence for atom H, 2": an H between a
carbon and the metal must lose one bond, and the excess rule is right there). But the whole-cohort
audit (`e_selfconsistency_capv1.jsonl`) moved **74 strings including 5 verified passes** and made
**11 molecules renumbering-fragile**. The AC diff over all 74 (`xyz2AC_obabel` under both rules)
separates the two populations cleanly by the excess of the metal contact at the atom whose cap
changed:

| | metal-contact excess | examples |
|---|---|---|
| repairs (a ligand bond kept) | **−0.30 … +0.07 Å** (a few to +0.16) | Pd–Se −0.21, Ti–Si −0.10, η-ring C +0.05 |
| regressions (a long contact kept) | **+0.10 … +0.44 Å** | JIXTES tBu *methyl* → Ni donor +0.37; GUSRAN C–F···Y 7th donor +0.34; IGOBOX B–H···Fe +0.20 |

The old rule was right for long agostic / C–F / B–H contacts — they *are* the longest neighbour —
and wrong only when the metal contact is a bond by the radii.

**Rule v2 — exempt only a SHORT metal contact** (`CAP_EXEMPT_EXCESS = +0.10 Å`; a longer one
counts and is cut as before; hydrogen unchanged). Probe over the 296 + the 25 other v1-moved
rows (`pop1d_report.txt`): `LIGAND_DIFF` 14/56, `HCOUNT` 3/84, `SPHERE_DIFF` 1/76 (the six v1
"repairs" there were long contacts), **the 5 verified passes byte-identical**, and 6 `P_DETACHED` /
`DATA_MULTI` rows now parse back ISO (a bridging Se/B ligand that was written slot-less). All
three levers together: 15 + 15 + 59 = 89 of 216 E1 strings describe their input.

**Rule v2 on the whole cohort — REFUTED by the offline re-score.** The audit itself looked fine
(no verified pass's string moved; fragile 502 → 495, NON-CANONICAL 120 → 119, NON-INJECTIVE
163 → 161; 9 became fragile — a hard +0.10 Å edge is noise-sensitive by construction). But
`honest_rescore.py` under v2 gave **195 self-consistent / 147 VERIFIED losses on rows whose input
string is unchanged** — the *generated* structures re-encode differently. Mechanism (AC diff on
the generated xyz): in a compressed generated structure two halides register a spurious
Br···Br / I···I contact at ≈ +0.3 Å; the shipped cap removed it *because it counted the metal
bond* (Br: M + Br = 2 > 1 → cut the longer). Exempting the short metal contact from the count
left the halide at 1 ≤ 1 and the spurious contact survived: `[Br]{0}[Br]{1}`. The count was doing
the coordination sphere's job for saturated ligand atoms.

**Rule v3 (current):** the count is the shipped one — every neighbour, the metal included — and
only the *victim* changes: a ligand bond with excess < +0.10 Å is never cut while a short metal
contact exists; that contact is dative, leaves the count, and nothing is cut. Smoke: KICSUM and
OBILAM repaired; the five verified passes and UMENEG byte-identical; 12/12 sampled v2 losses
re-encode to `smiles_1` again. **Offline re-score under v3 (4,760 structures): at most 5 passing
rows differ** (v2: 195), pending the audit's changed set to classify them.

**Rule v3 on the whole cohort (`e_selfconsistency_cap.jsonl`, `rescore_cap/`):** the string
moves on **31** molecules — `E1_GRAPH` 22 (12 of them byte-exact false passes), `P_DETACHED` 6,
`DATA_MULTI` 1, `E1_HCOUNT` 1, the timeout row — and on **no verified pass**. Canonicality: fragile
502 → 492 (2 became fragile, both failing rows; 12 became stable), mirror 3 unstable / 1 stable
(NON-CANONICAL 120 → 121, NON-INJECTIVE 163 → 162): a near-null, slightly positive. Offline
re-score of the 4,958 unchanged rows: self-consistent **5 losses** (4 are `G_CONSTRUCTION` false
passes turning honest; OHUTIV is a verified pass whose generated structure now re-encodes with a
different H count, `NONE → G_HCOUNT`) / 13 gains (11 false passes, 2 real); **VERIFIED −1 / +2**.
With all three levers (`fix3`): 109 strings move (80 + 31 − 2 overlap), same stability columns,
same offline transitions.

### 6a. Live A/B, rule v3 (`serializer_ab_report_cap.txt`, `serializer_ab_report_fix3.txt`)

| arm | cohort | controls | self-consistent | VERIFIED | runtime Σ OFF → ON |
|---|---:|---|---|---|---|
| `cap` alone | 31 | fault 31/31; `smiles_1` differs 31/31; reader coupling none (22/22) | 13 → 16 (+5/−2, the −2 false passes turning honest) | **0 → 6 (+6/−0)** | 2,300 → 2,314 s |
| all three (`fix3`) | 109 | 109/109; 109/109; H-reader coupling 93 vs 37 | 14 → 44 (+33/−3) | **0 → 30 (+30/−0)** | 26,234 → 8,198 s |

VERIFIED gains by census fault, `fix3`: `E1_HCOUNT` 20, `E1_GRAPH` 8, `P_DETACHED` 2 — exactly the
serializer lane's +24 plus the cap's +6: **the three levers are additive**. Still failing on the
ON arm: `E2_P_FRAGILE` 20, `E1_GRAPH` 14, `G_NOTHING` 12, `DETACHED` 11, `E1_NONINJECTIVE` 5,
`E1_NOT_ENCODED` 4, …

### 6b. Projection (live arm + offline re-score, unchanged rows exact by determinism)

| configuration | self-consistent | VERIFIED |
|---|---:|---:|
| v0.4.18 release sweep | 4,347 = 86.94% | 3,920 = 78.40% |
| `OIN_H_FAITHFUL` + `OIN_RC1_PROPAGATE` (§4) | 4,380 = 87.60% | 3,944 = 78.88% |
| `OIN_CAP_IGNORES_METAL` alone | 4,358 = 87.16% | 3,927 = 78.54% |
| **all three** | **4,378 = 87.56%** | **3,951 = 79.02%** |

(The all-three self-consistent number is below the pair's because the cap turns 4 + 2 false
passes into honest fails and OHUTIV into a real loss; VERIFIED, the number to quote, is
monotone.) Frozen `measurements/v0.4.19-serializer/` (42 files, `main` @ `c7b0be7d`);
`tools/v0419/freeze_stage.py --verify` re-derives all 26 cells. **Nothing promoted; no sweep run.**

**What the cap lever costs:** one verified pass (OHUTIV, `NONE → G_HCOUNT`: its generated
structure's H count re-encodes differently under v3), two failing rows newly renumbering-fragile,
and a hard +0.10 Å edge that 0.02 Å noise can cross (noise columns net +16 stable / −7). Against
+6 VERIFIED live, +2 offline. The owner's call; if promoted, sweep → `/refreeze-goldens` →
census `EXPECTED` → CHANGELOG.


## 7. Not done / not this lane

- H3 (adapter `is_haptic` on macrocycles), H4 (valence cap), H5 (η trims in perception): reader
  and perception defects; each moves the census's verdict and must be re-attributed, not assumed.
- No lever was promoted. No sweep was run.
- UVIWIG: shipped string encodes in the sweep (300.1 s), timed out at 300 s in the lever run (which
  shared the box with a 4-worker re-score); alone under the levers it encodes in 108.6 s. A slow
  encoder on the budget boundary, not a lever effect. Its lever string was not measured.
