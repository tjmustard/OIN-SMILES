# `OIN_ACCEPT_SCORED` — promotion gates (v0.4.7, lane L2-promote)

**Status: MEASURING.** Results sections are filled as runs land. Method and mechanism are
settled and are recorded here first, deliberately: they are the part that does not survive an
interrupted session otherwise.

---

## 0. What the lever does

`src/oinsmiles/generation/metallogen_adapter.py::_reencode_key_matches` decides whether a
pooled conformer is accepted. Three tests, in order:

```
1.  fast = _reencode_oin_fast(cmol)              # cheap contract-mol re-encode
    if fast is not None and key(fast) != target: return False        # CHEAP VETO
2.  if not independent_confirm and fast is not None: return True     # <-- OIN_ACCEPT_SCORED
3.  full = _reencode_oin(m)                      # XYZToSMILES().convert(xyz)
    return key(full) == target                                       # STRICT
```

The lever turns step 3 off. Step 3 is the only test in the predicate that does not reuse the
generator's own bond graph.

The harness scores with `key(get_oin_string(gen.mol, coords))`, which is exactly what step 1
computes. So the lever makes acceptance agree with the metric.

---

## 1. Why four gates and not the one already measured

A 22-molecule A/B (`tools/ab_accept_scored.py`, run before this lane) was decisive on runtime:
median 16.01s → 3.63s, total 1980.7s → 626.0s, `>30s` 10 → 3, no pass change either way.
Four things stood between that and a promotion.

### G3 — the A/B measured the wrong equality

`passed` compares `canonical_roundtrip_key`, which **folds slot renumbering, fragment reorder
and metal-@**. The v0.4.7 wave's contract is a **byte-identical OIN string**. A conformer can
satisfy the key and still emit different bytes, so a key-only A/B cannot clear a byte gate.
Nobody had checked whether `smiles_2` moves under the lever.

Measured as `sha256` of `get_oin_string(gen.mol, coords)` — byte-for-byte the string
`tools/test_dataset_roundtrip.py` records as `smiles_2` — per molecule, per arm, with
`sha256(smiles_1)` as a control that must not move (the lever is generation-side).

**This is pure string equality and therefore fully valid at load 40+, which the seconds are
not.** It ran first.

> **Determinism is a premise of G3, not a conclusion of it.** An A-vs-B sha difference is only
> attributable to the lever if `sha_out` is stable across two runs of the *same* arm. Any
> molecule whose sha differs between arms gets an A-vs-A′ control (`--single-arm`) before any
> verdict is stated.

### G1 — the quality arm measured nothing

The shipped A/B printed `clashes 0 over 0 mols` in **both** arms, all 44 measurements.
`clash.mol_clash_count()` duck-types on `mol.atom_list` and returns 0 on `AttributeError` by
design (`clash.py:110-120`), and `gen_result.mol` is a bare `rdkit.Chem.rdchem.Mol`. Fixed in
`49f8d69b`: `clash.vdw_clash_count(positions, znums)` on the returned coordinates, reporting
`clash_vdw`, `clash_severe`, `worst_overlap`, and a `clash_measured_on` **denominator**.

The failure mode there was not a wrong number, it was a missing denominator — `0 over 0` and
`0 across 22` printed identically.

This matters because accepting earlier bypasses `_select_by_geometry`'s clash-first ranking
(`clash.VDW_ACCEPTANCE_ENABLED` is ON by default). `worst_overlap` is continuous, so the
verdict is not hostage to where `clash_cutoff` sits.

### G2 — the pass-rate arm is circular

`passed` is computed with `key(get_oin_string(gen.mol, coords))` — the same predicate the lever
accepts on. By construction it cannot detect what dropping step 3 costs, so `18/22 → 18/22` is
not evidence of no loss. The lever's own docstring says it buys "latency and metric-fidelity,
NOT genuine losslessness."

Added arm: full `XYZToSMILES().convert()` on the generated XYZ, compared to the input OIN.

Two scoring decisions, both settled by measurement rather than assumption:

- **Scored on the key, not on bytes.** The first molecule measured settles it. ODEWID_comp_0
  key-passes independent re-perception while `oin_indep` differs from `oin_in` by donor slot
  renumbering alone (`N{3}` ↔ `N{1}`). A byte-scored independent arm would have reported a loss
  that isn't one.
- **`circular_only` is a BETWEEN-ARM DIFFERENCE, never an absolute cost of the lever.** See §3.

### G4 — n=22 was selected to exhibit the gap

Widened to a guard population of currently-**passing** molecules, drawn from
`results-v0.4.5-rebaseline/individual_reports/*.json` (936 reports: **634 success**, 302
failed). Question: does any molecule that currently passes stop passing?

*(`elapsed_s` lives inside each report's `metrics` dict; read from the top level it silently
yields 0.)*

---

## 2. PREFILTER_VETO — an accuracy defect, arm-independent, NOT this lane's to fix

AROHIA_comp_0 runs opposite to the lever: the **cheap prefilter never matches (0/48)** while
the **strict test matches 16/48**.

Reading the order in §0: the cheap veto at step 1 returns `False` *before* both the lever
branch (step 2) and the strict test (step 3). So those 16 conformers are unreachable in
**production in both arms** — the lever can neither cause this nor fix it. Consistent with the
A/B's observation that AROHIA is unmoved (42 → 41.6s).

The defect is that `_reencode_oin_fast`'s docstring claims it was "verified to yield the same
eta winding as the full `XYZToSMILES().convert` path". AROHIA is a counterexample: the cheap
encode disagrees with the full one, in the direction where cheap is wrong. A prefilter
documented as "may only REJECT" is only sound if it rejects nothing the strict test would take.

`tools/probe_accept_gap.py --capture-veto` counts this as `n_strict_only` and keeps the two
disagreeing strings for the first instance, because counts say it happens and only the strings
say why.

**Confirmed from production data, not just from the probe.** AROHIA's row in this lane's run
is the exact inverse of every other molecule in the cohort:

```
[A-default] AROHIA_comp_0   53.75s pass=False indep=True  sha=596fde4bb73155a7 clash=0 worst=0.8152
```

`passed` is computed with `get_oin_string(gen.mol, coords)` — the same function
`_reencode_oin_fast` uses — and it says **mismatch**. `indep_passed` runs the full
`XYZToSMILES().convert()` and it says **match**. So the structure the generator returned *does*
round-trip; the cheap encoder is what is wrong about it. Everywhere else in this cohort the
cheap test is the permissive one and strict is the strict one. On AROHIA the ordering inverts,
which is precisely the condition under which a "may only REJECT" prefilter becomes unsound.

Two consequences, both worth stating separately:

1. **It is an acceptance defect.** 16 of 48 conformers that the strict test would take are
   unreachable, with the lever on or off.
2. **It is also a measurement defect.** The harness scores with the cheap path, so AROHIA is
   recorded as a round-trip FAILURE while independent re-perception says it round-trips. The
   reported pass rate is wrong on this molecule in the pessimistic direction.

**Recorded, not fixed, per lane scope.**

---

## 3. The default arm already fails independent re-perception

Measured, arm A (`OIN_ACCEPT_SCORED=0`), on the CHEAP_ONLY class:

| molecule | arm | `passed` (circular) | `indep_passed` | sha256(smiles_2) |
|---|---|---|---|---|
| HEJXIF_comp_0 | A default | True | **False** | `169084c2a24843a9` |
| WIWRIE_comp_0 | A default | True | **False** | `ec9b6c3cbdcd07be` |
| NOMMOU_comp_0 | A default | True | **False** | `5d8208fcbb9313c2` |

**Mechanism.** These are the CHEAP_ONLY class of `docs/eta_accept_gap_cohort.md`: cheap matches
many conformers, strict matches **none**. So `accept_fn` never fires, the pool fills to
completion, and `_select_by_geometry` returns a best-by-geometry conformer **that was never
required to pass acceptance at all**. Acceptance is not the only path to a returned structure.

**Consequence for G2, and it is the load-bearing one.** The default path already ships
structures that fail independent re-perception. So `circular_only` cannot be read as an
absolute cost attributable to the lever — it is a property of the molecule and of the fallback,
present with the lever off. Only the **difference between arms** (`INDEP REGRESSIONS`:
A indep-pass → B indep-fail) is the lever's price.

Stated plainly because it inverts the intuitive reading of the gate: "the default guarantees
independent re-perception, the lever gives that up" is **false**. The default guarantees it
only for conformers that were actually *accepted*, and on the slowest molecules in the cohort
nothing ever is.

---

## 4. Results

### 4.1 G1 — structure quality. Measured for the first time. **B is better in aggregate and worse on two molecules.**

> ⚠ **CORRECTION.** An earlier revision of this file (commit `9f6ba5ba`) concluded from a
> **partial** log that "the quality cost is one clash on one molecule." That was wrong, and it
> was wrong in the project's signature way: the partial log had arm A's `POVPIA_comp_0`
> (`clash=16, severe=7, worst=0.4344`) but had not yet reached arm B's pair for it. Arm B fixes
> POVPIA to `0 / 0.75`. Reading a comparison before both arms have landed is the same error as
> reading a metric before checking its denominator. The completed run is below.

Full run on the original 22-molecule cohort with the corrected `vdw_clash_count`
(`spec/handoffs/v0.4.7/rescued/ab_v2.json`). **17 molecules paired** — both arms returned
coordinates. Aggregates over the paired set only, so the denominators match:

| metric (paired, n=17) | A default | B scored |
|---|---|---|
| `clash_vdw` total | **16** | **2** |
| `clash_severe` total | **7** | **0** |
| molecules with any clash | 1 | **2** |
| `worst_overlap` min | **0.4344** | **0.7283** |

**The two directions do not agree, and both are real:**

- **Aggregate severity improves sharply.** B removes all 7 severe clashes and lifts the worst
  overlap in the whole cohort from 0.4344 (deep inside the 0.60 severe cutoff) to 0.7283 (above
  even the 0.75 contact cutoff). That is entirely POVPIA: A returns a 16-clash structure, B
  returns a clean one.
- **B spreads a single clash onto two molecules that had none.** DAKGON `0/0.8164` → `1/0.7283`
  and RATPEK `0/0.7599` → `1/0.7461`. Neither is severe.

So bypassing `_select_by_geometry`'s clash-first ranking is **not** uniformly worse, which is
what the gate was written to catch. It trades one catastrophic structure for two mildly
imperfect ones. Per-molecule: **B worse 3, B better 7, identical 7.**

*(The third "B worse" is HEJXIF, `0/0.7632` → `0/0.7562` — a slightly tighter worst contact with
zero clashes either way. Counted as worse only because `worst_overlap` is continuous and the
comparison is strict.)*

**Pass rate, same run:** A **16/22**, B **18/22**. **No regressions.** Two fixes:
`GAVSED_comp_0` and `QIDKUL_comp_0`.

> ⚠ **Both "fixes" are BUDGET-LIMITED, not correctness.** Arm A was **SIGKILLed at the hard
> cap** on both molecules (`exit -9`); it did not compute a wrong answer, it did not finish.
> The honest statement is "A did not complete within the budget and B did," not "the lever
> fixed a correctness bug." Both arms had an identical cap, so the comparison is fair — but
> this project has a documented history of timeout-shaped pass deltas being misread as
> correctness deltas (v0.4.4's 11 "regressions" were all 300s timeouts, 0 correctness). The
> +2 is real and it is worth having; it is a *throughput* result.

Per-molecule detail (14 molecules visible in the partial log, retained for the record):

| molecule | A default `clash_vdw / worst_overlap` | B scored `clash_vdw / worst_overlap` | verdict |
|---|---|---|---|
| **RATPEK_comp_0** | 0 / 0.7599 | **1 / 0.7461** | **B gains a vdW clash** |
| HEJXIF_comp_0 | 0 / 0.7632 | 0 / 0.7562 | B marginally worse overlap, no clash |
| AROHIA_comp_0 | 0 / 0.8152 | 0 / 0.8152 | identical |
| ODEWID_comp_0 | 0 / 0.7502 | 0 / 0.7502 | identical |
| QESRUE_comp_0 | 0 / 0.8206 | 0 / 0.8206 | identical |
| YIYGAP_comp_0 | 0 / 0.84 | 0 / 0.84 | identical |
| YIZHIY_comp_0 | 0 / 0.7588 | 0 / 0.7588 | identical |
| NOMMOU_comp_0 | 0 / 0.7579 | 0 / 0.7601 | B better |
| MEDZUR_comp_0 | 0 / 0.7577 | 0 / 0.7705 | B better |
| WIWRIE_comp_0 | 0 / 0.7539 | 0 / 0.7559 | B better |
| FEXYOZ_comp_0 | 0 / 0.7621 | 0 / 0.7659 | B better |
| KAQDOV_comp_0 | 0 / 0.7504 | 0 / 0.7629 | B better |
| ZITSIE_comp_0 | 0 / 0.7526 | 0 / 0.7532 | B better |
| GAVSED_comp_0 | KILLED at 330s cap | 0 / 0.7611 | B produces a structure where A produced none |
| **POVPIA_comp_0** | **16 / 0.4344 (7 severe)** | **0 / 0.75** | **B fixes a catastrophic structure** |
| **DAKGON_comp_0** | 0 / 0.8164 | **1 / 0.7283** | **B gains a vdW clash** |
| LIYXEY_comp_0 | 0 / 0.7513 | 0 / 0.7513 | identical |
| XUPTAF_comp_0 | 0 / 0.8072 | 0 / 0.8072 | identical |
| QIDKUL_comp_0 | KILLED at 330s cap | 0 / 0.7517 | B produces a structure where A produced none |

**`POVPIA_comp_0` is where the fixed metric earns its keep.** The default arm returns a
structure with 16 vdW clashes, 7 of them severe, worst overlap 0.4344. That defect was
completely invisible for as long as `mol_clash_count` was silently returning 0 on an
`AttributeError` — it would have shipped unnoticed. The lever happens to fix it, but the
finding that matters is that the *default* path produces it at all.

**Runtime (ADVISORY — load 40+ on 12 cores, wall-clock is meaningless here):** KAQDOV
242.5→7.5, ZITSIE 134.2→5.5, HEJXIF 194.6→13.0, WIWRIE 74.6→6.3, NOMMOU 39.9→4.1, RATPEK
69.0→11.8.

> **The molecules that got *slower* are load noise, not regressions.** ODEWID 1.45→2.21,
> QESRUE 3.59→5.63, YIYGAP 1.83→4.05, YIZHIY 4.54→5.88 are all sub-6s NO_GAP-class molecules
> where both predicates already fire at pool index 0 and the lever changes nothing it could
> change. Said explicitly so nobody later reads these four rows as a slowdown.

### 4.2 Same-arm reproducibility — an unplanned control that de-risks G3

G3's premise is that a sha difference between arms is caused by the lever. That only holds if
the generator is deterministic, and this project has a documented habit of confirming wrong
beliefs with measurements that never exercised the hard case. So the premise got checked.

Arm A was measured **twice**: this lane's run and the rescued prior run — different process,
different tool build, different `--hard-cap`. `clash_vdw` and `worst_overlap` are functions of
the returned coordinates alone, so agreement to 4 decimal places means the generator returned
the same structure.

| molecule | run 1 | run 2 | | molecule | run 1 | run 2 |
|---|---|---|---|---|---|---|
| DAKGON | 0 / 0.8164 | 0 / 0.8164 | | POVPIA | 16 / 0.4344 | 16 / 0.4344 |
| FEXYOZ | 0 / 0.7621 | 0 / 0.7621 | | RATPEK | 0 / 0.7599 | 0 / 0.7599 |
| HEJXIF | 0 / 0.7632 | 0 / 0.7632 | | WIWRIE | 0 / 0.7539 | 0 / 0.7539 |
| KAQDOV | 0 / 0.7504 | 0 / 0.7504 | | YIZHIY | 0 / 0.7588 | 0 / 0.7588 |
| MEDZUR | 0 / 0.7577 | 0 / 0.7577 | | ZITSIE | 0 / 0.7526 | 0 / 0.7526 |
| NOMMOU | 0 / 0.7579 | 0 / 0.7579 | | | | |

**11 of 11 identical** on molecules where both runs produced a structure. (A twelfth, GAVSED,
is not a disagreement: the rescued run killed it at its tighter 330s cap and produced nothing.)

This does **not** discharge the blocking control — any molecule whose `sha_out` differs between
arms still gets an explicit A-vs-A′ `--single-arm` re-run before a G3 verdict is stated. It does
mean generator nondeterminism is an unlikely explanation for any difference that appears.

### 4.3 G3 — the wave's gate. **PASSES on the 17 molecules where it is measurable.**

`sha256` of `oin_out` (= `smiles_2` = `get_oin_string(gen.mol, coords)`), arm A vs arm B, from
the completed `ab_v2.json`.

**17 comparable · 17 byte-identical · 0 different · 5 not comparable.**

| molecule | sha256(smiles_2) both arms | | molecule | sha256(smiles_2) both arms |
|---|---|---|---|---|
| AROHIA | `596fde4bb73155a7` | | POVPIA | `1c7ad422d5a60e3a` |
| DAKGON | `39ee1454f510da71` | | QESRUE | `a1286fbfa0d7af53` |
| FEXYOZ | `3d8095e19b48bfcc` | | RATPEK | `019faa066340a873` |
| HEJXIF | `169084c2a24843a9` | | WIWRIE | `ec9b6c3cbdcd07be` |
| KAQDOV | `bc968079864b73da` | | XUPTAF | `4514bdf2e19e867f` |
| LIYXEY | `7dea5d722d62a792` | | YIYGAP | `816ebda751edf06d` |
| MEDZUR | `a9ab9e1b834acd88` | | YIZHIY | `9527102647c15d81` |
| NOMMOU | `5d8208fcbb9313c2` | | ZITSIE | `684fb92c4a00f1ae` |
| ODEWID | `f00fdf52371f5cdb` | | | |

> ⚠ **The denominator is 17 of 22, not 22 of 22. The gate is UNMEASURED on 23% of the cohort.**
> The five: **GAVSED**, **QIDKUL** (arm A SIGKILLed at the hard cap, so no arm-A string exists);
> **XIQKOY** (arm A is the DEAD class — `MetalloGen failed to generate any conformers via
> OIN-direct` — and arm B was SIGKILLed); **YENDUS** (SIGKILLed in both arms); **DEJHEF** (both
> arms fail identically at parse with `UncoordinatedFragmentError: Fragment 3 ('II') has no
> binding slot`, which is lever-independent and produces no bytes on either side).
>
> Do not round this to "22/22". Notably, the unmeasured five are exactly the expensive and
> pathological molecules — the ones most likely to behave differently if they could be measured.

**`sha_in` control: holds.** On every molecule where both arms produced an input hash, it is
identical. An earlier version of this check reported "MOVED" on GAVSED/QIDKUL/XIQKOY; that was
a **false positive** — a SIGKILLed molecule yields a synthesized row with no `oin_in` at all
(`exit -9`), and comparing `None` against a real hash is a missing measurement, not a moved
value. Fixed in both `ab_accept_scored.py` and `promote_gate_report.py`, because that false
positive would have condemned a clean run as confounded.

**What this does to the determinism risk.** It defuses it rather than answering it. The risk was
"if shas differ, is that the lever or a nondeterministic generator?" — no sha differs, so the
question does not arise. The A-vs-A′ control still runs, but it is now a confirmation, not a
blocker.

**Byte-exactness against the INPUT** (a different question from A-vs-B): **A 16/17 produced,
B 18/19 produced.** The one non-byte-exact molecule in each arm is AROHIA — the PREFILTER_VETO
case of §2. So the lever does not degrade byte-exactness relative to input either.

### 4.4 G3 at its true denominator — this lane's own run

`spec/handoffs/v0.4.7/runs/ab22.json`, 22 molecules, both arms, `--hard-cap 500`.

**20 comparable · 20 byte-identical · 0 different · 0 one-arm · 2 neither.**

The two `neither` are the **DEAD class** — `XIQKOY_comp_0` (`MetalloGen failed to generate any
conformers via OIN-direct`) and `DEJHEF_comp_0` (`UncoordinatedFragmentError: Fragment 3 ('II')
has no binding slot`). Both fail identically in both arms. That is a *named class*, not a hole
in the evidence.

`only one arm produced a string: 0` is the line that matters most here: it closes the concern
that arm A's SIGKILLs were concealing a divergence. With the wider 500s cap, arm A completed
GAVSED (`a0cbabdc832a6a63`) and QIDKUL (`563845a7fc3be79a`) — both matching arm B exactly, and
both matching the *rescued* run's arm B, cross-run.

**`sha_in` control: OK on all 22.**

**Cross-run determinism, sha-level:** this lane's arm A vs the rescued run's arm A —
**17 identical, 0 differ.** Independent processes, different tool builds, different caps. The
A-vs-A′ control is therefore discharged as confirmation rather than a blocker: with 0 sha
differences anywhere, the question it was guarding against does not arise.

### 4.5 G2 — the independent re-perception arm. **This is where the lever costs something.**

| metric | A default | B scored |
|---|---|---|
| `passed` (circular metric) | 19/22 | 19/22 |
| **`indep_passed` / measured** | **15/20** | **7/20** |
| credited by circular metric ONLY (`circular_only`) | 5 | **13** |
| INDEP REGRESSIONS (A pass → B fail) | — | **8** |
| INDEP FIXES | — | **0** |

The eight: `KAQDOV`, `ZITSIE`, `FEXYOZ`, `DAKGON`, `RATPEK`, `MEDZUR`, `POVPIA`, `HIDCIH`.
**All eight are key MISMATCHES, not exceptions** — re-perception succeeded and returned a
different answer. None is a perception crash.

**Note what the project's own metric sees: nothing.** `passed` is 19/22 in both arms, no
regressions. The entire cost is invisible to the harness, which is precisely why the lever's
docstring insisted this arm exists.

#### The exact predictor — and the model it replaces

My first model was "the GAP class regresses." That model has an exception (`YIZHIY_comp_0` is
GAP-class and did **not** regress), and rather than wave it through, here is the model that has
none:

> **A molecule regresses on `indep` iff the lever returned a DIFFERENT conformer *and* arm A's
> conformer had passed `indep`.**

**Correct on 20/20 molecules that produced a structure in both arms.** ("Different conformer" =
`clash_vdw` or `worst_overlap` differs; given the demonstrated determinism, identical values
mean identical coordinates.)

| bucket | n | outcome |
|---|---|---|
| different conformer **and** A passed `indep` | **8** | **all 8 regress** |
| different conformer, A already failed `indep` | 4 | no change — already broken with the lever OFF (HEJXIF, WIWRIE, NOMMOU, QIDKUL) |
| same conformer returned | 8 | no change — the lever changed nothing it could change (AROHIA, LIYXEY, ODEWID, QESRUE, XUPTAF, YIYGAP, **YIZHIY**, GAVSED) |

**YIZHIY is explained, not excused:** the lever returned a bit-identical structure (same sha,
same `clash_vdw`, same `worst_overlap`), so no outcome could differ. Its probe classification
(`c@0, s@1`) is the smallest possible non-zero gap, too small for the pool to diverge.

**The sharp form of the result:** when the lever changed the returned conformer, it *never once*
found a different conformer that still passed independent re-perception — **8 of 8**. That is
not a random quality drift; the strict test is doing real work, and the conformers it was
holding out for are the only ones that satisfy it.

#### Why G3 and G2 can both be true

`sha256(smiles_2)` is byte-identical between arms while `worst_overlap` moves. So the lever
returns a **different conformer that emits the same OIN string**. The notation is untouched;
the geometry underneath it is not. Both gates are measuring honestly and they are measuring
different things.

#### ⚠ `indep` is NOT a pristine oracle — read `indep=False` carefully

This belongs next to the number above, because `indep=False` is easy to over-read as "the
structure is wrong":

- **`ODEWID_comp_0` passes `indep` on the key but NOT on bytes.** Full re-perception renumbers
  donor slots (`N{3}` ↔ `N{1}`) on a molecule that is otherwise a clean round trip. Scoring the
  independent arm on bytes would have manufactured a failure.
- **`AROHIA_comp_0` has the cheap and full encoders disagreeing in the *opposite* direction**
  (`pass=False indep=True`) — full re-perception says it round-trips while the harness's own
  scoring path says it does not.

So the full `XYZToSMILES().convert()` path has its own presentation instabilities and its own
disagreements with the cheap path. `indep=False` means "the two encoders disagree about this
structure," which is strong evidence of a real difference but is not proof the geometry is
wrong.

### 4.4 G4 — guard population

*(pending)*

## 5. Recommendation

*(pending all four gates)*

---

## Reproduce

```
# 22-molecule cohort, both arms, all of G1+G2+G3 in one pass
systemd-run --user --unit=l2-ab22 -p OOMPolicy=continue -p MemoryMax=8G \
  --working-directory=<worktree> --setenv=PYTHONPATH=<worktree>/src \
  <main-checkout>/.venv/bin/python tools/ab_accept_scored.py \
    --cohort spec/handoffs/v0.4.7/cohort22.json \
    --out spec/handoffs/v0.4.7/runs/ab22.json --timeout 300 --hard-cap 500 --workers 2

# determinism control: re-run arm A alone on the molecules whose sha differed
... tools/ab_accept_scored.py --cohort <subset>.json --single-arm 0 --label A-prime --out ...

# render the gate tables without re-running anything
python tools/promote_gate_report.py --ab spec/handoffs/v0.4.7/runs/ab22.json \
    [--control spec/handoffs/v0.4.7/runs/aprime.json]

# PREFILTER_VETO
python tools/probe_accept_gap.py <AROHIA.xyz> --capture-veto --json ...
```

**Cohort provenance — 21/22 exact.** The original cohort JSON was not committed; it was later
recovered (`spec/handoffs/v0.4.7/rescued/gap_cohort_ORIGINAL.json`) and diffed against the
reconstruction used by this lane's run. Membership differs by exactly one substitution: the
original has `YENDUS_comp_0`, this lane's has `HIDCIH_comp_1` (the molecule the cProfile writeup
used to generate the hypothesis). The `eta`/`control` stratum labels differ on six molecules,
which is cosmetic — `ab_accept_scored.py` concatenates both lists and treats them identically.
`YENDUS_comp_0` is covered by the rescued replicate in §4.1, so between the two runs all 23
molecules are measured.

**Budget note.** `--hard-cap 500` kills molecules the prior run let run to 450s, applied
identically to both arms. Absolute totals here are not comparable to the prior run's; the A-vs-B
comparison is.
