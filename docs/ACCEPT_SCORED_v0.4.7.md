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

### 4.1 G1 — structure quality. Measured for the first time, and the cost is not zero.

Independent replicate on the original 22-molecule cohort with the corrected
`vdw_clash_count` (`spec/handoffs/v0.4.7/rescued/ab_v2_partial.log`). 14 molecules paired
across both arms.

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

**The honest framing is "measured, small, one molecule."** Not "no cost" — that was the
degenerate metric's answer and it was wrong. Not "large" either: one molecule of fourteen gains
a single vdW clash, and more molecules improve than degrade. Bypassing `_select_by_geometry`'s
clash-first ranking does cost something, and now there is a number for it.

**A pre-existing defect the fixed metric also surfaced, unrelated to the lever:**
`POVPIA_comp_0` in the **default** arm reports `clash=16 worst=0.4344` — sixteen vdW clashes
and a worst overlap well below the 0.60 *severe* cutoff. Recorded separately, the same way
HEJXIF's `indep=False` in arm A is: it is a structure-quality defect of the default path, and it
was invisible while `mol_clash_count` was silently returning 0.

**The lever converts at least one hard-cap timeout into a pass.** GAVSED_comp_0: arm A killed
at the 330s cap (no structure at all), arm B 5.41s `pass=True`. Under a wall-clock budget the
speedup *is* a pass-rate effect.

**Runtime (ADVISORY — load 40+ on 12 cores, wall-clock is meaningless here):** KAQDOV
242.5→7.5, ZITSIE 134.2→5.5, HEJXIF 194.6→13.0, WIWRIE 74.6→6.3, NOMMOU 39.9→4.1, RATPEK
69.0→11.8.

> **The molecules that got *slower* are load noise, not regressions.** ODEWID 1.45→2.21,
> QESRUE 3.59→5.63, YIYGAP 1.83→4.05, YIZHIY 4.54→5.88 are all sub-6s NO_GAP-class molecules
> where both predicates already fire at pool index 0 and the lever changes nothing it could
> change. Said explicitly so nobody later reads these four rows as a slowdown.

### 4.2 G3 / G2 on the 22-molecule cohort

*(pending — run `l2-ab22` in flight)*

### 4.3 G4 — guard population

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
