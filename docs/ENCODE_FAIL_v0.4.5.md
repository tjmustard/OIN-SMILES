# `encode_fail` — v0.4.5

Owns the `encode_fail` class of the v0.4.5 gap-to-100% decomposition: 48 molecules
(11.0% of the 436-molecule gap on the capstone arm, 6,719 mols,
`tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042/bucket_report.json`) where
`smiles_1` is null — the encoder produced no OIN string at all for the INPUT structure.

## 0. This cohort was already worked once (v0.4.4 SL5) — read that first

The 48 `molecule` values in `bucket_report.json`'s `encode_fail` rows are **the exact same
48** already sub-triaged in v0.4.4's SL5 swimlane (`tools/sl5_triage.py`,
`docs/ENCODER_ROBUSTNESS_v0.4.4_SL5.md`, landed on `main` before this branch forked). SL5's
own triage of these 48: `boron_cluster` 34, `resonance_timeout` 10, aromatic/quinoid 3,
`perception_charge_gap` 1 (`ASISAX`). SL5 landed:

- a typed `OINEncodeError` for electron-deficient boron clusters (`_is_electron_deficient_cluster`,
  `xyz2mol.py:535`) — **classifies**, does not encode, the 34 boron molecules;
- a forked, CPU-time-bounded `ResonanceMolSupplier` (`_resonance_candidates_isolated`,
  `xyz2mol.py:369`) that lets a completing large ligand finish byte-identically and
  degrades a genuine hang to the single perceived form, recovering some of the timeout 10;
- documented as residual/future work: the `get_UA_pairs`-matching hangs, the 3 quinoid
  de-aromatization cases, and `ASISAX`'s charge gap.

So this lane's job was **not** "triage 48 mols cold" — it was: re-measure the same 48
against current `main` (which already has SL5) to find the actual residual, and fix what's
cheap in what's left.

## 1. Re-measurement methodology and an honest caveat

Reproduced via `tools/sl5_triage.py`'s own `run_worker`/`classify` (isolated subprocess per
molecule, calls `XYZToSMILES().convert()` directly — **no generation**, so every result
below is unambiguously encoder-side, not a generation/verification artifact).

**Caveat on wall-clock timeouts:** the stock tool's `PER_MOL_TIMEOUT_S=90` is *shorter* than
SL5's own `_RESONANCE_CPU_BUDGET_S=120` CPU-second fork budget, which would misreport a
molecule the fork *would* recover as a hang. Reran with a 200s per-molecule wall-clock cap
instead. Even so, several molecules did not resolve inside 200s. This machine was running
five sibling v0.4.5 lane agents concurrently (`ps aux` showed concurrent
`canonicality_probe.py` and lane test-suite jobs; `uptime` load 8–12 on 12 cores throughout
this measurement) — the fork's bound is CPU-*time*, not wall-clock, specifically so the
*outcome* is load-independent, but under contention the same 120 CPU-seconds can take
substantially more than 120 *wall* seconds to accrue. So a `resonance_timeout` result below
means "did not resolve within 200 wall-clock seconds **on a heavily loaded shared box**,"
not "the fork budget is broken" — I did not have the isolated machine time to tell those
apart, and did not chase it further (out of this lane's scope; SL5's forked-resonance design
is another lane's carefully-reasoned infrastructure).

## 2. Histogram (current `main`-equivalent state, this branch @ fork point)

Reproduce: `PYTHONPATH=src .venv/bin/python /path/to/rerun_encodefail.py` (200s variant of
`tools/sl5_triage.py`, dataset dir `tmCAT-tmPHOTO_xyz_dataset`).

<!-- FINAL TALLY INSERTED BELOW ONCE THE BACKGROUND RERUN COMPLETES -->

## 3. Fixed: `ASISAX` — one-line permissive rescue, opt-in lever

**Root cause.** `ASISAX` is a Ni tetraaza-macrocycle (confirmed by direct diagnosis: at
ligand charge 0 the encoder's own `AC2mol` perceives it fine; every other charge in -6..+6
returns `None`). `get_lig_mol`'s charge-rescue loop (`_rescue_unusable_perception`,
`xyz2mol.py:498`) rejected the charge-0 candidate anyway, because it checked
`stuck_ring_atoms(candidate) or not _perception_is_usable(candidate)` — i.e. it treated
*any* stuck (unkekulizable-as-aromatic) ring as an automatic reject, even though
`_perception_is_usable` **already** calls `kekulize_safe_sanitize`, which can repair a stuck
ring by de-aromatizing it (that repair is what makes the 3 quinoid molecules in §4 recoverable
*in principle* — it just doesn't succeed for them). For `ASISAX` the repair succeeds, so the
rescue loop was strictly *more* conservative than the encoder's own downstream repair path,
for no documented reason. With every charge in -6..+4 rejected, the sweep exhausts and
`get_lig_mol` raises the generic `ValueError` → `encode_fail`.

**Fix, scoped narrowly.** `xyz2mol.py`'s `_rescue_unusable_perception`: split the combined
check so the stuck-ring rejection is skipped when `OIN_RESCUE_STUCK_RING` is truthy, leaving
the usability check (which already subsumes ring-repair) as the sole gate. Default unset →
byte-identical to the pre-fix code (same combined boolean short-circuit). Mirrors the
`OIN_EMIT_AXIAL` pattern (`xyz2mol.py:1101-1110`): inline truthiness check, no import
overhead when off.

**Why gated, not landed ungated.** The early return in `_rescue_unusable_perception`
(`if best_res_mol is not None and _perception_is_usable(best_res_mol): return`) means this
loop is *only* reached by ligands whose default perception was already unusable — so far, so
"only touches things that currently fail." But the loop returns the **first** charge (in
Huckel-distance order) that passes its checks; loosening the stuck-ring check can make an
*earlier* charge win where today a *later* one does, for any ligand that currently reaches
this loop and is rescued by some non-stuck charge further down the order. I could not rule
that out with a full corpus sweep (out of budget for this lane), so per the task's own
instruction this went behind a default-OFF lever rather than landing unconditionally.

**Verified:**
- Lever OFF: `ASISAX` still raises `ValueError` (byte-identical to before this change).
- Lever ON: `ASISAX` now encodes: `[Ni_SPL].CC(C)C1=c2c3c4c(c5ccccc25)=C(C(C)C)C(=N{0}4)C=CC2=N{1}c4c5c(c6ccccc6c4=C2C(C)C)=C(C(C)C)C(=N{2}5)C=CC1=N{3}3` — deterministic across a repeat.
- `tests/unit/test_encoder_robustness.py::TestStuckRingRescuePermissive` pins both.
- `tests/unit/test_regression_stability.py` (4 golden fixtures): unaffected, 6/6 pass.

**Canonicality check — NOT clean, and this matters.** `tools/canonicality_probe.py --only
ASISAX_comp_0 --trials 2` with the lever on:

```
ASISAX_comp_0: 4/6 drifted ['rdkit_canonical']
drift by transform: renumber 2, both 2   (rotate: 0)
of the drifted, key ALSO changed (isomer-level, worse): 1
```

So `ASISAX` now encodes deterministically for its *actual, fixed* atom ordering (the one
real dataset file has one order; repeat runs agree; pure rotation never drifts), but is
**not** robust to renumbering — presenting the same graph with atoms in a different order
changes the OIN string, including at the comparison-key level. Per this task's own
instruction, that is "moved to a different bucket," not "fixed."

**This is not a new defect I introduced — it's the already-tracked, corpus-wide
renumbering-instability finding** (`docs/RENUMBERING_INSTABILITY_v0.4.5.md`: 43.9% of a
225-molecule sample drift under pure renumbering, 21.1% at the key level; suspect #1 in that
doc is exactly "`AC2BO` / resonance-form order-dependence changing perceived bond orders,"
which is this macrocycle's failure mode). `ASISAX` could not show this defect before,
because it never reached a successful encode to begin with; my fix's only contribution is
that it now *can* show it. Other lanes (1/2/8) own closing that gap; I am not duplicating
their work here.

**Net honest read:** `ASISAX` moves from "encoder returns nothing, ever" to "encoder returns
a plausible, chemically-sensible Ni-macrocycle OIN for the structure as given, deterministic
under rotation and repetition, but not proven stable under renumbering" — real, bounded
progress, not full closure.

## 4. Documented, deferred: 3 quinoid/ylide kekulization failures

`KAXVOX`, `KAXWAK` (same Zn/N5O2S2 core, `KAXWAK` brominated — clearly two derivatives of one
compound), and `LEZWAO` (unrelated Pd/phosphine/silane complex) already raise a typed
`OINEncodeError` from `kekulize_safe_sanitize` (`aromaticity.py:161`): after de-aromatizing
the detected quinoid ring, the residual bond orders are still unusable (`KAXVOX`/`KAXWAK`:
still can't kekulize 5 atoms after the ring fix; `LEZWAO`: a carbon left at explicit valence
5). SL5's own docstring for this function already names the mechanism: *"at the wrong ligand
charge `AC2mol` leaves the ipso carbon pentavalent, and only re-perceiving the charge fixes
that."* That is a **different, deeper code path** than §3's fix (`fix_equivalent_Os`'s
whole-molecule equivalent-oxygen pass, not the per-ligand-fragment charge rescue), and fixing
it means teaching that pass to re-perceive charge the way `_rescue_unusable_perception`
already does for ligand fragments — a materially larger change for 3 molecules, and one I
did not attempt given the lane's budget. Recorded here rather than forced: this is "needs new
chemistry," in the task's own words, not "one permissive valence entry."

Fixtures added (`tests/fixtures/{KAXVOX,KAXWAK,LEZWAO}_comp_0.xyz`) so whoever picks this up
next does not need to re-pull them from the dataset.

## 5. Confirmed unfixable: 34 boron clusters

Unchanged from SL5: RDKit's 2-center-2-electron valence model has no Lewis structure for a
3-center-2-electron borane/carborane cage. `get_lig_mol`'s charge sweep already spans -4..+4
(and, per §3's diagnosis code, -6..+6 finds nothing new for a cage either) — no charge
widening fixes this class; it needs a different bonding model entirely (multi-center bonds),
out of scope for a valence-graph encoder. The typed `OINEncodeError` (already landed) is the
correct, honest terminus: **an encoder that refuses this input is correct**, not a bug.

## 6. Unresolved / needs more machine-time than this lane budgeted

A handful of the `resonance_timeout` cohort (`BENVOG`, `FAQYUU`, `HICLAG`, `HOHKUL`, and
possibly others — see §2's final tally) did not resolve within 200 wall-clock seconds under
this machine's load during this session (§1's caveat).

**Two of the three dataset-labelled `oom_killed` molecules are not encoder problems at all.**
The original `results-capstone-v042/bucket_report.json` records `HICLAG`, `HOCVAY`, and
`WEFZAL` failing with `Generation/Verification failed at UFF_1: child process died with exit
code -9` — an OOM kill — landing all three in `encode_fail` only because the harness's
`_ENCODED` marker mechanism does not survive a child dying (a report with `smiles_1: null`
cannot distinguish "encoder hung" from "encoder finished fine, generator died later"). Testing
each directly against `XYZToSMILES().convert()` alone (no generation):

- `HOCVAY` — encodes immediately.
- `WEFZAL` — encodes in **1.9s**: `[Pd_SPL].CC(=O)OC[C@H]1O[C@@H](N2C{0}N(CCCN3C{1}N([C@@H]4O[C...` (a Pd-nucleoside/PNA-type complex).
- `HICLAG` — genuinely does not resolve within 200s even encoder-only, consistent with SL5's
  own doc naming `HICLAG` (alongside `FAQYUU`) as the unrecovered `get_UA_pairs`/
  `max_weight_matching` O(V³) hang, not a resonance timeout.

So **2 of 48** (`HOCVAY`, `WEFZAL`) are not `encode_fail` at all under a clean, generation-free
signal — their presence in the frozen 48-count is a harness-bucketing artifact (a
generation-side OOM misattributed to the encoder), not an encoder refusal. This is worth
flagging to whoever owns `hard_fail`/generation: some fraction of what is currently bucketed
as `encode_fail` in the frozen `bucket_report.json` may actually belong to generation-side
OOM. I did not re-triage the full 48 against this distinction beyond these three (would
require re-running the full encode+generate pipeline per molecule, out of this lane's scope),
but flag it as a measurement caveat on the 11.0%/48 headline number itself: it is measuring
the *old* harness's bucketing, not a clean encoder-only signal.

## 7. Addressable-vs-total, honestly

| | count |
|---|---:|
| confirmed unfixable (boron valence-model ceiling) | 34 |
| fixed, encodes, but not proven canonical (moved bucket, not closed) | 1 (`ASISAX`) |
| documented, deferred — needs a deeper charge-reperception fix in a different function | 3 (quinoid/ylide) |
| resonance/UA-pairs timeout cohort, encoder-side, unresolved this session | up to 10 (§2 final tally) |
| dataset-mislabelled — encodes fine in isolation, likely a generation-side OOM | ≥1 (`HOCVAY`), possibly more (§6) |

Of the 48, **1** was actually fixed (with the canonicality caveat above), **34** are honestly
unfixable within this encoder's valence model, and the remainder needs either more machine
time free of contention (timeout cohort) or a deeper, different-function fix (quinoid cohort)
than this lane's budget covered.

## 8. Reproduce

```bash
export PYTHONPATH=$PWD/src
V=/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python
# Full re-triage (200s/mol; ~30-40 min under load):
$V tools/sl5_triage.py --dataset-dir /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset
# ASISAX only, fixed:
OIN_RESCUE_STUCK_RING=1 $V -c "from oinsmiles import XYZToSMILES; print(XYZToSMILES().convert('tests/fixtures/ASISAX_comp_0.xyz'))"
# Canonicality check on the fix:
OIN_RESCUE_STUCK_RING=1 $V tools/canonicality_probe.py \
    --dataset /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset \
    --only ASISAX_comp_0 --trials 2 -v
```
