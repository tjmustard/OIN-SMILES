# v0.4.14 Lane 2 — `PREFILTER_VETO`'s corpus prevalence

**Inherited from v0.4.13, already committed and green:** `OIN_PREFILTER_ADVISORY` (default OFF),
three telemetry counters, `tools/prefilter_prevalence.py`, 64 acceptance-path tests. v0.4.13 built
the instrument and measured `n = 1`, which this project forbids quoting. **This lane is the
measurement.**

```bash
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
V=$MAIN/.venv/bin/python
export PYTHONPATH=$PWD/src
```

---

## 1. Wiring gate — a zero here invalidates everything after it

```bash
$V tools/prefilter_prevalence.py --xyz $MAIN/tmCAT-tmPHOTO_xyz_dataset/cat/AROHIA_comp_0.xyz
```

| counter | v0.4.13 | v0.4.14 |
|---|---:|---:|
| cheap veto, strict **ACCEPTS** (`overridden`) | 2 | **2** |
| cheap veto, strict rejects (`confirmed`) | 1 | **1** |
| cheap veto total (`advisory` fired) | 3 | **3** |
| cheap PASSED, veto never consulted | 0 | **0** |
| verdict | `DEFECT_CONFIRMED` | **`DEFECT_CONFIRMED`** |

**The instrument is alive**, and reproduces v0.4.13's reading exactly. `overridden` is quoted only
beside `confirmed` and `cheap_pass`, because a lever that never fires and a lever that fires and
finds nothing print the same counter.

⚠ **`overridden = 2` is not "2 of 48".** The first override *accepts*, which stops the pool
filling, so the ON arm evaluates far fewer conformers than the pool holds. AROHIA's documented
`0/48` cheap vs `16/48` strict was measured with the pool **forced full**
(`tools/probe_accept_gap.py`). A count of 2 is the same defect observed *until it stopped
mattering*, not a smaller one.

⚠ Wall clock read **OFF 22.6 s → ON 2.5 s, delta −20.1 s**. That is **early exit, not speed**: the
ON arm accepts a conformer the OFF arm rejected and therefore stops generating. It is not a
performance claim and must not be quoted as one.

This gate was run under load (Lane 1's transition sim was executing). That is legitimate *for a
counter-based gate* — the counters are exact regardless of contention — and illegitimate for the
latency figures above, which is why the cohort run below is the one taken on a quiet machine.

---

## 2. The cohort — re-derived, never reused

⚠ **Any cohort frozen before v0.4.8 must be re-derived.** v0.4.12 pointed a pilot at the v0.4.6
accept-gap cohort and got a flat A/B because **all 8 of its molecules** now satisfy the key inside
`accept_fn`.

Population at risk = molecules the frozen sweep records as a genuine failure, since the defect's
direction is **pessimistic**: a wrongly-vetoed conformer can only turn a pass into a failure.

| bucket | population | returned a structure |
|---|---:|---:|
| `hard_fail` | 319 | 22 |
| `structural` | 417 | 417 |
| **total** | **736** | 439 |

Molecules producing **no** structure are kept in the population, not filtered out: if the cheap
prefilter vetoes every conformer, "produced nothing" is precisely the symptom.

### Why a sample, and how big

The full population costs **73.2 CPU-hours** across both arms (`hard_fail` median 300.3 s — the
timeout — and `structural` median 71.7 s, doubled). That is the same order as the 55 CPU-h sweep
v0.4.13 declined to run, for a measurement that does not need it.

**Stratified random sample, `n = 50`, seed 7**, proportional to the strata:

```
hard_fail    population  319 (43.3%)  sampled 22
structural   population  417 (56.7%)  sampled 28
cohort: tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.14-prefilter   50 files, 0 dangling symlinks ✓
```

Selection is recorded in `COHORT.json` (seed, strata sizes, the molecule list) so the draw is
reproducible rather than described.

```bash
$V tools/prefilter_prevalence.py \
    --cohort $MAIN/tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.14-prefilter \
    --out $MAIN/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane2/prefilter_prevalence.json
```

---

## 3. What this measurement can and cannot support

Stated **before** the numbers, so the framing is not chosen to fit them.

- **Prevalence — supported.** `overridden > 0` on a molecule is a direct observation inside the ON
  arm. It is not an A/B and no stochastic comparison enters it.
- **Recovery count (`pass_off → pass_on`) — CONFOUNDED.** It compares two independent runs of a
  **stochastic** generator, which is this project's standing trap: *never A/B by re-running a
  stochastic harness.* It is reported, and it is not causal evidence on its own.
- **Latency delta — direction only.** A negative delta is early exit. Only the quiet-machine run's
  numbers are quotable; v0.4.12 measured `UQUXAG_comp_0` at 17.93 s loaded vs 11.06 s clean (62%
  inflation, comparable to the whole effect) and discarded the loaded run.

**A `PREFILTER_VINDICATED` result is a real result and ships as one.** Four of this project's
releases have ended by refuting their own plan. If the population is small, this lane says so and
does **not** ship a fix.

---

## 4. Result

*(cohort run pending — see §5 of the release close-out)*
