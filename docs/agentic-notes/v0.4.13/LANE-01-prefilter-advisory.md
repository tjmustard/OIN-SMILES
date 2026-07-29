# v0.4.13 Lane 1 — `PREFILTER_VETO`: the defect is real, the prevalence is not measured

> **Status: INCOMPLETE, deliberately.** The lever, its telemetry and its harness are built and
> committed. The defect is **confirmed** on the one molecule where the answer was already known.
> **Corpus prevalence is n = 1**, which is precisely what this project forbids quoting — so it is
> handed to v0.4.14 rather than dressed up.

Instrument: `tools/prefilter_prevalence.py` · lever: `OIN_PREFILTER_ADVISORY` (default OFF)

---

## 1. What the lane was for

`_reencode_key_matches` (`generation/metallogen_adapter.py`) is two-stage. Step 1 is a **cheap
prefilter**: re-serialize the generated geometry through the generator's *own* contract-mol
connectivity and reject on a key mismatch, justified in-code as *"a MISMATCH here is a reliable
'geometry is wrong' signal"*.

`AROHIA_comp_0` falsifies that claim: the cheap test matches **0 of 48** conformers while the
strict independent test matches **16 of 48**. Because the cheap `return False` fires first, those
16 are unreachable **in both arms of every A/B ever run on this molecule**. The error runs in the
**pessimistic** direction — it makes the project look worse than it is — which is why it survived
five releases unexamined.

**Half the original charter was already dead on arrival.** The v0.4.13 sketch framed this as a
*scoring* defect. `OIN_INDEP_SCORE` went default-ON in v0.4.8 and the harness no longer scores
with the cheap path; on the frozen corpus cheap-fails-but-independent-passes is **28/5000** and
the honest metric already counts every one correctly. Only the *acceptance* half survives.

## 2. What was built

`OIN_PREFILTER_ADVISORY`, default OFF (so the default path is byte-identical): a cheap veto falls
through to the strict test instead of returning `False`. The decision is made **observable**,
which is the part that matters:

| counter | meaning |
|---|---|
| `adapter.prefilter_veto_overridden` | cheap NO, strict **YES** — a conformer the default discards |
| `adapter.prefilter_veto_confirmed` | cheap NO, strict NO — the prefilter was right |
| `adapter.prefilter_cheap_pass` | cheap YES — the veto never had an opinion |

`tools/prefilter_prevalence.py` runs both arms per molecule and emits one of four verdicts, which
exist because **they would otherwise all read "0"**: `INSTRUMENT_DEAD` · `NO_POPULATION` ·
`PREFILTER_VINDICATED` · `DEFECT_CONFIRMED`.

## 3. 🔴 A lever-interaction bug this change would have introduced

The `OIN_ACCEPT_SCORED` branch accepts *"on the predicate the SCORE uses"* — and the score **is**
the cheap re-encode. Before this lever, a cheap mismatch had already returned `False`, so reaching
that branch *implied* the cheap test passed. Falling through instead means the two levers
**together** would accept a conformer the score itself calls a failure.

Guarded with `and not cheap_vetoes`. **Neither lever's own A/B would have exercised the
combination** — this is the class of defect that only appears when two independently-gated changes
meet, and the reason each lever's charter has to name the others it can interact with.

## 4. What the fixture says

`AROHIA_comp_0`, the two-point gate (the answer is known independently, so a zero here means the
lever is not wired). **Scored honestly** — the distinction matters, see §5.3:

| | lever OFF | lever ON |
|---|---|---|
| `prefilter_veto_overridden` | — | **2** |
| `prefilter_veto_confirmed` | — | 1 |
| how acceptance ended | `early_exit_**miss**` | `early_exit_**hit**` |
| what was returned | `embed.best_rejected_returned: 11` — a **previously-rejected** conformer, via the geometry fallback | a conformer the strict test **accepted** |
| round-trips (honest) | **True** | **True** |
| elapsed | 22.73 s | **2.52 s** |

**Read the `recovered = 0` correctly: it is 0 because BOTH arms already pass, not because both
fail.** The accuracy delta on this molecule is genuinely nil.

**What does change is how the pass is obtained.** With the lever off, `accept_fn` accepts
*nothing* and `_select_by_geometry` hands back a conformer acceptance had already rejected — the
same unguarded-fallback path that Lane 2 measures at **280 molecules corpus-wide** (the GAVSED
class). With the lever on, one of the two overridden conformers is *accepted on its merits*, and
the molecule finishes **9× faster**.

So on n = 1 the lever converts *"accept nothing, return a reject and hope"* into *"accept a
conformer the independent test endorses"*, at the same honest verdict. Whether that is worth
anything at corpus scale is exactly the unmeasured question.

## 5. Three reasons these numbers are not what they look like

1. **The denominators are not pool sizes.** The first override *accepts*, which stops the pool
   filling, so the lever-ON arm evaluates far fewer conformers than the pool would hold. AROHIA's
   documented 0/48-vs-16/48 was measured with the pool **forced full**
   (`tools/probe_accept_gap.py`). "3 vetoes" here is the same defect observed *until it stopped
   mattering*, not a smaller one.
2. **The −20 s latency delta is not a speed claim.** It is early exit, and it was taken on a loaded
   machine besides. v0.4.12 measured `UQUXAG_comp_0` at 17.93 s loaded vs 11.06 s clean — a 62%
   inflation, comparable to the whole effect — and discarded the loaded run.
3. **🔴 The first version of this tool scored with the CIRCULAR predicate, and it INVERTED the
   answer.** It used `get_oin_string(res.mol, coords)` — the generator's own bond graph, the exact
   thing v0.4.8 replaced — copied from the older A/B tools in `tools/`. Uniquely wrong *here*,
   because this lane's entire subject is a cheap-vs-strict disagreement: scoring the outcome with
   the cheap predicate judges the lever by the very test it exists to override.

   **It did not merely blur the result — it reversed it:**

   | scored with | lever OFF | lever ON |
   |---|---|---|
   | `get_oin_string(res.mol, coords)` — circular | `passed=False` | `passed=False` |
   | `XYZToSMILES().convert(<written xyz>)` — honest | **`passed=True`** | **`passed=True`** |

   Identical telemetry, opposite verdict. The circular predicate calls `AROHIA_comp_0` a **double
   failure** when independent re-perception says it round-trips in both arms — a live, single-
   molecule instance of the 8 false negatives v0.4.8 measured, reproduced by accident while
   building something else. Anything written from the first run would have said "the molecule
   fails either way"; the truth is "it passes either way, by different routes".

## 6. Why this stops here

The acceptance bar says *"corpus prevalence stated with n and the command. **Not n = 1.**"* It is
n = 1. The measurement needs a **quiet machine** — its deliverable includes a latency cost, and
this session's machine was running mirror audits, a 62-fixture gate and an 11-molecule ARM 2
re-freeze.

**What v0.4.14 must run**, in this order:

```bash
# 1. the wiring gate — a zero here invalidates everything after it
OIN_PREFILTER_ADVISORY=1 $V tools/prefilter_prevalence.py \
    --xyz <main>/tmCAT-tmPHOTO_xyz_dataset/cat/AROHIA_comp_0.xyz

# 2. the cohort, RE-DERIVED from the frozen sweep — never a pre-v0.4.8 cohort
$V tools/prefilter_prevalence.py --cohort <re-derived> --out prefilter_prevalence.json
```

⚠ **Any cohort frozen before v0.4.8 must be re-derived, not reused.** v0.4.12 pointed a pilot at
the v0.4.6 accept-gap cohort and got a flat A/B because **all 8 of its molecules** now satisfy the
key inside `accept_fn`.

**And if the answer is small, say so and do not fix it.** `PREFILTER_VINDICATED` is a real verdict
the tool can emit. Four of this project's releases have ended by refuting their own plan; a fifth
would be unremarkable.
