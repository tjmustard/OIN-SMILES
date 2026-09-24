# v0.4.18 selection lane — `OIN_ETA_RETARGET`: ask each conformer whether its donors arrived

> **Harness A/B over all 1,146 η-bound molecules: self-consistent +49 / −4 → 86.10 → 87.00%;
> VERIFIED +43 / −5 → 77.66 → 78.42%.** Excluding rows on the 300 s budget boundary: +44 / −4 and
> +40 / −5. Runtime flat. **PROMOTED 2026-09-23 (owner delegation, §6); the release sweep follows.**

Branch `research/v0418-eta-selection` (off the promotion tip of `research/v0418-eta-detached`),
2026-09-20. Owner: *"Promote now and do the selection as the next lane."* Follows
`L2_ETA_DETACHED.md`, whose §5b ended on the finding this lane starts from: the two η target policies
(`target = scale × S`, shipped; `target = S`, `OIN_ETA_TARGET_UNSCALED`) **fix and break different
molecules** — a per-molecule oracle over the A/B's three arms is +222 verified with no losses, against
+150 for the better single configuration.

## 1. What the code does today (read before anything was written)

- **"First acceptable conformer wins" is the `accept_fn` early exit** (`generator3d/__init__.py`,
  `_file_and_maybe_stop`): the first conformer whose independent re-encode matches the round-trip key
  ends the pool, which is returned as a list of **one** — no sort, no selector. The test is the round
  trip and nothing else. 18 of the promotion's 39 verified losses pass it byte-exact.
- **The FF scan never checks that a binding group reached its target.** `abs_delta < d_converge` reads
  a *stalled* group as converged; a vdW-reverted step keeps `ff_success = True`. A conformer whose ring
  never arrived enters the pool as a success.
- The selector sees `(parsed, mols)` only — energy, contract mol, clash, fit, re-encode. Not the scale,
  not the target policy, not the input xyz. `OIN_ATTACH_RETURN`'s three exits are dead code (held off
  in v0.4.15: 0 gains of 289, *because every conformer in the pool was detached*).

So a selector has to judge coordination **from the string**, and the cheapest honest signal is one the
scan can compute about itself.

## 2. The bound — `tools/v0418/selection_sim.py` (offline, and it *can* express a loss)

The L2 A/B left three generated structures per η molecule, each with a known verdict. The simulated
predicate sees only what the generator sees: per element, the binding atoms `smiles_1` **declares**
(`{n}`, `{n>}`, `{n<}`), against the atoms of the candidate inside the encoder's contact cutoff
(covalent radii + 0.45 Å). Control: declared == the *input's* own contacts on 489 of 518 verified
passes (94.4%).

| selector (ties → the shipped default) | verified of 1,146 | vs shipped |
|---|---:|---|
| shipped (`on`) | 668 | — |
| **deficit only**, over on / on3 / off | **723** | **+57 / −2** |
| deficit only, over on / on3 | 712 | +45 / −1 |
| deficit, then **surplus** | 708 | +60 / **−20** |
| oracle | 740 | +72 / 0 |

**Surplus is a trap**: penalising contacts the string does not declare costs 20 losses — BOUNDARY
contacts are the modal state of a *passing* molecule (v0.4.16). Of the 72 molecules of oracle headroom
the deficit prefers the verified arm on 56, ties on 4, prefers the wrong one on 8. A bound, not the
A/B: one pool is not three separate runs.

## 3. The lever

`ff_clean` now ends every successful clean by counting, **per declared binding atom** (index-based —
the complex itself says which atoms bind), how many sit outside the contact cutoff:
`cleaner.last_arrival_deficit`. An observation; it changes nothing on its own.

With `OIN_ETA_RETARGET` on, in a complex that has an η group, at a pool scale ≠ 1.0: a conformer that
**failed to clean, or left a donor out**, is cleaned **once more from the same embedding** with the η
target at S (`clean_geometry(..., eta_unscaled=True)`), and the one with more donors in is kept.

- One conformer per attempt either way ⇒ the attempt / seed sequence is the shipped one.
- A conformer whose donors all arrived is never touched.
- A non-η molecule never reaches it: the *adapter* passes `ff_params["eta_retarget"]`, and only for a
  string with an η group.
- Lever OFF: 103 of 103 probe structures byte-identical to `results-v0.4.18-sweep`.

## 4. Probe (a sample — it chose what to run, nothing more)

63 oracle-headroom molecules + 40 verified controls, fresh process per (molecule, arm), honest round
trip: headroom **22 → 52** (+31 / −1), coordination intact 34 → 57, Σ time 956 → 759 s; controls
40 → 40, 38 of 40 structures byte-identical.

## 5. The harness A/B — all 1,146 η-bound molecules, commit `2f0b68d5`

`tools/v0418/run_selection_ab.sh`. **One arm was run.** The OFF arm is the v0.4.18 sweep of record's own
rows for the cohort (symlinks): shipped defaults *are* that sweep, and generation has been
byte-deterministic at this load three times over (L2's OFF arm == the v0.4.17 sweep 1,063 / 1,063; the
v0.4.18 sweep == L2's ON arm 1,068 / 1,068; this lane's lever-OFF probe == the v0.4.18 sweep 103 / 103).
⚠ That makes the report's "noise floor" line circular here — it prints 0 by construction.

| control | |
|---|---|
| completeness | 1,146 / 1,146 |
| input side | `smiles_1` differs OFF vs ON **0** |
| lever fired | **158** of 1,067 structures differ — it touches what it is aimed at and little else |

| | pass OFF → ON | gains | losses | net | projected on n = 5,000 |
|---|---:|---:|---:|---:|---|
| **self-consistent** | 825 → 870 | 49 | 4 | +45 | 86.10% → **87.00%** |
| **VERIFIED** | 669 → 707 | 43 | 5 | +38 | 77.66% → **78.42%** |
| excluding budget-boundary rows | | 44 · 40 | 4 · 5 | **+40 · +35** | +0.80 · +0.70 pts |

- **Budget boundary, stated rather than hidden.** 5 self-consistent gains (3 verified) were `hard_fail`
  at the 300 s budget in the sweep and finished here; 5 other molecules went *into* `hard_fail`
  (`structural → hard_fail`, not passes lost). With a built OFF arm those cannot be told from timing —
  the same rows (`TEQHOM`, `WOFGUT`) flipped between the sweep and L2's A/B. Book the lever at the
  excluded row: **+0.80 / +0.70 pts**.
- Every one of the 49 gains and 4 / 5 losses has a **different structure** in the two arms — none is a
  re-scoring of the same coordinates.
- Gains by census fault (verified): `DETACHED` 21 · `NONE` 13 (the promotion's own losses coming back)
  · `G_HCOUNT` 2 · `SPHERE_DIFF` 2. By metal: Ru 8 · Ir 7 · Rh 6 · Mo 5 · Zr 5 · Ti 3.
- Losses (verified, 5): `MERCIV`, `OGONOP`, `UKUMUJ` still round-trip byte-exact and fail only the
  ruler; `SARKEC`, `ZOSYEN` → `key_equal`. Ru 3, Pt 1, Mo 1.
- Structures: coordination intact 878 → **916**; DETACHED class 168 → 191 of 295, atoms in band
  60.5% → 68.7%.
- Runtime flat: Σ `elapsed_s` 11.63 → 11.81 h; > 30 s 238 → 243.

The offline bound said +45 / −1 for a choice between these two policies; the harness measured +43 / −5.

## 6. Promoted — owner delegation 2026-09-23

The owner delegated the rest of the v0.4.18 plan ("follow best practices … use as much compute as
needed"); the recommendation to promote had stood since 2026-09-20. `OIN_ETA_RETARGET` joined
`_DEFAULT_ON` at `cdb76eca` with its evidence block. Suite 1,076 green; ARM 1 byte-identical with
the lever on and at `"0"` (0 of 62 rows). The release sweep it owes is `results-v0.4.18-release-sweep`
(expected 87.00 / 78.42 up to budget-boundary rows); `/refreeze-goldens` follows it, and ARM 2 now
carries the honest observation column L2 §5e asked for — it will see this promotion where the gated
field cannot.

### What was on the table


| option | self-consistent | VERIFIED | owes |
|---|---|---|---|
| leave it OFF | 86.10% | 77.66% | nothing |
| **promote `OIN_ETA_RETARGET`** | **~87.0%** (+0.80 … +0.90) | **~78.4%** (+0.70 … +0.76) | full sweep, ARM 1 check, `/refreeze-goldens` |

What is left of the oracle's +72 after this: ~30 molecules, where the deficit ties or points the wrong
way — coordination judged per element count cannot see them. 🔴 And ARM 2 will not see this promotion
either (`L2_ETA_DETACHED.md` §5e): its re-encode goes through the generator's own bond graph.

## 7. The release sweep, and what the honest column found

**`results-v0.4.18-release-sweep` (shipped defaults, commit `cdb76eca`, 5,000 / 5,000, 6 h 29 min):
4,347 / 3,920 = 86.94% self-consistent / 78.40% VERIFIED** — the A/B predicted 4,350 / 3,921. Every
non-η structure is byte-identical to `results-v0.4.18-sweep` (3,713 / 3,713); every η structure is
byte-identical to the A/B's retarget arm (1,070 / 1,070); `smiles_1` moved on 0. All eleven verdicts
that differ from the projection are rows at the 300 s budget. Against the v0.4.17 sweep of record,
v0.4.18 is **+4.22 / +3.80 pts**.

`/refreeze-goldens`, full run: all 12 shards PASS, **0 of 425 rows owed**; audit 421 SAME + 4
sentinels; the lever fired on 41 of 198 η rows (structure differs from the previous full run) and
the gated field moved on 0 — as at the first promotion.

**The honest column** (`tools/gate_arm2_roundtrip_one.py`, field 9; `tools/v0418/arm2_honest_column.py`):

| 401 measured rows | honest pass | honest fail |
|---|---:|---:|
| gated pass | 243 | **135** |
| gated fail | 3 | 20 |

Validated: on the 370 rows whose gate structure is byte-identical to the release sweep's, the column
agrees with the sweep's honest verdict **370 / 370**. The 135 are `structural` 95, `key_equal` 21,
`hard_fail` 16, `facmer_divergent` 2 in that sweep, plus ULODUU, a golden-only row outside the
5,000-molecule cohort (the first draft of this table typed
134 against its own 401-row denominator; the instrument's output says 135). **ARM 2 passes 36% of its passing rows on the
generator's own bond graph alone** — the scored-vs-honest gap v0.4.8 measured on the whole cohort,
on cohorts chosen for being hard. Still an observation (the gate reads columns 1–3; the goldens do not
carry it); gating on it would re-freeze ~135 rows as honest fails or drop them — a change to what ARM 2
*means*, left to the owner with the evidence frozen.

## 8. Not done

- No two-policy *pool* (a third product dimension): the re-clean keeps one conformer per attempt, which
  is what made "shipped conformers are never touched" true by construction.
- No arrival term in the selector's sort or in `accept_fn`: a conformer that key-matches has its donors
  in by the encoder's own perception, so there was nothing for it to decide.
- `selection_sim.py`'s per-element predicate and the lever's per-atom one are different instruments;
  the second is stricter and index-based. The simulation bounded the lane; it did not tune the lever.
