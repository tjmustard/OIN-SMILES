# v0.4.8 · the honest round-trip baseline, N = 5000

> **`byte_exact` is 72.46%, not 82.80%.**
> The old number was measured with an instrument that asks the generator whether the generator
> was right. This is the same 5000 molecules, the same conformers, the same key and the same
> `status` gate — scored by re-perceiving each generated structure from its coordinates alone.
>
> **The drop is the deliverable, not a regression. Nothing in the library moved: the encoder is
> byte-identical on 4985/4985 inputs.**

---

## 1. The headline, both ways

| bucket | scored | % | honest | % | delta |
|---|---:|---:|---:|---:|---:|
| `byte_exact` | 4140 | 82.80% | **3623** | **72.46%** | **−517** |
| `key_equal` | 520 | 10.40% | 610 | 12.20% | +90 |
| `facmer_divergent` | 1 | 0.02% | 16 | 0.32% | +15 |
| `structural` | 9 | 0.18% | 417 | 8.34% | +408 |
| `hard_fail` | 315 | 6.30% | 319 | 6.38% | +4 |
| `encode_fail` | 15 | 0.30% | 15 | 0.30% | +0 |
| **total** | **5000** | | **5000** | | |

`key_equal` sub-split: `slot_renumber` 459 → 496, `rdkit_canonical` 61 → 114.

The two columns differ in **exactly one variable**: which round-trip string the verdict reads.

| column | string | what it is |
|---|---|---|
| scored | `smiles_2` | `get_oin_string(gen_result.mol, coords)` — the generator's own bond graph |
| honest | `smiles_2_indep` | a full `XYZToSMILES().convert()` of the generated XYZ — bonds *and* stereo re-derived from coordinates alone |

```bash
V=$PWD/.venv/bin/python; export PYTHONPATH=$PWD/src
D=$PWD/tmCAT-tmPHOTO_xyz_dataset
$V tools/honest_rescore.py --results-dir $D/results-v0.4.6-sweep \
    --output-dir $D/results-v0.4.8-honest --workers 10 --fill-coordination      # 334 s
$V tools/roundtrip_bucket_report.py --results-dir $D/results-v0.4.8-honest --score both
```

## 2. Where every corrected molecule went

"517 molecules moved" is not a reconciliation. A false positive can land in `key_equal` rather
than falling out of the passing set altogether, and those are different findings with different
owners. The full transition matrix:

| scored | → honest | n |
|---|---|---:|
| `byte_exact` | `structural` | **350** |
| `byte_exact` | `key_equal` | 180 |
| `key_equal` | `structural` | 62 |
| `key_equal` | **`byte_exact`** | **25** |
| `byte_exact` | `facmer_divergent` | 12 |
| `key_equal` | `facmer_divergent` | 4 |
| `byte_exact` | `hard_fail` | 3 |
| `structural` | **`byte_exact`** | **3** |
| `structural` | **`key_equal`** | **2** |
| `key_equal` | `hard_fail` | 1 |
| `facmer_divergent` | `structural` | 1 |
| *(unchanged)* | | 4357 |

**613 molecules degraded, 30 improved.** The corrections run in both directions, as predicted:
the shortcut asserted coordination that was not there *and* dropped stereo that was — the 5
`FAIL → pass` rows include `ZOYYUJ_comp_0`, whose `CC[S@]{2}(=O)` sulfoxide is the same
mechanism as the `YOSXIP` fixture.

### The haptic split, at corpus scale

| subset | scored `byte_exact` | lost under honest scoring | rate |
|---|---:|---:|---:|
| haptic (`{n>}` in `smiles_1`) | 897 | 329 | **36.7%** |
| non-haptic | 3243 | 216 | **6.7%** |

A 5.5× difference, same direction as the 936-molecule cohort's 28.1% / 2.8% and higher in both
arms. **More than a third of every haptic pass in this corpus was false.** Anyone validating
this metric on one subset alone would be wrong by a factor of five.

## 3. Is the honest arm itself trustworthy?

It has to be interrogated, not assumed — an encoder could simply mis-perceive generated
geometries and manufacture failures. Three independent checks say it does not.

### 3a. The control group

`report["coordination"]` is coordinate-only, ~2.2 ms, and consults **neither** bond graph. Over
the 3595 molecules that stayed `byte_exact` in both arms it flags a lost metal contact on
**48 — 1.3%**, comfortably inside its documented 3.7% false-alarm band. The honest arm does not
invent failures in the population where it agrees.

### 3b. The moved group

Over the 428 molecules that moved from a passing bucket into `structural` / `facmer_divergent`:

| `coordination` says | n | % |
|---|---:|---:|
| **FLAG — lost a metal contact** | **277** | **64.7%** |
| intact, but a contact within 0.1 Å of the cutoff | 52 | 12.1% |
| intact, clean | 99 | 23.1% |

**64.7% versus 1.3% is a 50× enrichment**, from an instrument that shares no code path, no bond
graph and no input with the one being checked.

### 3c. The two instruments agree on the *mechanism*, not just the count

Of the 277 doubly-flagged molecules, 95% show a **geometry-token change** in the honest string —
`TET→LIN` (66), `TPL→LIN` (54), `TET→TPL` (47) — i.e. the metal lost coordination sites. Ask
whether the size of that change matches the number of contacts `coordination` counted as lost:

| contacts lost ÷ coordination sites lost | n | reading |
|---:|---:|---|
| 1 | 20 | monodentate / σ-donor |
| 2 | 29 | η² alkene |
| 3 | 29 | η³ allyl |
| 4 | 40 | η⁴ diene |
| 5 | 35 | η⁵ Cp |
| 6 | 27 | η⁶ arene |
| 7 | 5 | η⁷ |

The ratio is not 1:1 — and it should not be, because a η⁵-Cp contributes five *contacts* while
occupying one coordination *site*. The distribution lands on the integers 1–7, and **95.2% of
the molecules with a ratio above 1 carry a haptic token in their input OIN** (against 40% of
those at ratio 1). Two instruments built from different data — one from the OIN geometry token,
one from raw interatomic distances — agree on which ligand left *and how many atoms it was
bound through*.

### 3d. The 99 that `coordination` calls clean

Not the honest arm misfiring — `coordination` compares per-element contact counts at the metal
and is blind to everything else, by construction:

| what actually differs | n | % |
|---|---:|---:|
| ligand body connectivity | 53 | 53.5% |
| bond order | 14 | 14.1% |
| geometry token, same CN (`OCT→TBP`, `SPL→SPY`, `TET→TPY`…) | 16 | 16.2% |
| stereo tags | 7 | 7.1% |
| fragment count | 6 | 6.1% |

Defects inside a ligand are outside a metal-contact probe's scope. The 16 same-CN geometry
changes are the documented `OGARAP` blindness — an η³→η² rearrangement preserves the count. This
is why **both instruments ship**: neither dominates, and the cheap one is the tripwire that
would catch a regression in the expensive one.

> **The acceptance bullet asked for disagreement within a ~3.7% band. What was measured is
> better than a band: every category of disagreement is named and attributed.**

## 4. The library did not move

A bucket report re-run over stored JSON proves only that the *classifier* did not move — it
never encodes anything. So the encoder was checked directly, at corpus scale: re-encode every
input XYZ under current `main` and diff against the `smiles_1` the v0.4.6 sweep recorded.

```
byte-identical : 4985/4985  (100.00%)
DRIFTED        : 0
encode errors  : 0
elapsed        : 582 s
```

`tools/encoder_identity_corpus.py`; 4985 is the correct denominator (15 molecules have
`smiles_1 = None`, the `encode_fail` bucket). **Zero drift** — so the correction reported here
is attributable to the scoring change and nothing else.

## 5. Why this was affordable, and why there is no A/B confound

The charter budgeted a **55 CPU-hour live re-sweep**. It was not needed. `save_artifacts`
writes `structures/<mol>_generated.xyz` from the *same* `gen_result.xyz` string that
`_attempt_generation` writes to the path `OIN_INDEP_SCORE` converts, so re-encoding the stored
file is **bit-identical to what the lever computes** — not an approximation of it. 4688 of 5000
molecules have such a structure (all 4660 successes plus the 28 failures that produced one); the
remaining 312 are failures in both arms with nothing to score.

Measured: **0.33 s/molecule median, 0.69 s mean, 334 s for the whole corpus.** The lever's own
rationale had priced this at "0.4–1.5 s/molecule, so a 5k sweep pays 1–2 CPU-hours" and held it
off on that basis. That figure was never measured; it is now, and it was never a reason to score
dishonestly. `levers.py` has been rewritten accordingly — a comment arguing against the code it
sits on is worse than no comment.

**The A/B confound the charter warns about does not apply to this measurement**, and it is worth
being precise about why rather than claiming immunity:

- Conformers are held fixed **by construction** — the same stored geometry is scored twice.
- `smiles_1` cannot move, because nothing is re-encoded from the input.
- No second encode enters `--mol-timeout 300`, so **no marginal pass converts into a timeout**.
- `metrics.elapsed_s` is therefore **unchanged**: `> 30 s` stays at 994 / 19.88%, median 7.19 s.
  The honest metric's own cost is recorded separately as `indep_encode_s` and is *not* inside
  the budget. **This release does not move runtime and does not claim to.**

The one thing offline cannot prove is that the *generator* still behaves — it never runs one.
That is what the bounded live arm covers (§7).

## 6. What the honest arm cannot score

4 molecules re-encode to nothing: `AFIROW_comp_0`, `MUZZUC_comp_0` (quinoid de-aromatization
leaves an invalid molecule), `HOBBUY_comp_0` (cannot kekulize, no quinoid ring to relax),
`RIPDEA_comp_0` (`N` valence 5). 0.08% of the corpus. They move `byte_exact → hard_fail` and are
counted as such rather than quietly dropped — an encoder that raises on a generated geometry is
a real limitation, not a scoring choice.

## 7. Scope: what changed, and what deliberately did not

`OIN_INDEP_SCORE` moved `_HELD_OFF → _DEFAULT_ON`. **That promotion changes what is *reported*,
not what is *accepted*.**

- `accept_fn` — untouched. The generator's pool acceptance is unchanged.
- The harness's **tier ladder** — untouched. It still escalates on `smiles_2`. Scoring the ladder
  honestly would change runtime *and* the failure mix in the very release that re-baselines the
  number, making both unmeasurable. That is the identical fault that let the `OIN_ACCEPT_SCORED`
  A/B report "zero regressions" while its honest arm read 8.
- **No molecule the honest metric newly fails was fixed here.** They are a worklist. A fix landed
  alongside a re-baseline makes the before/after unreadable, which is the entire cost this
  release is paying to avoid.

`src/oinsmiles/` behaviour is unchanged apart from the lever registry entry; §4 is the proof
rather than the assertion.

## 8. Artifacts

| what | where |
|---|---|
| honest re-scored twin (5000 reports) | `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/` |
| frozen side-by-side table | `…/bucket_report_both.{md,json}` |
| honest-only table | `…/bucket_report_honest.{md,json}` |
| per-molecule transition log + `#DONE` sentinel | `…/honest_rescore.jsonl` |
| encoder byte-identity log | `…/encoder_identity.jsonl` |
| the re-scorer | `tools/honest_rescore.py` |
| the corpus encoder gate | `tools/encoder_identity_corpus.py` |
| both-mode instrument | `tools/roundtrip_bucket_report.py --score {scored,honest,both}` |
| fixtures + tests | `tests/unit/test_honest_score.py`, `tests/fixtures/honest_score/` |

`FIYHUT_comp_0` scores a failure; `YOSXIP_comp_0` scores a pass; `OGARAP_comp_0` is caught by
the honest arm while `coordination` is blind to it. All three are pinned, driven from vendored
geometry so no generator run is required. **None of the three is in the 5000-molecule corpus** —
they live in `results-v0.4.5-rebaseline`, which is why the fixtures are vendored rather than read
from a sweep directory that a given machine may not have.
