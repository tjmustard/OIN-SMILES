# v0.4.14 — close-out

**`byte_exact` 75.88% → 77.30%, +1.42 points. 78 gains, 7 losses, net +71 molecules,
measured end-to-end over all 182 molecules the lever can affect.**

---

## 1. Predicted vs actual

| | predicted | actual | delta |
|---|---|---|---|
| `byte_exact` | **UP ~3.5–4.1 pts** | **+1.42** | **−2.1 to −2.7** |
| `> 30 s` count | FLAT | **+7 on the affected set** (29 → 36 of 182) | worse than predicted |
| `max(elapsed_s)` | FLAT | 323.1 s → 705.5 s **on the affected set** | worse than predicted |

### Why the accuracy miss — the charter's arithmetic was right and its attribution was wrong

The charter budgeted ~4.08 points from two blocks: `rdkit_canonical` (114 / 2.28) plus a
~90-molecule frozen-resonance class (1.80). Measured:

- **`rdkit_canonical` contributes 0.** It is **80.7% η-set denticity drift** — the generated ring
  slips and fewer carbons fall inside the metal-bonding cutoff. No string canonicalization can
  reach it. It was never a canonicality block; the bucket name was a hypothesis.
- **The resonance class is 103, not 90**, and **78 of it converted** — a 76% close rate on the part
  that was genuinely reachable.

So the lane took 78 of a reachable 103. The shortfall against the prediction is **114 molecules
that were never reachable by the kind of change the charter proposed**, plus the 7 losses.

### Why the runtime miss

The prediction assumed nothing default-path touched generation. **That assumption is this
release's central finding, and it is false**: the lever changes the emitted string, the generator
consumes the string, so it embeds differently — `GENERATED output moved` on **85 of 182**.

## 2. Runtime — measured on the affected set, NOT carried forward

Previous releases carried `994/5000 = 19.88% > 30 s` forward on the grounds that nothing
default-path changed in generation. **That justification does not hold here** and the number is not
re-asserted.

| affected population, n = 182, end-to-end (encode + generate + re-encode) | OFF | ON |
|---|---:|---:|
| median | 2.58 s | 3.71 s |
| mean | 22.89 s | 36.03 s |
| max | 323.1 s | **705.5 s** |
| total | 1.16 CPU-h | **1.82 CPU-h (+57.4%)** |
| `> 30 s` | 29 | **36 (+7)** |

Per-molecule delta: median **+0.70 s**, worst **+511.5 s** (`RAXJEH_comp_0`), best **−295.5 s**
(`RANBIT_comp_0`). The cost is bimodal by molecule, which is this project's standing finding about
every speed number it has ever taken.

⚠ **The corpus `994/5000` figure is NOT updated to `1001`, and the arithmetic that tempts you to do
that is wrong.** That figure comes from the harness's `metrics.elapsed_s`, which is a **SUM over up
to three separately SIGKILLed attempts**; the table above is a single in-process `generate()`. They
are not the same quantity and must not be added. The honest statement is: **the corpus figure is
carried forward unverified, and 7 molecules of the affected set newly cross 30 s.**

## 3. What this release actually delivered

**A smaller headline than chartered, and three findings worth more than the headline.**

1. **The generator-neutrality licence has a hole.** `fold_key_invariance.py` reading `0 keys
   changed` bounds the *acceptance* step — `accept_fn` decides by key — but the generator's input
   is the OIN **string**, so a slot relabeling changes `ParsedOIN`, the CoordMap, and the pool
   itself. An offline re-score holds the structure fixed, so **it reports `bad_direction = 0`
   whether or not losses exist**. v0.4.13 cancelled a chartered 55 CPU-h sweep on this licence.
2. **`key_equal` is not "benign canonicalization".** 183 of its 361 members (50.7%) are the
   generator building the **enantiomer**, invisible because `compare._parse_vertex_colors` folds
   reflection and `accept_fn` decides by that key.
3. **The generator's output depends on the slot labeling of its input.** That is what the 7 losses
   are — the encoder is canonical in both arms. Every future canonicalization lever pays this toll.

**And a correction to a published number:** v0.4.13's **+3.42 is ~+2.82** — 197 at-risk molecules,
never sampled, measured loss rate 6/40 = 15%.

## 4. No sweep was needed, for a better reason than v0.4.13's

v0.4.13 skipped its sweep on the neutrality licence, which was the wrong argument. The right one:
**the affected population is derivable.** A molecule whose `encode(input)` is byte-identical in both
arms hands the generator the same string, and generation is seeded, so it is unchanged **by
construction**. `tools/lever_string_movers.py` derives it from coordinates — 93 input-movers, 182
including the generated side — and an A/B over that set is *exact*.

⚠ **Derive it from coordinates, never from a frozen sweep's stored strings.** The stored-string
derivation read **179** and was wrong in both directions: 89 of those do not move `encode(input)`
under today's encoder, and **3 that do were missing**. Every coverage figure quoted against that
179 — including the mirror audit's "179/179 = 100%" — inherits the error.

## 5. Gates

| | result | coverage |
|---|---|---|
| ARM 1 | PASS byte-identical, `#DONE 62` | **0 of 62** — cannot see this change |
| arm2 goldens | 7/325 and 1/100 rows re-frozen, fields 1–6 only | 12/325, 2/100 |
| mirror audit, mover-enriched | **0 regressions**, 0 per-molecule verdict changes | 179/179 of the *stored-string* set |
| mirror audit, uniform `cat` | 0 regressions; reproduces v0.4.13's 157/92/1 | **1/250 = 0.4%** |
| `tests/unit` | see §6 | — |

## 6. Release hygiene

- `tests/unit` green · ruff 0.15.20 clean
- `pyproject` 0.4.13 → 0.4.14 · CHANGELOG stanza written · tagged `v0.4.14`
- `measurements/v0.4.14/` frozen (allowlist extended for 10 new instrument filenames)
- **Nothing pushed.**

## 7. Handed forward

- **v0.4.13's gain side is unsampled**, so its `~+2.82` is an upper bound.
- **The corpus `> 30 s` figure is unverified** since v0.4.8.
- **`OIN_PREFILTER_ADVISORY`** — prevalence measured at 1.5% of vetoes / 4.5% of molecules; stays
  off, but its latency objection is **refuted** (it is measured *faster*). Do not re-open it on cost.
- **The 25-molecule resonance residue and the 39 `NOT_A_MIRROR`** — the encoder ladder's remaining
  1.28 points.
