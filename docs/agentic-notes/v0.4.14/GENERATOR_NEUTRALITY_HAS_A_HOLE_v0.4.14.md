# v0.4.14 — the generator-neutrality licence has a hole, and two releases spent it

**`tools/fold_key_invariance.py` reading `0 keys changed` does NOT license an offline re-score.**
It bounds the *acceptance* step and says nothing about the *embedding* step. v0.4.13 cancelled a
chartered 55 CPU-h sweep on that licence and v0.4.14 skipped one. Both were **lucky, not right**.

This note supersedes the headline in the first draft of
[`LANE-01-resonance-fold.md`](LANE-01-resonance-fold.md). Corrected numbers in §4.

---

## 1. The argument, and where it breaks

The licence, as stated in `fold_key_invariance.py` and quoted in three charters:

> `accept_fn` decides by comparison key, so a change that never moves a key cannot alter what the
> generator returns. Therefore an offline re-score of the frozen corpus is **exact**.

The first sentence is true. The second does not follow, and the gap is one line of the call chain:

```
gen.generate(oin_in)  ->  generation/oin_parser.OINParser  ->  ParsedOIN
                      ->  metallogen_adapter: m-SMILES, CoordMap pins donors to template vertices
```

**The generator's input is the OIN _string_, not the key.** A slot relabeling — exactly what a
canonicalization lever does — produces a different `ParsedOIN`, a different CoordMap, and therefore
a different embedding. Key-invariance guarantees `accept_fn` cannot pick a *different conformer out
of a given pool*. It does not guarantee **the pool is the same pool**.

The offline re-score holds the generated structure fixed by construction, so it is structurally
incapable of observing this. Worse, for a molecule that already round-trips it can only ever print
*"still fine"* — **an offline re-score reports `bad_direction = 0` whether or not losses exist.**

> Ask what a BROKEN version of this instrument would print. For losses, the offline re-score prints
> the same thing as a working one.

## 2. Measured, not argued

`tools/generator_ab_honest.py` (new) runs the full shipped path in both arms and scores **honestly**
— `oin_out` is re-perceived from the written XYZ, never taken from `res.mol`:

```bash
$V tools/generator_ab_honest.py --cohort-dir <cohort-v0.4.5-5k> \
    --molecules-file <sample> --lever OIN_RESONANCE_DONOR_FOLD --out-json ab.json
```

| sample | n | result |
|---|---:|---|
| drawn from the **78 claimed gains** (seed 7) | 20 | **19 real gains, 0 losses** |
| drawn from the **39 movers already `byte_exact`** (seed 11) | 20 | **0 gains, 3 REAL LOSSES** |
| `HEKFEL_comp_0`, `--repeat 2` | 1 | **REAL LOSS**, `OFF=True → ON=False`, 0 non-determinism |

**`GENERATED output moved` — 10 of 20 and 8 of 20 respectively.** That counter is the direct
observation of the hole: the lever changed the string the generator consumes, so it embedded
differently. The offline sim cannot produce that number at all.

### The population the offline sim was blind to

The 179 string-movers split three ways, and only the first was ever measured end-to-end:

| | n | what the offline sim could see |
|---|---:|---|
| claimed gains | 78 | measured — and ~95% confirmed real |
| **already `byte_exact` at v0.4.13** | **39** | **nothing. These can only LOSE** |
| other (no bucket change either way) | 62 | no change either way |

## 3. This is not the stochastic-A/B trap, and that was checked

The standing rule is *never A/B by re-running a stochastic harness*. It does not bind here:
`MetalloGenAdapter` takes `seed=42`, and `HEKFEL_comp_0` reproduces its hashes exactly across
repeated runs in both arms (`--repeat 2`, `nondeterministic: []`). **A seeded harness can be
A/B'd.** What is *not* deterministic is the pool SIZE — `OIN3DGenerator(timeout=)` is advisory, so a
loaded box can produce fewer attempts. Hence `--repeat`, on every cohort, rather than inheriting
this one measurement.

## 4. Corrected headline

| | claimed (offline) | corrected (honest end-to-end) |
|---|---|---|
| gains | 78 | **~74** (19/20 sampled confirmed = 95%) |
| losses | **0** | **~6** (3/20 sampled of 39 at-risk = 15%) |
| net | +78 mol / **+1.56 pts** | **~+68 mol / ~+1.36 pts** |
| `byte_exact` | 77.44% | **~77.25%** |

Range over the plausible corners of both samples: **+0.94 to +1.51 points.** The lever stays
**default-ON** — it is net-positive across the entire range — but the release ships a **point
estimate with a stated method**, not a measured figure. **The exact number requires a sweep**, and
that sweep is now owed.

## 5. What the losses actually are — and why they are a v0.4.15 finding

The 4 confirmed losses are **not** the encoder emitting a wrong string. The encoder's output is
canonical and correct in both arms. They are molecules where the generator, handed an equally valid
but differently-*labelled* input string, builds a worse structure.

> **The generator's output depends on the slot labeling of its input. It should not.**

Two OIN strings differing only in which interchangeable donor got slot 0 describe the same complex,
and a generator that reproduces one but not the other is label-sensitive. That is the same *class*
of defect as the 183 `MIRROR_MATCH` molecules — the generator being sensitive to something the
notation says is not a distinction — and it belongs on the generator ladder, not the encoder one.

**It also means every future canonicalization lever carries this cost**, because every one of them
relabels slots. That is a standing argument for fixing the generator's label-sensitivity *before*
spending more releases on canonicalization.

## 6. What has to change in how this project measures

1. **`fold_key_invariance.py` is necessary, not sufficient.** Its `GENERATOR_NEUTRAL` verdict
   licenses "acceptance is unchanged". It must stop being read as "an offline re-score is exact".
   The tool's own output text is corrected to say so.
2. **Any lever that relabels slots owes a generator-side A/B**, on both the claimed gains *and* the
   already-passing movers. The second sample is the one nobody ran.
3. **`arm2` is a byte-identity gate, not an accuracy instrument.** `gate_arm2_roundtrip_one.py`
   scores with `get_oin_string(result.mol, coords)` — the generator's own bond graph, the circular
   predicate `OIN_INDEP_SCORE` replaced at a measured 9.6% false-positive rate. Reading its hashes
   as a round-trip verdict produced one wrong conclusion in this session before it was caught.
4. **v0.4.13's +3.42 rests on the same licence** and has the same unmeasured at-risk population. It
   is not re-derived here — that needs its own generator-side A/B — but it should not be quoted as
   exact until it is.
