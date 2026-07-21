# Round-trip failure elimination study (v0.4.3 research thrust)

**Status:** Complete. Phase 1 (zero-generation experiments) and Phase 2 (silent-degradation
telemetry, 861 molecules) both done. **Every table below is FF floor: `xtb` is not reachable in this
environment, so `optimizer_effective: "ff"` on 100% of rows.**

This study does not fix anything. Its purpose is **elimination**: to establish with evidence
what is *not* causing round-trip failures, so that a later development thrust spends its
effort where the cause actually is. Verdicts are strictly **REFUTED / SURVIVES / UNDETERMINED**.

The motivating question was posed as: is the cause (1) MetalloGen generation, (2) the OIN
logic and its conversion to MetalloGen format, (A) one of these, (B) both, (C) neither, or
(D) a combination plus something else?

**Answer: (D), but not in the proportions the question assumes.** See §5.

---

## 1. Executive summary — the denominator is wrong

### 1.1 "11.6% failure" is five different problems

The 2,917 failures in the 25,197-molecule `--quick` sweep partition cleanly by **how far down
the pipeline each one got**, which is recoverable from the stored strings alone:

| Bucket | What broke | n | % of failures |
|---|---|---:|---:|
| **A — encoder** | XYZ→OIN never emitted a string | **417** | 14.3% |
| **B — no conformers** | pool exhausted | 404 | 13.8% |
| **C — timeout** | hit the `--quick` 30 s wall | **1,144** | 39.2% |
| **D — other exception** | | 116 | 4.0% |
| **E — gate** | both strings exist; string/RMSD gate disagreed | 836 | 28.7% |

Bucket A splits into **240 encode timeouts + 177 `get_lig_mol` perception failures**.

**Bucket A never reached MetalloGen and never reached the OIN→m-SMILES handoff.** Both
proposed hypotheses are therefore refuted outright for 14.3% of the failure mass, before any
experiment is run.

### 1.2 A plurality of the rest is a budget artifact, not a defect

Re-judging the same molecules against the non-quick capstone sweep (1800 s budget) over the
6,719-molecule overlap:

- **352 of 652 adjudicated quick failures (54.0%) pass when given a real budget.**
- `C_timeout` collapses from **39.2% → 0.96%** of failures.
- Only 12 quick passes fail non-quick.

In the non-quick sweep the shape is entirely different: `B_no_conformers` becomes dominant
(40.7%) and timeouts nearly vanish.

**Consequence: the addressable "the generator makes bad structures" population is
roughly B+D+E ≈ 1,356 molecules ≈ 5.4% of corpus, not 11.6%.** Any analysis that uses the
headline failure rate as its dependent variable is measuring the timeout budget as much as
the generator.

### 1.3 The leading suspicion is bounded before it is tested

The prior suspicion was that failures originate in initial ligand placement around the metal.
That mechanism is **structurally incapable** of explaining buckets A + C (53.5% of failures),
and cannot produce a string-identity mismatch. It remains live for the distortion question —
where §4 shows it is in fact strongly supported.

### 1.4 The largest finding was invisible to the existing gate

> **ρ(coordination-sphere RMSD, full heavy-atom divergence) = 0.220** (n = 6,361 paired).
> ρ(coord_rmsd, clash count) = 0.165. ρ(coord_rmsd, MPO) = −0.345.

The pre-registered decision rule was: ρ > 0.8 means the cheap gate is a valid proxy;
**ρ < 0.5 means the pass/fail gate is measuring the wrong thing.** At 0.22 the gate is very
nearly uncorrelated with geometric quality.

| Cohort | n | % with ≥1 vdW clash | median MPO |
|---|---:|---:|---:|
| Passed both gates | 6,295 | **52.9%** | 0.700 |
| Failed (with geometry) | 109 | 63.3% | 0.589 |
| **Real crystal inputs (control)** | 6,404 | **5.3%** | 0.912 |

Passing structures clash at **ten times the rate of the crystal structures they were built
from**, and nearly as often as the failures. Distortion is population-wide, and the
round-trip gate cannot see it.

---

### 1.5 The "FF floor" is a venv split, not a missing dependency

Every row of every dataset carries `xtb_available: false` / `optimizer_effective: "ff"`, and
this study has captioned all its tables accordingly. The cause turns out to be configuration,
not absence:

| | |
|---|---|
| `uv run` interpreter | `/home/tjmustard/.venv/bin/python3` (**home** venv) |
| `xtb` binary | `…/OIN-SMILES/.venv/bin/xtb` (**project** venv), **xtb 6.7.1, works**, installed 2026-05-14 |
| `shutil.which("xtb")` under `uv run` | `None` |

Both `tools/test_dataset_roundtrip.py:84` and `generator3d/ml_optimizer.py:54` gate on
`shutil.which("xtb")`. The binary has been present since May but sits in a venv the harness
does not run under, so every sweep silently took the force-field path — exactly the
"deliberate, non-fatal degradation" that `ml_optimizer` documents.

**Consequence for the fix thrust:** the most-recommended next step (get a real optimizer into
the loop) is largely already done and needs a PATH change, not an install. But see Appendix A —
**the moment `xtb` becomes reachable, the multiplicity-forcing defect activates for 36.8% of
the corpus**, because that is when `.UHF` starts being consumed. These must land together.

---

## 2. Method

**Attribution set.** All causal claims use `results-capstone-v042`: non-quick, single clean
commit `58bba7ad`, `mol_timeout=1800`, 6,719 rows. `results-v0.4.0` is used only for
population framing and as a bulk geometry corpus.

That rule turned out to be more necessary than expected. Measured provenance:

| set | commit_ids | quick | statuses |
|---|---|---|---|
| `results-v0.4.0` | **five**: `c7edeeb6` (12,141), `5538b722-dirty` (10,574), `c7edeeb6-dirty` (2,394), `d832ab74` (49), `b12ab69e-dirty` (26) | true | 22,280 ok / 2,917 fail |
| `results-capstone-v042` | **one**: `58bba7ad` | false | 6,295 ok / 312 fail / 112 pending |

Three of the five v0.4.0 commits are dirty working trees. That dataset is not exactly
reproducible from any commit and cannot carry causal claims.

**Exclusions.** The 112 `pending_g-xtb` capstone rows are neither pass nor fail and are
excluded from every rate; the exclusion is reported rather than silently absorbed.

**Reproduction.** Every number below comes from a committed tool on branch
`research/v043-elimination`:

| Tool | Produces |
|---|---|
| `tools/elimination_corpus_stats.py` | §1.1 partition, §1.2 adjudication, flake floor, provenance |
| `tools/synthetic_oin_battery.py` | §3 handoff detection in isolation |
| `tools/elimination_static_sweeps.py` | §3 handoff sizing on real chemistry |
| `tools/structure_distortion_report.py` | §1.4 / §4 geometric quality (+ `docs/DISTORTION_v0.4.3_RESEARCH.md`) |

---

## 3. The handoff hypothesis was tested two ways and did not survive

Hypothesis (2) — that the OIN logic mis-specifies the SMILES and its conversion to MetalloGen
format — is the cheapest to test decisively, because `convert_parsed_to_msmiles` is a **pure
function**: no 3D input, no stochasticity. A mismatch is therefore a *proof* of corruption,
not a correlation.

### 3.1 Detection: 104 synthetic complexes of known composition

Hand-built OIN strings over every supported geometry × simple monodentate ligands
(Cl⁻, NH₃, H₂O, CO, CN⁻, OH⁻, CH₃⁻), plus adversarial probes aimed at each named hazard.
Composition, donor set and coordination number are known **by construction**, so no other
part of the pipeline can contribute error.

- **Slot-map bijection: 0 collisions across all 11 geometries.** Each OIN template slot maps
  to a unique MetalloGen slot, so no ligand can be silently overwritten.
- **Composition fidelity: 101 / 104 clean.**
- The 3 violations are all the same case, and it is **not an adapter bug** — see §3.3.

### 3.2 Sizing: 31,339 real OIN strings, no 3D generated

Running the real parser and the real handoff over every stored `.oin` in **both** corpora:

| Check | non-quick (capstone) | quick (v0.4.0) | Verdict |
|---|---|---|---|
| Slot drop (parser discards unresolvable slot) | **0 / 6,559** | **0 / 24,780** | REFUTED |
| Ligand slot collision | **0 / 6,559** | **0 / 24,780** | REFUTED |
| Heavy-atom formula preserved | **0 mismatches / 6,534** | **0 mismatches / 24,664** | REFUTED |
| Ligand count preserved | **0 mismatches / 6,534** | **0 mismatches / 24,664** | REFUTED |
| Handoff refuses outright | 25 (0.38%) | 116 (0.47%) | SURVIVES, small |
| Actinide silently rewritten | **0** | **0** | REFUTED |

**Across 31,223 molecules that reached m-SMILES, the handoff never once lost a heavy atom or a
ligand.** Heavy atoms are the deliberate comparison basis: unlike hydrogen they carry no
phantom-implicit convention, so a mismatch would be a real loss rather than a difference in
reading.

**Stereochemistry is preserved too.** Of 6,534 non-quick molecules, 831 (12.7%) carry a stereo
annotation and 221 show a textual `@`/`@@` count delta — but **0 of those 221 lost a perceived
stereocentre**. Every delta is RDKit re-expressing the same stereochemistry under a different
atom ordering.

The `@SP1/2/3` stripping at `metallogen_adapter.py:127-128` applies only to the *metal*
fragment, and **no metal fragment in the corpus carries an `@SP`/`@TB`/`@OH` tag at all** — the
regex never fires. (The 35 v0.4.0 molecules that match the pattern carry it on a ligand
sulfur, `[S@SP3]`, which that code path never touches.)

One number needs care: the textual-mismatch set fails at 10.0% against a 4.0% baseline. Since
**none** of those mismatches is a real stereo loss, this is confounded rather than causal —
molecules with more stereocentres are larger and more intricate, and so fail more often for
unrelated reasons. It is not evidence of a stereo defect.

The refusals (24 + 110 `UncoordinatedFragmentError`, 1 + 6 `ValueError`) all fail the round
trip, but they are **8% of the 312 non-quick failures / 4% of the 2,917 quick failures**, and
`UncoordinatedFragmentError` is a representation limit — outer-sphere counterions and free
solvent are unrepresentable in `metal|lig|geo` m-SMILES — not a defect in the conversion.

### 3.3 The one real information loss is in the format, not the code

The 3 synthetic violations are the bare-nitrogen case. `convert_parsed_to_msmiles` reads a
bare `N{n}` donor with no heavy neighbour as ammine, emitting `[NH3:1]` — the behaviour its
own comment documents ("the ammine reading wins").

The battery establishes *why* this is unfixable in the adapter: the nitride and ammine cases
**produce byte-identical OIN strings**. Two chemically different complexes serialize to one
string. The adapter has no information to discriminate on. This is a **non-injectivity in the
OIN grammar**, and fixing it requires a format change in `oin/inline.py`, not a fix in the
handoff.

### 3.4 Comparator bugs found and fixed mid-study

Both were caught because a result looked like a finding, and the finding turned out to be the
instrument. They are recorded because this is the failure mode the study is most exposed to:
**a comparator that encodes intuition rather than the representation under test will
manufacture findings.**

**(a) Charge asserted from chemistry instead of from the representation.** The first battery
run reported 59 charge violations. In the OIN/m-SMILES convention an X-type donor is written
neutral and the complex charge lives on the metal (`[Cl]`, not `[Cl-]`); asserting "chloride
is −1" tests the convention, not the implementation. Expected charge is now read back off the
OIN fragment itself, where formal charge is explicit and unambiguous.

**(b) Stereo counted textually instead of perceived.** The stereo screen initially reported
~3% of molecules changing their `@`/`@@` count through the handoff. Adjudicating the flagged
molecules with `FindMolChiralCenters` showed **identical perceived stereocentre counts on both
sides** — RDKit's canonical output simply re-expresses the same stereochemistry when atom
ordering changes. The textual delta was canonicalisation, not loss. The screen is retained
(it is cheap and it filters), but perception now decides.

The same reasoning forced a correction to the planned corpus-wide audit before it was written:
the "OIN-implied graph" cannot be a naive RDKit re-parse, because OIN's bare-donor convention
makes implicit H phantom, and that convention *is* the code under test. Heavy atoms — which
carry no such convention — became the comparison basis, with the crystal input as the only
independent oracle for hydrogen.

---

## 4. What the geometry actually shows

From `docs/DISTORTION_v0.4.3_RESEARCH.md` (6,404 non-quick structures, crystal inputs as
control cohort), plus the correlation analysis added by this study:

- **53% of generated structures carry ≥1 van-der-Waals clash vs 5% of the real inputs**
  (mean 2.4, worst 28). Bond-angle strain 6.3° vs 3.4°.
- **Bond *lengths* are not the problem** — the generator places near-ideal bond lengths and
  scores *better* than real crystals (0.038 vs 0.087 deviation). This metric is inverted and
  is excluded from the quality score.
- **MetalloGen's own acceptance gate is atomic-fusion-only** (covalent-radius ratio < 0.6).
  100% of shipped structures pass it while half carry vdW clashes.
- Distortion concentrates in **high-CN, low-symmetry geometries** (PBP 89% clashing, SQA 71%,
  TBP 68%, SPY 67%; LIN 27%) and **large early oxophilic metals** (Y 92%, Sc 97%, Ta, Nb, Hf,
  Zr, Ti), and is cleanest for Au, Hg, Ni, Zn.

That gradient — crowding worst where ligands are bulkiest and the coordination sphere most
crowded — is consistent with donors being placed onto ideal template axes without any
inter-ligand steric term. **The original suspicion about ligand placement is supported for the
distortion question, even though it was refuted as an explanation of the failure rate.**

---

## 5. Answer to (A) / (B) / (C) / (D)

**(D): a combination, plus something else — and the "something else" is the largest term.**

| Contributor | Attributable mass | Confidence |
|---|---|---|
| **Measurement**: the round-trip gate is blind to geometry (ρ=0.22); distortion is population-wide, affecting passes and failures alike | ~53% of *all* generated structures, passing ones included | High |
| **Generator (hypothesis 1)**: embed acceptance criteria are essentially never satisfied, so the pipeline ships the least-bad attempt; steric crowding from template-axis placement with no inter-ligand term; permissive fusion-only clash gate | dominant contributor to distortion; B+D+E ≈ 1,356 molecules for outright failure | High — mechanism now measured (§8b) |
| **Encoder** (neither hypothesis): `get_lig_mol` perception + encode timeouts | 417 (14.3% of failures) | High |
| **Budget** (neither hypothesis): `--quick` 30 s wall | up to 1,144 (39.2%), 54% of which are recoverable | High |
| **Handoff (hypothesis 2)** | ≤ 25 molecules (0.38%), all representation limits | High — REFUTED as a material cause |
| **Format non-injectivity** (nitride/ammine) | small, but unfixable in the adapter | High |

Hypothesis 2 as posed — that the OIN logic mis-specifies the conversion — is **refuted**.
Hypothesis 1 is **supported for structure quality but not for the failure rate**, and the
single largest effect is one that neither hypothesis named: **the instrument cannot see the
defect it is supposed to be measuring.**

---

## 6. Elimination table

| ID | Layer | Claim | Verdict | Evidence | Conf. | Scope |
|---|---|---|---|---|---|---|
| H1 | Encoder | XYZ→OIN perception/timeout failure precedes any generation | **SURVIVES** | 417/2,917 = 14.3% quick; 48/312 = 15.4% non-quick | High | A |
| H2b | Parser | Unresolvable slots silently dropped | **REFUTED** | 0 dropped slots / 31,339 | High | all |
| H3 | Handoff | OIN→m-SMILES loses/corrupts chemical information | **REFUTED** | 0 formula mismatches / 31,223; 0 ligand loss; 101/104 synthetic | High | all |
| H3a | Handoff | `@SP1/2/3` stripping loses square-planar stereo | **REFUTED** | 0 real stereocentre losses / 6,534; no metal fragment carries the tag at all | High | all |
| H3b | Handoff | Nearest-vector slot collision overwrites a ligand | **REFUTED** | 0 collisions across 11 geometries and 31,339 molecules | High | all |
| H3c | Handoff | H-count heuristics mis-protonate donors | **SURVIVES (narrow)** | only bare-N; a format limit, not an adapter defect (§3.3) | High | bare-N donors |
| H3d | Handoff | Kekulization failure + Cp charge guess corrupts aromatics | **REFUTED** (for composition) | 0 formula/charge corruption on real corpus | Medium | aromatic ligands |
| H3x | Handoff | Handoff refuses representable complexes | **SURVIVES (small)** | 0.38% non-quick / 0.47% quick; 90%+ `UncoordinatedFragmentError` | High | B/D |
| H4 | Generator | Ligand placement produces distorted structures | **SURVIVES (strong)** | 53% clashing vs 5% crystal control; strata track crowding | High | population-wide |
| H4a | Generator | "Best rejected candidate" ships known-bad geometry | **SURVIVES (reframed)** | fires on 86.3% of *clean* passes — it is the normal path, not an escape hatch; but enriched in distorted passes (OR 3.41, CI [1.80, 6.44]) and intensity rises 5.7→13.8→26.3→30.7 across strata | High | population-wide |
| H4b | Generator | PuLP+xyz2mol double failure flattens the graph | **REFUTED** | 1.7% of clean passes, 0.0% of bucket E; only elevation OR 3.16 CI [1.06, 9.41], below threshold | High | all |
| H4c | Generator | Dummy-center substitution distorts geometry | **UNDETERMINED** | not separately probed; no site isolates it | — | — |
| H4g | Generator | Blanket `except Exception` hides under-coordination failures | **SURVIVES (narrow)** | 4.6% in buckets B+D, 0.0% in all other strata, OR 30.98 (wide CI [1.76, 546]) | Medium | B/D |
| H4d | Generator | Disabled clash guards let clashes through | **SURVIVES (partial)** | 100% pass the fusion-only gate while 53% carry vdW clashes | Medium | population-wide |
| H4e | Generator | Multiplicity forced to singlet/doublet | **REFUTED as cause** / real latent defect | 35.4% of complexes are on high-spin-capable metals, all forced; but the FF path never consumes `.UHF` | High | latent |
| H4f | Generator | Actinide silently emitted as lanthanide | **REFUTED** | 0 actinides in corpus | High | — |
| H5b | Optimizer | Energy-parse failure poisons ranking | **UNDETERMINED** | no `xtb` on any row; untestable here | — | — |
| H6 | Measurement | Gate order hides geometry data | **CONFIRMED** | RMSD computed only after the string gate passes | High | E |
| H6a | Measurement | Coordination-sphere RMSD is blind to the defect | **CONFIRMED** | ρ = 0.220, far below the 0.5 invalidation threshold | High | all |
| H6b | Measurement | `--quick` is materially a different generator | **CONFIRMED** | 54% of quick failures pass non-quick | High | all |
| H5 | Selection | Geometry/winding sampled-and-filtered, not constrained | **SURVIVES** | QOSMER returns a *different SPY isomer* between identical runs (§7) | Medium | stochastic minority |
| H5a | Selection | Selection *fallthrough* silently ships the wrong conformer | **REFUTED** | geometry-select 0–0.9%, winding 0–0.9%, stereo-fallback 0–0.9% across all strata | High | all |
| H6c | Measurement | Single-run gates are noisy | **SURVIVES (narrow)** | ~98% of molecules reproduce byte-identically, but a minority does not; an earlier "0.00%" claim in this report was wrong (§7) | High | per-molecule claims |

---

## 7. Instrument validity

**Determinism is high but NOT universal — this corrects an earlier draft of this report.**

A 32-molecule stratified sample (16 passes, 8 bucket-E, 8 bucket-B/D) run twice under
identical invocation gave a perfect result: 0 gate-outcome flips, `smiles_2` identical for all
32, RMSD delta 0.0000 on all 19 measurable, **24/24 generated XYZ byte-identical**. That was
initially written up as "flake floor = 0.00%".

**That claim was too strong, and a later experiment refuted it.** 0 flips out of 32 has a 95%
upper confidence bound near 9% — it never supported "0%". The telemetry inertness run (§8),
on a *different* 20-molecule sample, found 1 structure differing between arms. A dedicated
control on that molecule settled it:

| run | `OIN_TELEMETRY` | generated XYZ sha256 |
|---|---|---|
| 1 | unset | `b05e7e92…` |
| 1 | **=1** | `b05e7e92…` — identical |
| 2 | unset | `b28144b7…` — **differs from run 1 under identical conditions** |

`QOSMER_comp_0` is genuinely nondeterministic run-to-run. The two outputs are not noise: they
are **different geometric isomers** — under SPY (square pyramidal) the apical slot 0 is taken
by the allenylidene carbene in one run and by the NHC in the other.

**Corrected statement:** the generator is deterministic for the large majority of molecules
(56 of 57 distinct molecules observed across both experiments, 98.2%), but a minority is
genuinely stochastic, and for those the *isomer* — not merely the coordinates — can change
between runs. This is direct supporting evidence for **H5**: coordination geometry is
sample-and-filter over a stochastic pool, not a constraint, so when the filter does not
discriminate the pool's randomness reaches the output.

**Consequences.** Effect sizes do not need wholesale discounting, but any single-molecule
claim must be replicated, and a per-molecule A/B is only meaningful with repeats. Byte-identity
remains a valid acceptance test for *most* molecules, which is why the inertness proof (§8)
needed a control rather than a pass/fail verdict.

**On the prior belief.** Project memory records TIPYEX flipping and the metric being
"conformer-flaky at the margin". TIPYEX is in the 32-molecule sample and did **not** flip; its
historical flip was *sharded vs. alone*, a different work partitioning. So there are plausibly
two distinct effects — partition sensitivity and per-molecule stochasticity — and this study
has now demonstrated the second directly. Shard-composition invariance remains untested.

**Metrics disqualified.** `bondlen_dev` is inverted (generator beats crystal) and the
ligands-only UFF relax proxy is confounded by its all-single-bond metal-free model; both are
informational and excluded from the quality score. A valid relaxation metric needs a
metal-capable engine that is not installed here.

---

## 8. Phase 2 instrument

The instrument for the four generator-internal / selection hypotheses. Its results are in §8b;
this section documents the instrument and its validation.

**Instrument.** `src/oinsmiles/generation/_telemetry.py` plus probes at six silent-degradation
sites: `embed.pulp_and_xyz2mol_both_failed`, `embed.best_rejected_returned`,
`pool.blanket_exception`, `pool.stereo_fallback_wrong_isomer`, `adapter.winding_fallthrough`,
`adapter.geometry_select_fallthrough`. It is off unless `OIN_TELEMETRY=1`, cannot raise,
consumes no randomness, and keeps state in a `ContextVar` so the harness's spawn workers do
not cross-talk. The diff is **insertion-only: 17 additions, 0 deletions**, so no existing
statement is modified or reordered.

**Validation status: PASSED, after a control.**

1. Imports confirmed non-circular; enable/disable semantics confirmed.
2. **Full unit suite green on the instrumented tree: 435 tests OK, 3 skipped** — exactly the
   known baseline count, so the probes changed no tested behaviour.
3. **Byte-identity, `OIN_TELEMETRY` unset vs `0`: 20/20 identical**, outcomes and re-encoded
   strings identical.
4. **Byte-identity, unset vs `1`: 19/20.** One molecule (`QOSMER_comp_0`) differed.

The one difference is *not* instrument perturbation, established two ways:

- **Mechanically.** Nothing outside `tools/telemetry_run.py` ever opens a collection, so under
  the harness `record()` reaches `_events.get() is None` and returns before doing anything.
  With no active collection the function is a pure no-op whatever `OIN_TELEMETRY` says, and no
  other code reads that variable.
- **Empirically.** A dedicated repeat control on that molecule (§7) showed telemetry-on
  reproducing the baseline byte-for-byte, while two telemetry-*off* runs disagreed with each
  other. The molecule is nondeterministic; the instrument is not implicated.

The probes are therefore inert. The caveat that matters for Phase 2 is not instrument bias but
**per-molecule stochasticity**: firing rates aggregated over hundreds of molecules are sound,
single-molecule attributions are not.

Note the proof's own report flags that it could not confirm probe *wiring*, because the
round-trip harness does not persist telemetry events. `tools/telemetry_run.py` closes that gap
by opening a collection itself and refusing to run if `OIN_TELEMETRY` did not take effect.

**Sample, and its power limit.** The plan called for 300 molecules per stratum. Two strata
cannot supply that — the non-quick failure population is simply too small:

| Stratum | planned | available | 95% CI half-width |
|---|---:|---:|---:|
| S1 passes, cleanest MPO quartile (0.799–1.000) | 300 | 300 | ±5.7% |
| S2 passes, most-distorted quartile (0.085–0.598) | 300 | 300 | ±5.7% |
| S3 bucket E (string/RMSD gate failures) | 300 | **109** | ±9.4% |
| S4 buckets B+D (no-conformers / other exception) | 300 | **152** | ±8.0% |
| **total** | 1,200 | **861** | |

At n=109 the S3 arm resolves a 5% vs 25% firing-rate difference but **cannot** resolve 5% vs
10%. Any site whose estimate lands close on these strata is UNDETERMINED, not refuted, and
needs escalation against the quick corpus where the failure classes are ~9× larger.

---

## 8b. Phase 2 results — which silent fallback actually fires

861 molecules, non-quick, sequential on an idle machine, generator kwargs copied from the
harness PASS-1 configuration. Firing rate = fraction of molecules where the site fired at
least once.

| site | S1 clean pass | S2 distorted pass | S3 bucket E | S4 buckets B+D |
|---|---:|---:|---:|---:|
| `embed.best_rejected_returned` | **86.3%** | **95.7%** | **96.3%** | 28.9% |
| `embed.pulp_and_xyz2mol_both_failed` | 1.7% | 0.3% | 0.0% | 5.3% |
| `pool.blanket_exception` | 0.0% | 0.0% | 0.0% | 4.6% |
| `pool.stereo_fallback_wrong_isomer` | 0.0% | 0.3% | 0.9% | 0.0% |
| `adapter.geometry_select_fallthrough` | 0.0% | 0.0% | 0.9% | 0.0% |
| `adapter.winding_fallthrough` | 0.0% | 0.0% | 0.9% | 0.0% |

Mean firings per molecule — intensity, not merely presence:

| site | S1 | S2 | S3 | S4 |
|---|---:|---:|---:|---:|
| `embed.best_rejected_returned` | **5.7** | **13.8** | **26.3** | **30.7** |
| `pool.blanket_exception` | 0.0 | 0.0 | 0.0 | 11.5 |

### 8b.1 Four of six probed fallbacks are REFUTED outright

`pool.stereo_fallback_wrong_isomer`, `adapter.geometry_select_fallthrough`,
`adapter.winding_fallthrough` and `embed.pulp_and_xyz2mol_both_failed` all fire at or below
the pre-registered 1% refutation threshold in the passing strata, with odds-ratio confidence
intervals spanning 1 everywhere.

This matters because three of them looked alarming on code inspection. The
all-bonds-order-1/all-charges-zero collapse at `embed.py:292-295` — flagged as the single most
dangerous silent degradation in the pipeline — fires on **1.7% of clean passes and 0.0% of
bucket-E failures**, and its only elevation (5.3% in S4, OR 3.16, CI [1.06, 9.41]) has a
lower bound below the 1.5 decision threshold. **It is real, rare, and not a driver of these
failures.** The same holds for the wrong-stereoisomer fallback and both adapter fallthroughs.

### 8b.2 `embed.best_rejected_returned` is the finding, but not in the obvious way

It fires on **86.3% of the cleanest structures in the corpus.** A path documented as the
last-resort escape hatch — return a candidate that *failed* `_finalize_positions` validation,
possibly scoring −100000 for overlapping atoms — is the **normal operating mode of the
generator**, not an error path.

So the naive reading ("this fallback is shipping bad structures, remove it") is wrong: removing
it would fail 86% of generations. The correct reading is that **`_finalize_positions` almost
never accepts a candidate outright**, and the generator is best-effort by default.

But it is *also* genuinely discriminating, on two axes:

- **Enrichment.** S2 vs S1 gives OR **3.41**, CI [1.80, 6.44] — lower bound above the 1.5
  threshold, so distorted passes fire it significantly more than clean ones. (S3 vs S1,
  OR 3.75, CI [1.38, 10.18], falls just under the threshold and is suggestive only.)
- **Intensity — the strongest signal in the study.** Mean firings per molecule rise
  monotonically with severity: **5.7 → 13.8 → 26.3 → 30.7**. A clean dose-response across four
  independently-defined strata.

The S4 rate (28.9%) is a survivorship artifact, not a contradiction: 145 of 152 S4 molecules
**raised** before completing, so most never reached that code path — and those that did fired
it hardest (30.7 mean).

**Interpretation.** The number of validation-failing embed attempts is a direct correlate of
final structure quality. This is the mechanistic bridge between §4's population-wide distortion
and the generator internals: structures are not distorted because a rare fallback misfires,
they are distorted because the acceptance criteria are essentially never met and the pipeline
ships the least-bad attempt. The lever is the **acceptance criteria and the placement that
feeds them**, not the fallback.

### 8b.3 `pool.blanket_exception` cleanly separates the generation failures

0.0% across all three other strata, 4.6% in S4, OR **30.98** — enriched, though the CI is very
wide ([1.76, 546]) at these counts. Mean 11.5 firings per affected molecule. This is the
blanket `except Exception` at `generator3d/__init__.py:332` that swallows under-coordination
`TypeError`s. It is specific to the no-conformers/exception bucket and invisible everywhere
else, which is exactly the signature of a real but narrow defect.

Consistent with that: S4 outcomes are 145 raised / 7 generated, while S3 is 108 generated /
1 raised — confirming buckets B+D fail *in generation* and bucket E fails *at the gate*.

### 8b.4 Reproducing

Every table in §8b is regenerated from the committed raw run:

```
uv run python tools/telemetry_analyze.py
```

which reads `docs/data/v0.4.3/telemetry_events.json` (the 861-molecule run,
`tools/telemetry_run.py` output) and `docs/data/v0.4.3/telemetry_strata.json` (the sample
definition). The odds ratios use a Haldane–Anscombe (+0.5) correction so an OR is defined even
when a cell is zero, which shrinks estimates toward 1 — the conservative direction for an
elimination study.

---

## 9. What this study could not test

| Gap | Minimum experiment |
|---|---|
| Which silent fallback actually fires (H4a/H4b/H4c, H5) | Phase 2: insertion-only telemetry probes at the ~17 degradation sites, 1,200-molecule stratified run |
| Whether distortion survives a real optimizer (H5b, the FF-floor caveat on every table) | Install `xtb`/g-xTB; re-run a control cohort with `optimizer_effective != "ff"` |
| Geometry for ~2,088 of the 2,917 quick failures | No 3D was ever produced; requires re-running, and for timeout/no-conformers may not yield one |
| Non-quick behaviour outside the 6,719-molecule capstone overlap | Extend the non-quick sweep toward full corpus coverage |
| Invariance to shard composition/ordering (§7) | Same molecules run alone vs. sharded, compared byte-wise |

---

## 10. What remains standing

Ranked by attributable mass, as input to the fix thrust:

1. **The measurement instrument (H6/H6a).** Until the gate can see geometry, "88% pass" is not
   a statement about structure quality. A physical-quality gate (clash/strain, or full-molecule
   divergence) beside the string/RMSD gates is a precondition for trusting any future A/B.
2. **Embed acceptance criteria and the placement that feeds them (H4/H4a/H4d).** Phase 2 makes
   this concrete: `_finalize_positions` accepts outright so rarely that the best-rejected path
   runs on 86% of the *cleanest* structures, and the count of validation-failing attempts rises
   monotonically with distortion (5.7 → 30.7). Donors are placed on ideal template axes; the FF
   clean constrains ligand atoms with the metal fixed and has no inter-ligand non-bonded term;
   the acceptance gate forbids only atomic fusion. **Do not remove the fallback — fix what it
   is compensating for.**
3. **The encoder (H1).** 417 molecules, entirely upstream of both proposed hypotheses.
4. **Generator budget/robustness (B+D).** `no_conformers` dominates non-quick failures (40.7%).
5. **Format non-injectivity (H3c)** and **representation limits (H3x)** — small, real, and
   requiring format work rather than bug fixes.

Explicitly **not** worth the fix thrust's time, on this evidence: the OIN→m-SMILES conversion
logic, slot resolution, the nearest-vector slot match, actinide handling, run-to-run
determinism, **the PuLP/xyz2mol double-failure collapse, the wrong-stereoisomer fallback, and
both adapter selection fallthroughs** — the last four measured directly in Phase 2 and found
at or below 1% with no enrichment, despite looking dangerous on inspection.

---

## Appendix A — benign and latent findings

**Real defects that are not causing these failures.** `om.py:440` forces multiplicity to
singlet or doublet for every complex. **35.4% of the non-quick corpus (2,324 molecules) and
36.8% of the quick corpus (9,109 — Ru 2,759, Ni 2,358, Fe 1,052, Cu 857, Mo 650, Co 480,
Mn 354, V 341, Cr 258)** sit on a metal that routinely carries S > ½. This is a genuine
mis-specification, but the FF-floor path never consumes `.UHF`, so it cannot be responsible
for any failure measured here. It becomes live the moment `xtb` is wired in — which is also
when it would start silently producing wrong-spin geometries.

That combination is worth stating plainly for the fix thrust: **the single most-recommended
next step (install a real optimizer) is also what activates the largest latent defect.** They
should land together.

The actinide rewrite (`om.py:445`) is the same shape of latent defect with **zero** exposure
in this corpus.
