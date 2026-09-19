# Fault-attribution census — what to measure before v0.4.17

> **Status:** APPROVED by the owner 2026-09-17. Planning artefact — no census number exists yet.
> Session handoff (gitignored, main checkout): `spec/handoffs/census/NEXT.md`.

## Context

Every release since v0.4.4 has chased **one suspected cause at a time** through one
instrument: the round trip `XYZ → E → OIN → G → XYZ' → E → OIN'`, scored by string compare.
A round-trip failure is a statement about the *composition* plus the comparator. It cannot say
which factor failed:

| | Component | A fault here looks like |
|---|---|---|
| **P**  | perception (XYZ → graph, charges, bond orders) | wrong graph, encoded faithfully |
| **E1** | first encode | lossy / non-canonical / non-injective string |
| **G**  | 3D generation | wrong graph, wrong hand, detached ring, or nothing |
| **E2** | re-encode of XYZ' | not invariant to the new atom order / conformer / H count |
| **C**  | comparator + metric | false pass, false fail |
| **H**  | harness | seed, load, timeout |

Note the structural confound: **G always emits a different atom order and a different
conformer than the input**, so every round trip is *implicitly* a renumbering + conformer
invariance test of E. A failure blamed on G may be E2 alone.

**Deliverable of this plan is not v0.4.17.** It is one per-molecule **attribution table** over
the frozen n=5000 cohort where every molecule — *passes included* — carries an independent
verdict per component. The next release is then chosen from a measured partition.

## What the survey found (3 read-only agents, this session)

**Already measured, keep:** `hard_fail` 266 = G produces nothing · `facmer` 15 = arrangement
only · `encode_fail` 12 root-caused · `slot_renumber` 201 = mirror image of the input ·
`DETACHED` 301 real (18.2× enriched), not selection-reachable (0/289) · E is rotation-invariant
(0 drift, n=225).

**The unknowns — nothing on disk answers these today:**

1. **How many of the 1,142 failures have a CORRECT generated structure?** No instrument judges
   G without routing through E. `coordination_report` (distances only) is the sole exception
   and covers the metal sphere only. No graph-isomorphism check exists. Generated XYZ is on
   disk for ~864 of the 1,142, so this is answerable **offline**.
2. **How many of the 3,858 passes are wrong?** FP rate was 9.6% (28.1% haptic) under the old
   scoring; under the honest metric it has **never been measured at corpus scale**. The key
   deliberately folds reflection, metal descriptors, H count. The un-run experiment is named in
   `docs/agentic-notes/v0.4.13/LANE-02-attach-classes.md` §6 (48 `byte_exact`+`DETACHED`).
3. **E's renumbering instability was last measured at v0.4.5: ~30% byte-drift, 5.4%
   key-broken** (n=298) — larger than the whole 22.84-pt gap — and has **never been crossed
   per-molecule against round-trip outcome**. The status doc itself says "100% round-trip
   accuracy is not currently well-defined". Coordinate-noise invariance: never tested.
   Conformer invariance: asserted; CREST tool exits 0 when CREST is absent; no numbers in tree.
4. **Generator "non-determinism" has never been sampled.** The harness never passes `seed=`;
   every tier runs `seed=42` (`generation/engine.py:32`). Nobody knows whether a G failure is
   *luck* (passes on seed 43) or *capability*. 77.16% has no seed error bar.
5. **`structural` is a residual bucket** ("none of the above" in
   `tools/roundtrip_bucket_report.py::classify`), not a diagnosis. Its 141 "PERCEPTION"
   molecules: perceiver-vs-geometry blame is explicitly unmeasured.
6. **No perception ground truth** beyond total `Charge:` in the XYZ comment line.
7. **No corpus-scale collision scan** (distinct molecules → same string; same molecule →
   different strings). Injectivity work covered 299 structures + 4 fixtures.
8. **Every headline since v0.4.6 is the same 5,000 molecules**; 11 levers were tuned on them.
   The other 20,197 have never seen a modern encoder/generator. Overfit is unmeasured.

## The census — four instruments, one join

All new code in `tools/census/`, on a `research/census` worktree (pre-push hook blocks it),
run with main checkout's `.venv/bin/python` + `PYTHONPATH=$PWD/src`. Outputs to
`tmCAT-tmPHOTO_xyz_dataset/results-census/` (gitignored) → frozen with `/freeze-measurements`
(extend the harvester ALLOWLIST). Notes → `docs/agentic-notes/v0.4.17/`.
Inputs: `results-v0.4.14-sweep/{summary_roundtrip.json, bucket_report_honest.json,
attach_class_audit.json, structures/*_generated.xyz}`. **Never re-score from
`structures/*.oin`** (falls back to `smiles_1`); use the JSON record.

**I1 — The neutral ruler + G-vs-input (offline, <1 CPU-h).** One perception-free,
encoder-free comparison applied identically to input XYZ and generated XYZ:
element-labelled distance adjacency (`perception_core.get_AC`/`xyz2AC_vdW` radii step only —
no bond orders, no charges) + metal contacts (`oin/coordination.py::metal_contacts`).
Per molecule: formula / H-count diff · heavy-atom graph isomorphic? · sphere verdict
(`coordination_report`, already stored) · chirality SAME / MIRROR / DIFFERENT
(`tools/injectivity/oracle.py::is_distinct_enantiomer`, `config_oracle.py::
configurational_signature`; ligand tetrahedral centres via metal-free fragments +
`AssignStereochemistryFrom3D`, which is automorphism-safe). Symmetric molecules where the
mapping is ambiguous are labelled `UNDECIDED`, never guessed.

**I2 — E self-consistency on the 5,000 inputs (no G; ~15 CPU-h detached, ~2.5 h on 6 cores —
forecast from `results-v0.4.8-honest/encoder_identity.jsonl`: 1.51 CPU-h per full encode pass,
median 0.24 s).** Extend `tools/canonicality_probe.py` (has identity/rotate/renumber, `--shard`,
`--n 0`) with `noise` (σ≈0.02 Å) and `mirror` modes and per-molecule JSONL. Per molecule:
deterministic? · renumber byte/key stable? · noise byte/key stable? · mirror: chiral-by-oracle
⇒ string MUST change; achiral ⇒ MUST NOT.

**I3 — String sufficiency, no 3D (offline, minutes).** (a) Parse-back: `smiles_1` →
`generation/oin_parser.OINParser` → `metallogen_adapter.build_contract_mol` → neutral graph;
isomorphic to the input's neutral graph? Isolates serializer/parser from embedding — matters
most for the 266 `hard_fail` that have no structure. (b) Collision scan with the I1 ruler:
same string / different molecule (injectivity) and same molecule / different string
(canonicality across independent crystal structures — `_comp_n` siblings and the 1,033
cat/photo duplicates are free natural twins). (c) Perception plausibility flags on the input:
charge sum vs comment line, metal oxidation state in the known set, valence anomalies.

**I4 — The join (`tools/census/attribution_table.py`).** One TSV, 5,000 rows: bucket,
attach class, I1/I2/I3 verdicts, features (metal, CN, hapticity, n_atoms, elapsed). Fault
assigned by first match:

| # | Condition | Fault |
|---|---|---|
| 1 | `smiles_1` None | P/E1 coverage |
| 2 | parse-back graph ≠ input | E1 serializer / parser |
| 3 | no generated structure | G-nothing (timeout vs no-conformer from `error`) |
| 4 | neutral graph ≠ | G-construction (DETACHED / wrong connectivity) |
| 5 | graph =, chirality MIRROR/DIFFERENT, and I2 says E *does* distinguish the mirror | G-handedness |
| 5b | same, but I2 says E(mirror)=E(x) | **E1 non-injective — G was never told** |
| 6 | graph = and chirality = but round trip FAILED | **E2 / P non-robust — not G's fault** (sub-split by I2 flags and by input-vs-generated perceived bond-order diff + `structure_distortion_report.py` → settles the 141 PERCEPTION blame) |
| 7 | round trip PASSED but graph ≠ or chirality ≠ | **metric false pass** |
| 8 | else | `UNATTRIBUTED` — reported as a number, never folded away |

**Every instrument ships with its broken-instrument control** (the project's own hard lesson):
I1 input-vs-input ⇒ 100% SAME; input-vs-mirrored-input ⇒ every oracle-chiral molecule MIRROR;
input-vs-ligand-deleted ⇒ graph ≠. I2 identity mode ⇒ 100% stable; a known v0.4.5 drifter must
drift. Every report prints its denominator and its `UNDECIDED`/`UNATTRIBUTED` count.

## Two small fresh-run probes (after the join, because the join picks their strata)

**S1 — Seed variance.** Add `--seed` to `tools/test_dataset_roundtrip.py` (research branch
only). ~120 molecules stratified by I4 fault class (pass / G-handedness / G-construction) ×
5 seeds, `--gen-timeout 60`, quiet box, BLAS=1. Yields P(pass | seed) per class = luck vs
capability, and the seed error bar on 77.16%.

**S2 — Held-out 1,000 (optional, ~8–10 CPU-h detached).** `tools/build_sweep_cohort.py`,
seed 43, excluding every molecule in any existing cohort. Shipped defaults. SE ≈ 1.3 pts.
If it reads well below 77%, the levers are fitted to the cohort.

## Second machine: the Mac Mini (available, not yet set up)

**Not needed for the core census.** C1/C3/C4 are offline and take minutes; I2 is ~2.5 h wall on
this box. Setting the Mac up costs a session, so **defer it until after C4** — it earns its
keep only for the fresh generator runs (S1, S2), where a quiet dedicated box matters because
*load biases accuracy* here (advisory timeout ⇒ starved pool ⇒ understated `byte_exact`).

**The trap: a Mac result is NOT like-for-like with the Linux baseline.** arm64 + Accelerate
BLAS vs x86 + OpenBLAS ⇒ the distance-geometry embed can differ at `seed=42`; and **g-xTB may
have no macOS build**, in which case the harness silently swaps PASS 2 to `FF_reroll`
(`_pass2_config`). So on the Mac only *within-machine* comparisons are valid:
- S1 seed variance ✅ (seed 42 re-run there as its own control arm).
- S2 held-out ✅ **only with a calibration arm**: ~300 molecules of the 5k cohort run on the
  Mac too, to measure the platform offset before the held-out number is compared to 77.16%.

**A free bonus unknown (#9): is the hash platform-independent?** For a canonical 1D hash this
is a requirement, and it has never been tested. One encode pass of the cohort on the Mac
(~1.5 CPU-h), diffed byte-for-byte against the Linux `smiles_1`. Any diff is a finding.

**M0 — readiness probe (one short session, only when we reach C5):** macOS Remote Login + key
auth so I can drive it over `ssh`; Mac **fetches from this box over SSH** (pull direction — the
`pre-push` hook blocks `research/*` and `main` is do-not-push, and `origin` is stale at
v0.4.13); `rsync` the gitignored `cat/`+`photo/` (125 MB) and the cohort manifest; `uv` +
rdkit `==2025.9.3` arm64 wheel; check xtb / g-xTB availability and record `optimizer_effective`;
no `systemd-run` on macOS ⇒ `nohup` + `caffeinate -i` + the same `#DONE` sentinel. Gate: unit
suite green + 62-fixture encoder byte-identity vs Linux, *then* the cohort encode diff.

## What each answer changes

- Rule 6 large ⇒ next release is **encoder invariance**, and v0.4.17's CONSTRUCTION charter is
  partly mis-aimed. Rule 5b large ⇒ the 201 enantiomers are an E problem, not a G problem.
- Rule 7 non-trivial ⇒ the headline is over-stated again; **gate the harness on I1 first**.
- I2-unstable molecules enriched among failures ⇒ G work is unmeasurable until E is canonical.
- S1 shows failures pass on other seeds ⇒ the lever is pool/budget, not construction.
- Rules 3–5 dominate and the rest are small ⇒ the current roadmap is confirmed — with evidence.

## Chunks (one Claude session each; CPU jobs run detached *between* sessions and cost no usage)

| Chunk | Work | Claude cost | CPU |
|---|---|---|---|
| **C0** | On approval: copy this plan to `docs/agentic-notes/v0.4.17/CENSUS_PLAN.md` + a `spec/handoffs/census/NEXT.md` handoff on `research/census`; save one memory | tiny | 0 |
| **C1** | I1 ruler + controls + run on 4,748 + cross-tab vs bucket | medium | <1 h |
| **C2** | I2: extend probe, smoke n=20, launch detached (`systemd-run` service unit, `#DONE` sentinel), end session | small | ~15 h |
| **C3** | I3 parse-back + collision scan + P flags | medium | minutes |
| **C4** | I4 join, fault partition, `CENSUS.md`, freeze, re-point `ROADMAP_100_100.md` | medium | minutes |
| **M0** | Mac Mini readiness probe + cross-platform hash diff (only if using the Mac for C5/C6) | small | ~2 h (Mac) |
| **C5** | S1 seed probe: add flag, launch; read next session | small | few h (Mac or here) |
| **C6** | S2 held-out + 300-mol calibration arm if on the Mac (optional) | tiny | ~10–13 h (Mac) |

Usage rules: each session opens from `NEXT.md`, not from the docs tree (this survey is not
repeated); no subagents or workflows needed — each chunk is one script; never poll a running
job, wait on the `#DONE` artifact; C1 and C3 are independent of C2's long run. C1 alone
answers unknown #1 and #2 and is the natural stopping point if usage is tight.

## Verification

- Controls above pass before any corpus number is read (a control failure voids the run).
- I4 row count = 5,000; fault classes sum to 1,142 failures + 3,858 passes; bucket column
  reproduces `bucket_report_honest` exactly (3858/365/15/484/266/12).
- Cross-checks against existing measurements: rule 4 ∩ `attach_class_audit` DETACHED ≈ 301;
  rule 5 ∩ `slot_renumber` ≈ the 201 `MIRROR_MATCH` from `veto_residue_chirality.py`. A
  disagreement is a finding about one of the two instruments, to be resolved before the join.
- `find tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k -xtype l | wc -l` = 0 before every run.
- Unit suite untouched (tools only); `uvx ruff@0.15.20` clean; commit normally, never push.
