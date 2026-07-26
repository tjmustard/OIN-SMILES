# v0.4.7 "Slow-100" Cohort and Byte-Identity Gate

Lane **L1-cohort**. This document freezes the target cohort for the v0.4.7 performance
wave and describes the two-arm byte-identity gate every other v0.4.7 lane gates its
changes against.

## 1. Selection method

`tools/select_slow_byte_exact.py` applies a fixed predicate over a sweep results
directory's `individual_reports/*.json`:

```
report["status"] == "success"
report["tier_passed"] == "UFF_1"
report["smiles_1"] and report["smiles_2"]
report["smiles_1"].strip() == report["smiles_2"].strip()      # byte-exact, NOT key-equal
isinstance(report["metrics"]["elapsed_s"], (int, float))
-> sort by metrics.elapsed_s DESC, take top N
```

`elapsed_s` lives **inside** `report["metrics"]`, not at the top level; reading it from
the top level silently returns a default (typically `0`) for every row and produces a
cohort that is an artifact of file-listing order, not slowness. The tool guards this
explicitly (excludes rather than defaults a malformed/missing value) and exits non-zero
on an empty input corpus or on 0 rows passing the predicate.

Bands are fixed on the *selection* `elapsed_s`: **A** > 200 s, **B** 100-200 s,
**C** < 100 s. `eta` flags a haptic OIN (`re.search(r"\{\d+[<>]\}", smiles_1)`).

**What this predicate certifies, and what it does not.** `smiles_1 == smiles_2` is a
**notation-level** guarantee: the two OIN *strings* are byte-identical. It does **not**
certify that the generator reproduced the same conformer, nor that a full re-perception
of the generated structure would still see the same coordination. A sibling v0.4.7 lane
measured exactly this gap: `OIN_ACCEPT_SCORED` can accept a *different* conformer that
happens to emit a byte-identical OIN string, and on 6 molecules a full re-perception of
that generated structure no longer sees the haptic ligand as coordinated at all. This
cohort's selection predicate — and the gate built on it — is deliberately notation-level
(matching this project's own round-trip contract, `oin/compare.py`), and that is the
correct contract for a *performance* wave whose job is "did the string change," not "is
the chemistry still right." But a future reader should not read "byte-exact round trip"
as a stronger, structure-level guarantee than it is.

## 2. Provisional cohort construction

Built from the two available v0.4.5 result dirs (the live `results-v0.4.6-sweep` was
explicitly excluded from cohort construction — it is incomplete and mid-config-change):

| label | dir | total reports | passed byte-exact predicate |
|---|---|---:|---:|
| A | `results-v0.4.5-sweep-partial-2697mols` | 2697 | 2253 |
| C | `results-v0.4.5-rebaseline` | 936 | 564 |

**Dataset-level overlap is small: only 100 of the 25,197-unique-basename dataset were
independently re-measured in both A and C.** This is not a bug — two uncorrelated random
draws of sizes 2697 and 936 from a 25,197-name universe are *expected* to share
`2697 * 936 / 25197 ≈ 100.1` names, and 100 is exactly what was found. `A` and `C` are
genuinely independent samples, not nested subsets of one frozen cohort.

Of those 100 shared names, **62 pass the byte-exact predicate in both A and C**
(the other 38 failed, errored, or landed on a different tier in at least one run). The
cohort is built from `min(elapsed_A, elapsed_C)` over those 62 — the "slow in both
independent runs" rule that screens out single-run contention flukes. Because only 62
candidates exist at all, every one of them clears a top-100 cutoff trivially: **there is
no 100th name to cut.**

**Cohort size actually achieved: 62, not 100.** This is reported honestly rather than
padded with single-run-only "slow" molecules pulled from outside the 100-name overlap —
doing so would reintroduce exactly the contention-fluke risk the both-runs rule exists to
remove, and would silently misrepresent what was actually double-confirmed.

**eta cross-check:** the haptic flag (computed from `smiles_1`) agreed between A and C
for all 62 overlapping molecules — 0 mismatches.

**Overlap with the live, still-growing `results-v0.4.6-sweep`** (1963/5000 done at time
of check): 50 of the 62 cohort molecules already have a report there. Informational only
— not used to build the cohort, since that sweep is mid-config-change.

### Band / eta split (on the cohort's selection `elapsed_s` = `min(elapsed_A, elapsed_C)`)

| band | definition | count |
|---|---|---:|
| A | > 200 s | 2 |
| B | 100-200 s | 3 |
| C | < 100 s | 57 |

eta fraction: **15 / 62 (24.2%)**.

**Caveat worth flagging plainly:** because the cohort is capped by the small
double-measured overlap rather than by an absolute speed floor, band C spans the full
range down to ~1.8 s — some cohort members are not "slow" by any absolute standard, only
slow-and-double-confirmed relative to what little overlap exists. See open questions.

### Full cohort (62 molecules), sorted by selection `elapsed_s` DESC

| molecule | elapsed_A (s) | elapsed_C (s) | selection elapsed_s = min | band | eta |
|---|---:|---:|---:|:-:|:-:|
| FUXMAN_comp_0 | 299.440 | 243.193 | 243.193 | A | yes |
| KELQOI_comp_0 | 227.559 | 276.877 | 227.559 | A | - |
| IJOSOP_comp_0 | 186.884 | 180.807 | 180.807 | B | - |
| JIVFOK_comp_0 | 167.199 | 209.020 | 167.199 | B | yes |
| FIBQAK_comp_0 | 184.385 | 127.646 | 127.646 | B | yes |
| GULPAG_comp_0 | 64.374 | 93.503 | 64.374 | C | yes |
| IJAXIB_comp_0 | 70.724 | 58.355 | 58.355 | C | - |
| AKOXII_comp_0 | 82.660 | 56.305 | 56.305 | C | yes |
| MIKDOA_comp_0 | 42.560 | 46.205 | 42.560 | C | yes |
| CIGLAI_comp_0 | 35.179 | 35.767 | 35.179 | C | - |
| IFIHAI_comp_0 | 30.787 | 42.441 | 30.787 | C | yes |
| FOJWOR_comp_0 | 26.575 | 29.695 | 26.575 | C | yes |
| GACTAG_comp_0 | 35.847 | 25.798 | 25.798 | C | yes |
| IHOCUE_comp_0 | 25.727 | 26.423 | 25.727 | C | yes |
| FILVON_comp_0 | 23.939 | 40.419 | 23.939 | C | - |
| DEZPUV_comp_0 | 23.939 | 22.094 | 22.094 | C | - |
| FOGCIN_comp_0 | 32.572 | 21.593 | 21.593 | C | - |
| DIJDAD_comp_0 | 19.932 | 19.864 | 19.864 | C | - |
| HORJAB_comp_0 | 14.835 | 15.014 | 14.835 | C | - |
| FIVKIG_comp_0 | 11.828 | 11.579 | 11.579 | C | - |
| HEBBAS_comp_0 | 10.846 | 23.037 | 10.846 | C | - |
| GUKDAS_comp_0 | 9.259 | 9.727 | 9.259 | C | - |
| BARVUO_comp_0 | 11.685 | 8.786 | 8.786 | C | - |
| FOJHIW_comp_0 | 12.397 | 8.708 | 8.708 | C | yes |
| MEKWEE_comp_0 | 8.558 | 12.781 | 8.558 | C | - |
| MUQLIV_comp_0 | 6.462 | 8.188 | 6.462 | C | - |
| AGELIJ_comp_0 | 15.422 | 6.442 | 6.442 | C | - |
| HINFOC_comp_0 | 8.621 | 6.259 | 6.259 | C | - |
| LOWXAW_comp_0 | 5.593 | 5.407 | 5.407 | C | yes |
| MULDAZ_comp_0 | 5.391 | 8.229 | 5.391 | C | - |
| CEYKOH_comp_0 | 5.920 | 5.373 | 5.373 | C | - |
| MAQYEL_comp_0 | 5.129 | 5.853 | 5.129 | C | - |
| MABKOQ_comp_0 | 5.035 | 4.835 | 4.835 | C | - |
| ETORIQ_comp_0 | 4.820 | 5.014 | 4.820 | C | - |
| FAHXUI_comp_0 | 6.966 | 4.730 | 4.730 | C | - |
| IKIJIX_comp_0 | 10.348 | 4.662 | 4.662 | C | - |
| BOLGIU_comp_0 | 5.564 | 4.595 | 4.595 | C | - |
| ECEXIU_comp_0 | 5.547 | 4.564 | 4.564 | C | - |
| JIMMAU_comp_0 | 4.761 | 4.558 | 4.558 | C | - |
| JEMVED_comp_0 | 4.485 | 4.816 | 4.485 | C | - |
| HAGKAD_comp_1 | 7.268 | 4.438 | 4.438 | C | - |
| LUVMIA_comp_0 | 5.779 | 4.161 | 4.161 | C | - |
| CUJJOH_comp_0 | 6.169 | 3.797 | 3.797 | C | - |
| CERQAT_comp_0 | 3.669 | 5.948 | 3.669 | C | - |
| JEGKEM_comp_1 | 3.288 | 3.358 | 3.288 | C | yes |
| LIMRIJ_comp_0 | 3.165 | 6.725 | 3.165 | C | - |
| HAZCIX_comp_0 | 3.119 | 4.044 | 3.119 | C | - |
| KOVXIE_comp_0 | 3.086 | 3.180 | 3.086 | C | - |
| NIKSUY_comp_0 | 5.907 | 3.066 | 3.066 | C | - |
| JANVUQ_comp_0 | 2.965 | 3.713 | 2.965 | C | yes |
| HIXJEF_comp_0 | 2.952 | 7.816 | 2.952 | C | - |
| IQEGEQ_comp_0 | 2.923 | 2.928 | 2.923 | C | - |
| COHJEQ_comp_0 | 5.374 | 2.771 | 2.771 | C | - |
| HIKJIW_comp_0 | 4.385 | 2.634 | 2.634 | C | - |
| BIMDOQ_comp_0 | 2.399 | 3.032 | 2.399 | C | - |
| DIPZAD_comp_0 | 2.379 | 2.885 | 2.379 | C | - |
| NOVMAO_comp_0 | 2.250 | 3.857 | 2.250 | C | - |
| JABDEW_comp_0 | 2.367 | 2.246 | 2.246 | C | - |
| LECSUJ_comp_0 | 2.240 | 2.216 | 2.216 | C | - |
| CAJROV_comp_0 | 2.057 | 2.605 | 2.057 | C | yes |
| FEJFAD_comp_0 | 2.009 | 2.019 | 2.009 | C | - |
| MUTYEG_comp_0 | 1.817 | 3.394 | 1.817 | C | - |

## 3. Cohort materialization

`tools/build_sweep_cohort.py` gained a `--names-file` mode (mutually exclusive with the
existing random `--n/--seed` sampling; both modes keep its refusal to overwrite an
existing `--out` dir). It resolves every requested name under `--subdirs` (default
`cat,photo`, dedup priority in that order) and exits non-zero naming exactly which names
were not found, rather than silently materializing a short cohort.

```
PYTHONPATH=src .venv/bin/python tools/build_sweep_cohort.py \
    --names-file cohort_names.txt \
    --dataset-dir tmCAT-tmPHOTO_xyz_dataset \
    --out tmCAT-tmPHOTO_xyz_dataset/cohort-v047-slow100 \
    --overlap-with tmCAT-tmPHOTO_xyz_dataset/results-v0.4.6-sweep \
    --overlap-with tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-sweep-partial-2697mols \
    --overlap-with tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-rebaseline
```

Result: **62/62 symlinks created**, 0 missing, at
`tmCAT-tmPHOTO_xyz_dataset/cohort-v047-slow100/` (gitignored, like all dataset dirs). A
symlink dir is used, not the raw tree, because the raw dataset has 1,033 basenames
duplicated across `cat/` and `photo/` — the established fix for the resulting
report-write race (same pattern as `cohort-v0.4.5-5k` and the v0.4.4
`regression_inputs/` cohort).

`tools/perf_attribute.py::find_xyz` previously searched only
`<dataset>/{cat,photo}` and could not see cohort symlink dirs. It now also globs
`<dataset>/cohort-*/` so other v0.4.7 lanes can point it at this cohort directly.

## 4. The byte-identity gate: `tools/gate_v047.sh`

```
tools/gate_v047.sh arm1 [--fixtures-dir DIR] [--golden PATH] [--out PATH]
tools/gate_v047.sh arm2 [--cohort-dir DIR] [--golden PATH] [--timeout S] [--out PATH]
tools/gate_v047.sh both ...
```

**The gate object is `sha256(smiles_1)` / `sha256(smiles_2)` — the OIN *string* — never
the generated XYZ.** `tools/perf_byte_identity_ab.py` shas the XYZ, which is strictly
stronger and would fail a lane for legitimately picking a different-but-equivalent
conformer. The generated XYZ's sha is recorded as an OBSERVATION column only in ARM 2's
output, never gated on.

**For other lanes reusing this cohort as a byte-gate directly** (e.g. an
`OIN_ACCEPT_SCORED`-style A/B): `tools/gate_v047_arm2_golden.tsv` records BOTH
`sha256(smiles_1)` AND `sha256(smiles_2)` per molecule, keyed by name, columns 2 and 3
respectively (`name<TAB>sha1<TAB>sha2<TAB>len1<TAB>len2<TAB>eta`) — a lane whose own
arm produces no output for some molecules can diff its own two arms' shas against this
manifest without re-deriving anything, and the frozen 62-name list
(`docs/COHORT_v0.4.7.md` §2 table) is exactly the population the byte-exact selection
predicate already vetted, so it doubles as a denominator check for "did my arm produce
output at all" as well as "is the string unchanged."

### ARM 1 — encoder, 61 fixtures

`tools/gate_arm1_encode.py` iterates every `tests/fixtures/*.xyz` file (asserted to be
exactly 61 at runtime — a fixture added/removed is itself gate-relevant drift and fails
loudly rather than silently comparing against a different corpus). For each: encode,
clear the AC2BO memo *between* molecules (copied discipline from
`tools/enc_byte_identity_ab.py` — a cross-molecule cache hit must never be what makes two
revisions agree), emit `name<TAB>sha256(oin)<TAB>len<TAB>eta`. Errors are part of the
contract too (`ERROR:<Type>:<msg>` in the sha column) — two revisions must raise the SAME
error. Output is sorted by name, followed by `# MANIFEST_SHA256=<...>` and `#DONE 61`.

**ARM 1 run to completion this session** (`PYTHONPATH=src .venv/bin/python
tools/gate_arm1_encode.py --fixtures-dir tests/fixtures`):

```
# molecules=61 fixtures_dir=.../tests/fixtures
# MANIFEST_SHA256=373cc387bae9a38c665bd8cbe4b5023682b933802b607499c645e78aa13aaf69
#DONE 61
```

Frozen as the committed golden at `tools/gate_v047_arm1_golden.tsv` (61 data rows,
sorted by name, plus the `# MANIFEST_SHA256=...` and `#DONE 61` lines). A future
`tools/gate_v047.sh arm1` run diffs against this file and fails loudly on any
byte difference, in either the manifest hash or any individual row.

### ARM 2 — round trip, the frozen 62-molecule cohort

`tools/gate_arm2_roundtrip_one.py` replicates the exact UFF_1-tier pipeline
(`optimizer=None, ensemble_size=1, ff_params=None`, the tier every cohort molecule is
known to have passed at, byte-exact, in two independent runs) for ONE molecule per
subprocess invocation. `gate_v047.sh arm2` loops over the cohort, invoking a fresh
interpreter per molecule — several `OIN_*` levers and module-level caches
(`generator3d/clash.py:VDW_ACCEPTANCE_ENABLED`, PuLP's topology memo, etc.) are frozen at
import time or live for the process's lifetime, so a fresh interpreter per molecule is
the only isolation guarantee that does not depend on enumerating every such cache
correctly.

**The ARM 2 golden manifest was NOT built by running the full round trip for all 62
molecules today** (`tools/gate_v047_build_arm2_golden.py`) — that would cost the exact
wall-clock this cohort was selected to be expensive at, which is the *quiet-phase sweep's*
job, not this lane's (env rule: keep any generation to a handful of molecules). Instead
the golden values are read directly from the two already-completed v0.4.5 sweep JSON
reports (`smiles_1`/`smiles_2`), re-asserting byte-exact agreement between A and C rather
than re-trusting the selection step blindly — 62/62 agreed, 0 mismatches.

**Mechanics smoke-tested on 3/62 molecules** (`MUTYEG_comp_0`, `CAJROV_comp_0`,
`FEJFAD_comp_0` — the three fastest in the cohort, ~2 s each): a fresh
`gate_arm2_roundtrip_one.py` invocation reproduced `sha256(smiles_1)` and
`sha256(smiles_2)` **byte-identical to the golden manifest** for all three. This proves
the gate mechanics (encode → generate → re-encode → hash) and the "golden sourced from
existing JSON" approach agree with an actual fresh re-run.

### Mandatory discipline built into both arms (each learned the hard way on this project)

- Every emitted line uses `print(..., flush=True)` (or a real `>>` file append at the
  shell level) — Python block-buffers stdout when not a tty, so a `timeout` kill mid-run
  can otherwise discard the buffer while a downstream `sort` still exits 0 on an empty
  file, silently looking like agreement.
- The `#DONE <n>` sentinel is required and its denominator is checked BEFORE any
  comparison is trusted, both in the encode arm (`#DONE 61`) and the round-trip arm
  (`#DONE <cohort size>`).
- `2>/dev/null` is never used anywhere in the gate.
- ARM 2 spawns one subprocess per molecule rather than looping in one interpreter.

## 5. Open questions / risks for the wave

- **62, not 100.** If a fuller "Slow-100" cohort is wanted, the only path is more
  independent-run overlap — re-run this selection once `results-v0.4.6-sweep` completes
  (5000 molecules; 50/62 of the current cohort already appear there), either as a third
  arm or in place of one of the two v0.4.5 arms.
- **Band C is not uniformly slow.** 57/62 cohort members are `< 100s` by the selection
  metric, and the fastest are only ~2s. The double-confirmation constraint, not an
  absolute speed floor, is what capped this cohort — worth a product decision on whether
  band C should be trimmed to only genuinely slow members (at the cost of a smaller
  cohort) once more overlap is available.
- **ARM 2's golden manifest certifies stability from today forward, not correctness
  against an independent oracle** — it was built from this project's own prior sweep
  output, not validated externally. A "gate passes" result means "unchanged since v0.4.5",
  not "was independently verified."
- **The gate is notation-level, not structure-level** (see §1) — a sibling lane found
  `OIN_ACCEPT_SCORED` can accept a different conformer under a byte-identical OIN string,
  and on 6 molecules a full re-perception of the generated structure no longer sees the
  haptic ligand as coordinated at all. A "PASS" on this gate means the string this lane
  emits is unchanged; it does not mean the 3D structure behind it is. That is the correct,
  deliberately-scoped contract for this performance wave, not a gap in this gate.
- `tests/fixtures/BENVOG_comp_0.xyz` (one of the 61 ARM 1 fixtures) is a known slow
  encode — resonance search bounded by a CPU-time `RLIMIT_CPU` fallback (not wall-clock),
  so ARM 1's wall-clock runtime is load-dependent even though its byte-identity outcome
  is not. Do not read a slow ARM 1 run under host contention as a regression signal by
  itself.
