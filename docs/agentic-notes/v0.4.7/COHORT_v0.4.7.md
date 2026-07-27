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

## 2. Cohort construction

### Revision note — the first build of this cohort was wrong, and here is why

The first pass at this section intersected two results dirs
(`results-v0.4.5-sweep-partial-2697mols`, 2697 mols, and `results-v0.4.5-rebaseline`, 936
mols) under a "slow in both independent runs screens out single-run contention flukes"
rule, on the assumption that both dirs were repeat measurements of the *same* population.
**They are not.** `results-v0.4.5-rebaseline` is the "gap ∪ guard" cohort from the v0.4.5
re-baseline; `results-v0.4.5-sweep-partial-2697mols` is a partial of the frozen seed-42 5k
cohort. They are different draws that happen to overlap by only ~100 names (consistent
with two uncorrelated random samples of sizes 2697 and 936 out of a 25,197-name universe:
`2697 * 936 / 25197 ≈ 100.1`, and 100 is exactly what was found). Intersecting them was
therefore a near-total sample destruction, not a robustness filter: it produced a
62-molecule cohort with 57/62 members under 100s and a fastest member at 1.8s — useless
for a wave whose entire purpose is the slow tail. The selection *arithmetic* was correct
and the overlap finding was correctly reported; the two-dir design itself was the mistake,
caught and corrected before this cohort was used by any other lane.

### Corrected method: top-N from ONE dir, corroboration only where it happens to exist

**The cohort is the top 100 molecules by `metrics.elapsed_s`, byte-exact selection
predicate, from `results-v0.4.5-sweep-partial-2697mols` ALONE:**

```
PYTHONPATH=src .venv/bin/python tools/select_slow_byte_exact.py \
    --results-dir tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-sweep-partial-2697mols \
    --n 100 \
    --corroborate-with tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-rebaseline
```

| source | total reports | passed byte-exact predicate | top-100 cutoff | max |
|---|---:|---:|---:|---:|
| `results-v0.4.5-sweep-partial-2697mols` (primary) | 2697 | 2253 | **93.06 s** | 299.44 s |

`results-v0.4.5-rebaseline` is used **only** as an optional corroboration column
(`--corroborate-with`, added to `tools/select_slow_byte_exact.py` in this revision) — for
a name that *happens* to also have a passing report there, its `elapsed_s` is recorded
alongside as a second observation; absence is expected and unremarkable, and it is never
a requirement for inclusion. Of the 100 selected, **5 have a corroborating observation**
in the rebaseline dir (`FUXMAN_comp_0`, `KELQOI_comp_0`, `IJOSOP_comp_0`, `FIBQAK_comp_0`,
`JIVFOK_comp_0`) — all 5 agree in relative order-of-magnitude with the primary
measurement, and all 5 have byte-identical `smiles_1` between the two sources (an
encoder-determinism sanity check the tool also performs: 0 mismatches, `eta` flag also
agrees on all 5).

### This is provisional — the live v0.4.6 sweep is the SAME frozen cohort as the primary source

`results-v0.4.6-sweep` (still running, 5000 mols) draws from the **same frozen seed-42 5k
cohort** as `results-v0.4.5-sweep-partial-2697mols` — not an independent draw. At 2440/5000
reports it already gives a top-100 cutoff of **97.4 s**, close to this build's 93.06 s, and
the two lists are expected to largely agree once it completes. **This build is therefore
explicitly provisional**: the quiet phase should re-derive the cohort from the completed
5000-molecule sweep and diff the resulting top-100 against this one (names and the 93.06 s
cutoff) rather than silently replacing it — a large divergence would itself be a finding
worth chasing, not just noise to overwrite.

### Band / eta split

| band | definition | count |
|---|---|---:|
| A | > 200 s | 27 |
| B | 100-200 s | 67 |
| C | < 100 s | 6 |

eta fraction: **60 / 100 (60%)**. This is the shape the wave was planned around: overwhelmingly
slow-tail molecules (94/100 at or above 100s), not a speed-agnostic sample.

### Full cohort (100 molecules), sorted by `elapsed_s` DESC

`elapsed_s` is from the primary source (`results-v0.4.5-sweep-partial-2697mols`) and is
the value the cohort was selected and banded on. `corroboration_s` is the optional second
observation from `results-v0.4.5-rebaseline` where one exists (`-` otherwise) — reported,
never filtered on.

| molecule | elapsed_s (s) | corroboration_s (s) | band | eta |
|---|---:|---:|:-:|:-:|
| FUXMAN_comp_0 | 299.440 | 243.193 | A | yes |
| CAHQEJ_comp_0 | 296.549 | - | A | yes |
| GAVSED_comp_0 | 288.116 | - | A | yes |
| BOMNOJ_comp_0 | 286.945 | - | A | yes |
| LUSBEI_comp_0 | 284.478 | - | A | - |
| EBUFOX_comp_0 | 274.977 | - | A | yes |
| KAQDEL_comp_0 | 264.796 | - | A | yes |
| DILZOQ_comp_2 | 264.126 | - | A | yes |
| MIQWEO_comp_0 | 257.492 | - | A | yes |
| CIMZAZ_comp_0 | 250.978 | - | A | yes |
| AHUKIZ_comp_0 | 247.294 | - | A | yes |
| KUJGEC_comp_0 | 243.326 | - | A | yes |
| FUWMAN_comp_0 | 235.510 | - | A | - |
| LUZDAL_comp_0 | 233.048 | - | A | yes |
| MECXEY_comp_0 | 231.608 | - | A | yes |
| DEMSOG_comp_0 | 228.478 | - | A | - |
| KELQOI_comp_0 | 227.559 | 276.877 | A | - |
| LIRFOJ_comp_0 | 223.991 | - | A | - |
| JOCWOP_comp_0 | 221.211 | - | A | - |
| KADVER_comp_1 | 214.444 | - | A | yes |
| AJOKIV_comp_0 | 213.951 | - | A | yes |
| ICUYUC_comp_0 | 213.361 | - | A | yes |
| OROJUB_comp_0 | 211.477 | - | A | - |
| EZURIB_comp_0 | 210.799 | - | A | - |
| DAWBOT_comp_0 | 207.944 | - | A | - |
| DIPPAW_comp_0 | 205.147 | - | A | - |
| MUHGIG_comp_0 | 202.427 | - | A | yes |
| LUMDON_comp_0 | 199.711 | - | B | - |
| HIJPUM_comp_0 | 198.136 | - | B | yes |
| CINWUT_comp_0 | 190.243 | - | B | - |
| DEFPOW_comp_0 | 189.616 | - | B | - |
| PIXGIP_comp_0 | 188.347 | - | B | yes |
| IJOSOP_comp_0 | 186.884 | 180.807 | B | - |
| JUXPID_comp_0 | 186.494 | - | B | yes |
| FIBQAK_comp_0 | 184.385 | 127.646 | B | yes |
| HEXFAS_comp_0 | 184.034 | - | B | yes |
| BIKRUL_comp_0 | 181.859 | - | B | yes |
| BOPZOX_comp_0 | 180.918 | - | B | - |
| KAHNEO_comp_0 | 180.404 | - | B | yes |
| IFICAD_comp_0 | 176.215 | - | B | yes |
| IBIZUP_comp_0 | 176.081 | - | B | - |
| COMXUY_comp_0 | 175.401 | - | B | yes |
| CISDUF_comp_0 | 174.046 | - | B | yes |
| MOCHAN_comp_0 | 172.745 | - | B | yes |
| MUBNAY_comp_0 | 171.242 | - | B | yes |
| DERMIX_comp_0 | 169.092 | - | B | yes |
| LOMFUQ_comp_0 | 168.169 | - | B | - |
| JIVFOK_comp_0 | 167.199 | 209.020 | B | yes |
| NEDCII_comp_0 | 161.650 | - | B | yes |
| ERATIA_comp_0 | 159.460 | - | B | - |
| FOJJUM_comp_0 | 159.435 | - | B | - |
| DOKROM_comp_0 | 156.323 | - | B | - |
| APAGOO_comp_0 | 156.114 | - | B | yes |
| FIXDAS_comp_0 | 155.923 | - | B | yes |
| OTIGIJ_comp_0 | 154.065 | - | B | - |
| IBAZIW_comp_0 | 153.445 | - | B | yes |
| DOGQEX_comp_0 | 152.639 | - | B | - |
| JEDPIU_comp_0 | 151.087 | - | B | yes |
| NICWUS_comp_0 | 150.539 | - | B | - |
| OCUKAA_comp_0 | 147.456 | - | B | yes |
| CISCUE_comp_0 | 144.896 | - | B | - |
| EYOCEC_comp_0 | 143.125 | - | B | - |
| NONGAA_comp_0 | 141.422 | - | B | yes |
| KIPBAP_comp_1 | 140.749 | - | B | - |
| HAMGEH_comp_0 | 137.409 | - | B | - |
| HIPQIH_comp_0 | 135.617 | - | B | yes |
| NOYREA_comp_0 | 134.689 | - | B | - |
| FOKQOM_comp_0 | 134.472 | - | B | yes |
| EHORID_comp_0 | 133.332 | - | B | - |
| HAMGAD_comp_0 | 127.861 | - | B | - |
| HOHTAB_comp_0 | 126.838 | - | B | yes |
| EWUQOC_comp_0 | 125.671 | - | B | yes |
| DIYWAM_comp_0 | 124.374 | - | B | yes |
| EZUROH_comp_0 | 124.321 | - | B | - |
| AGIKUW_comp_0 | 122.149 | - | B | yes |
| LARQED_comp_0 | 121.851 | - | B | yes |
| KIRVOY_comp_0 | 119.818 | - | B | yes |
| NOVYUU_comp_0 | 116.933 | - | B | yes |
| ENUGAX_comp_0 | 116.627 | - | B | - |
| DIDROZ_comp_0 | 116.237 | - | B | - |
| JAFFOL_comp_0 | 115.815 | - | B | yes |
| DERMET_comp_0 | 115.005 | - | B | yes |
| GIZQEN_comp_0 | 113.276 | - | B | - |
| IZIMOU_comp_0 | 112.933 | - | B | yes |
| IHOYAH_comp_0 | 112.389 | - | B | - |
| OVUFUI_comp_0 | 112.306 | - | B | yes |
| GIMWAA_comp_0 | 111.422 | - | B | yes |
| BITREE_comp_0 | 110.029 | - | B | - |
| GAMFOS_comp_0 | 106.674 | - | B | - |
| LUZGIY_comp_0 | 106.547 | - | B | yes |
| ADANEB_comp_0 | 105.927 | - | B | yes |
| KIKGAP_comp_0 | 105.890 | - | B | - |
| NISHAZ_comp_0 | 105.703 | - | B | yes |
| CAZXIM_comp_0 | 104.442 | - | B | - |
| MEVFAU_comp_0 | 99.925 | - | C | yes |
| CAHQAF_comp_0 | 99.676 | - | C | yes |
| INUWUL_comp_0 | 99.495 | - | C | yes |
| IXUTAW_comp_0 | 94.627 | - | C | yes |
| ILUDID_comp_0 | 93.848 | - | C | - |
| AJIJUY_comp_0 | 93.060 | - | C | yes |

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

Result: **100/100 symlinks created**, 0 missing, at
`tmCAT-tmPHOTO_xyz_dataset/cohort-v047-slow100/` (gitignored, like all dataset dirs). A
symlink dir is used, not the raw tree, because the raw dataset has 1,033 basenames
duplicated across `cat/` and `photo/` — the established fix for the resulting
report-write race (same pattern as `cohort-v0.4.5-5k` and the v0.4.4
`regression_inputs/` cohort).

Overlap report from the same command: **100/100** with the primary source (by
construction), **92/100** with the live `results-v0.4.6-sweep` (confirms it is drawing
from the same frozen cohort as the primary source, not an independent population), and
**5/100** with `results-v0.4.5-rebaseline` (the corroboration-only source, consistent
with the ~100-name total overlap between the two v0.4.5 dirs discussed in §2).

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
manifest without re-deriving anything, and the frozen 100-name list
(`docs/agentic-notes/v0.4.7/COHORT_v0.4.7.md` §2 table) is exactly the population the byte-exact selection
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

### ARM 2 — round trip, the frozen 100-molecule cohort

`tools/gate_arm2_roundtrip_one.py` replicates the exact UFF_1-tier pipeline
(`optimizer=None, ensemble_size=1, ff_params=None`, the tier every cohort molecule is
known to have passed at, byte-exact) for ONE molecule per subprocess invocation.
`gate_v047.sh arm2` loops over the cohort, invoking a fresh interpreter per molecule —
several `OIN_*` levers and module-level caches
(`generator3d/clash.py:VDW_ACCEPTANCE_ENABLED`, PuLP's topology memo, etc.) are frozen at
import time or live for the process's lifetime, so a fresh interpreter per molecule is
the only isolation guarantee that does not depend on enumerating every such cache
correctly.

**The ARM 2 golden manifest was NOT built by running the full round trip for all 100
molecules today** (`tools/gate_v047_build_arm2_golden.py`) — that would cost the exact
wall-clock this cohort was selected to be expensive at, which is the *quiet-phase sweep's*
job, not this lane's (env rule: keep any generation to a handful of molecules). Instead
the golden values are read directly from the primary source's already-completed v0.4.5
sweep JSON reports (`smiles_1`/`smiles_2`) — this tool now takes ONE required primary
source (must pass the full predicate for every cohort name, or the build fails loudly)
and an OPTIONAL corroboration source (informational only, never blocking): 100/100
resolved from the primary, 5/100 additionally cross-checked against the corroboration
source, 0 mismatches. Frozen at `tools/gate_v047_arm2_golden.tsv`
(`# MANIFEST_SHA256=6f61359bc61af48943e773031f913698005c9c02ef2d893511e36325fc4ab794`,
`#DONE 100`).

**Mechanics smoke-tested on 2/100 molecules** (`AJIJUY_comp_0`, `ILUDID_comp_0` — the two
fastest in the cohort, ~93 s each — kept to 2 rather than 3 this time since every cohort
member now genuinely costs ~93-300s, unlike the earlier mis-built cohort's ~2s fastest
members): a fresh `gate_arm2_roundtrip_one.py` invocation reproduced `sha256(smiles_1)`
and `sha256(smiles_2)` **byte-identical to the golden manifest for both**
(`AJIJUY_comp_0` -> `76aa9bde...0108d`; `ILUDID_comp_0` -> `3b88e003...9965`, both
matching sha1==sha2==golden exactly). This proves the gate mechanics (encode → generate →
re-encode → hash) and the "golden sourced from existing JSON" approach agree with an
actual fresh re-run.

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

- **This is a provisional cohort, and it says so explicitly (see §2).** The primary
  source (`results-v0.4.5-sweep-partial-2697mols`, 2697/5000 of the frozen seed-42
  cohort) is a partial view; the live `results-v0.4.6-sweep` is the SAME frozen cohort,
  further along (2440/5000 at time of writing, top-100 cutoff already 97.4s vs. this
  build's 93.06s). **Action for the quiet phase:** once `results-v0.4.6-sweep` completes,
  re-run `tools/select_slow_byte_exact.py --results-dir results-v0.4.6-sweep --n 100` and
  diff its name list and cutoff against this one (§2 table, 93.06s). Expected to largely
  agree; a large divergence would be a finding, not noise to silently overwrite this
  cohort with.
- **`results-v0.4.5-rebaseline` corroborates only 5/100 molecules.** This is expected and
  fine — it is a different, smaller, differently-drawn cohort (see the revision note in
  §2) and was never going to cover most of a top-100-by-elapsed_s drawn from a different,
  much larger population. The 5 that do overlap all agree (elapsed order-of-magnitude,
  byte-identical `smiles_1`, matching `eta`), which is the corroboration this column
  exists to provide — not a claim that the other 95 are somehow less trustworthy.
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
