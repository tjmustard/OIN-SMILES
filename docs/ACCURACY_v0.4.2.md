# OIN-SMILES v0.4.2 — round-trip accuracy & capabilities

This report accounts for what the v0.4.2 accuracy wave changed for XYZ ↔ OIN round-trip fidelity,
using per-molecule evidence rather than a single headline number. It is the companion to the
`CHANGELOG.md` `[0.4.2]` entry (mechanisms) and `docs/KNOWN_LIMITATIONS.md` (what is out of scope).

These are **pre-capstone, per-phase** results. The integrated regression gate — *no molecule that
passes on the baseline may regress on the merged tree* — is run by the capstone
(`spec/handoffs/v0.4.2/VALIDATION.md`) before v0.4.2 reaches `main`.

## Why there is no single pass-percentage

The continuous round-trip accumulator (`tmCAT-tmPHOTO_xyz_dataset/results-v0.4.0/`) is **not** a valid
floor for a headline percentage, for two independent reasons:

1. **Mixed provenance.** Its reports span several code commits, so a global rate averages current code
   on the hard cases with older code on the easy ones (the exact trap the v0.4.0 wave documented).
2. **`--quick` is a different generator.** The accumulator runs `--quick` (`uff_pool_size=2`,
   `max_attempts=10`, 30 s hard-kill) — a materially weaker generator than the non-quick default
   (`pool=5`, full budget). Most of this wave's "geometry" and "winding" failures **evaporate** under
   the default generator (see below), so a `--quick` percentage measures the harness, not the library.

Accuracy is therefore stated as **named per-molecule round-trip flips** measured on a single clean
commit, `c7edeeb6` (= tag `v0.4.1`). The floor is a **set of molecule IDs** — ~5,960 that pass on
`c7edeeb6` (`spec/handoffs/v0.4.2/baseline_pass_c7edeeb6.txt`) plus per-class failing goldens
(`baseline_fail_c7edeeb6.json`). Counts below are that single-commit floor; the mixed-provenance
backlog is larger and grows as the accumulator runs.

## Per-class before → after (on the `c7edeeb6` floor)

| class | floor | v0.4.2 outcome | phase |
|---|---:|---|---|
| `donor_H_atom_count` | 82 | **+6 flips** (`ARONEA`, `BOXJUU`, `CAHZIX`, `COTXAM`, `HINNOH`, `NORTAP`) via the haptic ring-carbon 0-H lock; the `+3`-atom rows are η-arene/`[CH]` H-fill, and the ammine/nitride notation ambiguity is **latent** — both routed to docs | S1 |
| `H_on_terminal_oxo_imido` | 2 | documented: downstream PuLP donor-charge protonation on the assembled complex, not fragment-fixable → docs | S1 |
| `encode_crash_other` | 4 | **+3 full wins** (`IROXET`, `SUNXAB`, `XEVMAN`) via the `get_lig_mol` charge-sweep fall-through; `ASISAX` (quinoid/stuck-ring) → docs | S3 |
| `kekulize_encode_crash` | 7 | **crash removed, canonical key matches** for `JOTJEK`/`TIYWUV`/`ZENZAW`×2 (residual is S1 `[CH]` H-fill); `KAXVOX`/`KAXWAK` (`=S=` ylide), `LEZWAO` (non-converging) → docs | S3 |
| `macrocycle_perception` | 9 | **re-triaged — 0 were aromatic-perception bugs**: 6 boron atom-stereo → S6b, 2 geometry → S5, 1 E/Z (`XIZXAG`) → S6a | S3 |
| `garbled_aromatic` | 2 | reduced-porphyrin (chlorin) generated-geometry saturation mismatch (`DIXXIS`, `ROJQOY`) → docs | S3 |
| `EZ_bond_stereo` | 23 | **+5 flips** (`AHAZOZ` nitrone, `AFECIZ`, `XIZXAG`, `AYUYIE`, `BOCGEH`); the latent force-DOUBLE crash is fixed and free monodentate-arm E/Z is now enforced. `AROHIA` (no input marker → encoder-side regen asymmetry) → docs/S5 | S6a |
| `atom_stereo` | 11 | **+6 fixes** (`JEKQAS`, `REPZUJ`, `ZORCOA` sulfonimidoyl S; `ORIHUU`, `XILZID` donor-C; `POYJIX` N⁺); re-triaged `KEBBUO` → S7, `SEMTOV`/`VEJXOZ` → S5; `JUCCUH`/`WEDYOU` deferred | S6b |
| `string_mismatch_other` (`[S@SP3]` subset) | 5 of 20 | **[S@SP3] closed**: `BAZMOH`/`HUGSEI`/`LUSKIV`/`YUMPIH` round-trip; `CIDDAU` `@` resolved (leaves to S1 on a `[SH]` count) | S6b |
| `geometry_or_fragment_change` | 29 | **`--quick` artifact**: under the non-quick default, CN4 **15/15**, CN5 **8/11**; no classifier change shipped (a TET/TPY hysteresis was designed then dropped — empty feasible margin) | S5 |
| `winding_flip` | 14 | **`--quick` artifact**: non-quick winding-face **12/14** correct; no winding code shipped; small irreducible residual → docs | S5 |
| `geometry_NON` | 1 | **CN-9 now encodable** — `XERTUK_comp_3` emits `[Y_TCT]` (was `[Y_NON]`); its 104-atom ligand still won't embed (generation limit) → docs | S5 |
| `high_rmsd` | 36 | **FF-floor artifact** — every row is an already-string-correct round-trip that fails only the 1.0 Å gate under FF-only; the `ENABLED_METALS` bond-length expansion is a **validated-negative** (ships pristine) | S7 |
| `timeout` | 339 | mostly **real valence `no_conformers` masked** by the `--quick` 30 s kill (not a pure budget artifact); harness now labels it honestly | S7 |
| `no_conformers` | 115 | genuine valence/perception gaps (not quick-starved) → routed to S1/S3/docs | S7 |
| `gen_exception_other` | 24 | **all `UncoordinatedFragmentError`** — outer-sphere counterions/solvent with no binding slot → docs (representation limit) | docs |
| `carborane_unsupported` | 36 | wontfix — 3c2e cage bonding is outside the two-centre `AC2mol` model → `spec/handoffs/v0.4.2/wontfix-carboranes.md` | docs |
| `rmsd_mapping_failed` | 1 | metric edge case | — |

## New capabilities in v0.4.2

- **CN-9 coordination geometry** — a tricapped-trigonal-prism template (`TCT`) makes 9-coordinate
  metals encodable (`[Y_TCT]`) instead of emitting `g:NON` and crashing generation. (Adds to the CN-8
  `SQA` support from v0.3.5; the geometry-code set is now `LIN TPL TET SPL SPY TBP OCT PBP TPY SQA TCT`.)
- **Free monodentate-arm E/Z is enforced** — the generator now reproduces the E/Z of a free C=N/C=C on
  a monodentate arm (previously only chelate-ring-locked bonds were handled), via a charge-aware
  double-bond promotion that no longer crashes the FF cleanup.
- **Cleaner heteroatom stereo** — spurious high-coordination donor-S and quaternary-N⁺ tags are cleared,
  and genuine sulfonimidoyl-S / donor-adjacent-C handedness is re-oriented against the metal-free
  fragment with a consistent fill-first CIP convention.
- **Two encoder crash paths recovered** — large saturated cages that defeated the charge ladder, and
  stale-aromatic-flag kekulization failures, now encode.
- **Honest failure accounting** — the round-trip harness stamps its run mode and separates FF-floor
  `high_rmsd` and `--quick` timeouts from real defects in its summary.

## The re-framing: the real accuracy-defect surface is small

The single most valuable output of this wave is not the ~25 named flips — it is **proving which
failures are not accuracy defects**. Of the ~750 failing rows on the `c7edeeb6` floor, the large
buckets are:

- **harness/FF artifacts** — `timeout` (339, mostly masked `no_conformers`), `high_rmsd` (36, FF
  geometric floor on string-correct round-trips);
- **representation / perception limits** — `carborane_unsupported` (36), `gen_exception_other` (24,
  all outer-sphere), the genuine part of `no_conformers` (115);
- **misfiled classes** — most `macrocycle_perception`/`garbled_aromatic` and much of
  `geometry_or_fragment_change`/`winding_flip`, which the non-quick default generator resolves.

The genuinely encoder/generator-fixable accuracy classes total on the order of a couple hundred rows,
of which v0.4.2 closes the tractable subset and routes the rest — with a named owner — to docs or a
future notation change. This is why the wave shipped **so little code** (8 source files, +346/−70) for
the size of the backlog: most of the backlog was never an accuracy bug.

## Reproducing / verifying

- Always evaluate conformer-yield-sensitive rows **non-`--quick`** (`--mol-timeout 1800`); `--quick`
  fabricates geometric distortion.
- Verify stereo only through `XYZToSMILES.convert()` (not `get_tmc_mol`, which defaults
  `with_stereo=False`), and run any `smiles_1` diff under **one** pinned rdkit (the two blessed
  versions disagree on `/`\` direction).
- Measure accuracy per-molecule against the `c7edeeb6` floor set, never as a percentage of the
  mixed-provenance accumulator.
- The integrated, cross-phase regression gate is the capstone's job: `spec/handoffs/v0.4.2/VALIDATION.md`.
