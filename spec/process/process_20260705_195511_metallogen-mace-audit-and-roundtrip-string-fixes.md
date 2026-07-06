# Process Document: MetalloGen MACE adoption, geometry-classification audit, and round-trip string fixes

**Generated:** 2026-07-05T19:55:11-07:00
**Session Focus:** Adopt the MACE MLIP optimizer for the MetalloGen 3D backend, audit the
geometry-classification criteria surfaced by the regenerated fixtures, and fix the OIN /
m-SMILES round-trip string errors that are not caused by geometry optimization.

## Problem Statement

OIN-SMILES converts 3D transition-metal-complex geometries to a 1D notation (OIN) and back.
The experimental **MetalloGen** OIN→XYZ backend (on `feature/metallogen-3d-generator`) generates
all 25 verification complexes — including the η/ansa metallocene family the legacy stitch engine
could never build — but the XYZ→OIN→XYZ **round-trip** only reproduced ~15/25 input strings. The
goal was to raise round-trip fidelity by fixing the *string/encoding* errors while treating
geometry-optimization quality as a separate concern.

## Starting State

- Branch `feature/metallogen-3d-generator` at `aa3cfc8` (CO-encoding fix); MetalloGen integrated
  behind the engine seam (`OIN3DGenerator(engine="metallogen")`), FF optimizer default.
- A **MACE MLIP optimizer** and regenerated MACE-geometry fixtures had been added to the working
  tree (uncommitted) by prior/parallel work; fixtures had been moved `tests/integration/*.xyz →
  tests/fixtures/*.xyz`. A 403 MB MACE weight file sat under `models/mace/` (gitignored).
- Two verification harnesses: `verify_xyz_to_oin.py` (pure encoder, hardcoded golden strings) and
  `verify_roundtrip.py` (generate → re-encode, compares normalized strings + coordination-sphere
  RMSD). With the new MACE fixtures the encoder test had regressed to 22/25 and the round-trip to
  15/25.
- Standing constraints: **no push to origin**; commit with `git commit --no-verify` (ruff not
  installed, pre-commit red on the vendored `generator3d` tree).

## Approach & Methodology

Iterative, evidence-first debugging with a strict separation between *string* and *geometry*
failures. For each failing complex we ran a **read-only diagnostic** to locate the exact root
cause before editing, then validated the fix in isolation (fast FF generation for string-only
checks; the slow MACE optimizer only to confirm geometry end-to-end). Work was committed in small,
single-concern commits. The session was framed against a plan file
(`~/.claude/plans/how-much-of-the-breezy-emerson.md`) that was revised as evidence changed the
difficulty estimates. This was not a HACF architect→execute cycle; it was direct diagnosis-driven
implementation, closed out with two forward-looking handoff briefs for the remaining geometry work.

## Steps Taken

1. **Reframed the plan around the MACE run.** The definitive `verification_artifacts_20260705_122853`
   (full 25, MACE) showed 15/25 — but with MACE the coordination-sphere RMSDs were all excellent,
   so almost every remaining failure was a *pure string* error. Categorized the 10 failures by
   scope (in-scope string vs out-of-scope geometry).

2. **Audited the geometry-classification criteria** (user redirect). The regenerated fixtures broke
   3 encoder classifications. A read-only diagnostic (`scratchpad/geo_diag.py`) dumped each
   complex's virtual-atom unit vectors, pairwise angles, and per-candidate Kabsch RMSD:
   - **VOacac2** classified `[V_TBP]` but is really square-pyramidal. Measured apex–basal 104–107°,
     trans-basal 146–152°, no ~180° pair — a *puckered* vanadyl. The idealized `SPY` template
     (basal atoms coplanar at 90°) lost the angular RMSD to `TBP` (0.470 vs 0.540). **Encoder was
     wrong.**
   - **TiCat3/4** classified `[Ti_TET]` (correct — bent metallocene, Cp–Ti–Cp 132°, Me–Ti–Me 90°);
     the hardcoded golden still said `[Ti_TPY]`. **Golden was stale.**
   - Answered the eta-vector question: it is the **haptic-group centroid** (`_reduce_hapticity`,
     one virtual atom per η ring), not per-atom.

3. **Fix A — puckered SPY template** (`oin_aligner.py` `TEMPLATE_SPECS["SPY"]`): re-derived the 4
   basal slots at ~105° apex–basal (z=−0.2588, xy·0.9659), keeping ±x/±y directions. VOacac2 SPY
   RMSD 0.540→0.055 → `[V_SPY]`; FeCO5 (TBP 0.003) and ReF7 (PBP) unaffected.

4. **Fix B — TiCat3/4 goldens** (`verify_xyz_to_oin.py`): updated the stale `expected_oin_string`/
   `expected_smiles` TPY→TET after confirming the encoder output is stable across 4 runs. Encoder
   test back to **25/25**.

5. **Committed the MACE work as its own commit** (user request). Cleared a stale `.git/index.lock`
   (0-byte, 3 h 20 m old, predating the running GitHub Desktop — a crashed-process leftover, not a
   live lock). Verified the 403 MB `*.model` was gitignored and excluded the `tests/backup/` scratch
   dir (added to `.gitignore`). To keep the MACE commit pure, temporarily reverted the 4 golden lines
   so they rode with the classification-audit commit instead.

6. **Fix C1 — stereo carry** (`metallogen_adapter.py` `build_contract_mol`): BDPP/BDNN backbone
   `@/@@` was nondeterministic. Diagnostics showed the embed picks a *random handedness* at sp3
   stereocentres (coord-sphere RMSD ~0.12 doesn't see backbone chirality), and — critically — a
   bare `SetChiralTag` does **not** survive `get_oin_string`'s fragment rebuild (proved by flipping
   a tag and observing no output change). The reliable lever is the *perceive-then-flip* pattern
   `ChiralityRecoveryUtility` already uses for Zone-A P: CIP-label the OIN fragment templates, map
   each sp3-C stereocentre template→contract via the existing substructure match, and flip the tag
   on CIP mismatch (bounded 3-pass). No-op for the 23 non-chiral-C complexes. Validated 8/8
   deterministic over stochastic FF embeds; BDPP/BDNN MACE round-trip PASS.

7. **Fix C2 — Ir binding-atom index** (`xyz2mol.py` `get_oin_string`): fac/mer-Ir re-encoded some
   cyclometalated carbons as `[cH]{N}` instead of `c{N}`. A cascade of diagnostics traced it past
   `generate_robust_smiles` (which emitted the binding C correctly as `[c]`) to the V2.4 sidecar
   w-tag: the binding **LocalIdx** pointed at the wrong ring carbon (atom 4, a CH) instead of the
   deprotonated `[c]` (atom 0). Root cause: `sanitized_mol.GetSubstructMatch(smiles_mol)` returns a
   *wrong automorphism* on the near-symmetric aryl (its two ortho carbons are identical to the
   matcher). Fixed by using RDKit's own **`_smilesAtomOutputOrder`** (symmetry-free), with the
   substructure match kept as a guarded fallback. Deterministic all-`c{N}` now; a **latent
   shared-encoder bug** the input path had only survived by luck. Encoder still 25/25.

8. **Definitive full MACE run** (`verification_artifacts_20260705_162501`): round-trip **19/25**.
   Confirmed every remaining failure is geometry, not string: TiCat1/3/4 + TiCp2Me2 have
   byte-identical strings; BDNN's stereo string is correct but the embed distorted the Pd square
   plane this run (`SPL`→`TPY`, RMSD 1.57, stochastic); TiCat2 is the lone residual string mismatch
   (kekulized Cp) but is also RMSD-blocked.

9. **Scoped the two remaining problems and wrote handoff briefs.** User confirmed `--ensemble-size >1`
   fixes BDNN. A 10-min characterization of the eta "999" (aided by the user adding distinct RMSD
   sentinels to `rmsd_utils.py`) found TiCat1 = **997 (element-key mismatch)**: the input
   distance-based coordination sphere counts the non-bonded ansa-Si (`{Ti,C:12,Si:1}`) while the
   generated bond-based sphere omits it (`{Ti,C:12}`) — a **test-harness** definition inconsistency,
   not eta geometry. Wrote `spec/handoff/HANDOFF_eta_rmsd999.md` and
   `spec/handoff/HANDOFF_BDNN_ensemble_stability.md`.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| Pucker the SPY template to ~105° | Add a topological "no ~180° pair ⇒ reject TBP" tie-break | Smallest change within the existing angular-RMSD framework; matches real square-pyramidal geometry; genuine TBP still wins via its true 180° axial pair |
| Update TiCat3/4 goldens (not the classifier) | "Fix" the classifier to output TPY | The structures *are* pseudo-tetrahedral; the encoder was correct and the golden was stale |
| MACE as its own commit; revert goldens temporarily | Bundle goldens into the MACE commit; partial-file staging | Keeps the large MACE fixture batch a clean, self-contained commit and the classification audit logically separate; avoids fragile `git add -p` |
| Stereo via perceive-then-flip vs template CIP | Enforce handedness in the MetalloGen embed; mirror coordinates | Bare tag flips don't survive `get_oin_string`; the embed fix is deep/geometry; carrying the *encoded* CIP decouples string fidelity from stochastic embed handedness |
| Ir fix via `_smilesAtomOutputOrder` | Strip H in `generate_robust_smiles`; force-aromatize | The binding C was already 0-H; the real bug was a wrong atom-index mapping, and the output order is the exact symmetry-free answer |
| Treat eta "999" as a test-harness fix | Assume deep eta hapticity geometry (Opus-hard) | Measured evidence: strings byte-identical, 997 = coordination-sphere *definition* mismatch (ansa-Si), so the geometry is fine |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| MACE optimizer + deps + fixtures | `pyproject.toml`, `src/oinsmiles/generator3d/*`, `tests/fixtures/*.xyz`, `models/mace/*` | committed `66a2725` |
| SPY template + TiCat goldens | `src/oinsmiles/utils/oin_aligner.py`, `tests/integration/verify_xyz_to_oin.py` | committed `c3d0cff` |
| Stereo carry | `src/oinsmiles/generation/metallogen_adapter.py` (`build_contract_mol`) | committed `56a543a` |
| Binding-atom index mapping | `src/oinsmiles/utils/xyz2mol.py` (`get_oin_string`) | committed `d415c16` |
| Worklog entries | `spec/worklog/NOTES.md` | committed `9004e32`, `1e90be2` |
| Handoff briefs | `spec/handoff/HANDOFF_eta_rmsd999.md`, `HANDOFF_BDNN_ensemble_stability.md` | committed `42f374c` |
| Distinct RMSD sentinels (user) | `tests/integration/rmsd_utils.py` | uncommitted (user's) |

## Results & Outcomes

- Encoder `verify_xyz_to_oin.py`: 22/25 → **25/25**.
- Round-trip `verify_roundtrip.py` (MACE): 15/25 → **19/25** — newly passing fac-Ir(ppy)₃,
  mer-Ir(ppy)₃, VOacac2, BDPP.
- All string-fidelity work complete: every remaining round-trip failure is geometry-dominated
  (TiCat1/3/4 eta coord-sphere definition, TiCp2Me2 rotamer, BDNN stochastic Pd distortion — fixed
  by `--ensemble-size >1`; TiCat2 a residual kekulized-Cp that is also RMSD-blocked).
- Unit suites green throughout: `discover tests/unit` 127 (skip 3), `discover tests` 55.
- Branch not pushed; 7 commits added this session (`66a2725`…`42f374c`).

## How to Reproduce

Prerequisites: `feature/metallogen-3d-generator`, `uv sync`, MACE weights present under
`models/mace/` (gitignored; see `models/mace/README.md`).

1. Encoder check: `uv run python tests/integration/verify_xyz_to_oin.py` → **25/25**.
2. Full round-trip (slow, MACE): `uv run bash tests/run_verification.sh
   --optimizer mace-omol-0-extra-large-1024 --ensemble-size 1` → round-trip **19/25**.
   Run convention: `--ensemble-size 1` first, then `5` only to rescue an ensemble-fixable failure
   (e.g. BDNN). Use `--only <name>` for a single complex; omit `--optimizer` for a fast FF
   string-only check.
3. Per-fix isolation (fast): generate a complex with FF (`optimizer=None`), build the contract mol,
   call `get_oin_string`, and compare `OIN(2)` to `OIN(1)` — string fixes are geometry-quality
   independent, so FF is sufficient to iterate; MACE is only needed to confirm geometry passes.
4. Gotchas: string-fidelity fixes must be validated over *several* stochastic embeds (the embed
   handedness/geometry varies run-to-run); a single lucky run hides nondeterminism.

## Patterns & Lessons

- **Chiral tags don't survive `get_oin_string`'s fragment rebuild** — you must re-perceive
  (`AssignStereochemistry(cleanIt=True, force=True)`) before flipping a tag for it to reach the
  emitted SMILES. This is the same pattern `ChiralityRecoveryUtility` uses for P.
- **Prefer `_smilesAtomOutputOrder` over `GetSubstructMatch`** to map a mol atom to its canonical
  SMILES position — substructure matching can return a wrong automorphism on symmetric substructures.
- **Validate string fixes over many stochastic embeds, not one run** — the MetalloGen embed picks a
  random handedness/geometry each time; determinism is the actual acceptance criterion.
- **Separate string failures from geometry failures early.** With MACE producing good geometry, the
  RMSD-99x sentinels and geo-code mismatches were the tell that remaining issues were geometry/harness,
  not encoding — which reshaped the plan and downgraded the eta work from "hard" to "moderate."
- **Distinct failure sentinels beat a single 999** — the user's per-branch RMSD codes turned a
  10-minute characterization into a precise root cause (element-key mismatch on the ansa-Si).
- **Commit hygiene under a messy tree:** keep a large mechanical batch (MACE fixtures) as its own
  commit and split conceptual fixes out, temporarily reverting a few lines rather than fighting
  partial-file staging.
