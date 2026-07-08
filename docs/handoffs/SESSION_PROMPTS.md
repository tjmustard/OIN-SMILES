# Session launch prompts — round-trip fix wave

Copy-paste prompts to start each of the 6 parallel worktree sessions. Every path
here is absolute so the prompts work regardless of which directory you launch
`claude` from.

- **Main checkout:** `/home/tjmustard/Documents/GitHub/OIN-SMILES`
- **Handoff docs:** `/home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/`
- **Dataset + case registry (gitignored — only in the main checkout):**
  `/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/`
- Unit-suite baseline to preserve: **233 OK / 4 skipped** (`uv run python -m unittest discover tests/unit`).

## Step 1 — create all six worktrees (run once, from the main checkout)

```bash
cd /home/tjmustard/Documents/GitHub/OIN-SMILES
git fetch origin && git checkout main && git pull
for s in donor-h eta-diene aromatic-perception eta-winding metrics stereo; do
  git worktree add "../OIN-SMILES-$s" -b "feature/roundtrip-$s" main
done
git worktree list
```

## Step 2 — in each worktree, sync deps and launch a session

```bash
cd /home/tjmustard/Documents/GitHub/OIN-SMILES-<slug>   # e.g. OIN-SMILES-donor-h
uv sync
claude
```

Then paste the matching prompt below. Each prompt assumes the session's working
directory is that worktree. The dataset/registry are read from the main checkout
by absolute path (they are gitignored, so not present in the worktree).

---

## S1 — donor-h  (worktree `OIN-SMILES-donor-h`, branch `feature/roundtrip-donor-h`)

```
You are fixing the "bare anionic donor protonation + unbound-fragment IndexError"
class of tmCAT/tmPHOTO round-trip failures (100 cases, the 2nd-largest bucket).

Read these first, in order:
  1. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/S1-donor-h.md   (your handoff — mission, proven root cause, evidence, acceptance)
  2. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/README.md        (shared protocol: env pins, private-output-dir rule, squash-PR flow)

Your full case list is the "S1-donor-h" rows of:
  /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/20260707-results/CASE_REGISTRY.md

You OWN only `convert_parsed_to_msmiles` in src/oinsmiles/generation/metallogen_adapter.py
(plus a new tests/unit/test_bare_donor_hydrogens.py). Do NOT edit the template
functions in that file — S2 owns them.

Reproduce a case (write to a PRIVATE output dir, never the shared results dir):
  uv run python tools/test_dataset_roundtrip.py \
    --dataset-dir /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset \
    --output-dir /tmp/rt-donor-h --quick \
    --only WAYHOW_comp_0,UDIVUY_comp_0,XADYAC_comp_0,FENMIX_comp_0,KIZQER_comp_0

The handoff proves the N-donor root cause by code-read, but VERIFY before coding
(the INENOF acetylide is a second, separate H leak). Then fix, add guard tests
for both the strip cases and the dative counter-cases that must keep H
([NH2]/[NH]/[OH2]), and confirm the full suite stays green (233 OK / 4 skipped).
Land via squash-PR to origin/main per docs/handoffs/README.md.
```

## S2 — eta-diene  (worktree `OIN-SMILES-eta-diene`, branch `feature/roundtrip-eta-diene`)

```
You are fixing the "eta-2 alkene/diene (COD) double-bond localization" class of
tmCAT/tmPHOTO round-trip failures (74 cases; drives Rh's 66% failure rate).

Read these first, in order:
  1. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/S2-eta-diene.md
  2. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/README.md

Your full case list is the "S2-eta-diene" rows of:
  /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/20260707-results/CASE_REGISTRY.md

You OWN only `_flatten_template`, `build_contract_mol`, and `_oin_fragment_templates`
in src/oinsmiles/generation/metallogen_adapter.py (plus a new
tests/unit/test_contract_mol_diene_transfer.py). Do NOT edit convert_parsed_to_msmiles
— S1 owns it. Regression floor: tests/unit/test_contract_mol_allyl_transfer.py must stay green.

Reproduce a case (PRIVATE output dir):
  uv run python tools/test_dataset_roundtrip.py \
    --dataset-dir /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset \
    --output-dir /tmp/rt-eta-diene --quick \
    --only GASBIN_comp_0,PENGAT_comp_0,ABIRIO_comp_0

This class burned two prior handoffs with wrong root-cause guesses — follow the
handoff's VERIFY-FIRST steps (instrument build_contract_mol's GetSubstructMatch)
before writing code. Add a COD-transfer guard test, keep the allyl tests green,
full suite 233 OK / 4 skipped. Land via squash-PR to origin/main.
```

## S3 — aromatic-perception  (worktree `OIN-SMILES-aromatic-perception`, branch `feature/roundtrip-aromatic-perception`)

```
You are fixing the "aromatic/quinoid/macrocycle perception" class of tmCAT/tmPHOTO
round-trip failures — the LARGEST bucket (158 cases), dominated by porphyrinoid
macrocycles that re-encode with inconsistent kekulization, plus mixed-aromatic
garbling and hard kekulize/NoneType encode crashes.

Read these first, in order:
  1. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/S3-aromatic-perception.md
  2. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/README.md

Your full case list is the "S3-aromatic-perception" rows of:
  /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/20260707-results/CASE_REGISTRY.md
(classes: garbled_aromatic, macrocycle_perception, kekulize_encode_crash, xyz2mol_none_crash)

You OWN src/oinsmiles/utils/xyz2mol.py and src/oinsmiles/generator3d/process.py
(plus a new tests/unit/test_aromatic_reencode.py). Read-only: oin_aligner.py (S4),
metallogen_adapter.py (S1/S2). Regression floor: tests/unit/test_quinoid_ligand_parse.py.

Reproduce a spread (PRIVATE output dir):
  uv run python tools/test_dataset_roundtrip.py \
    --dataset-dir /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset \
    --output-dir /tmp/rt-aromatic --quick \
    --only KIYWUM_comp_0,TIPDIG_comp_0,ABERIK_comp_0,NAXDOI_comp_0

Start with forward-encode STABILITY (encode the same input twice → identical
string) before chasing the re-encode. Porphyrins are the bulk — a consistent
macrocycle kekulization likely clears many at once. Add guard tests; crashes must
become specific errors, not bare tracebacks. Full suite 233 OK / 4 skipped.
Land via squash-PR to origin/main.
```

## S4 — eta-winding  (worktree `OIN-SMILES-eta-winding`, branch `feature/roundtrip-eta-winding`)

```
You are fixing the "automorphic-ring eta-winding flip + eta-slot placement" class
of tmCAT/tmPHOTO round-trip failures (52 cases). A fully symmetric ring (Cp*,
C6Me6) traversed CW vs CCW is the SAME structure, but the encoder emits {n>} vs
{n<} nondeterministically. Winding on SUBSTITUTED rings is load-bearing (TiCat
rac/meso) and must stay distinct.

Read these first, in order:
  1. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/S4-eta-winding.md
  2. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/README.md

Your full case list is the "S4-eta-winding" rows of:
  /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/20260707-results/CASE_REGISTRY.md

You OWN src/oinsmiles/utils/oin_aligner.py and src/oinsmiles/oin/compare.py
(plus a new tests/unit/test_automorphic_ring_winding.py). Read-only: xyz2mol.py (S3).
Regression floor: EVERY existing winding/eta test + tests/unit/test_roundtrip_canonical_key.py,
especially the TiCat3/4 rac/meso guards.

Reproduce (PRIVATE output dir):
  uv run python tools/test_dataset_roundtrip.py \
    --dataset-dir /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset \
    --output-dir /tmp/rt-winding --quick \
    --only CAHKEE_comp_0,SOJMIQ_comp_0

Prefer an ENCODER canonical tie-break (deterministic winding for automorphic
rings) over a comparator patch. Add a symmetric-ring stability test AND a
substituted-ring counter-case that must still distinguish </>. Verify on BOTH
rdkit 2025.09.3 and 2026.3.3. Full suite 233 OK / 4 skipped. Land via squash-PR.
```

## S5 — metrics  (worktree `OIN-SMILES-metrics`, branch `feature/roundtrip-metrics`)

```
You are fixing the "RMSD 996/999 sentinel + harness robustness" class of
tmCAT/tmPHOTO round-trip failures (62 cases). 996/999 are SENTINELS from
rmsd_utils.py (coordination-sphere per-element count mismatch / hard failure),
NOT bad geometry — a correct structure can hit 996. Also in scope: a hard
subprocess watchdog (SIGALRM can't kill native-code hangs like UGUHAH).

Read these first, in order:
  1. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/S5-metrics.md
  2. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/README.md

Your full case list is the "S5-metrics" rows of:
  /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/20260707-results/CASE_REGISTRY.md

You OWN tests/integration/rmsd_utils.py and tools/* (plus a new
tests/unit/test_rmsd_mapping.py). Read-only: everything under src/oinsmiles/.
NOTE: tools/test_dataset_roundtrip.py was just changed on main (commit d950f2a) —
rebase before opening your PR. A continuous runner executes it from the MAIN
checkout, so your worktree edits don't interfere until merge.

Reproduce (PRIVATE output dir):
  uv run python tools/test_dataset_roundtrip.py \
    --dataset-dir /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset \
    --output-dir /tmp/rt-metrics --quick \
    --only DAPZIF_comp_0,CAWYOR_comp_0,ABETIK_comp_0

First give sentinels an honest label (distinct from real high_rmsd); then improve
the mapping (element-agnostic eta slots, distance-anchored fallback). Use MEAN
RMSD, never max-per-atom. Full suite 233 OK / 4 skipped. Land via squash-PR.
```

## S6 — stereo  (worktree `OIN-SMILES-stereo`, branch `feature/roundtrip-stereo`)

```
You are fixing the "P/N donor @-stereo through generation + generator E/Z
honoring + no-conformer corners" class of tmCAT/tmPHOTO round-trip failures
(~82 cases incl. no-conformer triage). @-stereo is lost/flipped through the
stereo-blind embed; the generator ignores C=C E/Z (task #8).

Read these first, in order:
  1. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/S6-stereo.md
  2. /home/tjmustard/Documents/GitHub/OIN-SMILES/docs/handoffs/README.md

Your full case list is the "S6-stereo" and "S6-stereo (triage)" rows of:
  /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/20260707-results/CASE_REGISTRY.md

You OWN src/oinsmiles/core/chirality.py and src/oinsmiles/generator3d/ligand.py
(and, if needed, the embed path in generator3d/embed.py) + new tests
test_donor_stereo_roundtrip.py / test_generator_ez_honoring.py. If a fix needs
build_contract_mol stereo stamps, coordinate with S2 in your PR body rather than
editing that function. Regression floor: test_zone_a_p_donor_stereo.py,
test_backbone_heteroatom_stereo.py, test_double_bond_stereo_encoding.py.

Reproduce (PRIVATE output dir):
  uv run python tools/test_dataset_roundtrip.py \
    --dataset-dir /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset \
    --output-dir /tmp/rt-stereo --quick \
    --only QEXKOT_comp_0,SUXJOL_comp_0,XIZXAG_comp_0

CRITICAL convention (handoff §prior-art): stereo labels must come from
rdCIPLabeler on the metal-free template, NOT legacy Chem.AssignStereochemistry —
they disagree for 3-coordinate P and the wrong source round-trips the wrong
enantiomer. Verify @-stereo across ≥3 embed seeds and both rdkit versions. Full
suite 233 OK / 4 skipped. Land via squash-PR.
```

---

## Notes for whoever coordinates the wave

- **Sequencing:** S3 (158) and S1 (100) are the heavy hitters; start there if you
  are not running all six at once. S3's porphyrin fix likely clears many cases at once.
- **Shared-file pair:** S1 and S2 both edit `metallogen_adapter.py` (disjoint
  functions). Whichever lands first, the other rebases onto new main before its PR.
- **After any PR lands,** the other live sessions `git rebase origin/main`.
- The gitignored `CASE_REGISTRY.md` is regenerated any time by:
  `python tools/classify_failures.py --output-dir /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/20260707-results`
