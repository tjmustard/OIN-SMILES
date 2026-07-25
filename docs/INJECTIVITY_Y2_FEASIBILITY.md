# Injectivity Audit (Y2) — Feasibility & Fix-vs-Document Decision

**Status:** measurement / feasibility wave — **no encoder code changed yet.** Reproducer:
`tools/injectivity/config_oracle.py`. Parent: `docs/INJECTIVITY_Y1_OVERVIEW.md`.

## TL;DR — the headline reframe

Wave 1 confirmed three encoder blind spots (metal Δ/Λ, axial atropisomer, metal-bound 2°
amine) and filed them, per prior belief, as *documented-permanent* limitations. **Wave 2
overturns that framing.** For all three axes the configuration is **still present in the input
coordinates** — the encoder discards a signal it has, it does not lack one. A per-axis
descriptor recovered straight from the 3D geometry **flips between enantiomers** on every
confirmed fixture:

| axis | fixture | recovery route | base → mirror | recoverable? |
|---|---|---|---|---|
| **P1 metal Δ/Λ** | fac-Ir(ppy)₃ | RDKit `AssignStereochemistryFrom3D` → `@OH` permutation | **10 → 8** | ✅ yes (RDKit) |
| **P2 axial** | PdCl₂-R-BINAP | signed biaryl dihedral (ours) | **−75° → +75°** | ✅ yes (manual) |
| **P3 bound amine** | POJJOP | signed tetrahedral volume at locked N (ours) | **−9.4 → +9.4** | ✅ yes (manual) |

Control: cisplatin recovers **nothing** (achiral). So the question for each axis moves from
*"can the configuration be perceived?"* (**yes, all three**) to *"can it be emitted as a
**canonical** token the generator can **round-trip**?"* — a canonicalization + generator
problem, not a perception one.

## What the RDKit Book settled (perception routes)

- **Non-tetrahedral (metal) stereo from 3D is "Complete"** in RDKit ≥2022.09. The metal atom
  gets a `CHI_OCTAHEDRAL` / `CHI_SQUAREPLANAR` tag and a `_chiralPermutation` index; the index
  differs for Δ vs Λ (fac-Ir(ppy)₃: 10 vs 8). **P1 is perceived for free.** Caveats the Book
  flags as *"Totally missing"* for non-tet: **CIP assignment, canonicalization, SMILES-writing
  canonicalization** — i.e. RDKit hands us the raw permutation but not a canonical, orientation-
  invariant descriptor. That gap is ours to close.
- **Atropisomers are NOT perceived from pure 3D.** The Book (§"3D Coordinates for Input of
  Atropisomers"): *"the atropisomer bond is marked only if one of the neighbor bonds is
  wedged/hashed … the actual configuration is [then] determined by the 3D coordinates."* With no
  wedge, `FindPotentialStereo` returns nothing (confirmed on the BINAP ligand — 0 elements). So
  the **plan's assumed P2 fix mechanism ("RDKit atropisomer perception") does not fire from our
  geometry.** The recoverable signal is the raw **signed biaryl dihedral**, which we compute
  ourselves — no dependence on RDKit's wedge machinery.
- **Metal-bound amines are cleared by RDKit** as invertible (`CHI_UNSPECIFIED` even after
  `AssignStereochemistryFrom3D`), matching the encoder's own Zone-A clear. Recoverable only
  manually, as the **signed tetrahedral volume** of the metal-locked N's four neighbours.

## The configurational oracle (`tools/injectivity/config_oracle.py`)

Wave 1's oracle answered a scalar *"is the mirror distinct?"* (rigid RMSD). Wave 2 needs the
next thing up — *"what is the configuration?"* — a descriptor that flips for an enantiomer:

- `metal_stereo_descriptors(mol)` → `MetalStereo(shape, permutation)` from RDKit 3D perception.
- `axial_axes(mol)` → `AxialAxis(sign, dihedral, hindered)`; `hindered` is a cheap
  twist-plus-ortho-wall candidate filter (the **axis-selection** side), the sign is the
  **configuration** side. Planar chelated biaryls (ppy) are correctly rejected.
- `bound_amine_centers(mol)` → `AmineCenter(sign, volume)`; requires four **symmetry-distinct**
  neighbours, so NH₃ ammines and primary amines (equivalent H's) are correctly excluded.
- `mirror_flip_report(xyz)` → per-axis boolean that the descriptor flips for the z-mirror — the
  constructive proof of recoverability. Guarded by `tests/unit/test_config_oracle.py` (7 tests).

**These are distinguishers, not yet canonical tokens.** The metal permutation is relative to
RDKit's neighbour order; the axial/amine signs are relative to a canonical-rank reference
neighbour. Two *orientations of the same isomer* are not yet guaranteed one descriptor — that
orientation-invariance is exactly the canonicalization the fix must add (and overlaps v0.4.5).

## The central trade-off (why this is a decision, not just an implementation)

Per the Wave-1 user calibration, **a round-trip FAIL with a correct OIN is usually the
generator, not the notation.** That cuts directly against emitting these tokens:

> Today each blind spot is a **False Positive** — a lossy OIN that *passes* round-trip silently.
> The moment the encoder emits a distinguishing stereo token, the molecule can only round-trip
> if **the generator reproduces that stereochemistry**. MetalloGen does not build metal Δ/Λ,
> atropisomer, or locked-amine configuration today. So emitting a token **converts a silent
> False Positive into a loud False Negative** — the OIN becomes *more correct* (injective), but
> the headline round-trip pass-rate **drops**, and the new FAILs are generator-caused.

That is arguably the *right* trade for a losslessness project — but it is a product call, and it
couples the encoder change to generator work. It should not be made silently inside the encoder.

## Fix-vs-document recommendation (per axis)

| axis | perception | canonicalization cost | generator cost | recommendation |
|---|---|---|---|---|
| **P1 metal Δ/Λ** | free (RDKit) | **high** — needs canonical donor ordering (the v0.4.5 problem) | high — MetalloGen can't set OH permutation | **Document as *recoverable, deferred to v0.4.5-canonical*** — not permanent. Ship the oracle now; emit `@OHn`/`@SPn` only once canonical ordering + generator land together. |
| **P2 axial** | manual (dihedral) | **medium** — sign needs a canonical reference-neighbour rule; axis *selection* heuristic needs hardening | medium — generator must honour a biaryl dihedral sign | **Prototype the encoder emit** behind an opt-in flag (raw-string token, key folds it) as the lowest-risk fix candidate; gate ON only paired with generator support. Most self-contained of the three. |
| **P3 bound amine** | manual (signed volume) | **medium** — needs canonical neighbour order + a Zone-A override for the metal-locked case | high — generator can't set locked-N config | **Document as *recoverable, deferred*** alongside P1; the ~2.85 % motif blast-radius is real but the fix needs both a Zone-A carve-out and generator support. |

Net: **P2 is the lead fix candidate**, P1 and P3 are re-labelled *recoverable-but-deferred*
(not permanent), and all three now have an independent configurational measurement to build a
future fix against and to verify emitted tokens with.

### Decision taken (2026-07-24)

Per user direction, the recommendation above is adopted: **P2 opt-in emit implemented; P1 + P3
relabelled recoverable-but-deferred to v0.4.5** (their docs and `KNOWN_LIMITATIONS.md` updated).

The P2 emit ships behind the **`OIN_EMIT_AXIAL`** flag (default OFF → encoder output byte-identical,
zero regression). With the flag on, `get_oin_string` appends an axial-sign token computed from the
signed biaryl dihedral of each hindered axis (`src/oinsmiles/oin/axial.py`, the single source of
truth shared with the oracle): R-BINAP → ` |ax:-|`, S-BINAP → ` |ax:+|`. The round-trip key folds
the token (`oin/compare.py::_AXIAL_TOKEN_RE`), so the batch harness is unaffected whether or not the
flag is set — it stays a silent false-positive at the *key* level, but the *raw* string is now
injective on the axis. Molecules with no hindered biaryl axis are byte-identical even with the flag
on. Guard: `tests/unit/test_axial_emit.py`.

**Not yet done (gating ON):** the sign is still relative to a canonical-rank reference neighbour
(a distinguisher, not a fully canonical reference-free token), and MetalloGen does not reproduce the
axis — so turning the flag on by default would convert the silent false-positive into a
generator-caused false-negative. Gating ON waits on both a canonical token and generator support.

## Reproduce

```
# from the Wave-2 worktree root, main venv (no uv sync; rdkit pinned 2025.9.3):
PYTHONPATH=$PWD/src python -m tools.injectivity.config_oracle \
  tests/fixtures/fac-Ir\(ppy\)3.xyz tests/fixtures/PdCl2-R-BINAP.xyz \
  tests/fixtures/POJJOP.xyz tests/fixtures/CisPlatin.xyz
PYTHONPATH=$PWD/src python -m unittest tests.unit.test_config_oracle
```

Expected: metal flips only for fac-Ir(ppy)₃, axial only for BINAP, amine only for POJJOP,
nothing for cisplatin; tests `OK` (7).

## P2 close-out (2026-07-24) — the axis now round-trips

The deferral above rested on a premise that is **no longer true**: "the generator cannot reproduce
the axis, so emitting a token converts a silent false positive into a false negative." It can now.

**Measured A/B** (`tools/injectivity/axial_roundtrip_ab.py`, seed 42, FF optimizer), encoding each
BINAP enantiomer with `OIN_EMIT_AXIAL=1`, generating 3D from the resulting OIN, and reading the
generated structure's own axial token:

| arm | R (base) | S (mirror) | matched |
|---|---|---|---|
| baseline — no axial-aware pass | ✅ `-` | ❌ got `-`, wanted `+` | **1 / 2** |
| axial-aware selection | ✅ `-` | ✅ `+` | **2 / 2** |

The baseline returns the *same* handedness whichever enantiomer is requested — its single match is
luck, not fidelity. Note this is **n = 2** (one fixture, both enantiomers): a directional result,
not a rate.

### What made it work

1. **A canonical token.** The reference ortho neighbour is now chosen by *symmetry* rank
   (`CanonicalRankAtoms(breakTies=False)`), never tie-broken rank. The token is invariant under
   input atom renumbering and under any proper rotation, and flips only under reflection — which
   is what makes it comparable between the requested OIN and a freshly embedded conformer with
   entirely different atom numbering. (Honest note: the previous tie-broken code also survived
   permutation on BINAP; that hazard was latent, not live.)
2. **A stereogenicity gate — the real correctness win.** A ring end whose two ortho neighbours are
   symmetry-*equivalent* carries a local C2 through the axis, so the molecule is achiral however
   twisted it is. Without the gate the encoder emits a sign for achiral biaryls: *over-sensitivity*,
   the fourth cell of the Y1 confusion matrix. The independent geometric oracle agrees with the gate
   exactly — mirror RMSD **0.000 Å** (superimposable by proper rotation) on symmetric-end cases vs
   2.80 Å on the true atropisomer.
3. **Generator support by selection, not construction.** `_select_by_geometry` narrows the pool to
   conformers whose own token matches the request, following the project's thrice-confirmed
   "selection beats construction" pattern rather than trying to constrain a single-bond torsion in
   distance geometry. ETKDG has no atropisomer bias (21/19 on the free ligand), so the right
   conformer is always in the pool.

### Three integration bugs — and the hardening against their class

Each is repaired, but the durable problem was not the individual defect: two of the three
produced a *plausible wrong answer* with no error. A wrong atropisomer still satisfies the
round-trip key (the key folds the axial token), so nothing downstream complains — the
"self-consistent but wrong" failure this whole audit exists to hunt, reproduced inside the fix
for it. Three guards now make that class visible (`tests/unit/test_axial_failure_modes.py`):

- `_axial_narrow` returns `(kept, n_blind)`, separating *perceived-but-different* from
  *could not perceive at all*. `_axial_report_miss` **warns** when every candidate was blind,
  since that is a defect in the axial path rather than an absent conformer.
- `_verify_axial_honored` runs at the single exit of `_select_by_geometry` and warns when the
  returned structure is not the requested atropisomer — the one place that turns a silent
  wrong answer into an observable event (telemetry `adapter.axial_not_honored`).
- `mol_axial_token` logs when every perception strategy fails, so `None` is traceable rather
  than mute. `None` means *could not tell*, never *no axial token*.

- **The token broke the generation parser.** `generation/oin_parser.py` decides inline-vs-sidecar by
  the *absence* of `|`, so a trailing ` |ax:-|` misrouted inline OINs to the sidecar branch and lost
  the geometry code (`Geometry code '' not supported`). The token is now stripped before parsing and
  preserved on `original_oin` for the adapter.
- **Key-folding defeated acceptance.** The SL1 early-exit gate (default ON) accepts the first
  conformer matching the fac/mer key — and the key deliberately *folds* the axial token so the batch
  harness stays neutral. That made the key an unsound *acceptance* test: it stopped the pool at
  either atropisomer before selection could choose. Acceptance is now key-match **and** axial-match.
  **Generalisable lesson: a deliberately-lossy comparison key must never be reused as an acceptance
  predicate for an axis it folds.**
- **Perception failed silently on pool conformers.** Raw MetalloGen pool mols are unsanitized, so
  `CanonicalRankAtoms` raises; the `try/except` turned that into `None`, and the filter compared
  `None` against `-` and matched nothing. Axial perception now runs on the **contract mol** the
  adapter already builds for re-encoding, with sanitize fallbacks.

### Population-at-risk — 2.49 %, a sound rate

Motif counting is conformation-independent, so unlike the Wave-1 rigid-mirror scan this is a real
rate rather than an upper bound (`tools/injectivity/axial_population.py`, seeded sample of 1500,
1488 scanned):

| population | count | fraction |
|---|---:|---:|
| has any inter-ring aromatic single bond | 345 | 23.19 % |
| hindered (twisted + ortho-walled) | 56 | 3.76 % |
| **emitting (hindered + stereogenic)** | **37** | **2.49 %** |
| hindered but NOT stereogenic (gate suppresses) | 19 | 1.28 % |

Two things follow. First, the **blast radius is 2.49 %** — the structures whose OIN changes under
the flag, and equally those exposed to the P2 blind spot today. That is the same order as the P3
metal-bound-amine motif (2.85 %), so P2 is not a curiosity confined to BINAP. Second, the
stereogenicity gate suppresses **19 of 56 hindered axes (34 %)** as achiral; without it the encoder
would be over-sensitive on 1.28 % of the corpus — a defect of the same magnitude as the blind spot
it was meant to fix.

One observation worth following up rather than over-reading: among the emitting examples the `+`
sign dominates (roughly 3:1). The benign explanation is chemical — catalysis datasets are rich in
commercially dominant single-enantiomer ligands such as (R)-BINAP — but a systematic skew in the
sign convention has not been excluded.

### The sign convention had to be audited first — and it was broken

Before trusting any rate, the convention itself was audited: mirror **every** emitting structure
and require the token to flip. A sound convention flips every time regardless of any skew.

**34 of 37 flipped; 3 did not** — OJELAQ, YESKOZ, EBUHAN, all multi-axis, and all chiral by the
independent geometric oracle (mirror RMSD 3.2–7.1 Å). So the failures were real, not meso.

The cause was self-inflicted. Canonicalizing the token, ties between symmetry-equivalent axes were
broken with the sign itself:

```python
axes.sort(key=lambda ax: (sym_key, ax.sign))   # destroys the chirality
```

Equivalent axes tie on symmetry rank, so they sorted **by sign**, forcing signs into ascending
order — a molecule carrying `+-` then rendered identically to its mirror carrying `-+`. The token
silently stopped being a chirality descriptor for exactly the structures that need one, and the
single-axis BINAP fixture could never have caught it. Ties are now broken with the tie-broken
canonical rank (also graph-derived, so still renumbering-invariant, but it keeps each sign attached
to its own axis). **After the fix: 37/37 flip — convention SOUND.**

### Cohort A/B — the rate

With a sound convention, the round-trip A/B was run over a 12-structure slice of the emitting
cohort (24 enantiomers, both arms, seed 42, FF, 120 s cap):

| cohort | axial-aware selection | baseline (no axial pass) |
|---|---:|---:|
| **single-axis** | **22 / 22 (100 %)** | 8 / 22 (36.4 %) |
| multi-axis (EBUHAN, 2 axes) | 0 / 2 | 0 / 2 |
| overall | 22 / 24 (91.7 %) | 8 / 24 (33.3 %) |

Fourteen structures are fixed by the axial-aware pass. For single-axis atropisomers the axis now
round-trips **without exception**, against a baseline barely better than chance.

**The multi-axis gap is a generator limit, not a sign error.** The generated multi-axis structures
return an *empty* token — the biaryl torsions relax outside the hindered window, so no axis is
detected at all. The generator can hold one atropisomeric axis but not two, and the axial pass
cannot select a conformer the pool never contains.

### Status and recommendation

Capability and evidence are both in place for the single-axis case. **Recommendation: gate
`OIN_EMIT_AXIAL` ON by default**, accepting one known cost — multi-axis structures (2 of the 37
emitting, ~0.13 % of the corpus) will emit a token the generator cannot reproduce, converting a
silent false positive into a loud, correctly-attributed false negative. That trade is the one this
whole program argues for, and the failure is now instrumented (`adapter.axial_not_honored`) rather
than silent. Flipping the default changes emitted output, so it is left as an explicit product
call rather than taken unilaterally.

Follow-on, tracked separately: teach the generator to hold multiple hindered axes (a pool-diversity
or torsion-constraint problem, not a notation one).

## Roadmap into the rest of Wave 2 / Wave 3

- **P2 opt-in emit — DONE** (`swimlane/y1-w2-axial`; see "Decision taken" and the close-out above).
- **To gate P2 ON by default:** a broader A/B over an atropisomer-bearing cohort (the single fixture
  pair is directional only). Capability is in place; this is now an evidence question.
- **P1 / P3:** fold into the v0.4.5 canonical-string work (canonical donor / neighbour ordering);
  the config oracle is the verification instrument for any token emitted there.
- **Configurational oracle → Wave 3:** these three descriptors are the seed of the
  configuration-aware oracle the plan defers to Wave 3 for the dataset-scale UU collision hunt
  (they replace the conformation-inflated rigid RMSD for motif-bearing axes).
