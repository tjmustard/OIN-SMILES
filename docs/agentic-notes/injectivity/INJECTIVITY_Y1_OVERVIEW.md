# Injectivity Audit (Y1) — Blind Spots of the OIN Round-Trip Test

**Status:** measurement / falsification wave — no encoder code changed. Reproducer:
`tools/injectivity/`. Regenerable output: `results-injectivity-y1/` (gitignored).

## TL;DR

The round-trip test (`XYZ → OIN → 3D → OIN`, compared by a canonical key) is the *only*
current check that an OIN is correct, and it is **necessary but not sufficient**. It proves
the composition `E∘G` is a *retraction* — it can **never** prove the encoder `E` is
*injective* (that two distinct isomers get two distinct strings). We built an independent,
generator-free harness that constructs a mirror-image twin of a structure and asks whether
the encoder can still tell the two apart. On curated rigid fixtures it already exhibits
**two concrete collisions**, so the hypothesis *"a passing round-trip proves the OIN is
lossless"* is **REFUTED**.

## The Confusion-Matrix Frame

Judge every OIN outcome against ground-truth isomer identity `O(s)`, supplied by an oracle
that never touches the OIN string:

| | Round-trip PASS | Round-trip FAIL |
|---|---|---|
| **OIN truly correct** | True Positive | **False Negative** — *missed success* (mostly the generator, not the OIN) |
| **OIN wrong / lossy** | **False Positive** — *the dangerous cell* (KU + UU) | True Negative |

Round-trip cannot see either off-diagonal cell by construction. The dangerous cell (a
*passing* round-trip hiding a lossy encoding) is this wave's target.

**On the False-Negative cell:** a round-trip FAIL with a *correct* OIN is usually the
**generator** — MetalloGen not yet 100% accurate, or the job **timing out** — not the
notation being too strict. The Wave-3 missed-success audit must attribute FAILs into
`generator_inaccuracy` / `generator_timeout` (not OIN defects) vs `genuinely_different_isomer`
vs `canonicalization_noise`. This wave deliberately sidesteps that noise: **the probes never
invoke the generator**, so they isolate a pure property of `E`.

## Method — the Distinguishability (Collision) Probe

The encoder should be a *bijection* between isomer-classes and OIN strings:

- **Conformer-invariance** (already established — encoder PASSES): all conformers of one
  isomer → one string.
- **Injectivity** (this wave): distinct isomers → distinct strings.

To prove `E` is blind to an axis `d`, no generator is needed: mirror a structure (a
z-reflection is an enantiomer-generating improper isometry), then compare
`convert(base)` vs `convert(mirror)`:

- `raw_equal` — the raw OIN strings are byte-identical.
- `key_equal` — the round-trip equivalence key agrees (**what the batch harness gates on**).
- `oracle_distinct` — the independent oracle says the mirror is a genuinely different isomer.

A **collision** is `oracle_distinct ∧ key_equal`. If `raw_equal` too, the encoder is *totally*
blind. Every probe drives divergence through `convert()` on **real twin geometries**, fixing a
flaw in the pre-existing `tests/integration/test_isomer_divergence.py::test_metal_stereo_raw_only`,
which hand-feeds synthetic `@SP1`/`@SP2` strings the encoder has **no producer** for
(`oin/inline.py:110`) and so gives false confidence that metal chirality survives in the raw
string.

### The independent oracle

Primary signal is a **geometric enantiomer test** (`tools/injectivity/oracle.py`): a rigid
molecule is chiral iff its mirror cannot be superimposed on it by a *proper* rotation modulo a
graph automorphism. It uses only geometry + topology, so it certifies distinctness for metal
Δ/Λ **and** biaryl atropisomers — the very axes where RDKit / InChI descriptor perception is
itself blind (an incidental finding: `FindPotentialStereo` returns *nothing* on the BINAP
ligand and on fac-Ir(ppy)₃). Caveat: rigid superposition conflates conformation with
configuration, so the oracle is valid only for **rigid** species; the curated fixtures are
rigid by construction.

## Results — H0 REFUTED

| axis (probe) | fixture | oracle-distinct | raw_equal | key_equal | verdict |
|---|---|---|---|---|---|
| achiral control | CisPlatin | False (0.05 Å) | True | True | ⚪ invariant ok |
| **metal Δ/Λ (P1)** | fac-Ir(ppy)₃ | True (3.19 Å) | False | **True** | 🟠 **key-blind — batch false-positive** |
| **axial / atropisomer (P2)** | PdCl2-R-BINAP | True (4.02 Å) | **True** | **True** | 🔴 **encoder-blind (total)** |
| **metal-bound 2° amine (P3)** | POJJOP | True (2.35 Å) | **True** | **True** | 🔴 **encoder-blind (total)** |

- **P1 metal Δ/Λ — KEY-BLIND (confirmed):** fac-Ir(ppy)₃ and its enantiomer produce
  *different* raw strings, but the difference is only non-reproducible slot renumbering, which
  the round-trip key **deliberately folds** (`oin/compare.py` `_METAL_STEREO_RE` +
  `_polyhedron_signature`). The batch harness gates on the key, so two genuine Δ/Λ enantiomers
  round-trip as identical. This is the documented, deferred metal-stereo limitation — now
  demonstrated through `convert()` on real geometry rather than hand-fed strings.
- **P2 axial / atropisomer — ENCODER-BLIND (confirmed):** R-BINAP and S-BINAP encode to
  **byte-identical** OIN strings. The encoder has no single-bond atropisomeric axis
  perception at all; the skipped guard `tests/unit/test_axial_chiral.py` now has a paired
  enantiomer proof.
- **P3 metal-bound 2° amine — ENCODER-BLIND (confirmed):** the user's motivating case, pinned to
  a real dataset example (POJJOP), whose *sole* stereocentre is a metal-bound secondary amine.
  Its two enantiomers encode byte-identically — the metal-locked N configuration is dropped by
  the Zone-A `total_degree<4` clear (`core/chirality.py:722-727`). NB the documented
  `JUCCUH [N@@H]→[NH]` case is *not* this — there the encoder DOES capture a pendant amine's
  stereo; that loss is the generator (a False Negative), not an injectivity failure.

Per-axis detail: `docs/agentic-notes/injectivity/INJECTIVITY_Y1_P1_METAL.md`, `_P2_AXIAL.md`, `_P3_AMINE.md`.

## Population-at-risk — not yet a rate

A dataset population scan exists (`report.py --population N`) but its rigid oracle **conflates
conformation with configuration at dataset scale**, so its fractions are conformation-inflated
**upper bounds, not losslessness-failure rates** (a 30-structure sample reads ~97% "chiral",
which is the rigid test failing on flexible crystal conformers, not a real chirality rate). A
trustworthy per-axis population needs the configurational oracle deferred to **Wave 3**. The
curated rigid fixtures above are the sound Y1 result.

Where an axis has a **structural motif**, though, the population *is* soundly countable because
motif-matching is conformation-independent. For **P3**, a pure-geometry pre-filter (N within
2.6 Å of a metal, with 1 H + 2 C neighbours) finds a **metal-bound secondary amine in ~2.85 %**
of the corpus (171 / 6000 sampled) — a real blast-radius floor for that blind spot. The same
motif approach is how Wave 3 should quantify each named axis.

## Reproduce

```
# from the worktree root, main venv (no uv sync; rdkit pinned):
PYTHONPATH=$PWD/src python -m tools.injectivity.report --probes
PYTHONPATH=$PWD/src python -m tools.injectivity.twin_collision tests/fixtures/PdCl2-R-BINAP.xyz
```

Outputs `results-injectivity-y1/{injectivity_metrics.json,report.md}` (byte-stable; no
timestamps). Guard tests: `tests/unit/test_injectivity_probes.py`.

## Roadmap — all three waves complete

- **Wave 1 (this doc):** oracle + harness + the three named known-unknown probes. H0 REFUTED;
  P1/P2/P3 confirmed. ✅
- **Wave 2** (`INJECTIVITY_Y2_FEASIBILITY.md`): all three axes proved **recoverable from the
  input 3D** — the encoder discards a signal it has, so none is a permanent limitation. P2 was
  then closed out end-to-end: a canonical token, a stereogenicity gate (which fixed a real
  *over-sensitivity* defect), and generator support, so the axis now **round-trips**
  (2/2 vs 1/2 baseline). P1/P3 relabelled recoverable-but-deferred to v0.4.5. ✅
- **Wave 3** (`INJECTIVITY_Y3_UNKNOWN_UNKNOWNS.md`): the missed-success audit found that
  **77.8 % of round-trip failures never test the notation** (two thirds are timeouts), so the
  pass-rate largely measures generator throughput. The UU hunt, using InChI as a
  conformation-independent discriminator, found **no new blind spot in 299 structures**, while
  21.7 % of collapses trace to the already-named P1/P2/P3 axes. Linkage isomerism refuted as a
  blind spot; donor swap undetermined. ✅

### What the program established

The thesis holds: **the round-trip test cannot see either off-diagonal cell of its own
confusion matrix.** Both now have instruments that do not pass through the artefact under
test — generator-free collision probes for the false positives, cause attribution for the
false negatives. The practical corollary is that *pass-rate is not losslessness*: they are
different quantities, measured by different tools, and moving one does not imply moving the
other.

Confirmed blind spots are catalogued in `docs/KNOWN_LIMITATIONS.md`.
