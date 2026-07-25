# Injectivity Y1 · P3 — Metal-bound secondary-amine chirality (ENCODER-BLIND)

**Status:** measurement only. Verdict: **CONFIRMED encoder-blind (total)**. Parent:
`INJECTIVITY_Y1_OVERVIEW.md`. This is the user's motivating case ("a chiral amine, bound to a
metal, flips R→S but the system doesn't capture it") — pinned to a real dataset example.

## First, a correction to the hypothesis

The obvious candidate — the documented `JUCCUH` `[N@@H]→[NH]` residual — is **not** an encoder
blind spot. `convert(JUCCUH)` emits `...C[N@@H](C)Cc2ccc...` and its mirror emits `[N@H]`, so
the encoder **distinguishes** them (raw + key diverge). JUCCUH's stereogenic N is a *pendant*
backbone amine; the `[N@@H]→[NH]` loss happens on the **generated re-encode** (the generator
can't rebuild it), i.e. it is a False-Negative / generator effect, not an injectivity failure.

The true P3 blind spot is narrower: a **metal-bound** secondary amine, stereogenic *only*
because the metal locks its fourth position.

## The isomer pair (POJJOP)

`tests/fixtures/POJJOP.xyz` (CSD POJJOP, from the tmCAT set) is a square-planar Pd complex whose
**sole stereocentre is a metal-bound secondary amine** N — bonded to Pd, H, a tolyl, and a
CH₂-pyridyl (four distinct groups when the Pd bond is present); the rest of the complex (tolyl,
pyridyl, two Cl) is achiral. Its z-mirror inverts that N. Geometric oracle: distinct, mirror
RMSD **2.35 Å**.

## Measurement

`convert(base)` and `convert(mirror)` are **byte-identical**:

```
[Pd_SPL].Cc1ccc([NH]{0}Cc2ccccn{1}2)cc1.[Cl]{2}.[Cl]{3}
```

`raw_equal = True`, `key_equal = True` → **ENCODER-BLIND (total)**. The amine appears as bare
`[NH]{0}` with no `@` — the metal-locked configuration is dropped.

## Mechanism

Once the encoder strips the metal to fragment the complex, the amine N is trivalent
(`total_degree < 4`) and its chiral tag is cleared by the Zone-A rule
(`core/chirality.py:722-727`); nitrogen is also explicitly out of lone-pair-CIP scope
(`:33-36`, "trivalent [N@] is cleared by RDKit as non-stereogenic amine inversion"). The
information that the metal bond made the centre stereogenic is gone.

## Population-at-risk

Unlike the rigid-mirror chirality scan (conformation-inflated), the P3 motif is
**conformation-independent** and countable directly: a pure-geometry pre-filter (N within 2.6 Å
of a transition metal, with exactly 1 H and 2 C neighbours) found **171 metal-bound
secondary-amine motifs in 6000 sampled structures (~2.85 %)**. Of a 12-structure probe sample,
the majority (POJJOP, ABIFAV, CULGOF, ATEFIP, RIFGUJ, MOKCEV, …) were encoder- or key-blind;
a minority (TEGFET, KANYUW, LARYEI) were *distinguished* because they carry an **additional**
backbone stereocentre the encoder does capture. So the ~2.85 % is the motif prevalence; the
blind subset is those where the bound amine is the *only* stereo-determining element.

## Verdict & disposition

**Recoverable, deferred to v0.4.5 (NOT permanent).** Wave 2
(`docs/INJECTIVITY_Y2_FEASIBILITY.md`) showed the locked-N configuration is recoverable straight
from the 3D as the signed tetrahedral volume of the metal-bound N's four neighbours (POJJOP:
−9.4 vs +9.4). RDKit will not perceive it (it clears the N as an invertible amine), so recovery is
ours. Emitting it needs a canonical neighbour ordering (the v0.4.5 problem), a Zone-A carve-out for
the metal-locked case (`core/chirality.py:722-727`), and a generator that can set it — hence
deferred to v0.4.5 rather than filed permanent. Guard:
`tests/unit/test_injectivity_probes.py::test_metal_bound_amine_is_encoder_blind` (+ aspirational
`test_metal_bound_amine_should_diverge_in_raw_string`). Reproduce:
`PYTHONPATH=$PWD/src python -m tools.injectivity.twin_collision tests/fixtures/POJJOP.xyz`.
