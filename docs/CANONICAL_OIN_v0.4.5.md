# Canonical OIN-SMILES — v0.4.5

**Status:** released on local `main` (`0d165845`, tag `v0.4.5`). Not pushed.

## What changed

A canonical organic SMILES is byte-identical no matter which 3D orientation or conformer you
start from. Before v0.4.5 OIN-SMILES had not earned that: the *comparison key*
(`oin/compare.py`) was canonical, the *emitted string* was not. On the 6,719-molecule capstone
sweep, **12.32% of molecules round-tripped to the same isomer but a different string**.

v0.4.5 moves the machinery `compare.py` already had upstream into the encoder, and turns it on.

### The six levers now ON by default

Registered in `src/oinsmiles/oin/levers.py::_DEFAULT_ON`:

| lever | what it canonicalizes |
|---|---|
| `OIN_CANONICAL_BODY` | ligand body — reparse through RDKit, carrying donor identity by atom map number |
| `OIN_CANONICAL_PERCEPTION` | the half a reparse cannot fix: which resonance form / valence walk is chosen |
| `OIN_CANONICAL_SLOTS` | slot labels, by lex-min colored-vertex signature over the proper-rotation group |
| `OIN_CANONICAL_ETA_WINDING` | haptic winding heading atom and `>`/`<` sign |
| `OIN_STABLE_METAL_AC` | valence-capping order (highest-Z first) so perception stops depending on atom order |
| `OIN_STABLE_STEREO` | tetrahedral tags re-derived from parent geometry rather than fragment order |

**Why these six were safe to promote together, and the rule to keep:** each one *repairs a
renumbered presentation without rewriting the canonical answer*. That is why the corpus shows no
churn. Levers that **add information** to the string (`OIN_EMIT_AXIAL`,
`OIN_EMIT_LOCKED_DONOR`) are a different trade — the generator must then be able to reproduce
what they emit, so promoting one converts a silent false positive into a loud false negative.
Those stay opt-in, with reasons in `levers.py::_HELD_OFF`.

### `levers.py` also closes a live trap

Two spellings of "read a lever" were in use, and they disagree:

```python
os.environ.get("OIN_EMIT_AXIAL")              # truthy -> "0" ENABLES it
os.environ.get("OIN_EARLY_EXIT", "1") != "0"  # "0" disables
```

The first is a trap: `"0"` is a non-empty string, so opting out the obvious way switched the
lever **on**. Everything now routes through `lever_enabled()`, where `0/false/no/off/""`
disable. `OIN_BORON_CAGE` alone had five sites on the old spelling.

## Measured

* byte-stability under rotation/renumbering **58.1% → 69.6%** (300-molecule seed-42 sample)
* comparison-key instability **60 → 16 molecules** (1-in-5 to roughly 1-in-19)
* re-baseline over 936 molecules: **145 of 436** previously-failing molecules **FIXED (33.3%)**;
  of 500 previously-passing guards, 11 apparent regressions are **all** `TimeoutException
  exceeded 300s` against a capstone baseline that ran at 1800 s → **zero correctness
  regressions**
* vetoes green: fac/mer and cis/trans still distinct raw *and* at key level; mirror guard green;
  `geometry_tag_shift` 0/298 `[M_XXX]` changes; goldens byte-identical on the opt-out path
* suite: **837 tests OK** (3 skipped, 4 expected failures)

## Goldens that moved

Each was checked to be a **relabeling, not a new isomer**, before being re-pinned —
`canonical_roundtrip_key` identical is the assertion that would catch a canonicalization which
merged two isomers:

| fixture | change | key |
|---|---|---|
| CisPlatin | `[Cl]{0}.[Cl]{1}.N{2}.N{3}` → `N{0}.N{1}.[Cl]{2}.[Cl]{3}` | identical |
| TransPlatin | `[Cl]{0}.N{1}.[Cl]{2}.N{3}` → `N{0}.[Cl]{1}.N{2}.[Cl]{3}` | identical |
| fac-Ir(ppy)₃ | `n{3}…n{1}…n{4}` → `n{5}…n{1}…n{3}` | identical |
| PdCl₂-R-BINAP | P moves to slots 2,3 | identical |
| PdCl₂-RR-BDPP / BDNN | `@@/@` → `@/@@` | **changed — correctly** |

`OIN_CANONICAL_SLOTS` labels by lex-min vertex colour and `"N" < "[Cl]"` bytewise, so the amines
take the low slots. Cis is still cis: the chlorides land on 2,3, which are adjacent.

### The BDPP/BDNN case is worth reading

Those two goldens, and the "RDKit CIP oracle" that blessed them, were **wrong for four months**.
The test asserted `(S,S)`; the fixtures are named `(2R,4R)`; the geometry says `(R,R)`.

The oracle ran `rdCIPLabeler` on a SMILES **reparsed from the encoder's own output**.
`rdCIPLabeler` converts a parity tag into an R/S label — it does not check the tag against
anything. Hand it an inverted tag and it returns an inverted label with full confidence. So the
"oracle" was a snapshot of the encoder, and an inverted tag was self-consistent and passed.

`AssignStereochemistryFrom3D` on the parent complex is the arbiter, because it is the one thing
no encoder bug can rewrite. Both tests now derive truth from coordinates *and* cross-check that
the emitted string agrees — the loop the old form left open.

This is the same failure shape as the Y2 near-miss (a canonicalization that made the axial token
reflection-invariant, caught only by a corpus-wide mirror audit): **a measurement that only
exercises the easy case will confirm a wrong belief.**

## Known gaps — deferred deliberately, not overlooked

1. **P3 (metal-bound secondary amine) is not usable in the shipped default.** It is built,
   oracle-validated and green in `tests/unit/test_locked_donor.py` — but only with
   `OIN_CANONICAL_BODY=0`. The canonical body's reparse sanitizes a *metal-free* fragment, where
   RDKit clears `[N@]` on a 2-degree amine as a freely inverting amine — the exact behaviour the
   descriptor exists to work around. The fix is to stamp the tag inside `canonical_body_emit`
   before its final `MolToSmiles`; that call also re-derives donor marker positions, and a
   misplaced marker silently mislabels coordination, so it waits for v0.4.6 and its own A/B.

2. **`OIN_EMIT_AXIAL`'s promotion evidence needs re-measuring.** `_is_atropisomer_candidate`
   gates its steric wall on `not GetIsAromatic()`, so `OIN_CANONICAL_PERCEPTION` changes the
   hindered-axis count: YESKOZ measures **2 → 1**. No emitted string moves today (YESKOZ's axes
   are non-stereogenic, BINAP is unchanged at 1 hindered / 1 emitting), but the Y2 cohort numbers
   (single-axis 22/22, mirror audit 37/37) were taken with perception OFF. `axial.py`'s safety
   argument covers the *generator* reading fewer aromatic atoms — not the encoder reading more.

3. **`OIN_H_FAITHFUL` is blocked** by the same reparse that blocks P3.

4. **Lane 5 (metal Δ/Λ, P1) was not started.** 0/150 molecules emit a metal `@` tag today, so it
   would be *creating* a descriptor rather than fixing a collapse — a larger piece of work than
   the lane budget assumed, and it needs the Δ/Λ tris-bidentate fixture Lane 7 built.

5. **`OIN_STABLE_METAL_AC` on degenerate input.** On a deliberately broken geometry
   (`ticat3_generated_broken.xyz`) the new capping order lets the metal absorb contested bonds
   and perception *succeeds*, returning nonsense: 8 fragments, seven bare `[H+]`, `[Ti-14]`,
   where the old order failed loudly. Real data is clean (145 fixed, `geometry_tag_shift`
   0/298), so this is a degenerate-input concern. A sanity gate rejecting stranded bare-proton
   fragments is the obvious follow-up and was **not** added, because charged hydrides are
   legitimate and the gate needs its own corpus A/B.

## The confound that governs every accuracy number here

**77.8% of round-trip failures never test the notation** (67.4% are 300 s timeouts). So a
round-trip pass rate is substantially *generator throughput*, not notation quality, and any
change that alters runtime moves the rate for unrelated reasons.

That is why the canonicality verdict in this release comes from a **generator-free** instrument
(`tools/canonicality_probe.py` — hold the graph fixed, vary proper rotation and atom numbering)
and not from the round-trip sweep. It is also why all 11 apparent regressions above are
timeouts: the re-baseline ran at 300 s against a 1800 s baseline. Comparing pass rates across
different timeout budgets is the same config asymmetry that manufactured v0.4.4's 11 phantom
regressions.

`tools/v045_state.sh` prints live state computed from the repo rather than from a doc.
