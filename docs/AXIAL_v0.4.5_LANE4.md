# Axial atropisomerism, v0.4.5 Lane 4 — the multi-axis failure was never the generator

Lane 4 was chartered to "teach the 3D generator to hold two or more hindered biaryl axes", on the
basis that multi-axis structures round-tripped **0/2 in both A/B arms** because "the biaryl torsions
relax outside the hindered window (20–160°), so no axis is detected at all". Two prescribed
remedies followed from that reading: widen the embed pool, or guard the FF relaxation.

Instrumenting the pool first — the lane's own mandatory step 1 — refuted the premise. Neither
remedy would have helped, because neither addresses what is actually wrong.

---

## 1. What the pool histogram showed

`tools/injectivity/axial_pool_histogram.py` dumps, for every conformer the adapter actually sees,
the full per-axis detail behind its token: the signed dihedral, whether it passed the hindered
window, and whether it is stereogenic. On `YESKOZ` (the lane's primary multi-axis target, a
5,15-bis(2-(methylthio)phenyl) Zn porphyrin with an axial pyridine):

| rank | token | axes |
|---:|---|---|
| 0 | `""` | **(none detected)** |
| 1 | `""` | **(none detected)** |
| 2 | `""` | **(none detected)** |
| 3 | `""` | **(none detected)** |
| 4 | `""` | **(none detected)** |

"None detected" is not "relaxed flat". `detect_axial_axes` returns an axis for *every* qualifying
inter-ring single bond regardless of twist, and marks it `hindered=False` if the torsion is
near-planar. An empty list means no bond qualified as an axis **at all**.

Measuring the same bonds directly, ignoring the axis-selection criterion, shows the twists are
present and comfortably hindered:

| pool rank | axis 1 dihedral | axis 2 dihedral |
|---:|---:|---:|
| 0 | +87.7° | +122.1° |
| 1 | +128.0° | −14.7° |
| 2 | +157.2° | −124.5° |

**The generator holds both hindered axes, and the pool spans several sign combinations** — exactly
the raw material selection needs. Widening the pool would have added more of something already
present; guarding the relaxation would have protected a torsion that was never being flattened.

## 2. Why nothing could see them

The axis-selection test required **both** ends of the bond to be flagged aromatic. The encoder and
the generator do not agree about that flag on a metalloporphyrin:

| | macrocycle as perceived | aromatic atoms |
|---|---|---:|
| encoder — `get_tmc_mol`, bond orders from interatomic distances | aromatic pyrrolide core on Zn(II), two `[n-]` | 38 |
| generator — `build_contract_mol`, bond orders transferred from the OIN fragment SMILES | neutral localized tautomer, four dative N, written M→N | 18 |

A porphyrin *meso* carbon is therefore aromatic for the encoder and aliphatic for the generator, so
the meso-aryl axis exists on the input side and vanishes on the generated side. That is the entire
mechanism behind "multi-axis fails": the corpus's multi-axis structures are porphyrins, so
"multi-axis" was a **confound** for "macrocycle whose aromaticity perception is route-dependent".

Re-perceiving the generated coordinates through the encoder's own `get_tmc_mol` does not rescue it
(measured: still no axes), and neither would fixing the contract-mol transfer, because the OIN
fragment SMILES the encoder emitted does not itself survive RDKit sanitisation as aromatic:

```
CSc1ccccc1-c1c2nc(cc3ccc(n3)c(...)c3nc(cc4ccc1n4)C=C3)C=C2
  -> SanitizeMol -> CSc1ccccc1C1=C2C=CC(=N2)C=C2C=CC(=N2)... (12/40 aromatic)
```

So the descriptor had to stop depending on perception. It could not be worked around downstream.

## 3. The fix

`src/oinsmiles/oin/axial.py` now derives everything that selects an axis and fixes its sign from
properties both routes agree on:

| was | now | why |
|---|---|---|
| axis end is `GetIsAromatic()` | axis end is a **trigonal ring atom** — in a ring, 3 heavy neighbours, no H | a ring atom with 4 valences and 3 sigma bonds must hold a pi bond, whatever model perceived it. Strict superset: an aromatic end always has two ring neighbours plus the partner and no H, and an aromatic atom already holding three ring bonds has no valence left to bear an axis |
| reference neighbours are aromatic neighbours | reference neighbours are **ring** neighbours | identical set for an aromatic end (an axis end is never a fusion atom), but survives a localized ring |
| `CanonicalRankAtoms` on the mol as perceived | ranks on a **connectivity skeleton** — bond orders, charges, aromatic flags erased, metal bonds dropped | the reference neighbour sets the SIGN; ranks that differ between routes compare two differently-defined quantities and can report a match for the mirror image |
| stereogenicity + reference chosen on the intact graph | chosen with the **axis bond cut** | strictly finer, and removes a silent coin toss: on a tie `max()` took whichever neighbour the list yielded first, and the two candidates sit 180° apart |

Metal bonds are *dropped* from the skeleton rather than down-graded to single bonds: down-grading
closes every chelate ring, which buries BINAP's own axis inside the P–M–P ring and makes
`IsInRing()` true. It also sidesteps a second route disagreement — `DATIVE` direction is
begin-to-end and the two routes write it opposite ways round.

**Guard.** `tools/injectivity/axial_perception_sweep.py` applies a worst-case perception
perturbation (every non-metal bond flattened, aromaticity and charges cleared) to the same
coordinates and asserts the token is unchanged. Coordinates, elements and connectivity are
untouched, so handedness is untouched; only the perception is. Unit-level equivalents are in
`tests/unit/test_axial_emit.py::TestPerceptionInvariance`.

## 4. The porphyrin signs are COUPLED — which is why two lanes' results looked contradictory

With ranks taken from the connectivity rather than from an arbitrary resonance form, `YESKOZ`'s two
meso-aryl axes are reported **non-stereogenic**, and it emits nothing. Lane 7 meanwhile showed that
a `YESKOZ` **single-axis flip** *is* separated by the raw string with `OIN_EMIT_AXIAL=1`, and is
collapsed only by the key's `_AXIAL_TOKEN_RE` fold — evidence that the token carries real
information about exactly these structures.

Both results are correct, and measuring the raw signs (stereogenic gate ignored) on four variants
of `YESKOZ` shows precisely how they fit together:

| variant | raw signs | canonical up to global inversion |
|---|---|---|
| deposited (**anti**, α/β) | `--` | `++` |
| z-mirror of deposited | `++` | `++` |
| **single-axis flip** (**syn**, α/α) | `-+` | `+-` |
| z-mirror of the flip | `+-` | `+-` |

Two things fall out:

1. **Each configuration is reflection-invariant** — the mirror produces the exact complement, which
   canonicalizes to the same string. So `YESKOZ` is **achiral** about these axes, and the per-axis
   sign carries no handedness. This matches the symmetry argument: a 5,15-diarylporphyrin's meso
   carbon is flanked by two equivalent pyrroles, syn is fixed by the mirror plane through the *other*
   two meso positions, and anti by that plane composed with the C2 about the porphyrin normal.
2. **The pair still separates the diastereomers** — `++` (anti) ≠ `+-` (syn). The relative sign is
   real; only the absolute per-axis sign is not.

The porphyrin's own symmetry inverts **both** signs at once, so the sign vector is well-defined only
up to *global inversion*. That single fact explains every prior observation:

- the old descriptor's `|ax:+-|` **claimed chirality that is not there** — it survived the corpus
  mirror audit only because the mirror produces the complement, which the audit scores as "flipped
  correctly". The audit was confirming a false positive, not catching one;
- the absolute signs came from whichever arbitrary resonance form `AC2BO` returned
  (`xyz2mol_local.py:800` says so in as many words), so they were never reproducible;
- multi-axis structures scored 0/2 in *both* A/B arms. A generator limitation would have shown a
  difference between the arms. A request that is not well-posed shows none, which is what was
  observed.

**The fix this specifies (v0.4.6).** Do not drop coupled axes and do not emit their raw signs.
Canonicalize the sign vector over the orbit of the automorphism group's action on it — for a
coupled pair that is `min(token, complement)`, which is simultaneously reflection-invariant (correct:
achiral) and syn/anti-separating (correct: real diastereomers). The gate must distinguish *coupled*
ambiguity (one automorphism inverts every sign — quotient by global inversion, information survives)
from *independent* ambiguity (each axis has its own local C2 — quotient by independent flips,
nothing survives). Per-axis rank comparison cannot express that distinction; it needs the
automorphism group's action on the sign vector, enumerated with a blow-up guard and a conservative
fallback. Lane 7's `invert_stereocenter` operator is the right instrument to validate it, since
mirroring cannot isolate one axis of a multi-axis molecule.

**Until then the gate stays conservative** and emits nothing rather than a sign it cannot reproduce.
The cost is explicit and should be weighed by whoever schedules v0.4.6: the OIN currently encodes
syn and anti meso-arylporphyrins identically. That is a *smaller* loss than the status quo, which
emitted a non-reproducible sign and asserted a chirality that does not exist — but it is a real one.
Tracked in `docs/KNOWN_LIMITATIONS.md`.

## 5. Why the prescribed remedies were not applied

| prescribed remedy | why not |
|---|---|
| widen the pool when the token has length ≥2 (eta wide-pool precedent) | the pool already contained both hindered axes in several sign combinations; there was nothing to sample harder for |
| guard the FF relaxation so it cannot flatten a required torsion | the torsions were not being flattened (+87.7°/+122.1° at rank 0) |
| constrain the torsion at embed (`SetDihedralDeg` + constrained minimize) | this is *construction*, which the project has three negative results against — and it would have constructed a twist to satisfy a token that should never have been emitted |

Both cheap remedies target conformer *sampling*. The failure was in conformer *perception*, one
layer up, and it also disabled the instrument that was supposed to detect the failure —
`_verify_axial_honored` compares against the same blind `mol_axial_token`, so it reported nothing.

## 6. Reproduce

```bash
cd /home/tjmustard/Documents/GitHub/oin-v045-lane4
export PYTHONPATH=$PWD/src
V=/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python

# the diagnosis: what does the embed pool actually contain?
$V -m tools.injectivity.axial_pool_histogram tests/fixtures/YESKOZ.xyz

# perception invariance (the new guard), fixtures and/or corpus
$V -m tools.injectivity.axial_perception_sweep --fixtures
$V -m tools.injectivity.axial_perception_sweep --dataset default --n 400

# corpus population + sign-convention audit; --jobs is deterministic, --tag avoids
# overwriting a baseline scan
$V -m tools.injectivity.axial_population --n 1500 --mirror-check --jobs 4 --tag skeleton

# guards
$V -m unittest tests.unit.test_axial_emit tests.unit.test_axial_failure_modes
```
