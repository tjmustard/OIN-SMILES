# Lane 4 — Axial atropisomerism (Y1 blind spot **P2**)

**The blind spot:** a hindered biaryl axis (BINAP-type) has two non-superimposable
configurations, and the encoder wrote the same OIN string for both — `PdCl2-(R)-BINAP` and its
mirror image encoded **byte-identically**, so the notation could not say which enantiomer it
described.

---

## ELI5

Some molecules have two flat ring systems joined by a single bond that is too crowded to spin
freely, so the molecule is frozen with one ring twisted clockwise or anticlockwise relative to
the other. Those two frozen forms are **mirror images**: like your left and right hands, they
have exactly the same parts joined in exactly the same order, but no amount of turning will let
one sit exactly on top of the other. Drug chemistry cares intensely about this — one hand can
be a catalyst and the other useless — so a notation that writes the same text for both is
lying. OIN did exactly that: it wrote one string for both hands. This lane's job was the
*second half* of the fix: the string that tells the two hands apart already existed (added in
an earlier wave, behind an off-by-default switch), and Lane 4 was supposed to teach the 3D
structure *builder* to rebuild the harder cases. It turned out the builder was already doing
its job correctly, and the real bug was that nothing downstream could **see** the twists it
had built.

---

## The work, visually

```
 THE PIPELINE, AND WHERE P2 LIVES
 ================================

   input .xyz ──► get_tmc_mol()  ──► CIPAssigner ──► OINDiscreteAligner ──► OIN string
   (3D coords)    perception #1                                              │
                  (bond orders                                              │
                   from distances)                        OIN_EMIT_AXIAL ───┤ (default OFF)
                                                          appends " |ax:+-|"│
                                                                            ▼
                                                     [Pd_SPL]....[Cl]{3} |ax:-|
                                                                            │
   generated .xyz ◄── _select_by_geometry ◄── embed pool ◄── m-SMILES ◄──────┘
        ▲                     │
        │                     └── build_contract_mol()  ── perception #2
        │                         (bond orders transferred from OIN fragment SMILES)
        │
        └── axial-aware narrow: keep conformers whose OWN token == requested token

 ✗ THE DEFECT LANE 4 FOUND: perception #1 and perception #2 DISAGREE
 ==================================================================

   Zn-porphyrin, YESKOZ            aromatic atoms   meso carbon reads
   ------------------------------  --------------   ------------------
   perception #1 (encoder)              38          AROMATIC   -> axis FOUND
   perception #2 (generator)            18          aliphatic  -> axis GONE

   old axis test:  GetIsAromatic(end1) AND GetIsAromatic(end2)
                   └──────────── not stable across the two routes ────────────┘

 THE REFUTATION SEQUENCE
 =======================

   charter premise            "generator can't hold 2 hindered axes;
                               torsions relax out of the 20-160 deg window"
          │
          ▼  measure the embed pool FIRST (axial_pool_histogram.py, YESKOZ)
   pool token, ranks 0..4     ""  ""  ""  ""  ""        <- no axis DETECTED
   raw dihedrals, rank 0      +87.7 deg / +122.1 deg    <- both twists ARE there
   raw dihedrals, rank 1      +128.0 / -14.7            <- pool spans sign combos
   raw dihedrals, rank 2      +157.2 / -124.5
          │
          ✗ PREMISE REFUTED. Not "relaxed flat" — "never selected as an axis".
          │
          ├─► prescribed remedy 1: widen the pool         -> nothing to sample harder for
          ├─► prescribed remedy 2: guard FF relaxation    -> nothing was being flattened
          └─► prescribed remedy 3: constrain the torsion  -> construction; 3 prior negatives
          │
          ▼  make the descriptor perception-INDEPENDENT (axial.py)
   aromatic flag            -> trigonal ring atom (in ring, 3 heavy nbrs, 0 H)
   aromatic neighbours      -> ring neighbours
   ranks on perceived mol   -> ranks on a connectivity SKELETON (metal bonds dropped)
   ranks on intact graph    -> ranks with the AXIS BOND CUT
          │
          ▼  consequence
   YESKOZ now emits NOTHING: its 2 meso-aryl axes are non-stereogenic, and the
   4-variant sign table proves both configurations are ACHIRAL (mirror = complement)

 THE NEAR-MISS: CANONICALIZATION THAT DESTROYED THE CHIRALITY
 ============================================================

   axes.sort(key=lambda ax: (sym_key, ax.sign))     <-- WRONG
                                       ^^^^^^^
   equivalent axes tie on sym_key  ->  tie broken BY SIGN  ->  signs forced ascending

        molecule  "+-"  ─render─►  "-+"
        mirror    "-+"  ─render─►  "-+"        IDENTICAL. reflection-INVARIANT.

   caught by:  corpus mirror audit, 34 / 37 flip   (OJELAQ, YESKOZ, EBUHAN failed)
   fixed by:   tie-break on TIE-BROKEN canonical rank  ->  37 / 37 flip
   ⚠ the single-axis BINAP fixture could NEVER have caught it (one axis = no tie)

 LEGEND
 ======
   ──►      data flow            ✗       refuted / defect
   ▼        sequence in time     ⚠       trap worth remembering
   {n}      OIN slot marker      |ax:±|  the opt-in axial sidecar token
   "deg"    degrees              perception #1/#2  the two independent bond-order routes
```

---

## Initial assumptions and hypothesis

Lane 4's charter (`spec/handoffs/v0.4.5/Lane4-axial-gate-multiaxis.md`) had two parts and both
rested on prior-wave conclusions:

**Part A — flip `OIN_EMIT_AXIAL` to default ON.** The handoff states: *"The evidence to flip it
exists; the flip does not."* It prescribed changing the default in `xyz2mol.py`'s
`_axial_suffix` block, keeping `OIN_EMIT_AXIAL=0` as the opt-out, mirroring how `OIN_EARLY_EXIT`
was promoted in v0.4.4.

**Part B — teach the generator to hold ≥2 hindered axes.** The handoff's stated mechanism:
*"The generated structure returns an **empty** token: the biaryl torsions relax outside the
hindered window (20–160°), so no axis is detected at all. The generator holds one hindered axis
but not two."* Prescribed remedies, cheapest first: widen the embed pool when the requested
token has length ≥ 2 (precedent: the eta-winding wide pool), or guard the FF relaxation so it
cannot flatten a required torsion.

Both parts inherited their evidence from the Y2 wave (`docs/agentic-notes/injectivity/INJECTIVITY_Y2_FEASIBILITY.md`),
whose own summary reads: *"The multi-axis gap is a generator limit, not a sign error."*

### What the encoder side already was, before Lane 4

Landed in the prior Y2 wave (branch stack `y1-audit` → `y1-w2` → `y1-w2-axial`), so Lane 4 did
**not** build any of this:

- a canonical ` |ax:+|` / ` |ax:-|` / ` |ax:+-|` trailing sidecar token behind `OIN_EMIT_AXIAL`,
  computed from the **signed biaryl dihedral** in `src/oinsmiles/oin/axial.py`;
- a **stereogenicity gate**: an axis whose ring end has symmetry-equivalent ortho neighbours has
  a local C2 through the axis, so it is achiral however twisted, and emits nothing;
- the round-trip key **folding** the token (`oin/compare.py::_AXIAL_TOKEN_RE`);
- generator support by **selection**: `_select_by_geometry`'s axial-aware narrow plus an
  axial-aware acceptance predicate.

---

## What was actually found

### Confirmed (and unchanged by Lane 4)

| claim | measurement | source |
|---|---|---|
| P2 is total encoder blindness by default | `convert(R-BINAP)` and `convert(mirror)` byte-identical; `raw_equal = True`, `key_equal = True` | `docs/agentic-notes/injectivity/INJECTIVITY_Y1_P2_AXIAL.md` |
| the two atropisomers really are distinct | min proper-rotation mirror RMSD **4.02 Å** | same |
| RDKit does not perceive the axis from pure 3D | `FindPotentialStereo` returns 0 elements on the BINAP ligand; per the RDKit Book an atropisomer bond is marked only when a neighbour bond is *wedged* | `docs/agentic-notes/injectivity/INJECTIVITY_Y2_FEASIBILITY.md` |
| the configuration IS in the input coordinates | signed biaryl dihedral **−75° → +75°** base → mirror | same |
| single-axis round trip works with the token | cohort A/B **22/22 (100%)** vs baseline **8/22 (36.4%)** — 14 molecules fixed | `axial_cohort_ab.py --limit 12` |
| fixture-pair A/B | **2/2** with the axial pass vs **1/2** without (the baseline returns the same handedness whatever is requested — its single match is luck) | `axial_roundtrip_ab.py`, seed 42, FF |
| the stereogenicity gate is doing real work | suppresses **19 of 56** hindered axes (34%) as achiral; the independent geometric oracle agrees exactly — mirror RMSD **0.000 Å** on symmetric-end cases vs **2.80 Å** on the true atropisomer | `docs/agentic-notes/injectivity/INJECTIVITY_Y2_FEASIBILITY.md` |
| ETKDG has no atropisomer bias | 21 / 19 split on the free ligand, so the right conformer is always in the pool | same |
| multi-axis scored 0/2 in **both** A/B arms | EBUHAN, 2 axes | same |

Population-at-risk, measured by motif counting (conformation-independent, so a real rate rather
than an upper bound) — `tools/injectivity/axial_population.py`, seeded sample of 1500, 1488
scanned:

| population | count | fraction |
|---|---:|---:|
| has any inter-ring aromatic single bond | 345 | 23.19% |
| hindered (twisted + ortho-walled) | 56 | 3.76% |
| **emitting (hindered + stereogenic)** | **37** | **2.49%** |
| hindered but NOT stereogenic (gate suppresses) | 19 | 1.28% |

⚠ **These numbers were measured against the PRE-Lane-4 descriptor.** See
*Open questions* — no post-fix re-measurement is recorded anywhere in the repo.

### Refuted

**R1 — "the generator cannot hold two hindered axes."** REFUTED by
`tools/injectivity/axial_pool_histogram.py` on `YESKOZ` (5,15-bis(2-(methylthio)phenyl) Zn
porphyrin with an axial pyridine), the lane's primary multi-axis target. Every conformer the
adapter actually sees reports an **empty** token — and "no axis detected" is not "relaxed flat",
because `detect_axial_axes` returns an axis for *every* qualifying inter-ring single bond
regardless of twist and only marks it `hindered=False` when near-planar. An empty list means no
bond qualified as an axis **at all**. Measuring the same bonds directly:

| pool rank | axis 1 dihedral | axis 2 dihedral |
|---:|---:|---:|
| 0 | +87.7° | +122.1° |
| 1 | +128.0° | −14.7° |
| 2 | +157.2° | −124.5° |

The generator holds both hindered axes and the pool spans several sign combinations — exactly
the raw material selection needs. **A generator limitation would have separated the two A/B
arms; a request that is not well-posed shows no difference between them, which is what 0/2 vs
0/2 was actually reporting.**

**R2 — "aromatic flags are a safe basis for axis selection."** REFUTED by measuring the two
perception routes on the same metalloporphyrin:

| route | macrocycle as perceived | aromatic atoms |
|---|---|---:|
| encoder — `get_tmc_mol`, bond orders from interatomic distances | aromatic pyrrolide core on Zn(II), two `[n-]` | **38** |
| generator — `build_contract_mol`, bond orders transferred from the OIN fragment SMILES | neutral localized tautomer, four dative N written M→N | **18** |

So a porphyrin *meso* carbon is aromatic for the encoder and aliphatic for the generator. The
old test required **both** axis ends to be flagged aromatic, so the meso-aryl axis existed on
the input side and vanished on the generated side. "Multi-axis" was a **confound** for
"macrocycle whose aromaticity perception is route-dependent" — the corpus's multi-axis
structures are porphyrins.

Two downstream workarounds were also measured and rejected: re-perceiving the generated
coordinates through the encoder's own `get_tmc_mol` still finds no axes, and fixing the
contract-mol transfer would not help either, because the OIN fragment SMILES the encoder emitted
does not itself survive RDKit sanitisation as aromatic — **12 of 40** atoms aromatic after
`SanitizeMol`.

**R3 — "the emitted `|ax:+-|` for a porphyrin was real chirality."** REFUTED by measuring raw
signs (stereogenic gate ignored) across four `YESKOZ` variants:

| variant | raw signs | canonical up to global inversion |
|---|---|---|
| deposited (**anti**, α/β) | `--` | `++` |
| z-mirror of deposited | `++` | `++` |
| **single-axis flip** (**syn**, α/α) | `-+` | `+-` |
| z-mirror of the flip | `+-` | `+-` |

Each configuration's mirror is its exact **complement**, so each is reflection-invariant:
`YESKOZ` is **achiral** about these axes and the per-axis sign carries no handedness. The
symmetry argument agrees — a 5,15-diarylporphyrin's meso carbon is flanked by two equivalent
pyrroles; *syn* is fixed by the mirror plane through the *other* two meso positions, *anti* by
that plane composed with the C2 about the porphyrin normal. **The porphyrin's own symmetry
inverts both signs at once, so the sign vector is well-defined only up to global inversion.**

That single fact explains three prior observations at once: the old `|ax:+-|` claimed chirality
that is not there; the absolute signs came from whichever arbitrary resonance form `AC2BO`
returned (`xyz2mol_local.py:800` says so in as many words), so they were never reproducible; and
the corpus mirror audit **passed vacuously** on these structures, because the mirror produces
the complement, which the audit scores as "flipped correctly". **The audit was confirming the
false positive it was trusted to catch.**

The `++` (anti) vs `+-` (syn) *pair* still separates real diastereomers, so the **relative**
sign is genuine information. Only the **absolute** per-axis sign is ill-defined.

---

## What was done

Lane 4 is **two commits**, merged into `release/v0.4.5` as `43e461e0`:

| commit | what |
|---|---|
| `10214f3e` | `fix(axial): make the atropisomer descriptor perception-independent` |
| `3144bac0` | `docs(axial): record the Lane 4 diagnosis, the key-folding decision and the v0.4.6 recommendation` |

### The descriptor change — `src/oinsmiles/oin/axial.py`

Everything that selects an axis or fixes its sign now derives from properties both perception
routes agree on:

| was | now | why |
|---|---|---|
| axis end is `GetIsAromatic()` | axis end is a **trigonal ring atom** — `_is_trigonal_ring_atom()`: in a ring, 3 heavy neighbours, no H, not a metal | a ring atom with four valences and three sigma bonds must hold a pi bond whatever model perceived it. **Strict superset**: an aromatic axis end always has two ring neighbours plus the partner and no H, and an aromatic atom already carrying three ring bonds (a fusion atom) has no valence left to bear an axis — so no axis that used to be found is lost |
| reference neighbours are aromatic neighbours | reference neighbours are **ring** neighbours — `_ring_neighbors()` | identical set for an aromatic end (an axis end is never a fusion atom), but survives a localized ring |
| `CanonicalRankAtoms` on the mol as perceived | ranks over a **connectivity skeleton** — `_skeleton()`: bond orders, charges and aromatic flags erased, metal bonds **dropped** | the reference neighbour sets the **sign**; ranks that differ between routes compare two differently-defined quantities and can report a match for the *mirror image* |
| stereogenicity + reference chosen on the intact graph | chosen with the **axis bond cut** — `_axis_cut_ranks()` | strictly finer (cutting can only distinguish atoms the intact graph merged), and it removes a silent coin toss: on a tie `max()` took whichever neighbour the list yielded first, and the two candidates sit ~180° apart |

**Metal bonds are dropped, not down-graded to single.** Down-grading closes every chelate ring,
which buries BINAP's own axis inside the P–M–P ring and makes `IsInRing()` true, deleting the
axis. Dropping also sidesteps a second route disagreement: `DATIVE` direction is
begin-to-end and the two routes write it opposite ways round (encoder N→M,
`build_contract_mol` M→N).

Also in `10214f3e`: `OIN_EMIT_AXIAL=0` used to turn the lever **ON**, because the gate was a
bare truthiness test and `"0"` is a non-empty string. It now reads `"0"`-means-off, mirroring
`OIN_EARLY_EXIT`. Default stays OFF. (This was later centralized into
`src/oinsmiles/oin/levers.py::lever_enabled`.)

### Part A was deliberately NOT done

The confirmed v0.4.5 product call kept **all** injectivity levers default-OFF, because they ADD
information to the string and the generator must then be able to reproduce what they emit —
promoting one converts a silent **false positive** into a loud **false negative**.

⚠ **Flag the discrepancy explicitly:** this is a deliberate deviation from Lane 4's own handoff,
whose Part A said *"Change `OIN_EMIT_AXIAL` to default ON"* in step 1. The handoff was written
assuming the flip; the release did not take it. `spec/handoffs/v0.4.5/Lane4-axial-gate-multiaxis.md`
was **not** updated, so anyone reading the handoff will find an instruction the release
deliberately declined.

Consequences of holding the default OFF, all deliberate:

- **Three `@expectedFailure` tests stay xfail on the default path** in
  `tests/unit/test_injectivity_probes.py` — `test_metal_chirality_should_diverge_at_key`,
  `test_axial_should_diverge_in_raw_string`, `test_metal_bound_amine_should_diverge_in_raw_string`
  (verified: the module runs 11 tests, `OK (expected failures=3)`).
- An accidental promotion is now itself guarded:
  `tests/unit/test_axial_emit.py::TestDefaultOff::test_emit_gate_is_off_unless_the_env_var_is_set`
  asserts that an *unset* environment does not emit and that an explicit `"0"` stays off.

### Why each prescribed remedy was rejected

| prescribed remedy | why not |
|---|---|
| widen the pool when the token has length ≥2 (eta wide-pool precedent) | the pool already contained both hindered axes in several sign combinations; there was nothing to sample harder for |
| guard the FF relaxation so it cannot flatten a required torsion | the torsions were not being flattened (+87.7° / +122.1° at rank 0) |
| constrain the torsion at embed (`SetDihedralDeg` + constrained minimize) | this is **construction**, against which the project holds three recorded negative results — and it would have constructed a twist to satisfy a token that should never have been emitted |

Both cheap remedies target conformer **sampling**. The failure was in conformer **perception**,
one layer up — and it also disabled the instrument meant to detect it: `_verify_axial_honored`
compares against the same blind `mol_axial_token`, so it reported nothing.

### The gate stays conservative on coupled axes

Lane 4 did **not** ship a coupled-axis canonicalization. `YESKOZ` now emits nothing, which means
the OIN currently encodes *syn* and *anti* meso-arylporphyrins identically. That is a real loss,
recorded in `docs/KNOWN_LIMITATIONS.md`, and it is a **smaller** loss than the status quo it
replaced (a non-reproducible sign asserting a chirality that does not exist).

The fix is specified but not guessed at: canonicalize the sign vector over the orbit of the
automorphism group's action on it — for a coupled pair that is `min(token, complement)`, which is
simultaneously reflection-invariant (correct: achiral) and *syn/anti*-separating (correct: real
diastereomers). The gate must distinguish **coupled** ambiguity (one automorphism inverts every
sign ⇒ quotient by global inversion, information survives) from **independent** ambiguity (each
axis has its own local C2 ⇒ quotient by independent flips, nothing survives). Per-axis rank
comparison cannot express that distinction.

---

## Dead ends and refutations

### D1 — "Widen the embed pool for multi-axis tokens"

**Killed by:** `tools/injectivity/axial_pool_histogram.py` on `YESKOZ`. Every pooled conformer
already held both hindered twists (rank 0: +87.7° / +122.1°) and the pool already spanned
several sign combinations (rank 1: +128.0° / −14.7°; rank 2: +157.2° / −124.5°). Widening adds
more of something already present.

### D2 — "Guard the FF relaxation so it cannot flatten the torsion"

**Killed by:** the same histogram. The torsions were never being flattened. The premise —
"torsions relax outside the hindered window (20–160°)" — is a description of a mechanism nobody
had measured.

### D3 — "Constrain the biaryl torsion at embed time"

**Killed by:** prior art, not a new measurement. The project holds three recorded negative
results for construction-over-selection. Additionally, it would have constructed a twist to
satisfy a token that (per R3) should never have been emitted for those molecules.

### D4 — "Fix it downstream: re-perceive the generated coordinates, or fix the contract-mol transfer"

**Killed by:** two measurements. (a) Re-perceiving the generated XYZ through the encoder's own
`get_tmc_mol` still yields **no axes**. (b) The OIN fragment SMILES the encoder emitted does not
survive RDKit sanitisation as aromatic — **12/40** atoms aromatic afterwards:

```
CSc1ccccc1-c1c2nc(cc3ccc(n3)c(...)c3nc(cc4ccc1n4)C=C3)C=C2
  -> SanitizeMol -> CSc1ccccc1C1=C2C=CC(=N2)C=C2C=CC(=N2)...   (12/40 aromatic)
```

So the descriptor had to stop depending on perception; it could not be repaired downstream.

### D5 — "Canonicalize the multi-axis token by sorting equivalent axes on their sign"

This is the **near-miss**, and it is the most transferable lesson in the lane. To make the
multi-axis string order-independent, an earlier (Y2) version broke rank ties with the sign
itself:

```python
axes.sort(key=lambda ax: (sym_key, ax.sign))   # destroys the chirality
```

Symmetry-equivalent axes tie on symmetry rank, so they sorted **by sign**, forcing signs into
ascending order. A molecule carrying `+-` then rendered identically to its mirror carrying `-+`:
the token silently **stopped being a chirality descriptor** for exactly the structures that need
one.

**Killed by:** the corpus-wide sign-convention audit (`axial_population.py --mirror-check`) —
mirror **every** emitting structure and require the token to flip. **34 of 37 flipped; 3 did
not** — OJELAQ, YESKOZ, EBUHAN, all multi-axis, all called chiral by the independent geometric
oracle (mirror RMSD 3.2–7.1 Å), so the failures were real rather than meso. **Fix:** break ties
on the **tie-broken canonical rank** — also graph-derived, hence still renumbering-invariant,
but it keeps each sign attached to its own axis so reflection flips the token. **After the fix:
37/37 flip.**

⚠ **The single-axis BINAP fixture could never have caught this** — one axis means no tie to
break. Regression cover is now a *constructed* two-axis probe (`_two_axis_probe()` in
`tests/unit/test_axial_emit.py`: two disconnected copies of 2,2′-dimethylbiphenyl, whose twists
can be set independently), specifically so the guard cannot pass by luck on a molecule whose two
axes happen to agree. `YESKOZ` was retired from that role because it no longer emits.

### D6 — "The corpus mirror audit validates the sign convention"

**Killed by:** the four-variant `YESKOZ` sign table (R3). For **coupled** axes the mirror
produces the exact complement, which the audit scores as "flipped correctly". So the audit
returns a pass for a token that is reflection-invariant — the very defect it exists to detect.
The audit is sound for independent axes and **vacuous for coupled ones**.

### D7 (context, not Lane 4's own) — three Y2 integration bugs whose common shape matters

Recorded here because the class recurs: two of the three produced a *plausible wrong answer with
no error*, and a wrong atropisomer still satisfies the round-trip key because the key folds the
token.

- **The token broke the generation parser.** `generation/oin_parser.py` decides
  inline-vs-sidecar by the *absence* of `|`, so a trailing ` |ax:-|` misrouted inline OINs to the
  sidecar branch (`Geometry code '' not supported`). The token is now stripped before parsing and
  preserved on `original_oin`.
- **Key-folding defeated acceptance.** The SL1 early-exit gate accepts the first conformer
  matching the fac/mer key — and the key deliberately folds the axial token, so it stopped the
  pool at *either* atropisomer before selection could choose. Acceptance is now key-match **and**
  axial-match. **Generalisable lesson: a deliberately-lossy comparison key must never be reused
  as an acceptance predicate for an axis it folds.**
- **Perception failed silently on pool conformers.** Raw MetalloGen pool mols are unsanitized, so
  `CanonicalRankAtoms` raises; a `try/except` turned that into `None`, and the filter compared
  `None` against `-` and matched nothing. Perception now runs on the **contract mol** with
  sanitize fallbacks, and `mol_axial_token` logs when every strategy fails so `None` means
  *could not tell*, never *no axial token*.

---

## Where it landed

**Lever:** `OIN_EMIT_AXIAL` — **default OFF** (in `levers.py::_HELD_OFF`, not `_DEFAULT_ON`).
With the lever unset, encoder output is byte-identical to pre-Y2.

**Code:** `src/oinsmiles/oin/axial.py` (single source of truth, shared by the encoder emit and
`tools/injectivity/config_oracle.py`); emit site in
`src/oinsmiles/utils/xyz2mol.py` (`_axial_suffix`, computed **before** `_align_to_pai` because
that may reflect the coordinates and a reflection inverts a chirality descriptor; appended at the
single `return inline_oin + _axial_suffix + _metal_config_suffix`).
Key fold: `src/oinsmiles/oin/compare.py::_AXIAL_TOKEN_RE`.
Generator: `_select_by_geometry`'s axial-aware narrow plus the accept-first `target_axial` check
in `src/oinsmiles/generation/metallogen_adapter.py`.

**What is NOT usable in the shipped default:** the axial descriptor itself is fully functional
and correct, but **nothing in the default configuration emits it**, so on the shipped default the
encoder remains blind to P2 exactly as Y1 measured. That is a deliberate product call, not a
defect — but it means the three `@expectedFailure` probes stay xfail and no default-path
measurement can confirm the descriptor. Only a raw string with the lever ON can.

**Guard tests (all measured green on `main`):**

| module | result |
|---|---|
| `tests/unit/test_axial_emit.py` | **17 tests OK** |
| `tests/unit/test_axial_failure_modes.py` | **13 tests OK** |
| `tests/unit/test_injectivity_probes.py` | 11 tests, `OK (expected failures=3)` |

Named guards worth knowing by name:

- `TestDefaultOff::test_binap_blind_by_default` — default OFF stays byte-identical.
- `TestDefaultOff::test_emit_gate_is_off_unless_the_env_var_is_set` — catches an accidental
  promotion; asserts unset ⇒ no emit, and explicit `"0"` ⇒ no emit.
- `TestFlagOn::test_binap_atropisomers_diverge_in_raw_string` / `test_key_still_folds_the_token`
  / `test_non_atropisomer_unaffected`.
- `TestTokenIsCanonical` — invariant under atom renumbering, invariant under proper rotation,
  flips under reflection.
- `TestSymmetryEquivalentAxesStillFlip::test_two_equivalent_axes_of_opposite_sign_flip` and
  `::test_same_sign_pair_is_not_confused_with_opposite_sign_pair` — the D5 regression, on the
  constructed probe.
- `TestPerceptionInvariance::test_binap_token_survives_delocalization` /
  `test_two_axis_probe_survives_delocalization` — the Lane 4 fix, unit-level.
- `TestPorphyrinMesoAxesAreNotPerAxisStereogenic::test_yeskoz_emits_no_token` — asserts the
  conclusion under **both** aromaticity perceptions.
- `TestNotOverSensitive::test_gate_agrees_with_independent_oracle` — cross-validates the
  stereogenicity gate against `tools/injectivity/oracle.py::geometric_chirality`, which never
  sees the gate's logic.
- `TestNormalizerFolds::test_normalize_strips_axial_suffix` — the key fold.

**Corpus-scale guard:** `tools/injectivity/axial_perception_sweep.py` applies a worst-case
perception perturbation (every non-metal bond flattened, aromaticity and charges cleared) to the
same coordinates and asserts the token is unchanged. Coordinates, elements and connectivity are
untouched, so handedness is untouched — only the perception is.

**Reproduce:**

```bash
cd /home/tjmustard/Documents/GitHub/OIN-SMILES
export PYTHONPATH=$PWD/src:$PWD
V=$PWD/.venv/bin/python

# the diagnosis: what does the embed pool actually contain?
$V -m tools.injectivity.axial_pool_histogram tests/fixtures/YESKOZ.xyz

# perception invariance (the Lane 4 guard)
$V -m tools.injectivity.axial_perception_sweep --fixtures
$V -m tools.injectivity.axial_perception_sweep --dataset default --n 400

# corpus population + sign-convention audit (--jobs is deterministic; --tag avoids
# overwriting a baseline scan)
$V -m tools.injectivity.axial_population --n 1500 --mirror-check --jobs 4 --tag skeleton

# guards
$V -m unittest tests.unit.test_axial_emit tests.unit.test_axial_failure_modes
```

### ⚠ Source-document discrepancies to be aware of

1. **`docs/agentic-notes/injectivity/INJECTIVITY_Y2_FEASIBILITY.md` is stale on the multi-axis attribution.** It still
   reads *"The multi-axis gap is a generator limit, not a sign error… The generator can hold one
   atropisomeric axis but not two"* and *"Follow-on: teach the generator to hold multiple hindered
   axes (a pool-diversity or torsion-constraint problem, not a notation one)."* Lane 4 **refuted
   both** (R1, R3). The corrected account lives in `docs/agentic-notes/v0.4.5/AXIAL_v0.4.5_LANE4.md` and
   `docs/KNOWN_LIMITATIONS.md` §3 of the axial section. Prefer those.
2. **`docs/agentic-notes/injectivity/INJECTIVITY_Y2_FEASIBILITY.md` also still carries the "gate ON by default"
   recommendation** as its closing status. That recommendation was not taken in v0.4.5; the
   confirmed product call kept it OFF.
3. **`spec/handoffs/v0.4.5/Lane4-axial-gate-multiaxis.md` Part A ("flip the default") was
   deliberately declined** and the handoff was not amended.
4. **`docs/agentic-notes/injectivity/INJECTIVITY_Y1_P2_AXIAL.md` still describes the token's stereogenicity gate as
   operating "on the intact graph"** in spirit — it predates `_axial_cut_ranks`. It is not wrong
   about the *conclusion*, only about the mechanism.

---

## Open questions / for the next agent

### The promotion gate for `OIN_EMIT_AXIAL` — three conditions, not one

**Condition 1 — re-measure the cohort under `OIN_CANONICAL_PERCEPTION`.**
`OIN_CANONICAL_PERCEPTION` was promoted **default-ON** in v0.4.5, and it **interacts with the
axial descriptor**. `axial.py::_is_atropisomer_candidate` gates its steric wall on
`not oo.GetIsAromatic()` — the one part of the descriptor that still reads an aromatic flag, and
deliberately so (dropping the aromatic term makes every twisted biphenyl "walled" because the
ortho ring carbon itself becomes a wall; tightening it to "exocyclic substituent" alone would
stop a meso-arylporphyrin qualifying, since its porphyrin end is walled by fused pyrroles).

Measured on `YESKOZ`: **hindered axes go 2 → 1** under canonical perception, because more of the
macrocycle reads aromatic and one porphyrin wall stops qualifying. **No emitted string moves
today** — YESKOZ's two meso-aryl axes are non-stereogenic so its token is empty either way, and
`PdCl2-R-BINAP` is unaffected at 1 hindered / 1 emitting. **But the Y2 cohort numbers
(single-axis 22/22, corpus mirror audit 37/37) were all measured with perception OFF**, and
`axial.py`'s own safety argument only covers the **generator** reading *fewer* atoms as aromatic
(which can only ever *add* an axis, producing a loud token-length mismatch), **not the encoder
reading more**. Both cohorts must be re-measured with perception ON before promotion. This is
recorded at the gate itself in `levers.py::_HELD_OFF["OIN_EMIT_AXIAL"]` and pinned by the
two-perception subTest in `TestPorphyrinMesoAxesAreNotPerAxisStereogenic::test_yeskoz_emits_no_token`.

**Condition 2 — remove the key's fold in the same commit.** `compare.py::_AXIAL_TOKEN_RE` strips
` |ax:±|` before comparison. **A key that folds an axis is not a valid acceptance predicate for
that axis.** Today (lever OFF) the fold is a no-op that keeps the harness immune to a developer
running with the flag set — keep it. Once the lever is ON, folding makes the round trip
structurally unable to verify the one thing the token encodes: the generator would be free to
return the wrong atropisomer and the key would still match, leaving the failure only in
telemetry. This is no longer hypothetical — Lane 7's `invert_stereocenter` twin operator flips a
*single* axis of a multi-axis molecule (something mirroring cannot do) and on `YESKOZ` reports
`encoder_blind` by default but **`key_blind` with `OIN_EMIT_AXIAL=1`**: the raw strings differ and
`_AXIAL_TOKEN_RE` is the only thing collapsing them. **Decision recorded: unfold in the same
commit that promotes the lever.**

**Condition 3 — re-measure the population, which nobody has done.** The 2.49% emitting rate,
the 56/19 hindered/non-stereogenic split and the 37/37 mirror audit were all measured against
the **pre-`10214f3e` descriptor**, which used aromatic flags and ranked atoms on the mol as
perceived. `10214f3e` changed axis selection (trigonal ring atom), the reference-neighbour set
(ring neighbours), the rank basis (connectivity skeleton) and the stereogenicity scope (axis bond
cut) — and it demonstrably changed *which* structures emit (YESKOZ went from emitting `+-` to
emitting nothing). `docs/agentic-notes/v0.4.5/AXIAL_v0.4.5_LANE4.md` §6 documents the re-run command
(`axial_population --n 1500 --mirror-check --jobs 4 --tag skeleton`) but **no post-fix number is
recorded anywhere in the repo** — verified by grepping `docs/` for the `skeleton` tag. So the
population and the audit both need re-taking before any promotion claim rests on them. Run it
together with Condition 1 (perception ON) so one scan settles both.

### The genuinely open scientific question: coupled-axis canonicalization

The *relative* sign across a coupled pair is real information (`++` anti vs `+-` syn); the
absolute per-axis sign is not. The conservative gate shipped in v0.4.5 drops **both**, so
`YESKOZ`-type *syn* and *anti* meso-arylporphyrins currently encode identically.

The specified fix: canonicalize the sign vector over the orbit of the automorphism group's
action on it — `min(token, complement)` for a coupled pair. It requires distinguishing
**coupled** ambiguity (one automorphism inverts every sign; quotient by global inversion,
information survives) from **independent** ambiguity (each axis has its own local C2; quotient
by independent flips, nothing survives). Per-axis rank comparison cannot express that
distinction; it needs the automorphism group's action on the sign vector, enumerated with a
blow-up guard and a conservative fallback.

**The right validation instrument is Lane 7's `invert_stereocenter`, not mirroring** — mirroring
cannot isolate one axis of a multi-axis molecule, which is exactly why the mirror audit passed
vacuously (D6).

### Smaller residuals

- **Sign skew, not excluded.** Among the emitting examples `+` dominates roughly 3:1. The benign
  explanation is chemical (catalysis datasets are rich in commercially dominant
  single-enantiomer ligands such as (R)-BINAP), but a systematic skew in the sign convention has
  **not** been excluded. Worth one measurement, not an over-read.
- **`_verify_axial_honored` is only as good as `mol_axial_token`.** It compares against the same
  perception path it is supposed to police. It reported nothing throughout the multi-axis
  failure. Consider a second, independent signal (e.g. the raw dihedral) before trusting
  `adapter.axial_not_honored` as the sole promotion instrumentation.
