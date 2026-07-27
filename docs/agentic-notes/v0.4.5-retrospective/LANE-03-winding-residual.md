# Lane 3 — The eta winding residual

**What this lane was for:** drive the `winding_star_drift` bucket (13 molecules on the v0.4.2
capstone sweep) to zero by making eta-ring **heading-atom selection** and the **`>`/`<` winding
sign** embedding-invariant — replacing load-bearing *geometric* tiers with topological /
canonical-rank rules.

---

## ELI5

When a metal binds a whole flat ring (a "sandwich" complex like ferrocene), the OIN string records
*which face* of the ring the metal is looking at, as a single character: `>` or `<`. To compute it
you need a reference atom on the ring (the "star" or heading atom) and then you ask whether the ring
reads clockwise or counter-clockwise from the metal's point of view. Both halves could wobble: the
star was sometimes picked by geometry, so it landed on a different ring atom in a crystal structure
than in a computer-generated copy of the same molecule, and the character then moved with it. The
lane went in expecting to replace geometry with graph rules everywhere. It found that the sign math
was already rotation-proof, that only one molecule actually suffered from the geometric star, and
that most of the remaining "drift" was the encoder *correctly* reporting that the generated
structure was the **mirror image** — a genuinely different molecule. Folding those away would have
destroyed real chemistry, so the lane refused, and said so with numbers.

## The work, visually

```
   eta group (a ring bound face-on)              WHY THE MARKER CAN MOVE
   ┌───────────────────────────────┐             ─────────────────────────
   │  step 1: pick a HEADING atom  │             two independent halves:
   │  step 2: read the CIRCULATION │              (a) which atom is the star
   └───────────────────────────────┘              (b) which way it circulates

  ── FIX A (default ON) — oin_aligner.py step 4a, heading atom ───────────────────
     BEFORE                                    AFTER
     content-canonical tier                    both branches call
       ├─ orientation-free branch              _topological_heading_atom:
       │    → _canonical_heading_atom            ├─ tier 1  strict CanonicalRankAtoms
       └─ other branch                           ├─ tier 2  valence-tolerant symmetry
            → _canonical_heading_atom            │          graph (_orientation_symmetry_graph)
              └─ None ⇒ FELL THROUGH TO          └─ tier 3  lowest constituent index
                 the GEOMETRIC loop  ⚠         all three are functions of the GRAPH alone
                 (GIPDEQ: boron ylide
                  won't kekulize ⇒ None)

  ── the sign itself — ALREADY invariant, do not "fix" ──────────────────────────
     oin_aligner._determine_winding  ──►  oin/winding.py::signed_circulation
        sign = cross(v_star, v_next) · axis ,  axis = that ring's OWN metal→centroid
        all three factors rotate together ⇒ invariant by construction
        ⚠ signed_circulation is SHARED with the generator (haptic-face correction,
          Stereo Phase 3). Any change must keep both sides consistent or eta
          round-trips break.
     ⚠ orientation-free rings (Cp, arene, BPh4- phenyl: some automorphism REVERSES
       the cyclic order) are forced to a fixed '>'. That rule is LOAD-BEARING —
       it is what makes those rings' notation embedding-free. Do not simplify away.

  ── FIX C (OIN_CANONICAL_ETA_WINDING) — _canonical_eta_winding ─────────────────
     1. collect eta groups whose winding is load-bearing (skip orientation-free)
     2. colour every occupied slot by (chem_id, eta automorphism class)
              ⚠ winding EXCLUDED from the colour — we are deciding which slots may be
                relabelled, which must not depend on the characters being reassigned
     3. keep proper rotations preserving every slot colour (_brute_force_symmetries,
        all from Rotation.from_euler ⇒ no reflection can enter)
     4. for each 2-ORBIT of eta slots in one automorphism class:
              ε = _eta_swap_sense(a, b)          ← COMPUTED, never assumed
              ┌──────────────┬─────────────┬────────────────────────────────┐
              │      ε       │ signs equal │ verdict                        │
              ├──────────────┼─────────────┼────────────────────────────────┤
              │  +1 (two     │  differ     │ ACHIRAL ⇒ emit min(spelling,   │
              │  separate    │             │           mirror)              │
              │  copies)     │  same       │ CHIRAL  ⇒ LEAVE ALONE          │
              │  −1 (two     │  same       │ ACHIRAL ⇒ fold                 │
              │  rings in    │             │                                │
              │  ONE bridge) │  differ     │ CHIRAL  ⇒ LEAVE ALONE          │
              └──────────────┴─────────────┴────────────────────────────────┘
     5. FAIL SAFE on anything else — orbit ≠ 2, missing automorphism, ε = None
        ⇒ today's behaviour
```

Legend: `⚠` = a rule that is load-bearing and must not be "simplified". Fix A is unconditional;
Fix C is `OIN_CANONICAL_ETA_WINDING`, default-**ON** since v0.4.5 (shipped default-OFF in the lane).

## Initial assumptions and hypothesis

1. **The lane's stated premise:** the 13 `winding_star_drift` residuals come from the **geometric
   fallback tiers** — the heading atom and/or the sign being read off the 3D embedding, so a crystal
   structure and its regenerated twin disagree.
2. **The stated target:** `winding_star_drift` 13 → 0, by eliminating geometric fallbacks rather
   than making them deterministic.
3. **Assumed:** the drift class was an encoder canonicality defect throughout, so
   `tools/canonicality_probe.py` (rotate / renumber, graph held fixed) would be able to see it.

## What was actually found

**Refuted — the published list of 13 is stale.** `winding_star_drift = 13` comes from the **v0.4.2**
capstone sweep (`results-capstone-v042/bucket_report.json`). Re-encoding those same 13 stored
geometry pairs on current `main` reclassifies them:

| current bucket | n | molecules |
|---|---:|---|
| `key_equal` / **`winding_star_drift`** | **6** | GAMJAG, GIPDEQ, QIGZAJ, SIMDIE, TEYXEA, WUHRIB |
| `key_equal` / `rdkit_canonical` | 3 | BOJMAQ, MOHKOL, SATKIJ |
| `structural` | 3 | ABETIK, IFICAD, OVUBEO |
| `facmer_divergent` | 1 | FAHCIC |

v0.4.3/v0.4.4 had already moved 7 of 13 out of the class. **The true worklist was 6.** (BOJMAQ is
not a winding case at all — its two eta groups keep identical per-slot windings and merely trade
slot numbers, i.e. Lane 2. MOHKOL/SATKIJ show the sign pattern but with body drift on top, so
Lane 1 owns them first.)

**Refuted — the premise is wrong for 5 of the 6.** `tier` below is the tier that chose the heading
atom (`tools/winding_diag.py`), measured separately on the input and on the stored generated
structure; `rot-inv` = the emitted string is identical across N random **proper** rotations.

| molecule | symptom | deciding tier (input / generated) | rot-inv | `generated == mirror(input)` | root cause | verdict |
|---|---|---|:---:|:---:|---|---|
| **GIPDEQ** | star moves (B → CH₂), sign unchanged | **2-GEOMETRIC / 2-GEOMETRIC** | yes | no | `_canonical_heading_atom` returns `None` (boron ylide `[CH2]=c1cc(C)cc(C)c1=B(...)` will not kekulize) so Tier 1 fell through to the geometric heading, which tracks the embedding | **fixed** (Fix A) |
| **GAMJAG** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | two **separate, identical** benzylindenyl fragments; ε = +1, signs opposite ⇒ **achiral**, so the two spellings are one compound | **fixed** (Fix C) |
| **SIMDIE** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | Me₂Si-bridged bis(2-Me-4-Ar-7-OMe-indenyl); ε = −1, signs opposite ⇒ **chiral**: input and generated are **enantiomers** | **not an encoder bug** |
| **TEYXEA** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | as SIMDIE (ε = −1, opposite ⇒ chiral) | **not an encoder bug** |
| **WUHRIB** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | as SIMDIE (ε = −1, opposite ⇒ chiral) | **not an encoder bug** |
| **QIGZAJ** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | Cp ↔ fluorenyl are **not** automorphic (ε undefined) ⇒ **genuine enantiomers** | **not an encoder bug** |

Plus one **latent** case the characterization surfaced that was never in the drift list:
**MECJOU** (ε = −1, signs *same* ⇒ achiral) encodes differently from its own mirror image. Input
and generated happened to agree, so no sweep ever flagged it — same defect as GAMJAG, and Fix C
closes it too.

So the residual is **three distinct causes, not one**: geometric heading fallback (1 case);
un-canonicalized reflection freedom on an **achiral** eta arrangement (GAMJAG + latent MECJOU); and
a **generator stereochemistry error the key over-folds** (SIMDIE, TEYXEA, WUHRIB, QIGZAJ — 4 of 6),
where the encoder is *right* and the generated structure is the enantiomer. **The majority of this
lane's residual is not an encoder defect at all.**

**Confirmed — the sign math was already rotation-invariant.** `signed_circulation` computes
`cross(v_star, v_next) · axis` with `axis` the ring's *own* metal→centroid vector, so all three
factors rotate together and the sign is invariant by construction. Empirically: N random proper
rotations of every one of the 6 cases give **one** distinct string, and the independent corpus-wide
`canonicality_probe.py` run found **zero `rotate` drift**.

**Therefore — `winding_star_drift` is invisible to `canonicality_probe.py`.** That probe holds the
structure fixed and varies presentation; these residuals are *input structure vs a differently
generated structure of the same compound*, a different question. The acceptance number has to come
from the `reencode_ab` + bucket-report `key_equal` sub-split — the part of that instrument the
orchestrator's trust gate confirmed (12.04% vs the sweep's 12.32%).

**Refuted — RMSD cannot answer the question that matters.** Every one of the 6 satisfies
`generated == mirror(input)` byte for byte, so the question is not "did the encoder drift?" but
"**is `mirror(M)` the same compound as `M`?**". `tools/injectivity/oracle.py::is_distinct_enantiomer`
answers "chiral" for all 6 (RMSD 2.6–4.4 Å, automorphism cap hit) — but that is *conformational*
chirality: pendant aryls, tBu rotors and freely rotating rings mean a genuinely achiral meso complex
also fails to superimpose on its mirror. `tools/eta_core_chirality.py` was built to strip that away
(metal + coordinated ring systems + bridge + one substituent shell) and **does** calibrate on rigid
ferrocene (ACHIRAL, 0.103 Å, 200 core automorphisms — *after* flattening formal charge, which
perception localizes onto one arbitrary Cp carbon and which otherwise collapses the automorphism
group to 8 and reports ferrocene as chiral at 0.81 Å). It still cannot decide a metallocene whose
rings rotate freely about the metal axis, so it is **not** the oracle for this lane; it is kept, with
that limitation documented, as a correct instrument for rigid cores.

**Confirmed — the exact criterion, and the sense factor ε.** Two OIN spellings denote the same
compound exactly when related by a **proper rotation of the coordination polyhedron** composed with
a **graph automorphism of the ligand assembly**; both carry each ring's winding character along
unchanged, so `(ring identity, winding char)` travels as a unit and the only freedom is which
interchangeable slot each unit is written at. But the character is read off **each ring's own
canonical SMILES order**, so an automorphism carrying ring A onto ring B hands the character over
**possibly reversed**. Applying the slot-swapping rotation together with the automorphism maps
`(w₀, w₁)` to `(ε·w₁, ε·w₀)`, which equals the mirror spelling `(−w₀, −w₁)` exactly when:

| ε | the achiral arrangement is | typical structure |
|---:|---|---|
| **+1** | **opposite** signs `(>,<)` | two separate copies of one fragment — canonically ordered identically |
| **−1** | **same** signs `(>,>)` | two rings inside ONE bridged fragment — the canonical SMILES runs out along one ring and back along the other |

**ε is not a constant** — measured `+1` for the TiCat3/TiCat4 ligand, `−1` for SIMDIE and MECJOU.
Verdicts from the computed ε agree with the independent geometric oracle **6/6** (7 rows, QIGZAJ
being the not-automorphic "no fold" case):

| case | signs | ε | ε-verdict | oracle |
|---|---|---:|---|---|
| TiCat3 | `(>,>)` | +1 | CHIRAL | CHIRAL (2.14 Å) |
| TiCat4 | `(<,>)` | +1 | ACHIRAL | ACHIRAL (0.03 Å) |
| SIMDIE | `(<,>)` | −1 | CHIRAL | CHIRAL (2.78 Å) |
| MECJOU | `(>,>)` | −1 | ACHIRAL | ACHIRAL (0.04 Å) |
| TEYXEA | `(<,>)` | −1 | CHIRAL | CHIRAL (2.65 Å) |
| WUHRIB | `(>,<)` | −1 | CHIRAL | CHIRAL (1.45 Å) |
| QIGZAJ | `(>,<)` | — (not automorphic) | no fold | CHIRAL (2.70 Å) |

**Found — the key over-folds an enantiomer pair (a finding handed to Lane 2).** QIGZAJ buckets as
`key_equal`, i.e. the key says input and generated are the same isomer. They are not.
`compare.py::_parse_vertex_colors` colours every slot marker in a fragment with that **whole
fragment's** body, so for an ansa ligand the Cp slot and the fluorenyl slot receive the *same*
colour; the tetrahedral C2 that swaps slots 0↔1 then looks colour-preserving and
`_polyhedron_signature` folds the two spellings. Per-eta-ring identity (as in
`_eta_automorphism_class`) would fix it. **Lane 3 deliberately does not fold QIGZAJ**, so it remains
a visible `winding_star_drift` — correctly.

**Measured — Fix A (default ON), unmodified `main` vs the branch with the lever OFF.** 1176 common
re-encoded pairs. **2 molecules changed their emitted string**, both already `structural` in both
arms, so **no bucket moved**. GIPDEQ goes `winding_star_drift` → `byte_exact`. Scoped exactly as
predicted: only fragments where the strict canonical rank is unavailable.

**Measured — Fix C over its complete provable scope.** Fix C can only rewrite a winding character
inside a 2-orbit of eta slots, so every molecule it could touch has **≥ 2 eta winding markers**; all
**500** such corpus molecules were re-encoded in both arms:

| bucket | lever OFF | lever ON |
|---|---:|---:|
| `byte_exact` | 203 | **204** |
| `key_equal` / `winding_star_drift` | 5 | **4** |
| `key_equal` / `rdkit_canonical` | 54 | 54 |
| `key_equal` / `slot_renumber` | 55 | 55 |
| `facmer_divergent` | 52 | 52 |
| `structural` | 131 | 131 |

- **emitted string changed: 3 / 500** — all three `byte_exact` with the lever on (GAMJAG's fix plus
  two latent MECJOU-type recanonicalizations of pairs that already agreed);
- **bucket changed: 1** — GAMJAG, `winding_star_drift` → `byte_exact`;
- **`facmer_divergent` did not rise (52 → 52)**: no isomer over-folding.
- A broader 1082-molecule sample outside this scope showed **exactly one** string change (GAMJAG),
  consistent with the scope argument.

Suite at the time: `tests/unit/test_winding_canonical.py` (6 tests) plus the pre-existing
`tests/unit/test_winding_helper.py` green; full suite **605 OK / 3 skip / 3 xfail**; ruff clean.

## What was done

**Fix A — eliminate the geometric heading tier** (unconditional, `utils/oin_aligner.py` step 4a,
~line 1428). Both branches of the content-canonical tier now call
`OINDiscreteAligner._topological_heading_atom`, which has its own three sub-tiers — strict
`CanonicalRankAtoms` (`_canonical_heading_atom`) → valence-tolerant symmetry graph
(`_orientation_symmetry_graph`, which tolerates a BPh₄⁻ boron's four bonds where strict
sanitization returns `None`) → lowest constituent index — **all** functions of the graph alone.
Previously the non-orientation-free branch called `_canonical_heading_atom` directly and, on `None`,
fell through to the geometric loop. Why this cannot destabilize a pair that is byte-exact today: the
change is scoped to fragments where the strict rank is unavailable, a population that is by
construction already embedding-unstable, and identical strings imply identical fragment graphs while
the topological heading is a function of the graph. After Fix A the geometric heading loop is **dead
for eta groups**; it is retained only as a fail-safe. This follows the lane's instruction to
*eliminate* a geometric fallback rather than make it deterministic.

**Fix C — canonical winding across automorphic eta groups** (`OIN_CANONICAL_ETA_WINDING`, shipped
default-OFF by the lane), `OINDiscreteAligner._canonical_eta_winding` (~line 1616), hooked as step
4c and returning `{(rank, slot, heading_idx): winding_char}` overrides only. The five steps are in
the diagram above. Supporting helpers:

- `_eta_automorphism_class(smiles, constituent_indices)` — the canonical SMILES of the fragment with
  **every constituent atom given the same atom-map number**. Two eta groups match exactly when some
  automorphism of the fragment carries one constituent set onto the other. It works *across*
  fragments (two identical separate indenyls) and *within* one (a bridged pair), and it deliberately
  skips sanitization so an unkekulizable borate still gets an id. Returns `None` on failure, and a
  `None` id never joins any class.
- `_eta_swap_sense(a, b)` — **computes** ε. Different fragment SMILES ⇒ `None` (not automorphic, no
  fold permitted). Same SMILES on *different* fragments ⇒ `+1` (two copies of one ligand, canonically
  ordered identically). Same fragment ⇒ enumerate real self-matches with
  `GetSubstructMatches(mol, uniquify=False, useChirality=False, maxMatches=20000)` and test whether
  the image of A is a cyclic rotation of B (`+1`) or of `reversed(B)` (`−1`).

**Design choices, and the alternatives rejected.**

- **Winding is excluded from the slot colour** in step 2: the question being asked is *which slots
  may be relabelled*, and the answer must not depend on the very characters about to be reassigned.
- **Only proper rotations enter.** `_brute_force_symmetries` is built from `Rotation.from_euler`, so
  no reflection can appear; combined with graph automorphisms these are the two operations that
  preserve molecular identity by definition.
- **The fold is gated on a MEASURED achirality test, not an assumed one.** Concretely, on an
  ansa-metallocene the **rac** diastereomer carries the same character on both rings so sorting is a
  no-op and its mirror stays a different string; only **meso**, whose two mirror-related spellings
  really are one achiral compound, collapses. A Cp/fluorenyl pair is not automorphic, so each ring
  sits in a singleton orbit and nothing moves — correctly, because for an unsymmetrical bridge the
  two spellings are genuine enantiomers.
- **Scoped to 2-orbits, failing safe otherwise.** For a larger orbit the rotation group may realize
  only *some* rearrangements, so the reachability argument would need the full induced-group
  computation; 2 is every case the corpus presents.
- **`signed_circulation` (`oin/winding.py`) was left untouched.** It is the single source of truth
  **shared with the generator** (the OIN→XYZ haptic-face correction, Stereo Phase 3); changing the
  sign convention on one side breaks eta round-trips. The sign was measured invariant anyway, so
  there was nothing to fix.
- **The fixed-`>` rule for orientation-free rings is load-bearing and must not be "simplified
  away".** A ring that some automorphism carries onto itself in *reverse* cyclic order can be turned
  over by a proper rotation, so its geometric winding records only which face this particular
  embedding presented — it flips at random between an input structure and a regenerated one.
  Emitting the degenerate `>` is what makes those rings' notation a function of the structure alone.
  Note the test is strictly weaker than "every ring atom is in one symmetry class": mesitylene's
  ring has two classes and an arm-substituted Cp* has four, yet both are orientation-free, so
  testing symmetry classes would silently leave them broken.

## Dead ends and refutations

| tried | what killed it |
|---|---|
| "the 13 residuals come from the geometric fallback tiers" (the lane's charter premise) | measured per molecule: **5 of the 6** survivors never reach a geometric tier (`1a-strict / 1a-strict`). Only GIPDEQ does |
| the published `winding_star_drift = 13` as the worklist | stale (v0.4.2 sweep). Re-encoding the same 13 pairs on current `main` gives 6 winding / 3 `rdkit_canonical` / 3 `structural` / 1 `facmer_divergent` |
| target `winding_star_drift → 0` | **refused, with numbers.** 4 of the 6 are cases where the *generator* produced the enantiomer and the encoder correctly said so. Driving the class to 0 in the encoder means folding enantiomers |
| "make the sign embedding-invariant" | already was: `cross(v_star, v_next) · axis` with the ring's own axis; N random proper rotations give one string per case, and the corpus probe shows **zero** `rotate` drift |
| `canonicality_probe.py` as the acceptance instrument for this lane | it cannot see the class at all — it holds the structure fixed and varies presentation, whereas these residuals are input structure vs a differently *generated* structure |
| RMSD-based `is_distinct_enantiomer` as the chirality oracle | measures *conformational* chirality: says "chiral" for all 6 (RMSD 2.6–4.4 Å) including the achiral ones, because pendant aryls and rotors prevent superposition |
| `tools/eta_core_chirality.py` as the oracle | calibrates correctly on rigid ferrocene (ACHIRAL, 0.103 Å, 200 core automorphisms) but cannot decide a metallocene whose rings rotate freely about the metal axis. Kept for rigid cores only. Also note: without flattening formal charge it collapses ferrocene's automorphism group to 8 and calls it chiral at 0.81 Å |
| **first implementation of Fix C: just SORT the two winding characters** | made the encoder **reflection-invariant for every bridged case** and folded a rac pair — the v0.4.4 axial failure reproduced exactly. It passed every guard written against the easy fixture. What caught it was building the independent oracle and finding it reported ACHIRAL for the same-sign case and CHIRAL for the opposite-sign one — the **inverse** of the assumption. ε is now computed per fragment |
| folding QIGZAJ so the bucket would look better | its two spellings are genuine enantiomers (Cp ↔ fluorenyl are not automorphic). Left deliberately visible as `winding_star_drift` |
| tests spelling "lever off" by **deleting** the env var | after promotion that means ON (17 failures in v0.4.5, 6 more in v0.4.6). `test_winding_canonical.py` sidesteps it entirely by patching the module attribute `oin_aligner.CANONICAL_ETA_WINDING`, which is read at import |

## Where it landed

- **Fix A: unconditional**, `utils/oin_aligner.py` step 4a via
  `OINDiscreteAligner._topological_heading_atom` (line 1175, called at line 1428).
- **Fix C: lever `OIN_CANONICAL_ETA_WINDING`**, shipped default-OFF by the lane and **promoted to
  default-ON in v0.4.5** (`oin/levers.py::_DEFAULT_ON`). The module-level flag
  `oin_aligner.CANONICAL_ETA_WINDING = lever_enabled("OIN_CANONICAL_ETA_WINDING")` is read at
  import, so tests patch the attribute, not the environment.
- **Promotion evidence:** `docs/agentic-notes/v0.4.5/PROMOTION_GATE_v0.4.5.md` — all six canonicality levers together
  took byte-stability under rotation/renumbering from **58.1% (173/298) to 69.6% (208/299)** and
  comparison-key instability from **60 molecules to 16** on a 300-molecule seed-42 sample, with the
  over-folding veto (`test_facmer_key.py`, `tests/integration/test_isomer_divergence.py`) passing.
- **Result:** `winding_star_drift` **6 → 2** (GIPDEQ by Fix A, GAMJAG by Fix C; latent MECJOU also
  closed). Not 0, and deliberately so.
- **Guards:** `tests/unit/test_winding_canonical.py` —
  `TestWindingOrientationInvariance::test_eta_encoding_is_invariant_under_random_proper_rotations`
  (over `Ferrocene.xyz`, `Ferrocene-halide-face.xyz`, `TiCat1.xyz`, `TiCat3.xyz`, `TiCat4.xyz`, 3
  random proper rotations each, no hardcoded XYZ hash — a prior wave flaked CI that way);
  `TestEtaWindingCanonicalizationPreservesChirality::{test_rac_and_meso_bis_indenyl_stay_distinct_with_lever_on,
  test_chiral_eta_structure_still_differs_from_its_mirror_with_lever_on,
  test_achiral_eta_structure_folds_onto_its_mirror_only_with_lever_on,
  test_lever_leaves_orientation_free_metallocenes_untouched}`;
  `TestEtaSwapSense::test_sense_is_measured_per_fragment_not_assumed`. Supporting pre-existing
  suites: `tests/unit/test_winding_helper.py` (the `signed_circulation` contract shared with the
  generator), `tests/unit/test_automorphic_ring_winding.py` (the orientation-free classification —
  `test_mesitylene_is_free_despite_two_symmetry_classes`,
  `test_ansa_bis_indenyl_rings_are_bearing`, `test_borate_phenyl_is_free`),
  `tests/unit/test_eta_winding_generalization.py`.
  `TiCat3`/`TiCat4` are load-bearing fixtures precisely because they differ **only** in eta winding
  — same metal, same geometry, same ligand bodies.
- **Commits.** Branch `swimlane/v045-lane3`, tip `1d02ecf6`, merged into `main`: `3d90bf3a`
  (topological eta heading + measure the winding residual), `149fb15e` (chirality-safe eta winding
  canonicalization, gated OFF), `c1f644d7` (correct the winding-residual analysis to the measured
  result), `1d02ecf6` (WIP doc results — the agent halted mid-task). Shared tooling from Lane 1:
  `19d20042`, `7b85e123`, `20044883`. Promotion `1450b5ce`, release `0d165845`.

## Open questions / for the next agent

1. **The 4 remaining cases are a GENERATOR defect, not an encoder one.** SIMDIE, TEYXEA, WUHRIB and
   QIGZAJ: the generator does not reproduce the *coordinated face* of a bridged eta ring, so it
   returns the wrong diastereomer/enantiomer — the same shape of failure Lane 4 found for
   multi-axis atropisomers. **Next measurement:** for one of them, generate a pool and count how
   many pool conformers carry each `(w₀, w₁)` sign pair; if the correct pair is present but not
   selected, this is a `_select_by_geometry(honor_winding=True)` selection bug, not an embedding
   bug. That single measurement distinguishes the two repairs.
2. **The key should stop folding QIGZAJ** (Lane 2's territory).
   `compare.py::_parse_vertex_colors` colours every slot in a fragment with the *whole fragment's*
   body; per-eta-ring colour, as in `_eta_automorphism_class`, would separate an ansa ligand's two
   eta slots. **Next measurement before changing it:** re-run the `facmer_divergent` count over the
   500 ≥2-eta-marker corpus subset with per-ring colouring, since a finer colour can only *split*
   and splitting moves molecules out of `key_equal` into `facmer_divergent`/`structural` — quantify
   that before shipping.
3. **Orbits larger than 2 are unhandled** (`len(members) != 2` fails safe). No corpus molecule
   presents one today. **Next measurement:** count corpus molecules with ≥3 eta winding markers in
   one automorphism class before building the full induced-group computation.
4. **`tools/eta_core_chirality.py` cannot decide free-rotor metallocenes.** If a general achirality
   oracle is needed, the missing piece is a torsion-quotient superposition (allow ring rotation about
   the metal axis) — and it must be re-calibrated on ferrocene (expect ACHIRAL ≈ 0.10 Å with 200
   core automorphisms, *after* flattening formal charge) plus the TiCat3/TiCat4 rac/meso pair before
   it is trusted.
5. **Do not "simplify" two rules.** (a) The fixed `>` for rings with a reversing automorphism is
   load-bearing. (b) `oin/winding.py::signed_circulation` is shared with the generator; any change
   must land on both sides in one commit or eta round-trips break.
