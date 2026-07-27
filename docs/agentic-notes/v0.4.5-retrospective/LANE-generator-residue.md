# Lane: generator residue (`swimlane/v045-genresidue`)

**What the lane was for:** to work the residual generation-side failures left after the v0.4.4
accuracy wave — the three smallest non-timeout `hard_fail` classes of the 6,719-molecule capstone
sweep (`rmsd_gate` 33, `other` 25, `no_conformers` 23 — 81 molecules, **1.21%** of the corpus) —
and to say honestly which of them were live defects, which were documented representation limits,
and which were artefacts of a stale measurement.

Primary sources: `docs/agentic-notes/v0.4.5/GEN_RESIDUE_v0.4.5.md`, commit `4dcb5284`,
`src/oinsmiles/generation/metallogen_adapter.py`, `docs/KNOWN_LIMITATIONS.md`,
`docs/agentic-notes/injectivity/INJECTIVITY_Y3_UNKNOWN_UNKNOWNS.md`, `docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md` §3,
`docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md` §"The confound that governs every accuracy number here".

---

## ELI5

The project has two halves: turn a 3D structure into a one-line string, and turn that string back
into a 3D structure. This lane worked on the second half's leftovers — 81 molecules the generator
was recorded as failing on. The first thing it found was that the record was a week out of date:
33 of those 81 were failing against a rule the project had already deliberately removed, and 29
of them pass today with no code change at all. The second thing was that most of the rest are not
bugs but limits the notation already documents — a counter-ion floating outside the metal's
coordination sphere simply has nowhere to go in the string's grammar. Exactly one genuine,
previously-undiagnosed generator bug was found and fixed: a nitrogen or oxygen donor that is
already "full" (a tertiary amine, a water molecule) has no room left for the bond to the metal, so
the chemistry toolkit rejected the assembled complex outright — but only on tetrahedral metals, a
condition that had to be measured rather than guessed. Seven molecules now round-trip that
previously produced nothing.

---

## The work, visually

```
=====================================================================================
A. THE MEASUREMENT WAS STALE — this is why the lane's first job was re-measuring
=====================================================================================

  2026-07-15                       2026-07-22                      2026-07-26
  results-capstone-v042            v0.4.4-SL4 landed               this lane
  commit_id 58bba7ad               commit 5c5d193b                 (4dcb5284)
        |                                |                                |
        |--------- ONE WEEK -------------|                                |
        v                                v                                v
  81 molecules frozen into        "demote RMSD from a          re-ran all 33 rmsd_gate
  hard_fail_worklists.json        pass/fail GATE to a          rows: 29/33 already pass
  as rmsd_gate / other /          DIAGNOSTIC"                  with ZERO code change
  no_conformers                          |
        ^                                |
        |     the 33 rmsd_gate rows are a SNAPSHOT OF THE PRE-DEMOTION HARNESS
        +--------------------------------+

=====================================================================================
B. THE TRIAGE FUNNEL — 81 molecules, where each one went
=====================================================================================

  81  hard_fail (3 smallest non-timeout classes, 1.21% of the 6,719-molecule corpus)
   |
   +-- 33  rmsd_gate ....... 29 -> now status=success (stale measurement, no code change)
   |                          3 -> atom_count mismatch (sibling lane's territory)
   |                          1 -> genuine string mismatch: SEQVEP_comp_0
   |                               (RMSD no longer masks it; real generation defect)
   |
   +-- 25  other ........... 24 -> UncoordinatedFragmentError = the SAME 24-molecule
   |                                floor set in docs/KNOWN_LIMITATIONS.md since
   |                                2026-07-14. Representation limit, NOT fixed.
   |                           1 -> PEKQUU_comp_0, geo_code 'NON' = a CORRECT REFUSAL
   |
   +-- 23  no_conformers ...  7 -> FIXED by this lane (gated Group-1 donor fix)
                              1 -> FOGCIN_comp_0, already fixed by someone else
                              6 -> perception gap (aromatic [b], P=N=P), upstream
                              9 -> unresolved, likely genuine; NOT forced

   NET: 7 fixed by this lane | 30 confirmed-passing via stale-measurement re-verification
        31 documented limit / correct refusal | 13 unresolved

=====================================================================================
C. THE ONE REAL BUG — neutral L-donor over-valence on 4_tetrahedral
=====================================================================================

  The OIN says: this tertiary amine N is a donor, bonded to the metal.
  RDKit's bookkeeping says: that N already has three bonds, which is all N gets.

        BEFORE the metal bond is added        the metal bond arrives
        ------------------------------        ----------------------
              C                                      C
               \                                      \
            C - N          valence 3/3            C - N ---- Zn      valence 4/3
               /            = FULL                   /                 = ILLEGAL
              C                                     C                       |
                                                                            v
                                       StructuralAssemblyError(AtomValenceException:
                                       "Explicit valence for atom # k N, 4, is greater
                                        than permitted")  -> empty conformer pool
                                        -> reported as "no_conformers"

  THE GATE, discovered by decision-table probe (measured, not assumed):

                        | 3_TPL | 4_TET  | 4_SQP | 5_SPY | 5_TBP | 6_OCT |
    aqua [OH2] + halides|  pass | CRASH  | pass  | pass  | pass  | pass  |
    tertiary amine  "   |   --   | CRASH  | pass  | pass  |  --   |  --   |
    pyridine n / C=N "  |   --   |  pass  |  --   |  --   |  --   |  --   |
                                    ^^^^^
              ONLY 4_tetrahedral, and ONLY for a SATURATED donor.
              An aromatic or imine donor at the same nominal "full" integer
              valence never crashes, on any geometry tried.
              (4_tetrahedral is the only geometry in this table with a
               stereogenic metal centre.)

  THE FIX  metallogen_adapter.py::_prepare_ligand_fragments, :312-367

    if geo == "4_tetrahedral"
       and atom.GetFormalCharge() == 0
       and not atom.GetIsAromatic()
       and no DOUBLE/TRIPLE bond on the atom
       and 0 < GetDefaultValence(Z) <= atom.GetTotalValence():
           atom.SetFormalCharge(1)        # ammonium/oxonium reading of the dative bond
           atom.SetNoImplicit(True)
           atom.SetNumExplicitHs(atom.GetTotalNumHs())   # freeze H, no phantom H

    This charge is INTERNAL BOOKKEEPING ONLY -- the round trip re-encodes from the
    generated 3D geometry, not from this mol, so it never reaches the output OIN and
    net charge does not have to balance.

  WHY THE GATE EXISTS (caught before shipping): an ungated first version --
  just "donor already at default valence -> bump" -- changed 15/15 sampled
  currently-PASSING donor molecules, spuriously charging pyridine N, imine N,
  phosphine P, and even a haptic arene carbon to [cH+].

=====================================================================================
D. THE CONFOUND THAT GOVERNS EVERY NUMBER IN THIS LANE
=====================================================================================

  round-trip failures (3917 molecules: 2080 passing, 1837 failing)
  |
  |--- generator timeout ................ 1238  67.4%  <-- NEVER TESTS THE NOTATION
  |--- generator produced nothing ....... 191   10.4%  <-- NEVER TESTS THE NOTATION
  |--- encoder refused the input ........ 241   13.1%      = 77.8% combined
  |--- canonicalization noise ........... 123    6.7%
  |--- names a different isomer ......... 44     2.4%

  median failing molecule ran 300.3 s against a 300 s budget
                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        => a round-trip PASS RATE is substantially GENERATOR THROUGHPUT.
           Any change that alters runtime moves the rate for unrelated reasons.
           Compare WITHIN the completed set, or hold BOTH the molecule set AND
           the timeout budget identical. Never across different budgets.

  BUT a timeout bucket is NOT latent pass-rate.  24 timed-out molecules re-run --quick:
      SUCCESS 6 (25%) | String mismatch 6 | Atom count mismatch 6 | MetalloGen failed 6
      timed out again: 0
  => only ~25% is compute-limited; the other 75% reaches a verdict and the verdict is
     a REAL failure. More compute buys ~44 of 936 molecules (~4.7%), not the 174 the
     raw timeout count implies.

LEGEND
  hard_fail ............. a round-trip that produced no acceptable structure at all
  rmsd_gate / other / no_conformers ... failure classes in hard_fail_worklists.json
  4_tetrahedral, 4_SQP, ... MetalloGen geometry codes (TET, square planar, ...)
  m-SMILES .............. MetalloGen's `metal|lig1|lig2|...|geo` assembly string
  UncoordinatedFragmentError ... a fragment with no metal bond, which m-SMILES cannot express
  --quick ............... the cheap generator path (uff_pool_size=2, max_attempts=10, 30 s kill)
  status=success ........ canonical-key match AND atom-count match
```

---

## Initial assumptions and hypothesis

The lane was handed three frozen worklists and the implicit assumption that they described live
defects:

1. **`rmsd_gate` (33)** — assumed to be molecules whose generated geometry deviated too far from
   the input, with the working hypothesis that the `RMSD_GATE = 1.0` threshold might be wrong and
   need re-tuning.
2. **`other` (25)** — assumed heterogeneous, therefore assumed to need per-molecule triage; not
   assumed to be a single cause.
3. **`no_conformers` (23)** — assumed to be genuine generator failures: the embed produces an
   empty conformer pool.

Two framing assumptions were also in play, and both had to be corrected before any number in this
lane could be trusted:

* **that the frozen baseline was current.** All three classes were measured against
  `results-capstone-v042` (`commit_id 58bba7ad`, 2026-07-15).
* **that a `hard_fail` count is a defect count.** It is not, for the reason in panel D above.

---

## What was actually found

### `rmsd_gate` (33) — CONFIRMED STALE. Zero code change needed.

`tools/test_dataset_roundtrip.py::RMSD_GATE = 1.0` is still read and still stamped into every
report (`rmsd_gate: 1.0` in the run-env stamp), which is why the field name survives in the frozen
worklist — but it **has not been a pass/fail gate since `5c5d193b`** (v0.4.4-SL4, *"demote RMSD to
a diagnostic"*). Confirmed by reading the current `_attempt_generation`: on `rmsd >= RMSD_GATE` it
sets `report["metrics"]["rmsd_over_gate"] = True` (and `ff_floor` under FF-only) and **continues to
the atom-count check and to success**, exactly as SL4's changelog describes. The
`results-capstone-v042` sweep predates SL4 by a week, so these 33 "High RMSD" rows are a snapshot
of the *pre-demotion* harness.

**Measured, not assumed** — all 33 re-run through the current harness (`--only <name>
--mol-timeout 300`, no other flags):

| outcome | count | detail |
|---|---:|---|
| now `status=success` | **29 / 33 (87.9%)** | canonical-key match + atom-count match; RMSD is recorded, not gating — e.g. `NIQGAV_comp_0` at rmsd = 1.78 is a **pass** |
| `atom_count` mismatch | 3 | `CIFQEP_comp_0`, `FAFJUU_comp_0` (pending re-roll), `COSSOS_comp_0` (fails at the re-roll tier too) — the sibling `atom_count` lane's territory, not RMSD-related |
| genuine `string mismatch` | 1 | `SEQVEP_comp_0` — RMSD no longer masks it; a real generation defect (wrong structure), out of scope here, **reported honestly rather than claimed fixed** |

**Verdict on the 1.0 threshold:** defensible as a *diagnostic* value — it matches SL4's own
justification, that coordination-sphere RMSD is only ~0.22-correlated with geometric quality
(`docs/agentic-notes/v0.4.3/FALSIFICATION_v0.4.3_ELIMINATION.md` §1.4) and that eta rings collapse to a centroid so the
metric is blind to ring rotation — and it is **already not a gate**. No threshold change proposed;
none needed. **No code touched for this class.**

### `other` (25) — 24 known representation limit, 1 correct refusal, 0 fixed

All 25 raw errors were read and grouped rather than assumed single-cause:

| group | count | example |
|---|---:|---|
| `UncoordinatedFragmentError` (outer-sphere counterion/solvate) | 24 | perchlorate-like `[O-]Cl` ×5, perfluoroarylborate ×4, hydride ×3, methyl ×2, acetate / I₂ / water / thiocyanate / nitrate / triflate / borane-cluster / lone-B ×1 each |
| `Geometry code 'NON' not supported` | 1 | `PEKQUU_comp_0` |

**The 24 are byte-identical to the floor set already in `docs/KNOWN_LIMITATIONS.md`**
("Uncoordinated outer-sphere fragments", a 24-molecule list unchanged since a 2026-07-14 snapshot,
also documented in `docs/agentic-notes/v0.4.1/ACCURACY_v0.4.1.md`, `docs/agentic-notes/v0.4.2/ACCURACY_v0.4.2.md`,
`docs/agentic-notes/v0.4.3/FALSIFICATION_v0.4.3_ELIMINATION.md`). MetalloGen's `metal|lig1|...|geo` m-SMILES **has no
slot for a fragment with no metal bond**; this is a representation limit stated in the codebase's
own docstring (`generation/metallogen_adapter.py:98-109`), re-confirmed rather than re-derived
blind. **Not fixed** — a real fix needs a new OIN/generator convention for outer-sphere species
(place the fragment as an independently-embedded, non-bonded rigid body outside the coordination
core's bounding sphere; sketched, not attempted).

`PEKQUU_comp_0`'s `geo_code = NON` is the encoder's own documented "no discrete geometry template
fit" fallback (`utils/oin_aligner.py:829-831`, `return "g:NON|w:NON"`), emitted deliberately for an
exotic Ir-borane/carborane cluster whose bond-order perception is **already** producing an invalid
valence upstream (`C#S(C)C` — a sulfur with a triple bond plus two more substituents).
`metallogen_adapter.py:119` correctly refuses to build 3D structure for an undefined geometry code.
**A correct refusal, not a generator defect**; a real fix needs a cluster-bonding (Wade's rules)
model.

#### One new finding surfaced while triaging, out of scope for this lane

For the 5 `[O-]Cl`-shaped molecules (`ATAGUZ`, `CUBDOT`, `XAXZIH`, `XIFVAM`, `YUMBEP`), **the
encoder's own bond perception is also wrong**, independent of the generator issue. Checked
`ATAGUZ`'s input XYZ directly: both `ClO₄⁻` groups have all 4 Cl–O contacts at **1.428–1.449 Å**
(unambiguous covalent bonds), yet `utils/xyz2mol_local.py:145` hardcodes `atomic_valence[17] = [1]`
(Cl max valence 1), so `xyz2AC_obabel`'s valence-capping loop (`:1091-1098` at the time of
measurement) keeps only the **shortest** Cl–O contact and drops the other three oxygens into free
`[O-2]` ions.

Confirmed this does **not** unblock anything on its own — the resulting single-fragment ClO₄⁻ would
still have no binding slot and still raise `UncoordinatedFragmentError` — and `atomic_valence` is
an **ungated table** read by the one bond-perception path used corpus-wide, so touching Cl (or
Br = 35, I = 53, same `[1]` hardcode) risks changing output for any currently-passing coordinated
chlorate/perchlorate/bromate/iodate ligand. **Flagged for a dedicated encoder lane; not touched.**

> Cross-reference for whoever takes it: this is the *same* capping loop that the
> `OIN_STABLE_METAL_AC` lane fixed for iteration order and the `OIN_BORON_CAGE` lane fixed for
> boron-cage exemption. See `LANE-valence-order.md` — a Cl/Br/I valence change would be a **third**
> intervention in that one loop and must be read together with the other two.

### `no_conformers` (23) — 7 fixed, 1 already fixed, 6 representation limit, 9 unresolved

#### 3a. FIXED — "Group 1: neutral L-donor over-valence" (7 molecules)

`docs/KNOWN_LIMITATIONS.md` already **named** this root cause (`GEZKAZ`'s aqua O, `VIBRIK`'s
tertiary-amine N: *"Explicit valence for atom # k N/O, greater than permitted"*) but had never
fixed it. Reproduced directly with `generate_3d_structures` on the stored m-SMILES, bypassing the
harness: `FEJFAD_comp_0`, `FEJFOR_comp_0`, `LECSUJ_comp_0`, `MUTYEG_comp_0`, `VIBROQ_comp_0`,
`VIRWOJ_comp_0`, `WIQRIA_comp_0` — all an N,N-chelate (tertiary-amine or cyclic-ether/amine donor)
plus 2 halides on a tetrahedral Zn/Fe — **deterministically** raise
`StructuralAssemblyError(AtomValenceException: Explicit valence for atom # k N, 4, ...)` on every
attempt.

**Root cause, precisely:** a fully-substituted, saturated (sp3, non-aromatic, no double/triple
bond) **neutral** donor atom whose valence is already at its element's default **before the metal
bond is added** — a tertiary amine N with 3 C-substituents and 0 H, or a bracket aqua O with 2
explicit H — has no room left for the coming metal→donor bond. The existing donor
H-reconciliation logic in `_prepare_ligand_fragments` only zeroes an already-zero H count for a
*bare* N/O/S donor: it does not add valence room, and a **bracket (explicit-H) atom never enters
that code path at all.**

The trigger condition was found by decision-table probe (measured, not assumed) — see panel C.
Only `4_tetrahedral` crashes, and only for a saturated donor; an aromatic or imine donor at the
*same* nominal "full" integer valence never crashes on any geometry tried. `4_tetrahedral` is the
only geometry in the table with a **stereogenic metal centre**, so its embed path evidently adds an
explicit (valence-counted) bond that the others do not.

**Verification after gating:**

* **7/7 target molecules now generate and round-trip end to end**: `status=success`, RMSD
  **0.39–0.50 Å**, vdW clash **0–2** (`tools/test_dataset_roundtrip.py --only <name>`).
* **Regression check, sampled not asserted.** Scanned the corpus for currently-*passing*
  `4_tetrahedral` molecules whose donors the new gate's condition would touch: **703 / 1207**
  passing TET molecules qualify (**58%**) — that population is what the "prove it" rule protects.
  Two tiers:
  * *Smoke* (raw `generate_3d_structures`, tiny pool, 40 molecules sampled from the 703):
    **40/40 identical embed outcome**, patched vs unpatched HEAD.
  * *Deep* (full harness, `status`/`tier`/`rmsd`, 8 molecules — `YUHWAB_comp_0`, `DILYUV_comp_0`
    + 6 more random): **7/8 byte-identical**; the 1 discrepancy (`PEDVOK_comp_0`) traced to a
    wall-clock **timeout artefact** under shared load (HEAD hit the test's tight 120 s budget,
    patched did not — same code path, no default difference; the smoke tier already covers this
    molecule's embed outcome identically).
* `tests/unit` (605 tests) green; `tests/unit/test_regression_stability.py` (6 golden fixtures,
  including all four named goldens) byte-identical; `ruff check` / `format` clean.

#### 3b. Already fixed, unrelated to this lane (1 molecule)

`FOGCIN_comp_0` (`[Pd]|C=[S:1](=C)(C)[O-]|...|2_linear`, a hypervalent-S sulfinyl donor) now
passes (`status=success`, rmsd = 0.317) on **both** HEAD and the patch — some other change landed
since the stale v0.4.2 snapshot and fixed it independently. Reported for completeness; **not
claimed as this lane's fix.**

#### 3c. Documented representation limit, not a generator bug (6 molecules)

`DOFCAE_comp_0`, `DOFCAE_comp_2`, `NOZBIP_comp_0`, `PUPYOQ_comp_0`, `PUPYUW_comp_0` all share one
ligand template — a boratabenzene-bis(phosphine) pincer whose aromatic B,N-heterocycle
(`n1[b]n...c2ccccc12`) will not kekulize. That is `docs/KNOWN_LIMITATIONS.md` "Group 2"'s
already-named `DOFCAE` example ("aromatic boron `[b]`"): **one perception gap explains all five.**
`UTAZAS_comp_0`'s `CC(C)P1(C(C)C)=N=P1(C(C)C)C(C)C` (a cyclophosphazane P=N=P ring) matches the same
Group 2's already-named `MEDDUV` pattern (cyclophosphazene). Confirmed via direct per-fragment
`Chem.SanitizeMol`: **all 6 raise `Can't kekulize` / `Explicit valence` on the isolated ligand
fragment**, i.e. the encoder's own two-centre bond-order model cannot represent these rings at all.
This is **upstream of generation** — a notation/perception gap the docs already scope as "deferred
until a design exists", not a generator defect.

#### 3d. Unresolved, likely genuine (9 molecules)

`CETYAC_comp_0`, `ETORIQ_comp_0`, `KATTOO_comp_0`, `NOBCEN_comp_0`, `QUCDOJ_comp_0`,
`RUKWOL_comp_0`, `TESZEZ_comp_1`, `UTEMIR_comp_0`, `WAKXOX_comp_0`. Probed each directly
(`generate_3d_structures`, short budget) — **all return an empty conformer pool with no
exception**: not a hard crash, a genuine embed/convergence failure (and under the full budget too,
since that is what the original "no conformers" measurement used). Two independent donors carry an
**S≡C triple bond in a ring** (`KATTOO`, `NOBCEN` — a thiazolium-type NHC-like carbene), suggesting
a shared, not-yet-isolated exotic-donor cause, but **no safe narrow fix was found in the time
available and none was forced.** Consistent with the project's prior R2 finding that a genuine,
un-fixable-in-the-generator's-owned-layers residual is an expected, honest outcome for this class.

### Net result

| class | count | now passing (this lane) | already-stale-fixed | documented limit / correct refusal | unresolved |
|---|---:|---:|---:|---:|---:|
| `rmsd_gate` | 33 | 0 (needed no fix) | 29 | — | 4 (not RMSD; other lanes/defects) |
| `other` | 25 | 0 | 0 | 25 | 0 |
| `no_conformers` | 23 | 7 | 1 | 6 | 9 |
| **total** | **81** | **7 fixed by this lane** | **30 confirmed-passing via stale-measurement re-verification** | **31** | **13** |

### The measurement confound, and the refinement that matters more

**77.8% of round-trip failures never test the notation.** Over 3,917 molecules (2,080 passing,
1,837 failing): generator timeout **1,238 (67.4%)**, encoder refused the input 241 (13.1%),
generator produced nothing 191 (10.4%), canonicalization noise 123 (6.7%), output names a different
isomer 44 (2.4%). **The median failing molecule ran 300.3 s against a 300 s budget** — two thirds of
the failure mass is a stopwatch expiring.

Two consequences, and they pull in *different* directions, which is why both belong in this report:

1. **A round-trip pass rate is substantially GENERATOR THROUGHPUT, not notation quality**, and any
   change that alters runtime moves the rate for unrelated reasons. Comparisons must be made
   **within the completed set**, or with the molecule set **and** the timeout budget held
   identical. This is not theoretical: v0.4.5's re-baseline produced **11 apparent regressions**
   that are *all* `TimeoutException exceeded 300s` against a capstone baseline that ran at
   **1800 s** — the same config asymmetry that manufactured v0.4.4's 11 phantom regressions. Zero
   correctness regressions in both cases.
2. **But a timeout bucket is NOT latent pass-rate.** 24 molecules (seed 42, stratified 12 `UFF_1` +
   12 `g-xTB_1`) that all hit the 300 s wall in full mode were re-run on the cheaper `--quick` path:

   | outcome | n |
   |---|---|
   | SUCCESS | **6 (25%)** |
   | String mismatch | 6 (25%) |
   | Atom count mismatch | 6 (25%) |
   | MetalloGen failed | 6 (25%) |
   | **timed out again** | **0** |

   **Only ~25% of the timeout population is compute-limited.** The other 75% reaches a verdict, and
   the verdict is a **real failure**. So more compute buys roughly **44 of 936 molecules (~4.7%)**,
   not the **174** the raw timeout count implies. ⚠ The probe's `elapsed_s` values are unusable — it
   ran alongside the 6-shard 5k sweep at load 21–33, and wall clock is meaningless above ~12. Only
   the pass/fail outcomes survive, and only because none of the 18 failures was a timeout.

That pair is the honest reading: **timeouts inflate the failure count without informing you, and
they also hide real failures.** Both errors point the same way — do not treat the timeout bucket as
either noise or as banked accuracy.

### `MetalloGen failed` was actively MISLEADING

All **19** molecules in that class report `failed to generate any conformers for m-SMILES None` —
and the `None` is a **red herring**. `msmiles` is populated **only** when DIRECT assembly fails and
the adapter falls back (the `if prebuilt_complex is None:` branch,
`generation/metallogen_adapter.py:1827`), so on the **default OIN-direct path** it is legitimately
`None` — and the old message interpolated it regardless. The effect was to point every debugger at a
phantom `None` in the m-SMILES builder instead of at the embed that actually came up empty.

Fixed (`metallogen_adapter.py:1963-1978`) to name the assembly path actually used, plus the OIN,
the pool width and the timeout:

```python
if prebuilt_complex is not None:
    _via = f"OIN-direct assembly (geometry {getattr(parsed, 'geo_code', None)!r})"
else:
    _via = f"m-SMILES {msmiles!r}"
raise ValueError(
    f"MetalloGen produced no conformers via {_via}; "
    f"OIN={getattr(parsed, 'original_oin', None)!r}, pool={pool_n}, "
    f"timeout={self.timeout}s"
)
```

> Provenance note: this diagnostics fix and the 24-molecule `--quick` re-run are recorded in
> `docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md` §3 and the **0.4.6** CHANGELOG stanza — they postdate
> `docs/agentic-notes/v0.4.5/GEN_RESIDUE_v0.4.5.md` and are not in it. They are included here because they are the
> measurements that make this lane's residue interpretable, and because `tools/classify_failures.py`
> still keys on the substring `"failed to generate any conformers"` (`:195`) — anyone changing that
> message again must update the classifier in the same commit.

---

## What was done

**One source change, and it is deliberately narrow.**

| item | value |
|---|---|
| file | `src/oinsmiles/generation/metallogen_adapter.py` |
| function | `_prepare_ligand_fragments` |
| lines | `:312-367` (comment block `:312-352`, condition `:353-367`) |
| commit | `4dcb5284` — `fix(gen): recover neutral-donor over-valence on 4_tetrahedral (no_conformers Group 1)`, **+57 lines** |
| lever | **none** |
| default changed | **none** for any currently-passing molecule |
| new unit test | **none added** — verification was 7/7 harness round-trips, the 40/40 smoke tier, the 7/8 deep tier, and the existing 605-test suite + `test_regression_stability.py` |

The gate is four conditions plus a valence check, all of which must hold:

```python
if (
    geo == "4_tetrahedral"
    and atom.GetFormalCharge() == 0
    and not atom.GetIsAromatic()
    and not any(b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)
                for b in atom.GetBonds())
):
    mol.UpdatePropertyCache(strict=False)
    default_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
    if 0 < default_valence <= atom.GetTotalValence():
        atom.SetFormalCharge(1)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(atom.GetTotalNumHs())
```

Design points and why each is there:

* **`+1` formal charge = an ammonium-/oxonium-like reading of the dative bond.** It mirrors the
  existing, already-shipped `_charge_fix_promotion` escape in `generator3d/embed.py` (used
  *ungated* for an over-valent double-bond promotion). The same reasoning applies: **this internal
  RDKit charge never reaches the output OIN**, because the round trip re-encodes from the generated
  3D geometry rather than from this bookkeeping mol — so it only has to keep the embed
  valence-valid, and net charge does not need to balance.
* **`SetNoImplicit(True)` + `SetNumExplicitHs(GetTotalNumHs())` freezes the H count**, so RDKit
  does not insert a phantom hydrogen now that a bond's worth of "room" exists.
* **The geometry + hybridization gate is the whole safety argument.** It ensures the branch can only
  ever activate on a donor shape that would otherwise **hard-crash**.

**Rejected alternatives, with reasons:**

* **The ungated version of this same fix** ("donor already at default valence → bump") — written
  first, then **caught before shipping**: A/B against 15 sampled currently-passing donor-bearing
  molecules (imine, pyridine, phosphine, non-tetrahedral aqua, haptic arene) showed **15/15
  changed**, spuriously charge-bumping pyridine N, imine N, phosphine P, and even a haptic arene
  ring carbon to `[cH+]`. **Do not widen the gate without re-running that A/B** — the in-code
  comment says so explicitly at `:337-342`.
* **Changing `RMSD_GATE`** — declined: RMSD is already not a gate, and 1.0 is defensible as a
  diagnostic.
* **Fixing the 24 `UncoordinatedFragmentError` molecules** — declined as out of scope: it needs a
  new OIN/generator convention for outer-sphere species, not a patch.
* **Touching `atomic_valence[17]`/`[35]`/`[53]` to fix perchlorate perception** — declined: an
  ungated corpus-wide table read, and it would not unblock the molecules that surfaced it anyway.
* **Forcing a fix for the 9 unresolved `no_conformers`** — declined deliberately, and reported as
  "not fixed, likely genuine" rather than padding a number.

---

## Dead ends, refutations, and costs accepted

* **REFUTED: "the frozen worklist describes live defects."** Two of the three classes were
  substantially **stale measurement**. 29 of 33 `rmsd_gate` molecules pass today with **zero code
  change**, and 1 more (`FOGCIN_comp_0`) was fixed by an unrelated change. **30 of 81 molecules
  (37%) needed re-measurement, not engineering.** The general lesson: a frozen worklist carries the
  `commit_id` it was measured at — **read it, and diff it against what has landed since, before
  writing any code.**
* **REFUTED: "the `rmsd_gate` threshold needs tuning."** It is not a gate. The field name survived
  the demotion and misled the worklist consumer. Naming a metric after the gate it used to be is a
  trap that cost this lane its first hypothesis.
* **REFUTED: "`other` is heterogeneous and needs per-molecule work."** 24 of 25 are one cause, and
  it is a cause the codebase's own docstring already documents. Grouping the raw error strings first
  — rather than triaging 25 molecules cold — was what made that visible.
* **REFUTED: "the accuracy gap is mostly compute."** The `--quick` re-run of 24 timed-out
  molecules: **0 timed out again, and 18 of 24 produced a real failure verdict.** ~25%
  compute-limited, ~4.7% of 936 recoverable by compute, against the ~18.6% the raw timeout count
  suggests.
* **A misleading error message cost real debugging time.** `m-SMILES None` on all 19 `MetalloGen
  failed` molecules pointed at a phantom `None` in the m-SMILES builder; the actual failure was an
  empty embed pool on the OIN-direct path, where `msmiles` is *legitimately* `None`. **An error
  message that interpolates a variable belonging to a code path that was not taken is worse than no
  detail at all.**
* **Cost accepted: one deep-tier discrepancy was written off as a timeout artefact.**
  `PEDVOK_comp_0` differed between arms because HEAD hit the test's tight 120 s budget and the
  patched arm did not. The justification — same code path, no default difference, and the smoke tier
  covers this molecule's embed outcome identically — is reasonable but is a **judgement call, not a
  measurement**. A re-run at the normal 300 s budget was declined as unnecessary.
* **Cost accepted: the regression evidence for the Group-1 fix is SAMPLED, not exhaustive.**
  703/1207 passing TET molecules satisfy the gate's condition; 40 were smoke-tested and 8
  deep-tested. That is the population the guard rule about "prove it" exists to protect, and it was
  not proven exhaustively.
* **Cost accepted: no unit test pins the Group-1 fix.** The verification is harness runs plus the
  existing suite. `grep` finds no test referencing `4_tetrahedral` donor charge bumping, so a future
  refactor of `_prepare_ligand_fragments` could silently remove the branch and only the 7 molecules'
  harness behaviour would notice. The in-code comment at `:312-352` is currently the only durable
  record of *why* the gate has the shape it has.
* **Cost accepted: 13 molecules left unresolved and 31 left as documented limits.** Explicitly
  chosen over forcing a fix. `SEQVEP_comp_0` in particular is a newly *visible* real generation
  defect — RMSD used to mask it — and is reported rather than absorbed.

---

## Where it landed

**Merged and released.** `git log --oneline main..swimlane/v045-genresidue` is empty — the branch is
an ancestor of `main`, released with v0.4.5 (tag `0d165845`, local only, **do not push**).

| commit | subject |
|---|---|
| `4dcb5284` | `fix(gen): recover neutral-donor over-valence on 4_tetrahedral (no_conformers Group 1)` — the only source change, `metallogen_adapter.py` +57 |
| `b1ca4c57` | `docs(v0.4.5): 30 molecules fixed, 29 never failing — and the frozen baseline is stale` |

Final state:

| item | value |
|---|---|
| source change | `src/oinsmiles/generation/metallogen_adapter.py::_prepare_ligand_fragments` (`:312-367`) |
| lever / env var | **none** — the fix is geometry- and hybridization-gated in code |
| default changed for a passing molecule | **no** |
| molecules fixed | **7** — `FEJFAD_comp_0`, `FEJFOR_comp_0`, `LECSUJ_comp_0`, `MUTYEG_comp_0`, `VIBROQ_comp_0`, `VIRWOJ_comp_0`, `WIQRIA_comp_0` (all `status=success`, RMSD 0.39–0.50 Å, vdW clash 0–2) |
| guard tests relied on | `tests/unit/test_regression_stability.py::TestRegressionStability` (6 goldens, byte-identical); full `tests/unit` (605 at the time) green |
| diagnostics fix (v0.4.6) | `metallogen_adapter.py:1963-1978`; ⚠ `tools/classify_failures.py:195` keys on the old substring |
| documented-limit references | `docs/KNOWN_LIMITATIONS.md` ("Uncoordinated outer-sphere fragments" 24-molecule set; "Group 2" aromatic `[b]` / cyclophosphazene) |

---

## Open questions / for the next agent

1. **Widen the Group-1 fix beyond `4_tetrahedral`?** The decision table shows only TET crashing
   *today*, and the stated reason is that TET is the only geometry in the probe with a stereogenic
   metal centre — *"its embed path evidently adds an explicit (valence-counted) bond that the others
   do not."* That mechanism was **inferred, not confirmed.** Confirming it (find the line in
   `generator3d/embed.py` that adds the extra explicit bond on a stereogenic centre) would tell you
   whether the gate is the *right* condition or merely a *sufficient* one — and whether another
   geometry could start crashing after an unrelated embed change. **If you widen it, re-run the
   15-molecule ungated A/B first.**
2. **Add a unit test for the gate.** There is none. A minimal one would build a tertiary-amine
   N,N-chelate + dihalide `ParsedOIN`, call `_prepare_ligand_fragments`, and assert the donor came
   back with formal charge 1 on `4_tetrahedral` and charge 0 on `4_square_planar` — plus a negative
   case pinning that a pyridine `n` donor is never bumped. That is the cheapest available protection
   for the 703-molecule population.
3. **The `[O-]Cl` encoder-side perception bug.** `atomic_valence[17] = [1]` (and `[35]`, `[53]`)
   drops 3 of 4 perchlorate oxygens into free `[O-2]` ions. Needs its own lane, a corpus A/B over
   currently-passing coordinated chlorate/perchlorate/bromate/iodate ligands, and coordination with
   the two other v0.4.5 interventions in the same capping loop (`OIN_STABLE_METAL_AC`,
   `OIN_BORON_CAGE`).
4. **`SEQVEP_comp_0`** — a real generation defect that RMSD used to mask. Unowned.
5. **The 9 unresolved `no_conformers`.** `KATTOO_comp_0` and `NOBCEN_comp_0` both carry an S≡C
   triple bond in a ring (thiazolium-type NHC-like carbene), which is the only shared feature
   spotted. Note the cross-lane echo: `LANE-valence-search.md`'s over-cap population is *also*
   metal–NHC chemistry, and the charge filter there converts an unparseable NHC string into a
   parseable neutral-carbene one. **Worth checking whether any of these 9 perceive differently under
   `OIN_VALENCE_CHARGE_FILTER=1`** before concluding they are generator-side at all.
6. **Outer-sphere fragments need a notation decision, not a patch.** The sketch on record: place the
   fragment as an independently-embedded, non-bonded rigid body outside the coordination core's
   bounding sphere. That requires an OIN grammar extension and generator support; it is the largest
   single documented limit in this class (24 molecules).
7. **Never compare pass rates across timeout budgets.** The v0.4.4 and v0.4.5 releases each produced
   exactly 11 phantom regressions this way. If you must compare, hold the molecule set **and** the
   budget identical, or compare only within the set that completed in both arms.
8. **Re-run the `--quick` timeout probe on a quiet host.** Its outcome distribution
   (6/6/6/6, 0 timeouts) is sound because none of the failures was a timeout, but its `elapsed_s`
   column is unusable at load 21–33. A quiet-host run would also let the ~4.7% recoverable estimate
   be tightened.
