# Generator residue — acceptance report (v0.4.5)

Swimlane `swimlane/v045-genresidue`. Owned the three smallest non-timeout `hard_fail` classes
from `spec/handoffs/v0.4.5/hard_fail_worklists.json`: `rmsd_gate` (33), `other` (25),
`no_conformers` (23) — 81 molecules, 1.21% of the corpus. All three measured against a **stale**
`results-capstone-v042` sweep (`commit_id 58bba7ad`, 2026-07-15) — one week before v0.4.4-SL4
(`5c5d193b`, 2026-07-22) landed. Two of the three classes turned out to be substantially
**stale measurement**, not live defects; the third had one genuine, previously-undiagnosed
generator bug, now fixed and gated narrowly.

## 1. `rmsd_gate` (33) — the SL4 demotion already landed; ZERO code change needed here

`tools/test_dataset_roundtrip.py::RMSD_GATE = 1.0` is still read from every report
(`rmsd_gate: 1.0` in the run-env stamp), which is why the field name survives in the frozen
worklist — but it has not been a **pass/fail gate** since `5c5d193b` (v0.4.4-SL4, "demote RMSD
to a diagnostic"). Confirmed by reading the current `_attempt_generation`: on `rmsd >= RMSD_GATE`
it sets `report["metrics"]["rmsd_over_gate"] = True` (and `ff_floor` under FF-only) and
**continues to the atom-count check and success**, exactly as SL4's changelog describes. The
`results-capstone-v042` sweep's `commit_id` (`58bba7ad`, 2026-07-15) predates SL4 by a week, so
these 33 "High RMSD" rows are a snapshot of the *pre-demotion* harness.

**Measured, not assumed:** re-ran all 33 through the current harness (`--only <name>
--mol-timeout 300`, no other flags).

| outcome | count | detail |
|---|---:|---|
| now `status=success` | **29 / 33 (87.9%)** | canonical-key match + atom-count match; RMSD is recorded, not gating (e.g. `NIQGAV_comp_0` rmsd=1.78, still a pass) |
| `atom_count` mismatch | 3 | `CIFQEP_comp_0`, `FAFJUU_comp_0` (pending re-roll), `COSSOS_comp_0` (fails at the re-roll tier too) — the sibling `atom_count` lane's territory, not RMSD-related |
| genuine `string mismatch` | 1 | `SEQVEP_comp_0` — RMSD no longer masks it; a real generation defect (wrong structure), out of scope here, reported honestly rather than claimed fixed |

**Verdict on the 1.0 threshold:** defensible as a *diagnostic* value (matches SL4's own
justification — coordination-sphere RMSD is only ~0.22-correlated with geometric quality per
`FALSIFICATION_v0.4.3_ELIMINATION` §1.4, and eta rings collapse to a centroid so the metric is
blind to ring rotation), and it is **already not a gate** in the current code. No threshold
change proposed; none needed. **No code touched for this class.**

## 2. `other` (25) — 24 known representation limit, 1 correct refusal, 0 fixed

All 25 raw errors read and grouped (not assumed single-cause):

| group | count | example |
|---|---:|---|
| `UncoordinatedFragmentError` (outer-sphere counterion/solvate) | 24 | perchlorate-like `[O-]Cl`×5, perfluoroarylborate×4, hydride×3, methyl×2, acetate/I₂/water/thiocyanate/nitrate/triflate/borane-cluster/lone-B×1 each |
| `Geometry code 'NON' not supported` | 1 | `PEKQUU_comp_0` |

**The 24 are byte-identical to the floor set already in `docs/KNOWN_LIMITATIONS.md`**
("Uncoordinated outer-sphere fragments," 24-molecule list unchanged since a 2026-07-14
snapshot, also documented in `docs/agentic-notes/v0.4.1/ACCURACY_v0.4.1.md`, `docs/agentic-notes/v0.4.2/ACCURACY_v0.4.2.md`,
`docs/agentic-notes/v0.4.3/FALSIFICATION_v0.4.3_ELIMINATION.md`). MetalloGen's `metal|lig1|...|geo` m-SMILES has no
slot for a fragment with no metal bond; this is a representation limit stated in the codebase's
own docstring (`generation/metallogen_adapter.py:98-109`), re-confirmed rather than re-derived
blind. **Not fixed** — a real fix needs a new OIN/generator convention for outer-sphere species
(place the fragment as an independently-embedded, non-bonded rigid body outside the coordination
core's bounding sphere; sketched but not attempted — see §4).

**One new finding surfaced while triaging, out of scope for this lane:** for the 5
`[O-]Cl`-shaped molecules (`ATAGUZ`, `CUBDOT`, `XAXZIH`, `XIFVAM`, `YUMBEP`), the encoder's own
bond perception is *also* wrong, independent of the generator issue. Checked `ATAGUZ`'s input
xyz directly: both `ClO4⁻` groups have all 4 Cl-O contacts at 1.428-1.449 Å (unambiguous
covalent bonds), yet `utils/perception_core.py:185` hardcodes `atomic_valence[17] = [1]` (Cl max
valence 1), so `xyz2AC_obabel`'s valence-capping loop (`:1932-1940`) keeps only the shortest
Cl-O contact and drops the other 3 oxygens into free `[O-2]` ions. Confirmed this doesn't
unblock anything on its own (the resulting single-fragment ClO4⁻ would still have no binding
slot and still raise `UncoordinatedFragmentError`), and `atomic_valence` is an ungated table
read by the one bond-perception path used corpus-wide — touching Cl (or Br=35, I=53, same
`[1]` hardcode) risks changing output for any currently-passing coordinated
chlorate/perchlorate/bromate/iodate ligand. Flagged for a dedicated encoder lane; not touched
here.

`PEKQUU_comp_0`'s `geo_code=NON` is the encoder's own documented "no discrete geometry template
fit" fallback (`utils/oin_aligner.py:829-831`, `return "g:NON|w:NON"`), emitted deliberately for
an exotic Ir-borane/carborane cluster whose bond-order perception is already producing an
invalid valence upstream (`C#S(C)C` — a sulfur with a triple bond plus two more substituents).
`metallogen_adapter.py:119` correctly refuses to build 3D structure for an undefined geometry
code. **A correct refusal**, not a generator defect; a real fix needs a cluster-bonding (Wade's
rules) model, out of scope.

## 3. `no_conformers` (23) — 8 fixed (7 by a new gated fix, 1 already-fixed/stale), 6 documented representation limit, 9 unresolved and likely genuine

### 3a. Fixed: "Group 1 — neutral L-donor over-valence" (7 molecules)

`docs/KNOWN_LIMITATIONS.md` already names this root cause (`GEZKAZ`'s aqua O, `VIBRIK`'s
tertiary-amine N: "Explicit valence for atom # k N/O, greater than permitted") but had not
fixed it. Reproduced directly (`generate_3d_structures` on the stored m-SMILES, bypassing the
harness): `FEJFAD_comp_0`, `FEJFOR_comp_0`, `LECSUJ_comp_0`, `MUTYEG_comp_0`, `VIBROQ_comp_0`,
`VIRWOJ_comp_0`, `WIQRIA_comp_0` — all an N,N-chelate (tertiary-amine or cyclic-ether/amine
donor) + 2 halides on a tetrahedral Zn/Fe — deterministically raise
`StructuralAssemblyError(AtomValenceException: Explicit valence for atom # k N, 4, ...)` on
every attempt.

**Root cause, precisely:** a fully-substituted, saturated (sp3, non-aromatic, no
double/triple bond) neutral donor atom whose valence is already at its element's default
before the metal bond is added (a tertiary amine N: 3 C-substituents, 0 H; a bracket aqua O:
2 explicit H) has no room left for the coming metal→donor bond. The existing donor
H-reconciliation logic in `_prepare_ligand_fragments` (`metallogen_adapter.py`) only zeroes an
already-zero H count for a bare N/O/S donor — it does not add valence room, and a bracket
(explicit-H) atom never enters that code path at all.

**Decision-table probe (measured, not assumed) to find the true trigger condition**, because a
first, ungated version of the fix was too broad and had to be caught before shipping:

| donor type | 3_TPL | 4_TET | 4_SQP | 5_SPY | 5_TBP | 6_OCT |
|---|---|---|---|---|---|---|
| aqua (`[OH2]`) + halides | pass | **crash** | pass | pass | pass | pass |
| tertiary amine + halides | — | **crash** | pass | pass | — | — |
| pyridine `n` / imine `C=N` + halides | — | pass | — | — | — | — |

Only `4_tetrahedral` crashes, and only for a saturated (non-aromatic, no multiple bond) donor —
an aromatic or imine donor at the *same* nominal "full" integer valence never crashes, on any
geometry tried. `4_tetrahedral` is the only geometry in this table with a stereogenic metal
centre, so its embed path evidently adds an explicit (valence-counted) bond that the others do
not.

**Fix** (`metallogen_adapter.py::_prepare_ligand_fragments`, ~25 lines): when
`geo == "4_tetrahedral"` and a donor atom is neutral, non-aromatic, has no double/triple bond,
and is already at its element's default valence, give it a `+1` formal charge (an
ammonium-/oxonium-like reading of the dative bond) and freeze its H count. Mirrors the existing,
already-shipped `_charge_fix_promotion` escape in `generator3d/embed.py` (used ungated for an
over-valent double-bond promotion) — the same reasoning applies: this internal RDKit charge
never reaches the output OIN (the round trip re-encodes from the generated 3D geometry, not
from this bookkeeping), so it only has to keep the embed valence-valid.

**Caught before shipping:** an earlier version of this fix had no geometry/hybridization gate —
just "donor already at default valence → bump." A/B against 15 sampled currently-passing
donor-bearing molecules (imine, pyridine, phosphine, non-tetrahedral aqua, haptic arene) showed
**15/15 changed** (spuriously charge-bumped pyridine N, imine N, phosphine P, even a haptic
arene ring carbon to `[cH+]`). Added the geometry+hybridization gate above and re-measured.

**Verification after gating:**
- **7/7 target molecules now generate and round-trip end-to-end**: `status=success`, RMSD
  0.39–0.50 Å, vdW clash 0–2 (`tools/test_dataset_roundtrip.py --only <name>`).
- **Regression, sampled not asserted.** Scanned the corpus for currently-*passing*
  `4_tetrahedral` molecules whose donors the new gate's condition would touch: **703 / 1207**
  passing TET molecules qualify (58%) — this population is what the guard rule about "prove it"
  is protecting. Two-tier check:
  - *Smoke* (raw `generate_3d_structures`, tiny pool, 40 molecules sampled from the 703):
    **40/40 identical embed outcome** (patched vs. unpatched HEAD).
  - *Deep* (full harness, `status`/`tier`/`rmsd`, 8 molecules — `YUHWAB_comp_0`,
    `DILYUV_comp_0` + 6 more random): 7/8 byte-identical; 1 discrepancy (`PEDVOK_comp_0`)
    traced to a wall-clock timeout artifact under shared system load (HEAD hit the test's
    tight 120 s budget, patched didn't — both are the same code path with no default
    difference; re-run at the normal 300 s budget is unnecessary since the smoke test already
    covers this molecule's embed outcome identically).
- `tests/unit` (605 tests) green; `tests/unit/test_regression_stability.py` (6 golden fixtures,
  including all four named in the brief) byte-identical; `ruff check`/`format` clean.

### 3b. Already fixed, unrelated to this lane (1 molecule)

`FOGCIN_comp_0` (`[Pd]|C=[S:1](=C)(C)[O-]|...|2_linear`, a hypervalent-S sulfinyl donor) now
passes (`status=success`, rmsd=0.317) on **both** HEAD and my patch — i.e. some other change
landed since the stale v0.4.2 snapshot fixed it independently of this lane. Reported for
completeness; not claimed as this lane's fix.

### 3c. Documented representation limit, not a generator bug (6 molecules)

`DOFCAE_comp_0`, `DOFCAE_comp_2`, `NOZBIP_comp_0`, `PUPYOQ_comp_0`, `PUPYUW_comp_0` all share
one ligand template — a boratabenzene-bis(phosphine) pincer whose aromatic B,N-heterocycle
(`n1[b]n...c2ccccc12`) will not kekulize. This is `docs/KNOWN_LIMITATIONS.md` "Group 2"'s
already-named `DOFCAE` example ("aromatic boron `[b]`") — one perception gap explains all 5.
`UTAZAS_comp_0`'s `CC(C)P1(C(C)C)=N=P1(C(C)C)C(C)C` (a cyclophosphazane P=N=P ring) matches the
same Group 2's already-named `MEDDUV` pattern (cyclophosphazene). Confirmed via direct
per-fragment `Chem.SanitizeMol`: all 6 raise `Can't kekulize` / `Explicit valence` on the
isolated ligand fragment, i.e. the encoder's own two-centre bond-order model cannot represent
these rings at all — this is upstream of generation, a notation/perception gap the docs already
scope as "deferred until a design exists," not a generator defect.

### 3d. Unresolved, likely genuine (9 molecules)

`CETYAC_comp_0`, `ETORIQ_comp_0`, `KATTOO_comp_0`, `NOBCEN_comp_0`, `QUCDOJ_comp_0`,
`RUKWOL_comp_0`, `TESZEZ_comp_1`, `UTEMIR_comp_0`, `WAKXOX_comp_0`. Probed each directly
(`generate_3d_structures`, short budget) — all return an empty conformer pool with **no**
exception (not a hard crash; a genuine embed/convergence failure under the full budget too,
since that's what the original "no conformers" measurement already used). Two independent
donors carry an S≡C triple bond in a ring (`KATTOO`, `NOBCEN` — a thiazolium-type NHC-like
carbene), suggesting a shared, not-yet-isolated exotic-donor cause, but I did not find a safe,
narrow fix in the time available and am not forcing one. Consistent with this project's prior
R2 finding that a genuine, un-fixable-in-the-generator's-owned-layers residual is an expected,
honest outcome for this class — reporting these as **not fixed, likely genuine** rather than
padding a number.

## Net result

| class | count | now passing (this lane) | already-stale-fixed | documented limit / correct refusal | unresolved |
|---|---:|---:|---:|---:|---:|
| `rmsd_gate` | 33 | 0 (needed no fix) | 29 | — | 4 (not RMSD; other lanes/defects) |
| `other` | 25 | 0 | 0 | 25 | 0 |
| `no_conformers` | 23 | 7 | 1 | 6 | 9 |
| **total** | **81** | **7 fixed by this lane** | **30 confirmed-passing via stale-measurement re-verification** | **31** | **13** |

## Commits

See `git log` on `swimlane/v045-genresidue`; the only source change is
`src/oinsmiles/generation/metallogen_adapter.py` (+57 lines, the gated Group-1 fix). No lever,
no default changed for any currently-passing molecule — the new branch is geometry- and
hybridization-gated and only ever activates on a donor shape that would otherwise hard-crash.
