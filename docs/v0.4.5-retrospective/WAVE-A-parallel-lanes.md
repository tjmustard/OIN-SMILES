# Wave A — the parallel lanes (1, 3, 4, 7) and the eight-plus lanes nobody planned

**Purpose of the wave:** run every v0.4.5 lane that has *no* dependency on another lane
concurrently, in separate git worktrees, so the release's critical path is only as long as its one
genuine chain (Lane 1 → Lane 2 → Lanes 5/6). Lane 1 was the blocking lane; Lanes 3, 4 and 7 were
scheduled alongside it because none of them reads another lane's output.

---

## ELI5

An OIN string is a one-line text spelling of a 3D metal complex. v0.4.5's job was to make that
spelling *canonical* — the same complex must always produce the same text, no matter how the input
file happened to be rotated or how its atoms happened to be numbered. Four separate defects were
known to break that, and none of them touched the others, so four agents worked on them at the same
time in four private copies of the repository: the ligand *body* text (Lane 1), the ring-winding
`>`/`<` character (Lane 3), the twisted-biaryl "axial" descriptor (Lane 4), and the measurement
fixtures and oracles the later lanes needed (Lane 7). Every change shipped behind an environment
switch that defaulted to *off*, so nothing in the wave could change what users got until a single
measured decision later in the release. What actually happened is that three of the four lanes
found their own stated target was wrong, and ten *more* lanes were opened along the way — one of
which found that 13% of molecules were silently emitting a different stereochemistry depending on
the input file's atom order.

## The wave, visually

```
                        WAVE W0 — build and TRUST the instrument
                     (tools/canonicality_probe.py, reencode_ab.py,
                      build_sweep_cohort.py, run_sweep.sh)   [WAVE-W0-instrument.md]
                                        │
        ┌───────────────┬───────────────┼───────────────┬────────────────┐
        │               │               │               │                │
   ┌────▼────┐    ┌─────▼────┐    ┌─────▼────┐    ┌─────▼────┐           │
   │ LANE 1  │    │  LANE 3  │    │  LANE 4  │    │  LANE 7  │           │
   │ canon.  │    │ winding  │    │  axial   │    │ research │           │
   │  body   │    │ residual │    │  P2      │    │residuals │           │
   ├─────────┤    ├──────────┤    ├──────────┤    ├──────────┤           │
   │target   │    │target    │    │target    │    │target    │           │
   │rdkit_   │    │winding_  │    │generator │    │close Y3  │           │
   │canonical│    │star_drift│    │multi-axis│    │residue + │           │
   │500 → ~0 │    │  13 → 0  │    │  case    │    │build the │           │
   ├─────────┤    ├──────────┤    ├──────────┤    │2 missing │           │
   │◆ TARGET │    │◆ TARGET  │    │◆ CHARTER │    │fixtures  │           │
   │REFUTED  │    │ STALE +  │    │ PREMISE  │    ├──────────┤           │
   │396/500  │    │REFUSED   │    │ REFUTED  │    │✔ 2 fixt. │           │
   │are Lane │    │13→6 real,│    │generator │    │✔ torsion │           │
   │2's      │    │6→2 not 0 │    │HOLDS both│    │  oracle  │           │
   │✔ 2 levs │    │✔ 1 lever │    │axes      │    │✔ twin ops│           │
   │  built  │    │ + 1 uncd.│    │✚ DEVIATED│    │✖ NO      │           │
   │         │    │  fix     │    │from own  │    │ encoder  │           │
   │         │    │          │    │handoff   │    │ change   │           │
   └────┬────┘    └──────────┘    └──────────┘    └────┬─────┘           │
        │                                              │                 │
        │ BLOCKS (canonical body = the vertex COLOUR)   │ unblocks        │
        ▼                                              ▼                 │
   ══ WAVE B: LANE 2 canonical slots ══           ══ WAVE C: LANE 5 ══   │
                                                                          │
   ─────────────── ten lanes opened DURING Wave A, unplanned ─────────────┘
   ‼ lane8  stable stereo ......... LOAD-BEARING: 13% wrong absolute stereo
   ★ lane9  inequivalent donors ... user-requested; REFUTED a "soundness class"
   ★ boron  cage perception ....... RETRACTED a documented permanent ceiling
   ★ atomcount hydrogen ........... 3 inverted premises; 23 molecules pass
   ★ valsearch valence budget ..... 20,000 → 200 candidates, ~100× floor cut
   ★ valorder valence ordering .... found OIN strings that DO NOT RE-PARSE
   ★ encodefail triage ............ 14/48 addressable, 34/48 correct refusal
   ★ genresidue rmsd/other/noconf .. baseline was STALE by one week
   ★ perf   generation ............ a 48–57 s call ran twice per rejection
   ★ encspeed encoder profile ..... AC2BO is 99.8% of a slow encode

Legend
  ◆ the lane's own stated target was measured WRONG
  ✔ delivered      ✖ deliberately delivered nothing
  ✚ deliberate documented deviation from the lane's handoff
  ‼ unplanned and load-bearing      ★ unplanned
  ══ downstream wave (see WAVE-B / WAVE-C)
  Every lever in this wave shipped default-OFF. Nothing here changed user-visible
  output until the Wave D promotion gate.
```

## Initial assumptions and hypothesis

The plan entering Wave A believed the following, all of it inherited from the v0.4.4 capstone
`bucket_report.json` (6,719 molecules, `byte_exact` 81.19% against a 93.51% notation ceiling, with
12.32% of molecules round-tripping to the *same isomer under a different string*):

1. **The 12.32% decomposes cleanly into three owned buckets:** `rdkit_canonical` 500,
   `slot_renumber` 315, `winding_star_drift` 13. Lane 1 owns the first, Lane 2 the second, Lane 3
   the third, and the three are disjoint.
2. **Lane 1 is the big lane and the blocking lane.** `oin/compare.py::canonical_fragment_body`
   already folded body drift at *comparison* time — that is precisely why the key was canonical and
   the string was not — so promoting the same function into the emit path was expected to be
   plumbing. Lane 2 could not start until it existed, because Lane 2's lex-min runs over ligand
   bodies as vertex **colours**.
3. **Lane 3 is small and mechanical:** 13 molecules, cause assumed to be the geometric heading-atom
   fallback tiers, so making those tiers topological should drive the class to 0.
4. **Lane 4 is a generator problem.** Its handoff asserted that multi-axis atropisomers round-tripped
   0/2 in both A/B arms because "the biaryl torsions relax outside the hindered window (20–160°), so
   no axis is detected at all", and prescribed widening the embed pool or guarding the FF relaxation.
   Its Part A required `OIN_EMIT_AXIAL` to be promoted to default-ON.
5. **Lane 7 is pure measurement.** It changes no encoder output, so its acceptance criterion is
   different in kind from every other lane: the suite floor has to hold **exactly**, not merely stay
   above a number, because any movement at all would mean it had changed behaviour it promised not to.
6. **Seven lanes, one plan.** Lanes 1–7 were the whole release.

## What was actually found

### Confirmed

- **`rotate` drift is zero, corpus-wide.** `tools/canonicality_probe.py` holds the molecular graph
  fixed and varies only presentation. Across every arm of every lane's A/B, random *proper rotations*
  of the input produced **0** drifted strings. The encoder was already fully orientation-invariant.
  Every canonicality defect in the release is an **atom-numbering** dependence. This reframed all
  four lanes at once: the enemy is the input file's atom order, not its orientation.
- **`CanonicalRankAtoms(breakTies=True)` is not an invariant.** Over 20 random renumberings of
  `CC(N)=NC` it returned a different ranking **18 times in 20**. What *is* invariant is the graph
  induced by relabelling with `_smilesAtomOutputOrder` — adjacency bytes plus symbol vector, always
  exactly 1 distinct value. Several pre-existing helpers were built on the wrong invariant.
- **A canonical fixed point exists for almost every ligand body.** Of 6,062 distinct capstone corpus
  ligand bodies, **6,056 converge in one reparse pass**; **6 oscillate with period two** (RDKit flips
  `@`/`@@` on adamantane-cage carbons every parse/write cycle), so no fixed point exists for those and
  they bail.
- **Lane 4's descriptor could be made perception-independent**, and the pool did contain the raw
  material selection needed.

### Refuted — three of four lanes' own targets

| lane | its stated target | measured reality |
|---|---|---|
| 1 | `rdkit_canonical` 500 are ligand-**body** drift | **396 of 500 (79.20%) are `slot_or_order`** — body multisets byte-identical, only slots/fragment order differ, i.e. Lane 2's. Only **104 (20.80%)** are reparse-fixable. `rdkit_canonical` is the *fallthrough* branch of `_key_equal_subclass`'s cascade, so any pair mixing two kinds of drift lands there. The headline was an upper bound **off by ~5×**. Over all 828 `key_equal` rows the body population is **104 (12.56%)**, of which 101 go byte-exact from the body fix alone. |
| 1 | promoting the compare-layer reparse closes the gap | The reparse alone: **6 byte-stability fixes, 0 regressions, 0 key-level defects fixed** (an earlier read of the same shape was drift **17 → 17**). It folds two *serializations* of one graph; renumbering hands the serializer a genuinely *different graph*. The key repair came from the **perception** lever, not the reparse. |
| 3 | `winding_star_drift = 13`, cause = geometric heading tiers | The 13 is **stale** — from the v0.4.2 capstone. Re-encoding the same 13 stored pairs on current `main` reclassifies them: **6** still `winding_star_drift`, 3 → `rdkit_canonical`, 3 → `structural`, 1 → `facmer_divergent`. v0.4.3/v0.4.4 had already moved 7 of 13 out. Of the true 6, the geometric tier is load-bearing for **exactly one** (GIPDEQ). |
| 3 | drive the class to 0 | **Refused, on measurement.** All 6 satisfy `generated == mirror(input)` byte for byte, so the question is "is `mirror(M)` the same compound as `M`?" and it differs per molecule. GAMJAG and (latent) MECJOU are achiral → fixable. SIMDIE, TEYXEA, WUHRIB, QIGZAJ are **chiral** — input and generated are genuine **enantiomers**, correctly reported. So `winding_star_drift` goes **6 → 2, not → 0**; **driving it to 0 would fold enantiomers.** 4 of 6 are the *generator* building the wrong enantiomer and the encoder faithfully saying so. |
| 4 | the generator cannot hold two hindered biaryl axes | **The generator holds both.** `tools/injectivity/axial_pool_histogram.py` on `YESKOZ`: rank 0 axis dihedrals **+87.7° / +122.1°**, rank 1 **+128.0° / −14.7°**, rank 2 **+157.2° / −124.5°** — comfortably hindered, spanning several sign combinations. The descriptor could not *perceive* them. Axis selection required both bond ends flagged `GetIsAromatic()`, and the two routes disagree on a metalloporphyrin: encoder `get_tmc_mol` sees **38** aromatic atoms (aromatic pyrrolide core on Zn(II)), generator `build_contract_mol` sees **18** (neutral localized tautomer). "Multi-axis" was a **confound** for "macrocycle whose aromaticity perception is route-dependent". |
| 4 | `\|ax:+-\|` on YESKOZ encodes real chirality | **It does not.** Each configuration's mirror is its exact **complement** (`--`↔`++`, `-+`↔`+-`), so both syn and anti are reflection-invariant and therefore **achiral** about those axes. The old token claimed chirality that does not exist. Worse: **the 37/37 corpus mirror audit was confirming that false positive, not catching it** — a complement satisfies the audit's flip test vacuously. Also explains the 0/2-in-*both*-arms result: a generator limitation separates arms; an ill-posed request does not. |

### Where reality diverged hardest from the plan: ten unplanned lanes

The plan had seven lanes. **Sixteen** lanes were merged into `release/v0.4.5`. By branch count, ten
of the sixteen were not in the plan at all (contemporaneous notes variously say "eight" or "nine"
unplanned — quote the branch list, not a count):

| branch | why it was opened | the finding that made it worth a lane |
|---|---|---|
| `swimlane/v045-lane8` | opened in response to the renumbering finding | **LOAD-BEARING.** 29 of 223 molecules (**13.0%**) emit a *different absolute stereochemistry* under pure atom renumbering; `rotate` gives **0** drift. A soundness bug no prior instrument could see. Became the sixth promoted lever. [Lane 8](LANE-08-stable-stereo-renumbering.md) |
| `swimlane/v045-lane9` | user-requested, after Lane 2 reported "7 molecules with slots on inequivalent donors" | **Refuted the class.** For all 7 the competing strings place donors on vertices related by a **proper rotation** (`identity 4 / in_rotation_group 4 / NOT_A_ROTATION 0`; `\|Δrssd\|` 0.0 … 1.2e-14). Shipped no encoder change, because a fix would have frozen one of two equally faithful labelings. [Lane 9](LANE-09-inequivalent-donors.md) |
| `swimlane/v045-boron` | to confirm a documented "permanent ceiling" | **Retracted it.** 48/48 round-trip behind `OIN_BORON_CAGE`. And found that **14 molecules are scored as passes while describing the wrong graph.** [boron](LANE-boron-cage.md) |
| `swimlane/v045-atomcount` | 74-molecule atom-count class | **Three inverted premises.** 60 of 74 *gain* atoms; the generator is completely innocent (59/59 reproduced from the OIN alone); hydrogen-only, 27/27, `heavy_delta == 0`. [atom count](LANE-atom-count-hydrogen.md) |
| `swimlane/v045-valsearch` | the encoder's 20,000-candidate valence search | The search **saturates at ~100 of 20,000** candidates. `QIDKUL` 148.10 s → 0.92 s with identical `best_BO` and identical OIN. [valence search](LANE-valence-search.md) |
| `swimlane/v045-valorder` | test "is `found_valid = 0` an ordering artefact?" | **Hypothesis refuted**, and a losslessness defect found underneath: the default emits NHC carbenes as aromatic imidazolium and **that string does not re-read** (`KekulizeException`, 3 of 4 fragments parse). [valence order](LANE-valence-order.md) |
| `swimlane/v045-encodefail` | triage 48 `encode_fail` | 14/48 addressable, 34/48 a **correct refusal**. Also: the frozen baseline contains misattributed rows. [encode fail](LANE-encode-fail.md) |
| `swimlane/v045-genresidue` | `rmsd_gate` / `other` / `no_conformers` (81) | **The baseline was stale by a week** — 29 of 33 `rmsd_gate` molecules already passed. [generator residue](LANE-generator-residue.md) |
| `swimlane/v045-perf` | generation cost | A **48–57 s** call ran **twice per rejected conformer**; counter 4 → 2, exact halving, byte-identical. [perf](LANE-perf-generation.md) |
| `swimlane/v045-encspeed` | nobody had ever profiled the *encoder* | **`AC2BO` is 99.8%** of a slow encode, and it is **not** an eta phenomenon — a non-eta 39-atom molecule is slower than an eta 59-atom one. [encoder speed](LANE-encoder-speed.md) |

Lane 5 was **never started in v0.4.5** — see [Wave C](WAVE-C-injectivity-descriptors.md).

## What was done

### Lane 1 — canonical ligand body + canonical perception → [LANE-01-canonical-body.md](LANE-01-canonical-body.md)

**Levers:** `OIN_CANONICAL_BODY`, `OIN_CANONICAL_PERCEPTION` (both default-OFF in the lane, both
promoted in Wave D).

New module `src/oinsmiles/oin/canonical_body.py` — `canonical_body()` (a re-export of
`compare.canonical_fragment_body`, deliberately *not* a copy, so the comparison key, the encoder and
Lane 2's vertex colours cannot drift apart) and `canonical_body_emit(mol, donor_indices)`, hooked in
`src/oinsmiles/utils/xyz2mol.py` at the `lever_enabled("OIN_CANONICAL_BODY")` branch (~line 1746).
Slot identity is carried through the reparse **by atom map number, never re-derived** — the obvious
`GetSubstructMatch` alternative can return a wrong automorphism and put the marker on a CH instead of
the deprotonated X-type carbon (`c{N}` → `[cH]{N}`). Three load-bearing guards (composition, donor
identity, idempotence) plus an aromaticity guard; every failure returns `None` and the caller keeps
the un-reparsed body for the **whole** fragment.

The two perception seams: `lig_checks` candidate ranking now sorts on
`(-N_aromatic, N_pos + N_neg, canonical SMILES)` instead of settling ties on
`ResonanceMolSupplier` order; `AC2BO` is closed **by conjugation** — relabel into canonical order,
perceive there, map the bond-order matrix back. Regression fixture `tests/fixtures/NAXDOI.xyz`:
drift over 6 presentations **3/6 → 3/6 (body lever) → 0/6 (perception lever)**.

A/B, 250 molecules (248 encoded), 3 trials × 3 transforms: levers off 146 byte-stable (58.87%), body
lever 152 (61.29%), both levers 156 (62.90%); key-broken **49 → 49 → 41**.

### Lane 3 — eta winding residual → [LANE-03-winding-residual.md](LANE-03-winding-residual.md)

**Lever:** `OIN_CANONICAL_ETA_WINDING` (default-OFF in the lane, promoted in Wave D). Two fixes:

- **Fix A, unconditional** (`utils/oin_aligner.py` step 4a): both branches of the content-canonical
  tier now call `_topological_heading_atom`, so the geometric heading loop is dead for eta groups
  (retained as a fail-safe). Measured over 1,176 common re-encoded pairs: **2 strings changed, 0
  buckets moved**; GIPDEQ `winding_star_drift` → `byte_exact`.
- **Fix C, behind the lever** (`OINDiscreteAligner._canonical_eta_winding`): colour occupied slots by
  `(chem_id, eta automorphism class)` with **winding excluded** (you cannot decide which slots may be
  relabelled using the characters you are about to reassign), keep only the proper rotations
  preserving every colour, and for each 2-orbit compute the **sense factor ε** with
  `_eta_swap_sense` and apply a measured achirality test. Fail-safe on anything larger.

**ε is not a constant** — measured `+1` for TiCat3/TiCat4, `−1` for SIMDIE and MECJOU — and this is
the lane's near-miss. Its first implementation simply **sorted** the two winding characters, which
made the encoder reflection-invariant for every bridged case and folded a rac pair: **the Y2 axial
bug, reproduced exactly**, passing every guard written against the easy fixture. An independently
built oracle (`tools/eta_core_chirality.py`) caught it by reporting ACHIRAL for the same-sign case and
CHIRAL for the opposite-sign one — the inverse of the assumption. Computed-ε verdicts then agree with
the geometric oracle **6/6**.

Fix C's provable scope is molecules with **≥ 2 eta winding markers**; all 500 such corpus molecules
re-encoded in both arms: `byte_exact` 203 → 204, `winding_star_drift` 5 → 4, **`facmer_divergent`
52 → 52** (no over-folding), 3/500 strings changed, 1 bucket moved. Suite **605 OK / 3 skip /
3 xfail**.

### Lane 4 — axial atropisomerism, Y1 blind spot P2 → [LANE-04-axial-atropisomer-P2.md](LANE-04-axial-atropisomer-P2.md)

**Lever:** `OIN_EMIT_AXIAL`, **kept default-OFF** — a deliberate, documented deviation from the
lane's own handoff. That handoff's Part A required promoting the lever; the confirmed product call
was that **injectivity levers stay opt-in for v0.4.5**. So the lane delivered Part B (the descriptor
fix), recorded the promotion evidence for a later release, and added a guard test asserting the
default is off:
`tests/unit/test_axial_emit.py::TestDefaultOff::test_emit_gate_is_off_unless_the_env_var_is_set`.

`src/oinsmiles/oin/axial.py` now derives axis selection and sign from properties both perception
routes agree on: an axis end is a **trigonal ring atom** (in a ring, 3 heavy neighbours, no H) rather
than `GetIsAromatic()`; reference neighbours are **ring** neighbours; ranks are taken on a
**connectivity skeleton** with bond orders, charges and aromatic flags erased and metal bonds
*dropped* (down-grading them to single bonds closes every chelate ring and buries BINAP's own axis
inside the P–M–P ring); and stereogenicity plus reference are chosen with the **axis bond cut**,
which removes a silent coin toss between two candidates 180° apart.

Recorded promotion evidence, unchanged and still valid for the single-axis cohort: **22/22 (100%)
vs 8/22 (36.4%)** baseline (`axial_cohort_ab.py --limit 12`); fixture pair 2/2 vs 1/2. Recorded
key-folding decision: `_AXIAL_TOKEN_RE` (`src/oinsmiles/oin/compare.py:104`) must **stop** folding
in the same commit that promotes the lever, because folding makes the round trip structurally unable
to verify the one thing the token encodes. The conservative gate now emits nothing for coupled axes,
which is a **real, stated loss**: syn and anti meso-arylporphyrins currently encode identically.
The fix is specified but deliberately not guessed — canonicalize the sign vector over the
automorphism group's orbit, `min(token, complement)` for a coupled pair.

### Lane 7 — research residuals, fixtures and oracles → [LANE-07-research-residuals.md](LANE-07-research-residuals.md)

**Changed no encoder output**, by charter. Its acceptance was therefore that the suite floor hold
*exactly*.

**The two fixtures Lane 5 was blocked on, and which had never been built**, selected by a census of
the 26,230-structure corpus:

| fixture | what | mirror per the oracle | metal descriptor | vdW clashes |
|---|---|---|---|---|
| `tests/fixtures/ZUMNEC.xyz` | tris(catecholato)Mo, 37 atoms | **distinct, 1.34 Å** | `@OH` permutation **11 → 9** | 0 |
| `tests/fixtures/JEGKOW.xyz` | Rh(I) square planar, donors N/P/C(carbonyl)/I, 31 atoms | **not distinct, 0.099 Å** | `@SP` permutation 2 | 0 |

**Why this matters and is not bookkeeping:** `fac-Ir(ppy)₃` is the *easy* case for metal Δ/Λ,
because its three chelates are unsymmetric (C,N) — so a Lane-5 descriptor that secretly encoded
**fac/mer** rather than **helicity** would still appear to work on it. ZUMNEC is the **only**
homoleptic tris-bidentate in the corpus whose two donors per chelate are symmetry-equivalent
(`CanonicalRankAtoms(breakTies=False)` gives them one rank), so it has no fac/mer distinction at all
and metal helicity is its **sole** stereogenic element; it carries no axial or bound-amine axis
either. Its 1.34 Å mirror RMSD is 2.7× the 0.5 Å threshold and 13–27× the achiral controls
(0.05–0.10 Å). And JEGKOW's mirror is correctly *not* a new isomer: a square-planar coordination
plane **is** a mirror plane, so four different donors give **diastereomers, not enantiomers** —
demanding oracle-distinctness there is a category error, and the right operator for `@SP` is a donor
swap.

Also delivered: the **torsion-aware configurational oracle**
(`tools/injectivity/torsion_oracle.py`), 10/10 on cross-validation (EDOQIZ conformational 0.076;
BINAP configurational 3.448 against control 0.032). Building it caught two bugs that made it
confidently wrong — rotatability cut the *atom* rather than the *bond* (a metal is a cut vertex, so
every chelate's metal–donor bond looked rotatable), and automorphisms were enumerated on the
H-explicit graph where methyls starve the budget (EDOQIZ: 4,000 matches → **6** distinct heavy
images, against **864** on the heavy skeleton). The same starvation is a **defect in the pre-existing
`tools/injectivity/oracle.py`**: it inflates rigid RMSD on methyl-rich species (EDOQIZ 2.55 → 0.50 Å
once complete), so the rigid oracle **over-reports chirality at dataset scale** and `uu_hunt`'s
"chiral" gate inherits it. Left unchanged and flagged — it moves no curated fixture verdict.

Twin operators `invert_stereocenter` and `swap_donor`
(`tools/injectivity/twin_operators.py`) plus `report.py --operators`. Donor swap on real geometry:
all four *cis* swaps of PtMeNH₃ClBr and all six of FeH₂(CO)₄ are distinct isomers the key
**separates**; the two *trans* swaps correctly read not-distinct — a trans exchange is a 180°
rotation, i.e. a built-in negative control.

### The unplanned lanes, in one line each

[Lane 8](LANE-08-stable-stereo-renumbering.md) · `src/oinsmiles/oin/stable_stereo.py` plus a 14-line
`xyz2mol.py` hook: do not translate a chiral tag through the fragment rebuild, **re-derive** it from
the parent conformer's coordinates with `AssignAtomChiralTagsFromStructure`, because geometry is not
a function of atom numbering. Measured over 10 of the 29 flip molecules: byte-stable **0/10 → 8/10**,
key-level defects **10 → 1**; permanent guard `tests/unit/test_stable_stereo_mirror.py` (10/10 mirrors
differ — the "stable because constant" veto).
[Lane 9](LANE-09-inequivalent-donors.md) · [boron](LANE-boron-cage.md) ·
[atom count](LANE-atom-count-hydrogen.md) · [valence search](LANE-valence-search.md) ·
[valence order](LANE-valence-order.md) · [encode fail](LANE-encode-fail.md) ·
[generator residue](LANE-generator-residue.md) · [perf](LANE-perf-generation.md) ·
[encoder speed](LANE-encoder-speed.md) — see the table above and each lane's own file.

## Dead ends and refutations

| tried / believed | what killed it |
|---|---|
| "the 500 `rdkit_canonical` molecules are ligand-body drift" | `tools/diagnose_body_drift.py`: **396/500 (79.2%) `slot_or_order`**, 104 reparse-fixable |
| "promote `compare.py`'s reparse into the encoder and the gap closes" | canonicality probe: 6 byte-stability fixes, **0 key defects fixed** |
| build the canonical relabelling on `CanonicalRankAtoms(breakTies=True)` | not invariant — a different ranking in **18 of 20** renumberings of `CC(N)=NC` |
| re-derive the donor's post-reparse position with `GetSubstructMatch` | measured failure mode `c{N}` → `[cH]{N}` on a near-symmetric cyclometalated aryl — a wrong automorphism silently corrupts the coordination sphere |
| naive reparse with no aromaticity guard | Ni-porphyrin: `c1c2nc(…` → `C1=C2C=CC(=N2)…`, broke `tests/unit/test_aromatic_reencode.py` |
| iterate the reparse to a fixed point unconditionally | 6 of 6,062 corpus bodies oscillate with period two (RDKit adamantyl `@`/`@@`) |
| perception retry chosen by "keep whichever scored a higher total bond order" | silently re-imported the order dependence the lever removes; caught only because it broke the `NAXDOI` invariance guard |
| "`winding_star_drift` is 13 and the cause is the geometric heading tiers" | the 13 is a v0.4.2 number; re-encoding gives 6, and the geometric tier is load-bearing for **1** of them |
| drive `winding_star_drift` to 0 | 4 of 6 are genuine **enantiomers**; folding them would destroy stereochemistry. Target corrected to 6 → 2 |
| **sort the two eta winding characters** to canonicalize them | made the encoder reflection-invariant for every bridged case and folded a rac pair — the Y2 bug verbatim. Caught by an independently built oracle, **not** by any guard written against the easy fixture |
| "the generator cannot hold two hindered biaryl axes" | pool histogram: both axes hindered at every rank, several sign combinations present |
| widen the embed pool / guard the FF relaxation / constrain the torsion at embed | nothing to sample harder for; the torsions were never being flattened; and constraining would have *constructed* a twist to satisfy a token that should not be emitted (the project carries three negative results for construction over selection) |
| "`\|ax:+-\|` on YESKOZ records real axial chirality" | its mirror is the exact complement ⇒ reflection-invariant ⇒ **achiral**. The 37/37 mirror audit was **confirming** the false positive |
| trust `tools/injectivity/oracle.py`'s rigid RMSD at dataset scale | automorphism starvation on the H-explicit graph: `ROGYAO_comp_0` called "distinct, ENCODER-BLIND (total)" at 2.586 Å is **achiral** (0.423 Å) |
| demand that JEGKOW's mirror be a distinct isomer | category error — a square-planar coordination plane is a mirror plane; four different donors give diastereomers |

## Where it landed

All Wave A branches merged into `release/v0.4.5` and from there into local `main` via `0d165845`
(tag `v0.4.5`). **Not pushed**, per the standing instruction.

| lane | branch | tip | commits | merge into `release/v0.4.5` |
|---|---|---|---:|---|
| 1 | `swimlane/v045-lane1` | `12569f03` | 5 | `c9ebac35` |
| 3 | `swimlane/v045-lane3` | `1d02ecf6` | 4 | `4c05237d` |
| 4 | `swimlane/v045-lane4` | `3144bac0` | 2 | `43e461e0` |
| 7 | `swimlane/v045-lane7` | `d42194ef` | 8 | `4d92d828` (merged **first** — it changes no encoder output) |
| 8 | `swimlane/v045-lane8` | `b7355cfa` | 4 (incl. `220c191b`, taking Lane 1's perception lever as a shared root) | `6eb82071` |
| 9 | `swimlane/v045-lane9` | `c804044f` | 2 | `bbbfb3f8` |
| perf | `swimlane/v045-perf` | `1e5399c8` | 2 | `9e1fe6aa` |
| encspeed | `swimlane/v045-encspeed` | `3b15e26c` | 4 | `ddbc9fbd` |
| genresidue | `swimlane/v045-genresidue` | `4dcb5284` | 1 | `504f8158` |
| encodefail | `swimlane/v045-encodefail` | `e7d58ca8` | 2 | `2579bfbb` |
| valsearch | `swimlane/v045-valsearch` | `1497cc90` | 6 | `d075f0d6` |
| valorder | `swimlane/v045-valorder` | `7f54454e` | 7 | `e4661843` |
| boron | `swimlane/v045-boron` | `924727c1` | 11 | `f4c3525a` |
| atomcount | `swimlane/v045-atomcount` | `c1ae2759` | — | folded in as the **second parent of the promotion commit** `1450b5ce` |

**Levers created in this wave:** `OIN_CANONICAL_BODY`, `OIN_CANONICAL_PERCEPTION` (Lane 1),
`OIN_CANONICAL_ETA_WINDING` (Lane 3), `OIN_EMIT_AXIAL` semantics fix (Lane 4),
`OIN_STABLE_STEREO` (Lane 8), `OIN_STABLE_METAL_AC` (landed on the Lane 2 branch, see
[Wave B](WAVE-B-canonical-slots.md)), `OIN_BORON_CAGE`, `OIN_H_FAITHFUL`,
`OIN_RESCUE_STUCK_RING`, plus the two `valsearch`/`valorder` levers. **Every one shipped
default-OFF within the wave.** Four were promoted at the Wave D gate.

**Suite,** as each lane reported it on its own branch: Lane 3 **605 OK / 3 skip / 3 xfail**;
`atomcount` **611 OK / 0 failures** on its final code; Lane 6 **620 / 3 / 3** against a 605
baseline. The wave's cumulative suite was verified only at integration —
**837 OK / 3 skipped / 4 expected failures** on `release/v0.4.5`, re-verified on the merged `main`.

**Docs written in the wave:** `docs/CANONICAL_BODY_v0.4.5.md`,
`docs/WINDING_RESIDUAL_v0.4.5.md`, `docs/AXIAL_v0.4.5_LANE4.md`,
`docs/INJECTIVITY_Y3_RESIDUALS.md`, `docs/RENUMBERING_INSTABILITY_v0.4.5.md`,
`docs/WRONG_DONOR_v0.4.5.md`, `docs/BORON_CAGE_v0.4.5.md`, `docs/ATOM_COUNT_v0.4.5.md`,
`docs/VALENCE_SEARCH_v0.4.5.md`, `docs/VALENCE_ORDER_v0.4.5.md`,
`docs/ENCODE_FAIL_v0.4.5.md`, `docs/GEN_RESIDUE_v0.4.5.md`, `docs/PERF_v0.4.5.md`,
`docs/ENCODER_PERF_v0.4.5.md`.

**A hard stop landed mid-wave.** All eight then-running lane agents were terminated by an account
monthly spend limit (`docs/V045_STATUS_2026-07-25.md`). Nothing was lost — every lane was committed
to its own `swimlane/` branch, and the three with uncommitted edits were committed as clearly-marked
`WIP(...) INCOMPLETE` (`3722b18e` Lane 2, `8fdccb55` Lane 8, `1d02ecf6` Lane 3). At that moment
Lanes 5 and 6 had never started, and `swimlane/v045-perf` had zero commits.

## Open questions / for the next agent

1. **`OIN_EMIT_AXIAL`'s promotion evidence is now stale in one specific way, and must be
   re-measured before the lever is flipped.** `_is_atropisomer_candidate` gates its steric wall on
   `not GetIsAromatic()`, and `OIN_CANONICAL_PERCEPTION` is now default-ON — measured on `YESKOZ`,
   hindered axes go **2 → 1** because more of the macrocycle reads aromatic. No emitted string moves
   today (YESKOZ's axes are non-stereogenic so its token is empty either way; BINAP is unchanged at
   1 hindered / 1 emitting), but the Y2 cohort numbers backing the evidence (single-axis 22/22,
   corpus mirror audit 37/37) were taken with perception **OFF**, and `axial.py`'s safety argument
   covers the *generator* reading **fewer** aromatic atoms, not the encoder reading **more**.
   Re-measure both cohorts with perception ON. Reason recorded in
   `src/oinsmiles/oin/levers.py::_HELD_OFF["OIN_EMIT_AXIAL"]`.
2. **The coupled-sign canonicalization Lane 4 specified and did not build.** Distinguish *coupled*
   ambiguity (one automorphism inverts every sign → quotient by global inversion, information
   survives) from *independent* ambiguity (each axis has its own local C2 → nothing survives).
   Per-axis rank comparison cannot express that; it needs the automorphism group's action on the
   sign vector, with a blow-up guard and a conservative fallback. Validate with Lane 7's
   `invert_stereocenter`, since mirroring cannot isolate one axis of a multi-axis molecule. Until
   then the OIN encodes syn and anti meso-arylporphyrins identically.
3. **Lane 3's residual is owned by two other places, not by the encoder.** SIMDIE, TEYXEA, WUHRIB
   need the **generator** to reproduce the coordinated face of a bridged eta ring; QIGZAJ needs the
   **key** to stop folding it — `compare.py::_parse_vertex_colors` colours every slot in a fragment
   with the *whole fragment's* body, so an ansa ligand's Cp and fluorenyl slots get the same colour
   and the tetrahedral C₂ looks colour-preserving. Per-eta-ring colour, as in
   `_eta_automorphism_class`, would separate them.
4. **The rigid-oracle over-report is still live.** `tools/injectivity/oracle.py` enumerates
   automorphisms on the H-explicit graph and starves its budget on methyl-rich species. It changes
   no curated fixture verdict, so `BASELINE.md` §3 does not move — but any *dataset-scale* chirality
   claim taken from it, including `uu_hunt`'s "chiral" gate, is inflated. Port
   `torsion_oracle.py`'s heavy-skeleton enumeration into it, then re-run `uu_hunt`.
5. **`YOYBIY_comp_0`** is the one molecule `OIN_CANONICAL_PERCEPTION` makes worse (its canonical
   valence walk lands on an all-single-bond perception of a bis(pyridylamidine) ligand). It is
   *usable*, so the `get_tmc_mol` retry does not fire. Decide whether the ranking should prefer the
   more-aromatic *assembled* perception; judge it on the amidine/amidinate population, not on this
   one molecule.
6. **The 6 oscillating ligand bodies** will canonicalize if RDKit ever stabilizes the adamantyl
   `@`/`@@` flip. Re-run `canonical_body(canonical_body(cage)) == canonical_body(cage)` after any
   rdkit bump (currently pinned `==2025.9.3`).
7. **Do not re-derive the fixture rationale.** If a future metal-stereo descriptor passes on
   `fac-Ir(ppy)₃` and is not tested on `ZUMNEC`, it has not been tested — that is the whole reason
   Lane 7 built ZUMNEC. Likewise, do not "fix" JEGKOW's non-distinct mirror.
