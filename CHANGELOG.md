# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [0.4.14] - 2026-07-28

> ### The lever worked, the charter didn't — and the charter being wrong is the bigger result.
>
> **`byte_exact` 75.88% → 77.44%, +1.56 points, 78 molecules, 0 moved in a bad direction.**
> Predicted +3.5–4.1. **The miss is the finding:** 114 of the 204 molecules this release was
> chartered to close were never reachable by the kind of change it proposed.
>
> `OIN_RESONANCE_DONOR_FOLD` ships **default-ON**. It ranks a fragment's *constitutional
> skeleton* — bond orders, aromatic flags, charges and hydrogens erased, connectivity, element and
> chiral tag kept — so acac written ketone/enol, carboxylates and sulfonates stop being read as two
> inequivalent donors. Ester `-O-`/`=O`, ether/ketone and amide N/O still do not merge.
>
> **🔴 Two blocks were re-filed off the encoder ladder. Neither number moved; who owns them did.**
>
> | | n | pts | was | is |
> |---|---:|---:|---|---|
> | `rdkit_canonical` | 114 | 2.28 | "RDKit canonical write order" | **80.7% η-set denticity drift** — the generated ring slips, fewer carbons fall inside the bonding cutoff. No string change can reach it |
> | veto-reverted `slot_renumber` | 183 | 3.66 | "benign canonicalization" | **the generator built the ENANTIOMER.** `byte_exact` failing is CORRECT |
>
> The second one is the serious one. `key_equal` is documented as *"benign canonicalization — the
> win reclaimed"*; for **183 of its 361 remaining members (50.7%)** that is false. They are
> enantiomer pairs, invisible because `compare._parse_vertex_colors` folds reflection deliberately
> — **and `accept_fn` decides by that same key**, so the generator accepted a mirror-image
> structure and the harness recorded it as a same-isomer string difference. Measured, not inferred:
> `tools/veto_residue_chirality.py` re-encodes the input, the stored round trip, and the *mirrored*
> input, and reads **222/222 classified, 0 excluded, 183 `MIRROR_MATCH`**.
>
> That is the **third** time this project's headline has rested on a metric that folds the axis
> under test, after v0.4.8 (scored vs honest) and v0.4.11 (the fold itself). The transferable form:
> **a bucket name that asserts a cause is a hypothesis, not a measurement.**
>
> **🔴 The gate that passed was blind, and saying so is the point.** A uniform 250-molecule mirror
> audit contains **1 of the 179 moved molecules — 0.4% coverage** — so its identical before/after
> tally proves only that the lever does not damage molecules it never touches. That is v0.4.13's
> ARM 1 failure reproduced exactly, and it would have shipped as clean evidence if the coverage
> had not been computed. The gate that *can* see the change runs on a mover-enriched cohort at
> **179/179 = 100%**: 0 regressions, **0 per-molecule verdict changes**, `achiral_or_preexisting_fold`
> unmoved at 71, and — the number that separates "safe" from "never fired on anything chiral" —
> **33 of the 78 gains are on molecules the encoder resolves as chiral and still resolves after.**

### Added
- **`OIN_RESONANCE_DONOR_FOLD`, promoted default-ON.** Widens `OIN_CANONICAL_DONOR_FOLD`'s
  donor-equivalence test to the fragment's constitutional skeleton
  (`canonical_slots._skeleton_ranks`). Grouping is union-find (`_merge_classes`), never a composite
  key, so the partition can only get **coarser** — a composite key could move a slot into a
  different bucket and lose a labeling the shipped encoder already reaches. Subject to the same
  coupling invariant as the fold it widens: only safe with `OIN_FOLD_PARITY_VETO` on, pinned by
  `test_resonance_donor_fold::TestResonanceFoldInheritsTheVetoCoupling`.
- **`tools/veto_outcome_audit.py`** — separates the parity veto's five outcomes, which no bucket
  report can, since every reverted molecule lands in the same bucket regardless of why. Over
  **393/393 movers, 0 excluded: 222/222 `vetoed_collapse`, 0 `no_evidence`.** `levers.py` asserted
  the veto was alive at corpus scale; this measures it.
- **`tools/veto_residue_chirality.py`** — settles whether the veto-reverted residue is a generator
  or an encoder problem, rather than inferring it from `vetoed_collapse` (which is a statement
  about one structure and its mirror, not about input vs round trip).
- **`tools/resonance_transition_sim.py`** — both arms re-encode from coordinates, because the state
  being compared against runs the parity veto and therefore needs a conformer.
- `tools/fold_key_invariance.py` grew `--lever` / `--holding`. `--holding` is not a convenience:
  measuring the widening against a fold-OFF baseline would report the *fold's* movement as the
  widening's.

### Changed
- `docs/agentic-notes/ROADMAP_100_100.md`: gap re-derived to **`100 − 77.44 = 22.56`**, and the
  decomposition re-filed by owner — **5.94 points move off the encoder ladder**. The encoder ladder
  has **1.28 points** of reachable work left (39 `NOT_A_MIRROR` + 25 resonance residue); the rest of
  the distance to 100% is generator work.
- arm2 goldens re-frozen for the promotion — **7 of 325 rows** in `gate_v049_arm2_golden.tsv` and
  **1 of 100** in `gate_v047_arm2_golden.tsv`, `MANIFEST_SHA256` recomputed (arm2 does not verify
  it, so a stale one goes unseen). ARM 1 is byte-identical, `#DONE 62`. *Coverage of the moved
  population is a different count and is 12/325 and 2/100 — the v0.4.9 golden predates the donor
  fold and some of its labelings coincide with today's.*
  Two traps caught in the re-freeze: only fields 1–6 are taken from the fresh run, because a v0.4.9
  golden's field 7 is the **band** `--band` filters on while a fresh row carries `xyz_sha` there;
  and **field 3 is a fresh stochastic generation**, so the two rows whose round-trip status changed
  (`HEKFEL`, `FOJJUM`) are **not** lever-caused — the deterministic control on the frozen sweep
  structures reads `True→True` and `False→False` respectively. Both facts are recorded as comments
  inside the goldens.

### Fixed
- Nothing. This release changes what the encoder emits for 78 molecules and what four numbers on
  the roadmap mean; it fixes no reported defect.

## [0.4.13] - 2026-07-28

> ### The donor fold ships — and the release that was supposed to cost 55 CPU-hours cost none.
>
> **`byte_exact` 72.46% → 75.88%, +3.42 points, 171 molecules, 0 moved in a bad direction.**
> The first time this project's headline has moved **up**; the only other time it moved at all was
> v0.4.8, which took it down 10.34 points on purpose to make it honest.
>
> `OIN_CANONICAL_DONOR_FOLD` + `OIN_FOLD_PARITY_VETO` go **default-ON together**. The gate ran on
> **two** populations, not one, and both zero with exact accounting:
>
> | | cat/ | cat+photo |
> |---|---:|---:|
> | collapses, veto OFF → promoted | **19 → 0** | **33 → 0** |
> | `achiral_or_preexisting_fold`, unmoved | 157 | 134 |
> | accounting | 73 + 19 = 92 ✓ | 83 + 33 = 116 ✓ |
>
> Every collapsing molecule moved into `distinct_both_arms` and **nothing else moved at all** — so
> the zero is bought by *separating*, not by *abstaining*, which is how v0.4.12's own first veto
> implementation failed (it declined on 18 of 18 while all three fixture tests passed).
>
> **🔴 The 55 CPU-h re-sweep v0.4.12 made a precondition was the WRONG INSTRUMENT, and was not
> run.** `results-v0.4.8-honest` — the 72.46% baseline — is itself an offline re-score, not a
> sweep; no generator sweep has run in this project since v0.4.6. A fresh run would have
> contaminated a +3.42-point signal with the stochastic generator's run-to-run variation and would
> not even have been like-for-like with the table it was compared against. What made the offline
> route *exact* was measured, not assumed: over all **9669** corpus strings the fold moves **1019**
> and changes **0** comparison keys, so `accept_fn` — which decides by key — returns bit-identical
> conformers either way.

### Added
- **`OIN_CANONICAL_DONOR_FOLD` and `OIN_FOLD_PARITY_VETO` promoted to default-ON**, together and
  never separately. `test_levers::TestDonorFoldAndParityVetoAreCoupled` pins the pair: the fold ON
  with the veto OFF is the configuration v0.4.11 refuted, and **neither `byte_exact` nor the
  round-trip key can see the damage** — only `tools/mirror_audit_donor_fold.py` can.
- `OIN_PREFILTER_ADVISORY` (**default OFF**) — makes the cheap acceptance prefilter advisory
  instead of dispositive, with `overridden` / `confirmed` / `cheap_pass` telemetry so the decision
  is observable. On `AROHIA_comp_0` the cheap test rejects **2** conformers the strict test accepts.
- `tools/attach_class_audit.py` — sizes the MEDZUR and GAVSED classes, both handed forward on
  **n = 1**. Over 767 genuine failures: **GAVSED 280, MEDZUR 99**, with a `byte_exact` control at
  1.32% against 24.11% (**18.2× enrichment**) and `UNKNOWN` = 0.
- `tools/fold_key_invariance.py` — the generator-neutrality proof above.
- `tools/prefilter_prevalence.py` — Lane 1's two-arm harness; distinguishes `INSTRUMENT_DEAD`,
  `NO_POPULATION`, `PREFILTER_VINDICATED` and `DEFECT_CONFIRMED`, which would otherwise all read 0.
- `tools/freeze_v0413_table.py` — refuses to freeze a table from an arm that excluded anything.

### Changed
- **`structural` is re-mechanised, not re-labelled.** 417 molecules / 8.34 pts, scheduled v0.4.17
  as *"bounded by what the generator can assemble"* — **266 of them (63.8%) are `DETACHED`**: the
  generator assembled a structure and *returned* it with ligands off the metal, because
  `_select_by_geometry`'s fallback ranking is not attachment-aware. A one-site return-path guard,
  not a capability floor.
- ARM 2's v0.4.9 golden: **11 of 325 rows re-frozen** (the fold-movers), re-run individually rather
  than bulk-accepting a regenerated manifest. ARM 1 and the v0.4.7 ARM 2 golden are unchanged —
  **provably**, not as an economy: 0 of their molecules are fold-movers.

### Fixed
- 🔴 **`tools/fold_transition_sim.py` printed the REFUTED number after measuring nothing.** Run
  from a worktree its `--dataset` default is *relative*, so all 393 movers landed in `unavailable`
  and it reported the bare fold's **+7.86** — the figure v0.4.11 refuted for collapsing
  enantiomers — under a heading saying "veto", and exited 0. It now refuses when it measures 0 of
  N movers.
- 🔴 **`tools/run_sweep.sh` resolved its interpreter by globbing sibling checkouts** — the trap
  `gate_v047.sh` was hardened against in v0.4.9, left live in the one place it costs 55 CPU-h.
  From the v0.4.13 worktree it selected `EtaCatalysis/.venv` (no rdkit); the next candidate carries
  rdkit 2025.09.2 against the pinned 2025.9.3. Now resolves via `--git-common-dir` and refuses on
  version drift.
- A lever-interaction bug introduced by `OIN_PREFILTER_ADVISORY` and caught before it shipped:
  without `and not cheap_vetoes`, that lever plus `OIN_ACCEPT_SCORED` would together accept a
  conformer the score itself calls a failure — a combination neither lever's own A/B exercises.

### Known limitations
- **Lane 1's corpus prevalence is NOT measured.** The lever, telemetry and harness are built and
  the defect is confirmed on `AROHIA_comp_0`, but n = 1 is exactly what this project forbids
  quoting. Handed to v0.4.14 with its instrument gated.
- **ARM 1 is structurally blind to this release**: 0 of its 62 fixtures are fold-movers, so its
  PASS means "no regression", not "the promotion works". Any future canonicality promotion must
  state its gate's *coverage of the moved population*, not just its verdict.
- The **48 `byte_exact` molecules that read `DETACHED`** are unexplained; two readings survive and
  neither was tested.

## [0.4.12] - 2026-07-28

> ### Reflection parity — the filter v0.4.11's refutation demanded, and it works.
>
> v0.4.11 built the donor fold, measured **+7.86 `byte_exact` points across 393 molecules**, and
> refuted it: the fold **collapses enantiomers in 221 of those same gains**. This release builds
> the reflection-parity filter its close-out specified.
>
> **The gate reads 0.** A uniform 250-molecule mirror audit goes **19 → 0**
> `REGRESSION_raw_collapsed`, with those 19 moving into `distinct_both_arms` (73 → 92) and
> `achiral_or_preexisting_fold` **unmoved at 157** — the direct evidence that the veto did not
> simply refuse to fold everything.
>
> **The surviving gain is +3.42 points** (171 of 393). v0.4.11 bounded the safe set at *"at most
> ~172 (~3.44 pts)"* by counting collapses; this filter counts survivors through the shipped
> predicate and lands on 171. Two independent routes agreeing to one molecule.
>
> **The headline is FLAT at 72.46% and that is deliberate** — both levers ship default-OFF, so
> the default path is byte-identical (ARM 1 PASS, 62/62; suite 988 OK).

### Added
- `OIN_FOLD_PARITY_VETO` (**default OFF**) — declines the donor fold on any molecule where
  folding would make the structure's mirror encode identically. Lives in `get_oin_string` on the
  **pristine** conformer, because **reflection parity is not a property of the emitted string**:
  a donor swap is a transposition fixing every other vertex, so the obvious `det > 0` test on
  the polyhedron rejects *every* swap and degenerates the fold to the identity.
- `OIN_ETA_ACCEPT_EXIT` (**default OFF**) — the eta winding criterion relocated from
  `_select_by_geometry_impl` into `accept_fn`, the only site consulted *during* pool filling.
  A conjunction with geometry classification and ligand attachment, because winding alone would
  bypass clash-first ranking — the defect `OIN_ACCEPT_SCORED` has.
- `tools/fold_transition_sim.py` — commits the `+393` measurement v0.4.11 made with an
  uncommitted script. Its veto arm re-encodes from coordinates and carries a drift control.
- `tools/ab_accept_scored.py` gains `--lever` / `--extra-env` and a **G5 metal-configuration
  arm**. G1–G4 are structurally blind to Δ/Λ because `compare.py` folds `|mc:|`, so an arm can
  return the opposite enantiomer and every gate reports "identical".

### Changed
- `OIN_ETA_EARLY_EXIT`'s `_HELD_OFF` entry now records that it is **runtime-inert as sited** and
  that its promotion gate is **void, not unrun** — it runs downstream of a fully-filled pool.
- The carry-forward licence is now **measured**: re-encoding the frozen corpus's inputs *and*
  stored generated structures reproduces the v0.4.8 strings on all 393 movers (0 drift), so
  v0.4.9/v0.4.10/v0.4.11's byte-identity claims are confirmed rather than inherited.

### Fixed
- Two silent defects in the veto's own construction, both caught only by disbelieving a clean
  result. `tmc_mol`'s atom order is **not** the coordinate order (`__origIdx`) — zipping
  positionally encoded `BIWDIV` as `[Co_TBP]` with invented bonds. And the mirror was encoded
  with the fold *inherited*, which disarmed the achiral guard and made the self-check decline on
  **18 of 18** movers **while all three fixture tests passed** — because declining to fold also
  separates a mirror pair. `resolve()` now records *why* it decided, and tests assert the
  outcome rather than the string.

### Known limitations
- `OIN_ETA_ACCEPT_EXIT` ships with its correctness arms and **no runtime claim**. Its timing A/B
  was stopped rather than banked at load 26 (*never interleave timing runs with gate runs*), so
  its predicted tail reduction is **unverified**.
- The chartered v0.4.6 accept-gap cohort is **stale on 8 of 8** molecules; any cohort frozen
  before v0.4.8 must be re-derived. The real eta target population is **405 molecules whose key
  never matches, 378 of them > 30 s**.

## [0.4.11] - 2026-07-27

> ### `slot_renumber` — the fold works, buys 7.86 points, and must not ship.
>
> This release built the fix v0.4.5 Lane 2 specified in writing for the **largest single block in
> the gap** (496 molecules / 9.92 points), measured it, and **refuted it**. The within-fragment
> donor fold does exactly what it was designed to do: **393 molecules move
> `key_equal/slot_renumber → byte_exact`, none in any other direction**, `facmer_divergent` holds
> at 16, the comparison key moves on **0 of 992** strings, and both gate arms stay byte-identical
> with the lever off.
>
> It also **collapses enantiomers.** A corpus mirror audit on a **uniform** draw found 19 of 250
> structures (7.6%) whose mirror image encodes identically once the fold is on — and run directly
> on the 393 molecules the fold *claims as wins*, **221 of 393 (56.2%)** collapse. An independent
> geometric oracle (`tools/injectivity/oracle.py`, which shares no machinery with the encoder)
> confirms **18 of 19** uniform-draw collapses and **26 of a 30-sample** of the 221 are genuinely
> chiral; three are cap-free and unambiguous. **More than half the gain is the damage**, so at most
> ~172 of 393 (~3.44 of the 7.86 points) are safe.
>
> **Why the safety argument failed.** It assumed two donors in one `breakTies=False` symmetry class
> are interchangeable. `CanonicalRankAtoms` computes the symmetry of the **isolated ligand graph**,
> but those donors sit at distinct vertices whose relation to the other ligands is chirality-bearing
> — so the vertex permutation their exchange induces can be **improper**. *A fragment's automorphism
> says nothing about the parity of the vertex permutation it induces.* v0.4.5's restriction to proper
> rotations was not conservatism; it was the load-bearing correctness condition.
>
> **The finding that outranks the points: `byte_exact` can be raised by deleting information, and
> the comparison key will agree.** Both are blind to reflection, because `_parse_vertex_colors`
> folds that axis deliberately. A one-directional transition matrix is *not* evidence of safety for
> anything touching canonicalization. Mirror-audit every future canonicality lever, on a uniform
> draw, before quoting its points.

**Accuracy: FLAT — 72.46%.** `OIN_CANONICAL_DONOR_FOLD` ships **default OFF** and no default answer
moves. ARM 1 **62/62 byte-identical**, ARM 2 **90/90 gated**. The `BASELINE.md` §1 carry-forward
licence therefore **survives**: v0.4.12 does not owe a re-sweep.

> ⚠ **Two independent cohorts find it** — uniform 19/250 (7.6%) and the runtime-stratified cohort
> 31/300 (10.3%) — so the damage is not an artifact of one sample. The named Δ/Λ fixtures
> (`ZUMNEC`, `fac-Ir(ppy)₃`) **both pass with the lever on**: neither carries the vulnerable motif,
> and fixtures alone could never have caught this.
>
> ⚠ **A methodology error made during this release, recorded rather than buried.** The stratified
> audit was read mid-run as "0 regressions in the first 200" and briefly written up as a *"wrong
> stratum"* finding. The tool prints a verdict only every 50th molecule, so four clean progress
> lines had been mistaken for 200 clean molecules; the completed run reports 31 collapses. **A
> partial run is not a result, and a progress line is not a tally.**

### Added
- `tools/mirror_audit_donor_fold.py` — the instrument that caught the collapse. Reports a four-way
  table and defines a regression as an implication (`OFF_distinct and not ON_distinct`), so a
  pre-existing fold is never miscounted against the lever under test. **Should gate every future
  canonicality lever.**
- `tools/slot_drift_mechanism.py`: `--roundtrip` (classify round-trip pairs, not just
  re-presentation pairs), `--expect N` (abort unless exactly N pairs are selected), and
  `--explain-distinct` (re-test `distinct_donors_LOCAL` resonance-insensitively).
- `OIN_CANONICAL_DONOR_FOLD`, registered in `_HELD_OFF` with the full refutation and the condition
  for promotion (a reflection-parity filter).

### Measured
- **All 496 `slot_renumber` molecules classified for the first time** — the taxonomy had only ever
  run on re-presentation pairs. `same_vcolor_identical` **496/496**; `diff_occupancy`,
  `diff_geometry`, `diff_colors` and `postpass_BUG_diverges` all **0**. The charter's headline risk
  is refuted *in the release's favour*: the block is entirely reachable in principle.
- **90 of the 118 `distinct_donors_LOCAL` are frozen-resonance-form artifacts** — acac binds through
  two equivalent oxygens but is written ketone/enol, so `CanonicalRankAtoms` separates them. That is
  a ligand-**body** gap, so **v0.4.14 is re-sized from 2.28 to ~4.08 points**.
- Two stale v0.4.5 priors corrected: the canonicality probe reads **35 / 28 / 7 / 0** today, not
  32 / 23 / 7 / 2 — the 2 unparsable are gone because `OIN_BORON_CAGE` was promoted in v0.4.6, after
  the prior was recorded.

### Changed
- `TestResidualClassIsOutOfReachByDesign` **inverted, not deleted**, as v0.4.5 asked: the residual
  class *is* reachable. What is not true is that reaching it this way is safe.
- The fixture-based mirror test is renamed to say what it actually proves — nothing — because
  passing fixtures are how this nearly shipped.

### Added (tests)
- `TestDonorFoldCollapsesEnantiomers` pins the defect at string level (`BIWDIV_comp_0` and its
  mirror), to be **inverted, not deleted**, when a parity filter lands.
- `TestDonorFoldScope` tests the three scope conditions by removing one at a time — asserting on
  `_donor_swap_permutations` **directly**, because the rotation group already converges most small
  cases and an equal-output assertion passes whether or not the fold over-reached.

## [0.4.10] - 2026-07-27

> ### Cost per attempt — byte-identical by construction, and the arbiter was broken before we started.
>
> **The first thing this release measured was its own gate, and the gate was dead.**
> `tools/gate_v047.sh arm1` had been exiting 1 *before comparing anything* since `dd51a515`:
> `ULODUU_comp_0.xyz` was added as a fixture and the golden was never extended, so the
> `EXPECTED_FIXTURE_COUNT = 61` guard hard-refused every run. v0.4.9 froze a 328-molecule benchmark
> and measured a 0.28% noise floor while the encoder arm of its own gate was non-runnable.
>
> Behind that refusal sat a second, real drift: `ASISAX_comp_0`'s frozen row carries
> `ERROR:ValueError:xyz2mol failed:` and the tree emits `perception_tmc failed:` — the v0.4.7 rename,
> recorded as *behaviour-neutral*. It was, of behaviour. It was not neutral for the error **string**,
> and ARM 1 hashes error strings deliberately. **The other 60 rows are byte-identical**, so the
> encoder itself has not moved; both rows are re-frozen at 62 with the diff recorded.
>
> **The transferable finding: a gate that fails before it compares is indistinguishable from a gate
> that is merely inconvenient to run, and it silently stops covering everything else it was watching.**

**Accuracy: FLAT, by construction.** No answer changes. Byte-identity is the primary acceptance gate,
not a secondary check, and both arms pass for both lanes: **ARM 1 62/62 byte-identical, ARM 2 90/90
gated**, plus identical generated-structure fingerprints across every A/B run.

> ⚠ **Every speed number here is bimodal by molecule and none of them is a corpus number.** The same
> two changes measure −32.9% and −86.7% on the molecules they target, and −2.4% and *nothing* one
> molecule over. v0.4.9 shipped a lever aimed at whichever function profiled expensive last and
> measured nothing; this release reports both halves.
>
> ⚠ **Lane A's figures were corrected before release.** A first A/B taken while four gate processes
> competed for the box (load average **35**) reported **−50.2%** and **+9.6%**; quiet-box
> re-measurement gives **−32.9%** and **+0.3%**. The gain was over-stated by 17 points and the null
> was buried in 30% within-arm spread. **Byte-identity gates are load-immune and can be parallelised
> freely; wall-clock is neither** — that is the method note this release paid for.

### Changed
- **The discarded `.index()` scan is gone** (`generator3d/embed.py`). `get_embedding`'s outer loop
  called `alternative_ace_mol_list.index(alternative_ace_mol)` and **threw the result away**. Not a
  no-op: `list.index` compares with `Molecule.__eq__` → `is_same_molecule` → `get_c_eig_list` →
  `numpy.linalg.eig`, so it eigendecomposed a Coulomb matrix per candidate per outer iteration,
  O(n²) in the candidate list. **`CAHQEJ_comp_0` 86.35 s → 57.96 s (−32.9%)**;
  **`FOSNEI_comp_0` +0.3% — no effect**, and attributed rather than excused: that molecule makes
  **3** such comparisons costing 0.03 s, against `CAHQEJ`'s **99** costing 38.52 s. Ships **on by
  default, with no lever** — gating provably dead code would permanently ship
  `if not lever: <discarded computation>`.
- **ARM 1 golden re-frozen at 62 rows**, `EXPECTED_FIXTURE_COUNT` and the `#DONE` sentinel 61 → 62.
  The two counts stay hardcoded and independently asserted: deriving one from the other would make
  this class of drift self-heal, which deletes the guard rather than fixing it.

### Added
- **`OIN_MEMO_CIP_REPARSE` (default OFF)** — memoises `chirality._reparse_cip_label_once`, a pure
  function of three immutable scalars that costs **2.43 s a call**. **`VAFMIA_comp_0` 81.89 s →
  10.87 s (−86.7%)**, ~310× the noise floor. That molecule set v0.4.9's budget ε: the bound can
  decline to *start* the next `accept_fn` call but cannot interrupt the running one, which is why it
  held to 2.09× rather than 1.0×. At 10.87 s VAFMIA now clears 30 s — **the worst of the eleven
  molecules that exceeded budget in both of v0.4.9's arms, and only that one.**
  ⚠ The helper is also on the **encoder** path (`chirality.recover()` →
  `_reparse_aromatic_cip_label`), which the charter's generator-side framing hides, so this lane is
  gated on ARM 1 as well as ARM 2. Both pass.
- `tests/unit/test_embed_dead_scan.py` (5) — pins the *purity* properties the deletion rests on, so
  adding the "obvious" cache to `get_c_eig_list` fails there with an explanation instead of
  surfacing as an unattributable string diff in a corpus run. Plus a lint against reintroduction.
- `tests/unit/test_cip_reparse_memo.py` (11) — warm equals cold on every branch; the memo key is
  **complete** (`O[C@:99]([O])(C)CC` returns `None` with `fill_deficit` and `"R"` without, checked in
  both orders on a cleared cache); default OFF produces zero cache traffic; `=0` disables.

### Not done, deliberately
- **`OIN_MEMO_CIP_REPARSE` is not promoted.** The charter permits same-release promotion *"only if
  byte-identity holds on the whole benchmark"*. This release ran the **fast band — 90 of 328**. The
  promotion gate is the full cohort with the lever on, ~10 CPU-h sharded 6-way, against the existing
  frozen golden. Fast-band-only evidence is not that.
- **`max(elapsed_s) < 30 s` is still not delivered**, and this release does not claim to approach it.
  v0.4.9 measured the same **11** molecules over budget in both arms; one of them is now fast.

## [0.4.9] - 2026-07-27

> ### Speed becomes measurable — and this release refutes its own premise.
>
> v0.4.9 was chartered on one number: *"759.9 s against a 300 s budget. That single number is this
> release's justification: the budget is not a budget."* **That number is arithmetic on a sum.**
> `metrics.elapsed_s` accumulates up to three separately SIGKILLed harness attempts. Split by
> `tier_passed`, all **4658** single-attempt rows in the 5000-molecule sweep finish within
> **0.2 s** of their 300 s cap — the harness enforces to ε ≈ 0.2 s.
>
> The advisory-timeout defect is real; it is in the code and two direct-call probes measure it
> (60 s asked, 137.9 s and 172.8 s spent). The corpus number was simply never evidence for it.
>
> **What replaces it is a stronger argument.** The `≥ 300 s` band is 291 molecules burning
> **27.3 CPU-h — 49.8% of the entire 54.8 CPU-h sweep — for three honest passes.** And 93.1% of
> all honest passes already finish under 30 s, so a per-molecule cap is far cheaper than assumed:
> 30 s recovers **37.8 CPU-h** and costs **251 passes (5.02 points of `byte_exact`)**; 300 s
> recovers 3.1 CPU-h for 3 passes.

**Accuracy: unchanged.** `byte_exact` is **72.46%**, the honest v0.4.8 baseline, confirmed here by
a live 300-molecule run: **296/300 per-molecule agreement** with the offline re-score and
**0/300 encoder drift**. No encoder, generator, or notation behaviour changed by default.

### Added
- **`OIN_ENFORCE_BUDGET` (default OFF)** — makes `OIN3DGenerator(timeout=)` a bound instead of a
  hint. Checked in three places, because there is no single cost sink: inside `get_embedding`'s
  nested loops, before each `accept_fn` re-encode, and before `_select_by_geometry`.
- **`BudgetExhaustedError`** — budget exhaustion is now typed and distinguishable from an assembly
  failure, so the next release can tell its own regressions from this release's intent.
- **A frozen, stratified runtime benchmark** — `tools/gate_v049_arm2_golden.tsv`, 328 molecules
  across four runtime bands × eta/non-eta plus a fast control, drawn from one source by a
  deterministic rule, every row labelled with its honest round-trip class. **Reproduces to 0.28%**
  (277.01 s vs 277.79 s, byte-identical rows) — the noise floor every later runtime claim must clear.
- `tools/select_runtime_strata.py`, `tools/budget_bound_ab.py`, `tools/budget_bound_report.py`.
- `--gen-timeout` on the sweep harness, and `--hard-timeout` / `--shard` / `--band` on the gate.

### Fixed
- **`--mol-timeout` never reached the generator.** The harness hardcoded the generator's budget, so
  the two were fully decoupled and a bound could not be A/B-ed through the harness at all.
  Behaviour-identical for every invocation this project has run.
- **The v0.4.7 gate silently resolved to an arbitrary interpreter** — from a worktree it selected an
  unrelated project's venv, once with rdkit 2025.09.2 against the pinned 2025.9.3. A byte-identity
  gate on a different rdkit reports MISMATCHes that read as code regressions.
- **The gate itself was unbounded**, the same defect it was measuring.

### Measured, not delivered
- **`max(elapsed_s) < 30 s` is NOT delivered.** With the bound on, ε = **+32.8 s** on a 30 s
  budget and the same **11** molecules exceed it in every arm. The bound compresses the tail
  (2.63× → 2.09×); it does not remove it. ε is one *in-flight* `accept_fn` re-encode, ~24 s of
  `chirality._reparse_cip_label_once`.
- **The honest gap reorders the roadmap.** `structural` was 9 molecules / 0.18 points when scored
  dishonestly; honestly it is **417 / 8.34 points — the second-largest block in the gap**, ahead of
  `hard_fail`, and the ladder parks it at v0.4.16 as "knowledge, not points".
- **A 22% win, found and deliberately not taken** — `get_embedding`'s outer loop calls
  `alternative_ace_mol_list.index(...)` and discards the result: 3711 calls, 198
  eigendecompositions, 22% of an eta generation. It belongs to v0.4.10.

## [0.4.8] - 2026-07-27

> ### ⚠ The headline accuracy figure goes DOWN, on purpose. `byte_exact` is **72.46%**, not 82.80%.
>
> **This is a re-baseline, not a regression.** No encoder, generator or notation behaviour changed
> in this release — the *measurement* changed. Every accuracy number this project has published
> before now is over-stated, and the correct comparison for anything after this point is 72.46%.
>
> The harness scored a round trip with `get_oin_string(gen_result.mol, coords)` — **the generator's
> own bond graph**. That is not merely inaccurate, it is circular: `gen_result.mol` is exactly the
> artifact that would have to be wrong for the test to fail. A ligand could drift a full ångström
> off the metal and the graph would still call it bonded. `FIYHUT_comp_0` ships both Cp rings 0.85 Å
> off the iron — 10 bonded carbons to 0 — and scored a byte-exact pass.
>
> The same shortcut erred in the other direction too, discarding stereo the coordinates *do*
> support: `YOSXIP_comp_0`'s `[S@]{5}` sulfoxide flattened to `S{5}` and was scored a mismatch.
> Both directions close with one call — a full `XYZToSMILES().convert()` of the generated XYZ.

### Changed
- **`OIN_INDEP_SCORE` promoted to default-ON.** The round-trip verdict is now re-derived from the
  generated coordinates alone. Measured over the same 5000-molecule corpus, same conformers, same
  key and same `status` gate:

  | bucket | scored | honest | delta |
  |---|---:|---:|---:|
  | `byte_exact` | 4140 (82.80%) | **3623 (72.46%)** | **−517** |
  | `key_equal` | 520 (10.40%) | 610 (12.20%) | +90 |
  | `facmer_divergent` | 1 | 16 | +15 |
  | `structural` | 9 (0.18%) | 417 (8.34%) | +408 |
  | `hard_fail` | 315 | 319 | +4 |
  | `encode_fail` | 15 | 15 | ±0 |

  613 molecules degraded, **30 improved** — the correction runs in both directions. At corpus
  scale **36.7% of haptic `byte_exact` passes were false**, against 6.7% non-haptic.

- **Scope is deliberately narrow: this changes what is *reported*, not what is *accepted*.**
  `accept_fn` is untouched and the harness's tier ladder still escalates on the old predicate.
  Scoring the ladder honestly would move runtime *and* the failure mix in the same release that
  re-baselines the number, making both unmeasurable. Runtime did not move: nothing was
  re-generated, so `metrics.elapsed_s` is unchanged (`> 30 s` stays 994 / 19.88%).

- **`levers.py`'s rationale for holding the lever off was rewritten.** It priced the second encode
  at "0.4–1.5 s/molecule" and cited the cost as the reason to wait. That figure had never been
  measured. It is **0.33 s/molecule** — a low single-digit percentage of the sweep it corrects.

### Added
- **`tools/honest_rescore.py`** — re-scores a completed sweep offline, with no generator run.
  `save_artifacts` stores the same `gen_result.xyz` string the lever converts, so re-encoding it is
  bit-identical to the lever rather than an approximation. The full corpus re-scores in **334 s**
  against the ~55 CPU-hours a live re-sweep costs, and because conformers are held fixed by
  construction it also removes the A/B confound: `smiles_1` cannot move and no second encode enters
  the `--mol-timeout` budget. Resumable, `#DONE` sentinel, `--fill-coordination` backfill.
- **`tools/encoder_identity_corpus.py`** — corpus-scale encoder byte-identity gate. A bucket report
  re-run over stored JSON proves only that the *classifier* did not move; it never encodes anything.
  Result for this release: **4985/4985 byte-identical, zero drift.**
- **`tools/roundtrip_bucket_report.py --score {scored,honest,both}`** — `both` emits the two tables
  and the per-molecule transition matrix. A single number replacing another with no reconciliation
  is precisely what makes a re-baseline read as a regression.
- **`tools/atom_count_provenance.py`** and the Lane 2 verdict: the `Atom count mismatch` gate is
  **load-bearing, not a third error direction**. 18 molecules at corpus scale (not the 27 assumed),
  hydrogen-only in 18/18 — and **8 of them re-encode byte-identically to their input**, so no string
  comparison of any kind can separate a structure carrying two extra hydrogens from the original.
  `XAKCAP_comp_0` defeats four instruments at once: scored string PASS, honest string PASS (the same
  string), key EQUAL, `coordination` INTACT, atom count 61 ≠ 63.
- Vendored metric fixtures + tests: `tests/unit/test_honest_score.py`,
  `tests/unit/test_atom_count_gate.py`.

### Why the drop is trustworthy
The honest arm could itself be wrong — an encoder that mis-perceives generated geometries would
manufacture exactly this result. It was checked against `report["coordination"]`, which reads
distances only and consults neither bond graph:

- **control** (3595 `byte_exact` in both arms): 1.3% flagged — inside its 3.7% false-alarm band;
- **moved** (428 pass → `structural`/`facmer_divergent`): **64.7% flagged — a 50× enrichment**;
- the two agree on **mechanism**, not just count: `contacts_lost ÷ sites_lost` lands on the integers
  1–7, and 95.2% of the ratio-above-1 population carries a haptic token — because an η⁵-Cp is five
  *contacts* but one coordination *site*. Two instruments built from different data agree on which
  ligand left and how many atoms it was bound through.
- The 99 molecules `coordination` calls clean are its documented blindness (ligand-interior
  connectivity, bond order, same-CN rearrangement — the `OGARAP` class), not the honest arm
  misfiring.

This refutes the previously recorded suspicion that a ~19× `structural` inflation on re-encode was
an artifact of the method. It is real, and two thirds of it is independently confirmed.

Full analysis: `docs/agentic-notes/v0.4.8/HONEST_BASELINE_v0.4.8.md` and
`docs/agentic-notes/v0.4.8/ATOM_COUNT_GATE_v0.4.8.md`.

## [0.4.7] - 2026-07-26

> **A release of measured negatives.** Five swimlanes ran; four of them ended by refuting the thing
> they were built to ship, and the fifth refuted its own sibling. Every lever here is **default-OFF**
> and the `results-v0.4.6-sweep` bucket report is **byte-for-byte unchanged** — `byte_exact`
> 4140 / 82.80 %, `> 30 s` 994 / 19.88 % — which is the whole point: nothing shipped to the default
> path, and the evidence that nothing did is on the record.

### Added
- **`OIN_ATTACH_CHECK`** (default **OFF**) — `OIN_ACCEPT_SCORED`'s missing safety condition:
  *accept the first conformer the score credits **that still has its ligands attached***.
  Coordinate-only, 7–81 ms/conformer against the 48–57 s strict test it replaces. It **never reads a
  bond object** as evidence of attachment — a detached ligand keeps its bond, so a `GetBonds()`-based
  check certifies exactly the defect it exists to catch. Falsified before it was built: separates
  **7 of the 8** known independent regressions with **0 false positives** over 22 round-tripping
  conformers, and both predicates the promote lane proposed fail (count-based wrongly rejects 11 of
  22, set-based 3). Known residual shipped honestly: `POVPIA` is not caught, so this is 7/8, never
  8/8. **Pairing rule:** against the bare lever it is strictly one-directional — 17 fixes, **zero**
  regressions, severe clashes 14 → 5 — so nobody should run `OIN_ACCEPT_SCORED` without it.
- **The two-arm v0.4.7 regression gate** (`tools/gate_v047.sh`) — ARM 1 encode-only over 61 fixtures,
  ARM 2 a full round trip over a frozen 100-molecule slow cohort, both against committed goldens
  carrying a `MANIFEST_SHA256`. It can tell a notation change from a compute change, which no
  existing instrument could do. The cohort was rebuilt mid-lane: the first construction intersected
  two result dirs and drew from a population that was not the frozen seed-42 5 k.
- **`docs/agentic-notes/v0.4.7/`** — five lane reports: `ACCEPT_SCORED`, `ATTACH_CHECK`, `COHORT`,
  `ENCODE_FLOOR`, `BORON_GEN_CEILING`.

### Changed
- **`docs/` is split by audience and the root is closed.** Of 76 tracked files, 72 were session
  artifacts with four genuine product docs lost among them. Everything measured, tried or refuted
  moved to `docs/agentic-notes/<release>/` (via `git mv`, so `log --follow` works); the root now
  holds exactly `README`, `OPTIMIZERS`, `GENERATION_PIPELINE`, `KNOWN_LIMITATIONS`. Enforced by
  `tools/check_docs_layout.sh` from `pre-commit`, with the rule in `.agents/rules/docs-layout.md`.
  ⚠ The guard does **not** run on merge commits — `pre-commit` is skipped for merges — so a merge is
  the one path into the repo that walks past it. All five v0.4.7 lane docs arrived that way and were
  routed by hand.
- **`OIN_ACCEPT_SCORED`: DO NOT PROMOTE**, superseding this lever's own earlier "promote" reading.
  The gate that recommended it was **circular** — `passed` is computed with
  `get_oin_string(gen.mol, coords)`, the same predicate the lever accepts on, so "18/22 both arms,
  zero regressions" could not detect what dropping the strict step costs. Measured with a genuinely
  independent arm: indep **15/20 → 7/20, 8 regressions, 0 fixes**, one-way; **6 of the 8 lose haptic
  coordination outright**, with the metal geometry tag degrading in lockstep
  (`[Ru_TET]`→`[Ru_TPL]`, `[Zr_TET]`→`[Zr_LIN]`). At n=100 it costs +28 vdW clashes and takes severe
  clashes 5 → 14. **Byte-identical notation, changed geometry** — the cost is invisible to the very
  metric that would police it.
- **The encode floor is three cost regimes, not one.** `XIRMER`'s 20-minute encode is **3 forked
  resonance timeouts at 95.8 % of the fork budget each** — the CPU-limit backstop firing, not slow
  arithmetic, so tuning the solver would never have touched it. The `AC2BO` memo does hit now
  (LEAD 3); LEAD 2 does not matter. Sub-cap valence search no longer materialises a product it never
  reads, pinned by `tests/unit/test_valence_order.py`.
- **The boron generation ceiling, priced at the generator.** Promoting `OIN_BORON_CAGE` in v0.4.6
  moved most of that class from failing *instantly* to failing *slowly*: 40 of 48 burn a large,
  uncapped amount of compute — a full-complex PuLP/CBC bond-order-and-charge solve, re-primed per
  dummy-metal option, on a cage vertex that has no 2c-2e Lewis structure. ≈ **2.1 CPU-h** per 5,000-
  molecule sweep. The internal generation budget is again confirmed **advisory**: molecules ran to
  **3.7×** the requested cap.

### Fixed
- **`OIN_BORON_GEN_FASTFAIL`'s discriminator was refuted, and the lever corrected before it could do
  harm.** The lane shipped a `{LIN}` geometry exclusion and a docstring naming TET among geometries
  that "failed 100 % of the time"; `ULODUU_comp_0` is `[Zr_TET]` and **does** generate (61.8 s). The
  lane's sweep capped at 30 s, so a 61.8 s success read there as a cap-burner — **the class boundary
  moves with the compute you give it**. Cross-tabulated, **every geometry with a success also has
  failures** (LIN 1/3, TET 1/4), so geometry separates nothing — and it was the last discriminator
  standing after hapticity, size and denticity. Safe set widened to `{LIN, TET}`, priced exactly: it
  gives up **4** of 25 recoverable cap-burners (all TET, 60–64 s), while the LIN exclusion costs
  **nothing** (its 3 failures already fail in 0.00–0.01 s). Promotion now requires a *mechanism*, not
  more correlation. A real `ULODUU_comp_0` fixture pins it.
- **16 function-level imports across the merged lanes still named `utils/xyz2mol{,_local}.py`.**
  Because they are function-level, nothing failed at import time — the modules loaded clean and would
  have raised only when the code ran.

### Added
- **`docs/agentic-notes/v0.4.4/DIRECT_DG_VALIDATION.md`**: scale-validation section recording the v0.4.0-vs-v0.4.4
  regression sweep (3,917 molecules = 2,917 v0.4.0 failures + 1,000 seed-42 successes, full
  quality) — **11 regressions / 1,092 fixes / net +1,081** round-trip-OK, zero correctness
  regressions on the 1,000-success guard, confirming direct-DG-as-default at dataset scale.
- **`docs/KNOWN_LIMITATIONS.md`**: note that full-quality generation can time out (300 s) on 11
  medium-large molecules that quick mode round-trips in <30 s — a generation compute-time regime
  (full pool + direct-DG), not a wrong-answer/notation regression; root cause not yet isolated.
- **`docs/agentic-notes/v0.4.4/ACCURACY_v0.4.4.md`**: round-trip accuracy report — the measured 3,917-molecule
  regression sweep (11 regressions / 1,092 fixes) and a full-set pass-rate **projection** of
  **~92% (±~1)**, up from the v0.4.0 baseline of 88.4%, with the projection caveats stated.
- **`src/oinsmiles/oin/coordination.py`** — coordination-integrity check, and the first thing in a
  report that can see a **detached ligand**. The round-trip metric scores via
  `get_oin_string(gen_result.mol, …)`, i.e. the *generator's own bond graph*, so a Cp ring that
  drifted 0.85 Å off the metal is still "bonded" there and the molecule scores a pass. Geometric
  only (distances, ~10 ms, consults neither bond graph); validated **90.2 % recall / 3.7 % false
  alarm** against 633 molecules with independent ground truth. Wired into
  `tools/test_dataset_roundtrip.py` as `report["coordination"]` — **always on, never a gate.**
- **`docs/agentic-notes/v0.4.6/METRIC_FALSE_POSITIVES.md`** — the reported accuracy is not the accuracy, measured in both
  directions over 633 scored successes and 302 reported failures of `results-v0.4.5-rebaseline`.
  One root cause, two opposite errors: `gen_result.mol` **asserts bonds the geometry lacks**
  (61/633 = **9.6 %** false positives, **28.1 % of haptic** molecules — FIYHUT ships both Cp rings
  off the Fe and passes) and **lacks stereo the geometry has** (8/302 = 2.6 % false negatives —
  YOSXIP's `[S@]{5}` sulfoxide flattens to `S{5}` and fails). Net **~5.7 points over-stated**.
- **`docs/agentic-notes/v0.4.6/SWEEP_v0.4.6_5K.md`** — the seed-42 5,000-molecule sweep: **round-trip 4660/5000 = 93.2 %**,
  `byte_exact` **82.80 %**. Failure attribution shows **78.8 % of failures never test the notation**
  (median failing molecule 300.2 s against a 300 s budget), so the notation-attributable gap is
  **~1.1 %** and the `<30 s` and `100 %` goals are one goal: generator compute. PASS 2's recovery
  ladder recovered **zero**.
- **`OIN_INDEP_SCORE`** (default **OFF**) — records the honest round-trip verdict beside the scored
  one as `smiles_2_indep` / `indep_key_match`, via a full `XYZToSMILES().convert()` of the generated
  geometry. Fixes both error directions with one call. Off because it costs a second encode
  (~0.4–1.5 s/molecule) and because switching the *scored* predicate would move ~53 molecules in one
  step — indistinguishable from a regression.
- **`OIN_ACCEPT_SCORED`** (default **OFF**) — pool acceptance on the predicate the score uses,
  skipping the stricter independent re-perception. Median **16.0 s → 3.6 s**, molecules over 30 s
  **10 → 3**. Held off on measurement, not caution: re-running the sweep's 340 failures with it on
  "recovered" 90, but `report["coordination"]` shows those 90 are **60 DEGRADED / 21 BOUNDARY /
  9 INTACT** — **+90 → 95.0 % is really +9 → 93.4 %**, so a global flip manufactures 60 phantom
  passes no existing gate could detect. Recommendation: promote **with scope**, throughput runs only.
- Tools: `haptic_false_positive.py`, `probe_accept_gap.py`, `summarize_accept_gap.py`,
  `ab_accept_scored.py`, `profile_eta.py`, `boron_gen_time.py`.

### Changed
- **`docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md`** — its 34/34 "round-trip" is corrected to **notation-level**: all
  nine checks are encoder-side and none invokes the generator. The pipeline arm, 33 of 34 measured:
  only **2 generate** a structure, 25 burn the whole cap (**~2.1 CPU-h wasted per 5k sweep**), 6 fail
  instantly. So `OIN_BORON_CAGE` moved most of this class from failing *instantly* to failing
  *slowly* — a cost missing from its original pricing, now recorded in `levers.py`.
- **`OIN3DGenerator(timeout=…)` is documented as ADVISORY, not a bound** — 60 s requested,
  60.0–172.8 s spent. Only the harness's per-molecule SIGKILL enforces a budget, so any timing taken
  without it understates the tail.
- **Renamed `utils/xyz2mol_local.py` → `utils/perception_core.py` and `utils/xyz2mol.py` →
  `utils/perception_tmc.py`** (no behavior change): the `_local` suffix read as "unmodified vendored
  copy," but the file is a hard fork — roughly a third of its functions are project-owned, and it
  imports `oin.levers.lever_enabled` back from the package, something a genuinely-vendored file
  wouldn't do. New names separate the distance→bond-order core (`perception_core`) from the
  TMC-specific driver built on top of it (`perception_tmc`). The upstream `xyz2mol()` function name
  is retained unchanged inside `perception_core.py`.

### Fixed
- `clash.mol_clash_count()` returns 0 for any bare RDKit `Mol` (no `.atom_list`), so a
  structure-quality gate built on it reads 0 everywhere and certifies without measuring. Quality is
  now computed with `vdw_clash_count(positions, Z)` plus the continuous `worst_overlap`.
- The coordination check's boundary band was **one-sided** — it saw contacts just past the cutoff but
  not those *held* by a hair. Two-sided, the count of false positives passing with no signal drops
  **6 → 2**, with the primary verdict's recall and false-alarm rate unchanged.

## [0.4.6] - 2026-07-26

### Added
- **`OIN_BORON_CAGE` promoted to default-ON** — deltahedral borane/carborane cages now encode.
  Measured on the 936-molecule re-baseline: 34 of the 36 `XYZToSMILES failed` molecules are
  electron-deficient boron clusters, and the lever takes that population from **0/36 encoding to
  34/36** (0.2–4.2 s each). The boron lane separately measured 48/48 round-tripping.
  ⚠ It moves **14 molecules from scored-passing to failing**, which is correct — they passed while
  describing the wrong graph (`VEJXOZ` invents a C=B double bond) — but it trades 14 silent false
  positives for 14 loud honest failures, so a headline pass rate can move either way.
- **Metal-centred Δ/Λ descriptor (Y1 blind spot P1)** — `src/oinsmiles/oin/metal_config.py`, a
  complete pipeline where none existed: the lane began with 0/150 molecules emitting any metal
  stereo tag because no descriptor existed. Chelate-aware symmetry test → trailing `|mc:±|` sidecar
  behind `OIN_EMIT_METAL_CONFIG` (**default OFF**) → key fold with the un-fold obligation recorded →
  generator reproduction of the requested helicity. ZUMNEC emits and inverts under reflection;
  square-planar JEGKOW correctly emits nothing.
- **`test_levers::TestNoTestUnsetsAPromotedLever`** — a lint for the "unset means off" trap, which
  cost 23 test failures across two promotions. It found three further instances on its first run,
  one of which was passing vacuously and invisible to 838 green tests.

### Fixed
- **`OIN_BORON_CAGE` valence-bypass scope.** `_parse_fragment`'s cage rung skips
  `SANITIZE_PROPERTIES`; gated on the lever alone it applied that bypass to *every* fragment, so
  `C#O` parsed instead of reaching the `RAW:` fallback — and CO is among the commonest ligands in
  transition-metal chemistry. Now scoped to boron-containing fragments. Verified over all **1,194**
  distinct fragment bodies the corpus emits: 56 differ, all 56 contain boron, **0 boron-free
  affected**.
- **`canonical_body_emit` is H-faithful.** Both of its `MolToSmiles` writes now route through
  `h_faithful_smiles`, so `OIN_CANONICAL_BODY` no longer silently discards the `OIN_H_FAITHFUL`
  repair. Byte-identical on all 61 fixtures with that lever off. ⚠ No measured accuracy benefit —
  the `atom_count` class is unmoved (8/45 both arms) and is not a serialization defect.
- **`MetalloGen failed` diagnostics** — the message reported `m-SMILES None` for every OIN-direct
  failure, because `msmiles` is only populated on the fallback path. It now names the assembly path
  actually used, plus the OIN, pool width and timeout.

### Documented (measured negative results)
- `docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md` records **eight refuted hypotheses**, each with the measurement
  that killed it: P3 tag restoration through the reparse (rewrites RIFGUJ's *relative* ring-carbon
  stereo); "the accuracy gap is mostly compute" (75% of timeouts hide real failures, not latent
  passes); a donor-cut rule for hydrogen (`dH` spans −36…+14, matches the bare-donor count in 4/45);
  eta incremental pool widening (the fill loop already short-circuits via `accept_fn`); and three
  successive Δ/Λ formulations.
- **The `<30 s` tail is attributed**: eta-ring winding appears in **63.5%** of the >30 s tail vs
  19.8% of the fast set. The cost is a low *acceptance rate*, not pool bookkeeping, so the remedy is
  raising the probability of generating the requested ring face — construction over selection, which
  this project carries three prior negative results for.

## [0.4.5] - 2026-07-26

### Added
- **Canonical OIN-SMILES**: the *emitted* string is now canonical, not just the comparison key. Six
  levers ship ON via the new single-source registry `src/oinsmiles/oin/levers.py`:
  `OIN_CANONICAL_BODY`, `OIN_CANONICAL_PERCEPTION`, `OIN_CANONICAL_SLOTS`,
  `OIN_CANONICAL_ETA_WINDING`, `OIN_STABLE_METAL_AC`, `OIN_STABLE_STEREO`.
  What made them safe to promote together: each **repairs a renumbered presentation without
  rewriting the canonical answer**, which is why the corpus shows no churn.
- `src/oinsmiles/oin/canonical_slots.py` — lex-min colored-vertex signature over the proper-rotation
  group, exporting `canonical_slot_permutation()`. Also unifies the rotation group, fixing PBP from
  2 of 10 proper rotations to 10.
- `docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md`, plus 15 per-lane measurement documents.

### Fixed
- **Encoder instability under pure atom renumbering** (unplanned Lane 8): 13% of molecules emitted a
  *different absolute stereochemistry* when the input atoms were reordered. No prior instrument
  could see it. Fixed by `OIN_STABLE_STEREO`, which re-derives fragment tags from the parent
  geometry rather than translating a parity relative to a destroyed neighbour order.
- **Inverted CIP goldens** for `PdCl2-RR-BDPP` / `PdCl2-RR-BDNN`, wrong for four months. The test
  that "verified" them ran `rdCIPLabeler` on a SMILES reparsed from the encoder's *own* output —
  `rdCIPLabeler` converts a parity tag into a label without checking it, so an inverted tag was
  self-consistent and passed. Ground truth from `AssignStereochemistryFrom3D` is (R,R), agreeing
  with the fixtures' own `(2R,4R)` names.
- The lever registry closes a live trap: `bool(os.environ.get("X"))` made `X=0` *enable* X.
  `OIN_BORON_CAGE` alone had five sites on that spelling.

### Measured
- Byte-stability under rotation/renumbering **58.1% → 69.6%**; comparison-key instability
  **60 → 16 molecules**.
- Re-baseline over 936 molecules: **145 of 436** previously-failing molecules fixed (33.3%);
  of 500 previously-passing guards, all 11 apparent regressions are `TimeoutException exceeded 300s`
  against an 1800 s baseline ⇒ **zero correctness regressions**.
- Suite: 837 tests OK.

### Known limitations
- **P3 (metal-bound 2° amine) is not usable in the shipped default** — `OIN_CANONICAL_BODY`'s
  reparse clears the `[N@]` it stamps. The obvious fix was tried and is *measurably wrong*; see
  `levers.py::_HELD_OFF`.
- `OIN_EMIT_AXIAL`'s promotion evidence needs re-measuring under canonical perception, which changes
  the hindered-axis count (YESKOZ 2 → 1).

## [0.4.4] - 2026-07-23

An **accuracy + measurement** release, developed as six parallel worktree swimlanes (SL0–SL5)
squash-merged to `main`, then a promote pass. Where v0.4.3 attacked structure *quality*, v0.4.4
attacks round-trip *string accuracy* and the instrument used to measure it. The centerpiece is a
**fac/mer-aware canonical round-trip key** (SL0) that distinguishes a genuine *fac*↔*mer* miss
from benign slot-relabeling — turning "did it round-trip" from a brittle byte-compare into a
symmetry-aware hash, and decomposing the accuracy gap into byte_exact / key_equal /
facmer_divergent / structural / hard_fail / encode_fail buckets
(`tools/roundtrip_bucket_report.py`). Two levers earned promotion to default on the worst-cohort
A/B: **early-exit conformer acceptance** (accept the first conformer that re-encodes to the key —
byte-exact 44.7% → 60.5%, ~5× faster, zero regressions) and **OIN-direct assembly** (build the
internal complex straight from the parsed OIN, retiring the lossy m-SMILES bridge to a fallback —
accuracy-neutral, and metal `@SPn` chirality now rides through the representation). Three other
levers (rigid η-winding construction, greedy placement, stretched-bond metric) were measured,
found net-negative or neutral, and **kept opt-in** — a third confirmation that *selection beats
construction*. No OIN format-version change (still v3.7 inline). Full unit suite **488 → 551 OK /
3 skip**; FF-path goldens round-trip byte-identically under the new defaults. See
`docs/GENERATION_PIPELINE.md` for the full default pipeline and `docs/agentic-notes/v0.4.4/DIRECT_DG_VALIDATION.md` for
the direct-assembly decision record.

### Added
- **fac/mer-aware canonical round-trip key** (`oin/compare.py::canonical_roundtrip_key`): a
  symmetry-aware slot canonicalization that distinguishes *fac* from *mer* (and cis/trans) while
  folding benign slot-relabeling, so a round-trip "pass" measures isomer identity rather than byte
  equality. Paired with **`tools/roundtrip_bucket_report.py`**, which classifies every round-trip
  into byte_exact / key_equal / facmer_divergent / structural / hard_fail / encode_fail. (SL0.)
- **`docs/GENERATION_PIPELINE.md`**: a sectioned description of the full default OIN → 3D pipeline
  (parse → direct assembly → distance-geometry embed → winding search → early-exit acceptance →
  geometry selection) and how each stage is validated.
- **Opt-in generation levers** (all OFF by default; success path byte-identical): direct-assembly
  rigid η-winding construction (`OIN_DIRECT_ASSEMBLY` / `ff_params={"oin_direct": True}`, SL2),
  difficulty-ordered collision-aware Kabsch greedy placement (SL3), a stretched-bond acceptance
  metric (SL1), and a no-acceptance-progress embed cutoff (`OIN_EMBED_NO_PROGRESS`, SL4).

### Changed
- **Early-exit conformer acceptance is now ON by default** (`generator3d/__init__.py`,
  `generation/metallogen_adapter.py`; opt out with `OIN_EARLY_EXIT=0` or
  `ff_params={"early_exit": False}`). Generation stops at the first conformer that independently
  re-encodes to the fac/mer key instead of exhausting the pool. On the worst-cohort A/B this lifted
  byte-exact round-trip **44.7% → 60.5%** and key-match **55.3% → 73.7%** with **zero regressions**,
  ~5× faster. (SL1 promote.)
- **OIN-direct assembly is now the default generation path** (`generation/metallogen_adapter.py`;
  opt out with `OIN_DIRECT_DG=0` or `ff_params={"direct_dg": False}`). The internal `MetalComplex`
  is built directly from the parsed OIN via `om.get_om_from_parsed`, retiring the winding-lossy
  m-SMILES bridge (`convert_parsed_to_msmiles`) to a **fallback** used only if direct assembly
  raises. Proven accuracy-neutral twice (0 regressions / 0 gains / identical buckets on a
  38-molecule stratified A/B and a confirmation A/B); metal `@SPn` chirality now survives into 3D
  generation. (direct-dg promote; `docs/agentic-notes/v0.4.4/DIRECT_DG_VALIDATION.md`.)
- **The round-trip harness demotes RMSD to a diagnostic** (`tools/test_dataset_roundtrip.py`): a
  string-exact round trip that only exceeds the coordination-sphere RMSD gate is now a **success**
  carrying an `rmsd_over_gate` diagnostic (RMSD is only ~0.22-correlated with geometric quality),
  reclaiming 51 of 306 hard-fail molecules (28 directly attributable to the demote) with no
  regression to the passing middle. (SL4; `docs/agentic-notes/v0.4.4/RELIABILITY_v0.4.4_SL4.md`.)

### Fixed
- **Electron-deficient boron clusters now fail with a typed, classified error** (`utils/perception_tmc.py`,
  `core/translator.py`): a carborane / closo-nido borane cage (≥3 borons with a B–B bond) that RDKit
  cannot perceive into a Lewis structure now raises a typed **`OINEncodeError`** naming the cause
  instead of returning `None`, so a known encode ceiling is distinguishable from an unexpected
  failure. `OINEncodeError` subclasses `ValueError`, so existing handlers are unaffected. Classifies
  34 of the 48 capstone `encode_fail` molecules. (SL5; `docs/agentic-notes/v0.4.4/ENCODER_ROBUSTNESS_v0.4.4_SL5.md`.)
- **xyz2mol perception hangs on large conjugated ligands bounded and recovered**
  (`utils/perception_core.py`, `utils/perception_tmc.py`): the `AC2BO` valence-combination sort is now capped
  (`_VALENCE_COMBO_CAP`), and the `ResonanceMolSupplier` enumeration runs in a **forked,
  CPU-time-bounded child** (`RLIMIT_CPU`) — a completer returns its true resonance form
  byte-identically, while a genuine hang is killed and degrades to the single form, recovering
  previously-unencodable molecules (e.g. `BENVOG`, `HUCNAU`) without changing any currently-encodable
  OIN. (SL5.)

## [0.4.3] - 2026-07-21

A **structure-quality + conformer-invariance** release, developed as nine parallel worktree
sessions (A0–A5, B1, B2, C1), each squash-merged to `main`. Where v0.4.2 attacked round-trip
*string* accuracy, v0.4.3 attacks the *geometry the round-trip gate cannot see*: the elimination
study (`docs/agentic-notes/v0.4.3/FALSIFICATION_v0.4.3_ELIMINATION.md`) found the round-trip gate is geometry-blind
(`ρ(coord_rmsd, full_divergence)=0.22`) and the generator was routinely shipping
validation-failing structures — **53% of generated structures carried a van der Waals clash vs
5% for real crystals**. The headline fix is a whole-complex vdW clash acceptance term, now on by
default. On a pinned worst-cohort A/B the generated clash fraction dropped **92.5% → 5.1%** (near
the 2.5% real-crystal reference) with round-trip fidelity held (92.5% → 90.0%, zero donor
ejections). An early 100-molecule validation sample (80:20 split) gained **31+ round-trip
successes over v0.4.2 with zero regressions**; a full ~2,900-molecule sweep is in progress and
will be reported when complete. No OIN format-version change (still v3.7 inline). Full unit suite
**435 → 488 OK**; FF-path goldens still round-trip under the new defaults.

### Added
- **Whole-complex vdW clash acceptance term** (`generator3d/clash.py`, `generator3d/embed.py`,
  `generator3d/clean_geometry.py`, `generator3d/__init__.py`): a van der Waals steric-overlap
  gate and candidate score applied at the embed acceptance gate, the UFF pool ranking, geometry
  selection, and the FF-clean scan (previously only atomic *fusion* was rejected). **On by
  default** (`clash.VDW_ACCEPTANCE_ENABLED`; opt out per run with `OIN_VDW_ACCEPTANCE=0`). Drops
  the generated clash fraction 92.5% → 5.1% on the worst-cohort sample with no donor ejections.
  (Sessions A3 + A5.)
- **Kabsch/Umeyama rigid-placement embed (`option=3`)** (`generator3d/embed.py`): builds each
  ligand independently and rigidly places it onto the ideal coordination vectors, reflection-
  guarded so chelate handedness never flips. Opt-in via `ff_params={"use_kabsch": True}` (or
  `{"kabsch_only": True}` to isolate it for A/B). Kept opt-in — it gave no clash benefit over the
  vdW term alone. (Session A4.)
- **Electronic (ligand-field) geometry prior** (`utils/oin_aligner.py`, `utils/perception_tmc.py`):
  resolves ambiguous low-coordination spheres (square-planar vs tetrahedral vs square-pyramidal)
  using the metal d-electron count when the geometric RMSD is a near-tie (e.g. d8 Ni), so
  distinct conformers of one isomer converge to a single OIN string. (Session B1.)
- **Conformer-invariance integration tests** (`tests/integration/test_conformer_convergence.py`,
  `test_isomer_divergence.py`, `test_conformer_invariance.py`) and a 30-molecule stratified
  conformer test set (`tests/fixtures/conformer_set/`): assert that all conformers of one isomer
  encode to the same canonical OIN string, and that a *different* isomer (SPL↔TET, `{n>}`↔`{n<}`,
  E/Z, cis↔trans slot order) encodes to a *different* string. (Sessions B1 + B2.)
- **Optional CREST conformer cross-check** (`tools/conformer_invariance_crest.py`,
  `tools/install_crest.sh`, `tools/run_conformer_crest_sweep.sh`): builds a real conformer
  ensemble per structure (g-xTB pre-opt + per-conformer re-opt) and reports whether they collapse
  to one OIN string. CREST is an optional external binary; the workflow skips gracefully when it
  is absent. (Session C1.)
- **Structure-distortion tooling** (`tools/structure_distortion_report.py`,
  `tools/compare_failures.py`): scores generated structures on a geometry/graph/reference MPO with
  a vdW-clash headline — the primary quality metric for this release, since round-trip pass is
  geometry-blind. (Session A0.)

### Changed
- **Default 3D-generation optimizer is now `ff` (fast FF-only + vdW acceptance), not g-xTB.** The
  A5 A/B found FF + the vdW term gives the lowest clash fraction and the best round-trip fidelity,
  deterministically and fast; g-xTB (`optimizer="xtb"`, most geometry-accurate but slower) and
  MACE remain opt-in for callers who want physical refinement. Applies to `OIN3DGenerator`,
  `SMILESToXYZ`, and the `oin-smiles oin2xyz` CLI (`--optimizer` default `g-xtb` → `ff`).
  (Session A5.)
- **Weak-field high-spin multiplicity** (`generation/om.py`): the generator no longer forces
  every complex to singlet/doublet — the correct spin state is assigned for high-spin-capable
  metals (inert on the FF path, consumed by the opt-in g-xTB/MACE relax). (Session A1.)
- **Project `xtb` binary now reachable under the run environment** (`generator3d/ml_optimizer.py`):
  the `.venv/bin/xtb` build is honored via `OIN_XTB_BIN` / path resolution, so opt-in
  `optimizer="xtb"` invokes g-xTB instead of silently falling back to FF. (Session A1.)

### Fixed
- **Silent "no conformers" on structurally-impossible complexes** (`generator3d/__init__.py`): a
  blanket `except Exception` that swallowed structural embed failures (e.g. over-valent dative
  donors) now surfaces a typed `StructuralAssemblyError`, so the failure is diagnosable instead of
  reported as a generic conformer miss. (Session A2.)

## [0.4.2] - 2026-07-15

The tmCAT/tmPHOTO round-trip **accuracy** wave (parallel worktree phases S1/S3/S5/S6a/S6b/S7 plus a
docs phase), staged on the `release/v0.4.2` integration branch off v0.4.1. Where v0.3.7 closed the
residual failure classes on a ~2,600-complex sample, the continuous accumulator has since swept the
**full 25,197-complex** corpus (the sweep is now complete); this wave attacks what that larger sweep
re-surfaced. It delivers three things: **107 real encoder/generator accuracy fixes** (confirmed
against a matched v0.4.1 `--quick` control — see the quick-mode A/B below); a **measurement
re-framing** — most
of the largest "failure" buckets (`geometry_or_fragment_change`, `winding_flip`, `timeout`,
`high_rmsd`) turn out to be `--quick`/FF-only **harness artifacts or misfiled classes**, not accuracy
defects; and **honesty tooling + documentation** that stops the backlog conflating the two.

No single headline pass-% is quoted, deliberately: the accumulator is mixed-provenance (multiple
commits) and runs `--quick` (a different, pool-of-2 generator), so a global percentage is not a valid
floor. Accuracy is stated as **named per-molecule round-trip flips** against the single-commit
`c7edeeb6` floor (a set of ~5,960 molecule IDs that pass on the baseline, `spec/handoffs/v0.4.2/`).
No OIN format-version change (still v3.7 inline); one new CN-9 geometry code is added. Full unit
suite **409 → 435 OK**; FF-path geometry is byte-identical to v0.4.1 for every non-targeted molecule.
See `docs/agentic-notes/v0.4.2/ACCURACY_v0.4.2.md` for the per-class before→after and `docs/KNOWN_LIMITATIONS.md` for what
remains out of scope.

### Fixed
- **Bare 0-H haptic ring carbons re-protonated on generation (S1)** (`generation/metallogen_adapter.py`):
  a bare η ipso / ring-fusion aromatic carbon (`c{n}`, 0 H) was over-protonated because
  `get_ligand_from_smiles`'s `Chem.AddHs(explicitOnly=False)` adds a phantom H to a carbon already
  bonded to three heavy neighbours — the H count is lost crossing into MetalloGen (`get_ace_mol_from_rd_mol`
  copies atomic number / charge / bond order, **not** H counts or `NoImplicit`). `convert_parsed_to_msmiles`
  now locks a haptic ring carbon whose perceived H count is 0 (`SetNumExplicitHs(0)` + `SetNoImplicit(True)`),
  H==0-only so genuine `[cH]` haptic carbons and correctly-perceived substituted carbons are byte-identical.
  Flips `ARONEA`, `BOXJUU`, `CAHZIX`, `COTXAM`, `HINNOH`, `NORTAP`; Cp*/bis-indenyl/ferrocene passers
  unchanged. Guarded by `tests/unit/test_haptic_carbon_hcount.py`.
- **Two XYZ→OIN encoder crash paths recovered (S3)** (`utils/perception_tmc.py`, `utils/aromaticity.py`):
  (a) `get_lig_mol` early-returned `None` for large saturated polyamine/phosphine cages whose extended
  Hückel charge is several electrons off — it now falls through to the existing `_rescue_unusable_perception`
  charge sweep, so `IROXET`, `SUNXAB`, `XEVMAN` encode and round-trip. (b) `kekulize_safe_sanitize`'s
  "no quinoid ring" branch now retries a fresh aromaticity re-perception before raising, clearing the
  stale-aromatic-flag kekulize crash on `JOTJEK`, `TIYWUV`, `ZENZAW` (their canonical key now matches;
  a residual `[CH]`-radical H-fill count is S1's domain). Both changes are additive — they only alter
  a code path that is *already* crashing, so no passer can regress. Guarded by
  `tests/unit/test_encode_crash_recovery.py`.
- **Free monodentate-arm C=N/C=C E/Z dropped by the generator (S6a)** (`generator3d/embed.py`,
  `generator3d/ligand.py`): a documented latent bug blocked the fix — `_apply_double_bond_stereo`
  restored a PuLP-demoted bond to `DOUBLE` to enforce its E/Z but never adjusted formal charge, so
  enforcing a free arm (e.g. `AFECIZ`'s C=N, which PuLP wants single-with-charged-N) produced a
  4-valent neutral N, `SanitizeMol` rejected it, and *every* `ff_clean` raised. A new
  `_charge_fix_promotion` bumps an over-filled endpoint's formal charge, gated by the existing valence
  probe (safe because the round trip re-encodes from geometry — the generator's internal charges never
  reach the output). With the crash gone, `near_donor` is narrowed from the broad donor-neighbour proxy
  to the encoder's chelate-ring predicate (a new `_chelate_locked_atoms` virtual-metal ring test), so
  free-arm E/Z is enforced while ring-locked bonds stay suppressed. Flips `AHAZOZ` (nitrone
  `/C=[N+](\[O-])`), `AFECIZ`, `XIZXAG`, `AYUYIE`, `BOCGEH`; `AFECIZ`'s `ff_clean` no longer exhausts
  the 250-attempt budget. Guarded by `tests/unit/test_chelate_locked_ez.py` (the pending-charge-fix
  guard is flipped to assert enforcement; +2 charge-aware promotion tests).
- **sp3 / heteroatom atom-stereo `@`/`@@` disagreements (S6b)** (`core/chirality.py`,
  `generation/metallogen_adapter.py`): both encode paths funnel through
  `ChiralityRecoveryUtility.recover`, so the fixes are symmetric. (a) A spurious high-coordination
  donor-S tag (`[S@SP3]`/`[S@SP1]`/`[S@TB9H]`) from `AssignStereochemistryFrom3D` is now cleared when
  the S is not a genuine stereocentre on the metal-free fragment (`FindMolChiralCenters(useLegacyImplementation=False)`
  — a divalent thioether S is absent, a genuine sulfonimidoyl S(VI) is kept): `BAZMOH`, `HUGSEI`,
  `LUSKIV`, `YUMPIH` round-trip and `CIDDAU`'s `@` is resolved. (b) A metal-stripped donor becomes a
  radical (`[O]`, carbene/alkene C) that `rdCIPLabeler` refuses to rank; the reparse now **fills the
  open valence with H** (skipping aromatic atoms, so η-Cp/arene donors are untouched), and both the
  `_OIN_CIPCode_SP3` stamp and the comparison read the label through this one fill-first convention —
  re-orienting sulfonimidoyl S (`JEKQAS`, `REPZUJ`, `ZORCOA`), quaternary N⁺ (`POYJIX`, routed through
  the same degree-4 machinery), and donor-adjacent C (`ORIHUU`, `XILZID`). Guarded by
  `tests/unit/test_heteroatom_atom_chirality.py` (+4). Re-triaged out: `KEBBUO` → S7 (a `--quick`
  timeout), `SEMTOV`/`VEJXOZ` → S5 (borane winding).

### Added
- **CN-9 tricapped-trigonal-prism geometry (`TCT`) (S5)** (`utils/oin_aligner.py`,
  `generation/oin_parser.py`, `generation/metallogen_adapter.py`): a 9-coordinate metal previously had
  no template, so the encoder emitted `g:NON` and generation raised. A `TCT` template added across the
  same 5-site lockstep as v0.3.5's CN-8 `SQA` (encoder `TEMPLATE_SPECS` + candidate branch, the
  load-bearing `oin_parser.py::TEMPLATES` mirror, the generator geo map, and `tests/unit/test_cn9_geometry.py`)
  now classifies it as `[Y_TCT]`. `XERTUK_comp_3` encodes `[Y_TCT]` (was `[Y_NON]`); it still cannot
  embed its 104-atom ligand — a generation-embedding limit, documented, not a geometry gap.
- **Round-trip harness honesty (S7)** (`tools/test_dataset_roundtrip.py`): a `RMSD_GATE` constant;
  every report now stamps `mol_timeout` and `rmsd_gate`, and marks `ff_floor=True` on a `high_rmsd`
  gate reached under FF-only; the run summary prints a `_honesty_breakdown` that separates real defects
  from FF-floor `high_rmsd` and `--quick` timeouts. So the backlog can no longer silently conflate an
  FF/quick artifact with an accuracy failure. Guarded by `tests/unit/test_harness_honesty.py`.
- **Two-result-set A/B comparator (`tools/ab_compare.py`)**: joins a control and a candidate
  round-trip result set on `molecule` and reports per-class `fail→pass` / `pass→fail` / `fail→fail` /
  `pass→pass` transitions with a `fixes = N, regressions = R` headline and explicit ID lists — the
  two-directory comparator the harness lacked. Accepts a results dir or a `summary_roundtrip.json`,
  groups by an optional `molecule→class` map, stdlib-only. Produced the quick-mode A/B below.
- **Milestone backlog snapshots (`tools/milestone_report.py`)**: freezes a full tiered backlog
  (`reports/REPORT_<M>.md` + `case_registry_<M>.json`) each time the accumulator crosses a clean
  1,000-molecule multiple; idempotent via `.v041_milestones.json`, read-only on the harness outputs,
  safe to run alongside a live `--continue` sweep.

### Changed
- `.gitignore` now ignores `/scratchpad/` — throwaway analysis state should never be linted or
  committed (S7).
- **`tools/job_dashboard.py`**: the live FastAPI results dashboard's `get_stats` now deduplicates by
  molecule (keeping the latest attempt) so re-runs no longer double-count, keeping the served totals in
  step with `case_registry.json`.

### Notes — the measurement re-framing (the load-bearing finding of this wave)
- **`geometry_or_fragment_change` (29 on the floor) and `winding_flip` (14) are overwhelmingly
  `--quick` conformer-pool artifacts, not classifier or generator defects (S5).** Re-run under the
  non-quick default generator (seed 42, pool 5): CN4 geometry **15/15**, CN5 **8/11**, winding-face
  **12/14** correct — with the CN4/CN5 code path byte-identical to pristine. A calibrated TET/TPY
  encoder hysteresis was designed and then **dropped**: the feasible margin interval is empty (a 0.05 Å
  margin flips 13 passing molecules; the passer-gap floor is 0.0022 Å). **No classifier, pool, or
  winding code shipped** — the CN-9 template is S5's only code change.
- **The curated bond-length `ENABLED_METALS` expansion is a validated-negative;
  `generator3d/bond_lengths.py` ships pristine (S7).** A paired median-of-deltas coordination-sphere
  RMSD A/B (≥10 seeds, **full** conformer pool — collapsing it fabricates results) with a Zn positive
  control (−0.025 Å, reproducing P4's landed win, so the harness is trusted) showed **Ru +0.009 Å**
  (confirming P4's deliberate exclusion) and W/Mn/Co flat — the swept scale factor already compensates
  the covalent-sum overestimate. No candidate metal cleared the per-metal gate.
- **`--quick`'s 30 s hard-kill mislabels real chemistry (S7).** A full-budget confirmatory sample found
  `timeout` is mostly real valence `no_conformers` masked by the kill (not a pure budget artifact),
  `no_conformers` is essentially all genuine, and `high_rmsd` is the FF-only geometric floor on an
  already-string-correct round-trip.
- **Many class labels were symptoms, not root causes.** S3 re-triaged its 22 goldens and refuted all
  four of its class labels (the `macrocycle_perception`/`garbled_aromatic` rows were boron atom-stereo,
  geometry, E/Z, and generator H-fill); the docs phase refuted the S1 handoff's `AJIJUY` nitride/ammine
  example (it contains no nitrogen and re-encodes byte-identically), found the nitride/ammine notation
  ambiguity currently **latent** (no clean dataset trigger), and found every `gen_exception_other` row
  is an `UncoordinatedFragmentError` (outer-sphere counterions/solvent).

### Triage / documentation
- `docs/KNOWN_LIMITATIONS.md` reconciled and extended: the nitride/ammine notation ambiguity (latent),
  FF-floor `high_rmsd` + `--quick` timeouts, `UncoordinatedFragmentError` outer-sphere fragments,
  irreducible generator stereochemistry (`QOFTOU` rac/meso, `winding_flip` residual), the S6b
  atom-stereo fixes, and the deferred residuals `JUCCUH` (trivalent-N inversion RDKit clears) and
  `WEDYOU` (macrocyclic multi-P relative configuration). `spec/handoffs/v0.4.2/wontfix-carboranes.md`
  records the carborane class (3-centre-2-electron bonding outside the two-centre `AC2mol` model).
- New `docs/agentic-notes/v0.4.2/ACCURACY_v0.4.2.md`: the per-class before→after table, the measurement-integrity
  rationale, and a **quick-mode A/B confirmation** — re-running all **2,917** v0.4.1 `--quick`
  failures against a matched v0.4.1 control (via `tools/ab_compare.py`, each arm on its own
  checkout's harness+venv) shows v0.4.2 delivers **107 code-attributable round-trip fixes** across the
  deterministic accuracy classes (E/Z 33, string-mismatch 25, atom-stereo 22, donor-H 14, encode-crash
  8, macrocycle 2, winding 2, geometry 1) with **zero deterministic regressions** (the 3 raw
  `pass→fail` are all stochastic `--quick` timeouts).
- New `docs/agentic-notes/v0.4.1/ACCURACY_v0.4.1.md`: the full-corpus failure-mode distribution over the now-complete
  25,197-molecule sweep (88.4% `--quick` round-trip pass / 95.8% accuracy-clean, stated as a screening
  floor).
- Version bumped to 0.4.2 in `pyproject.toml`; the `v0.4.2` git tag is applied when the wave is pushed
  to `origin` (not yet pushed).

## [0.4.1] - 2026-07-14

Consolidation release: three parallel worktrees — a code-quality remediation, a dead-code
purge, and a second performance wave — merged onto the v0.4.0 line. The headline is a **real
`SMILESToXYZ` public reverse API** (previously a dummy-atom stub), plus continued generator
speedups, CI/type-check hardening, `print()` → `logging`, and removal of genuinely dead code.
No OIN format change (still v3.7); FF-path golden geometry is byte-identical to v0.4.0, and a
deterministic 100-molecule dataset A/B reproduced v0.4.0 round-trip outcomes with zero regressions.

### Added
- **Real `SMILESToXYZ` public reverse API** (`core/translator.py`): the exported `SMILESToXYZ`
  — previously a stub that emitted dummy `"X"` atoms and never parsed the OIN — now delegates to
  `OIN3DGenerator`. `convert(oin) -> str` returns the XYZ block (symmetric with
  `XYZToSMILES.convert`); `generate(oin)` returns the full `GeneratedStructure` (XYZ + bonded mol).
  The engine is imported lazily so the forward direction doesn't pay for the generation backend.
  Guarded by a rewritten `tests/unit/test_translator.py` (mock delegation + a real FF end-to-end
  test that fails if dummy `"X"` atoms ever return).
- **CI round-trip smoke guard** (`tests/integration/test_roundtrip_smoke.py`,
  `.github/workflows/ci.yml`): a deterministic OIN → XYZ → OIN guard over the canonical goldens
  (+ transplatin) under `optimizer="ff"`, `seed=42`, so a generation/encoding regression fails CI
  instead of only surfacing in the out-of-band dataset harness.
- **mypy type-checking in CI** (`pyproject.toml`, `.github/workflows/ci.yml`): a first-party
  `[tool.mypy]` scope with vendored-code carve-outs, run as a `Type check (mypy)` CI job.
- **Opt-in embedding/optimization parallelism CLI flags** (`cli.py`): `--embed-threads N` engages
  batched parallel conformer embedding (off by default — it samples conformers differently, so it
  is not byte-identical), and `--optimize-workers N` sets the parallel g-xTB optimize-loop worker count.

### Removed
- **Dead `SMILESToXYZ` support cluster** — `oin/parser.py` (the incomplete `OINParser`),
  `oin/writer.py` (`OINWriter`), `core/graph.py` (`TMCGraph`/`Atom`), and their unit tests. These
  backed the old dummy-atom stub and were unreachable from any working path. The public
  `SMILESToXYZ` name is **kept and now works** (see Added); `engine.py` docstrings that pointed
  users at the old broken path were corrected.
- **Orphaned one-off scripts** from completed waves: `tools/verify_metal_first.py`
  (DirectParser gate spike; feature never shipped), `tools/test_uff_pool_size.py` (P0-era
  UFF pool experiment), `tests/integration/{debug_welrow,reproduce_issue,
  verify_prd_compliance,verify_ir_complexes}.py` and the committed regenerated output dir
  `tests/integration/verification_artifacts_IR_TEST/`.
- **Committed generated artifacts** in `tests/candidate_outputs/`: `molassembler_cisplatin.xyz`,
  `spike_cisplatin.xyz`, `cli_cisplatin.xyz`, `AnsaMetallocene_TiCat1_before_baseline.xyz`.
- **Inert dead code**: stray `pass` + dead string literal in `oin/inline.py`
  (`generate_inline_string`), no-op `__main__` block in `generator3d/om.py`.
- `docs/ETKDG_AROMATIC_FIX.md` (historical doc for the Molassembler backend removed in
  v0.3.7; preserved in git history).

### Fixed
- **CI "Top-level tests" step was a silent no-op** — `unittest discover tests` collected
  0 tests (the only module directly under `tests/` was pytest-style; subdirectories are not
  packages). The step is removed; `tests/test_generator3d.py` was converted to
  `unittest.TestCase` and moved to `tests/unit/` so it actually runs in CI.
- `docs/KNOWN_LIMITATIONS.md` pointed at the removed `generation/molassembler_adapter.py`
  for the bond-length table; now points at `generator3d/bond_lengths.py`.
- **`cli.py` `oin2xyz` `ff_params` typing** — the `--embed-threads`/`--optimize-workers` block
  reassigned `ff_params` from `dict` to `dict | None`, which the new mypy job flagged once both
  the perf CLI flags and the type-check landed together; `ff_params` stays a `dict` and `or None`
  is inlined at the call site (runtime and generated geometry unchanged).

### Changed
- **Second performance wave on the optimizer-refinement path**
  (`generator3d/{clean_geometry,embed,ml_optimizer,__init__}.py`, `cli.py`): the g-xTB optimize
  loop is single-threaded per call and parallelized across conformers (deterministic, ~8× on that
  stage); the MACE calculator is reused across conformers instead of rebuilt per structure;
  per-iteration FF construction in the `ff_clean` scan is now lazy; and `get_alternative_molecule`
  is memoized per generation. FF-path golden geometry stays byte-identical to v0.4.0.
- **Library `print()` routed to `logging`** (`core/translator.py`, `generation/metallogen_adapter.py`,
  `generator3d/*`, `oin/*`, `utils/*`): diagnostics go through module loggers, and the package
  attaches a `NullHandler` so importing `oinsmiles` never writes to the consumer's stderr on its own.
- HACF-toolchain framework tests (`test_autonomous_resolution`, `test_dynamic_orchestrator`,
  `test_provenance_integration` + 7 manual-plan `.md` files) moved from `tests/integration/`
  to `.agents/tests/` — they test `.agents/scripts/`, not the chemistry library. Deleted
  25 empty `.claude/commands/*.md` stubs and the stale `.claude/backup.commands/` snapshot.

## [0.4.0] - 2026-07-12

The MetalloGen **performance** wave (parallel worktree phases P0–P11), landed on top of the
v0.3.7 accuracy line. It makes 1D → OIN → 3D generation measurably faster **without changing
the generated chemistry** — at a fixed `--seed`, every generated XYZ is byte-identical to
v0.3.7. Golden-complex speedups (load-fair A/B, seed 42, FF-only): **cisplatin 12.0×,
ferrocene 5.6×, fac-Ir(ppy)₃ 6.9×, PdCl₂-BINAP 1.9×**. The wave confirmed two cost regimes:
small/medium complexes are bound by the PuLP/CBC bond-order solver (attacked by P2), large
ones by RDKit's ETKDG embedding (attacked by P3/P9). No OIN format change (still v3.7).

### Added
- **Deterministic generation + `--seed` (P1)** (`generation/engine.py`, `cli.py`, `generator3d/__init__.py`, `generation/metallogen_adapter.py`): generation is seeded end-to-end (default `seed=42`, per-attempt `seed + i*1009`), so the same OIN yields byte-identical XYZ across runs. `oin-smiles oin2xyz --seed N` and `OIN3DGenerator(seed=N)` sample a different but reproducible conformer. The dead `pool_size` parameter was dropped.
- **`tools/benchmark_generation.py` (P0)**: a serial, per-stage generation profiler for the four golden OINs and a CN-stratified sample, reporting median + IQR over N seeded runs. This is the ground-truth benchmark the wave's speed claims are measured against.
- **Curated metal–ligand bond-length table (`generator3d/bond_lengths.py`, P4/P8)**: metal-gated σ-donor distances (`ENABLED_METALS = {Ni, Pd, Pt, Zn, Cd, Hg, Ag}`) used as the FF-clean scan target so σ donors are pinned at physically realistic separations instead of the covalent-radius sum; unlisted pairs fall back to the covalent sum. Pd–P is exempted (`SHORT_PIN_EXEMPT_PAIRS`) to avoid over-tightening BINAP.

### Changed
- **Topology-keyed PuLP/CBC memoization (P2)** (`generator3d/utils/compute_chg_and_bo_pulp.py`, `generator3d/__init__.py`): the bond-order/charge solver dominates wall time for small and medium complexes. A per-generation memo keyed on fragment topology collapses redundant CBC subprocess solves within a single generation (cisplatin **432 → 21** solves). Load-fair marginal gain: cisplatin **9.7×**, ferrocene **5.3×**, BINAP **1.74×**. The cache is cleared per generation, so the speedup is a genuine per-invocation win, not a warm-cache artifact.
- **Removed the dead ETKDG rebuild-retry ladder (P3)** (`generator3d/embed.py`, `generator3d/__init__.py`): a redundant embed retry loop forced ~4× wasted ETKDG embeds per ligand on failure. Removing it cuts fac-Ir(ppy)₃'s ETKDG embed time from **38.6 s → 3.2 s** (~32.9 s of a 57 s baseline) with byte-identical output. An opt-in batched parallel embed (`embed_num_threads`) is available but off by default (it samples conformers differently, so it is validated by the accuracy gate rather than byte-identity).
- **Geometry-matcher candidate prefilter (P11)** (`utils/oin_aligner.py`): the O(n!) coordination-geometry matcher is gated behind a batched-numpy candidate prefilter that discards non-viable donor↔template permutations before the expensive `Rotation.align_vectors` step. Byte-identical winning permutation, RMSD, and rotation matrix vs the exhaustive sweep; `_map_to_template` **1.30 s → 0.07 s (17.8×)** on a CN-6 complex.
- **Deterministic embed rescue for octahedral loose-scale failures (P9)** (`generator3d/embed.py`): the primary ETKDG embed returns `-1` (not an exception) for some OCT loose-scale conformers; a deterministic `useRandomCoords=False` retry rescues each. fac-Ir(ppy)₃ ETKDG failures **16 → 1**; the success path is byte-identical.
- **Single-pass coordination perception (P5)** (`generation/metallogen_adapter.py`, `utils/oin_aligner.py`): `_select_by_geometry` classified and fit each candidate's coordination geometry twice; a `classify_and_fit` dedup does it once (Ir perception calls 20 → 10). Byte-identical.
- **Honest optimizer fallback + PASS-2 tier relabel (P6)** (`generator3d/ml_optimizer.py`, `tools/test_dataset_roundtrip.py`, `manuscript/figures/plot_timing.py`): when no `xtb`/g-xTB binary is present the optimizer now warns loudly and returns the FF geometry, and the round-trip harness relabels the FF-fallback tier (`g-xTB_N` → `FF_reroll_N`) so reports no longer silently claim g-xTB refinement that did not happen.

### Notes
- **P10 (in-process CBC/HiGHS solver) shipped no code** — an in-process solver was neither byte-identical to nor faster than the existing CBC subprocess call; recorded as a negative result.
- **P4/P8 are net-neutral on accuracy.** The curated bond lengths did not move coordination-sphere mean RMSD beyond noise (median +0.016 Å) and every distance stays well under the `rmsd ≥ 1.0` fail gate, so no molecule flips pass/fail; P8 erases the FF-clean speed regression P4 introduced on BINAP. Kept as harmless; a v0.4.1 follow-up will decide whether the pair earns its complexity.
- **Accuracy is unchanged from v0.3.7.** Because generation is byte-identical at a fixed seed (P2/P3/P5/P11) or only turns a failed embed into a success (P9), the ≈89% tmCAT/tmPHOTO round-trip pass rate carries over; this wave adds speed, not coverage.

## [0.3.7] - 2026-07-12

The tmCAT/tmPHOTO round-trip **residual** fix wave (parallel worktree sessions R1–R5),
building on the v0.3.6 fidelity wave. Each session owned a disjoint set of files and
closed one class of residual round-trip failure. Dataset FF-only pass rate: **≈89%**
(`--quick` sampled estimate, 95% CI 87.7–90.3%; up from the 85.4% pre-wave baseline).
This release also removes the legacy SCINE Molassembler 3D-generation backend and its
`scine-molassembler` dependency (see **Removed** — MetalloGen has been the default
engine since v0.3.3).

### Removed
- **Legacy SCINE Molassembler 3D-generation backend** (`generation/molassembler_adapter.py`, ~3,500 lines): MetalloGen has been the default engine since v0.3.3, and the Molassembler `engine="legacy"` path is now removed along with the `scine-molassembler>=2.0.0` dependency. `OIN3DGenerator` keeps its signature; `engine` must be `"metallogen"` (the default) and any other value raises `ValueError`. The CLI `--engine` flag and the `MolassemblerTimeoutError` exception are gone. This drops the project's only heavy C++ dependency, which was previously imported eagerly at module load even on the default path.
- **Deferred "direct parser" dead code** (`generation/engine.py`, `generation/oin_parser.py`): the never-integrated `parse_oin_direct` pipeline and its masm-coupled helpers (`construct_molassembler_mol`, `convert_bond_type`, `_extract_oin_constraints`, `tokenize_unsanitized_smiles`, `_dg_worker`, `_translate_eta_vertex_to_atoms`, and the unused `normalize_template`) were removed. The live `OINParser.parse()` path (inline + sidecar) is unchanged.
- **Orphaned `core.chirality._attach_dummy_metal` helper**: only ever called from the removed Molassembler stitch path, now dead and deleted.

### Changed
- **Pinned RDKit to an exact version** (`pyproject.toml`): `rdkit>=2025.9.3` → `rdkit==2025.9.3`, so every install resolves to a single, tested RDKit release rather than floating up to whatever is newest (RDKit's version affects ETKDG embedding determinism).

### Internal
- The `GeneratedStructure` result type moved to a lightweight `generation/structure.py` module (previously defined inside the deleted `molassembler_adapter.py`). Its legacy-only `haptic_face_decisions` field was dropped.
- Removed the tests dedicated to the legacy engine / direct parser, and reduced two mixed test modules to their engine-agnostic (encoder / parser) coverage.
- Reconciled the `spec/compiled/architecture.yml` hypergraph: deleted the 6 nodes for the removed engine + direct parser and their edges, relocated the `GeneratedStructure` node to `structure.py`, and moved the affected surviving nodes back to `clean`.

### Fixed
- **Metalloporphyrin meso E/Z spuriously stamped on the fast re-encode (R1)** (`generation/metallogen_adapter.py`): `build_contract_mol` feeds `get_oin_string` (the generator's fast re-encode), where `AssignStereochemistryFrom3D` stamps an E/Z marker on every localized ring double bond — including a metalloporphyrin's meso C=C/C=N bridges, which are ring-locked through the metal and carry no free E/Z. The forward XYZ→OIN encode strips those (v0.3.6.1 `_clear_chelate_locked_bond_stereo`), but the fast path skipped it, so the generated OIN carried slashes the input never had. The same clear now runs on the contract mol. `macrocycle_perception` **25/46 → 34/46** round-trip (+9 E/Z-only rows, zero regressions); the 2 remaining E/Z-only rows (VAVRAN, XIZXAG) carry E/Z on pendant, non-ring-locked bonds and are correctly left alone. Guarded by `tests/unit/test_contract_mol_chelate_ez.py`.
- **`_apply_double_bond_stereo` forced an over-valent double bond, so some molecules generated no conformer at all (R2)** (`generator3d/embed.py`): a carried C=C/C=N stereo bond was restored to `DOUBLE` unconditionally to enforce its E/Z. When the dummy-metal PuLP re-perception had relocated the double bond, that promotion made an endpoint over-valent (`FIXYER_comp_0`: a 5-valent carbon), and every downstream `SanitizeMol`/`MolToSmiles` — including a debug `print` inside the embed loop — raised, so **generation returned nothing**. The promotion is now guarded by a valence check (`_promotion_keeps_valence`, a throwaway `SANITIZE_PROPERTIES` copy) and only applied when it keeps the molecule valence-valid, degrading to the documented "leave it and skip the constraint" behavior otherwise. Recovers `FIXYER`, `EDOFUB`, `EDOGEM`, `ZIHGEE` (all carry a `/C=C/`/`/C=N/` whose forced promotion over-valenced a carbon) to **clean round-trips** (rmsd 0.10–0.58); `PILWUC` (same mechanism) now generates but reclassifies to `string_mismatch`. Attribution A/B'd against pristine `embed.py`: without the guard `FIXYER` exhausts the attempt budget and generates no conformer (273 s, `MetalloGen failed to generate any conformers`); with it, a clean round-trip in 132 s. Guarded by `tests/unit/test_embed_budget.py`; regression floor `tests/unit/test_generator_double_bond_stereo.py` (a valence-valid pendant alkene is still promoted and its E/Z enforced).
- **Chelate-locked E/Z lost to the comparator's `RAW:` fallback + nitride/ammine notation collision (R3)** (`oin/compare.py`, `oin/inline.py`): (Part C) a slot-stripped eta-Cp/Cp* or bare-`n` azole/pyridyl donor ring raises `KekulizeException` on a full sanitize, so the chelate-lock E/Z clearing never reached it and the raw `/`…`\` slashes survived as a spurious mismatch; `_parse_fragment` now retries with a partial sanitize that skips kekulization. Closes 23 dataset rows with 0 regressions (6 `EZ_bond_stereo`, 9 `macrocycle_perception`, 4 `garbled_aromatic`, 3 `string_mismatch_other`, 1 carborane). The 4 exocyclic-azine E/Z residuals (PEDPOG/PEDPUM/PEDQEX/VUQDUI) are left distinct on purpose — topologically identical to the genuine flips AFECIZ/RIQFON, so masking them would false-pass real diastereomers. (Part B) a heavy==0 terminal nitride now serializes as a bracketed `[N]{n}`, distinguishable from an ammine `N{n}`; WAYHOW/AFIZEV/FUVNER/OPUYAB round-trip 0-H while the ammines AFAVIO/RIZVAY/OQIHUT/XILBIF still round-trip as NH₃. (Part A) the nitro `[N+](=O)[O-]` vs `N(=O)=O` case was verified already resolved on current main — no comparator change. Guarded by `tests/unit/test_chelate_ez_comparator.py` (extended) and `tests/unit/test_nitride_ammine_notation.py`.
- **Non-binding 0-H chalcogen donors gained a phantom hydrogen (R4)** (`utils/oin_aligner.py`): `OINSanitizer.generate_robust_smiles` locked H counts only for metal-binding atoms, so a non-binding heteroatom with a valence deficit (croconate/oxo ring O, nitrito –O) serialized bare (`c(O)`, `ON=O`) and the adapter's `MolFromSmiles` re-added the implicit H — a phantom the input never had (COLWIK croconate 55→58, ACOXEX nitrito 75→77). Each non-binding O/S deficit is now charged by its deficit (a valence-1 O → `[O-]`) — 0-H, closed-shell, and embeddable (a neutral 0-H O would be a radical UFF cannot type). 21/21 targeted O-deficit rows flip failed→success; aqua/hydroxo/carbonyl/ether O and real O–H are untouched, and the S1 ammines stay NH₃. N is excluded (it overlaps R3's nitride/ammine notation). Guarded by `tests/unit/test_non_binding_donor_hydrogens.py`.
- **sp3 / heteroatom atom-stereo (@/@@) disagreement between encode and re-encode (R5)** (`core/chirality.py`, `generation/metallogen_adapter.py`): the forward encode and the generator's contract-mol re-encode disagreed on `@`/`@@`, taking `atom_stereo` from **0/25 to 18/25** full round-trips (the remaining 7 fail on a different class; no `@`-disagreement residual remains). Four mechanisms: (1) stop inventing sp3 stereo the OIN left unspecified — `AssignStereochemistryFrom3D` stamped a tag on every chiral-looking embed centre (KAPCEM 0→4); (2) drop the spurious `CHI_OCTAHEDRAL` tag on achiral –SF₅ sulfur (MEDHUB); (3) orient a metal/eta-adjacent specified sp3 centre against the **metal-free** fragment's `rdCIPLabeler` label (`_OIN_CIPCode_SP3` stamp + a `recover()` verify-and-flip, mirroring the Zone-A P lone-pair branch), since a metal-adjacent centre's CIP label flips between the metal-present contract mol and the emitted metal-free fragment (AHEBEV); (4) take that label aromatic-preserving (`SANITIZE_ALL ^ SANITIZE_KEKULIZE`) on a **fresh re-parse** of the template SMILES, since a `RemoveHs(sanitize=False)` corrupts the aromatic state of a fused indenyl/fluorenyl ring and flips the label (BABWAD, KAGXUM, NOSGAD; the same re-parse fixes the aromatic-arm P-donor GUXPIA). Guarded by `tests/unit/test_heteroatom_atom_chirality.py`.

### Added
- **Generation-internal wall-clock budget for the embed loop (R2)** (`generator3d/__init__.py`, `generation/metallogen_adapter.py`): the FF-only attempt loop had no time bound — `timeout` was consumed only by the ASE optimizer — so a molecule whose embed never validated ran the full `max_attempts` (250) budget before returning nothing (`ZIHGEE_comp_0`: ~1696 s). `generate_3d_structures` takes a new `embed_time_budget` (wired to the existing per-molecule `timeout`, 300 s full / 60 s quick) and stops the loop at the deadline — checked *between* attempts, so a molecule that does embed is never interrupted, and the pool built so far is returned; an empty pool becomes the same `[]` as before, only fast. Turns a pathological non-terminating case into a fast, honest failure without changing the outcome of a molecule that embeds. `None` (the default) preserves the prior unbounded behavior for direct callers. Guarded by `tests/unit/test_embed_budget.py`.

### Triage / documentation
- **`no_conformers` class (R2):** all 36 rows the post-S6 registry filed under `no_conformers` were re-run serially (the registry over-counts the class — a concurrent sweep fabricates the error, so several rows are contention flakes that generate on a serial retry). After the two fixes above, **27 of 36 now generate a conformer**, and every row has a verdict: **16 round-trip cleanly** (net-new passes — `ABERUW CITGAO DURSOZ FIXYER EDOFUB EDOGEM ZIHGEE AGOGEJ RUBPAH ZUJNAT YEPXID IJAXIB IMELIW VEZYOQ TEZTAV YOMDAH`, rmsd 0.10–0.75); **11 generate but reclassify** into smaller-defect classes owned by other sessions (6 `atom_count`, 3 `string_mismatch`, 2 `high_rmsd`); and **9 are genuine and documented** in `docs/KNOWN_LIMITATIONS.md` (neutral L-donor over-valence `FUVNER`/`GEZKAZ`/`VIBRIK`; exotic bond orders the two-centre perception can't build `DAHXOB`/`MEDDUV`/`IREPAX`/`DOFCAE`/`HURGOS`; and one geometry-realization gap `BOBJIM`). The generation seed is fixed (42), so these verdicts are deterministic.
- **Residual-tail triage (R4)** (`docs/agentic-notes/v0.3.7/roundtrip_residual_triage_R4.md`, `tools/triage_overrides.json`): the remaining unowned buckets were RDKit-diffed and routed — `string_mismatch_other` (28) is entirely stereo-only under per-fragment canonicalization (13 → `EZ_bond_stereo`, 13 → `atom_stereo`, 2 R3-resolved); `atom_stereo` (25) → the R5 session; singletons documented (DEKQAN `geometry_NON` → unperceived Y–O bonds; atom-count-decrease rows → generator atom loss, not encoder phantom-H).

## [0.3.6.1] - 2026-07-10

### Fixed
- **Chelate-locked C=N/C=C E/Z regression in the round-trip comparator** (`oin/compare.py`): v0.3.6's S6 encoder drops E/Z on a double bond a metal ring holds rigid, but only where the *input* structure's donor bonds make that ring perceivable — a generated structure whose donor is bonded differently keeps the marker, so `smiles_1` (no slash) and `smiles_2` (slash) described the same chelate yet failed the string gate (a re-run flagged 62 salicylaldimine/hydrazone rows, `EZ_bond_stereo` 4→64). `canonical_roundtrip_key` now reconstructs the chelate rings (a dummy metal bonded to every slot atom) and clears E/Z on the double bonds those rings lock, on both sides, so they compare equal — while a pendant, freely-rotatable alkene/imine in no metal ring keeps its E/Z and still distinguishes a real diastereomer. Fixes 52 of the 62 (the other 10 are `RAW:`-fallback parse failures deferred to the v0.3.7 R3 session); genuine flips (AFECIZ, RIQFON) still fail correctly. Guarded by `tests/unit/test_chelate_ez_comparator.py`.

## [0.3.6] - 2026-07-10

The tmCAT/tmPHOTO round-trip fidelity wave (parallel worktree sessions S1–S6). Each
session owned a disjoint set of files and fixed one defect class found in the dataset
baseline; the round-trip harness compares by the v0.3.5 `canonical_roundtrip_key`.

### Fixed
- **Bare anionic terminal donors re-protonated on generation (S1)** (`generation/metallogen_adapter.py`): `convert_parsed_to_msmiles` stripped a bound donor's implicit H only when it had ≥2 heavy neighbours, so a *terminal* anionic N donor (silylamide, anilide, azide, phosphinimide `P=N`) — which has exactly one heavy neighbour — was rebuilt as NH₂/NH₃ and the generated atom count no longer matched (e.g. XADYAC, UDIVUY, FENMIX, QIPKOS, XOSSEE). Gate widened to `heavy >= 1`. Separately, an outer-sphere/uncoordinated fragment with no binding-slot vectors raised a raw `IndexError` at `frag_vectors[0]`; it now raises a typed `UncoordinatedFragmentError` (an honest "MetalloGen cannot place a free counter-ion", e.g. BEYHEU, CUBDOT, NECCIH, NASZOY). Guarded by `tests/unit/test_bare_donor_hydrogens.py`.
- **η²-diene / COD double bond localized onto the ring backbone (S2)** (`generation/metallogen_adapter.py`): `build_contract_mol` transferred bond orders from the OIN template to the generated fragment by an unconstrained automorphism search (both graphs are equal-size, all-single, heavy-atom-only), so a valid map could land the `C=C` off the metal-bound carbons (e.g. 1,5-COD `[CH2]=[CH2]` on GASBIN, PENGAT; drove ~66% of Rh failures). `_flatten_template` now colours every atom with its OIN coordination `_oinSlot`, a bitmask-DP assigns generated donors to slots, and `_transfer_score` breaks intra-slot ties against the embedded 3D geometry — the legacy map is kept only when slot-valid within tolerance (a repair, not a re-pick, so charges and `_CIPCode` ride along). Cleared the `eta_diene_localization` bucket (78→0 of the sampled rows). Guarded by `tests/unit/test_contract_mol_diene_transfer.py`; regression floor `test_contract_mol_allyl_transfer.py`.
- **Aromatic / charge perception garbled in the XYZ→OIN encoder (S3)** (`utils/perception_tmc.py`, new `utils/aromaticity.py`, `generator3d/process.py`): `get_oin_string` rebuilt each fragment copying bond *type* but with `IsAromatic=False`, so an aromatic ring re-serialized as an unparseable mixed single/double `c=c` (`RAW:` token, key never matched); a blanket `SetFormalCharge(0)` turned `[N+](=O)[O-]` into unparseable `N(O)=O`; `lig_checks` spent a `ResonanceMolSupplier` cursor with `len()` then iterated it → `None`; and `=P` ylide / `[CH]` radical ligands crashed with `KeyError(BondType.AROMATIC)`. Now normalizes aromatic-flagged bonds to `BondType.AROMATIC`, restores formal charges in pairs, indexes the resonance supplier, and runs a `kekulize_safe_sanitize` + `_rescue_unusable_perception` charge sweep so all 17 kekulize cases at least encode. Guarded by `tests/unit/test_aromatic_reencode.py`, `test_encoder_perception.py`. (Porphyrinoid macrocycles and the ylide/radical *generation* remain open — see KNOWN_LIMITATIONS.)
- **Eta-ring winding flipped on symmetric rings (S4)** (`utils/oin_aligner.py`, `utils/perception_tmc.py`): a fully symmetric η-ring (Cp*, mesitylene, borate phenyl) emitted `{n>}` vs `{n<}` depending on the embedding, because winding is meaningful only when *no* fragment automorphism reverses the eta group's cyclic order (turning a ring over is a proper rotation that leaves the metal fixed). The encoder now computes this via a labelled canonical-SMILES symmetry graph and emits a fixed marker for orientation-free rings; substituted rings that carry a real rac/meso distinction (TiCat3/4/5) are untouched. Also fixed the set-valued eta-slot placement on tetraphenylborate (SOJMIQ) and a latent `None` fallback on 4-coordinate neutral boron. Guarded by `tests/unit/test_automorphic_ring_winding.py`.
- **Coordination-sphere RMSD reported false 996/999 sentinels (S5)** (`tests/integration/rmsd_utils.py`): 62 rows reported `High RMSD: 996/999` despite correct chemistry — the metric never ran. `999` rows were all Y or Sc, the two transition metals missing from a hand-copied `METAL_ATOMIC_NUMBERS` set (now imports `core/constants.py::TRANSITION_METALS_NUM`); `996`/`997` rows chose the input coordination sphere by a covalent-radius cutoff that erred both ways (a real 2.57 Å apical Pd–N excluded; a 2.19 Å non-donor C admitted), now replaced by a per-element nearest-k match with a ceiling safety net. The metric returns an honest `(rmsd, None)` / `(None, reason)`; the harness emits `RMSD mapping failed: <reason>` instead of a magic number. The `--mol-timeout` watchdog (a `signal.alarm` that could not interrupt native code and only wrapped generation) is rewritten to run each pass in a `spawn` subprocess SIGKILLed on expiry, with the encode moved into the child — so an encoder hang (UGUHAH, >150 s inside `XYZToSMILES.convert()`) no longer wedges a run. Guarded by `tests/unit/test_rmsd_mapping.py`, `test_roundtrip_watchdog.py`.
- **Stereo lost through the MetalloGen embed (S6)** (`generator3d/embed.py`, `generator3d/ligand.py`, `generator3d/__init__.py`, `core/translator.py`): the embed drew its RDKit seed from an *unseeded* `random`, so generation was non-deterministic (and the "task #8" E/Z test was flaky for that reason); the seed is now threaded so three runs are byte-identical. sp3 atom chirality is captured from the parsed ligand and re-asserted after sanitize (`enforceChirality` + a global remap), so a `[C@H]`/`[C@@H]` donor backbone returns the correct enantiomer. And chelate-locked C=C/C=N E/Z is preserved: RDKit ring perception ignores `DATIVE` bonds, so a ring-locked double bond looked acyclic and was dropped — `_clear_chelate_locked_bond_stereo` now perceives rings on a scratch copy with DATIVE upgraded to SINGLE (fixes VOacac2; `verify_xyz_to_oin.py` 26/27→27/27). Guarded by `tests/unit/test_generator_atom_chirality.py`, `test_chelate_locked_ez.py`, and a de-flaked `test_generator_double_bond_stereo.py`.

### Known limitations (documented, not fixed this wave — see `docs/KNOWN_LIMITATIONS.md`)
- **Porphyrinoid macrocycles** re-encode localized because the OIN string carries no formal charge, so MetalloGen builds the pyrrolide N neutral and the contract mol will not sanitize (a donor-charge-layer problem, not encoder aromatic perception).
- **Ylide (`=P`) / radical (`[CH]`) ligands** now encode without crashing but do not round-trip — the 3D generator cannot reproduce those bond orders.
- **Borane / carborane clusters** remain unsupported (multi-centre bonding outside the two-centre perception model).

## [0.3.5] - 2026-07-08

### Added
- **Structure-level canonical round-trip comparator** (`oin/compare.py`): `canonical_roundtrip_key(oin)` = (metal + geometry, sorted multiset of RDKit-canonical ligand-fragment SMILES, winding multiset). It collapses notation drift that made chemically-identical structures compare unequal (implicit-vs-explicit donor H, NHC carbene bare-C vs `[CH2]`, which symmetric carboxylate O carries the slot, fragment ordering) while still distinguishing genuine connectivity/metal/geometry/winding differences. Lightweight (rdkit + `OINInlineHandler` only, no 3D stack); `normalize_oin_for_comparison`/`winding_canonical_key` moved here with `verify_roundtrip.py` re-exporting for back-compat, and the dataset harness now compares by this key.
- **Canonical symmetric-donor binding slot** (`utils/oin_aligner.py`, `utils/perception_tmc.py`): a monodentate ligand that binds through one of two resonance-equivalent atoms (e.g. a carboxylate's two oxygens, which differ as `=O`/`-O` in any single Kekulé structure) now always carries the `{slot}` marker on a canonically-chosen atom. Structures that differ only in which atom 3D bond perception happened to pick now encode identically, fixing a spurious round-trip "String mismatch". Guarded by `tests/unit/test_canonical_donor_binding.py`.
- **`tools/recalculate_oin_smiles.py`**: A utility to recalculate OIN SMILES strings from both the input XYZ and the generated XYZ structures for previously processed datasets. Updates the `summary_roundtrip.json` and `individual_reports` statuses if a codebase change causes a previously failed mismatch to now perfectly round-trip.
- **`tools/rebuild_summary.py`**: Rebuilds `summary_roundtrip.json` from the per-molecule `individual_reports` already on disk.
- **Dataset Roundtrip Tools enhancements**: Added `--quick`, `--continue`, `--rerun-failed`, `--random`, and `--mol-timeout` options to `tools/test_dataset_roundtrip.py` to allow robust, resumed, and time-bounded background processing of large datasets without hanging on pathological UFF geometries.

### Fixed
- **Carbene / dative-amine H miscount in m-SMILES** (`generation/oin_parser.py`): `convert_parsed_to_msmiles` used an ad-hoc heuristic to strip a binding atom's implicit H — it built an NHC carbene carbon as `CH2` (+2 H) and stripped a dative secondary amine's N–H (misclassifying any N with ≥2 heavy neighbours as anionic amido), producing wrong generated atom counts (ACAWOR 61→65, ABESAD 91→90). Replaced with a principled rule: an explicit bracket-H count (`[NH]`/`[OH2]`/`[CH2]`) is authoritative and kept (neutral L-type donor); only a *bare* binding atom is reinterpreted (bare chalcogen → 0 H alkoxide/thiolate/oxo, bare N ≥2 heavy → amido, bare non-aromatic C ≥2 heavy → NHC carbene/carbanion). Guarded by `tests/unit/test_msmiles_donor_hydrogens.py`.
- **Backbone P/S/Si stereocentre lost on round-trip** (`generation/metallogen_adapter.py`): a generated backbone P came back achiral (e.g. ABOPOY's iminophosphorane `N=[P@]`) because the round-trip path skips `CIPAssigner`, so `recover()` treated the emitted `[P@]` as a stray and cleared it. `build_contract_mol` now stamps `_OIN_CIPCode` on a backbone (non-donor) P so `recover()` keeps and orients it, and the perceive-then-flip carry set is extended to Si/S so their handedness is matched to the template. (Zone-A chiral P donors are handled separately below.) Guarded by `tests/unit/test_backbone_heteroatom_stereo.py`.
- **C=C cis/trans (E/Z) stereo dropped end-to-end** (`utils/perception_tmc.py` encoder, `generator3d/ligand.py`/`chem.py`/`embed.py` generator): `get_oin_string` rebuilt each fragment copying bond *type* but not the stereo reference atoms, so both E and Z re-encoded as a bare `C=C`; and the stereo-blind MetalloGen embed re-perceived bond orders with a cis bias, so an E input could regenerate as Z. The encoder now carries double-bond stereo across the rebuild and materializes `/`\`\` directions before the canonical SMILES, and the generator recovers E/Z from the parsed bond directions and enforces it on the embed (skipping metal-coordinating/donor-adjacent double bonds, which chelation already locks). Because the MetalloGen embed uses a random seed and newer RDKit honors the distance-geometry stereo constraint less reliably, `generate_3d_structures` now also verifies each conformer's alkene dihedral and rejects any that embedded the wrong (or a distorted, ambiguous) side, so the returned structure reproduces the requested E/Z deterministically across RDKit 2025.09.3 and 2026.03.3. Guarded by `tests/unit/test_double_bond_stereo_encoding.py` and `tests/unit/test_generator_double_bond_stereo.py`.
- **CN-8 (square-antiprismatic) crash + quinoid-ligand parsing** (`utils/oin_aligner.py`, `generation/oin_parser.py`, `generator3d/process.py`): 8-coordinate complexes (e.g. AFEPIM, Hf + 4 bidentate salicylaldimine) had no encoder template, so the aligner emitted `g:NON` and generation crashed with "Geometry code 'NON' not supported" — added a square-antiprismatic `SQA` template (encoder `TEMPLATE_SPECS` + `n==8` match, parser `TEMPLATES`, and the `OIN_TO_METALLOGEN_GEO` mapping). Separately, amidinate / 2-iminopyridine ligands (a quinoid aromatic ring with an exocyclic C=N and no valid Kekulé) crashed conformer generation with `KeyError(BondType.AROMATIC)`; a targeted `_dearomatize_stuck_rings` de-aromatizes only the offending ring. Guarded by `tests/unit/test_cn8_geometry.py` and `tests/unit/test_quinoid_ligand_parse.py`.
- **η³-allyl double-bond loss on round-trip** (`generation/metallogen_adapter.py`): `_flatten_template` built the connectivity-only query for `build_contract_mol`'s per-fragment substructure match but cleared only aromaticity and charge — not radical electrons. A ligand atom that binds the (stripped) metal is under-valent, so its template atom carried a radical, and `GetSubstructMatch` treats radical count as a match constraint — so bond-order/aromaticity transfer silently failed and the ligand was emitted all-single and de-aromatized (the η³-allyl double-bond loss, e.g. ABAZEK). Now clears radicals and normalizes H valence so the match succeeds and the allyl `=` is preserved by transfer; five dataset allyl cases (ABAZEK/ABETIK/ABETOQ/ACALOI/AGOVOK) now key-match. Guarded by `tests/unit/test_contract_mol_allyl_transfer.py`.
- **Zone-A chiral P donor stereo lost on round-trip** (`generation/metallogen_adapter.py`): a phosphorus binding the metal directly with a stereogenic lone pair encoded as `[P@]{0}` on XYZ→OIN but re-encoded achiral `P{0}` after OIN→3D→OIN (e.g. ACUWUT). `build_contract_mol` never populated the lone-pair path `recover()` needs for the donor, and the dative metal→P bond makes `AssignStereochemistryFrom3D` return `CHI_UNSPECIFIED`. The generated donor is now primed with `_OIN_CIPCode_LP` and a seeded chiral tag from `rdCIPLabeler` on the metal-free template (not the legacy `Chem.AssignStereochemistry` label, which disagrees for 3-coordinate P and would round-trip the wrong enantiomer), so `recover()`'s lone-pair verify-and-flip branch re-asserts the encoded handedness. Guarded by `tests/unit/test_zone_a_p_donor_stereo.py`.
- **`--quick` roundtrip crash from unsupported `max_attempts`** (`generator3d/__init__.py`): `--quick` mode passes `ff_params={"max_attempts": 10}`, which `generate_3d_structures` forwarded verbatim to `TMCOptimizer(**ff_params)` — but `TMCOptimizer.__init__` has no such parameter, so every molecule raised `TypeError: unexpected keyword argument 'max_attempts'` and crashed before any geometry was produced. `max_attempts` is now filtered out before constructing the optimizer while still capping the embedding retry loop.
- **UFF loop hang during 3D generation**: Fixed a bug in `generate_3d_structures` where an unrecognized UFF atom type (which immediately fails FF cleaning) caused the generator to blindly brute-force 250 random embeddings before giving up. Added support for a `max_attempts` override in `ff_params` (10 under `--quick`) that caps the embedding loop.

### Changed
- **g-xTB optimizer renamed `xtb` → `g-xtb`** across the user-facing surface: the `oin-smiles oin2xyz --optimizer` default, the `ASEOptimizer` method (which still accepts both spellings), and the dataset roundtrip harness tiers. The dataset harness now also short-circuits hard generation/verification failures and timeouts instead of escalating them to the slow g-xTB pass.

## [0.3.4] - 2026-07-07

### Fixed
- **Default g-xTB optimizer crash.** `ASEOptimizer.optimize()` unconditionally called `atoms.get_potential_energy()` after refinement, but the g-xTB path (a subprocess, not an ASE calculator) never attaches one — so every *successful* g-xTB optimization raised `Atoms object has no calculator`. The energy is now read from the ASE calculator only on the MACE path; the g-xTB path uses the energy parsed from `xtb` stdout. Guarded by a new regression test. (`generator3d/ml_optimizer.py`)
- **Non-functional charge/bond-order code paths.** `process.get_chg_and_bo` / `get_bo_matrix_from_adj_matrix` / `get_chg_list_from_bo_matrix` referenced a `frag`/`compute_scipy` module not vendored in this subset (would `NameError` if reached); they now route through the active PuLP solver (`compute_chg_and_bo_pulp`). Also fixed a missing `self` parameter on `Molecule.get_screening_result`, a `print_x` typo in the `Ligand` debug print, and a duplicate `get_electron_count` definition.

### Changed
- **MACE / PyTorch are now an opt-in extra.** `mace-torch` and the pinned CUDA-11.8 `torch` moved from hard dependencies to `pip install/uv sync --extra mace`. The default install is lightweight (FF + g-xTB, no `torch`), matching the documented "fast FF-only path". MACE optimizers (`mace-omol-*`) require the extra and still fail loudly if unavailable.
- Removed dead, unreachable, partially-ported code from the vendored MetalloGen engine (fragment-based charge/BO heuristics, `get_kth_neighbor_atom_list`, `scatter_molecules`, and the `frag`-based `Detect_EZ`/`Detect_RS`/`Detect_stereocenter` stereo path — OIN-SMILES handles stereo upstream).

### Added
- **Configurable g-xTB timeout with FF fallback.** A `timeout` (default 300s) threads through `OIN3DGenerator` -> `MetalloGenAdapter` -> `generate_3d_structures` -> `ASEOptimizer`; when `xtb` exceeds it, `subprocess.TimeoutExpired` is caught and generation falls back to the FF geometry instead of hanging (the old fixed 60s cap was too short for large complexes).
- **`tools/install_mace_weights.sh`**: idempotent downloader for the public `MACE-omol-0-extra-large-1024` checkpoint (ACEsuit GitHub release) that also registers `MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH` in `.env`. `models/mace/README.md` reworked accordingly (the extra-large model is freely downloadable; OMol25 remains Hugging-Face-gated).
- **`tests/unit/test_generator3d_units.py`** (32 tests): unit coverage for the vendored `generator3d` engine — `chem.Atom`/`Molecule`, `process` helpers, the frag→PuLP reroute, package helpers, and `ASEOptimizer` (incl. the g-xTB regression test).
- **Lint tooling & CI**: `ruff` added as a dev dependency (`[dependency-groups]`); the entire repo is now `ruff check`/`ruff format` clean (was 1817 findings). Added `.github/workflows/ci.yml` running lint + the unit and encoder suites on push/PR.
- **Packaging metadata**: `pyproject.toml` gains a proper description, `authors`, `license`, `keywords`, `classifiers`, and `[project.urls]`.
- **Docs**: `docs/OPTIMIZERS.md` (FF vs g-xTB vs MACE selection/install) and a `docs/README.md` index; README installation steps corrected (numbering + light-vs-MACE install).

### Removed
- Stray repository clutter: `os` (empty), `main.py` (hello-world stub), `REPORT.md`, `OpenSourceTMCBuilderReport.md`, a stray `verify_xyz_to_oin.py.fragment`, and four unrelated HACF-framework docs under `docs/`. The root HACF installer was renamed `install.sh` → `HACF-install.sh` to distinguish it from project installation (`uv sync`).

## [0.3.3] - 2026-07-06

### Changed
- **Default generation engine → MetalloGen; default optimizer → g-xTB.** `OIN3DGenerator.__init__` now defaults to `engine="metallogen"` and `optimizer="xtb"` (was `engine="legacy"`, `optimizer=None`), and `cli.py oin2xyz` rides that default. The MetalloGen backend is the better-validated generator (full MACE round-trip 25/25; eta-winding rac/meso, BDNN square-plane, and TiCat eta TET/TPY all fixed). g-xTB is fast and, unlike MACE, **degrades gracefully to FF** when the `xtb` binary is not on `PATH` (so the default never hard-fails). MACE (`mace-omol-0-extra-large-1024` / `mace-omol25`, higher accuracy, requires `mace-torch` + weights and fails loudly if absent) and the legacy SCINE Molassembler backend (via `engine="legacy"`, still the reference for Zone-A P stereo enforcement) remain opt-in.
- **`verify_roundtrip.py --optimizer` default → `xtb`** (was MACE) for fast iteration; pass `mace-omol-0-extra-large-1024` for the accurate sign-off.

### Added
- **g-xTB optimizer** (`generator3d/ml_optimizer.py`): a subprocess wrapper around the `xtb` binary (`xtb <struc> --gxtb --opt`) selected via `optimizer="xtb"`. Fast semi-empirical refinement of the FF pool; warns and falls back to FF if the binary is missing. Helper `tools/install_gxtb.sh` installs the binary; `tests/integration/run_optimization_grid.py` benchmarks FF/g-xTB/MACE across the fixtures.
- **`oin-smiles oin2xyz --engine {metallogen,legacy}` and `--optimizer`** (default `xtb`; accepts `ff`/`none`/`mace-omol-*`), giving a fast, no-torch, or legacy path from the CLI.

### Fixed
- **Legacy-specific real-generation unit tests pinned to `engine="legacy"`** (`test_winding_inertness.py`, `test_zone_a_p_genenforce.py`, `test_stereo_roundtrip_diagnostics.py`) so the default flip keeps the fast unit suite deterministic and heavy-optimizer-free while those Molassembler-only behaviors stay under test.

## [0.3.2] - 2026-07-06

### Added
- **Eta-ligand winding (rac/meso) round-trip fidelity**: the encoder now emits a winding marker per haptic slot — not just the first ring — measured against each ring's *actual* metal→centroid axis (`oin_aligner.py` `_permute_and_serialize` / `_determine_winding`, was the idealized template slot axis that flipped the 2nd ring under a distorted ansa bite). The MetalloGen generator honors it via winding-**multiset** conformer selection over a widened eta pool (`metallogen_adapter.py` `_eta_winding_multiset` / `_reencode_oin_fast` / `ETA_SELECT_POOL`), fixing the TiCat3/TiCat4 rac↔meso diastereomer swap (both now round-trip to the correct isomer). Generalized to N eta ligands and variable hapticity (η³ definite / η² degenerate); `verify_roundtrip.py` compares via a winding-canonical key (winding-stripped string + sorted multiset). Adds `TiCat5`/`TiCat6` fixtures and `test_eta_winding_generalization.py` (8 tests).
- **`oin_aligner.py`**: Added `classify_coordination_geometry()` (best-matching OIN geo code for a set of metal-centred donor vectors) and `coordination_geometry_fit()` (RMSD of the best donor-to-template assignment — the fit quality the classifier itself discards). Both wrap the existing discrete-geometry matcher.
- **`MetalloGenAdapter`**: Added `_select_by_geometry()` geometry-code-aware conformer selection. From the energy-ranked pool it keeps only conformers whose coordination sphere classifies as the requested geometry, then returns the tightest template fit (energy breaks ties). Haptic/η donors are gated out (donor count ≠ coordination number), making selection a deliberate no-op there and strictly non-regressive versus lowest-energy.
- **`generate_3d_structures`**: Added conformer deduplication over the FF pool via new `uff_pool_size`, `rmsd_threshold`, and `energy_threshold` parameters, plus a `calculate_heavy_atom_rmsd()` helper. `MetalloGenAdapter` surfaces these through `ff_params`.
- **`rmsd_utils.py`**: Added `_compute_robust_rmsd()` (anchor-pair candidate rotations → Hungarian assignment → Kabsch refine → ICP polish, floored by the greedy estimate) for the >5-atoms-per-element branch, fixing bent ansa-metallocene mis-pairing.
- **`test_geometry_selection.py`**: Added 18 unit tests covering the classifier, template-fit ranking, coordination perception, the haptic gate, and lowest-energy fallbacks.
- **`run_verification.sh`**: Added a `--limit N` pass-through to `verify_xyz_to_oin.py` and `verify_roundtrip.py`.
- **`tools/test_dataset_roundtrip.py`, `tools/test_uff_pool_size.py`**: Added dataset round-trip and UFF-pool-size sweep scripts.

### Fixed
- **`MetalloGenAdapter`**: Fixed the stochastic PdCl2-RR-BDNN failure where the generated Pd distorted from square-planar (`SPL`) toward trigonal-pyramidal, giving RMSD ~1.4 and a geo-code mismatch. Geometry-fit-ranked selection now returns the cleanest square-plane from the pool (BDNN: 5/5 round-trip PASS, RMSD 0.10–0.20; was intermittent). The prior `--ensemble-size` lever was a no-op under an optimizer — the pool is fixed at `pool_size` and only the lowest-energy conformer was used.
- **`generator3d/__init__.py`**: Fixed a `TypeError: '<' not supported between instances of 'NoneType'` crash that intermittently aborted generation when pool energies were unset, by making the energy sort None-safe and computing a final FF energy for each conformer.
- **`rmsd_utils.py`**: Fixed a `997` coordination-sphere false positive on ansa-metallocenes (TiCat1–4). The non-bonded ansa-bridge Si fell inside the distance cutoff and broke element-set equality; `calculate_tmc_rmsd` now drops from the distance-based input sphere any element absent from the bond-based generated sphere (bond sphere is donor ground truth).
- **`MetalloGenAdapter`**: Fixed TiCat2 Cp radical aromaticity loss by calling `Chem.RemoveHs(t, sanitize=False)` on sanitize-failed fragments, so the re-encoded OIN keeps aromatic `c1[cH]…` instead of kekulized `C1[CH]=…`.
- **`MetalloGenAdapter`**: Closed the stochastic TiCat1/3 `[Ti_TET]`↔`[Ti_TPY]` round-trip string drift by extending geometry-fit-ranked selection to haptic ligands. `_coordination_vectors` now reduces hapticity to centroid donors (new `_reduce_haptic_positions`: <1.6 Å transitive clustering, only when the group count equals the expected coordination number), so bent metallocenes (TiCat 14→4, ferrocene 10→2, η²-alkene 5→4) become eligible for selection instead of falling back to the lowest-energy conformer (which sometimes lands a TPY-ish embed). Strictly non-regressive (falls back to lowest-energy unless a conformer both classifies as the target *and* fits tighter); TiCat1/3 now hold `[Ti_TET]` deterministically (8/8 FF, RMSD 0.05–0.23; DEBUG confirms selection actively picks a non-lowest-energy rank), with Ferrocene/TiCp2Me2/Zeise non-regressive.

### Changed
- **`generate_3d_structures`**: Signature gained `uff_pool_size=50, rmsd_threshold=0.5, energy_threshold=2.0`; the conformer pool is now energy-sorted and deduplicated before selection.
- **`TMCOptimizer` (`clean_geometry.py`)**: `clean_geometry()`/`ff_clean()` now return `(success, final_energy)` and stamp `.energy` on the molecule (propagated through `MetalComplex.get_molecule()`), so conformers carry a rankable energy.

## [0.3.1] - 2026-07-05

### Added
- **`MetalloGenAdapter`**: Added support for MACE MLIP optimizer (`mace-omol-0-extra-large-1024`) to refine structures.
- **Force Field configuration**: Surfaced FF convergence knobs via presets, environment variables, and CLI (`--ff-preset`).

### Changed
- **`rmsd_utils.py`**: Replaced generic `999.0` return codes with distinct descending error codes (`998.0`, `997.0`, etc.) for easier debugging.
- **`MetalloGenAdapter`**: Mapped `FF` and `none` (case-insensitive) to `None` for the optimizer flag to default to FF-relaxed geometry.

### Fixed
- **OIN encoder**: Mapped binding atoms to SMILES index via canonical output order, fixing Ir `[cH]` drift.
- **`MetalloGenAdapter`**: Carried encoded sp3-carbon stereo into contract mol, fixing BDPP/BDNN round-trip failures.
- **`MetalloGenAdapter`**: Implemented manual bond-order transfer, fixing CO encoding for FeCO5/FeH2CO4.
- **`MetalloGenAdapter`**: Completed classification audit for MACE geometries, specifically handling SPY pucker and TiCat3/4 TET goldens.

## [0.3.0] - 2026-07-04

### Fixed
- **OIN v3.7: descriptor-free metal token** (`[Pt_SPL]`, was `[Pt@SP1_SPL]`). The `@desc` was an RDKit non-tetrahedral stereo leak via a stale `is_metal` variable in perception_tmc.py; isomer information was and remains fully encoded by slot ordering. Parsers continue to accept legacy `@desc` strings.
- **Eta-ring canonicalization** (`utils/oin_aligner.py`): multi-substituted η-rings now round-trip byte-stably. Fragment order for same-mass haptic ligands is keyed on a heading-independent canonical ring SMILES (was perception_tmc arrival order), and the heading/marker atom of a substituted η-ring not in `SYMMETRIC_LIGANDS` is chosen by lowest `Chem.CanonicalRankAtoms` rank (was 3D geometric alignment to the template slot vector, which varied between a hand-built and a generated structure). Winding-sign computation is untouched; symmetric-η and non-η fragments are unchanged by construction. Dead `base_sort_key` retired.
- **Square-planar Zone-A P enforcement was one-sided** (`generation/molassembler_adapter.py`, `core/chirality.py`): on SPL complexes a metal-bound P stereocenter could only ever generate one enantiomer (the other emitted the wrong 3D and a "could not be enforced" warning), because the metal-present CIP is fixed by placement geometry while the enforcement loop only re-seeded the embed. Fixed by embedding the fragment with a Z=0 dummy metal so `[P@]`/`[P@@]` embed as true 4-coordinate mirror images (symmetric with the encode side). Both enantiomers now generate correctly with no warning.
- **Bidentate incompatible-bite chelates route to the DG fallback** (`generation/molassembler_adapter.py`): DIPAMP-class ligands whose isolated conformation cannot span the chelate bite were placed on the template path and collided with the metal (non-binding H atoms landing ~1.4–1.65 Å from the metal, later misread as hydrides). A non-binding-H proximity guard now routes them to distance geometry, which round-trips them byte-identically.

### Added
- **Winding round-trip preservation (Stereo Phase 1)** (`oin/inline.py`, `generation/oin_parser.py`): the slot-tag parser now captures η-ligand winding markers (`{n>}` CW / `{n<}` CCW) and threads them through to `ParsedOIN.winding_by_slot`, so winding survives XYZ→OIN parsing into the 3D generator instead of being silently dropped.
- **Haptic-face control on 3D generation (Stereo Phase 3)** (`generation/molassembler_adapter.py`, `oin/winding.py`): the winding marker now steers which face of an η-ring the metal binds. A signed-circulation check per ring mirrors the fragment across the ring plane when its embedded winding disagrees with the marker (a proper, CIP-invariant correction).
- **Zone-A P stereocenter encoding (Stereo Phase 4)** (`core/chirality.py`): phosphorus stereocenters bonded directly to the metal are encoded as `[P@]`/`[P@@]` using a lone-pair CIP convention derived from a dummy-metal copy, and verified/enforced on regeneration. (Zone-A **N** encoding remains deferred — RDKit clears trivalent `[N@]` amine tags, so it needs an out-of-band marker.)
- **Direct Parser Fragment Mapping (v0.2.2 Blocker #1) audit completion** (2026-05-10): `_extract_oin_constraints()` audited and verified. Returns 3-tuple `(stripped_smiles, constraints_dict, fragment_to_atom_mapping)` for downstream eta-bond and polydentate-ligand processing. Fragment mapping associates OIN fragment ranks to atom indices in the connected SMILES. Renamed from public `extract_oin_constraints` to private `_extract_oin_constraints` (never a public API; 31 total call sites updated, zero unprefixed references remaining). Verification spike `tools/verify_metal_first.py` confirms metal-first invariant on 6 baseline fixtures; 3 Pd chirality test fixtures documented in `tests/fixtures/_exclusions.yml` (all verified geometrically valid via round-trip RMSD < 1.0 Å). Audit tool `tools/audit_extract_calls.py` confirms rename completeness. 55/55 tests passing (5 new fragment mapping tests verify determinism, cisplatin/polydentate correctness, contiguous atom indices, and metal-at-fragment-zero invariant). Hypergraph node `atom_direct_parser_regex` updated with new output type and status set to `clean`. MiniPRD archived.
- **Direct Parser Molassembler Instantiation audit completion** (2026-05-06): `MiniPRD_DirectParser_MolassemblerInstantiation.md` audited and verified. 20/20 unit tests passing (deterministic Cisplatin/TiCat1 construction, shape assignment, eta bond handling, error cases, all-or-nothing semantics). Implementation in `src/oinsmiles/generation/oin_parser.py` includes `construct_molassembler_mol()` (all-or-nothing transaction wrap), `convert_bond_type()` (RDKit→SCINE mapping), and `extract_oin_constraints()` (OIN v3.6 annotation extraction). SCINE shape mapping covers 10 geometries (SQP, SPL, OCT, TBP, LIN, TPL, TET, TPY, SPY, PBP). New hypergraph node `atom_direct_parser_masm` added to `architecture.yml`. MiniPRD archived.
- **Direct Parser AST Tokenization audit completion** (2026-05-06): `MiniPRD_DirectParser_ASTTokenization.md` audited and verified. 16/16 unit tests passing (deterministic atom/bond extraction, aromatic preservation, implicit H handling, error cases). Implementation in `src/oinsmiles/generation/oin_parser.py::tokenize_unsanitized_smiles()` confirmed to parse unsanitized SMILES with RDKit atom maps, preserve aromatic flags, and defer validation to Molassembler. Hypergraph node `atom_direct_parser_ast` status verified as `clean`. MiniPRD archived.
- **Direct Parser MiniPRD audit completion** (2026-05-06): `MiniPRD_DirectParser_RegexPreprocessor.md` audited and spec-aligned to OIN v3.6 inline format. All 14 unit+integration tests passing. New hypergraph node `atom_direct_parser_regex` added to `architecture.yml`. Updated 4 related MiniPRDs (AST Tokenization, Molassembler Instantiation, Integration, Verification) for consistent constraint dict keys and format examples.
- **MiniPRD audit completion** (all v0.2.0 release specs audited, 2026-05-05): All 5 core feature MiniPRDs now audited and archived — Molassembler Spike, Molassembler Adapter, Chiral Encoding, Chiral Tests, CLI. Updated MiniPRD_MolassemblerAdapter Test 5 to reflect core baseline (5 Pt/Fe/Ir complexes); v0.2.1 eta-ligand regressions documented as known limitations.
- **`GeneratedStructure` dataclass** (`generation/molassembler_adapter.py`, re-exported from `engine.py`): `OIN3DGenerator.generate()` now returns `GeneratedStructure(xyz: str, mol: Optional[Chem.Mol])` instead of a plain string. `mol` carries full RDKit bond connectivity and a 3D conformer with the template-placed positions, enabling callers to write MOL/SDF files with proper bond tables.
- **Bond-preserving MOL/SDF output** in QA test scripts: `verify_roundtrip.py` and `compare_dg_strategies.py` now use `gen_result.mol` for MOL/SDF file output so generated structures include bond connectivity (N–H, M–Cl, M–N dative bonds, etc.). XYZ-only mols are still used for RMSD calculation where matching topology is required.
- **`--include-tmqm` flag** (`verify_xyz_to_oin.py`): tmQM examples are now opt-in. The fast script (`run_verification_fast.sh`) and default roundtrip script exclude the ~103-example tmQM dataset; `run_verification_ALL.sh` includes it.
- **DG strategy comparison script** (`compare_dg_strategies.py`): benchmarks `single`, `ensemble`, and `directed` conformer strategies side-by-side on all curated examples with RMSD, min-distance, and timing metrics. Integrated into all three verification scripts.
- **`Ex{N}_{Name}_` prefixed output files**: verification scripts write named output artifacts (e.g. `Ex1_CisPlatinXYZ-OIN-SMILES_original.xyz`, `…_single.mol`, `…_generated.sdf`) for human QA.
- **Molassembler input diagnostics** in `verify_roundtrip.py` Step 2: logs parsed OIN geometry code, fragment/slot assignments, connected SMILES, permutation index, trans-sym pairs, and expected binding atoms before generation.
- **P/N stereocenter test fixtures** (`tests/integration/`): Added three Pd complex fixtures to verify chirality encoding — PdCl2-R-BINAP (axial-chiral BINAP), PdCl2-RR-BDNN (N-chiral diphosphine), PdCl2-RR-BDPP (P-chiral diphosphine). All pass round-trip verification and extend integration test coverage to 25 examples.

### Fixed
- **`cli.py` `oin2xyz` command**: Fixed `_cmd_oin2xyz` to access `.xyz` attribute from `GeneratedStructure` return value. Prior to fix, the function treated the return value as a plain string, causing `AttributeError` after `OIN3DGenerator.generate()` return type changed in v0.2.0.
- **TiCat1/3/4 3D structure generation** (`generation/molassembler_adapter.py`): `_stitch_multi_eta_fragment` was failing to generate 3D coordinates for ansa-metallocenes with aromatic η5 ligands (Cp, indenyl). Root cause: Phase 4 attempted to kekulize extracted ring SMILES (`[cH]1[cH][cH][cH][cH]1`), which fails for 5-membered all-carbon aromatic rings (5π electrons violates Hückel's 4n+2 rule). Solution: replaced with ETKDG embedding on the **full bridged fragment** (both rings + Si + methyls) with de-aromatization (aromatic bonds→SINGLE, clear aromatic flags). Phase 5 was rewritten to extract ring positions directly from the ETKDG conformer and transform them via centroid/plane alignment. Phase 7 methyl placement corrected: removed spurious `*2.0` scaling and fixed H direction sign (`-cos(tet_angle)` not `cos(tet_angle)`). Result: TiCat1/3/4 now generate 3D with correct atom counts, Si–C bonds (1.87 Å), and tetrahedral methyls. Known trade-off: de-aromatization causes round-trip bonding inference to fail (SINGLE instead of AROMATIC), but geometry quality is good (RMSD ~1.6 Å vs prior ~999 Å failures). See `docs/ETKDG_AROMATIC_FIX.md` for full technical details.

### Changed
- `OIN3DGenerator.generate()` return type changed from `str` to `GeneratedStructure`. Callers that previously used the return value as a string should access `.xyz` for the XYZ block.
- `_stitch_fragment()` and `_stitch_eta_fragment()` now return a 3-tuple `(positions, symbols, mol)` instead of a 2-tuple. `mol` is the RDKit mol with bond topology; for `_stitch_eta_fragment` it is `None` when the analytic geometry fallback is used (e.g. Cp anion ligands in ferrocene).
- `_template_generate()` return type changed from `str | None` to `tuple[str, Chem.Mol | None] | None`. Builds a combined RDKit mol by `CombineMols`-ing the metal atom and each fragment mol, adding dative metal–ligand bonds, and setting a conformer from the final `all_pos` array.

## [0.2.0] - 2026-03-07

### Added
- **SCINE Molassembler backend** (`generation/molassembler_adapter.py`): template-based 3D placement for all ligand types; DG fallback for remaining conformers. Replaces Architector entirely.
- **P/N stereocenter encoding** (`core/chirality.py`): `CIPAssigner` reads the full-TMC 3D conformer (pre-fragmentation) and stores CIP codes on P/N atoms; `ChiralityRecoveryUtility` verifies/corrects chiral tags post-fragmentation; `PseudoAtomStrategy` provides fallback for uncomputable stereocenters.
- **CLI** (`oin-smiles`): two subcommands — `xyz2oin <path>` and `oin2xyz <oin>` — registered as a package entry point.
- **`MolassemblerTimeoutError`** exported from `generation/engine.py`; `OIN3DGenerator` accepts a `timeout` parameter (default 60 s).
- **OIN v3.6 inline format** as canonical output of `XYZToSMILES.convert()` (e.g. `[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}`; the `@SP1` descriptor was a stale-variable bug, not a v3.6 design element — see v3.7 fix above).
- Integration round-trip tests (`verify_roundtrip.py`): unified XYZ → OIN → XYZ → OIN flow with RMSD < 1.0 Å and string-identity checks.
- Unit tests for chirality encoding, Molassembler adapter, regression stability, and axial chirality.

### Changed
- `OIN3DGenerator.generate()` now returns an **XYZ block string** (was an Architector `Molecule` object).
- `XYZToSMILES.convert()` now runs `CIPAssigner.assign_all()` on the full TMC mol before fragmentation.
- `pyproject.toml`: replaced Architector dependency with `scine-molassembler>=2.0.0`; added `oin-smiles` CLI entry point.
- OIN format examples in README updated from V2.4 sidecar to V3.6 inline.

### Removed
- `generation/architector_adapter.py` — Architector integration removed.
- `generation/wrapper.py` — Architector wrapper removed.
- `tests/unit/test_architector.py` — superseded by Molassembler adapter tests.

## [0.1.0] - 2025-xx-xx

Initial release with OIN v2.4 sidecar format, Architector backend, and `XYZToSMILES`/`OIN3DGenerator` APIs.
