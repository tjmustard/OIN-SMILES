# Boron cages — v0.4.5

Spike verdict on the largest declared "permanent ceiling" of the `encode_fail` triage:
**34 of 48 `encode_fail` molecules are boron/carborane cages, documented as permanently
unfixable.** They are not. All 34 now encode, and re-encode consistently, behind
`OIN_BORON_CAGE` (promoted default-ON in v0.4.6).

> ⚠ **"round-trip" in this document means the NOTATION round trip, not the pipeline one.**
> Every check in §5's 34/34 table is encoder-side — encodes, byte-identical repeat encode,
> fragments re-parse, atom and bond multisets, key computable, key stable. **Not one invokes
> the 3D generator.** The original summary line here read "encode and round-trip", which a
> reader takes as OIN → XYZ → OIN; that was never measured, so the wording is corrected.
> §9 measures the pipeline arm. It is a runtime **cost**, not a win.

---

## 1. The claim, and the specific thing wrong with it

`docs/agentic-notes/v0.4.5/ENCODE_FAIL_v0.4.5.md` §5 reads:

> RDKit's 2-centre-2-electron valence model has no Lewis structure for a 3-centre-2-electron
> borane/carborane cage. `get_lig_mol`'s charge sweep already spans −4..+4 (and a diagnostic
> widening to −6..+6 finds nothing new), so no charge widening fixes this class; it needs a
> different bonding model entirely. The typed `OINEncodeError` is the correct, honest terminus.

Every sentence in that is true. The conclusion still does not follow, because **the encoder
never gets a cage to reason about.** The failure happens two stages *before* any bond-order
search, in adjacency perception:

`xyz2AC_obabel` (`utils/xyz2mol_local.py`) builds the adjacency matrix from a covalent-radius
distance criterion, and then runs this:

```python
# filter adjacency matrix if max valence is exceeded
for i in range(num_atoms):
    a_i = mol.GetAtomWithIdx(i)
    N_con = np.sum(AC[i, :])
    while N_con > max(atomic_valence[a_i.GetAtomicNum()]):
        AC = remove_weakest_bond(mol, i, AC, dMat, pt)
        N_con = np.sum(AC[i, :])
```

`atomic_valence[5] == [3, 4]`, so the cap for boron is **4**. A closo/nido deltahedral vertex
has **5 or 6** neighbours. So on every carborane this loop deletes cage bonds until each boron
is 4-coordinate. The cage shatters, the debris fragments, `get_lig_mol` fails on the debris, and
`_is_electron_deficient_cluster` (≥3 B + ≥1 B–B bond) labels the result an irreducible ceiling.

**The ceiling was a valence *rule* applied during adjacency perception.** Nothing about 3c-2e
bonding was ever tested by the failure that was attributed to it.

## 2. The measurement (`tools/boron_ac_probe.py`, 34/34)

Raw distance criterion (tolerance 0.5, no pruning) vs. the encoder's shipped AC:

| | raw (geometry only) | shipped (after pruning) |
|---|---|---|
| max boron degree | **6 or 7** (34/34) | **4** (34/34) |
| total B–B edges over the cohort | **993** | 587 |
| B–B edges deleted | — | **406 (40.9%)** |
| per molecule | — | 7–19 deleted |
| length of the deleted edges | — | **1.712–2.105 Å** |

Two things make this decisive rather than suggestive.

**The deleted edges are real bonds.** 1.71–2.11 Å is the textbook carborane B–B range
(~1.75–1.80 Å typical). The pruning loop is not trimming marginal long-range contacts; it is
deleting chemically correct cage bonds because a table says boron cannot have five neighbours.

**The raw edge counts are textbook-exact,** which is what proves the distance criterion is
recovering real cages rather than noise:

| nB | B–B edges, raw | cage this is | molecules |
|---:|---:|---|---|
| 10 | 21 | one *o*-carborane `C2B10` icosahedron (30 edges − 1 C–C − 8 C–B) | CAKBEW CAKBOG HAXJAS ICEZIC JABGAX PAQBOZ PAQCAM RANMUR RAWJEG RIWKAK RIWKEO ULODUU XUKRIF YIVLAQ |
| 20 | 42 | **two** *o*-carborane cages (2 × 21) | GANYEZ HAXJOG JAFMIP JAFTAO JAFTES MAFSIY RONQET RONQOD YIBZIV |
| 9 | 18 | one nido `C2B9` dicarbollide (25 − 1 − 6) | COZCEZ GOHWOQ PAYTUH |
| 18 | 36 | **two** dicarbollides — a COSAN-type sandwich | AVOFIB BEKLUA BEKMIP |
| 12 | 30 | closo `B12H12(2−)`, **all 30 icosahedral edges B–B** | OZAREO |
| 11 | 25 | `B11` / `C2B10` nido | MODZUA RANCIU |
| 19 | 39 | mixed dicarbollide + carborane | RONPES |
| 20 | 40 | two cages, one nido | RULBUV |

### Failure-site histogram over the 34

The spike asked where each one fails — `xyz2AC_obabel`, `AC2BO`, `SanitizeMol`, or
`MolToSmiles`. The honest answer is that the *reported* site and the *causal* site are
different stages, which is why the class was misdiagnosed:

| stage | count | what happens |
|---|---:|---|
| `xyz2AC_obabel` pruning loop | **34/34** | **causal site.** 7–19 cage bonds silently deleted. Raises nothing. |
| `MetalDisconnector` + `GetMolFrags` | 34/34 | the amputated cage now falls apart into sub-cages plus loose `[H]B` fragments, each treated as its own ligand |
| `get_lig_mol` → `OINEncodeError` | **34/34** | **reported site** (`xyz2mol.py:775`). The charge sweep is asked to find a charge for debris. |
| `AC2BO` | 0 | never reached with an intact cage — and see §5, with an intact cage it does not raise, it calls `sys.exit()` |
| `SanitizeMol` / `MolToSmiles` | 0 | never reached |

So the histogram is degenerate: **one site, 34/34.** That is the good news — one representation
change serves the whole class.

**Is the class mixed?** All 34 carry a genuine deltahedral cage (verified independently of the
`nB >= 3` bucket label, which was never evidence of one — see §7). But they span at least 5
topologies and include metallacarborane sandwiches, exo-substituted cages, and cages that are
merely a spectator ligand alongside ordinary phosphines (MODZUA is Ag(PPh₃)₂ + carborane). No
subset is an "ordinary polyhedral borane a targeted charge would fix": the charge sweep is
irrelevant because its input is already topologically wrong.

## 3. Representations tried (`tools/boron_repr_bench.py`)

Cage graph built from the unpruned AC, then benched for sanitize → serialize → **re-parse to
the same graph** → canonical-SMILES idempotence. Round-trippability was the bar, not elegance.

| representation | sanitize | serialize | re-parse | graph identical | verdict |
|---|---|---|---|---|---|
| single bonds, full sanitize | ✗ `AtomValenceException` | — | — | — | the status quo |
| **single bonds, `^ SANITIZE_PROPERTIES`** | **OK** | **OK** | **OK** | **yes** | **shipped** |
| single bonds, `^ PROPERTIES ^ KEKULIZE` | OK | OK | OK | yes | equivalent; redundant |
| B–B `DATIVE`, full sanitize | ✗ `AtomValenceException` | — | — | — | dative still counts toward explicit valence |
| B–B `UNSPECIFIED`, full sanitize | mixed | OK | ✗ | **no** | **not round-trippable** |
| B–B `ZERO`, full sanitize | mixed | OK | ✗ | **no** | **not round-trippable** |
| single + pinned explicit Hs, full sanitize | ✗ | — | — | — | valence check is on the graph, not the Hs |

The two negative results are the load-bearing ones, because dative/zero-order bonds were the
obvious candidates:

- **`DATIVE` does not help at all.** RDKit counts a dative bond toward the *end* atom's explicit
  valence, so a 6-connected cage boron trips the same check.
- **Zero-order bonds move the failure downstream instead of fixing it.** They sanitize, and
  `MolToSmiles` writes them as `~`. But re-parsing `~` yields a *mix* of `SINGLE` and
  `UNSPECIFIED` — on OZAREO, 30 uniform `ZERO` B–B bonds came back as 19 `SINGLE` + 11
  `UNSPECIFIED`. SMILES cannot carry a zero-order bond. A representation that encodes but does
  not read back has not solved anything.

What works is the boring option: **plain single bonds, and skip the one check that objects.**
That is the sense in which the graph was always fine — the cage needs no exotic bond type, only
for RDKit's valence *table* to be stepped around on the atoms that violate it.

## 4. What shipped, and how it is scoped

Everything is behind `OIN_BORON_CAGE`, default OFF, and every gate additionally requires the
molecule to *contain the motif*.

**The motif is a B–B–B triangle** (`boron_cage_vertices`, `utils/xyz2mol_local.py`) — the
deltahedral face signature. Chosen over "contains ≥3 boron" so the relaxation cannot reach
ordinary boron chemistry. Pinned by unit test:

| species | B–B–B triangle? | touched? |
|---|---|---|
| `BPh4-`, `BF4-` borates | no B–B bond at all | no |
| diboron / diboryl `B–B` | no third boron | no |
| boroxine (`B-O-B-O-B-O`) | zero B–B bonds | no |
| linear `B–B–B` chain | 2 B–B bonds, no triangle | no |
| any closo/nido cage vertex | yes | yes |

This is **stricter** than the pre-existing `_is_electron_deficient_cluster` (≥3 B and ≥1 B–B
bond), which would also match the linear chain.

Four changes:

1. **`xyz2AC_obabel` pruning exemption** — cage vertices skip the connectivity cap, computed
   from the *pre*-pruning AC so the exemption cannot be triggered by pruning itself.
2. **`_cage_frag_mol` (`utils/xyz2mol.py`)** — perceives a cage fragment *directly* rather than
   searching for bond orders: every cage edge single, hydrogens as the geometry gives them,
   formal charges zero, `SANITIZE_ALL ^ SANITIZE_PROPERTIES`. It deliberately does **not** claim
   a Lewis structure, and does **not** claim the chemically-correct cage charge (a dicarbollide
   is really 2−, closo-B₁₂H₁₂ is 2−); it claims only a canonical, round-trippable graph. See §6.
3. **`sanitize_allowing_boron_cage` (`utils/aromaticity.py`)** — the three downstream full
   sanitizes a cage trips (`GetMolFrags(asMols=True)`, `kekulize_safe_sanitize`,
   `CIPAssigner.assign_all`) route through it. With the lever unset, or on a mol without the
   motif, it is **exactly `Chem.SanitizeMol(mol)`** — same call, same exceptions. If the relaxed
   sanitize also fails, the strict error is re-raised, so nothing is silently swallowed.
4. **`compare._parse_fragment` rung** — without it a cage fragment falls to the `RAW:` fallback
   and contributes its *literal input SMILES* to `canonical_roundtrip_key`, making the key
   atom-order dependent. With it, the cage fragment canonicalizes properly.

### Why `AC2BO` is bypassed rather than taught about boron

The tempting one-liner is `atomic_valence[5] = [3, 4, 5, 6]`. Do not. With an over-connected
boron, `AC2BO` hits:

```python
if not possible_valence:
    logger.debug(...)
    sys.exit()
```

A bare **`sys.exit()` in a library perception path** — `SystemExit` is a `BaseException`, so
`except Exception` cannot catch it and the whole process dies. Widening the table changes which
molecules reach that line; it does not remove the landmine. `_cage_frag_mol` routes cage
fragments around `AC2BO` entirely. (The `sys.exit()` is a latent hazard for any other
over-coordinated element and is worth filing separately; it is not touched here.)

## 5. Round-trip proof — the acid test (`tools/boron_roundtrip.py`)

`OIN_BORON_CAGE=1`, 34/34, isolated subprocess per molecule:

| check | result |
|---|---|
| encodes | **34/34** |
| byte-identical on a repeat encode | **34/34** |
| every OIN fragment re-parses | **34/34** |
| heavy-atom multiset == encoder's own `tmc_mol` | **34/34** |
| heavy-**bond** multiset (element-pair) == encoder's | **34/34** |
| hydrogens on boron conserved | **34/34** |
| `canonical_roundtrip_key` computable | **34/34** |
| key stable across a repeat | **34/34** |
| key free of the `RAW:` fallback | **34/34** |

Worked example, OZAREO (closo-B₁₂H₁₂ amide on Rh, 62 atoms):

```
[Rh_TPY].CN(C)C(O{2})NB1234[BH]{0}567[BH]89%10[BH]%11%12%13[BH]58%14[BH]%1158
[BH]%12%11%15[BH]9%13%12[BH]{3}16%10[BH]2%11%12[BH]35%15[BH]47%148
.Cc{1}1c{1>}(C)c{1}(C)c{1}(C)c{1}1C
```

All 12 borons present, the icosahedron carried as ring-closure digits, 11 × `[BH]` plus the one
exo-substituted vertex written bare. Re-parsing gives back 29/29 heavy atoms, identical
heavy-bond multiset, 11/11 boron hydrogens.

### One honest caveat on "same atom count"

Total-H is **not** used as the pass criterion, and the reason matters. Stripping a slot marker
from an OIN fragment turns a coordinated donor into a free ligand, which then legitimately fills
implicit hydrogens (`C(O{2})` → `C(O)`, i.e. an −OH). Over the 34 the total-H delta ranges 0–25.
That drift is a property of the OIN format's slot notation, **not** of the cage: it appears on
ordinary passing molecules too (see `tools/boron_regression_ab.py`'s donor-H histogram) and the
per-boron hydrogen count is exact 34/34. Charging it to the cage work would be wrong. The
primary criteria are the heavy-atom multiset, the heavy-bond multiset, and boron-H — the three
things a shattered cage would fail immediately.

## 5a. The blast radius is not 34 — 14 "passing" molecules are silently wrong

This came out of the regression A/B, and it is the most consequential finding in the spike.

The A/B (`tools/boron_regression_ab.py`, 120 passing molecules, seed 0) returned:

| | result |
|---|---|
| lever OFF vs the frozen capstone OIN | **120/120 byte-identical** — the change is inert when off |
| lever ON vs lever OFF | **119/120 identical**, **1 differs** |

The one difference is `VEJXOZ_comp_0` — and it is not a regression, it is an **additional fix**
that exposes a second, worse failure mode:

| | B–B cage bonds | spurious bonds | key |
|---|---|---|---|
| geometry (truth) | 12 | — | — |
| lever OFF | **6** (50% deleted) | invents a **C=B double bond** to balance valences | falls back to `RAW:` |
| lever ON | **12** | none | canonical |

`VEJXOZ` was scored as a *pass*. Its OIN round-trips. It describes the wrong molecule. The
round-trip key cannot see this, because the corrupted encode is compared against **its own
corrupted mol** — the same "a lossy key must never be the acceptance predicate for an axis it
folds" trap this release has hit before.

So the 34 `encode_fail` molecules are only the subset where amputating the cage happened to
produce something `get_lig_mol` could not perceive *at all*. Where the debris happens to remain
perceivable, the encoder silently emits a plausible, self-consistent, wrong graph.

Corpus scan for the actual population (`tools/boron_blast_radius.py`, no encoding, no
generation — a text filter to ≥3 boron, then adjacency on those 192 files):

| | count |
|---|---:|
| xyz files in `cat/` + `photo/` | 26,230 |
| with ≥3 boron | 192 |
| carrying a real deltahedral cage motif | **186** |
| of those, cage bonds **deleted** by the pruning loop | **186 — every single one** |
| ├─ known `encode_fail` 34 | 34 |
| ├─ **counted as PASSING in the frozen capstone reports** | **14** |
| └─ not covered by the capstone arm (unmeasured) | 138 |
| boron but no cage motif (borates etc., correctly untouched) | 6 |

The 14 silently-corrupted passers lose **133 of 269 cage bonds (49.4%)**:

`PEKQUU` (17/34 deleted) · `RAJNEY` (12/21) · `ULOFIK` (11/21) · `DUDTIG` (10/18) ·
`KIXXOF` (10/18) · `RAJNOI` (10/21) · `XIQKOY` (10/18) · `UYEJAK` (9/21) · `XIQLAL` (9/18) ·
`PEKQII` (8/16) · `VOFHUW` (8/21) · `CIDHAY` (7/18) · `SEMTOV` (6/12) · `VEJXOZ` (6/12)

All 14 also **round-trip cleanly with the lever on** (`tools/boron_roundtrip_14passing.json`,
14/14 `ROUNDTRIP_OK`: deterministic, all fragments re-parse, heavy-atom + heavy-bond multisets
and boron-H exact, key canonical with no `RAW:`). Two of them — `KIXXOF` and `DUDTIG`, both Rh
thiaboranes — needed §5b's crash fix first; before it they killed the worker process outright.

Read plainly: **the pruning defect reaches 186 corpus molecules. 34 fail loudly, 14 fail
silently while being scored correct, and 138 were never measured.** The "34 permanent ceiling"
number was both a misdiagnosis and an undercount, and the accuracy metric was reporting 14
wrong answers as right.

## 5b. The one genuinely nasty obstacle: a cage chiral tag is a native crash

Worth its own section because it is the only thing in this spike that could not have been
predicted from reading code, and because "it raises an exception" would have been the wrong
answer.

`AssignAtomChiralTagsFromStructure` stamps a permutation tag on a 5-/6-connected cage vertex,
because its 3D neighbourhood *is* asymmetric. RDKit's stereo machinery has no permutation table
for that shape, and it does not report the problem — it corrupts the heap. Observed, all with
`OIN_BORON_CAGE=1` before the fix:

| molecule | symptom |
|---|---|
| `KIXXOF` (Rh thiaborane) | `RuntimeError: basic_string::_M_create` from `Chem.AssignStereochemistry` |
| `KIXXOF`, encoded twice in one process | `free(): invalid pointer`, `Fatal Python error: Aborted`, inside `FindPotentialStereo` |
| `DUDTIG` (Rh thiaborane) | `free(): invalid size` → SIGABRT; separately SIGSEGV |

Three things about this that matter:

1. **It manifests on the *second* encode.** The first encode of `KIXXOF` succeeded and returned a
   correct OIN. The abort came on the next one. That is latent heap corruption, so a
   single-molecule test would have passed and the corpus sweep would have died at a random point
   with no attributable cause.
2. **`except Exception` cannot help.** A `SIGABRT` is not catchable, so the only correct fix is to
   never set the tag. `clear_boron_cage_stereo` (`core/chirality.py`) clears it in *both*
   `CIPAssigner.assign_all` and `ChiralityRecoveryUtility.recover` — `recover` is entered from
   `get_oin_string` on a mol `CIPAssigner` never touched, so neither can rely on the other.
3. **It is not a boron problem, it is a *cage vertex* problem.** The tag that aborted
   `FindPotentialStereo` on `KIXXOF` was on the thiaborane's **cage sulfur**. Clearing only boron
   would have left the crash in place. Cage heteroatom vertices (carborane C, thiaborane S) are
   identified as atoms bonded to ≥3 cage-vertex borons — a bound an exocyclic substituent
   (bonded to one) cannot reach.

Clearing loses nothing: a cage vertex's "handedness" is the polyhedron, and the polyhedron is
already carried by the cage's bond graph. Verified by encoding four cage molecules five times
each in one process — no crash, all deterministic, and no `[B@` anywhere in the output.

## 6. What this does *not* claim

Being precise, because the encode is real but bounded:

- **Not a Lewis structure.** The cage is carried as a plain single-bonded graph. It asserts
  connectivity, nothing about electron count. That is the honest reading of an object whose
  bonding a 2c-2e model genuinely cannot express — and it is enough for a canonical 1D hash,
  which is what OIN is for.
- **Not the chemically-correct cage charge.** `_cage_frag_mol` returns charge 0, so the derived
  metal oxidation state for these complexes is the "neutral cage" reading, not the −2 a chemist
  would write for a dicarbollide. Fixing that means deciding a charge convention for cages; it is
  a separate call and would change the metal token.
- **Not proven stable under renumbering.** Same limitation the `ASISAX` fix hit: these encode
  deterministically for their actual atom ordering, and the key is now canonical rather than
  `RAW:`, but the corpus-wide renumbering instability
  (`docs/agentic-notes/v0.4.5/RENUMBERING_INSTABILITY_v0.4.5.md`) has not been re-measured for cages. Lanes 1/2/8 own
  that gap; this spike did not duplicate it.
- **Zone-A CIP degrades on cage molecules.** `_build_dummy_metal_copy` (`core/chirality.py`)
  builds its probe with a full sanitize and so fails on a cage, emitting the existing
  `OINStereoWarning` and falling back to today's clearing behaviour. Graceful, already-designed
  degradation, not a crash — but on MODZUA it means the phosphine lone-pair CIP is not computed.
  A further call site to route through `sanitize_allowing_boron_cage` if anyone wants it closed.
- **No stereo is carried on the cage itself.** Per §5b every cage-vertex chiral tag is cleared,
  which is correct (RDKit cannot represent that shape, and trying corrupts memory) but does mean
  a genuine cage-substitution diastereomer would not be distinguished by an `@` on the vertex.
  Whether the cage bond graph plus the slot markers already separate such isomers is untested —
  it would need a mirror-twin collision probe of the kind Y1 built, and this spike did not do one.

## 7. Ruled out, with evidence

- **"obabel's distance criterion misses the long cage B–B bonds"** — refuted. It finds every one
  of them, at 1.71–2.11 Å, in textbook-exact counts. The loss is entirely downstream in the
  pruning loop.
- **"No charge widening fixes this class"** — confirmed, and irrelevant. The fragments handed to
  the charge sweep are already topologically wrong, so no charge on them could be right. Widening
  −4..+4 to −6..+6 was measuring the wrong stage.
- **"It needs a different bonding model entirely (multi-centre bonds)"** — refuted. It needs a
  plain single-bonded graph and one skipped sanitize flag. No multi-centre bond is involved.
- **DATIVE and zero-order bonds** — tried, both fail, §3. Recorded so nobody re-chases them.
- **The `34` was never earned as a cage count.** `classify()` in `tools/sl5_triage.py` buckets on
  `nB >= 3` read off the raw xyz, with no B–B check — unlike the runtime detector. It happens to
  be right for all 34 (independently verified here), but the bucket label was not evidence.

## 8. Reproduce

```bash
export PYTHONPATH=$PWD/src
V=/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python
D=$PWD/tmCAT-tmPHOTO_xyz_dataset

# the reframe: geometry perceives the cage, the pruning loop deletes it (34/34, ~1 min)
$V tools/boron_ac_probe.py --dataset-dir $D

# representation A/B: single/^PROPERTIES round-trips, DATIVE and zero-order do not
$V tools/boron_repr_bench.py --dataset-dir $D --mols AVOFIB,OZAREO,CAKBEW

# the acid test, 34/34 ROUNDTRIP_OK (~1 min CPU)
OIN_BORON_CAGE=1 $V tools/boron_roundtrip.py --dataset-dir $D

# where each one fails, stage by stage
$V tools/boron_characterize.py --dataset-dir $D

# regression A/B: lever OFF vs frozen, and lever ON vs OFF, on passing molecules
$V tools/boron_regression_ab.py --dataset-dir $D \
   --reports-dir $D/results-capstone-v042/individual_reports --sample 120 --seed 0

$V -m unittest tests.unit.test_boron_cage tests.unit.test_regression_stability
```

## 8a. Suite, lint, and exactly which commit each number came from

| gate | result | measured at |
|---|---|---|
| `uvx ruff@0.15.20 check` + `format --check`, whole tree | clean | final (`262f83da`) |
| `tests/unit/test_regression_stability.py`, 4 goldens | **6/6 OK**, lever OFF *and* ON, byte-identical across the lever | final |
| `tests/unit/test_boron_cage.py` | **19/19 OK** | final |
| 4 cage molecules x 5 encodes in one process | no crash, all deterministic, no `[B@` in output | final |
| all 6 corpus molecules with >=3 B and no cage motif | **6/6 byte-identical** with the lever ON | final |
| **full `discover tests/unit`** | **624 tests, OK (skipped=3, expected failures=3)** — the full collected count, one clean invocation | **final tree** (`tools/suite_clean.txt`) |
| ├─ pre-existing tests (non-boron) | **605, all OK** — exactly the pre-change baseline | final |
| └─ new lever tests | **19, all OK**; 0 load failures | final |
| 120-molecule regression A/B (§5a) | 120/120 OFF==frozen, 119/120 OFF==ON | **final src** (re-run; bit-identical to the earlier arm) |

### A measurement of mine that was not clean, kept on the record

Before the clean run above, the suite reported **`Ran 623 tests, OK`** against the same src — one
short of the 624 the loader collects. The cause was mine: **I edited
`tests/unit/test_boron_cage.py` while that run was in flight**, so it executed 605 pre-existing
tests plus 18 of what are now 19 boron tests. `discover` imports at collection time, so an
in-flight edit yields a count matching neither tree, and the discrepancy is a *single digit* —
exactly the size that gets rounded away rather than chased.

It is superseded by the clean 624/OK run and is recorded only for the lesson: **do not edit a test
file while a suite run is in flight in the same worktree.** Check the reported count against
`TestLoader().discover(...)` rather than trusting it.

The **120-molecule A/B was never affected by this** — it imports `src`, not `tests`. It was re-run
against the final src and came back bit-for-bit identical to the earlier arm: 120/120 OFF==frozen,
119/120 OFF==ON (`VEJXOZ` again), same donor-H histogram
`{-4: 2, 0: 29, 1: 22, 2: 37, 3: 13, 4: 8, 5: 5, 6: 4}`. Confirmed by
`git diff <crash-fix commit> HEAD -- src` being empty, so the src the A/B imported is the src
shipped here. That the numbers did not move is expected — `clear_boron_cage_stereo` is gated on
both the lever and the B-B-B motif, so it cannot reach the OFF arm — but it is now measured rather
than argued from the gate.

The **120-molecule A/B is not affected by any of this** — it was re-run against the final src and
came back bit-for-bit identical to the earlier arm: 120/120 OFF==frozen, 119/120 OFF==ON
(`VEJXOZ` again), same donor-H histogram `{-4: 2, 0: 29, 1: 22, 2: 37, 3: 13, 4: 8, 5: 5, 6: 4}`.
Confirmed by `git diff <crash-fix commit> HEAD -- src` being empty, so the src the A/B imported is
the src shipped here. That the numbers did not move is the expected result —
`clear_boron_cage_stereo` is gated on both the lever and the B-B-B motif, so it cannot reach the
OFF arm at all — but it is now measured rather than argued from the gate.

## 9. Verdict on the ceiling

Sharper than "the model cannot do it":

> **The 34 boron cages were never blocked by RDKit's bonding model. They were blocked by
> `xyz2AC_obabel`'s valence-capped bond-pruning loop, which deletes 406 of 993
> geometrically-correct B–B cage bonds (40.9%, at 1.71–2.11 Å) because
> `atomic_valence[5]` tops out at 4. With cage vertices exempted from that cap and the cage
> fragment carried as a plain single-bonded graph sanitized with `SANITIZE_PROPERTIES` skipped,
> all 34 encode and round-trip: deterministic, every fragment re-parsing, heavy-atom and
> heavy-bond multisets and boron-H counts exact, and `canonical_roundtrip_key` canonical rather
> than falling back to a raw string.**

The irreducible part is much smaller than a ceiling, and it is a *convention* question rather
than a representational one: a single-bonded cage graph is not a Lewis structure, and its formal
charge (currently 0) is a decision nobody has made. Neither blocks a canonical, lossless 1D
hash, which is what OIN actually requires.

**Addressable count for the `encode_fail` cohort: 34 of 48 move from "permanent ceiling" to
"encodes, and re-encodes consistently, behind a lever."** Combined with §7 of
`ENCODE_FAIL_v0.4.5.md` (4 already work, 1 fixed, 3 deferred, 7 timeout-bound), the
"confirmed unfixable" row of that table is now **0** — *for encoding*.

And the correct total is **48 molecules, not 34**: the 34 loud failures plus the 14 silent ones
from §5a, all 48 now encoding and notation-round-tripping (34/34 and 14/14 `ROUNDTRIP_OK`).
The 14 are the more interesting half, because they were being counted as *successes*.

⚠ Read §10 before quoting any of this as a round-trip result. Encoding was the ceiling this
spike set out to break, and it broke it. The **generator** is a separate ceiling, it was not
measured here, and it is not passing.

---

## 10. The pipeline arm, measured (2026-07-26) — the lever's unpriced runtime cost

§5's table is nine encoder-side checks. This is the arm nobody ran: does a boron molecule that
now *encodes* also **generate a 3D structure**? Sample of 10 of the 34, `optimizer=None`,
`ensemble_size=1`, generation cap 60 s (`scratchpad/boron_gen_times.jsonl`):

> ⚠ **CORRECTED 2026-07-26 — the first version of this section reported a 10-molecule sample as
> "0 of 10 produce a 3D structure". With 33 of the 34 measured, it is 2/33: two DO assemble.** The
> sample was not representative, and the error mattered — it was the basis of a proposed "boron
> fast-fail", which this measurement REFUTES, because a blanket fast-fail would cost those 2 real
> passes. Same lesson as the rest of this wave: a sample that only exercises the common case
> confirms whatever you already believed.

| outcome | n | detail |
|---|---|---|
| **produced a 3D structure** | **2** | `RAWJEG` (LIN, 2 slots), `ULODUU` (TET, 4 slots) |
| **burned the whole cap, produced nothing** | 25 | 60.0 – 172.8 s, all `failed to generate any conformers via OIN-direct assembly` |
| instant, loud failure | 3 | `Geometry code 'NON' not supported by MetalloGen mapping`, 0.00–0.01 s |
| instant, loud failure | 3 | `UncoordinatedFragmentError` |

**2 of 33 assemble.** So the promotion moves most — not all — of this class from failing instantly
to failing slowly. The 25 cap-burners cost **~2.1 CPU-hours per full 5,000-molecule sweep at the
300 s budget, for zero passes**, which is a real waste; but it cannot be reclaimed by a blanket
fast-fail without losing the 2 that work, and no clean discriminator between the two groups was
found (both span small and large cages; the 2 successes are low-denticity LIN/TET, the failures
run to OCT/SPY at 5-6 slots, which is suggestive but n=2).

A precise gate — fast-fail only above some denticity or cage size — needs the mechanism, not this
correlation. Deferred to v0.4.6+.

`XIQKOY_comp_0` is the clean two-point demonstration, same molecule, one lever:

| `OIN_BORON_CAGE` | encode | generate | total |
|---|---|---|---|
| `0` | 0.86 s → OIN keeps the cage as a **disconnected** fragment | `UncoordinatedFragmentError` in **0.01 s** | **0.87 s** |
| `1` | 1.12 s → OIN is a correct, fully-coordinated B₁₀ cage (5 slots) | never returns | **>340 s** |

The lever is doing its job perfectly: with it on, the encoder emits a *better* string — a
genuinely coordinated closo-borane cage instead of an amputated fragment. The generator then
cannot assemble that cage, and it does not fail fast; it spends the entire per-molecule budget
discovering this.

### Why this matters beyond boron

1. **It works directly against the <30 s per-molecule goal.** At the sweep's 300 s budget, 34
   molecules that previously cost ~1 s now cost up to 300 s each — roughly **2.8 CPU-hours
   added to a full sweep for zero additional passes.** `levers.py` prices the promotion at "14
   molecules move from scored-passing to failing"; the runtime line was missing, and it is now
   in that entry.
2. **`timeout` is not a hard bound.** Every capped molecule overran: 60 s requested, 60.7–137.9 s
   spent (GOHWOQ 2.3×). `embed_time_budget=self.timeout` bounds the embed attempt loop, not the
   OIN-direct assembly path around it. What actually enforces the budget in a sweep is the
   harness's per-molecule SIGKILL subprocess — so **any timing measured without that watchdog
   understates the tail**, and any caller relying on `timeout` as a wall-clock guarantee is
   wrong.
3. **This does not argue for reverting the promotion.** A notation that describes the right
   molecule and fails loudly in the generator is still better than one that silently describes
   the wrong graph (§5a's 14). The honest framing is: `OIN_BORON_CAGE` fixed the **encoder**
   ceiling and exposed a **generator** ceiling that was previously hidden behind an encode
   failure. Assembling a polyhedral borane cage from m-SMILES is the open problem, and it is a
   generator3d problem, not a notation one.

**Reproduce:**

```bash
V=.venv/bin/python; export PYTHONPATH=$PWD/src
# the two-point demonstration
for v in 0 1; do OIN_BORON_CAGE=$v timeout 340 $V -u -c "
import sys,time; from oinsmiles import XYZToSMILES
from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen as G
t=time.monotonic(); o=XYZToSMILES().convert(sys.argv[1]); print('encode',round(time.monotonic()-t,2),o)
t=time.monotonic()
try: G(optimizer=None,ensemble_size=1,timeout=300,ff_params=None).generate(o); print('gen ok')
except Exception as e: print('gen fail',round(time.monotonic()-t,2),type(e).__name__)
" tmCAT-tmPHOTO_xyz_dataset/cat/XIQKOY_comp_0.xyz; done
```
