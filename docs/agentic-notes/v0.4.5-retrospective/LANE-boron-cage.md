# Lane — Boron cages (`OIN_BORON_CAGE`)

**The failure population:** the **34 `boron_cluster` molecules** of the 48-molecule `encode_fail`
cohort — crystal structures whose XYZ produced no OIN string at all — which had been written up
twice as a *permanent representational ceiling*. They are not a ceiling. And the population is
not 34: another **14 molecules were being scored as PASSING while the encoder described the wrong
graph**, plus **138** more that carry the same defect and were never measured.

This lane produced the single largest directly-measured accuracy gain of the release.

---

## ELI5

In ordinary chemistry a bond is two atoms sharing two electrons, which you can draw as one line
between two dots. A borane "cage" (a carborane, a dodecaborate) does not work that way: it shares
two electrons across **three** atoms at once, so every corner of the cage touches five or six
neighbours instead of the three or four a boron atom is "allowed". Our structure reader looks at
the raw 3D coordinates, correctly notices all those close contacts, and then consults a lookup
table that says *boron may have at most 4 neighbours* — and to make that true it **deletes** the
longest bonds until each boron is down to 4. The cage falls apart into rubble, the rubble cannot
be interpreted as a molecule, and the resulting error was blamed on "chemistry we can't
represent". The actual fix is to stop deleting: a cage is a perfectly ordinary *graph* (dots and
lines), it just breaks one bookkeeping rule, so we skip that one rule for cage corners only. No
new chemistry was needed at any point.

## The work, visually

```
  crystal XYZ (62 atoms, closo-B12H12 amide on Rh)
        │
        ▼
┌──────────────────────── utils/perception_core.py::xyz2AC_obabel ─────────────────────────┐
│                                                                                        │
│  STEP 1  covalent-radius distance criterion   d(i,j) <= Rcov_i + Rcov_j + tolerance     │
│          ───────────────────────────────────────────────────────────────────────────    │
│          ✔ PERCEIVES THE CAGE CORRECTLY.  max boron degree 6 or 7 (34/34 molecules);   │
│            993 B-B edges over the cohort, topologies textbook-exact                     │
│            (o-carborane 10 B / 21 edges · closo-B12 12 / 30 · dicarbollide 9 / 18)      │
│                                                                                        │
│  STEP 2  the pruning loop      ◄── THE CAUSAL SITE, 34/34.  Raises nothing.            │
│          ─────────────────                                                             │
│            for i in cap_order:                                                         │
│         ┌─►   while N_con > max(atomic_valence[Z_i]):   # atomic_valence[5] == [3, 4]   │
│         │        AC = remove_weakest_bond(...)          # deletes the LONGEST bond      │
│         │                                                                              │
│    ●────┘   LEVER INTERVENES HERE:                                                     │
│    ● if lever_enabled("OIN_BORON_CAGE"):                                               │
│    ●     exempt = boron_cage_vertices(atomic_nums, AC)   # from the PRE-pruning AC      │
│    ● for i in cap_order:  if i in exempt: continue                                     │
│                                                                                        │
│          ✘ WITHOUT IT: 406 of 993 B-B edges deleted (40.9%), 7-19 per molecule,        │
│            at bond lengths 1.712 - 2.105 A  (textbook carborane B-B is ~1.75-1.80 A)   │
└────────────────────────────────────────────────────────────────────────────────────────┘
        │
        ▼   ┌─────────────────────────── two outcomes ───────────────────────────┐
        │   │                                                                    │
   SHATTERED CAGE                                                        INTACT CAGE
        │                                                                        │
        ├─ MetalDisconnector + GetMolFrags (34/34)                               │
        │    sub-cages + loose [H]B, each treated as its own ligand              │
        │                                                                        │
        ├─► PATH A (34 mols) get_lig_mol charge sweep -4..+4 on DEBRIS           │
        │      → OINEncodeError "electron-deficient boron cluster"               │
        │        REPORTED SITE perception_tmc.py:959 · NO OIN STRING                    │
        │                                                                        │
        └─► PATH B (14 mols) the debris happens to stay perceivable              │
               → encoder INVENTS a C=B double bond to balance valences           │
               → emits a plausible, self-consistent, WRONG graph                 │
               → round-trip key compares it against its OWN corrupted mol        │
               → ★ SCORED AS A PASS ★                                            │
                                                                                 │
                              ●  LEVER ON  ────────────────────────────────────────┘
                              ●  get_lig_mol → _has_boron_cage → _cage_frag_mol
                              ●    every cage edge SINGLE · charges 0 · H as geometry gives
                              ●    sanitize with SANITIZE_ALL ^ SANITIZE_PROPERTIES
                              ●  aromaticity.py::sanitize_allowing_boron_cage  (3 call sites)
                              ●  chirality.py::clear_boron_cage_stereo  (the native-crash pin)
                              ●  compare.py::_parse_fragment cage rung  (key ≠ "RAW:")
                                    │
                                    ▼
                     OIN string · 48/48 ROUNDTRIP_OK · key canonical
```

Legend — `●` = code the `OIN_BORON_CAGE` lever adds or gates. `✔` = works today and always did.
`✘` = the defect. `★` = the finding that changed the product decision. `AC` = adjacency matrix.
`cap_order` = the iteration order of the pruning loop (input order by default; heaviest-Z-first
under the *separate* `OIN_STABLE_METAL_AC` lever, which lives in the same loop).

## Initial assumptions and hypothesis

**The headline assumption, and it was wrong: this class was a permanent ceiling.** It had been
recorded that way twice, in prose that is individually true and collectively a non-sequitur:

- `docs/agentic-notes/v0.4.4/ENCODER_ROBUSTNESS_v0.4.4_SL5.md` §W1 — *"RDKit cannot perceive a 3c-2e boron cage into a
  sanitizable Lewis structure, and `get_lig_mol`'s charge sweep already spans −4..+4, so no charge
  widening will ever encode these… This classifies 34/48 of the cohort (it does not *encode* boron
  cages — out of scope)."*
- `docs/agentic-notes/v0.4.5/ENCODE_FAIL_v0.4.5.md` §5, titled **"Confirmed unfixable: 34 boron clusters"** — *"it
  needs a different bonding model entirely (multi-center bonds), out of scope for a valence-graph
  encoder. The typed `OINEncodeError` (already landed) is the correct, honest terminus: **an
  encoder that refuses this input is correct**, not a bug."*
- The same belief is still in code, in `perception_tmc.py::_is_electron_deficient_cluster`'s docstring:
  *"this is a permanent representational ceiling of the RDKit valence model, not a missed charge
  guess."*

Every sentence in those is true about 3c-2e bonding. **None of it was ever exercised by the
failure attributed to it**, because the encoder never received a cage to reason about. The
ceiling was a valence *rule* applied during adjacency perception, two stages before any
bond-order search.

Four further hypotheses were held going in, and all four are answered under "Dead ends":

1. *"obabel's distance criterion misses the long cage B–B bonds"* — the natural explanation for
   why perception would produce a broken cage.
2. *"the charge sweep just needs widening"* — the fix the previous two write-ups had already
   tried and rejected (−4..+4 widened diagnostically to −6..+6).
3. *"a cage needs an exotic bond type"* — `DATIVE`, or zero-order (`ZERO` / `UNSPECIFIED`) B–B
   bonds, were the obvious candidates for "a bond RDKit will not count against valence".
4. *"the blast radius is the 34 loud failures"* — i.e. this is a pure encode-coverage lane with no
   correctness dimension.

One assumption that held: **the `34` itself.** All 34 do carry a genuine deltahedral cage,
verified independently — but see the last dead end, because the *bucket label* was never the
evidence for it.

## What was actually found

### Confirmed — the cause is the pruning loop, one site, 34/34

`tools/boron_ac_probe.py` (artifact `tools/boron_ac_probe.json`, 34 rows) compares the raw
distance criterion (tolerance 0.5, no pruning) against the encoder's shipped AC. Re-derived from
the artifact for this write-up:

| | raw (geometry only) | shipped (after pruning) |
|---|---|---|
| max boron degree | **6 or 7**, 34/34 | **4**, 34/34 |
| total B–B edges over the cohort | **993** | 587 |
| B–B edges deleted | — | **406 (40.9%)** |
| per molecule deleted | — | **7 – 19** |
| length of the deleted edges | — | **1.712 – 2.105 Å** (406 edges) |

Two things make this decisive rather than suggestive. The deleted edges are at textbook carborane
B–B distances, so the loop is not trimming marginal long-range contacts. And the raw edge counts
are topologically exact — 10 B / 21 B–B edges is exactly one *o*-carborane C₂B₁₀ icosahedron
(30 edges − 1 C–C − 8 C–B); 12 B / 30 edges is closo-B₁₂H₁₂ with all 30 icosahedral edges B–B;
9 B / 18 edges is one nido C₂B₉ dicarbollide. The distance criterion is recovering real cages,
not noise.

Failure-site histogram over the 34 (`tools/boron_characterize.py`): the *reported* site and the
*causal* site are different stages, which is exactly why the class was misdiagnosed.

| stage | count | what happens |
|---|---:|---|
| `xyz2AC_obabel` pruning loop | **34/34** | **causal site.** 7–19 cage bonds silently deleted. Raises nothing. |
| `MetalDisconnector` + `GetMolFrags` | 34/34 | the amputated cage falls into sub-cages plus loose `[H]B`, each treated as its own ligand |
| `get_lig_mol` → `OINEncodeError` | **34/34** | **reported site**, `perception_tmc.py:959`. The charge sweep is asked to find a charge for debris. |
| `AC2BO` | 0 | never reached with an intact cage — and with one it does not raise, it calls `sys.exit()` |
| `SanitizeMol` / `MolToSmiles` | 0 | never reached |

A degenerate histogram is the good news: one representation change serves the whole class.

### Confirmed — 48/48 encode and round-trip with the lever on

`tools/boron_roundtrip.py`, isolated subprocess per molecule. Artifacts:
`tools/boron_roundtrip_34.json` (34 rows, **34/34 `ROUNDTRIP_OK`**, `deterministic: true` in all
34) and `tools/boron_roundtrip_14passing.json` (14 rows, **14/14 `ROUNDTRIP_OK`**).

Per-molecule checks, all 34/34: encodes; byte-identical on a repeat encode; every OIN fragment
re-parses; heavy-atom multiset equals the encoder's own `tmc_mol`; heavy-**bond** multiset
(element-pair) equals the encoder's; hydrogens on boron conserved; `canonical_roundtrip_key`
computable; key stable across a repeat; key free of the `RAW:` fallback.

Worked example, OZAREO (closo-B₁₂H₁₂ amide on Rh, 62 atoms) — all 12 borons present, the
icosahedron carried as ring-closure digits, 11 × `[BH]` plus the one exo-substituted vertex
written bare:

```
[Rh_TPY].CN(C)C(O{2})NB1234[BH]{0}567[BH]89%10[BH]%11%12%13[BH]58%14[BH]%1158
[BH]%12%11%15[BH]9%13%12[BH]{3}16%10[BH]2%11%12[BH]35%15[BH]47%148
.Cc{1}1c{1>}(C)c{1}(C)c{1}(C)c{1}1C
```

**Total-H is deliberately not a pass criterion**, and the reason is not the cage: stripping a slot
marker from an OIN fragment turns a coordinated donor into a free ligand which then legitimately
fills implicit hydrogens (`C(O{2})` → `C(O)`, an −OH). Over the 34 the total-H delta ranges 0–25,
and the same drift appears on ordinary passing molecules — the 120-molecule A/B's donor-H
histogram is `{-4: 2, 0: 29, 1: 22, 2: 37, 3: 13, 4: 8, 5: 5, 6: 4}`. The primary criteria are the
heavy-atom multiset, the heavy-bond multiset and boron-H: the three things a shattered cage fails
immediately.

### ★ Confirmed — 14 molecules were scored as PASSING while describing the WRONG GRAPH

This is the most consequential finding in the lane, and it came out of the regression A/B rather
than the coverage work.

`tools/boron_regression_ab.py` (120 passing molecules, seed 0; artifact
`tools/boron_regression_ab.json`):

| arm | result |
|---|---|
| lever OFF vs the frozen capstone OIN | **120/120 byte-identical** — the change is inert when off |
| lever ON vs lever OFF | **119/120 identical**, **1 differs** |

The one difference is `VEJXOZ_comp_0`, a nido-C₂B₇ on Ru, and it is not a regression — it is an
additional fix that exposes a second, worse failure mode:

| | B–B cage bonds | spurious bonds | key |
|---|---|---|---|
| geometry (truth) | 12 | — | — |
| lever OFF | **6** (50% deleted) | **invents a C=B double bond** to balance valences | falls back to `RAW:` |
| lever ON | **12** | none | canonical |

`VEJXOZ` was scored a *pass*. Its OIN round-trips. It describes the wrong molecule. The round-trip
key cannot see this, because the corrupted encode is compared against **its own corrupted mol** —
the same "a lossy key must never be the acceptance predicate for an axis it folds" trap this
release hit elsewhere.

So the 34 `encode_fail` molecules are only the subset where amputating the cage happened to produce
something `get_lig_mol` could not perceive *at all*.

Corpus scan for the real population (`tools/boron_blast_radius.py`, artifact
`tools/boron_blast_radius.json` — a text filter to ≥3 boron, then adjacency on those 192 files;
no encoding, no generation):

| | count |
|---|---:|
| xyz files in `cat/` + `photo/` | 26,230 |
| with ≥3 boron | 192 |
| carrying a real deltahedral cage motif | **186** |
| of those, cage bonds **deleted** by the pruning loop | **186 — every single one** (verified row-by-row: `BB_deleted > 0` in 186/186) |
| ├─ known `encode_fail` 34 | 34 |
| ├─ **counted as PASSING in the frozen capstone reports** | **14** |
| └─ not covered by the capstone arm (unmeasured) | 138 |
| boron but no cage motif (borates etc., correctly untouched) | 6 |

The 14 silently-corrupted passers lose **133 of 269 cage bonds (49.4%)**:

`PEKQUU` (17/34 deleted) · `RAJNEY` (12/21) · `ULOFIK` (11/21) · `DUDTIG` (10/18) ·
`KIXXOF` (10/18) · `RAJNOI` (10/21) · `XIQKOY` (10/18) · `UYEJAK` (9/21) · `XIQLAL` (9/18) ·
`PEKQII` (8/16) · `VOFHUW` (8/21) · `CIDHAY` (7/18) · `SEMTOV` (6/12) · `VEJXOZ` (6/12)

All 14 round-trip cleanly with the lever on. Read plainly: **the pruning defect reaches 186 corpus
molecules — 34 fail loudly, 14 fail silently while being scored correct, and 138 were never
measured.** The "34 permanent ceiling" number was both a misdiagnosis and an undercount, and the
accuracy metric was reporting 14 wrong answers as right.

**This is the explicit product argument that promoted the lever, and it is not a free win.**
Promoting `OIN_BORON_CAGE` moves 14 molecules from scored-passing to **failing**. That is correct —
it trades 14 silent false positives for 14 loud honest failures — but it means a headline pass rate
can move either way, and nobody reading a single aggregate number would be able to tell which
happened. The decision recorded at `oin/levers.py::_DEFAULT_ON`: *a lossless notation that silently
emits a wrong graph is worse than one that fails audibly.*

### Confirmed — a chiral tag on a cage vertex is a NATIVE CRASH, not a bad descriptor

The only obstacle in the lane that could not have been predicted from reading code, and "it raises
an exception" would have been the wrong answer.

`AssignAtomChiralTagsFromStructure` stamps a permutation tag on a 5-/6-connected cage vertex,
because its 3D neighbourhood genuinely *is* asymmetric. RDKit has no stereo permutation table for
that shape and does not report the problem — it corrupts the heap. Observed, all with
`OIN_BORON_CAGE=1` before the fix:

| molecule | symptom |
|---|---|
| `KIXXOF` (Rh thiaborane) | `RuntimeError: basic_string::_M_create` from `Chem.AssignStereochemistry` |
| `KIXXOF`, encoded twice in one process | `free(): invalid pointer`, `Fatal Python error: Aborted`, inside `FindPotentialStereo` |
| `DUDTIG` (Rh thiaborane) | `free(): invalid size` → SIGABRT; separately SIGSEGV |

Three properties that matter:

1. **It manifests on the *second* encode.** `KIXXOF`'s first encode succeeded and returned a
   correct OIN; the abort came on the next one. That is latent heap corruption, so a
   single-molecule test passes and a corpus sweep dies at a random point with no attributable
   cause.
2. **`except Exception` cannot help.** `SIGABRT` is not catchable, so the only correct fix is to
   never set the tag.
3. **It is not a boron problem, it is a *cage vertex* problem.** The tag that aborted
   `FindPotentialStereo` on `KIXXOF` was on the thiaborane's **cage sulfur**. Clearing only boron
   would have left the crash in place.

Verified by encoding four cage molecules five times each in one process: no crash, all
deterministic, no `[B@` anywhere in the output.

### Confirmed — the shipped representation is the boring one

`tools/boron_repr_bench.py` built the cage graph from the unpruned AC and benched
sanitize → serialize → **re-parse to the same graph** → canonical-SMILES idempotence.
Round-trippability was the bar, not elegance.

| representation | sanitize | serialize | re-parse | graph identical | verdict |
|---|---|---|---|---|---|
| single bonds, full sanitize | ✗ `AtomValenceException` | — | — | — | the status quo |
| **single bonds, `^ SANITIZE_PROPERTIES`** | **OK** | **OK** | **OK** | **yes** | **shipped** |
| single bonds, `^ PROPERTIES ^ KEKULIZE` | OK | OK | OK | yes | equivalent; redundant |
| B–B `DATIVE`, full sanitize | ✗ `AtomValenceException` | — | — | — | dative counts toward explicit valence |
| B–B `UNSPECIFIED`, full sanitize | mixed | OK | ✗ | **no** | not round-trippable |
| B–B `ZERO`, full sanitize | mixed | OK | ✗ | **no** | not round-trippable |
| single + pinned explicit Hs, full sanitize | ✗ | — | — | — | the valence check is on the graph, not the Hs |

## What was done

Everything is behind `OIN_BORON_CAGE`, and **every gate additionally requires the molecule to
contain the motif** — so the lever is not "relax valence checking", it is "relax valence checking
on deltahedral cage vertices".

**The motif is a B–B–B triangle** — `utils/perception_core.py::boron_cage_vertices(atoms, AC) ->
set[int]`, the deltahedral face signature. Every closo/nido vertex sits on at least one triangular
face; nothing else in this corpus does. Chosen over "contains ≥3 boron" so the relaxation cannot
reach ordinary boron chemistry, and it is **stricter than the pre-existing
`_is_electron_deficient_cluster`** (≥3 B and ≥1 B–B bond), which also matches a linear chain.
Negative controls, each a unit test in `tests/unit/test_boron_cage.py::TestCageMotifDetection`:

| species | B–B–B triangle? | detector output |
|---|---|---|
| `BPh4-` borate | no B–B bond at all | `set()` |
| `BF4-` | no B–B bond at all | `set()` |
| diboron / diboryl `B–B` (`OB(O)B(O)O`) | no third boron | `set()` |
| boroxine `B1OB(O)OB(O)O1` (B-O-B-O-B-O) | three borons, **zero** B–B bonds | `set()` |
| linear three-boron chain (`OB(O)B(O)B(O)O`) | 2 B–B bonds, **no triangle** | `set()` |
| any closo/nido cage vertex (`B1B(O)B1O`) | yes | 3 vertices |

Six code changes:

1. **`utils/perception_core.py::xyz2AC_obabel` — the pruning exemption.**
   `exempt = boron_cage_vertices(atomic_nums, AC)` is computed from the **pre-pruning** AC, so the
   exemption cannot be triggered by pruning itself, and the loop does `if i in exempt: continue`.
   This is the causal fix; the other five make the rest of the pipeline survive an intact cage.

2. **`utils/perception_tmc.py::_cage_frag_mol(frag_mol)` — perceive the cage directly.**
   Called from `get_lig_mol` *before* `_select_lig_mol`, gated on
   `lever_enabled("OIN_BORON_CAGE") and _has_boron_cage(mol)`. Every cage edge single, hydrogens as
   the geometry gives them, formal charges zero, `SANITIZE_ALL ^ SANITIZE_PROPERTIES`.
   `_has_boron_cage` (`perception_tmc.py:624`) is the fragment-level wrapper around the same
   `boron_cage_vertices` motif, so both halves of the lever agree on what a cage is.

   **Why `AC2BO` is bypassed rather than taught about boron.** The tempting one-liner is
   `atomic_valence[5] = [3, 4, 5, 6]`. Do not. With an over-connected boron, `AC2BO` reaches
   `if not possible_valence: … sys.exit()` — a bare `sys.exit()` in a library perception path.
   `SystemExit` is a `BaseException`, so `except Exception` cannot catch it and the whole process
   dies. Widening the table changes *which* molecules reach that line; it does not remove the
   landmine. (That `sys.exit()` is a latent hazard for any other over-coordinated element and is
   still unfixed — see open questions.)

3. **`utils/aromaticity.py::sanitize_allowing_boron_cage(mol)`** — the three downstream full
   sanitizes a cage trips (`GetMolFrags(asMols=True)`, `kekulize_safe_sanitize`,
   `CIPAssigner.assign_all`) route through it. With the lever off, or on a mol without the motif,
   it is **exactly `Chem.SanitizeMol(mol)`** — same call, same exceptions. If the relaxed sanitize
   also fails, the strict error is re-raised, so nothing is silently swallowed. The predicate is
   `_boron_cage_relaxation_applies(mol)`: lever set **and** ≥3 boron **and** a B–B–B triangle.

4. **`core/chirality.py::clear_boron_cage_stereo(mol)`** — the native-crash pin. Clears the chiral
   tag on every cage vertex, in **both** `CIPAssigner.assign_all` and
   `ChiralityRecoveryUtility.recover`, because `recover` is entered from `get_oin_string` on a mol
   `CIPAssigner` never touched — neither can rely on the other. Cage **heteroatom** vertices
   (carborane C, thiaborane S) are identified as atoms bonded to **≥3 cage-vertex borons**, a bound
   an exocyclic substituent (bonded to one) cannot reach. Clearing loses nothing: a cage vertex's
   "handedness" is the polyhedron, and the polyhedron is already carried by the cage's bond graph.

5. **`oin/compare.py::_parse_fragment` — the cage rung.** Without it a cage fragment falls to the
   `RAW:` fallback and contributes its *literal input SMILES* to `canonical_roundtrip_key`, making
   the key atom-order dependent. The rung retries with `_NO_VALENCE` (`SANITIZE_ALL ^
   SANITIZE_PROPERTIES`) and then `_NO_VALENCE_NO_KEKULIZE`. **See the blast-radius leak below —
   this rung is where it lived.**

6. **Lever plumbing.** `OIN_BORON_CAGE` originally used the bare-truthiness read
   `os.environ.get("OIN_BORON_CAGE")` at **five separate sites**, which means `OIN_BORON_CAGE=0`
   *ENABLED* it — anyone opting out the obvious way got the opposite of what they asked for.
   Migrated to the central registry `oin/levers.py::lever_enabled`, where `0` / `false` / `no` /
   `off` / `""` disable.

### ⚠ The blast-radius leak, found and closed at promotion time

`compare.py::_parse_fragment`'s cage rung skips `SANITIZE_PROPERTIES` — that is a valence-**RULE**
bypass, not a cosmetic relaxation. Gated on the lever **alone**, it applied that bypass to **every
fragment**. `C#O` fails the valence check and nothing else, so **carbon monoxide started PARSING**
instead of reaching the `RAW:` fallback — and CO is among the commonest ligands in
transition-metal chemistry. The promotion was silently changing chemistry for a large population
that has nothing to do with boron.

Caught by exactly one suite failure out of 840:
`tests/unit/test_canonical_body.py::test_unparseable_body_gets_stable_raw_token`, failing with
`AssertionError: 'C#O' != 'RAW:C#O'`.

Fix: also require boron in the fragment —
`if _lever_enabled("OIN_BORON_CAGE") and ("B" in smiles or "b" in smiles)`. Deliberately a cheap
substring test rather than a parse: this runs before any mol exists, boron's presence in the SMILES
is a *necessary* condition for a cage, a false positive (the `B` inside `Br`) merely restores the
previous unscoped behaviour for that one fragment, and a false negative is impossible because a
cage cannot be written without a boron symbol.

Verified over **all 1,194 distinct fragment bodies the corpus emits** (harvested from
`smiles_1` / `smiles_2` across the 936-molecule re-baseline, slot markers and `_GEO` stripped),
comparing the parse result lever-ON vs lever-OFF:

| | |
|---|---:|
| distinct emitted fragment bodies | **1,194** |
| fragments whose parse result differs ON vs OFF | **56** |
| …of those, containing boron | **56** — the intended scope |
| …of those, boron-FREE | **0** — the leak is closed |

All 56 go `None → parsed` (an 8-atom `[BH]1B[BH][B@H]2…` cage, a 46-atom carborane thioether),
which is the recovery the lever exists for. Guard:
`tests/unit/test_boron_cage.py::TestValenceBypassIsScopedToBoron::test_carbon_monoxide_still_falls_back_with_the_lever_ON`.

Two process notes worth keeping. First, this instrument **replaced** a killed 61-fixture
whole-string byte-identity run and is better on every axis: n=1194 not 61, it exercises the one
predicate that changed rather than a downstream proxy, it runs in seconds rather than ~80 minutes,
and it is load-independent. Second, the A/B had been deferred as "needs an idle machine once the
5k sweep frees the cores" — **which was wrong**: string equality is deterministic, only wall-clock
is load-sensitive. Deferring a load-independent measurement on load-dependent grounds cost real
time.

## Dead ends and refutations

### "obabel's distance criterion misses the long cage B–B bonds" — REFUTED

The natural explanation for a broken cage, and the reason the ceiling story was believable.
**Killed by `tools/boron_ac_probe.py`:** the raw criterion finds **every** cage edge — 993 over the
cohort, at 1.712–2.105 Å, in topologically exact counts (o-carborane 21, closo-B₁₂ 30,
dicarbollide 18). Max boron degree raw is 6 or 7 in 34/34 and 4 in 34/34 after pruning. The entire
loss is downstream, in the pruning loop.

### "No charge widening fixes this class" — CONFIRMED, and irrelevant

The one prior hypothesis that was *true*. `get_lig_mol`'s sweep spans −4..+4 and a diagnostic
widening to −6..+6 finds nothing new. **But it was measuring the wrong stage:** the fragments handed
to the charge sweep are already topologically wrong, so no charge on them could be right. A
confirmed negative on the wrong stage is how a misdiagnosis acquires evidence.

### "It needs a different bonding model entirely (multi-centre bonds)" — REFUTED

**Killed by the shipped fix existing.** It needs a plain single-bonded graph and one skipped
sanitize flag. No multi-centre bond is involved anywhere in the implementation.

### "A dative bond will dodge the valence check" — REFUTED

**Killed by `tools/boron_repr_bench.py`:** RDKit counts a dative bond toward the *end* atom's
explicit valence, so a 6-connected cage boron trips the identical `AtomValenceException`. `DATIVE`
does not help at all.

### "A zero-order bond will dodge the valence check" — REFUTED, and this one is subtler

Zero-order (`ZERO`) and `UNSPECIFIED` B–B bonds *do* sanitize, and `MolToSmiles` writes them as
`~`. They move the failure downstream instead of fixing it. **Killed by the re-parse arm:**
re-parsing `~` yields a *mix* — on OZAREO, 30 uniform `ZERO` B–B bonds came back as **19 `SINGLE` +
11 `UNSPECIFIED`**. SMILES cannot carry a zero-order bond. A representation that encodes but does
not read back has not solved anything, and for a *lossless notation* project that is the whole bar.

### "Just widen `atomic_valence[5]`" — REFUTED, and it is worse than not working

**Killed by reading what `AC2BO` does when `possible_valence` is empty:** a bare `sys.exit()`.
`SystemExit` is a `BaseException`, so widening the table does not make the encoder handle cages —
it changes which molecules walk into an uncatchable process kill.

### "The blast radius is the 34 loud failures" — REFUTED, and this is the lane's main result

**Killed by `tools/boron_regression_ab.py` then `tools/boron_blast_radius.py`:** 120/120
byte-identical OFF-vs-frozen but 119/120 OFF-vs-ON, and the one difference (`VEJXOZ`) turned out
to be a *scored pass* with 6 of 12 cage bonds deleted and an invented C=B double bond. The corpus
scan then found 186 affected molecules: 34 loud, **14 silently wrong while scored correct**, 138
unmeasured.

### "The `34` was a cage count" — REFUTED as *evidence*, though the number is right

`classify()` in `tools/sl5_triage.py` buckets on `nB >= 3` read off the raw xyz, with **no B–B
check** — unlike the runtime detector. It happens to be right for all 34 (verified independently
here), but the bucket label was never evidence that any of them contained a cage. A count that is
correct for the wrong reason still cannot be reasoned from.

### A hand-written cage SMILES as a scoping test — ABANDONED

Tried as the positive half of `TestValenceBypassIsScopedToBoron` and withdrawn; per the test
module's own note it *"only proved that hand-written cage SMILES are easy to get wrong."* The
positive half is instead carried by real fixtures:
`TestPruningExemption::test_lever_on_keeps_the_full_icosahedron` (30/30 B–B edges survive on
OZAREO) and `TestCageEncodes::test_key_does_not_degrade_to_the_raw_fallback`.

### A measurement of my own that was not clean, kept on the record

Before the final clean run, `discover tests/unit` reported **`Ran 623 tests, OK`** against the same
`src` — one short of the 624 the loader collects. The cause was mine: **`tests/unit/test_boron_cage.py`
was edited while that run was in flight**, so it executed 605 pre-existing tests plus 18 of what
are now 19 boron tests. `discover` imports at collection time, so an in-flight edit yields a count
matching neither tree — and the discrepancy is a *single digit*, exactly the size that gets rounded
away rather than chased. Superseded by the clean 624/OK run. Lesson: do not edit a test file while
a suite run is in flight in the same worktree, and check a reported count against
`TestLoader().discover(...)` rather than trusting it. The 120-molecule A/B was never affected — it
imports `src`, not `tests`, and was re-run against final `src` bit-for-bit identically.

## Where it landed

**Lever:** `OIN_BORON_CAGE`. **`PROMOTED TO DEFAULT-ON IN v0.4.6`, not v0.4.5** — merge commit
**`d799de1f`** (`release(v0.4.6): boron cage promotion + Lane 5 metal Delta/Lambda (P1)`), listed
in `oin/levers.py::_DEFAULT_ON` alongside the six v0.4.5 canonicality levers.

⚠ **Through v0.4.5 it shipped OFF *by omission*.** It appeared in **neither** `_DEFAULT_ON` nor
`_HELD_OFF` in `src/oinsmiles/oin/levers.py` — the one state that leaves no recorded reason. A
lever that is off because someone wrote down why is a decision; a lever that is off because nobody
listed it is an accident that looks like a decision.

**The promotion's measured justification** (recorded at the `_DEFAULT_ON` definition site and in
`CHANGELOG.md` §0.4.6): on the **936-molecule re-baseline**, **34 of the 36 `XYZToSMILES failed`
rows are electron-deficient boron clusters**, and the lever takes that population from **0/36
encoding to 34/36**, at **0.2–4.2 s each**. The two non-boron holdouts are a quinoid-ring case and
an `ASISAX` `get_lig_mol` case (both owned by the `encode_fail` lane). The boron lane separately
measured **48/48 round-tripping**.

**The promotion's cost, recorded in the same place:** 14 molecules move from scored-passing to
failing. Correct, and it means a headline pass rate can move either way.

**Lane commits** (`swimlane/v045-boron`, merged to `release/v0.4.5` at **`f4c3525a`**, then into
the release integration merge **`1450b5ce`**):

| commit | what |
|---|---|
| `5ca3f6e9` | `boron(v0.4.5): the cage ceiling is an AC pruning rule, not the valence model` — the pruning exemption, `_cage_frag_mol`, `sanitize_allowing_boron_cage`, the `_parse_fragment` rung |
| `2008026c` | docs: the cage ceiling, measured — 34/34 encode and round-trip |
| `ea4a6ceb` | `boron: a cage chiral tag is a NATIVE CRASH, and 14 "passing" molecules are wrong` — `clear_boron_cage_stereo` |
| `72d8cef8` | docs: the silent-corruption blast radius and the native-crash obstacle |
| `ae96aafb` | 48/48 round-trip on final code + the frozen measurement artifacts |
| `8172bdd9` | `test(boron): the sharpest scoping control -- boron-rich, cage-free, must not move` |
| `262f83da` | chore: ignore the run logs the measurement tools drop in `tools/` (the commit the final gate numbers were measured at) |
| `eb5b7fb0`, `1d8514c9`, `30b34edb`, `924727c1` | docs corrections: which commit each number came from; the in-flight-edit suite miscount; the A/B re-run; the clean 624/OK |
| `d799de1f` (v0.4.6) | promotion to default-ON **+ the `_parse_fragment` boron scoping fix** |

**Guard tests** — `tests/unit/test_boron_cage.py`, 19 tests at v0.4.5, extended at promotion:

| class | pins |
|---|---|
| `TestCageMotifDetection` | the six motif cases above — borate, `BF4-`, diboron, boroxine, linear B–B–B chain emit `set()`; a triangle emits 3 |
| `TestPruningExemption` | `test_lever_off_amputates_the_cage` (max B degree 4, B–B < 30 on OZAREO) / `test_lever_on_keeps_the_full_icosahedron` (max degree 6, **exactly 30** B–B) |
| `TestCageEncodes` | typed `OINEncodeError` with the lever off; encodes + deterministic with it on; all 12 borons present; `test_key_does_not_degrade_to_the_raw_fallback`; `test_mixed_complex_encodes_cage_and_ordinary_ligands` (MODZUA = Ag(PPh₃)₂ + carborane — the phosphines must still come out as ordinary aromatic fragments) |
| `TestCageStereoMustNeverBeTagged` | `KIXXOF` encodes without crashing; no `[B@` survives |
| `TestSilentCorruptionOfAPassingMolecule` | `VEJXOZ`: 6 B–B + ≥1 spurious C=B with the lever off, **12 B–B + 0 spurious** with it on |
| `TestNonCageMoleculesUnaffected` | four goldens byte-identical across the lever; **plus the sharpest control** — `ASUVIV` (Ir boroxine, 3 B, 0 B–B) and `AROTAE` (Fe with `[BH3-]`/`[BH-]`/borane groups) byte-identical across the lever. A golden has no boron at all, so it cannot show whether the gate is the *motif* or merely the *element*; these two can. |
| `TestValenceBypassIsScopedToBoron` | `C#O` still returns `None` from `_parse_fragment` with the lever ON |

`_LeverMixin` sets `OIN_BORON_CAGE` **explicitly in both directions**, because deleting the
variable used to mean "off" and stopped meaning that at promotion — every `test_lever_off_*`
silently became a second lever-ON test asserting amputated-cage behaviour against the fixed path.
That was the third and fourth occurrence of the trap in this release; it is now a lint,
`tests/unit/test_levers.py::TestNoTestUnsetsAPromotedLever`.

**Other gates, all at the final tree:** `uvx ruff@0.15.20 check` + `format --check` clean;
`tests/unit/test_regression_stability.py` (4 goldens) **6/6 OK lever OFF *and* ON**, byte-identical
across the lever; 4 cage molecules × 5 encodes in one process — no crash, deterministic, no `[B@`;
all 6 corpus molecules with ≥3 B and no cage motif **6/6 byte-identical** with the lever ON; full
`discover tests/unit` **624 tests OK (skipped=3, expected failures=3)** — 605 pre-existing (exactly
the pre-change baseline) + 19 new, 0 load failures.

### ⚠ Superseded text still in the tree — flag for whoever edits these next

| location | stale claim |
|---|---|
| `docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md` | says "default OFF" throughout, and §4's opening *"Everything is behind `OIN_BORON_CAGE`, default OFF"*. True for v0.4.5, **false since `d799de1f`**. |
| `docs/agentic-notes/v0.4.5/ENCODE_FAIL_v0.4.5.md` §5 + §7 | *"Confirmed unfixable: 34 boron clusters"*, *"needs a different bonding model entirely"*, *"34 of 48 (70.8%) are an honest, correct ceiling"*. **Refuted by this lane** — the doc's own §7 header row "confirmed unfixable" should read 0. |
| `docs/agentic-notes/v0.4.4/ENCODER_ROBUSTNESS_v0.4.4_SL5.md` §W1 | the "irreducible ceiling" framing. Its test module docstring was corrected at promotion; the doc was not. |
| `src/oinsmiles/utils/perception_tmc.py::_is_electron_deficient_cluster` docstring | *"this is a permanent representational ceiling of the RDKit valence model"*. False for the default configuration. |
| `src/oinsmiles/utils/perception_core.py::xyz2AC_obabel` comment | closes with *"Both default OFF; with neither set this loop is byte-identical to pre-v0.4.5"* — both `OIN_BORON_CAGE` and `OIN_STABLE_METAL_AC` are now default-ON. |
| `src/oinsmiles/utils/perception_tmc.py::_cage_frag_mol` docstring | *"Returns: (mol, 0) on success"* — it returns the bare mol; the caller adds the `, 0`. |
| checked-in artifacts | `tools/boron_roundtrip.json` (2 rows) and `tools/boron_characterize.json` (2 rows: AVOFIB, BEKLUA) are **small samples**, not the 34-molecule results. The 34-molecule round-trip artifact is `tools/boron_roundtrip_34.json`; the 14 passers are in `tools/boron_roundtrip_14passing.json`; the AC probe's 34 rows are in `tools/boron_ac_probe.json`. Do not cite the 2-row files as cohort evidence. |

## Open questions / for the next agent

1. **The 138 unmeasured cage molecules.** They carry the identical defect (cage bonds deleted in
   186/186 rows) but were outside the capstone arm, so nobody knows how many are silent
   false-positives like the 14 and how many are loud failures like the 34. This is the single
   biggest unknown the lane leaves, and it is a cheap sweep now that the lever is default-ON.
2. **Cage formal charge is an unmade decision, not a bug.** `_cage_frag_mol` returns charge 0, so
   the derived metal oxidation state for these complexes is the "neutral cage" reading, not the
   −2 a chemist would write for a dicarbollide or closo-B₁₂H₁₂. Fixing it means choosing a charge
   convention for cages — and it would change the `[M_XXX]` metal token, so it needs its own
   geometry-tag veto (`tools/geometry_tag_shift.py`).
3. **Not proven stable under renumbering.** These encode deterministically for their actual atom
   ordering and the key is canonical rather than `RAW:`, but the corpus-wide renumbering
   instability (`docs/agentic-notes/v0.4.5/RENUMBERING_INSTABILITY_v0.4.5.md`) has **not** been re-measured for cages.
   Run `tools/canonicality_probe.py` over the 186.
4. **Zone-A CIP degrades on cage molecules.** `core/chirality.py::_build_dummy_metal_copy` builds
   its probe with a full sanitize, so it fails on a cage, emits the existing `OINStereoWarning` and
   falls back to clearing. Graceful, already-designed degradation — but on MODZUA it means the
   phosphine lone-pair CIP is not computed. One more call site to route through
   `sanitize_allowing_boron_cage` if anyone wants it closed.
5. **No stereo is carried on the cage itself.** Every cage-vertex chiral tag is cleared, which is
   correct (RDKit cannot represent that shape and trying corrupts memory), but a genuine
   cage-substitution diastereomer would not be distinguished by an `@` on the vertex. Whether the
   cage bond graph plus the slot markers already separate such isomers is **untested** — it needs a
   mirror-twin collision probe of the kind Y1 built (`tools/injectivity/`).
6. **The `sys.exit()` in `AC2BO` is still there** and is a latent hazard for any over-coordinated
   element, not just boron. `_cage_frag_mol` routes cages around it; nothing else is protected.
   Worth filing separately.
7. **The stale text in the table above.** In particular, `ENCODE_FAIL_v0.4.5.md`'s "confirmed
   unfixable: 34" is now the most misleading sentence in the docs tree, because it is the row a
   planner would read when deciding what to work on next.
