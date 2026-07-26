# Boron cages — v0.4.5

Spike verdict on the largest declared "permanent ceiling" of the `encode_fail` triage:
**34 of 48 `encode_fail` molecules are boron/carborane cages, documented as permanently
unfixable.** They are not. All 34 now encode and round-trip behind `OIN_BORON_CAGE`
(default OFF).

---

## 1. The claim, and the specific thing wrong with it

`docs/ENCODE_FAIL_v0.4.5.md` §5 reads:

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
  (`docs/RENUMBERING_INSTABILITY_v0.4.5.md`) has not been re-measured for cages. Lanes 1/2/8 own
  that gap; this spike did not duplicate it.
- **Zone-A CIP degrades on cage molecules.** `_build_dummy_metal_copy` (`core/chirality.py`)
  builds its probe with a full sanitize and so fails on a cage, emitting the existing
  `OINStereoWarning` and falling back to today's clearing behaviour. Graceful, already-designed
  degradation, not a crash — but on MODZUA it means the phosphine lone-pair CIP is not computed.
  A fourth call site to route through `sanitize_allowing_boron_cage` if anyone wants it closed.

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
"encodes and round-trips behind a default-OFF lever."** Combined with §7 of
`ENCODE_FAIL_v0.4.5.md` (4 already work, 1 fixed, 3 deferred, 7 timeout-bound), the
"confirmed unfixable" row of that table is now **0**.
