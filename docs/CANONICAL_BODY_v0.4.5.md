# Canonical ligand body — v0.4.5 Lane 1

Two opt-in levers, both default **OFF** (unset → byte-identical output):

| lever | what it canonicalizes | code |
|---|---|---|
| `OIN_CANONICAL_BODY` | the **serialization** of a ligand body: round-trip it through `MolFromSmiles`/`MolToSmiles` and clear chelate-locked E/Z | `src/oinsmiles/oin/canonical_body.py`, hooked at `utils/xyz2mol.py` |
| `OIN_CANONICAL_PERCEPTION` | the **graph** the serializer is handed: resonance-form tie-breaks and the `AC2BO` valence walk / Kekulé matching | `utils/xyz2mol.py::lig_checks`, `utils/xyz2mol_local.py::AC2BO` |

The reusable entry point Lane 2 should call for its vertex colors is
**`oinsmiles.oin.compare.canonical_fragment_body`** (re-exported as
`oinsmiles.oin.canonical_body.canonical_body`) — the canonical form of one slot-stripped
ligand body. One implementation, three consumers (comparison key, encoder, slot
canonicalization), so they cannot drift apart.

---

## 1. Step 0 — attribution, measured before any encoder code was touched

`tools/diagnose_body_drift.py` attributes each capstone `key_equal` row to exactly one
root cause. Run over the 500 rows labelled `rdkit_canonical`:

| cause | count | share |
|---|---:|---:|
| `slot_or_order` — ligand-body multisets **identical**, only slots/order differ | **396** | 79.20% |
| `reparse_fixable` — aromatic-vs-Kekulé / implicit-vs-explicit H | **104** | 20.80% |
| `ez_chelate_locked` | 0 | 0.00% |
| `resonance_charge` | 0 | 0.00% |
| `formal_charge_placement` | 0 | 0.00% |
| `connectivity` (out of scope) | 0 | 0.00% |
| `unattributed` | 0 | 0.00% |

**`rdkit_canonical` is the fallthrough branch of `_key_equal_subclass`'s cascade.** Its
three positive tests (`fragment_reorder`, `slot_renumber`, `winding_star_drift`) each
require the drift to be *exclusively* of one kind, so a pair that mixes a renumbered slot
with a reordered fragment matches none of them and lands in `rdkit_canonical` even though
no ligand body changed at all. The headline "500 body-drift molecules" is therefore an
upper bound off by about 5×. Over all 828 `key_equal` rows the body-drift population is
**104 (12.56%)**, of which **101** would go byte-exact from the body fix alone; the other
3 also need Lane 2.

Reproduce:

```bash
PYTHONPATH=src .venv/bin/python tools/diagnose_body_drift.py \
    --bucket-report tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042/bucket_report.json
```

## 2. The measurement that changed the plan

The capstone bucket report compares a crystal geometry against a *generated* geometry, so
it can only see drift the generator happened to expose. `tools/canonicality_probe.py`
instead holds the molecular **graph fixed** and varies only presentation — random proper
rotation, random atom renumbering, both — so the expected answer is byte-identical and any
difference is a canonicality defect.

Under that instrument — 250 dataset molecules (248 encoded), 3 trials × 3 transforms each:

| arm | byte-stable | drifted | `rdkit_canonical` | `slot_renumber` | key broken |
|---|---:|---:|---:|---:|---:|
| levers off | 146 (58.87%) | 102 | 71 | 42 | **49** |
| `OIN_CANONICAL_BODY` | 152 (61.29%) | 96 | 65 | 39 | **49** |

Molecule-by-molecule: **6 fixed, 0 regressed**. The reparse is strictly non-regressive but
it fixes **zero** key-level defects.

- **`rotate` alone causes no drift.** The encoder is already orientation-invariant.
- **All drift comes from renumbering** — permuting the atom lines in the XYZ file.

That is a diagnosis rather than a disappointment. The reparse folds two *serializations*
of one graph; renumbering hands the serializer a **different graph**, and `MolToSmiles` is
faithful to whichever resonance form it is given. Only the perception lever can close it.

## 3. Step 2 — where the graph actually moves

**`lig_checks` tie-breaks (the one that mattered).** Every consumer of its candidate list —
`_select_lig_mol`'s three accumulate loops and `_rescue_unusable_perception`'s `max` —
ranks resonance forms by *(most aromatic, fewest formal charges)* and settles a tie by
taking whichever came **first**. First is `ResonanceMolSupplier` enumeration order, which
is a function of the input atom numbering. Renumbering therefore silently picks a different
resonance form: an amidinate flips `N=C(N-)` to `N-C(=N)`, a 2-iminopyridine flips its
C=N. Sorting the candidate list on `(-N_aromatic, N_pos + N_neg, canonical SMILES)` fixes
every consumer at once and changes no selection *logic* — only which member of an exact tie
wins, and now that is a property of the molecule rather than of the file.

**`AC2BO`.** It returns, in its own words, "an arbitrary resonance form": the valence walk
(`_ordered_valences` → `itertools.product` over per-atom lists in input index order) and
the Kekulé double-bond placement inside it (`get_UA_pairs` → `nx.max_weight_matching`,
whose result depends on edge insertion order) both key off atom numbering. Closed by
conjugation — relabel the atoms into a canonical order, perceive there, map the bond-order
matrix back — so every index-order dependence inside the core becomes a function of the
canonical labelling, in one place, rather than hardening each of them separately.

### ⚠ The canonical order is the SMILES write order, not `CanonicalRankAtoms`

Measured over 20 random renumberings of `CC(N)=NC`:

| source of order | invariant under renumbering? |
|---|---|
| `Chem.CanonicalRankAtoms(mol, breakTies=True)` | **no** — a different ranking 18 times in 20 |
| `Chem.MolToSmiles(mol)` → `_smilesAtomOutputOrder` | **yes** — one string every time |

`breakTies=True` settles ties between symmetry-equivalent atoms on the **input index**. The
canonical *string* is the invariant RDKit actually guarantees; the rank assignment among
symmetric atoms is not. Building the permutation on the wrong one would have produced a
"canonicalization" that is not canonical — and every single-fixture guard would still have
passed. (This is the Y2 lesson: a measurement that only exercises the easy case confirms a
wrong belief.)

## 4. Regression fixture

`tests/fixtures/NAXDOI.xyz` is the smallest in-repo structure that reproduces the defect.
Permuting its atom lines changes the emitted string **and** the comparison key:

| levers | drift over 6 presentations | key stable |
|---|---:|---|
| none | 3/6 (`rdkit_canonical`) | **no** |
| `OIN_CANONICAL_BODY` only | 3/6 — unchanged | **no** |
| `OIN_CANONICAL_PERCEPTION` | **0/6** | yes |

## 5. Known residue

- **Stereo perception is order-dependent.** With both levers on, the remaining
  key-breaking drift is dominated by `@`/`@@` flips under renumbering (QUPWUT: `[S@]` →
  `[S@@]`; OJOXAM: `[C@@H]` → `[C@H]`). That is a chirality-perception defect, not a
  body-serialization one, and it is out of Lane 1's scope — it belongs with the metal- and
  centre-stereo lanes.
- **Six of 6062 distinct corpus ligand bodies oscillate** under `canonical_fragment_body`:
  RDKit flips `@`/`@@` on adamantane-cage carbons every parse/write cycle, so there is no
  fixed point. `canonical_body_emit` bails on those and keeps the un-reparsed body. The
  same instability affects the comparison key, so it is pre-existing.
- **Donor brackets come from the inline handler, not the body.** `OINInlineHandler`
  reparses each fragment, stamps a map number on the donor (which forces a bracket) and
  keeps the bracket content, so `[cH]{0}` / `[NH2]{0}` survive whatever the body reparse
  did. This is deterministic, so it does not harm canonicality, but it does mean an emitted
  body is generally *not* equal to its own `canonical_fragment_body` — the acceptance
  predicate has to be applied to the pre-inline body, which is what
  `canonical_body_emit` guarantees by construction.
