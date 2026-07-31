# v0.4.15 Lane 2 — `OIN_ACCEPT_STRING_EXACT` · **works, outside the class it was aimed at**

**Verdict: default OFF pending an owner call. +48 molecules (+0.96 pts), zero losses, at 4.00×
runtime.**

## The defect

`compare._parse_vertex_colors` folds reflection **deliberately**: it colours every donor atom of a
ligand with that ligand's *whole* canonical body, so two same-coloured vertices are
interchangeable and a transposition between them — an **odd** permutation, i.e. a reflection — is
invisible to `_polyhedron_signature`. And `accept_fn` decides by that key. So the generator builds
the mirror image, acceptance takes it, and the harness files the result as a benign same-isomer
string difference, in a bucket whose name says the difference is benign.

> **A lossy key must never be reused as an acceptance predicate for an axis it folds.**
> Third instance, after v0.4.8 (scored vs honest) and v0.4.11 (the donor fold).

## Why the chartered fix could not work

The charter said: use `oin/metal_config.py` as the acceptance predicate. **An acceptance test
needs a reference handedness, and the generator's only input is the OIN string.**
`_select_by_geometry_impl` already has a helicity-aware branch and it is dead by construction —
`parse_metal_config_token(parsed.original_oin)` is `None` for the whole corpus with
`OIN_EMIT_METAL_CONFIG` off, and emitting the token is v0.4.16 by the charter's own scoping.

## The replacement, chosen by measurement

Over the 183 known `MIRROR_MATCH` molecules, re-read from the baseline sweep's own
`smiles_1`/`smiles_2_indep`:

```
norm_differs = 183/183      key_same = 183/183
```

The handedness **survives `normalize_oin_for_comparison` and is folded only by the key**:

```
AFADOC_comp_0  input …c2O{5})c(O{4})c…   generated …c2O{4})c(O{5})c…
AGAVIQ_comp_0  input …P{3}(…OP{4}…)…     generated …P{4}(…OP{3}…)…
```

So the predicate is **normalized-string equality** — no new descriptor, nothing emitted.
Normalized rather than raw because raw carries the metal `@OH`/`@SP` labels the encoder documents
as atom-order-dependent and irreproducible.

**Non-regressive by construction.** An accepted conformer is returned as the *sole* pool member
(`return [early_hit]`), so plainly rejecting a key-equal conformer would let the energy-sorted pool
hand back a *different* mol. Instead the predicate returns `generator3d.ACCEPT_INCUMBENT`: the pool
keeps filling, and if nothing string-exact appears the first incumbent is returned — byte-identical
to the pre-lever answer. The whole cost is latency.

## Measured

| arm | n | gains | losses | runtime | `>30 s` |
|---|---:|---:|---:|---:|---:|
| `key_equal` | 365 | **48** | **0** | **4.00×** | 30→**122** |
| `MIRROR_MATCH` | 201 | 1 | 0 | 4.50× | 5→48 |
| control `byte_exact`/INTACT | 200 | 0 | **0** | 1.00× | 3→3 |

### 🔴 The decomposition is the finding

| subset | n | gains | rate |
|---|---:|---:|---:|
| **non-enantiomer** (`slot_renumber` / `rdkit_canonical`) | 164 | **47** | **28.7%** |
| enantiomer (`MIRROR_MATCH`) | 201 | 1 | **0.5%** |

**57× apart.** Everything the lane recovers is *outside* the class the charter aimed it at.

- **`slot_renumber` is a genuine selection bug**: the pool *does* hold a string-exact conformer and
  acceptance was stopping on a merely-key-equal one first. The lever fixes it.
- **The 201 enantiomers are not reachable by acceptance**: telemetry shows
  `pool.accept_incumbent_recorded = 1` — exactly **one** key-matching conformer in the whole pool,
  and it is the mirror. Construction must fix these, not selection.

⚠ **Scoped to the chartered 201, this release would have shipped +1 molecule.** The
owner-accepted widening to all 365 `key_equal` (2026-07-29) is the entire lane. The charter's
framing — enantiomers are the target, `slot_renumber` is a later lane — was measured backwards.

⚠ **1/201, not 0/201.** An earlier partial read 0 of 45. A single counter-example matters: the
correct handedness is not *categorically absent* from the pool, just overwhelmingly rare. That is
weaker and more accurate than "the generator cannot build it", and it changes what v0.4.16 should
attempt.

## The open decision

**+0.96 pts `byte_exact`** against **~+92 molecules over 30 s** corpus-wide (678 → ~770) and ~+3.6
CPU-h per sweep. The roadmap targets `byte_exact` 100% **and** `max(elapsed_s) < 30 s`, so this is
close to a wash between the two halves. Left default-OFF so the trade is taken deliberately.

If promoted, a full 5k sweep becomes mandatory to establish the new headline — it was skipped this
release precisely because the shipped default did not change.
