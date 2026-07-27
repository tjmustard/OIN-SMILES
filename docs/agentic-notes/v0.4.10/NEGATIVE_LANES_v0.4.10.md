# v0.4.10 · The two chartered lanes, and why neither of them shipped code

> **Both of v0.4.10's *original* lanes are negatives.** Lane C's premise does not reproduce; Lane D's
> mechanism was already fixed by an earlier release. The two changes this release actually shipped
> were named by v0.4.9's close-out and were not in the charter at all.
>
> This is the fourth release in this project to refute part of its own plan. Writing it down is the
> process working.

---

## Lane C — SVD in `_finalize_positions` (charter Lane 1, "start here")

### What the charter claimed

> **15.6 s over 1189 SVD calls** on `HIDCIH_comp_1`, via `_finalize_positions` → `ic.update_xyz` →
> `pinv`/`svd`. Recorded at `v0.4.5-retrospective/LANE-eta-runtime-30s.md:313-314` as open item 4,
> *"the obvious next target"*, and left undone.
>
> 1189 SVD calls across 33 attempts is ~36 per attempt — **consistent with the iteration running to
> or near its cap every time.**

The charter asked which of three things was happening, since each has a different fix and only one is
"make the SVD faster":

- **(a)** converges early and keeps iterating anyway → check convergence sooner;
- **(b)** runs to the cap without converging → stop earlier, or precondition;
- **(c)** genuinely needs its iterations → cheaper decomposition, or nothing.

### What is actually happening: **none of the three**

`update_xyz` runs with `criteria=1e-4, max_iteration=30`. Instrumented by wrapping it and
`solve_equation` (one `pinv` per iteration, so the `pinv` count **is** the iteration count) and
running a full generation on three molecules spanning eta, non-eta and the CIP-bound class:

| molecule | generate | `update_xyz` calls | `pinv`/svd | % of generate | converged | hit cap | **mean iters** | B shape |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `HIDCIH_comp_1` | 20.49 s | 74 | 1.75 s | 8.5% | **74/74** | **0** | **1.00** | 55 × 171 |
| `CAHQEJ_comp_0` | 54.15 s | 90 | 8.45 s | 15.6% | **90/90** | **0** | **1.00** | 114 × 333 |
| `VAFMIA_comp_0` | 79.64 s | 1 | 0.05 s | 0.1% | **1/1** | **0** | **1.00** | 63 × 177 |

**165 calls. 165 converged. Zero reached the cap. Every single one converges on the first
iteration** — the histogram is `{1: N}` for all three molecules, with 29 of the 30 allowed
iterations never used.

So:

- **(a) is false** — it does not keep iterating; it returns immediately after the first update.
- **(b) is false** — the cap is never approached, let alone hit.
- **(c) is false in the form that matters** — it does not *need* 36 iterations, it needs **one**.

**The chartered fix is not merely unnecessary, it is impossible.** "The same result from fewer
iterations" cannot be delivered when the iteration count is already the minimum a loop can execute.

### The charter's own number does not reproduce

On `HIDCIH_comp_1`, the exact molecule the figure was taken from:

| | charter (v0.4.5-era) | measured 2026-07-27 |
|---|---:|---:|
| SVD calls | **1189** | **74** |
| attributed time | **15.6 s** | **1.75 s** |
| generation | 96.2 s | 20.49 s |

A **16× discrepancy in call count.** The figure is two releases stale: v0.4.5's duplicate-re-encode
memo, v0.4.9's changes and this release's own Lane A all sit between it and the current tree. The
charter's `⚠ this is n = 1` warning was right to be there, but the problem turned out to be that the
single measurement was **stale**, not that it was unrepresentative.

### What is left, and why it is not this release's to take

A single `pinv` per call remains, and it is not free — **15.6% of `CAHQEJ_comp_0`**. The only way to
reduce it is a **different decomposition** (`lstsq`, Cholesky where the matrix admits it). That is
precisely the class of change v0.4.10's rule forbids:

> A different decomposition can converge to coordinates that differ in the last bits, which then
> propagates through perception into a different string.

The charter names this as the lane's likeliest way to break its own rule, and it is right. A release
whose deliverable is *"only the wall-clock changes"* cannot ship a numerical substitution on the
strength of an argument about conditioning. **Handed on** to a release permitted to change numerics,
with the target correctly sized: not "1189 wasted iterations" but "one `pinv` of a 114 × 333 matrix,
8.45 s, 15.6% of one molecule and 0.1% of another."

### Ships as

**Nothing.** No lever, no `OIN_FAST_FINALIZE`, no code. The charter states plainly that *"a clean
negative is a shipped result here"*, and this is one.

---

## Lane D — redundant per-attempt work (charter Lane 2)

### What the charter claimed

`v047-boronfast` attributed the boron generation class to two unbounded mechanisms — a per-option
PuLP/CBC solve, and **a nested `EmbedMolecule` sweep that reruns every attempt** — pinned both call
sites, and deliberately implemented neither, because the paths are shared with non-boron molecules
and the fix needed a corpus A/B it had no time to run. The charter's framing: *"The blocking
condition has been removed; the work has not been done."*

Its stated deliverable was **the measurement**, not necessarily a fix: *"what is re-run per attempt,
and what it costs — stated as a measurement, not as an inference from reading the code."*

### The expensive half was already hoisted, two releases ago

`get_alternative_molecule` — the function that contains the PuLP/CBC solve — is **already memoized**:

```python
# generator3d/embed.py
def _alt_mol_cached(new_complex, option, cache):
    """Memoized get_alternative_molecule via a caller-supplied per-generation dict."""
    if cache is not None and option in cache:
        return cache[option]
```

The `alt_cache` dict is created **once per generation** (`generator3d/__init__.py:298`) and threaded
to all three call sites (`:562`, `:632`, `:656`). So the attempt-invariant preparation the charter
was aimed at is not re-run per attempt; it is computed once per `option` and cached. v0.4.9
independently corroborates this from the other direction — it measured CBC at **1.74 s of 82.44 s
(2.1%)** on `FOSNEI_comp_0` and noted *"the topology memo already collapsed it"*.

### What is genuinely re-run per attempt

Reading the loop narrows the candidates to exactly one hoistable item:

| work | attempt-invariant? | hoistable? |
|---|---|---|
| `_alt_mol_cached(...)` | yes | **already cached** |
| `get_rd_mol()` — rebuild + `SanitizeMol` from the ace_mol | **yes** | only with a defensive copy — the mol is then mutated by `_apply_double_bond_stereo` / `_apply_atom_chirality` and written into by `EmbedMolecule`, so a shared instance would leak state across attempts |
| cmap construction | **no** — depends on `haptic_scale`, which varies across the inner sweep | — |
| `AllChem.EmbedMolecule` | **no** — it *is* the attempt | not by definition |
| `alternative_ace_mol_list.index(...)` | yes, and its result was discarded | **removed — this release's Lane A** |

### The measurement — the charter's actual deliverable

Wrapping `Molecule.get_rd_mol` and `AllChem.EmbedMolecule` through a full generation on the
post-Lane-A tree:

| molecule | generate | `get_rd_mol` — **attempt-invariant**, the only hoistable candidate | `AllChem.EmbedMolecule` — per-attempt, **is** the attempt |
|---|---:|---|---|
| `CAHQEJ_comp_0` | 54.45 s | 192 calls, **0.19 s — 0.35%** | 90 calls, 25.75 s — **47.3%** |
| `FOSNEI_comp_0` | 84.41 s | 23 calls, **0.03 s — 0.03%** | 19 calls, 63.99 s — **75.8%** |

**The entire remaining attempt-invariant surface is 0.03–0.35% of a generation.** And it is not free
to take: `get_rd_mol()`'s output is immediately mutated by `_apply_double_bond_stereo` /
`_apply_atom_chirality` and then written into by `EmbedMolecule`, so hoisting it requires a
`Chem.Mol()` copy per attempt — which costs a large fraction of the 0.19 s it saves, in exchange for
a live aliasing risk in code that `2a43ce09` records as already having been broken once by a dedent
(*"restore `return []` into its except block — I dedented it and broke the embed"*).

Meanwhile **47.3% and 75.8% is `EmbedMolecule` itself**, which is per-attempt by definition. **There
is no "nested sweep rerun every attempt" left to remove: what reruns *is* the attempt.**

Note the two molecules disagree by an order of magnitude on both columns (192 vs 23 `get_rd_mol`
calls; 47.3% vs 75.8% embed share) — **the same bimodality the shipped lanes show.** It does not
change the verdict here only because both ends of the range are negligible.

### Ships as

**Nothing new.** The one measured instance of "redundant per-attempt work" in this codebase was the
discarded `.index()` scan, and Lane A removed it — 38.8% of `CAHQEJ_comp_0` against this lane's
best remaining candidate at 0.35%. The mechanism `boronfast` named was already fixed.

> **The corpus A/B `boronfast` could not run was not run here either — because there is no longer a
> change to run it on.** That is a different outcome from "still blocked", and the distinction is the
> lane's contribution.
