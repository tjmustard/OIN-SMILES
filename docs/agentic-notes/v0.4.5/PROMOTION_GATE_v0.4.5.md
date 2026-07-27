# v0.4.5 promotion gate — measured result and recommendation

Measured on `trial/v045-integration` (all six lane branches merged into local `main`), with
`tools/canonicality_probe.py --n 300 --trials 2`, seed 42 fixed so **every arm samples the same
molecules**. Generator-free: the probe holds the molecular graph **fixed** and varies only proper
rotation, atom renumbering, and both — so the correct answer is byte-identical, a known ground
truth rather than an inferred baseline.

---

## 1. Headline: all six canonicality levers together

**Final, all 3 shards** (298 / 299 molecules encoded):

| arm | byte-stable | comparison **key** broken |
|---|---|---|
| all levers OFF | 173/298 — **58.1%** | 60 — **20.1%** |
| all levers ON | 208/299 — **69.6%** | 16 — **5.4%** |
| **delta** | **+11.5 pts (+35 molecules)** | **60 → 16, a 73% reduction** |

Drift by subclass:

| subclass | OFF | ON | |
|---|---:|---:|---|
| `rdkit_canonical` | 91 | **18** | −73; the perception and body levers do the heavy lifting |
| `slot_renumber` | 42 | 74 | +32 — **reclassification, not regression** (see below) |
| `encode_fail` | 1 | 0 | |

Drift by transform: `renumber` 176 → 107, `both` 165 → 119, and **`rotate` is 0 in both arms** —
orientation-invariance was already sound and is preserved.

> **Correction to an earlier draft of this document.** I first published these figures from 2 of 3
> shards, as 62.0% → 73.5% and key-broken 39 → 8 (−79%). The third shard contained harder
> molecules, so the absolute levels are lower and the key reduction is −73% rather than −79%. The
> **delta held at exactly +11.5 points** across both reads, which is the quantity the promotion
> decision rests on — but the partial absolutes should not stand as if final, so they are replaced
> here rather than left in place.

**The key-instability number is the one that matters most.** The comparison key is the harness's
acceptance predicate and the basis of every accuracy figure this project reports. It goes from
unstable on **1 molecule in 5** to **1 in 19**.  
*(Corrected: an earlier revision said "1 in 25". 16 unstable of 298 is 1 in 18.6, so ~1 in 19. `docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md` already carried the correct figure.)*

### Why `slot_renumber` rising is not a regression

Lane 2 established this independently, with per-molecule accounting: **0 molecules broken** in
either of its arms. A molecule that previously drifted in *both* its ligand body and its slot
numbers is counted under `rdkit_canonical` (the first-matching subclass); once the body is
canonical, its remaining slot drift reclassifies it into `slot_renumber`. Byte-stability rising
while `slot_renumber` rises is exactly the expected signature.

The residual 45 `slot_renumber` and 8 key-broken molecules are the structurally harder class Lane 2
characterized and declined to force: **32/32 of its residual pairs are `same_vcolor_identical`**, so
no relabeling at that seam can close them (`docs/agentic-notes/v0.4.5/CANONICAL_SLOTS_v0.4.5.md` §7a), plus the 7
wrong-donor molecules now owned by Lane 9.

---

## 2. Veto checks — all passed

| veto | instrument | result |
|---|---|---|
| **`facmer_divergent` must not rise** (over-folding — the standing risk) | `test_facmer_key.py` + `test_isomer_divergence.py`, re-checked directly with `compare.py` | **PASS.** fac≠mer and cis/trans stay distinct **raw and at key level** with all six levers ON. The one non-green run (`OIN_CANONICAL_SLOTS` alone) is a **stale hardcoded golden string** — the lever deliberately relabels slots — not an isomer merge; verified directly. |
| Levers-OFF byte-identical | `test_regression_stability.py` (goldens) | **PASS** |
| **`OIN_STABLE_METAL_AC` geometry-tag shift** | `tools/geometry_tag_shift.py --n 300` | **PASS, 298 molecules: 0 string changes, 0 `[M_XXX]` changes, 0 coordination-number changes, no transitions.** |
| **`OIN_STABLE_STEREO` must not be stable-because-constant** | `test_stable_stereo_mirror.py`, re-run under all six levers ON | **PASS.** 10/10 mirrors differ; nothing collapsed. |
| Integrated suite | `discover tests/unit` | 729 tests, 717 OK, 3 skip, 4 xfail; the 5 errors were **one** missing fixture, since fixed. |

The geometry veto deserves emphasis because it was set *against* my own fix. Capping the metal
first can only **add** metal bonds the old atom-order-dependent iteration discarded, so coordination
numbers could have risen and the template fit could have reclassified polyhedra corpus-wide — and
`[M_XXX]` is not cosmetic: it selects the vertex table, hence the rotation group, hence the
canonical slot labelling and the key's entire vertex signature. **0/298 refutes that concern.**

---

## 3. Recommendation

**Promote all six to default-ON**, using the `OIN_EARLY_EXIT` template
(`metallogen_adapter.py:1636-1653`): membership test on `ff_params` so `False` can opt out, default
string `"1"` with `!= "0"`, and an in-code comment naming this document and its measured deltas.

| lever | lane | recommend |
|---|---|---|
| `OIN_CANONICAL_BODY` | 1 | **ON** |
| `OIN_CANONICAL_PERCEPTION` | 1 | **ON** |
| `OIN_CANONICAL_SLOTS` | 2 | **ON** |
| `OIN_CANONICAL_ETA_WINDING` | 3 | **ON** |
| `OIN_STABLE_METAL_AC` | 2 | **ON** |
| `OIN_STABLE_STEREO` | 8 | **ON** |

**Keep OFF** — these emit *new descriptors* the generator must reproduce, which converts a silent
false positive into a loud false negative. That is the right direction but it is a separate product
call, already made: injectivity levers stay opt-in for v0.4.5.

| lever | lane | recommend |
|---|---|---|
| `OIN_EMIT_AXIAL` | 4 | **OFF** (staged for v0.4.6, evidence recorded) |
| `OIN_EMIT_BOUND_AMINE` | 6 | **OFF** |
| `OIN_RESCUE_STUCK_RING` | encode_fail | **OFF** — its one molecule encodes but is not renumbering-stable |

The distinction that makes this coherent: the six promoted levers **repair renumbered presentations
without rewriting the canonical answer** — which is exactly why the corpus shows no churn
(`geometry_tag_shift` 0/298; goldens byte-identical). The three held levers **add information** to
the string. Determinism first, new descriptors second.

---

## 4. What this does and does not buy

It does **not** move the round-trip headline much, and that is structural rather than
disappointing: a canonicality defect lands in `key_equal`, which **already counts as passing**. Of
the 332 closeable molecules in the gap to the ~98.45% ceiling, **`hard_fail` is 306 — 92%** — so
generator throughput, not the notation, is the dominant accuracy lever.

What it buys is that the number becomes **meaningful**. Before these levers, 13% of molecules
encoded a different absolute stereochemistry depending on the input file's atom order, and the key
moved for 1 in 5 — so "round-trip success" was partly a property of how the XYZ happened to be
numbered. That is now 1 in 19 (16 of 298).

---

## 5. Reproduce

```bash
cd /home/tjmustard/Documents/GitHub/oin-v045-trial
export PYTHONPATH=$PWD/src
V=/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python
DS=/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset

# all-OFF arm
env -u OIN_CANONICAL_BODY -u OIN_CANONICAL_PERCEPTION -u OIN_CANONICAL_SLOTS \
    -u OIN_CANONICAL_ETA_WINDING -u OIN_STABLE_METAL_AC -u OIN_STABLE_STEREO \
  $V tools/canonicality_probe.py --dataset "$DS" --n 300 --trials 2 --shard 1:3 --out <dir>

# all-ON arm
env OIN_CANONICAL_BODY=1 OIN_CANONICAL_PERCEPTION=1 OIN_CANONICAL_SLOTS=1 \
    OIN_CANONICAL_ETA_WINDING=1 OIN_STABLE_METAL_AC=1 OIN_STABLE_STEREO=1 \
  $V tools/canonicality_probe.py --dataset "$DS" --n 300 --trials 2 --shard 1:3 --out <dir>

# the geometry veto
$V tools/geometry_tag_shift.py --lever OIN_STABLE_METAL_AC --n 300 --dataset "$DS" --out <dir>
```
