# Lane — `encode_fail` (the encoder returns nothing at all)

**The failure population:** the **48 molecules** of the v0.4.2 capstone whose crystal XYZ produced
**no OIN string at all** — `smiles_1 is None`. 0.71% of the 6,719-molecule corpus, and 11.0% of the
436-molecule gap-to-100% on the capstone arm
(`tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042/bucket_report.json`). They can never round-trip,
so they are a hard ceiling on the accuracy goal by construction. Worked twice: **v0.4.4 SL5**
classified and bounded, **v0.4.5** re-measured the residual against `main`.

⚠ **The single largest sub-population — 34 boron clusters — was declared a permanent ceiling by
this lane and by SL5 before it. That was wrong.** See `LANE-boron-cage.md`; this document's §5-equivalent
is superseded and the discrepancy is flagged explicitly below.

---

## ELI5

To turn 3D coordinates into a string, the encoder first has to *guess the chemistry*: which atoms
are bonded, and which of those bonds are double. For most molecules that guess is fast. For a few
it is a combinatorial nightmare — the code was building the full multiplication table of every
atom's possible bonding options (millions of rows) just to sort them, or asking RDKit to enumerate
every way of shuffling double bonds around a big conjugated ring, which it does inside a C++ call
that Python cannot interrupt. Those molecules did not fail; they **hung**, sometimes forever. The
work here was to put hard ceilings on both steps so an oversized ligand either finishes or gives up
quickly, while every ordinary ligand produces byte-for-byte the same string it always did — and to
give the genuinely impossible cases a *named* error instead of a generic crash, so a known
limitation can be told apart from a bug.

## The work, visually

```
  crystal XYZ  ──►  XYZToSMILES().convert()  ──►  get_tmc_mol  ──►  perception
                                                                         │
   ┌─────────────────────────────────────────────────────────────────────┴──────────┐
   │            THREE INDEPENDENT SUPER-POLYNOMIAL STAGES, each revealed            │
   │            only after bounding the previous one                                │
   │                                                                                │
   │  ① AC2BO valence-order sort            utils/xyz2mol_local.py                  │
   │     materialised the FULL Cartesian product of per-atom valences just to        │
   │     sort candidate assignments — exponential in the number of multivalent       │
   │     atoms.                                                                     │
   │     ● BOUNDED: _VALENCE_COMBO_CAP = 500_000  (test pins >= 100_000)             │
   │       above the cap: skip the sort, iterate the LAZY product bounded to         │
   │       _VALENCE_FALLBACK_TRIES = 20_000  (env OIN_VALENCE_FALLBACK_TRIES)        │
   │       main loop early-returns on the first valid assignment                    │
   │     ✔ PROVABLY byte-identical: a >cap ligand HANGS on unbounded main, so no     │
   │       currently-encodable molecule can reach the fallback                       │
   │                                                                                │
   │  ② ResonanceMolSupplier                utils/xyz2mol.py::lig_checks             │
   │     builds the conjugation-electron GROUPS in a C++ call BEFORE any             │
   │     enumeration, so maxStructs cannot bound it (even maxStructs=2 hangs).       │
   │     ⚠ IT HOLDS THE GIL ⇒ a watchdog THREAD cannot interrupt it.                 │
   │     ● BOUNDED BY THE KERNEL, not by Python:                                    │
   │          _resonance_candidates_isolated()                                      │
   │            os.fork  +  pipe  +  select                                         │
   │            child: RLIMIT_CPU = _RESONANCE_CPU_BUDGET_S = 120 CPU-seconds        │
   │            parent: select timeout  _RESONANCE_WALL_SAFETY_S = 900 s backstop    │
   │                                                                                │
   │       finishes in budget ──► child returns forms (property-preserving          │
   │                              ToBinary), rebuilt in parent ⇒ BYTE-IDENTICAL     │
   │       burns past budget  ──► kernel SIGXCPU kills the child (enforced even      │
   │                              mid-C-call) ⇒ degrade to the SINGLE perceived      │
   │                              form ⇒ the molecule RECOVERS                      │
   │       gate: _resonance_needs_isolation — >=50 heavy OR >=35 aromatic atoms      │
   │                                                                                │
   │  ③ get_UA_pairs → networkx max_weight_matching   xyz2mol_local.py               │
   │     O(V³) in the unsaturated-atom graph.                                        │
   │     ✘ NOT BOUNDED.  Documented residual (FAQYUU, HICLAG).                       │
   └────────────────────────────────────────────────────────────────────────────────┘
                                          │
                                          ▼
   ┌──────────────── the 48, re-triaged against current main (200 s/mol) ───────────┐
   │  boron_cluster          34  ← declared PERMANENT here.  ✘ WRONG — see the      │
   │                                boron lane: 0/36 → 34/36 encoding              │
   │  resonance_timeout       7   BENVOG FAQYUU HICLAG HOHKUL KEMTED KESWUB NAKLET   │
   │  encodes_now             3   HOCVAY HUCNAU WEFZAL — need NOTHING from this lane │
   │  other (quinoid/ylide)   3   KAXVOX KAXWAK LEZWAO — typed error, deferred       │
   │  perception_charge_gap   1   ASISAX — fixed, OIN_RESCUE_STUCK_RING (default OFF)│
   └────────────────────────────────────────────────────────────────────────────────┘

  TYPED TERMINUS (all buckets):  aromaticity.py::OINEncodeError(ValueError)
     names the specific limitation · core/translator.py::convert re-raises rather
     than flattening · subclasses ValueError so `except ValueError` still catches it
```

Legend — `●` = the bound that was added. `✔` = a byte-identity property that is *provable*, not
sampled. `✘` = not fixed / refuted. `①②③` = the three hang stages, in the order they were
discovered (the handoff hypothesis named only ②).

## Initial assumptions and hypothesis

1. **"The 34 boron clusters are a permanent representational ceiling."** Inherited from SL5 W1 and
   restated by this lane as *"Confirmed unfixable"*: RDKit's 2c-2e valence model has no Lewis
   structure for a 3c-2e cage, the charge sweep already spans −4..+4 (widened diagnostically to
   −6..+6 with nothing new), so it needs multi-centre bonding — out of scope for a valence-graph
   encoder. **This is the lane's biggest wrong answer.** It is refuted in full by
   `LANE-boron-cage.md`: the cages were destroyed by a valence *rule* in adjacency perception, two
   stages before any bond-order search, and `OIN_BORON_CAGE` now encodes them.
2. **"The timeout cohort hangs in `ResonanceMolSupplier`."** The handoff named one stage. There are
   **three**, and each was only visible after the previous one was bounded.
3. **"A wall-clock timeout is the obvious way to bound a hang."** It is the wrong instrument here —
   see the dead ends.
4. **"A size-based skip is the cheap version of a time bound."** Tried first, withdrawn.
5. **"The reported failure site is the causal site."** False for the boron cohort (reported
   `get_lig_mol`, caused in `xyz2AC_obabel`) and it is the general reason this cohort was
   mis-triaged.
6. **"`tests/unit/test_xyz2mol_errors.py` tests the error path."** It had stopped doing so
   entirely — see below. Nobody had assumed to check.

## What was actually found

### Confirmed — the 48 split into five buckets, and 3 of them need no fix

Full 48/48 re-triage against current `main`-equivalent state (200 s/mol wall-clock cap,
generation-free — `XYZToSMILES().convert()` in an isolated subprocess per molecule, so every result
is unambiguously encoder-side):

| bucket | count | molecules |
|---|---:|---|
| `boron_cluster` | 34 | ⚠ **not a ceiling** — see the boron lane |
| `resonance_timeout` | 7 | `BENVOG`, `FAQYUU`, `HICLAG`, `HOHKUL`, `KEMTED`, `KESWUB`, `NAKLET` |
| `encodes_now` | 3 | `HOCVAY`, `HUCNAU`, `WEFZAL` — **already succeed on unmodified `main`** |
| `other` (quinoid/ylide kekulize) | 3 | `KAXVOX`, `KAXWAK`, `LEZWAO` |
| `perception_charge_gap` | 1 | `ASISAX` — fixed this lane |
| **total** | **48** | |

SL5's own earlier triage of the identical 48 was: `boron_cluster` 34, `resonance_timeout` **10**,
aromatic/quinoid 3, `perception_charge_gap` 1. The 10 → 7 movement is the `encodes_now` 3.

### Confirmed — the frozen 48-count overstates the *encoder's* population

**3 of the 48 already encode with no fix from this lane at all.** And only one of the three is an
encoder story:

- **`HUCNAU`** — genuinely recovered by SL5's own forked, CPU-time-bounded resonance. This
  confirms that fix works on real data, not only in its own artificially-shortened-budget unit
  test.
- **`HOCVAY`, `WEFZAL`** — **not an encoder story at all.** The frozen
  `results-capstone-v042/bucket_report.json` records `HICLAG`, `HOCVAY` *and* `WEFZAL` failing with
  `Generation/Verification failed at UFF_1: child process died with exit code -9` — an **OOM kill**.
  They land in `encode_fail` only because the harness's `_ENCODED` marker mechanism does not survive
  a child process dying: a report with `smiles_1: null` cannot distinguish *"encoder hung"* from
  *"encoder finished fine, generator died later"*. Tested directly against
  `XYZToSMILES().convert()` alone: `HOCVAY` encodes immediately; `WEFZAL` encodes in **1.9 s** —
  `[Pd_SPL].CC(=O)OC[C@H]1O[C@@H](N2C{0}N(CCCN3C{1}N([C@@H]4O[C…` (a Pd-nucleoside/PNA-type
  complex).
- `HICLAG`, the third member of that OOM-labelled trio, genuinely does **not** resolve within 200 s
  even encoder-only — consistent with SL5 naming it (alongside `FAQYUU`) as the unrecovered
  `get_UA_pairs` / `max_weight_matching` O(V³) hang. A real encoder-side cost, not a harness
  artifact.

So **2 of 48 are not `encode_fail` at all** under a clean generation-free signal. This is a real
caveat on the 11.0% / 48 headline: it measures the *old harness's bucketing*, not a clean
encoder-only signal, and some further fraction of the frozen `encode_fail` rows may belong to
generation-side OOM. Flagged for whoever owns `hard_fail`/generation; a full re-triage would need
re-running encode+generate per molecule and was out of this lane's scope.

### Confirmed — `ASISAX` is a one-line over-conservatism, and the fix is honest about what it buys

`ASISAX` is a Ni tetraaza-macrocycle. Direct diagnosis: at ligand charge **0** the encoder's own
`AC2mol` perceives it fine; **every other charge in −6..+6 returns `None`**.
`get_lig_mol`'s charge-rescue loop (`_rescue_unusable_perception`, `xyz2mol.py:498`) rejected the
charge-0 candidate anyway, because it tested
`stuck_ring_atoms(candidate) or not _perception_is_usable(candidate)` — treating *any* stuck
(unkekulizable-as-aromatic) ring as an automatic reject, even though `_perception_is_usable`
**already** calls `kekulize_safe_sanitize`, which can repair a stuck ring by de-aromatizing it. For
`ASISAX` that repair succeeds. **The rescue loop was strictly more conservative than the encoder's
own downstream repair path, for no documented reason.** With every charge rejected, the sweep
exhausts and `get_lig_mol` raises a generic `ValueError` → `encode_fail`.

With `OIN_RESCUE_STUCK_RING=1`, `ASISAX` encodes, deterministically across a repeat:

```
[Ni_SPL].CC(C)C1=c2c3c4c(c5ccccc25)=C(C(C)C)C(=N{0}4)C=CC2=N{1}c4c5c(c6ccccc6c4=C2C(C)C)
=C(C(C)C)C(=N{2}5)C=CC1=N{3}3
```

**But the canonicality check is NOT clean, and it was reported that way rather than as a win.**
`tools/canonicality_probe.py --only ASISAX_comp_0 --trials 2` with the lever on:

```
ASISAX_comp_0: 4/6 drifted ['rdkit_canonical']
drift by transform: renumber 2, both 2   (rotate: 0)
of the drifted, key ALSO changed (isomer-level, worse): 1
```

So `ASISAX` encodes deterministically for its *actual* atom ordering (the one real dataset file has
one order; repeat runs agree; pure rotation never drifts) but is **not robust to renumbering** —
presenting the same graph with atoms reordered changes the OIN string, including at the
comparison-key level. **That is "moved to a different bucket", not "fixed."**

This is not a new defect: it is the already-tracked corpus-wide renumbering instability
(`docs/RENUMBERING_INSTABILITY_v0.4.5.md` — 43.9% of a 225-molecule sample drift under pure
renumbering, 21.1% at key level), whose suspect #1 is exactly *"`AC2BO` / resonance-form
order-dependence changing perceived bond orders"*, which is this macrocycle's failure mode.
`ASISAX` could not *show* the defect before, because it never reached a successful encode; the fix's
only contribution is that it now can.

### Confirmed — the encoder was also emitting UNPARSEABLE strings, and the same molecules were 100× too slow

This is adjacent to the lane rather than inside it, and it is worth a subsection because it changes
what "the encoder succeeded" means. Two metal–NHC complexes — `QIDKUL_comp_0` (Rh(COD)Cl) and
`QIDKIZ_comp_0` (AuCl), both of a bis(2,4-dinitrophenyl)imidazol-2-ylidene — were emitting a string
that **does not read back**, and paying two minutes of CPU to do it.
Measured (`docs/VALENCE_ORDER_v0.4.5.md`):

| | default | `OIN_VALENCE_CHARGE_FILTER=1` |
|---|---|---|
| `QIDKUL_comp_0` whole encode | **124.44 s**, `found_valid=0` on all 3 over-cap calls | **0.87 s**, `found_valid` on all of them |
| emitted OIN | `…-n2c{0}n(…)cc2…` (aromatic imidazolium) | `…N2C{0}N(…)C=C2…` (neutral carbene) |
| does the string re-read? | **3 of 4 fragments — the ligand fails `KekulizeException`** | **4 of 4** |
| `QIDKIZ_comp_0` whole encode | **129.53 s**, `found_valid=0` | **0.32 s** |
| does the string re-read? | **2 of 3 fragments — same `KekulizeException`** | **3 of 3** |

The shipped default emits the carbene as a lowercase aromatic imidazolium that RDKit cannot
kekulize; the filter emits it as the neutral carbene ring, which does re-read. **This is not "a
different string" — it is the difference between a string that parses and one that does not.**

⚠ **Attribution correction.** These numbers are often paraphrased as "one case went from 124 s to
0.87 s once valence search was bounded." That is not quite what happened, and the distinction
matters for anyone reasoning about the cap. The **bound** (`_VALENCE_COMBO_CAP` /
`_VALENCE_FALLBACK_TRIES`) is what makes an over-cap ligand *terminate at all* instead of hanging.
The **124.44 s → 0.87 s** improvement comes from `OIN_VALENCE_CHARGE_FILTER` (default **OFF**,
`swimlane/v045-valorder`), which enumerates only the provably-feasible valence assignments: for
QIDKUL's 37-atom ligand exactly **16 of 1,259,712** candidates can satisfy `AC2BO`'s own acceptance
predicate, and the earliest sits at **rank 209,858** — so no prefix of 20,000, *in any order*, could
ever have reached one. The bound stops the hang; the filter is what finds the answer.

### Confirmed — `test_xyz2mol_errors` was SILENTLY TESTING NOTHING

The lane's cleanest process finding, and it has nothing to do with the 48.

The original test drove `get_tmc_mol`'s error path with a captured TiCat3 *generated* structure
(`tests/fixtures/ticat3_generated_broken.xyz`) whose clashing geometry over-connects a ligand
fragment. Promoting **`OIN_STABLE_METAL_AC`** to default-ON changed that. That lever caps valences
**highest-Z-first**, so perception no longer depends on input atom order — and on this degenerate
geometry it lets the **titanium absorb the contested bonds** instead of the walk dead-ending.
Perception then **succeeds**, and returns nonsense: **48 atoms in 8 fragments, seven bare `[H+]`
ions, and a `[Ti-14]` centre.** The error-path contract test therefore no longer fired at all.

Rewritten as **fault injection**:

```python
# (None, charge): get_lig_mol returns a 2-TUPLE that the call site unpacks BEFORE the
# `if not lig_mol` guard, so injecting a bare None would fail in the unpack -- which is
# the very TypeError this contract exists to prevent, raised from the wrong place.
with mock.patch.object(xyz2mol_module, "get_lig_mol", return_value=(None, 0)):
    with self.assertRaises(ValueError) as ctx:
        get_tmc_mol(_FIXTURE, 0, with_stereo=False)
self.assertIn("get_lig_mol failed", str(ctx.exception))
```

That parenthetical is the whole lesson in miniature: the *shape* of the injected fault has to match
the call site, or the test passes for the wrong reason. And a second test now **pins the degenerate
behaviour** rather than pretending it is fine
(`test_broken_fixture_perceives_a_degenerate_graph_under_stable_metal_ac`: >1 fragment, >0 bare
protons), so if a future sanity gate makes this raise again, the docstring gets revisited instead of
the behaviour drifting unobserved.

**On real data the lever is clean:** capstone A/B **145 molecules fixed, zero correctness
regressions**, and `tools/geometry_tag_shift.py` shows **0 of 298 `[M_XXX]` changes**. So this is a
degenerate-input concern, not a shipped-accuracy one. The follow-up worth considering — a sanity
gate rejecting a perceived molecule that contains isolated bare-proton fragments — was deliberately
**not** added, because charged hydrides are legitimate and the gate needs its own corpus A/B first.

### Not resolved, and honestly labelled

The remaining **7** `resonance_timeout` molecules did not resolve within 200 wall-clock seconds. The
caveat on that number is recorded in full because it changes what the number means: the box was
running **five sibling v0.4.5 lane agents concurrently**, load **8–12 on 12 cores** throughout
(`uptime`, plus concurrent `canonicality_probe.py` and lane test-suite jobs visible in `ps aux`).
The fork's bound is CPU-*time* specifically so the *outcome* is load-independent — but under
contention the same 120 CPU-seconds can take substantially more than 120 *wall* seconds to accrue.
So a `resonance_timeout` result here means *"did not resolve within 200 wall-clock seconds on a
heavily loaded shared box"*, **not** "the fork budget is broken."

`faulthandler`-based profiling of `BENVOG` (`tools/sl5_profile_hang.py`) confirmed it is genuinely
inside `_resonance_candidates_isolated` (blocked in the parent's `select()`) at 150 s — the
mechanism is doing what it was designed to do, just not within the reproduction's budget. Whether
these 7 would resolve on an isolated machine within the existing 900 s `_RESONANCE_WALL_SAFETY_S`
backstop is **not measured**.

## What was done

### `OIN_RESCUE_STUCK_RING` — one-line permissive rescue, opt-in, default OFF

`utils/xyz2mol.py::_rescue_unusable_perception`. Split the combined boolean so the stuck-ring
rejection is skipped when the lever is set, leaving `_perception_is_usable` (which already subsumes
ring repair) as the sole gate. Lever unset ⇒ byte-identical to the pre-fix code, by the same
short-circuit.

**Why gated rather than landed unconditionally**, recorded because the reasoning is not obvious: the
early return `if best_res_mol is not None and _perception_is_usable(best_res_mol): return` means the
loop is *only* reached by ligands whose default perception was already unusable — so far, so "only
touches things that currently fail". **But** the loop returns the **first** charge in
Hückel-distance order that passes its checks, so loosening the stuck-ring test can make an *earlier*
charge win where today a *later* one does, for any ligand that currently reaches this loop and is
rescued by some non-stuck charge further down the order. That could not be ruled out without a full
corpus sweep, so it went behind a default-OFF lever.

### `_VALENCE_COMBO_CAP` — the AC2BO valence-order bound (SL5, default behaviour)

`utils/xyz2mol_local.py`. `_VALENCE_COMBO_CAP = 500_000`; above it, skip the sort and iterate the
lazy product bounded to `_VALENCE_FALLBACK_TRIES = 20_000` (overridable with
`OIN_VALENCE_FALLBACK_TRIES=<int>`); the main loop early-returns on the first valid assignment.

**The cap must stay far above ordinary-ligand scale.** A cisplatin/ferrocene-scale ligand's valence
product is a few dozen at most, so pinning it is what stops a future edit silently lowering it into
the range that would perturb real ligands. The guard asserts a floor, not the value:
`tests/unit/test_encoder_robustness.py::TestAc2boCapIsByteIdentical::test_cap_is_large` asserts
`_VALENCE_COMBO_CAP >= 100_000`. **Byte-identity here is provable rather than sampled:** a
>cap-combo ligand *hangs* on unbounded `main`, so no currently-encodable molecule can reach the
fallback.

### Forked, CPU-time-bounded `ResonanceMolSupplier` (SL5, default behaviour)

`utils/xyz2mol.py::_resonance_candidates_isolated`. `os.fork` + pipe + `select`; child sets
`RLIMIT_CPU = _RESONANCE_CPU_BUDGET_S = 120` CPU-seconds; parent's `select` has a
`_RESONANCE_WALL_SAFETY_S = 900` wall-clock backstop for a starved child that never even runs. Only
large ligands take this path — `_resonance_needs_isolation`: `_RESONANCE_ISOLATION_HEAVY = 50` heavy
atoms **or** `_RESONANCE_ISOLATION_AROMATIC = 35` aromatic atoms. Ordinary ligands run inline,
unchanged.

Three design choices, each with a reason that is easy to lose:

- **The GIL forces the process boundary.** `ResonanceMolSupplier` builds the conjugation-electron
  groups in a C++ call *before* any enumeration, so `maxStructs` does not bound it (even
  `maxStructs=2` hangs), and it holds the GIL inconsistently — a watchdog **thread** cannot reliably
  interrupt it. Only the kernel can: `SIGXCPU` is enforced even mid-C-call.
- **CPU-time, not wall-clock.** A wall-clock timeout would make the *outcome* depend on machine
  load — a starved *completer* could wrongly fall back and change its OIN. A given ligand's
  resonance burns roughly constant CPU seconds whether the box is idle or saturated, so CPU-time
  keeps the encode deterministic. A higher budget only lets *more* completers finish (hangs burn CPU
  without bound and are always killed), so the budget is set generously above every observed
  completer.
- **`os.fork`, not `multiprocessing`.** It works inside the dataset harness's *daemon* workers
  (which forbid child processes) and inherits the ligand copy-on-write, so there is no input
  pickling. Forms come back as property-preserving `ToBinary` and are reconstructed in the parent —
  **byte-identical to the inline path**, verified on e.g. `EHADAV` (Co macrocycle, ~21 CPU-s) and
  `RUTJEW` (72 aromatic atoms), both reproducing `main`'s OIN exactly.

`BENVOG_comp_0` is the fixture: its macrocycle burns CPU in `ResonanceMolSupplier` without ever
finishing.

### `OINEncodeError` — the typed terminus

`utils/aromaticity.py:24`, `class OINEncodeError(ValueError)`. It **subclasses `ValueError` for
back-compat**, so existing `except ValueError` handlers (and `core/translator.py`, which re-raises)
keep working — and `core/translator.py::convert` re-raises the typed error rather than flattening
it, so callers can distinguish a known ceiling from an unexpected failure. Raised by
`_is_electron_deficient_cluster` (boron), by `kekulize_safe_sanitize` (`aromaticity.py:161`) for the
quinoid/ylide cases, and named per-limitation rather than generically.

### Fixtures added, so the next person does not re-pull from the dataset

`tests/fixtures/{ASISAX,KAXVOX,KAXWAK,LEZWAO}_comp_0.xyz`.

## Dead ends and refutations

### "The 34 boron clusters are permanently unfixable" — REFUTED (and it was this lane's own conclusion)

**Killed by `swimlane/v045-boron`.** The cages were perceived *correctly* by `xyz2AC_obabel`'s
distance criterion (993 B–B edges, textbook-exact topologies) and then destroyed by a pruning loop
that deletes an atom's longest bonds while its connectivity exceeds `max(atomic_valence[Z])` — 4 for
boron, against a cage vertex's 5–6. **406 of 993 B–B bonds deleted, at 1.712–2.105 Å.** The 2c-2e
vs 3c-2e argument is individually true and was never exercised by the failure attributed to it.
`OIN_BORON_CAGE` takes the population from **0/36 to 34/36** encoding on the 936-molecule
re-baseline and was promoted to default-ON in v0.4.6 (`d799de1f`).

**Why the wrong answer was believable:** the *reported* failure site (`get_lig_mol` →
`OINEncodeError`, `xyz2mol.py:775`) and the *causal* site (`xyz2AC_obabel`'s pruning loop, which
raises nothing) are different stages. Confirming that no charge in −6..+6 helps is a true
measurement — of the wrong stage. **A confirmed negative on the wrong stage is how a misdiagnosis
acquires evidence.**

### "Widening the charge sweep will fix it" — CONFIRMED, and irrelevant

−4..+4, widened diagnostically to −6..+6, finds nothing new. True. Also measuring debris.

### "A wall-clock timeout is the right way to bound the resonance hang" — REFUTED

**Killed by the determinism requirement, not by a benchmark:** wall-clock makes the *outcome*
load-dependent, so a starved completer would fall back and its OIN would change. Encoder output must
not be a function of machine load. Hence `RLIMIT_CPU`.

### "A watchdog thread can bound it" — REFUTED

**Killed by the GIL.** `ResonanceMolSupplier` holds it inconsistently across a C++ call, so a Python
thread cannot reliably interrupt. This is why the design pays for a `fork`.

### "A size-based skip is the cheap version" — TRIED, WITHDRAWN, with a named regression

Skipping resonance for large ligands **changed the perceived form of passing molecules whose
resonance completes quickly**. `EHADAV` regressed: resonance finds the aromatic form, the skip kept
the localized tautomer. The forked bound keeps a completer's real resonance result and falls back
only on a *genuine* hang. Recorded so nobody re-proposes the skip as a simplification.

### "`sl5_triage.py`'s stock timeout is a valid signal" — REFUTED, and it is a live trap

The stock tool's `PER_MOL_TIMEOUT_S = 90` is **shorter than SL5's own
`_RESONANCE_CPU_BUDGET_S = 120`** CPU-second fork budget, so it would misreport as a hang a molecule
the fork *would* recover. Re-run with a 200 s wall-clock cap instead. Anyone re-running this triage
must raise the cap or they will reproduce the wrong histogram.

### "`test_xyz2mol_errors` tests the error path" — REFUTED

**Killed by reading what the fixture now does under `OIN_STABLE_METAL_AC`:** perception *succeeds*
on the deliberately-broken input, returning 48 atoms in 8 fragments with seven bare `[H+]` and a
`[Ti-14]`. The contract test had stopped firing. Generalisable lesson, in the module's own words: *a
contract test that depends on a fixture staying unbuildable is one perception improvement away from
silently testing nothing* — and here it would have hidden a regression back to the bare-`None`
`TypeError` the contract exists to prevent.

### "Bounding the valence search is what made the NHC cases fast" — IMPRECISE

The bound makes them *terminate*; `OIN_VALENCE_CHARGE_FILTER` is what makes them *fast and correct*
(124.44 s → 0.87 s, 3/4 → 4/4 fragments re-reading). See the attribution correction above. Cited
because the two are routinely conflated and the conflation makes the cap look like an accuracy fix,
which it is not.

### The 3 quinoid/ylide cases — DEFERRED, with the reason, not silently dropped

`KAXVOX`, `KAXWAK` (the same Zn/N5O2S2 core, `KAXWAK` brominated — clearly two derivatives of one
compound) and `LEZWAO` (an unrelated Pd/phosphine/silane complex) already raise a typed
`OINEncodeError` from `kekulize_safe_sanitize`: after de-aromatizing the detected quinoid ring the
residual bond orders are still unusable (`KAXVOX`/`KAXWAK` still cannot kekulize 5 atoms;
`LEZWAO` leaves a carbon at explicit valence 5). The mechanism is already named in that function's
docstring — *"at the wrong ligand charge `AC2mol` leaves the ipso carbon pentavalent, and only
re-perceiving the charge fixes that"* — but it lives in a **different code path** from the
`ASISAX` fix (`fix_equivalent_Os`'s whole-molecule equivalent-oxygen pass, not the per-ligand-fragment
charge rescue). Fixing it means teaching that pass to re-perceive charge the way
`_rescue_unusable_perception` already does: a materially larger change for 3 molecules. **"Needs new
chemistry", not "one permissive valence entry."**

## Where it landed

| change | lever | default | where |
|---|---|---|---|
| stuck-ring rescue permissiveness | `OIN_RESCUE_STUCK_RING` | **OFF** (held) | `utils/xyz2mol.py::_rescue_unusable_perception` |
| AC2BO valence-combo bound | none (SL5) | **ON** | `utils/xyz2mol_local.py`, `_VALENCE_COMBO_CAP` / `_VALENCE_FALLBACK_TRIES` |
| forked CPU-bounded resonance | none (SL5) | **ON** | `utils/xyz2mol.py::_resonance_candidates_isolated` |
| typed encode ceiling | none (SL5) | **ON** | `utils/aromaticity.py::OINEncodeError` |
| boron cage encoding | `OIN_BORON_CAGE` | **ON since v0.4.6** | see `LANE-boron-cage.md` |
| valence charge filter (the NHC fix) | `OIN_VALENCE_CHARGE_FILTER` | **OFF** | `utils/xyz2mol_local.py` (`swimlane/v045-valorder`) |
| `get_UA_pairs` / `max_weight_matching` O(V³) | — | **unbounded** | documented residual |

`OIN_RESCUE_STUCK_RING`'s held-off reason is recorded verbatim in
`src/oinsmiles/oin/levers.py::_HELD_OFF`: *"its one molecule (ASISAX) encodes but is not
renumbering-stable, so promoting moves it between buckets rather than fixing it."*

**Lane commits** — `swimlane/v045-encodefail`, merged into `release/v0.4.5` at **`2579bfbb`**, then
into the integration merge **`1450b5ce`** and the v0.4.5 release commit **`0d165845`**:

| commit | what |
|---|---|
| `b9143ac1` | `WIP(encodefail): encode_fail triage + ASISAX rescue -- histogram tally PENDING` — the substantive commit: `xyz2mol.py` lever split, 4 fixtures, `test_encoder_robustness.py` +48 lines. ⚠ Its message carries the "34 boron are PERMANENTLY unfixable" claim, now refuted. |
| `e7d58ca8` | `docs(v0.4.5): finalize the encode_fail histogram -- 14/48 addressable, 34/48 a correct ceiling` — the final 48/48 tally, the `encodes_now` 3, the OOM-misbucketing finding, the `BENVOG` profile |

SL5's own work (`_VALENCE_COMBO_CAP`, the forked resonance, `OINEncodeError`) landed on `main` in
v0.4.4, **before** this branch forked.

**Guard tests** — `tests/unit/test_encoder_robustness.py`:

| class / test | pins |
|---|---|
| `TestBoronClusterTypedCeiling` | ⚠ **now pins the OPT-OUT contract, not shipped behaviour.** `setUp` sets `OIN_BORON_CAGE=0` explicitly; asserts a typed `OINEncodeError` naming `"boron cluster"` from `get_tmc_mol` *and* propagated through `convert`. It previously set nothing and relied on the default — which stopped meaning "off" at promotion, at which point both tests were asserting a ceiling the shipped encoder no longer has. `RAWJEG_comp_0.xyz` now **encodes** in the default configuration. |
| `TestBoronClusterTypedCeiling::test_typed_error_is_a_value_error` | back-compat: `issubclass(OINEncodeError, ValueError)` |
| `TestElectronDeficientClusterDetector` | `_is_electron_deficient_cluster` fires on a B–B–B chain, not on `BPh4-`, not on scattered borons with no B–B bond, not on pyridine |
| `TestAc2boCapIsByteIdentical::test_cap_is_large` | `_VALENCE_COMBO_CAP >= 100_000` |
| `TestAc2boCapIsByteIdentical::test_ordered_valences_matches_unsorted_content` | the extracted sorter is a **permutation** of the raw product — same members, so the sort cannot change *which* candidates exist |
| `TestForkedResonanceRecovery::test_benvog_recovers_via_cpu_budget_fallback` | monkeypatches `_RESONANCE_CPU_BUDGET_S = 2` so `SIGXCPU` fires fast; `BENVOG` encodes via the fallback, leads with `[Ni`, and is deterministic across a repeat |
| `TestStuckRingRescuePermissive` | lever unset ⇒ `ASISAX` still raises `ValueError`; lever on ⇒ encodes, leads with `[Ni`, deterministic |
| `tests/unit/test_xyz2mol_errors.py::test_get_tmc_mol_raises_valueerror_when_get_lig_mol_fails` | the fault-injected error-path contract, with the 2-tuple caveat in a comment |
| `tests/unit/test_xyz2mol_errors.py::test_broken_fixture_perceives_a_degenerate_graph_under_stable_metal_ac` | pins the degenerate perception (>1 fragment, >0 bare protons) so it cannot drift unobserved |

Assertions are deliberately rdkit-version robust — a non-empty encode led by the correct metal
token plus forward-encode stability, no exact bond-direction strings.

**Suite and lint at the lane's close:** full `tests/unit` discovery **607 tests, OK (skipped=3,
expected failures=3)** — no new failures from the `OIN_RESCUE_STUCK_RING` change or the 4 added
fixtures; `test_regression_stability.py` (4 goldens) **6/6 pass, byte-identical**;
`test_encoder_robustness.py` **12/12**; `uvx ruff@0.15.20 check` and `format --check` clean on every
changed file.

### ⚠ Superseded text still in the tree

| location | stale claim |
|---|---|
| `docs/ENCODE_FAIL_v0.4.5.md` §5 | **"Confirmed unfixable: 34 boron clusters"**, *"needs a different bonding model entirely (multi-center bonds)"*, *"an encoder that refuses this input is correct, not a bug"*. **Refuted.** |
| `docs/ENCODE_FAIL_v0.4.5.md` §7 | the summary table's *"confirmed unfixable — 34"* row, and *"34 of 48 (70.8%) are an honest, correct ceiling"*. That row is now **0**; the correct total addressed is **48**, not 14 (34 loud + 14 silent, all round-tripping). This is the row a planner reads when deciding what to work on, which makes it the most costly stale sentence in the docs tree. |
| `docs/ENCODER_ROBUSTNESS_v0.4.4_SL5.md` §W1 and its "Net effect" bullet | *"34/48 encode-fails now fail with a typed, classified `OINEncodeError`"* as the terminal state. |
| `tests/unit/test_encoder_robustness.py` module docstring | **already corrected** — it now carries `⚠ NO LONGER A CEILING IN THE DEFAULT CONFIGURATION` and explains that W1 pins the opt-out contract. Use this as the template for the doc fixes above. |
| `src/oinsmiles/utils/xyz2mol.py::_is_electron_deficient_cluster` docstring | *"a permanent representational ceiling of the RDKit valence model"*. |
| `docs/ENCODE_FAIL_v0.4.5.md` §9 reproduce block | `tools/sl5_triage.py` still ships `PER_MOL_TIMEOUT_S = 90`, shorter than the 120 CPU-second fork budget. Running it as written reproduces a wrong histogram. |

## Open questions / for the next agent

1. **How much of the frozen `encode_fail` bucket is actually generation-side OOM?** Two of 48
   (`HOCVAY`, `WEFZAL`) were proven to be, because the harness's `_ENCODED` marker cannot survive a
   dying child. Nobody has re-triaged the rest of `bucket_report.json` with a clean encoder-only
   signal. Any future `encode_fail` number should be produced generation-free, per this lane's
   methodology, or it inherits the same ambiguity.
2. **Do the 7 `resonance_timeout` molecules resolve on an idle machine within the 900 s
   `_RESONANCE_WALL_SAFETY_S` backstop?** Not measured. The measurement requires *uncontended*
   machine time — and note the load-independence rule from the boron lane: outcomes are
   load-independent, wall-clock is not, so design the probe to report pass/fail rather than seconds.
3. **Bound stage ③.** `get_UA_pairs` → networkx `max_weight_matching` is O(V³) and unbounded;
   `FAQYUU` and `HICLAG` die there, not in resonance. It is a candidate for the same fork +
   `RLIMIT_CPU` treatment the resonance stage got, and the design is already written down — the only
   open question is whether a degraded matching yields a *usable* perception or just a different
   wrong answer.
4. **The 3 quinoid/ylide cases need charge re-perception inside `fix_equivalent_Os`.** Fixtures are
   already committed. The scope is "teach that whole-molecule pass to re-perceive charge the way
   `_rescue_unusable_perception` does for fragments".
5. **`OIN_RESCUE_STUCK_RING` cannot be promoted on `ASISAX` alone.** Two things gate it: a corpus
   sweep answering *"does any ligand that currently reaches `_rescue_unusable_perception` and is
   rescued by a later charge now land on an earlier one?"*, and `ASISAX` becoming
   renumbering-stable, which is owned by the canonical-perception lanes rather than here.
6. **`AC2BO`'s bare `sys.exit()` is still live** and is not boron-specific: `SystemExit` is a
   `BaseException`, so no `except Exception` in the encoder can contain it. Only cage fragments are
   currently routed around it.
7. **A sanity gate on perceived molecules containing isolated bare-proton fragments** is the
   documented follow-up to the `OIN_STABLE_METAL_AC` degenerate-input finding. Deliberately not
   added: charged hydrides are legitimate, so the gate needs its own corpus A/B first.
