# The success metric's false-positive rate — measured

**9.6 % of harness-scored successes do not survive independent re-perception with their
coordination intact. For haptic molecules it is 28.1 %.**

Measured 2026-07-26 on `results-v0.4.5-rebaseline` (633 scored successes, generator-free).
Tool: `tools/haptic_false_positive.py`. Data: `docs/agentic-notes/v0.4.6/metric_false_positives.json`.

---

## ELI5

The project grades itself by turning a 3D structure back into a string and checking the string
matches. But the grader asks the *generator* which atoms are bonded, instead of working it out
from the coordinates. So if the generator builds a ferrocene whose two rings have drifted almost
an ångström away from the iron, the generator still believes they are attached — and the grader
believes the generator. The molecule is scored as a perfect round trip.

We measured how often that happens. About 1 in 10 of all "successes", and **more than 1 in 4** of
the sandwich-type (haptic) molecules, are structures whose coordination has actually fallen apart.

This does not mean the notation is broken. It means the **ruler** is short, and every accuracy
number the project has reported is measured with it.

---

## The shape of the measurement

```
   report says status == "success"
              │
              ├─ re-encode the INPUT xyz        with the CURRENT encoder → oin_in
              └─ re-encode the STORED GENERATED xyz, same encoder        → oin_indep
                    (full XYZToSMILES().convert -- bonds perceived from
                     COORDINATES ALONE, blind to the generator's beliefs)
              │
              ▼
        compare, and classify by WHAT changed -- not by whether the string changed
              │
     ┌────────┴─────────────────────────────────────────────┐
     │ KEY_MATCH            562   independent re-perception agrees → genuine pass
     │ KEY_DIFF_COORD_OK     10   string drifted, coordination intact → NOT charged
     ├──────────────────────────────────────────────────────┤
     │ GEO_DEGRADED          54 ┐                            │
     │ HAPTIC_LOST            4 │ coordination NOT supported │ → FALSE POSITIVES
     │ HAPTICITY_REDUCED      2 │ by the geometry            │    61 / 633 = 9.6%
     │ DENTICITY_LOST         1 ┘                            │
     └──────────────────────────────────────────────────────┘
              (54 + 4 + 2 + 1 = 61; all seven classes sum to 633)
              │
     concentrated exactly where predicted:
       haptic inputs      48/171 = 28.1%
       non-haptic inputs  13/462 =  2.8%     ← a 10x concentration
```

---

## 1. Why this was worth measuring

`docs/agentic-notes/v0.4.7/ACCEPT_SCORED_v0.4.7.md` §4.7 established the mechanism: `tools/test_dataset_roundtrip.py`
scores a round trip with

```python
oin2_string = get_oin_string(mol_gen_bonded, xyz_coords)   # mol_gen_bonded = gen_result.mol
```

`gen_result.mol` carries the **generator's own bond graph**. A haptic ring that has drifted off the
metal is still bonded *in that graph*, so the re-encode reproduces the input's coordination and the
key matches. That document named the consequence and left it open:

> *"Nobody currently knows the false-positive rate with the lever OFF."*

That question outranks every lever in the release. A lever changes a few percent of runtime; this
sets whether the reported accuracy **is** the accuracy. Chasing "100 % round-trip" on an instrument
with an unmeasured false-positive rate is chasing a number that may not be real.

## 2. Method, and the two confounds it avoids

Both sides are re-encoded with the **current** encoder. The report's stored `smiles_1` is *not*
used, so a difference cannot be a version artifact between the code that produced the sweep and the
code reading it.

Classification is by **what changed structurally**, not by whether the string changed. This matters
because `[[reencode-vs-harness-smiles2]]` already showed this path inflates `structural` ~19× via
harmless presentation drift (slot renumbering, fragment order). Charging that to the metric would
have produced a large and completely wrong number. `KEY_DIFF_COORD_OK` isolates it: 10 molecules,
not counted as false positives.

The four charged classes are coordination failures:

| class | test | n |
|---|---|---|
| `GEO_DEGRADED` | the `[El_GEO]` metal geometry tag changed | 54 |
| `HAPTIC_LOST` | a haptic slot present in the input is absent | 4 |
| `HAPTICITY_REDUCED` | a haptic slot survives but binds fewer atoms (η⁵→η²) | 2 |
| `DENTICITY_LOST` | fewer distinct slots overall — a donor detached | 1 |

> **Counting note, because this project's rule is "measure counts, don't trust quotes."** An
> earlier revision of this file quoted 57/6/4/2 and `KEY_MATCH` 565. Those came from grepping the
> shard **logs** for class names — which also matched the class names printed in each shard's own
> summary block, inflating every class by a few. The counts above are recomputed from
> `docs/agentic-notes/v0.4.6/metric_false_positives.json` and sum exactly to 633. The headline rates (9.6 % / 28.1 % /
> 2.8 %) were computed from the archived rows and never moved.

`HAPTICITY_REDUCED` exists because the winding head `{n>}` and the slot number both survive a ring
slip, so the slot-set checks cannot see it. Without it, four genuine failures read as clean.

## 3. Results

```
scored successes measured : 633
FALSE POSITIVES           : 61/633 = 9.6%
  of HAPTIC inputs        : 48/171 = 28.1%
  of NON-haptic inputs    : 13/462 =  2.8%
```

### 3.1 The `GEO_DEGRADED` transitions are coordination-number LOSS, not classifier noise

The obvious objection is that the geometry classifier flipped between two similar polyhedra. It did
not. Almost every transition **lowers the coordination number**:

*(top 14 transitions of 54; the remaining 5 are a one-each tail, all also CN-lowering.)*

| transition | CN | n |
|---|---|---|
| TET → TPL | 4 → 3 | 13 |
| TET → LIN | 4 → 2 | 11 |
| TPL → LIN | 3 → 2 | 10 |
| TPY → LIN | 5 → 2 | 3 |
| SPY → LIN | 5 → 2 | 3 |
| SPY → TPL | 5 → 3 | 2 |
| OCT → TBP / TPL / SPL | 6 → 5 / 3 / 4 | 3 |
| SPL → LIN, TBP → LIN, TBP → TPL, TPY → TPL, SPY → SPL | all lower | 5 |

There is no CN-preserving reassignment in the whole table. Ligands are leaving.

### 3.2 The haptic cases, with distances — the verification that settles it

Perception near the covalent-radius + 0.45 Å cutoff can be marginal, so "the encoder no longer sees
a bond" is not by itself proof of detachment. Measured metal–carbon distances:

| molecule | metal | input M–C (nearest 10) | generated M–C | C within cutoff |
|---|---|---|---|---|
| **FIYHUT_comp_0** | Fe | **2.02 – 2.05 Å** | **2.84 – 2.96 Å** | **10 → 0** |
| ATUROX_comp_0 | Zr | 2.48 – 2.61 Å | 2.91 – 3.47 Å | 10 → 2 |
| FOJWOR_comp_0 | Hf | 2.45 – 2.53 Å | 2.46 – 3.27 Å | 10 → 6 |
| XEVPEU_comp_0 | Ti | 2.13 – 2.48 Å | 2.53 – 2.94 Å | 7 → 4 |

FIYHUT is decisive. Textbook ferrocene Fe–C is ≈ 2.05 Å, and the input matches it to 0.03 Å. The
generated structure puts all ten ring carbons at **2.84–2.96 Å — roughly 0.85 Å too far**. That is
not a threshold judgement; at 2.9 Å there is no Fe–C bond. **Both cyclopentadienyl rings have come
off the iron, and the molecule is recorded as a successful round trip.**

## 4. What this means

1. **Reported accuracy is inflated, and the inflation is structured, not random.** It concentrates
   10× on haptic molecules — which are ~23 % of the corpus and 35.6 % of generator CPU. Any figure
   quoted for the haptic population should be read as **~72 % of its face value**.
2. **The instrument cannot detect its own failure mode.** No existing sweep, bucket report, or
   promotion gate will ever surface these, because they all consume `smiles_2`, which is computed
   from the generator's bond graph. This is why the rate had to be measured out-of-band.
3. **It is independent of `OIN_ACCEPT_SCORED`.** These 61 molecules are the **default** path.
   The lever makes the population larger; it did not create it. The v0.4.7 lane's G2 result (indep
   re-perception 15/20 → 7/20 with the lever on) and this 9.6 % baseline are the same phenomenon
   measured from two directions.
4. **It reframes the "100 % round-trip" goal.** The gap to 100 % is not only molecules that fail;
   it includes molecules that *pass and should not*. Fixing the metric must come before, or at
   least alongside, closing the remaining failure classes — otherwise the closing work is graded
   by a ruler that credits detached ligands.

## 5. BUILT — `src/oinsmiles/oin/coordination.py` (commit `5118f1ef`)

The recommendation in §5.1 below has been **implemented and validated**. `coordination_report()`
compares metal-contact multisets computed from **distances only** — no perception, no re-encode,
and it consults **neither** bond graph, so it cannot be fooled the way the metric is.

| | result |
|---|---|
| FLAG on the 61 known false positives | **55/61 = 90.2 %** |
| false alarm on the 572 genuine passes | **21/572 = 3.7 %** |
| cost | **2.2 ms/molecule** |

Wired into `tools/test_dataset_roundtrip.py` as `report["coordination"]`, recorded *before* the key
comparison returns so it is present on mismatches too. **A diagnostic, never a gate.** End-to-end
through the real harness:

```
ADAMAT_comp_0    status=success  coordination.intact=True
FIYHUT_comp_0    status=success  coordination.intact=False  Fe: lost {'C': 10} (10 -> 0)
```

FIYHUT still scores `success` — no gate changed — but the report now carries the evidence.

**The `MARGINAL_BAND` refinement was forced by measurement.** A raw loss verdict flagged 45 genuine
passes; 36 of those had a contact within 0.10 Å of the cutoff, the worst being η⁶ arenes (Ru 9→3,
Cr 9→4) where a per-atom covalent criterion is stricter than the encoder's ring perception.
Reporting a loss made *entirely* of boundary contacts as `boundary_only` rather than degraded cut
false alarms 7.9 % → 3.7 % for one point of recall. FIYHUT is unaffected — 0.33 Å beyond cutoff is
over 3× the band.

**Update (`5f565b6a`, `5057c630`) — the two limits were attacked; one closed, one refuted.**

*Refuted as a verdict:* `denticity_signature()` groups metal contacts by ligand (connected
components of the non-metal covalent graph — ferrocene is `(5,5)`), needing **no** slot→atom
correspondence, so the "needs correspondence" claim below was wrong. It catches 5 of the 6 misses,
lifting recall to 98.4 % — and fires on **55.8 % of genuine passes**. Eight recall points for a 15×
worse false-alarm rate is a different, useless instrument. It is **recorded** (`denticity_in` /
`denticity_gen`, informative per-molecule) but `intact` stays loss-based, guarded by a test.

*Closed:* the boundary band was **one-sided** — it only saw contacts just *past* the cutoff.
OGARAP_comp_0, the last unexplained miss, retains 3 Pd–C contacts but at margins **+0.084 / +0.028
/ +0.007 Å**; a contact held by 0.007 Å is as ambiguous as one lost by 0.007 Å. With `at_boundary_gen`
the band is two-sided:

| | one-sided | two-sided |
|---|---|---|
| FLAG recall / false alarm | 90.2 % / 3.7 % | **90.2 % / 3.7 %** (unchanged) |
| FLAG-or-BOUNDARY recall | 91.8 % | **96.7 %** |
| **false positives passing with no signal** | 6 | **2** |

Both remaining misses are gain-driven over-coordination. A side finding worth keeping: `BOUNDARY`
covers 210/572 genuine passes, i.e. **better than a third of generated structures hold a ligand
within 0.1 Å of the perception cutoff.**

**Two scope limits, as originally measured:**
1. The verdict is **loss**-based, so gain-driven over-coordination is invisible — 4 of the 61
   *gained* contacts (6→11, 7→12) and changed geometry tag without losing any. No threshold is
   asserted, because a genuine pass in the same corpus gained 2 (Mo 6→8).
2. **Same-count hapticity rearrangement** is invisible — OGARAP goes η³→η² with the carbon count
   unchanged, which an aggregate per-element multiset cannot see.

Together these account for 5 of the 6 misses. Closing them needs slot→atom correspondence, which is
deliberately out of scope here (see the module docstring).

---

## 5.1 The recommendation, as originally written

Add a **coordination-integrity check** to the harness, next to the existing key comparison, and
record it per molecule rather than gating on it at first:

- for each slot in the input OIN, confirm the corresponding atom(s) in the **generated
  coordinates** are within the covalent-radius contact criterion of the metal;
- report hapticity per haptic slot (η-order in, η-order out).

This is the same "cheap attachment check" assessed as feasible in `ACCEPT_SCORED_v0.4.7.md` §6
(11–81 ms per molecule, 6–7 of 8 recoverable). This measurement is independent corroboration that
it is worth building, and it supplies the corpus-wide baseline that proposal lacked: **9.6 %
overall, 28.1 % haptic.**

Record it as a diagnostic first, exactly as v0.4.4 did with RMSD. Gating on it immediately would
move ~61 molecules from pass to fail in one step and make the change indistinguishable from a
regression — the confound this project has already been caught by twice.

## 6. The OTHER half — false NEGATIVES, and the single root cause behind both

Measured the same way with `--status failed` over the 302 reported failures
(`docs/agentic-notes/v0.4.6/metric_false_negatives.json`). Only **51 produced a stored structure**; the other 251 have
nothing to check (timeout / no-conformers / encode failure) and are genuine failures.

Of those 51, **35 key-match under independent re-perception**. That number is *not* the
false-negative rate, and checking why is what keeps it honest:

| why the harness marked it failed | n | reading |
|---|---|---|
| `Atom count mismatch` | 27 | a **separate deliberate gate**, not the string comparison — the known atom-count class |
| `String mismatch` | **8** | **genuine false negatives** — the harness says the string differs, full re-perception says it matches |

So the false-negative population is **8**, not 35: 8/302 reported failures = 2.6 %.

### The two error directions have ONE cause

`YOSXIP_comp_0` shows it exactly:

```
input             ... CS{1}CC[S@]{5}(C)=O ...     chiral sulfoxide sulfur
harness smiles_2  ... CS{1}CCS{5}(C)=O ...        tag LOST  -> scored a mismatch
full re-encode    ... CS{1}CC[S@]{5}(C)=O ...     tag PRESENT -> round-trips
```

The generated geometry carries the correct sulfoxide chirality. The harness cannot see it because
`gen_result.mol` was never run through stereo-perception-from-structure.

That is the **same shortcut** as §1, producing the opposite error:

| | what `gen_result.mol` does | consequence |
|---|---|---|
| §1 | **asserts bonds** the geometry does not support | 61 **false positives** — detached ligands scored as passes |
| §6 | **lacks stereo** the geometry does support | 8 **false negatives** — correct structures scored as failures |

One root cause: scoring through the generator's own molecule instead of re-perceiving the
coordinates. It over-credits connectivity and under-credits stereochemistry at the same time.

### Net effect on the headline number

On this 936-molecule cohort the reported pass count is **inflated by 61 and deflated by 8** — a net
**~53 molecules, 5.7 points**, of over-statement. Both halves must be fixed together, and both are
fixed by the same change: score from a full re-perception of the generated coordinates. That is
precisely what `accept_fn`'s strict step already does, and what `docs/agentic-notes/v0.4.7/ACCEPT_SCORED_v0.4.7.md`
measured as costing runtime — so the cost of an honest metric is now quantified from both sides.

## 7. Reproduce

```bash
export PYTHONPATH=$PWD/src; V=.venv/bin/python
D=tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-rebaseline
for i in 1 2 3; do
  $V tools/haptic_false_positive.py --results-dir "$D" --shard $i:3 --json fp$i.json &
done; wait
# haptic subpopulation only
$V tools/haptic_false_positive.py --results-dir "$D" --haptic-only
# the other half: reported failures that actually round-trip
$V tools/haptic_false_positive.py --results-dir "$D" --status failed
```

Generator-free: it reads stored `structures/*_generated.xyz` and runs no 3D generation, so it is
unaffected by machine load and by the timeout confound that dogs every pass-rate comparison in
this project.
