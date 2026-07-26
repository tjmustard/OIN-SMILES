# The valence-search ORDER — v0.4.5 (`swimlane/v045-valorder`)

Owner: valorder swimlane, based on `swimlane/v045-valsearch` @ `1497cc90`. The git-durable
record; the gitignored `spec/handoffs/v0.4.5/PROGRESS-valorder.md` has the working narrative.

## The lead this lane was given

`docs/VALENCE_SEARCH_v0.4.5.md` closed with a hypothesis it called "one lead that looks better
than what this lane shipped":

> The over-cap branch does not merely cap the search; it disables the heuristic that makes the
> search succeed. […] sorted ligands succeed in under 10 candidates; unsorted ones fail in
> 20 000. `found_valid = 0` on the over-cap arms is then not obviously a statement about the
> ligand being unperceivable — it may be a statement about search *order*.

**The order hypothesis is wrong.** The measurement that refutes it also hands over a fix that
is much better than the one the hypothesis proposed, and it is an accuracy fix.

## Headline

**`found_valid = 0` is not about order and not about budget — it is about DENSITY.** For
QIDKUL's 37-atom ligand exactly **16 of 1 259 712** candidate valence assignments can satisfy
`AC2BO`'s own acceptance predicate, and the earliest sits at **rank 209 858**. No prefix of
20 000, in any order, could ever reach one. But the feasible set is *decidable in closed form*,
because the predicate's charge condition is **additive over atoms** — so you can enumerate the
16 directly instead of grinding 20 000 of the 1 259 696 that cannot work.

Doing that turns a fallback guess into a real Lewis structure:

| | default (raw prefix, 20 000) | `OIN_VALENCE_CHARGE_FILTER=1` |
|---|---|---|
| `QIDKUL_comp_0` whole encode | 124.44 s, `found_valid=0` on all 3 over-cap calls | **0.87 s**, `found_valid` on all of them |
| emitted OIN | `…-n2c{0}n(…)cc2…` | `…N2C{0}N(…)C=C2…` |
| does the string re-read? | **3/4 fragments — the ligand fails `KekulizeException`** | **4/4 fragments** |
| `QIDKIZ_comp_0` whole encode | 129.53 s, `found_valid=0` | **0.32 s** |
| does the string re-read? | **2/3 fragments — same `KekulizeException`** | **3/3 fragments** |

Both molecules are metal–NHC complexes (Rh(COD)Cl and AuCl of a
bis(2,4-dinitrophenyl)imidazol-2-ylidene). The shipped default emits the carbene as a lowercase
**aromatic imidazolium** that RDKit cannot kekulize; the filter emits it as the **neutral
carbene** ring that does re-read. So this is not "a different string" — it is the difference
between a string that parses and one that does not.

## How the question was made decidable

`AC2BO` accepts a candidate only via `BO_is_OK`, and a valid return forces
`BO.sum(axis=1) == valences` **exactly**: either `UA` is empty and `BO = AC`, or
`valences_not_too_large` bounds every atom above by `valences` while `(BO - AC).sum() == sum(DU)`
fixes the total, so every per-atom slack is zero. `charge_is_OK` therefore evaluates
`get_atomic_charge` on the *candidate's own* valences — a pure per-atom function. Write
`Q0(valences) = Σ_i get_atomic_charge(z_i, v_i)`. Its only correction is `Q += 2` per trivalent
single-bonded carbon, and only while running below the target. Hence two necessary conditions:

* **C1** `Q0 ≤ charge` and `charge − Q0` even;
* **C2** `sum(valences) − sum(AC_valence)` even, since every added bond raises `BO.sum()` by 2.

Both are additive, so a suffix DP over `(Q0, parity)` answers, for the **whole** space and
without running the search:

* how many candidates can possibly be valid;
* the exact **rank** of the first one in any lexicographic enumeration order (mixed radix, so
  exact for spaces far beyond 2**64);
* and therefore whether a prefix budget could ever have found it.

`tools/valorder_feasibility.py` does this for all eight known over-cap molecules in **8.4 s of
CPU**. The instrument is validated three ways rather than trusted:

* **brute force.** Enumerating QIDKUL's full 1 259 712-candidate space by hand agrees exactly:
  `survivors=16 first_rank=209858`, in both orders.
* **the real code path.** Every rank the DP predicted was then hit by production `AC2BO` on the
  nose — BENVOG's first valid candidate predicted at rank 14 and found as candidate **15**
  (1-based), ZAZREZ's at 648 found as **649**, BENVOG under the heuristic order at 224 found as
  **225**.
* **the sibling lane.** The probe reproduces `docs/VALENCE_SEARCH_v0.4.5.md`'s independently
  recorded `best_BO` sha `9d234c9103a8` and both whole-encode OIN shas (`0d428d9dfc56`,
  `796221298c6c`).

## Q1 — is `found_valid = 0` an ordering artefact?

**No. Where the two orders differ, the heuristic order is WORSE, and for 3 of 8 molecules it is
not a different order at all.**

The reason is structural. `_ordered_valences` sorts by `(order_idx, tuple)` where `order_idx` is
the position in `product(O_sums, N_sums, C_sums, P_sums, S_sums)` — so **O varies slowest** and S
fastest. In these ligands the atoms with a choice are overwhelmingly **O and N**, so the
heuristic defers changing exactly the atoms that must change.

| molecule | ligand | charge | space | `Q0 == charge` | **C1 (what the code enumerates)** | first feasible rank, raw → heuristic |
|---|---|---|---|---|---|---|
| `QIDKUL_comp_0` | 37 | -8 | 1 259 712 | **16** | **16** | 209 858 → 209 858 *(same order)* |
| `QIDKIZ_comp_0` | 37 | -8 | 1 259 712 | **16** | **16** | 209 858 → 209 858 *(same order)* |
| `LIYFAA_comp_0` | 92 | -10 | 1 679 616 | **0** | **0** | — |
| `HICLAG_comp_0` | 147 | -2 | 53 747 712 | 97 868 | 105 620 | **4 → 32 768** |
| `BENVOG_comp_0` | 148 | -2 | 26 873 856 | 74 108 | 80 548 | 14 → 224 |
| `ZAZREZ_comp_0` | 144 | -2 | 8 503 056 | 1 424 304 | 1 734 696 | 648 → 648 *(same order)* |
| `KESWUB_comp_0` | 188 | -2 | 6 718 464 | 41 020 | 45 358 | 14 → 896 |
| `HOHKUL_comp_0` | 220 | -2 | 26 873 856 | 74 108 | 80 548 | 14 → 224 |

**Two counts, and the wider one is the operative one.** `Q0 == charge` is the strict condition;
**C1** additionally admits `Q0 < charge` with matching parity, because `charge_is_OK` can climb by
`+2` per trivalent single-bonded carbon. `iter_charge_feasible_valences` enumerates **C1** — it
must, or the filter would not be necessary-only — so every "candidate *n*" index below is an index
into the C1 sequence, and C1 is the count a budget has to cover. The ranks are identical under
either condition on these molecules (the first C1 member is also the first strict member).

"Same order" is literal: when every multi-choice atom belongs to one element group, and a group
is traversed in atom order, the heuristic sequence *is* the raw sequence. Single-choice atoms
contribute nothing. So for QIDKUL, QIDKIZ and ZAZREZ the hypothesis is not merely false, it is
**vacuous** — there was no reordering to be had.

`BENVOG` and `HOHKUL` share every number in that table, which is exactly what a cache-keying bug
looks like — the sibling lane shipped one — so it was checked rather than shrugged at. The two
captures are distinct ligands (148 atoms `C76H52N4O16`, `AC.sum` 320, sha `fea3a4d18cac`; 220
atoms `C100H100N4O16`, `AC.sum` 464, sha `4a37c79878f2`) that happen to have the **same
choice-bearing census** — 4 N and 16 O with a choice, everything else pinned. The DP depends only
on that census and the target charge, so identical profiles are the expected result, not a
collision.

Measured on the real code path, at the shipped budget of 20 000:

| molecule | raw | heuristic order | charge filter |
|---|---|---|---|
| `QIDKUL` ligand | 20 000 cands, no valid, 46.0 s | 20 000 cands, no valid, 41.2 s | **VALID at candidate 10, 0.02 s** |
| `LIYFAA` ligand | timeout (19 312 cands / 90 s) | timeout (16 354 cands / 90 s) | 0 feasible — see below |
| `HICLAG` ligand | timeout (5 942 cands / 90 s) | timeout (6 153 cands / 90 s) | **VALID at candidate 1 129, 18.6 s** |
| `BENVOG` ligand | VALID at 15, 0.28 s | VALID at 225, 4.56 s | **VALID at 1, 0.01 s** |
| `ZAZREZ` ligand | VALID at 649, 8.10 s | VALID at 649, 8.80 s | **VALID at 1, 0.02 s** |
| `KESWUB` ligand | timeout (6 033 cands / 90 s) | timeout (6 188 cands / 90 s) | **VALID at candidate 1 216, 22.1 s** |
| `HOHKUL` ligand | timeout (3 356 cands / 90 s) | timeout (3 947 cands / 90 s) | **80 000 of its 80 548 feasible cands examined, none valid**, 2 572 s |

`KESWUB` is the row that shows why "how deep does the budget reach" matters more than "how many
candidates does it examine". The 20 000-candidate budget spends its examinations on *feasible*
candidates once the filter is on, and the 1 216th feasible candidate lies far beyond raw rank
20 000 — so the filter does not merely reorder the same reachable set, it **extends the reach of
the same budget** by the reciprocal of the feasible density.

Note the third and fourth rows: on BENVOG and ZAZREZ all three arms return the **same
`best_BO` sha** (`f4e60eba2807`, `fddae6872a0b`). That is the subsequence property doing its
job — the filter skips candidates, it never reorders them, so where the default already
succeeds the filter agrees with it exactly and merely arrives sooner.

The ordering lever is shipped anyway, default OFF (`OIN_VALENCE_ORDERED_FALLBACK`), so this
refutation is reproducible rather than asserted. `iter_ordered_valences` reproduces
`_ordered_valences`' order **element for element** in O(1) memory — asserted as an equality over
150+ random configurations — which is also what makes "sub-cap is untouched" provable for it.

## Q2 — a second correction to the sibling lane's reading

`docs/VALENCE_SEARCH_v0.4.5.md` recorded, from its fine budget sweep on QIDKUL's ligand:

> **1–25 candidates: genuinely worse.** Sum 92, a different matrix. Cutting this far would
> silently degrade perception.

That inference is **inverted**, and the sha proves it: the sum-92 matrix it flagged as worse is
`53205ad7836d` — which is exactly the matrix the charge filter returns with **`found_valid=1`**.
The sum-94 matrix (`9d234c9103a8`) that the full 20 000-candidate search returns is the
*unvalidated guess*.

Checked atom by atom, both arms of QIDKUL's ligand:

| | raw, 20 000 tries | charge filter |
|---|---|---|
| `found_valid` / `over_cap_exhausted` | 0 / 1 | **1 / 0** |
| candidates examined | 20 000 | **10** |
| `BO.sum()` | 94 | 92 |
| formal charges sum to the target `-8` | yes | yes |
| atoms above their max tabulated valence | none | none |
| SanitizeMol | OK | OK |
| perceived form | aromatic imidazolium ylide, `[c-]` + `[n+]`, 10 formal charges | **neutral NHC carbene**, `[C]`, 8 formal charges |

Both are drawable in isolation. The difference is that only one of them satisfies the
algorithm's own predicate — and only one of them survives being written into an OIN and read
back. **`BO.sum()` is not a quality metric. Validity is.**

### Why the default's string does not re-read, exactly

The guess is a *charge-separated* perception, and the emitted OIN drops the charges. Taking the
emitted ligand fragment with slot markers stripped:

| the fragment | result |
|---|---|
| as emitted, `…-n2cn(…)cc2…` | **`KekulizeException`, unkekulized atoms 8 22 23** |
| with the ylide charges the perception actually had, `…-[n+]2[c-]n(…)cc2…` | **OK** |
| as an imidazolium cation instead, `…-n2c[nH+](…)cc2…` | `KekulizeException` |

So the perception was valid *as perceived* and became unreadable *as written*: an aromatic
5-ring whose carbene carbon has three heavy neighbours, no hydrogen and no charge has no
electron to give the ring. Restoring the metal bond that `{0}` denotes cannot fix it either —
that makes the carbon 4-connected, which is still not kekulizable. The filter's perception needs
no ring charges at all, so the question does not arise.

A caveat on the instrument, because it also reports failures that are **not** defects: the
re-parse check strips slot markers and metal tags, so it is naive about some OIN spellings
(`[cH]` with an explicit hydrogen, `C#O` carbonyls) and flags a few sub-cap molecules in **both**
arms identically. It is therefore evidence of a *difference between arms on one molecule*, which
is what is claimed here, rather than an absolute grammar check.

### The rescued structures, checked rather than assumed

`HICLAG` and `KESWUB` produce **no** string on the default path (both order arms time out), so
there is no before/after diff to show. What can be shown is that the structure the filter finds
is sound on its own terms:

| | `HICLAG` (147) | `KESWUB` (188) | `BENVOG` (148) |
|---|---|---|---|
| `found_valid` at candidate | 1 129 | 1 216 | 1 |
| `BO.sum()` / `best_BO` sha | 392 / `7ff80f862e9d` | 450 / `3deb2077e1c0` | 398 / `f4e60eba2807` |
| per-atom formal charges sum to the target `-2` | **yes** | **yes** | **yes** |
| atoms above their max tabulated valence | none | none | none |
| `SanitizeMol` | OK | OK | OK |
| nonzero formal charges | 2 × `O⁻` | 2 × `O⁻` | 2 × `O⁻` |

Each sha matches the one the independent ligand probe recorded, so the verification and the
measurement agree.

One nuance worth recording, because it is the sort of thing that looks like a charge bug and is
not. `charge_is_OK` (and `set_atomic_charges` with it) applies a `+2` correction per trivalent
single-bonded carbon while the running total is below the target. On the four validated
structures above the per-atom charges hit the target **exactly**, so that correction never
fires. On the 200-candidate `best_BO` guesses captured for `HICLAG` and `KESWUB` they sum to
`-4` against a target of `-2`, and only reach the target *through* the correction — which writes
`+1` onto a carbon `get_atomic_charge` calls `-1`. The net charge ends up right either way; the
difference is that the guess needs a heuristic patch on two carbons and the validated structure
needs none.

## Q3 — does it help the large-ligand class?

**Partly, and the split is not the one the ligand-size table predicted.**

`docs/VALENCE_SEARCH_v0.4.5.md` grouped the population by ligand atom count and concluded the
148–220-atom class was bound by *per-candidate cost*. On this lane's measurements the binding
constraint is **feasible-candidate density**, and it cuts across size:

| class | molecules | what the filter does |
|---|---|---|
| dense feasible set, valid early | `BENVOG` (148), `ZAZREZ` (144) | already succeeded on the default path; filter reaches the **same** answer (identical `best_BO` sha) in 1 candidate instead of 15 / 649 |
| sparse feasible set, valid reachable | `QIDKUL`, `QIDKIZ` (37), `HICLAG` (147), `KESWUB` (188) | **timeout or guess → real Lewis structure.** 10, 10, 1 129 and 1 216 candidates |
| **no feasible candidate at all** | `LIYFAA` (92) | provably hopeless at this charge — see below |
| effectively refuted | `HOHKUL` (220) | **99.3% of its feasible set searched (80 000 / 80 548), none valid**; `best_BO` never left `AC.copy()` |

The ligand-size split in `docs/VALENCE_SEARCH_v0.4.5.md` does not survive: **`KESWUB` at 188
atoms is rescued and `BENVOG` at 148 was never failing on this axis**, while the 37-atom pair is
among the *hardest* by feasible density (16 candidates in 1.26 M). Size predicts per-candidate
cost; it does not predict whether the search can succeed. Density does.

Note also that per-candidate cost measured here is **0.012–0.027 s**, not the ~0.36 s
`docs/VALENCE_SEARCH_v0.4.5.md` reported for `HICLAG`. That figure came from whole-encode wall
divided by candidates and so carried the rest of the encode; on the ligand alone the largest
over-cap ligand costs about 25 ms per candidate. The 70x spread between ligand sizes it inferred
is really about 5x.

## Q4 — `LIYFAA`, and the finding that closes a question permanently

**At charge -10 with carbenes allowed, `LIYFAA_comp_0`'s 92-atom ligand has ZERO candidates that
can satisfy the predicate — out of 1 679 616.** Not "none in the first 20 000": none at all. No
budget, no ordering and no matcher can perceive it there. The default path nevertheless spends
its full budget discovering this (19 312 candidates in 90 s before the probe's cap, in both
order arms).

That relocates the problem. `LIYFAA` is not a search failure; it is a **charge-proposal**
failure. The extended-Hückel proposal handed `AC2BO` a target the ligand cannot reach, and the
ladder's job is to move off it. Worth its own lane; the DP makes "is this charge even
reachable?" a millisecond question, which is a cheap gate the ladder does not currently have.

## What is shipped, and the gating

**Shipped:** `iter_ordered_valences`, `iter_charge_feasible_valences`, three new counters
(`over_cap_ordered_calls`, `over_cap_filtered_calls`, `over_cap_infeasible`), and three tools
(`valorder_probe.py`, `valorder_feasibility.py`, `valorder_encode_ab.py`). Two levers, **both
default OFF**:

* `OIN_VALENCE_ORDERED_FALLBACK` — take the bounded prefix in `_ordered_valences`' order.
  Measured **worse**; kept only so the refutation above is reproducible. Nothing should ship
  behind it.
* `OIN_VALENCE_CHARGE_FILTER` — enumerate only candidates that can possibly be valid.

Both are read **only inside `if over_cap:`**, so sub-cap ligands — 99.8% of the corpus per the
sibling lane's scan — cannot reach them and are byte-identical **by construction**, the same
structural argument `OIN_VALENCE_FALLBACK_TRIES` rests on. Lever reads prefer
`oinsmiles/oin/levers.py` when it is present (it is **not** on this branch's base; it arrives with
`trial/v045-merge2`) and otherwise use identical local semantics. Either way
`OIN_VALENCE_CHARGE_FILTER=0` **disables**, closing the sense-inversion trap that registry exists
for.

That forward compatibility was checked rather than hoped for: dropping the registry file in and
importing both from the package root and from the submodule resolves `_lever_enabled` to
`oinsmiles.oin.levers` with **no circular import** (`oinsmiles/oin/` is a namespace package, so
nothing re-enters `oinsmiles/__init__.py`), and unset still reads OFF because neither lever is in
`_DEFAULT_ON`. When the registry lands, the local fallback can be deleted; until then it is not
a duplicate of the registry, only of its two-line semantics.

### The filter is designed to be strictly dominant

Two properties, both tested:

1. **It is a subsequence of the raw product, in the same relative order.** It skips candidates;
   it never reorders them. Every candidate a valid perception could use is still yielded at its
   original relative position, so *the first valid candidate found is the one an unbounded raw
   search would have found*. BENVOG and ZAZREZ confirm this empirically — identical `best_BO`
   sha across all three arms.
2. **When the feasible set is empty it falls back to the historical enumeration.** `best_BO` is
   assembled from candidates the filter drops (its own charge test is on the *BO*, not on the
   candidate), so without this the lever could turn one guess into a *different* guess.
   With it, the lever's entire blast radius is "a guess becomes a real Lewis structure".
   `LIYFAA` is the real molecule this covers.

Property 2 costs `LIYFAA` its speedup: it keeps grinding a budget the DP has already proved
cannot succeed. Short-circuiting that is a **separate, provably safe optimisation** — the same
DP that proves the grind is pointless would authorise skipping it — but it changes the fallback
answer, so it is left as an explicit decision rather than folded in silently.

### Correctness claims and how each is pinned

`tests/unit/test_valence_order.py`, 14 tests:

| claim | how it is pinned |
|---|---|
| `iter_ordered_valences` == `_ordered_valences` | order equality over 150+ random configurations, with a guard that the comparison is not vacuous |
| the charge filter never drops an acceptable candidate | **brute force** over 11 whole small candidate spaces, running the real `BO_is_OK`/`get_BO` predicate on every member, with a guard that some candidate was actually accepted |
| output is a subsequence, not a reordering | positions in the raw product are checked to be increasing |
| the filter actually prunes | a filter that keeps everything would pass every test above |
| an infeasible charge reproduces the default answer | `best_BO` and candidate count compared arm to arm |
| both levers default OFF, and `"0"` disables | counters checked on a real over-cap `AC2BO` call |
| sub-cap perception cannot change | ethene encoded under every lever combination |

## Sub-cap byte-identity, measured as well as argued

Structural argument first: the lever reads are inside `if over_cap:`, so a sub-cap ligand cannot
observe them. Measured anyway, filter OFF vs ON, whole encodes in one process with the OFF arm
repeated as a determinism self-check:

| result | count |
|---|---|
| OIN identical | **16 / 16** |
| candidates examined identical | **16 / 16** |
| `REPEAT-OK` | **16 / 16** |
| `over_cap_calls` | **0** on every molecule and every arm |
| OIN changed | **0** |

Twelve dataset molecules (a seed-3 sample confirmed sub-cap by `tools/valsearch_scan.py`, 0 of 53
ligand fragments over cap) plus all **four goldens** — CisPlatin, Ferrocene, fac-Ir(ppy)3,
PdCl2-R-BINAP. Sub-cap ligands examined **4 to 6** candidates each, so they never come close to
the over-cap branch, exactly as the code structure implies.

```
$V tools/valorder_encode_ab.py --dataset <dir> --mols JIQBET_comp_0,… \
   --files "tests/fixtures/CisPlatin.xyz,tests/fixtures/Ferrocene.xyz,…" --cap 150
```

## Suite and lint

Two full `discover tests/unit` runs, one per tree state, because the tree moved mid-run:

| tree | tests | result | wall |
|---|---|---|---|
| `0790946d` (14 lane tests) | **632** | OK, 3 skipped, 3 expected failures | 890 s |
| final, `7fd10555` (16 lane tests) | **634** | OK, 3 skipped, 3 expected failures | 975 s |

632 = the `valsearch` baseline of **618** + 14; 634 = 618 + 16, the transition-metal guard having
added two. **Skipped and expected-failure counts are identical in both runs and identical to the
baseline**, which is the check that matters: a new xfail or skip is how a regression hides in a
green suite. `tests/unit/test_regression_stability.py` is inside both runs and green on the
default path. `uvx ruff@0.15.20 check` and `format --check` clean across `src`, `tools` and
`tests`.

Honest gap, same one the sibling lane recorded: **no pre-change baseline was run on this host.**
The 618 figure is `valsearch`'s measurement on the branch this one is based on, and the evidence
that nothing regressed is that arithmetic together with the unchanged skip/xfail counts — not a
paired before/after run. The release sweep is holding six cores throughout.

## Verdict on the ordering hypothesis

**`found_valid = 0` is not an ordering artefact. It is a density artefact, and for one molecule
it is a genuine impossibility.**

* Ordering: refuted. Vacuous for 3 of 8 molecules, and strictly worse for the other 5.
* Budget: refuted as well. QIDKUL's earliest possibly-valid candidate is at rank 209 858, so no
  affordable prefix could work; and `LIYFAA` has none at any rank.
* **A valid Lewis structure does exist, and the encoder can now reach it,** for QIDKUL, QIDKIZ,
  HICLAG, BENVOG, ZAZREZ and KESWUB — six of the eight, four of which the default path cannot
  perceive at all today. For `HOHKUL` the answer is very nearly **"no"**: 80 000 of its 80 548
  possibly-valid candidates were examined (2 572 s) and **not one** is accepted — and `best_BO`
  never left `AC.copy()`, so no candidate even qualified for the fallback update. The last **548
  (0.7%)** are unexamined, and closing them costs a fresh 43-minute run for 0.7% of a molecule, so
  it is left stated rather than spent. `HOHKUL`'s defect is therefore **not** the enumeration; it
  is that `get_BO`'s greedy iterated matching saturates none of its feasible assignments.
* For `LIYFAA` the answer is **no, provably, at that charge** — the "closes the question
  permanently" outcome the brief asked for. It relocates that molecule to the charge proposal.

### Recommended, not taken

Promote `OIN_VALENCE_CHARGE_FILTER` to default ON. The case:

* sub-cap ligands cannot reach it (read only inside `if over_cap:`), so 99.8% of the corpus is
  byte-identical by construction, and 16/16 molecules measured agree;
* where the default already succeeds, the filter returns the **identical** `best_BO`
  (BENVOG, ZAZREZ), because it is a subsequence and not a reordering;
* where the default fails, it produces a structure that satisfies the algorithm's own predicate
  instead of a guess — and on QIDKUL/QIDKIZ the guess emits an OIN that does not re-read;
* when nothing is feasible it falls back to the historical enumeration, so it cannot turn one
  guess into a different guess.

The case for caution, which is why the flip is left as a product call following the Y2 axial
precedent:

* it **changes emitted strings** for the over-cap population (2 measured, 4 more that go from no
  string at all to a string), so any hardcoded expectation for those refcodes moves;
* the population is ~0.2% of the corpus, and the counterfactual for `HICLAG`/`KESWUB`/`HOHKUL`
  cannot be diffed — there is no baseline string to compare against;
* the `get_BO` matching defect the sibling lane and this one both brushed against is untouched,
  so a promotion should not be read as "the over-cap class is solved".

The lead the sibling lane recorded was worth chasing and was wrong. What made the difference was
not implementing it — that took an hour — but noticing that the acceptance predicate is additive
and can therefore be *counted* instead of *sampled*. Every claim above is a closed-form count or
a counter, not a wall clock, on a host running at load 30+.
