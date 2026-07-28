# v0.4.12 — reflection parity, and an acceptance lever moved to where it can pay

> **Both levers ship default-OFF. The default path is byte-identical and the headline is
> FLAT at 72.46% by construction.** This release buys two things: the filter v0.4.11's
> refutation demanded, and a correction to where the eta runtime lever was installed.

Lane detail: [`LANE-01-parity-veto.md`](LANE-01-parity-veto.md) ·
[`LANE-02-eta-accept-exit.md`](LANE-02-eta-accept-exit.md)

---

## 1. What this release set out to do

The `LADDER DECISION 2026-07-28` re-pointed v0.4.12 at two lanes, one each:

- **L1 — reflection parity.** v0.4.11 built the donor fold, measured **+7.86 `byte_exact`
  points across 393 molecules**, and refuted it: it **collapses enantiomers in 221 of those
  same 393 gains**. Its close-out named the next step, and this is it.
- **L2 — honest acceptance.** Chartered as *"run the `OIN_ETA_EARLY_EXIT` corpus A/B its
  promotion gate demands"*.

Both lanes ended somewhere other than where the charter pointed. That is now four consecutive
releases; it is the process working.

---

## 2. The finding that generalizes beyond this release

v0.4.11 established that **`byte_exact` can be raised by deleting information, and the
comparison key will agree**. v0.4.12 adds the companion, and it bit *inside this lane's own
implementation*:

> ### 🔴 A safety check that declines and a safety check that works produce the same output.
>
> The first veto encoded the mirror with the fold **inherited** rather than forced off. That
> made `s_rot_m` identical to `s_fold_m`, reduced the left conjunct to *"did the fold fire?"* —
> true by construction — and **disarmed the achiral guard**. The self-check then tripped on
> **18 of 18** movers, so the lever emitted the rotation-only labeling every single time.
>
> **All three fixture tests passed.** Declining to fold *also* separates a mirror pair.

The remedy is not "be more careful". It is that the decision must be **observable**:
`fold_parity.resolve` records *why* it decided — `fold_inactive`, five distinct `declined_*`
reasons, `vetoed_collapse`, `allowed_preexisting_fold`, `allowed_separation_survives` — and the
tests assert the **outcome**, not the string. The corpus number that proves the veto is alive
is **`declined_*` = 0**, not the separation count.

This is the third instance of the same shape in three releases (v0.4.7's attachment check was a
silent no-op whose complete A/B "reported what a genuine null result looks like"; v0.4.11's
four instruments were each blind to reflection). The invariant to carry forward:

> **Before quoting an instrument, ask what a broken version of it would print. If that is the
> same thing, you have not measured anything yet.**

The same trap appeared one layer out, in the *measuring tool*: `fold_transition_sim.py`'s first
draft implemented its veto arm as a string operation, which — because the veto lives one level
up — would have reported a clean, plausible, false *"the veto costs nothing"*.

---

## 3. Reflection parity is not a property of the string

The tempting fix is to require each donor swap to be a proper rotation of the polyhedron,
reusing `derive_rotation_group`'s `det > 0` test. It is wrong **before it is built**: a donor
swap is a transposition fixing every *other* vertex, which for a rank-3 polyhedron does not
preserve the Gram matrix at all — so that test rejects **every** swap and the fold degenerates
to the identity.

The fold is justified by a symmetry of the **ligand**, realized as a proper rotation of the
whole complex. That is not a property of the coordination graph and cannot be read off the
emitted string. Hence the veto lives in `get_oin_string`, on the **pristine** conformer
captured before `_align_to_pai` — which can itself reflect, and mirroring an already-reflected
conformer composes two improper maps into a proper one, so the veto would never fire while
looking perfectly healthy.

A direct consequence, recorded because v0.4.11's close-out predicted otherwise:
`TestDonorFoldCollapsesEnantiomers` was to be *"inverted, not deleted"*. **It is not inverted.**
`canonicalize_oin_slots` still collapses its two strings and always will. **The fix cannot live
where the defect is visible.**

---

## 4. The charter's Lane 2 was measuring a lever that cannot work

`OIN_ETA_EARLY_EXIT` sits in `_select_by_geometry_impl`, downstream of a fully-filled pool. The
only site consulted per conformer *during* filling is `accept_fn`. Its own in-code A/B already
read *fires, attempts 32 → 32*. **Its promotion gate is void, not unrun** — and running the
chartered corpus A/B would have spent real CPU confirming a documented structural null.

The criterion moved to `accept_fn` as `OIN_ETA_ACCEPT_EXIT`, conjoined with geometry
classification and the attachment check — because winding alone would bypass clash-first
ranking, which is exactly the defect `OIN_ACCEPT_SCORED` has.

Then the *cohort* turned out to be stale too: all 8 molecules of the v0.4.6 accept-gap cohort
now satisfy the key inside `accept_fn`, so a pilot A/B on them was flat by construction
(telemetry: `adapter.early_exit_hit` in both arms, `adapter.eta_accept_*` never). Re-derived
from the frozen sweep, the real population is **405 eta molecules whose key never matches, 378
of them over 30 s**.

---

## 5. A licence that was assumed for three releases is now measured

The carry-forward licence lets a release skip a 55 CPU-h sweep when it moves no default answer.
v0.4.9, v0.4.10 and v0.4.11 each *claimed* byte-identity; nothing re-checked it against the
frozen corpus.

Lane 1's veto arm carries a **drift control**: it re-encodes the frozen corpus's inputs *and*
its stored generated structures with today's encoder and requires them to reproduce the v0.4.8
strings before using them. They do. So the licence now rests on measurement rather than on
three releases of inherited assertion — and both of v0.4.11's headline anchors (393 / +7.86 pts,
and 19/250 mirror regressions) reproduce **exactly**.

---

## 6. Predicted vs actual

The charter's predictions were written **before** either lever was built.

| | predicted | actual | |
|---|---|---|---|
| `byte_exact`, default path | FLAT 72.46% | **FLAT 72.46%** | ✓ both levers OFF |
| `byte_exact`, L1 lever ON | **+3.0 to +3.5 pts** (≈150–175 of 393 survive) | **+3.42 pts, 171 of 393** | ✓ |
| mirror audit, uniform 250 | 19 → **0** | **19 → 0 at seed 7 AND seed 11** | ✓ |
| `> 30 s`, default path | FLAT | **FLAT** | ✓ nothing default-path changed |
| eta `> 30 s`, L2 lever ON | down; magnitude unpredicted | **11 → 5; median −69.5%** | ✓ |
| suite | ≥ v0.4.11's floor | **988 OK** (skipped 3, xfail 5) | ✓ |

Two of these deserve more than a tick.

**The L1 prediction was not a lucky guess, and it is corroborated by an independent route.**
v0.4.11 bounded the safe set at *"at most ~172 of 393 (~3.44 of the 7.86 points)"* by counting
**collapses** in a mirror audit. This release counted **survivors** through the shipped
predicate and got **171 / +3.42**. Two different mechanisms agreeing to a single molecule is
the strongest evidence here that the veto separates the right set — stronger than either number
alone.

**L2's runtime thesis is confirmed and L2 still does not ship.** On the re-derived population
the lever takes the median **61.8 s → 18.82 s** and `> 30 s` **11 → 5**, with `RIRYOJ_comp_0`
going 119.04 s → **4.4 s** for a byte-identical string. Then the fifth gate arm failed:
`KIHHUG_comp_0`'s arms returned structures whose metal-configuration descriptors differ
(`''` → `|mc:-|`) with **byte-identical output**, and it is also the run's only `indep`
regression. Four of the five gates call that molecule unchanged.

> A large, real speedup that costs something **four of five instruments cannot see** is exactly
> the trade this release exists to be able to detect. The lever stays OFF.

The first attempt at this A/B was run at load 26 and **thrown away** rather than banked — the
standing trap is *never interleave timing runs with gate runs*. The two runs disagree by enough
to have wrecked the conclusion: `UQUXAG_comp_0` reads **17.93 s** under load and **11.06 s**
clean, a 62% inflation comparable to the whole effect being measured.

---

## 7. What v0.4.13 inherits

1. **The GAVSED class.** `_select_by_geometry`'s fallback ranking is not attachment-aware, so
   the check guards *acceptance*, not *return*. Closing it changes arm A's behaviour and needs
   its own gate.
2. **`PREFILTER_VETO`, acceptance-side only.** Its scoring half was closed by `OIN_INDEP_SCORE`
   in v0.4.8 — on the frozen corpus, cheap-fails-but-independent-passes is **28/5000** and the
   honest metric already counts those correctly. What remains is the generator-side veto, with
   the MEDZUR class.
3. **A stale-cohort audit.** The v0.4.6 accept-gap cohort was wrong about 8 of 8 molecules.
   Any cohort frozen before v0.4.8 should be re-derived, not reused.
