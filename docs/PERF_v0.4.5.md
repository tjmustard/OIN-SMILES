# v0.4.5 perf lane — findings

Owner: perf swimlane (`swimlane/v045-perf`, based on `e8b603d5`). Goal: per-molecule round-trip
under 30s without moving accuracy. This file is the git-durable record; the gitignored
`spec/handoffs/v0.4.5/PROGRESS-perf.md` has the incremental working narrative.

## Summary

The largest measured cost in the slow tail is **not** what the release handoff's "prior art"
section named (haptic scale-sweep, CBC, ETKDG dead retries — all already optimized on this
branch's base). It is a **duplicate expensive re-encode**: `metallogen_adapter.py`'s SL1
generate-until-key-exact early-exit (`OIN_EARLY_EXIT`, default ON since v0.4.4) tests every
newly-embedded conformer against the requested OIN's canonical key via
`_reencode_key_matches`, whose confirming step (`_reencode_oin`) is a full
`XYZToSMILES().convert()` re-perception of the generated geometry — measured at **48-57 seconds
per call** on an eta/haptic test case (`QIDKUL_comp_0`, 59 atoms). That test runs once inside
`generate_3d_structures`'s pool-fill loop (`accept_fn`) for every conformer as it's accepted, and
then — for every conformer that reaches it — **a second time**, redundantly, inside
`_select_by_geometry_impl`'s own early-exit re-scan over the same returned pool right after.

Every mol in the pool `_select_by_geometry` sees (other than ones drawn from the untested
`stereo_rejects` fallback) was already tested by `accept_fn` when `generate_3d_structures`
accepted it — a mol that HAD matched would have short-circuited the whole pool-fill and been
returned alone (`return [early_hit]`, `generator3d/__init__.py:590`), so anything that survives to
`_select_by_geometry` is provably already a confirmed non-match. The second test recomputes an
already-known answer.

## The fix

`metallogen_adapter.py`:
- `_reencode_key_matches(..., cache=None)` — new optional param. When provided, the expensive
  `full = _reencode_oin(m)` step is memoized in `cache` keyed on `id(m)`.
- `_select_by_geometry_impl` / `_select_by_geometry(..., reencode_cache=None)` — thread the cache
  through to their own `_reencode_key_matches` call.
- `MetalloGenAdapter.generate()` — creates one `_reencode_cache = {}` per `generate()` call,
  passes it into the `accept_fn` closure (so the pool-fill loop populates it) and into the
  `_select_by_geometry(..., reencode_cache=_reencode_cache)` call right after (so the re-scan
  reads it back).

### Why this is safe / byte-identical

- `_try_accept` (`generator3d/__init__.py`) returns the *same* mol object it appends to
  `successful_mols`, and `accept_fn` is called on that exact object — object identity survives
  from the fill loop into the `mols` list `_select_by_geometry` receives, for the default FF-only
  path (`optimizer=None`).
- Default `cache=None` at every call site that doesn't opt in (there are none left un-threaded in
  this codebase, but any future direct caller of `_reencode_key_matches` gets the pre-existing
  behavior for free) — behavior is unchanged whenever caching isn't wired up.
- A cache **miss** recomputes exactly as before (same function, same inputs) — the memo can only
  ever make an answer arrive faster, never different. This is the same pattern as the existing
  PuLP topology memo (`compute_chg_and_bo_pulp.py`) and the `alt_cache` memo in
  `generator3d/embed.py` — this codebase's established, previously-validated way to remove
  provably-redundant work.
- Degrades gracefully (falls back to a clean miss, i.e. no speedup but no wrong answer) whenever
  object identity does *not* survive between the two sites — e.g. the `optimizer` (xtb/g-xTB)
  path constructs new mol objects via `ASEOptimizer.optimize`, so `id(m)` won't match; also the
  `stereo_rejects` fallback path, whose mols were never tested by `accept_fn` in the first place.
  Both cases just recompute, same as pre-fix.
- Ungated: this is a pure dead-work removal, not a search-behavior change, so per the release's
  own rule ("a pure dead-work removal that is *provably* byte-identical may land ungated") it does
  not need an env lever — unlike, say, reordering `scales_for_haptic`, which WOULD change which
  conformer wins and must be gated if ever attempted.

## Measured

| what | before | after | note |
|---|---|---|---|
| `reencode_oin_full_calls`, `QIDKUL_comp_0`, `max_attempts=2` | 4 | 2 | exact halving, deterministic counter |
| cost of the 2 eliminated calls | ~56s each | 0.04-0.05s (cache hit) | same session, same code path |
| `generate_s`, same bounded scenario | 212.8s | 116.0s | contended host (load ~7-13); same-code-version before/after comparison, not cross-host |

Full unit suite (`discover tests/unit`, 605 tests) green; `tests/integration/test_roundtrip_smoke.py`
green; `ruff check` / `format` clean.

## Scope note: this is bigger than, and different from, the handoff's eta hypothesis

The release handoff attributed eta's slowness to the haptic machinery (4x `scales_for_haptic`
sweep, `ETA_SELECT_POOL` pool widening, no C++ batch path for haptic) and estimated an upper bound
of ~23% off total runtime. Direct per-call tracing on `QIDKUL_comp_0` shows the embed/scale-sweep
itself costs ~1-2s per pool-fill attempt — cheap. The ~50-210s per attempt instead comes almost
entirely from the acceptance-check re-encode described above, which fires on **every** rejected
conformer regardless of whether it came from a haptic scale sweep or a plain DG embed. This
duplication is not eta-specific in mechanism (it fires for any molecule using the default
`early_exit=True` path with more than a trivial rejected-conformer count); eta molecules simply
generate more rejected conformers before a winding-correct one appears, so they pay the cost most
often and most visibly. A non-eta molecule with an unusually large rejected-conformer count before
its geometry/energy pick would show the same pattern — see the still-open `QIDKIZ_comp_0` question
in the progress file.

## What's still open (see PROGRESS-perf.md for the live version)

1. Confirm byte-identity (fix vs `git show HEAD:...` pre-fix) across eta + non-eta molecules —
   in flight as this doc is written.
2. Re-measure the stratified slow-tail sample with the fix in place, attempt-count-bounded (not
   wall-clock-bounded — that approach wasted real session time chasing a moving contention
   target), to quantify the aggregate win.
3. The handoff's own named eta costs (scale sweep, pool width, no batching) are still real and
   still unclaimed; revisit once this bigger confound is out of the measurement.
4. Encode-side cost (`XYZToSMILES().convert()` itself measured 48-71s on this same molecule) is
   out of this lane's scope but flagged: an eta molecule whose bare encode alone exceeds 30s
   cannot hit the 30s round-trip target regardless of any generation-side fix.
