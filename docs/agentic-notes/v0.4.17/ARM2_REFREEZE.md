# v0.4.17 — ARM 2 golden re-freeze

> **51 of 425 rows re-frozen, from a FULL gate run, each with one recorded reason.**
> 40 moved because of `OIN_EXACT_DONOR_FOLD`. **15 were already stale when v0.4.17 started** (4 are
> in both sets). Re-running only the predicted movers — the v0.4.13 / v0.4.14 procedure — would
> have missed **19 of the 51 (37%)**.

Branch `research/v0417-encoder-canonicality`, 2026-09-19. Owner request: *"Go ahead with the ARM 2
golden re-freeze"*. Lane record: `L1_EXACT_DONOR_FOLD.md`.

| golden | rows | re-frozen | `LEVER` | `STALE+LEVER` | `STALE` | killed, not re-frozen |
|---|---|---|---|---|---|---|
| `tools/gate_v047_arm2_golden.tsv` | 100 | **6** | 4 | 1 | 1 | 0 |
| `tools/gate_v049_arm2_golden.tsv` | 325 | **45** | 32 | 3 | 10 | 2 |
| | 425 | **51** | 36 | 4 | 11 | 2 |

`MANIFEST_SHA256`: v047 `59fe3710…dd038`, v049 `be0c2398…b619` (recomputed; ARM 2 never verifies it).

## 1. Why the full gate and not "the movers"

The precedent re-runs the rows a sweep predicts will move (v0.4.13: 11 of 325; v0.4.14: 7 of 325,
1 of 100). This lane had such a list — 6 + 30 rows. Before spending generator time, a free check:
does `sha256(smiles_1)` from the **record sweep** reproduce the goldens it supposedly matches?
It did not — 2 + 13 rows differed — and `KELQOI` was on no mover list at all. Either those golden
rows were stale or the sweep's encode path differs from ARM 2's, and the two readings lead to
different *recorded reasons* for the same row. So the question had to be settled first.

**Step 1 — every row's field 2, encode-only** (`tools/v0417/arm2_field2_audit.py`). Field 2 is a
pure function of the input file and the encoder, so all 425 rows cost under 20 minutes at 6 jobs. Each molecule
is encoded twice in a **fresh process** (perception memos are process-lifetime; levers are read at
import): shipped (no `OIN_*` set) and `OIN_EXACT_DONOR_FOLD=0`.

| | v047 (100) | v049 (325) |
|---|---|---|
| `SAME` golden = off = shipped | 94 | 288 |
| `MOVED_BY_LEVER` golden = off ≠ shipped | 4 | 20 |
| **`STALE`** golden ≠ off = shipped | **1** | **11** |
| **`STALE_AND_MOVED`** | **1** `JUXPID` | **1** `JUXPID` |
| **`HEALED`** golden ≠ off, shipped = golden | 0 | 1 `NASZOY` |
| `SENTINEL` (`NO_ENCODE@300s`) | 0 | 4 |
| `ENCODE_FAILED` | 0 | 0 |

*What a broken audit would print:* a worker importing the main checkout's `oinsmiles` reports
shipped = off on every row — a flawless "nothing moved". Every worker reports the package path it
imported and the driver refuses a foreign one; the run aborts if the lever fired on 0 rows (it fired
on 5 and on 22). The audit also agrees **row-for-row** with `sha256(smiles_1)` from both stored sweeps, which
were produced by the harness, not by this tool.

**Step 2 — the gate itself, all 425 rows** (`run_arm2_refreeze.sh on`; commit `f5f507ef`; shipped
defaults, `--timeout 300`, 6 shards, BLAS=1; 75 min wall).

| | v047 | v049 |
|---|---|---|
| fully compared, reproduces | 94 | 234 |
| sha1-only (`NO_STRUCTURE@300s`), reproduces | — | 35 |
| deterministic no-structure, reproduces | — | 5 |
| observed only (`NO_ENCODE@300s`) | — | 4 |
| **field-2 mismatch** | **6** | **32** |
| **field-3-only mismatch** | 0 | **13** |
| killed at `--hard-timeout`, scored as a field-2 mismatch | 0 | 2 |

The gate's field-2 mismatches and the audit's `shipped ≠ golden` are **the same 6 and the same 32
rows** — two code paths, one answer; `arm2_refreeze.py diff` aborts if they ever differ.

**Step 3 — the control** (`run_arm2_refreeze.sh off`; commit `6cf9b185`; 6 + 47 rows,
`OIN_EXACT_DONOR_FOLD=0`; 21 min wall). Generation is **seeded** (`MetalloGenAdapter seed=42`), so a healthy old
row comes back byte-for-byte with the lever off:

- `LEVER` — the control reproduces the old row, shipped does not. **36 rows; all 36 came back.**
- `STALE` — the control does *not* reproduce the old row, and equals shipped. **11 rows.**
- `STALE+LEVER` — not reproduced, and differs from shipped. **4 rows** (`JUXPID` ×2, `HIKCAI`,
  `HIKCIQ`; the last two are stale on field 2 and lever-moved on field 3).

The determinism this rests on was measured, not assumed: the shipped run reproduced **all 368
untouched gated rows** (328 on both fields), six shards wide.

## 2. The 15 rows that were stale before this lane

`KELQOI JUXPID` (v047) · `HIKCAI HIKCIQ HIMFAM IMEPAS JUXPID KUPVEX MEXWUH PAKWAB QIDHEQ UNEGEB
WAQXAP XENZAS` on field 2 and `KIHHUG` on field 3 (v049). `NASZOY` was a 16th until this lever moved
it back onto the golden's string.

For every field-2 one: **the last stored sweep whose `smiles_1` hashes to the golden is a v0.4.8
sweep (2026-07-26); the first that does not is the v0.4.14 sweep (2026-07-29); none since does.**
The sweeps only bracket it (v0.4.9 – v0.4.14); §6 names it — **v0.4.13's donor fold** — and they were
left behind by re-freezes that re-ran only predicted movers. ARM 2 is run in full only at a release (`--band fast` is the routine control), so
a red full gate is invisible between releases. §6 pins which lever.

## 3. Predicted vs measured

| | v047 | v049 |
|---|---|---|
| predicted movers (A/B union list ∩ golden) | 6 | 30 |
| lever-moved, measured | 5 | 35 |
| predicted **and** moved | 5 | 27 |
| predicted, row did **not** move | 1 `DOGQEX` | 3 `GUXDUZ HAMKOV NASZOY` |
| lever-moved, **not** predicted | 0 | 8 `FOSNEI HACYEQ HIKCAI HIKCIQ KEHNOA ULORAQ ZELKUY ZOYYUJ` |
| stale-only, on nobody's list | 1 | 10 |

All 8 unpredicted lever rows move on **field 3 alone**: the input string is unchanged, the seeded
generator builds the same structure, and it is the *re-encode of that structure* the fold relabels.
The generated-side mover list was derived from the sweep's structures through the harness's scorer;
ARM 2 scores `get_oin_string(result.mol, …)` on the UFF_1-tier structure. Close, not identical —
the sweep proxy also flagged `PEDPEW` and `PUVWEK`, which the gate reproduces.

⚠ This does not reopen the A/B or the sweep: both *measured* every molecule. It says a mover list is
a **sample frame**, good enough to choose an A/B population and not good enough to decide which
golden rows are correct.

## 4. Round-trip status of the 51 (ARM 2's predicate)

| old golden → re-frozen | rows |
|---|---|
| `in == out` → `in == out` | 18 |
| `in ≠ out` → **`in == out`** | **21** |
| `in ≠ out` → `in ≠ out` | 5 |
| `in == out` → `in ≠ out` | **0** |
| budget rows (`NO_STRUCTURE@300s` kept) | 7 |

⚠ **ARM 2's predicate is circular** — it scores the generator's own bond graph, the path
`OIN_INDEP_SCORE` replaced in v0.4.8 at a 9.6% false-positive rate. This table says the byte-identity
gate lost nothing; it is not an accuracy reading. That is the sweep: 82.72% / 74.60% verified.

4 of the 7 budget rows (`ECIGAZ PODZEO TOGKOR XENZAS`) assemble a structure on this box today —
with the lever **on and off**. A fact about the box; the sentinel stays.

## 5. 🔴 Gate finding: a killed budget-row reads as a string mismatch

`EQEROI` and `MUKGUW` were SIGKILLed at `--hard-timeout` (450 s) in **both** arms. The gate
synthesises a `HARD_TIMEOUT@450s` row and then compares that token with the golden's `sha_in` — a
**field-2 MISMATCH for a string it never saw**; the kill landed in the generator, after a healthy
encode. Encode-only, both strings are unchanged under both lever settings, and both goldens already
carry `NO_STRUCTURE@300s`, so **nothing is owed and nothing was re-frozen**. But the full gate
**fails on these two rows on a loaded box**, whatever the golden says.

Not fixed here — it changes the gate's contract, which is the owner's call. The one-line repair is
for `gate_arm2_roundtrip_one.py` to flush `sha_in` *before* it starts generating, so a kill leaves
the half of the row that was computed; the gate would then report these as sha1-only rows that still
produce no structure, which is what they are.

## 6. Which lever left the field-2 rows stale: `OIN_CANONICAL_DONOR_FOLD`, i.e. v0.4.13

Settled by intervention (`tools/v0417/arm2_stale_cause.py`): base = the pre-lane encoder
(`OIN_EXACT_DONOR_FOLD=0`), then each of the 11 other default-ON levers set to `"0"` **alone**, plus
the coupled pair donor-fold + parity-veto; fresh process per arm; 14 molecules × 13 arms.

| switching this off brings the golden's hash back | molecules |
|---|---|
| `OIN_CANONICAL_DONOR_FOLD` (alone, and with the veto) | **14 of 14** |
| …and also `OIN_CANONICAL_SLOTS` | 1 (`KELQOI`) |
| …and also `OIN_BORON_CAGE` | 1 (`PAKWAB`) |
| nothing restores it | 0 |

The 14 are the 13 distinct field-2-stale molecules (`JUXPID` sits in both goldens) plus `NASZOY`.
*Control:* the base arm reproduced the audit's lever-off hash on 14/14 — a worker importing the
wrong package would have printed "nothing restores it" fourteen times, and aborted here instead.

So these rows have been wrong since **v0.4.13**, which promoted the donor fold on an *offline
re-score* — no sweep ran with the fold on — and re-froze 11 of 325 rows from a simulated mover list.
`OIN_RESONANCE_DONOR_FOLD` (v0.4.14) restores none of them. `KIHHUG` (field 3) is not in this
trace: its input string was never stale, and what changed its generated structure's encoding was
not traced.

## 7. The real gate against the re-frozen goldens

`run_arm2_refreeze.sh verify` — commit `1a3a1075`, shipped defaults, the 6 + 47 affected rows,
14 min wall. This is `tools/gate_v047.sh arm2` itself reading the new files, not this lane's
re-implementation of its comparison.

| | v047 | v049 |
|---|---|---|
| re-frozen rows fully compared, **PASS** | 6 | 38 |
| re-frozen budget rows, sha1-only, **PASS** | — | 7 |
| **MISMATCH** | 0 | 2 — `EQEROI`, `MUKGUW` (§5: killed, never re-frozen) |

**51 of 51 re-frozen rows pass the real gate.** All 53 rows carry the same `(sha_in, sha_out)` as in
the first shipped run — a third independent run agreeing byte-for-byte.

⚠ Coverage, stated plainly: after the splice the real gate has read the 51 changed rows. The other
374 rows (368 gated and reproduced, 4 observe-only, 2 killed) were read by the real gate *before*
the splice, in the full run, and their golden text is byte-identical (`changed rows == owed rows`,
checked). No single post-splice invocation covered all
425; the next release's full run will be the first.

## 8. What was built

| | |
|---|---|
| `tools/v0417/arm2_field2_audit.py` | step 1; `--reclassify` re-derives classes from stored shas |
| `tools/v0417/run_arm2_refreeze.sh` | `on` / `off` / `verify`, detached, 6 shards, refuses a dirty tree, an `OIN_*` in the user manager, and an existing result directory |
| `tools/v0417/arm2_stale_cause.py` | §6: one lever off at a time over the stale rows; aborts if the base arm does not reproduce the audit |
| `tools/v0417/arm2_refreeze.py` | `diff` (the gate's comparison; field 2 must agree with the audit; builds the control cohorts) and `splice` (fields 1–6 fresh, 7+ preserved, sentinels kept, layout kept, manifest recomputed; aborts on a lost structure or a `NO_STRUCTURE_DET` row) |
| `src/oinsmiles/oin/levers.py` | comment-only: the v0.4.14 block's "field 3 is a fresh stochastic generation" marked WRONG (the goldens' header retracted it on 2026-07-28; this file never heard); the v0.4.17 block records the debt as paid |

Results: `results-v0.4.17-exactfold/arm2_refreeze/` (`on/`, `off/`, `verify/`, the two audits,
`arm2_refreeze_rows.tsv`, the `.BEFORE.tsv` goldens).

## 9. Not done / caveats

- **Field 8** of a v0.4.9 golden (the v0.4.8 honest-class transition, e.g. `key->FAIL`) is preserved
  and is now out of date for most re-frozen rows. Nothing compares it; v0.4.13/14 did not update it
  either. Rebuilding it means rebuilding the golden from the v0.4.17 sweep — a different job.
- The `STALE` reason says *the pre-v0.4.17 encoder on this box does not reproduce the row*. For
  field 2 that is a code fact with a named cause (§6). For `KIHHUG` (field 3) it is not traced.
- v0.4.13's and v0.4.14's own re-frozen rows were not re-audited for *why* they were chosen; only
  that the full gate now passes.
- The 2 killed rows were not re-run alone on an idle box.

## 10. Traps

- **A mover list is a sample frame, not a verdict.** It missed 37% of what this re-freeze owed.
- **`sys.path.insert` beats `PYTHONPATH`** — `gate_arm2_roundtrip_one.py` inserts its own `../src`,
  so the copy of the gate you run decides whose code is tested. Run the worktree's.
- **`chmod +x` dirties a tree.** The launcher's dirty-tree refusal caught it.
- **Re-running `splice --write` on a written golden duplicates the comment block** — `git checkout`
  the goldens first.
- A comment in `src/` outlives its retraction elsewhere. The "stochastic field 3" line shaped this
  lane's first plan (an "out-hash may be noise" hedge that was simply false).
