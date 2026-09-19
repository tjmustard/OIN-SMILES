---
description: "Re-freeze the ARM 2 gate goldens after a shipped DEFAULT changes — from a FULL gate run, one recorded reason per row"
---

Re-freeze `tools/gate_v047_arm2_golden.tsv` and `tools/gate_v049_arm2_golden.tsv`.
Self-contained — do not read any other skill or handoff to run this.

## The rule (owner decision, 2026-09-19)

**Whenever a lever is promoted to default-ON — or any change moves a shipped string — the ARM 2
goldens are re-frozen from a FULL run of the gate: every row of both files. Never from a list of
rows somebody predicted would move.**

## Why this exists

v0.4.13 and v0.4.14 re-froze only the rows a sweep (or an offline re-score) predicted. v0.4.17 ran
the whole gate instead and found:

- **15 of 425 rows had been wrong since v0.4.13** — for two months. ARM 2 is only run in full at a
  release (`--band fast` is the routine control), so a red full gate is invisible between releases.
- **Its own 36-row mover list would have missed 19 of the 51 rows it owed (37%)** and "re-frozen"
  4 that had not moved. 8 of the missed rows move on field 3 alone: the input string is unchanged,
  the seeded generator builds the same structure, and it is the *re-encode of that structure* the
  lever relabels. No input-side prediction can see that.

A mover list is a **sample frame** — good enough to choose an A/B population, not to decide which
golden rows are correct. The full run costs **75 minutes of wall time**, 6 shards.

Record: `docs/agentic-notes/v0.4.17/ARM2_REFREEZE.md`. Evidence: `measurements/v0.4.17/v0417_arm2_*`.

## Argument

`$ARGUMENTS` is the lever being promoted (e.g. `OIN_EXACT_DONOR_FOLD`) and the release
(e.g. `v0.4.18`). If either is missing, ask. If the change is not a lever (no `=0` spelling brings
the old behaviour back), say so and stop: the control in step 4 needs the **previous commit's**
checkout as its `off` arm, and that is a judgement call for the owner, not a default.

## What makes this possible: ARM 2 is DETERMINISTIC

Generation is seeded (`MetalloGenAdapter seed=42`). v0.4.17 reproduced all 368 untouched gated rows
six shards wide, and a third run matched the first on 53/53. So a moved row **is attributable**:
re-run it with the lever at `"0"`. If the old row comes back, the lever moved it. If it does not,
the golden was already stale and the recorded reason must say so.

⚠ `levers.py` and `CHANGELOG.md` (v0.4.14) once said field 3 was "a fresh stochastic generation".
That was retracted in the goldens' own header and is false. Do not plan around it.

## Steps

Work in the lane's **worktree**, with the MAIN checkout's python (rdkit is pinned there; never
`uv sync` in a worktree). Nothing else CPU-heavy may run beside steps 2, 4 and 6 — **load biases
the budget-dependent rows**.

```bash
V=<main checkout>/.venv/bin/python
D=<main checkout>/tmCAT-tmPHOTO_xyz_dataset
export REFREEZE_LEVER=<OIN_THE_LEVER>
export REFREEZE_TAG=<v0.4.NN>
export REFREEZE_OUT=$D/results-$REFREEZE_TAG-arm2-refreeze     # inside the dataset tree, never /tmp
A="--out-dir $REFREEZE_OUT --cohort-tag $REFREEZE_TAG"
mkdir -p $REFREEZE_OUT
```

0. **ARM 1 first — it is cheap and the gate verifies its manifest.**
   `bash tools/gate_v047.sh arm1`, then again with `$REFREEZE_LEVER=0` in the environment. The
   `=0` run must reproduce the old manifest byte-for-byte; re-freeze the rows the shipped run
   moves and write the reason into the golden's header.

1. **Audit field 2 of EVERY row, encode-only (~15 min).**

   ```bash
   for g in v047:cohort-v047-slow100 v049:cohort-v049-strata; do
     $V tools/v0417/arm2_field2_audit.py --lever $REFREEZE_LEVER \
        --golden tools/gate_${g%%:*}_arm2_golden.tsv --cohort-dir $D/${g#*:} \
        --out $REFREEZE_OUT/field2_audit_${g%%:*}.jsonl | tee $REFREEZE_OUT/field2_audit_${g%%:*}.txt
   done
   ```

   Read the classes. **`STALE` / `STALE_AND_MOVED` / `HEALED` mean the golden was wrong before you
   started** — a finding in its own right. The tool aborts if a worker imported a foreign
   `oinsmiles`, and if the lever fired on 0 rows (a dead lever and a clean audit print the same).

2. **Run the FULL gate, shipped defaults.** Commit first — the launcher refuses a dirty tree, and a
   `chmod +x` dirties it.

   ```bash
   bash tools/v0417/run_arm2_refreeze.sh on        # 6 detached units; done when DONE_1..6 exist
   ```

3. **`$V tools/v0417/arm2_refreeze.py $A diff`** — the gate's own comparison over every golden row.
   It **aborts unless the gate's field-2 mismatches are exactly the audit's**; if it aborts, look
   at the rows it names before trusting either. It builds the control cohorts (the unreproduced
   rows) and separates `KILLED` rows — see *Things that will catch you*.

4. **The control:** `bash tools/v0417/run_arm2_refreeze.sh off` — the same rows with the lever at `0`.

5. **Splice.** Dry run, read it, then write:

   ```bash
   $V tools/v0417/arm2_refreeze.py $A splice
   $V tools/v0417/arm2_refreeze.py $A splice --write --comment-file <block.txt>
   ```

   Every re-frozen row gets exactly one reason — `LEVER`, `STALE`, or `STALE+LEVER` — written
   into the golden. The comment block is yours to write: say how the run was done, what the counts
   were, and how many rows changed round-trip status (**and that ARM 2's predicate is circular —
   byte identity, not accuracy**). Commit the goldens.

6. **Let the REAL gate read them:** `bash tools/v0417/run_arm2_refreeze.sh verify`. Every
   re-frozen row must pass. Anything else that fails is either a `KILLED` row or a real problem.

7. **If step 1 found stale rows, name the cause:**
   `$V tools/v0417/arm2_stale_cause.py --base-lever $REFREEZE_LEVER --out-dir $REFREEZE_OUT --out $REFREEZE_OUT/stale_cause.jsonl`
   — each other default-ON lever set to `0` alone; the one that brings the golden's hash back
   dates the staleness.

8. **Freeze the evidence** with `/freeze-measurements` (`arm2_refreeze.py $A rows` merges the shards
   first). ⚠ The harvester rebuilds a release's `README.md` index from the **current picks only** —
   stage these files *with* the release's other files and re-harvest the whole set.

9. **Report**, in the lane note and to the owner: rows re-frozen per golden, the split by reason,
   predicted-vs-measured if a mover list existed, round-trip transitions, and the coverage of the
   final verify run stated plainly (it reads the changed rows; the untouched rows were read by the
   full run before the splice).

## Things that will catch you

- **`sys.path.insert` beats `PYTHONPATH`.** `gate_arm2_roundtrip_one.py` inserts its own `../src`,
  so *which copy of the gate you run* decides whose code is tested. Run the worktree's.
- **A SIGKILLed row reads as a string mismatch.** The gate synthesises `HARD_TIMEOUT@450s` and
  compares that token with `sha_in`. `EQEROI` and `MUKGUW` die in the generator in every run, lever
  on and off. `diff` files them as `KILLED`, settles field 2 from the audit, and re-freezes nothing;
  they go in the control anyway — **if `off` is NOT killed too, the lever is involved**, and
  `splice` aborts. (Open with the owner: have the runner flush `sha_in` before it generates.)
- **Never splice whole rows.** Fields 1–6 come from the fresh run; **fields 7+ are preserved** —
  field 7 of a v0.4.9 golden is the band `--band` filters on, and a fresh row carries `xyz_sha`
  there. A `NO_STRUCTURE@300s` sentinel is **kept** and only `sha_in`/`len_in` move: whether this
  box assembles a structure inside the budget today is a fact about the box.
- `splice` **aborts** on a row that lost its structure and on any `NO_STRUCTURE_DET` row. That is a
  judgement — budget or defect — and a person makes it.
- **`splice --write` twice duplicates the comment block.** `git checkout` the goldens first.
- `# MANIFEST_SHA256` is recomputed by `splice`; ARM 2 never verifies it, so a stale one goes unseen.
- Field 8 of the v0.4.9 golden (the v0.4.8 honest-class transition) is a stale observation column.
  Nothing reads it. Do not "fix" it row by row; rebuilding it means rebuilding the golden.
- `--shard` is **1-based** (`1:6` … `6:6`). The launcher handles it; a hand-rolled loop over `0..5`
  silently drops a sixth of the cohort.
