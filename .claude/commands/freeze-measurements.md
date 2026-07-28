---
description: "Freeze a release's comparison numbers into the tracked measurements/ tree before they expire"
---

Freeze measurement artifacts into `measurements/`. Self-contained — do not read any other skill
or handoff to run this.

## Why this exists

`tmCAT-tmPHOTO_xyz_dataset/` is **gitignored in its entirety** and so is `spec/handoffs/`, so a
release's numbers live on exactly one disk until they are harvested. **v0.4.11's mirror-audit
JSON — the 19/250 measurement that chartered all of v0.4.12 — is gone.** It was recoverable only
because the tool was committed and deterministic, at ~50 minutes of CPU.

Session scratchpads under `/tmp` are worse: they expire. v0.4.12's *entire* output lived there.

## Argument

`$ARGUMENTS` is a release (`v0.4.13`), or `--backfill` to sweep the historical `results-*` dirs.
If empty, ask which release, and ask **where the numbers actually are** — a results directory, a
scratchpad, or both. Do not guess; a wrong `--from` silently harvests nothing.

## Steps

1. **Locate the sources.** Check `tmCAT-tmPHOTO_xyz_dataset/results-<release>*` *and* this
   session's scratchpad. Most releases have both, and the scratchpad is the one that expires.

2. **Dry-run first, always.**

   ```bash
   V=$PWD/.venv/bin/python
   $V tools/harvest_measurements.py --release <rel> --from <dir> [--from <dir2>] --dry-run
   ```

   Show the user the selection and the totals. Expect tens of KB to a few hundred KB. If it runs
   to megabytes something is wrong with the source — investigate rather than passing `--force`.

3. **Read the output for two things the tool prints deliberately:**
   - `[provenance UNKNOWN]` — the file will be recorded that way rather than given a
     plausible-looking command. Add a `PROVENANCE` pattern to the tool if you know the real one.
   - `🔴 REFUSED` — a local path survived scrubbing. The tree is public; find out why before
     forcing anything.

4. **Write it**, then verify:

   ```bash
   $V tools/harvest_measurements.py --release <rel> --from <dir>
   du -sh measurements/
   grep -rl "/home/\|/tmp/claude" measurements/ && echo "LEAK" || echo "clean"
   ```

5. **Commit.** End the message with `Claude-Session: <url>` — the `commit-msg` hook only
   *rewrites* that line, it never inserts one, so omitting it means no trailer at all.

## Things that will catch you

- **It writes to the MAIN checkout, never your worktree** — by design, resolved via
  `git rev-parse --git-common-dir`. Worktree files vanish on `git worktree remove` and untracked
  files do not survive `git clean -fd`. Both nearly ate this data already. Use `--link` if you
  want read access from a worktree.
- **`measurements/` is PUBLIC.** `origin` is a public GitHub repo and `main` is pushed regularly.
- **An absent tally key means zero, not "not measured."** `tally.get('REGRESSION_raw_collapsed')`
  returns `None` on a perfectly clean audit. Use `.get(k, 0)`.
- **Do not hand-edit anything under `measurements/`.** Rerun the tool; it regenerates each
  release's `README.md` index with sizes, sha256s and provenance.

Full rule: `measurements/README.md`.
