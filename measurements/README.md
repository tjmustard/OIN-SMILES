# `measurements/` — the numbers, frozen so a later release can diff them

**Written by `tools/harvest_measurements.py`. Do not hand-edit — rerun the tool.**

Every sweep this project has ever run lives in `tmCAT-tmPHOTO_xyz_dataset/`, which is
**gitignored in its entirety** (`.gitignore:85`) — 2.1 GB of dataset plus ~20 `results-*/`
directories, none of which survive a clone. `spec/handoffs/` is gitignored too. So for six
consecutive releases the project wrote careful prose about numbers that existed on exactly one
disk. Then:

> 🔴 **v0.4.11's mirror-audit JSON is gone.** The 19/250 measurement that chartered the whole of
> v0.4.12 exists nowhere. It was recoverable only because the *tool* was committed and
> deterministic — at ~50 minutes of CPU to re-run. The next one might not be.

This tree is the fix. It holds the small, derived artifacts a **later release actually diffs**.

## What goes here, and what does not

| data | goes in | the test |
|---|---|---|
| **Cross-release comparison** — frozen baselines, gate/audit tallies, transition matrices, A/B rows | **`measurements/<release>/`** | *"Will a future release diff this to know whether it regressed?"* Operationally: **if a CLOSEOUT step produced it, it belongs here.** |
| **Report-backing** — the evidence for one paragraph of one write-up | `docs/agentic-notes/<release>/` | *"Is this read alongside one specific note?"* |

**Never here**, and the tool enforces it: `individual_reports/` (~20 MB per sweep),
`structures/`, the `cat/`+`photo/` inputs, `*.log` (v0.4.12's logs alone were 29 MB of rdkit
noise), and anything over the 512 KB per-file cap. A harvest totalling more than 5 MB refuses
outright. This tree is a comparison index, not a data dump.

## ⚠ This tree is PUBLIC

`origin` is a public GitHub repo and `main` is pushed regularly. Everything here is world-
readable and permanent. Two consequences:

1. **Off-machine durability comes free** — the data is backed up by virtue of being pushed. It
   is *not* private, and nothing here should ever be treated as if it were.
2. The harvester **scrubs local absolute paths** on the way in (`<repo>/` stripped, `/home/<user>/`
   → `<HOME>/`, scratchpad paths → `<SCRATCH>`) and drops any file where a path survives
   scrubbing. That guard covers the index files the tool generates as well, because the first
   run leaked the full scratchpad path into exactly one of them.

## ⚠ Reading a tally: an absent key means ZERO, not "not measured"

`mirror_audit_donor_fold.py` omits a verdict key entirely when its count is zero. So a clean
run reads:

```python
{'achiral_or_preexisting_fold': 157, 'distinct_both_arms': 92, 'encode_failed': 1}
#  ^ no REGRESSION_raw_collapsed key at all -- that IS the pass
tally.get('REGRESSION_raw_collapsed')      # -> None, NOT 0
```

**Use `tally.get(k, 0)`.** This project has been bitten repeatedly by the same shape — *"0 fail"*
and *"0 fail over 0 measured"* printing identically — and a `None` here reads as *"the audit
didn't run"* when it actually means *"the audit was perfectly clean"*.

## Known gaps — real history, not oversight

| release | state |
|---|---|
| v0.4.0 – v0.4.8 | backfilled by `--backfill` from the surviving `results-*/` dirs |
| **v0.4.9, v0.4.10, v0.4.11** | **no artifacts preserved.** v0.4.10/v0.4.11 ran no sweep (byte-identical by construction / carry-forward licence), and **v0.4.11's mirror-audit JSON is confirmed lost** |
| v0.4.12 | complete — harvested from the session scratchpad before it expired |
| v0.4.13-honest | present, harvested from a concurrently-running release; may be superseded |

`results-v0.4.9-bound/*.jsonl` still exists on disk and is recoverable if anyone wants it. It is
not harvested because it is per-attempt telemetry rather than a comparison table.

## Usage

```bash
# manual, and always dry-run first
python tools/harvest_measurements.py --release v0.4.13 --from <results dir> --dry-run
python tools/harvest_measurements.py --backfill --dry-run

# from a worktree: the destination is still the MAIN checkout, by design
python tools/harvest_measurements.py --link      # optional read-access symlink
```

Or `/freeze-measurements` for the guided path. It runs automatically from `CLOSEOUT.md` §4b.

**The tool always writes to the main checkout**, resolved via `git rev-parse --git-common-dir`,
never the current worktree — worktree files vanish on `git worktree remove`, and untracked files
do not survive `git clean -fd`. That is precisely how the data this tree exists to protect was
nearly lost twice.
