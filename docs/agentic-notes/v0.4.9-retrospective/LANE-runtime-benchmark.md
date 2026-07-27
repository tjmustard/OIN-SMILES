# Lane 2 retrospective — freeze the runtime benchmark

Findings and provenance: `docs/agentic-notes/v0.4.9/RUNTIME_BENCHMARK_v0.4.9.md`.
This records what was *tried*, what was *abandoned*, and what the lane got wrong on the way.

## What the charter asked for vs what shipped

| charter | shipped | why the difference |
|---|---|---|
| "re-derive from the finished 5k sweep" | drew from `results-v0.4.8-honest` | it **is** the 5k sweep plus the honest fields — the only dir carrying both `elapsed_s` and `honest_class`. Drawing from `results-v0.4.6-sweep` would have thrown away the labelling that makes the cohort honest. |
| "extend `gate_v047_build_arm2_golden.py`; do not write a second one" | extended it (`--predicate`, `--strata-file`, `--allow-unbaselined`) **and** wrote a second *selector* | the golden builder was genuinely reusable. The **selector** was not: `select_slow_byte_exact.py` is top-N-by-elapsed under a byte-exact predicate, and stratified sampling with no predicate is a different algorithm. Bolting it on would have made the v0.4.7 cohort non-re-derivable. |
| `docs/RUNTIME_BENCHMARK_v0.4.9.md` | `docs/agentic-notes/v0.4.9/…` | `docs/` root is closed to four product docs and a `pre-commit` guard enforces it. The charter's path would have been rejected at commit time. |

## Negative results and course corrections

**Integer stride under-delivered and cut the top off every cell.** The first draw used
`cell[::ceil(n/quota)]`. On the 151-member `b30_60/eta` cell that gives `k=4`, **38** picks
instead of 40, and a last pick at rank 148 — losing the three slowest molecules in a band, in a
benchmark whose entire purpose is a **p100** target. Replaced with even rank spacing
(`round(i·(n−1)/(quota−1))`), which hits both endpoints and the quota exactly. *A sampling rule
that silently drops the extreme is worse than no stratification, because it looks stratified.*

**Three pins are not in the source sweep at all.** `GOHWOQ`, `ULODUU` and `XIQKOY` — the
advisory-overrun fixture and two of the three boron fixtures — are absent from the frozen
seed-42 5k draw. The first instinct was to drop them from the cohort. That is backwards: it
would remove the release's own named evidence from the release's own benchmark. They are in the
cohort dir, measured, and carry no frozen golden, with `--allow-unbaselined` saying so out loud.

**`NO_STRUCTURE` was initially applied to all 50 no-structure rows.** Splitting them by the
source report's error text showed **45 timed out and 5 raised**. Only the timeouts are
budget-dependent; the 5 that raised are code-determined and are now gated strictly
(`NO_STRUCTURE_DET`), so a lane that accidentally makes `GURKUA` generate is caught rather than
shrugged at. *Weakening a gate uniformly, because part of it must be weak, throws away the part
that could have stayed strong.*

**Four rows carry no byte-identity signal at all.** `LUJJUX`, `NOGWOX`, `NOXREZ`, `UROFUD` were
hard-killed **inside the encoder** (`exceeded 300s while encoding`). Corpus-wide that is 15
molecules and 1.25 CPU-h, and it is exactly the `encode_fail` bucket. **No generator-side bound
can reach them** — they never get as far as the generator. That is a limit on Lane 1's claim,
found by Lane 2, and it is an argument for the per-molecule bound over a generator-only one.

## Two defects found in the tool doing the measuring

Both were shipped in v0.4.7 and both would have corrupted this lane's own results.

1. **The gate silently ran on an arbitrary interpreter.** Its fallback globbed
   `$(dirname $REPO)/*/.venv/bin/python`. From a worktree it selected `EtaCatalysis/.venv` (no
   rdkit → every molecule died → short `#DONE` → reads as "truncated run"), then
   `EtaTMCSMILES/.venv` (rdkit **2025.09.2** vs the pinned **2025.9.3**). The second is the
   dangerous one: **a byte-identity gate on a different rdkit reports MISMATCHes that look like
   code regressions.** An `import rdkit` check was not enough — the fix is to resolve the main
   checkout deterministically via `git rev-parse --git-common-dir` and refuse on version drift.
2. **The gate itself was unbounded.** `--timeout` is only the generator's *advisory* budget — the
   very defect Lane 1 exists to fix, sitting inside the tool measuring it. `--hard-timeout` now
   applies a real `SIGKILL` and synthesizes a `HARD_TIMEOUT@Ns` row, so a kill is a **result**
   rather than a missing line that makes `#DONE` short.

## What this lane deliberately did not do

- **Did not re-run the sweep.** The goldens come from the existing primary source; re-generating
  100+ slow molecules to freeze values already recorded would have cost the exact wall-clock the
  cohort was selected to be expensive at.
- **Did not gate correctness.** 71 of 325 rows are known `byte->FAIL`. ARM 2 detects *change*;
  `roundtrip_bucket_report.py --score honest` detects *correctness*. Conflating them is how the
  v0.4.7 cohort came to be read as a quality statement.
- **Did not touch ARM 1.** 61 fixtures, same golden, same verdict.
