# v0.4.2 — session launch & landing prompts

Copy-paste prompts for driving the wave. All paths are the main checkout
`/home/tjmustard/Documents/GitHub/OIN-SMILES`.

## The phases

| Phase | Open this file | Branch / worktree | Gate |
|---|---|---|---|
| P0 | `spec/handoffs/v0.4.2/P0-baseline.md` | `feature/roundtrip-baseline` | clean floor set + per-class goldens committed; no headline % |
| S1 | `spec/handoffs/v0.4.2/S1-donor-h.md` | `feature/roundtrip-donor-h` | fixable donor-H/oxo-imido subset round-trips; notation residual → docs |
| S3 | `spec/handoffs/v0.4.2/S3-aromatic.md` | `feature/roundtrip-aromatic` | the 4 aromatic-perception classes encode; 0%→pass on goldens |
| S5 | `spec/handoffs/v0.4.2/S5-geometry.md` | `feature/roundtrip-geometry` | geometry-change + NON + winding goldens round-trip; no CN regressions |
| S6a | `spec/handoffs/v0.4.2/S6a-ez-stereo.md` | `feature/roundtrip-ez-stereo` | free-arm C=N E/Z round-trips; AFECIZ/XIZXAG unblocked |
| S6b | `spec/handoffs/v0.4.2/S6b-atom-stereo.md` | `feature/roundtrip-atom-stereo` | atom_stereo + `[S@SP3]` goldens round-trip |
| S7 | `spec/handoffs/v0.4.2/S7-metrics.md` | `feature/roundtrip-metrics` | ENABLED_METALS RMSD win; honest triage of high_rmsd/timeout/no_conformers |
| docs | `spec/handoffs/v0.4.2/docs.md` | `feature/roundtrip-docs` | carborane + notation/artifact residuals documented |
| capstone | `spec/handoffs/v0.4.2/VALIDATION.md` | `feature/roundtrip-validation` | per-molecule set-inclusion gate; then one squash → main |

## To launch a phase

Open a fresh Claude Code session **in the main checkout** and hand it the doc:

```
@spec/handoffs/v0.4.2/S1-donor-h.md
```

The doc's `▶ START HERE` bootstrap creates the worktree and proceeds. Launch order:
**P0 first**, then S3 / docs / S5 / S6b / S7 concurrently, and the serial chain **S1 → S6a**.

## To land a finished phase (paste into the finished phase session)

> Your phase is complete and green. Land it into the integration branch **`release/v0.4.2`** by
> squash-merge — do **not** open a PR to `main`, and do **not** push. Steps:
>
> 1. In your worktree `../OIN-SMILES-<slug>`, commit everything on `feature/roundtrip-<slug>`.
>    Then `git rebase release/v0.4.2` (absorbs anything landed since you branched) and re-run your
>    gate **and** `uv run python -m unittest discover tests/unit` to confirm still-green.
> 2. Squash-merge from the staging worktree:
>    ```bash
>    cd ../OIN-SMILES-v0.4.2
>    git merge --squash feature/roundtrip-<slug>
>    git commit    # subject: "v0.4.2-<phase> <slug>: <one-line summary>"
>                  # body: named molecule flips (per-class) + regression spot-check result
>    ```
> 3. Confirm content landed: `git diff --quiet release/v0.4.2..feature/roundtrip-<slug> -- <owned files>` exits 0.
> 4. Tag + clean up: `git tag archive/roundtrip-<slug> feature/roundtrip-<slug>`, then
>    `git worktree remove ../OIN-SMILES-<slug>` and `git branch -D feature/roundtrip-<slug>`.
> 5. Announce to live phases that they should `git rebase release/v0.4.2`.
>
> If the rebase conflicts with another phase's squash, resolve it **in your own functions only**; if
> it reaches a function you do not own, stop and report rather than editing it.

## Pausing / resuming the live accumulator (before any large sweep)

P0 and the capstone run a **large** sweep and MUST pause the `--quick` accumulator first (only one
sweep at a time, across all worktrees):

```bash
# find it
pgrep -af "test_dataset_roundtrip.*--continue"
# pause: stop the --continue loop (kill the wrapper/loop, not mid-molecule if avoidable)
# ... run the large sweep into a PRIVATE --output-dir ...
# resume: restart the accumulator loop on results-v0.4.0 exactly as it was
uv run python tools/test_dataset_roundtrip.py \
  --dataset-dir tmCAT-tmPHOTO_xyz_dataset --quick --output-dir tmCAT-tmPHOTO_xyz_dataset/results-v0.4.0 \
  --limit 1 --random --continue --mol-timeout 30   # (its original invocation)
```

Per-phase checks do **not** pause it — they use `--only <goldens>` into `/tmp/rt-<slug>`.

## Final merge to main (capstone only, after sign-off)

> The capstone's per-molecule gate is green on `release/v0.4.2` (no molecule that passed on
> `c7edeeb6` regressed; named fixes listed). Merge the wave to `main` as **one** squash:
> ```bash
> cd /home/tjmustard/Documents/GitHub/OIN-SMILES   # main checkout, on main
> git merge --squash release/v0.4.2
> git commit    # subject: "v0.4.2: round-trip accuracy wave (S1/S3/S5/S6a/S6b/S7 + docs)"
> ```
> Leave **unpushed** unless the user says otherwise. Then bump `pyproject` to 0.4.2 and tag `v0.4.2`
> per project release convention, on the user's go-ahead.
