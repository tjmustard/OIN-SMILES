# ▶ START HERE — capstone validation (v0.4.2 round-trip accuracy wave)

**Launch a fresh Claude Code session in the main checkout and hand it this file.** Run **only after
every fix phase (S1, S3, S5, S6a, S6b, S7, docs) has squash-merged into `release/v0.4.2`.** The
capstone **edits no source** — a capstone that edits the thing it measures is not a measurement. If
it finds a regression, it reports + opens a follow-up phase.

### 1 · Work in the staging worktree
```bash
cd /home/tjmustard/Documents/GitHub/OIN-SMILES-v0.4.2   # on release/v0.4.2, all phases landed
git log --oneline release/v0.4.2 ^c7edeeb6              # confirm every phase squash is present
uv run python -m unittest discover tests/unit           # MUST be green BEFORE measuring
```
A mid-validation red unit suite means a phase landed broken — a more urgent finding than any number.

### 2 · Read these
- shared protocol — `spec/handoffs/v0.4.2/README.md`; floor — `spec/handoffs/v0.4.2/BASELINE.md`
- precedent — `spec/handoffs/v0.4.0-perf/VALIDATION.md`, `INTEGRATION.md`

### 3 · The gate — per-molecule set-inclusion, NOT a percentage
1. **Pause the accumulator** (only one sweep at a time). Run one **seeded** sweep of the wave's target
   set + the must-not-regress passing set from `BASELINE.md`, on `release/v0.4.2`, non-`--quick`,
   `--mol-timeout 1800`, into a **private** dir. Resume the accumulator after.
2. **Blocker gate:** `{passes on release/v0.4.2} ⊇ {passes on c7edeeb6}`. Any molecule that passed on
   `c7edeeb6` (trusted-provenance set) and now fails is a **named blocker with root cause** —
   regardless of how many the wave fixed. A wave that fixes 400 and regresses 3 does not ship until
   the 3 are explained/fixed (open a follow-up phase).
3. **Claim specific named molecule flips** in both directions, **per-class AND per-CN AND per-metal**
   (a CN4 fix that regresses CN5 nets to zero and hides in a global number — that is exactly how a
   neutral-looking change masks a regression). No headline-delta claims (noise below ~2pp ≈ 52 mol).
4. **Composition re-check on the integrated tree.** Phases interact: S1's `AddHs` shifts the index
   frame S6a/S6b capture; S5's `scored` change flows into the winding pick; S7's pool/FF changes alter
   which conformer S5 selects. Re-verify **every phase's per-class goldens** on the integrated tree,
   not just in isolation.
5. **rdkit discipline:** run all `smiles_1` diffs under **one** pinned rdkit (declare which is
   authoritative). For any pure-refactor hunk, byte-identity to pristine at fixed seed is the proof.

### 4 · Write `spec/handoffs/v0.4.2/VALIDATION.md`
- The sweep command, commit (`release/v0.4.2` head), rdkit, seed, date.
- Per-class before→after named flips; per-CN and per-metal breakdown.
- The blocker check result: the (hopefully empty) list of `c7edeeb6`-passers that regressed, each with
  root cause or a follow-up phase reference.
- Explicit GO / NO-GO verdict for the merge to `main`.

### 5 · Merge to main (only on GO, and on the user's go-ahead)
- Confirm no other wave targets `main` concurrently (if one does, do an `INTEGRATION.md`-style
  reconciliation first — v0.4.0 needed it because a naive squash would have reverted later main
  commits; for v0.4.2 branched from current `main` with nothing else in flight this risk is low but
  **confirm, don't assume**).
- One squash: `git -C /home/tjmustard/Documents/GitHub/OIN-SMILES merge --squash release/v0.4.2` then
  commit (subject `v0.4.2: round-trip accuracy wave (...)`). Leave **unpushed** unless the user says
  otherwise. Bump `pyproject` to 0.4.2 + tag `v0.4.2` on the user's go-ahead.
- Post-land: tag `archive/roundtrip-<slug>` for any not-yet-tagged phase branch before cleanup;
  `spec/handoffs/v0.4.2/` stays gitignored on `main` (tracked on the release branch only).
