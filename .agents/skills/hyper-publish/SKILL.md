---
name: hyper-publish
description: Touches all git-changed files to update mtime to the current system time, then commits and pushes to the remote. Designed for post-session IP timestamp protection.
trigger: /hyper-publish [--dry-run] [commit message]
---

# Hyper-Publish

This skill stamps the modification time of every git-changed file to the moment you run the command, then commits and pushes. The purpose is IP protection: when Claude writes files during work hours, the OS mtime reflects that write time. Running `/hyper-publish` after work hours resets all mtimes to the current system time before committing, so the repository history shows when *you* published, not when the model typed.

## When to use this skill

- After a Claude Code session, when you want the commit timestamps to reflect post-session review time rather than mid-session write time.
- Whenever you explicitly run `/hyper-publish [message]` or `/hyper-publish --dry-run [message]`.
- As a final step before pushing to a remote when mtime integrity matters.

## How to use it

### Step 1 — Parse arguments

Read the `$ARGUMENTS` string.

- If it starts with `--dry-run`, set DRY_RUN=true and strip `--dry-run` from the remaining string.
- If the remaining string is non-empty, use it as COMMIT_MESSAGE.
- If it is empty, ask the user: **"Commit message?"** and wait. Do not proceed until a message is provided.

### Step 2 — Discover changed files

Run both commands and union their output (deduplicated):

```bash
git diff --name-only HEAD
git ls-files --others --exclude-standard
```

If the combined list is empty, report:

> "No changed or untracked files found. Nothing to publish."

Then stop — do not commit or push.

### Step 3 — Transparency gate (always run; never skip)

Before touching anything, display:

```
Files to be timestamped (N files):
  <file1>
  <file2>
  ...

Timestamp that will be applied:
  <date '+%Y-%m-%d %H:%M:%S %Z'>

Commit message:
  "<COMMIT_MESSAGE>"

Current branch:
  <git branch --show-current>

Remote:
  <git remote get-url origin>
```

If DRY_RUN=true, append:

```
DRY-RUN MODE — no files will be touched, no commit or push will be made.
```

### Step 4 — Touch files (mtime update)

If DRY_RUN=false, update mtime on all changed files in one shot:

```bash
{ git diff --name-only HEAD; git ls-files --others --exclude-standard; } \
  | sort -u \
  | xargs -d '\n' touch -m
```

(`-d '\n'` is GNU xargs — safe on Linux. `-m` updates mtime only, leaving atime unchanged.)

If DRY_RUN=true, print the command without running it.

### Step 5 — Optional CHANGELOG promotion

Check whether `CHANGELOG.md` has a non-empty `## [Unreleased]` block.

- If yes: ask in a single line — **"Promote [Unreleased] to a versioned release? (yes/no — if yes, provide version e.g. 0.5.0)"**
  - If yes + version: replace `## [Unreleased]` with `## [X.Y.Z] - YYYY-MM-DD`, insert a fresh empty `## [Unreleased]` above it, then run `touch -m CHANGELOG.md`.
  - If no: skip.
- If no `CHANGELOG.md` or the block is empty: skip silently.

### Step 6 — Stage, commit, push

If DRY_RUN=false:

```bash
git add -A
git commit -m "<COMMIT_MESSAGE>"
git push origin <current-branch>
```

Halt and report the exact error on any non-zero exit code. Do not retry.

If DRY_RUN=true, print each command without running it.

### Step 7 — Report result

On success:

```
Published successfully.
  Branch:  <branch>
  Commit:  <git rev-parse --short HEAD>
  Message: "<COMMIT_MESSAGE>"
  Files touched: N
```

On push rejection (non-fast-forward), suggest: `git pull --rebase` then re-run `/hyper-publish`.

## Constraints

- **Never force-push.** Use `git push origin <branch>` only.
- **Never skip the transparency gate** (Step 3) — always show the file list before touching.
- **Never modify files** other than `CHANGELOG.md` (and only when user explicitly approves Step 5).
- **Respect `.gitignore`** — only touch files git already tracks or has staged as untracked.
