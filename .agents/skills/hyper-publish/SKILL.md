---
name: hyper-publish
description: "Enhanced commit-push-PR automation: AI-suggests commit messages from CHANGELOG, uses HITL approval gates at each step, optionally creates PRs with auto-drafted descriptions. No system-level git/gh approval prompts required."
trigger: /hyper-publish [--dry-run] [--skip-pr] [commit message]
---

# Hyper-Publish: Commit, Push, and PR Creation

This enhanced skill automates the full git commit → push → (optional) PR creation workflow with human-in-the-loop (HITL) approval at each critical step.

**Key improvements over basic `git commit && git push`:**
- AI-proposes commit message from CHANGELOG `[Unreleased]` + git diff (no typing required)
- Timestamps all file mtimes to post-session time
- Displays transparency gate before any changes
- HITL approval gates before each file modification and git command
- Optionally creates a GitHub PR with auto-drafted description
- No system-level approval prompts for git/gh commands (pre-approved via settings.json)

---

## When to use this skill

- After a Claude Code session, when you want to commit changes with an AI-suggested message based on your CHANGELOG
- When you want to commit, push, and optionally create a PR in one unified workflow
- Whenever you run `/hyper-publish [message]` or `/hyper-publish --dry-run`

---

## How to use it

### Step 0 — Initialize Permissions (automatic)

Run `python .agents/skills/hyper-publish/scripts/setup_permissions.py` to patch `.claude/settings.json` with pre-approved patterns for `Bash(git *)` and `Bash(gh *)`. This prevents Claude Code from prompting for shell command approval on every invocation.

If settings.json is already configured, this step skips silently.

---

### Step 1 — Parse Arguments

Read the `$ARGUMENTS` string.

- If it starts with `--dry-run`, set `DRY_RUN=true` and strip `--dry-run` from the remaining string.
- If the remaining string starts with `--skip-pr`, set `SKIP_PR=true` and strip it.
- If remaining string is non-empty, use it as `COMMIT_MESSAGE` (skip Step 3 and proceed to Step 4).
- If remaining string is empty, proceed to Step 3 (generate message suggestion).

**Examples:**
- `/hyper-publish` → Suggest message, offer PR
- `/hyper-publish "Update README"` → Use provided message, offer PR
- `/hyper-publish --dry-run` → Dry-run mode, suggest message
- `/hyper-publish --skip-pr "Release v1.0"` → Use provided message, skip PR step

---

### Step 2 — Discover Changed Files

Run both commands and union their output (deduplicated):

```bash
git diff --name-only HEAD
git ls-files --others --exclude-standard
```

If the combined list is empty, report:

> "No changed or untracked files found. Nothing to publish."

Then stop — do not proceed.

---

### Step 3 — NEW: Propose Commit Message

**Only if no commit message was provided in `$ARGUMENTS`.**

Run the helper script:

```bash
python .agents/skills/hyper-publish/scripts/generate_commit_message.py
```

The script reads:
- `CHANGELOG.md` → extracts `## [Unreleased]` bullet points
- `git log --oneline -5` → infers style (Release vX.Y.Z, feat: prefix, or imperative)
- `git diff --name-only HEAD` → infers scope from changed file paths
- Outputs JSON: `{"subject": "...", "body": "...", "style": "...", "source_bullets": [...]}`

Parse the JSON and use **AskUserQuestion** to present options:

```
Proposed commit message:
  <subject line>

Options:
- Option A: Accept this message
- Option B: Edit it (enter your own in the text field)
- Option C: Enter a completely custom message
```

Bind the result to `COMMIT_MESSAGE` (the variable used in later steps).

**Fallback:** If the script fails or CHANGELOG is empty, prompt the user directly:

> "Commit message (required)?"

---

### Step 4 — Transparency Gate + HITL Approval Gate #1

Before any modifications, display a preview panel showing:

```
Files to be timestamped (N files):
  <file1>
  <file2>
  ...

Timestamp to be applied:
  <date '+%Y-%m-%d %H:%M:%S %Z'>

Commit message:
  "<COMMIT_MESSAGE>"

Current branch:
  <git branch --show-current>

Remote:
  <git remote get-url origin>

PR creation: <Enabled / Skipped (--skip-pr flag)>
```

If `DRY_RUN=true`, append:

```
DRY-RUN MODE — no files will be touched, no commit or push will be made.
```

Then use **AskUserQuestion** to request explicit approval:

```
Review the changes above. Proceed?

- Option A: Proceed — touch files, commit, and push
- Option B: Abort
```

**Do not proceed past this point without explicit user approval.**

---

### Step 5 — Touch Files (mtime update) + HITL Approval Gate #2

If `DRY_RUN=false`, display the exact touch command that will be executed:

```bash
{ git diff --name-only HEAD; git ls-files --others --exclude-standard; } \
  | sort -u \
  | xargs -d '\n' touch -m
```

Then use **AskUserQuestion**:

```
Touch N files to update mtime?

- Option A: Touch files (update mtime to current time)
- Option B: Skip mtime update (proceed directly to commit)
- Option C: Abort
```

If the user selects:
- **Option A**: Run the touch command
- **Option B**: Skip and proceed to Step 6
- **Option C**: Stop immediately

If `DRY_RUN=true`, print the command without running it.

---

### Step 6 — CHANGELOG Promotion + HITL Approval Gate #3

Check whether `CHANGELOG.md` has a non-empty `## [Unreleased]` block.

- If yes: use **AskUserQuestion**:
  ```
  Promote [Unreleased] to a versioned release?

  - Option A: Promote to vX.Y.Z (enter version, e.g., 1.2.3)
  - Option B: Leave CHANGELOG as-is
  ```

  If promoted: replace `## [Unreleased]` with `## [X.Y.Z] - YYYY-MM-DD`, insert a fresh empty `## [Unreleased]` above it, then touch CHANGELOG.md.

- If no `CHANGELOG.md` or the block is empty: skip silently.

---

### Step 7 — Stage, Commit, Push + HITL Approval Gate #4

If `DRY_RUN=false`, display the exact commands that will run:

```
git add -A
git commit -m "<COMMIT_MESSAGE>"
git push origin <current-branch>
```

Then use **AskUserQuestion**:

```
Ready to stage, commit, and push?

- Option A: Run these commands
- Option B: Abort
```

**Important:** These commands execute without system-level approval prompts (pre-approved via `setup_permissions.py`), but the user has already approved them explicitly at the skill level via this HITL gate.

Halt and report the exact error on any non-zero exit code. Do not retry.

If `DRY_RUN=true`, print each command without running it.

---

### Step 8 — NEW: Offer PR Creation + HITL Approval Gate #5

After successful push, unless `--skip-pr` was set, use **AskUserQuestion**:

```
Create a GitHub pull request?

- Option A: Create PR
- Option B: Skip PR creation
```

If user selects **Option B** or `--skip-pr` was set: jump to Step 10.

If user selects **Option A**: proceed to Step 9.

---

### Step 9 — NEW: Draft and Create PR + HITL Gates #6 & #7

Run the helper script:

```bash
python .agents/skills/hyper-publish/scripts/generate_pr_description.py
```

The script reads:
- `CHANGELOG.md` → extracts `## [Unreleased]` or most recent versioned block
- Formats as markdown:
  ```markdown
  ## Summary
  - <bullet from ### Added>
  - <bullet from ### Changed>

  ## Test plan
  - [ ] <auto-generated test item 1>
  - [ ] <auto-generated test item 2>
  ```
- Outputs to stdout

Parse the output as `PR_BODY`.

**HITL Gate #6 — Description Approval:**

Use **AskUserQuestion** to present the drafted description:

```
Proposed PR description:
  <formatted markdown>

Options:
- Option A: Accept this description
- Option B: Edit it (provide your own markdown in the text field)
- Option C: Use minimal description (just the commit message as title)
```

Bind the result to `PR_BODY`.

**HITL Gate #7 — Final PR Gate:**

Display the exact `gh pr create` command that will run:

```bash
gh pr create --title "<COMMIT_MESSAGE>" --body "<PR_BODY>"
```

Then use **AskUserQuestion**:

```
Create the pull request?

- Option A: Create PR now
- Option B: Abort
```

If approved, run:

```bash
gh pr create --title "$COMMIT_MESSAGE" --body "$PR_BODY"
```

Capture the PR URL from the output and proceed to Step 10.

---

### Step 10 — Report Result

On success:

```
✅ Published successfully.
  Branch:   <branch>
  Commit:   <git rev-parse --short HEAD>
  Message:  "<COMMIT_MESSAGE>"
  Files:    <N files touched>
  Remote:   <git remote get-url origin>
  PR:       <URL>  (or "PR skipped" if no PR was created)
```

On push rejection (non-fast-forward), suggest:

> "Push rejected. Pull latest changes with `git pull --rebase` and try again."

---

## Constraints

- **Never force-push.** Use `git push origin <branch>` only.
- **Never skip HITL gates.** Each step touching files or running git commands requires explicit user approval via AskUserQuestion.
- **Never modify files** other than CHANGELOG.md (and only when user explicitly approves Step 6).
- **Respect `.gitignore`** — only touch files git already tracks or has staged as untracked.
- **No system-level approval prompts** — all git/gh commands are pre-approved in settings.json, but skill-level HITL gates ensure user control.

---

## Examples and Templates

Refer to:
- `.agents/skills/hyper-publish/examples/example_commit_messages.md` — commit message examples from this repo
- `.agents/skills/hyper-publish/examples/example_pr_description.md` — sample PR descriptions
- `.agents/skills/hyper-publish/resources/commit_template.md` — blank commit message template
- `.agents/skills/hyper-publish/resources/pr_description_template.md` — blank PR description template

---

## Troubleshooting

**"No changed files found":**
- Run `git status` to see what's staged vs. unstaged
- Stage changes with `git add <files>` and re-run `/hyper-publish`

**"Push rejected":**
- Pull latest: `git pull --rebase`
- Re-run `/hyper-publish` with your commit message

**"gh pr create failed":**
- Ensure you're authenticated: `gh auth login`
- Ensure your branch is pushed and visible on GitHub

**Permissions still prompting:**
- Run `/hyper-publish` once to auto-patch `.claude/settings.json`
- Or run manually: `python .agents/skills/hyper-publish/scripts/setup_permissions.py`
