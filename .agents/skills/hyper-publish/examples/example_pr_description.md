# Example PR Description

This shows the format that `/hyper-publish` generates when creating a pull request after pushing a commit.

---

## Example: From CHANGELOG [0.4.1] Block

**Title (from commit message):**
```
Release v0.4.1: Framework upgrades with /hyper-update and /hyper-recover
```

**Body:**

```markdown
## Summary
- **`hyper-publish` skill**: New `/hyper-publish` command that touches all git-changed files to update mtime to the current system time, then commits and pushes to the remote. Designed for post-session IP timestamp protection.
- **`hyper-update` skill**: Smart framework upgrade command that fetches latest upstream framework, auto-updates non-sensitive files, and guides user through section-by-section collaborative merge of sensitive customizable files (CLAUDE.md, AGENTS.md, IDE rules).
- **`hyper-recover` command**: Backup management and interrupt recovery. Creates timestamped backup directories with lazy initialization, detects interrupted merges, and provides restore operations.

## Test plan
- [ ] Run `/hyper-publish` with a commit message — verify mtime is updated and commit is pushed
- [ ] Run `/hyper-update` — verify upstream changes are fetched and merged correctly
- [ ] Run `/hyper-recover --file=CLAUDE.md --date=<date>` — verify backup restoration works
```

---

## Structure

The generated PR description always has:

1. **## Summary** — bulleted list of changes from CHANGELOG `### Added` / `### Changed` sections
2. **## Test plan** — auto-generated checklist (one [ ] item per added feature, with descriptive action)

---

## How It's Generated

The `generate_pr_description.py` script:

1. Reads `CHANGELOG.md`
2. Extracts the `## [Unreleased]` block (or most recent versioned block if Unreleased is empty)
3. Parses `### Added`, `### Changed`, `### Removed` subsections
4. Formats as markdown with `## Summary` + `## Test plan` sections
5. Outputs to stdout as a single string

If the CHANGELOG is missing or has no content, falls back to a minimal template:
```markdown
## Summary
Update to reflect recent changes.

## Test plan
- [ ] Manual verification of changes
```

---

## Template Reference

See `resources/pr_description_template.md` for the blank template you can customize.
