# Commit Message Template

Use this as a reference when editing a commit message manually in `/hyper-publish`.

---

## Format

```
<imperative verb> <what>: <brief context or reason>

<Optional body — 1-3 sentences explaining WHY, not what>
<Reference any related issues, PRs, or decisions>

<Optional footer: Co-Authored-By: Name <email>>
```

---

## Guidelines

### Subject Line (required)
- **Imperative mood**: "Add", "Fix", "Update", "Release" — as if giving a command
- **No period** at the end
- **≤72 characters** (aim for ≤50 for best GitHub display)
- **Specific**: Say what changed, not "fix bug" but "fix authentication timeout in login flow"

### Body (optional but recommended for significant changes)
- Explain **WHY** the change was made, not what was changed (the diff shows that)
- Keep to 1-3 sentences
- Wrap at 72 chars
- Blank line between subject and body

### Footer (optional)
- Co-authored credit: `Co-Authored-By: Claude Haiku 4.5 <noreply@anthropic.com>`
- Issue reference: `Fixes #123` or `Related: #456`

---

## Examples

### Release
```
Release v0.4.1: Framework upgrades with /hyper-update and /hyper-recover
```

### Feature
```
Add hyper-publish skill for timestamp-safe git publishing

Implements post-session IP protection by updating file mtimes
before commit, ensuring the repository history reflects when
you published, not when the model wrote.
```

### Simple fix
```
Fix typo in README.md
```

### Update with context
```
Update CHANGELOG to reflect recent bug fixes

These fixes were discovered during integration testing and
must be documented before the next release.
```

---

## Anti-patterns

❌ "fix stuff"  
❌ "WIP"  
❌ "asdf" (placeholder)  
❌ All caps or run-on sentences  
❌ Commits without verbs ("The new feature")  

✅ "Fix race condition in token counter"  
✅ "Add streaming support to API responses"  
✅ "Release v0.5.0: Performance optimizations"
