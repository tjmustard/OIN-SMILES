# PR Description Template

Use this template when `/hyper-publish` asks you to manually edit the PR description.

---

## Format

```markdown
## Summary
- <What changed: one change per bullet>
- <Why it matters: brief context>

## Test plan
- [ ] <How to verify this works: concrete action>
- [ ] <Edge case: what if X happens>
- [ ] <Integration: does this break anything else>
```

---

## Guidelines

### Summary Section
- **Bullet list** of what changed (1-3 bullets, unless significant feature)
- Be **specific**: "Fixed race condition in token counter" not "Fixed a bug"
- Include **why** it matters: performance improvement, security fix, UX improvement

### Test Plan Section
- **Checklist format**: `- [ ] <action>`
- **Concrete actions** (not vague): "Run `npm test` and verify all 42 tests pass" not "Test it"
- **Cover main path + edge cases**: happy path, error handling, integration points
- Reviewers will use this to validate your code

---

## Real Example

```markdown
## Summary
- Added `/hyper-publish` skill for timestamp-safe git publishing
- Updated mtime on all changed files before commit
- Implements IP protection by reflecting post-session time in repo history

## Test plan
- [ ] Create a test branch and make a file change
- [ ] Run `/hyper-publish "Test commit"` and verify mtime is updated to current time
- [ ] Verify commit message matches what was provided
- [ ] Confirm push completed successfully to origin
- [ ] Create test with `--dry-run` flag and verify no files were touched
```

---

## Anti-patterns

❌ No summary, just "Fix stuff"  
❌ Test plan with vague items: "Test it" or "Make sure it works"  
❌ Test plan that requires special knowledge (describe steps, don't assume)  
❌ Empty sections or copy-paste boilerplate  

✅ Clear bullet points with "What" and "Why"  
✅ Concrete, reproducible test steps  
✅ Covers both happy path and edge cases  
✅ Concise but specific enough for a reviewer to validate
