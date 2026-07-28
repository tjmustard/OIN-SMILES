# Documentation Layout (OIN-SMILES)

`docs/` holds two kinds of file with two different audiences. Keeping them apart is the
difference between a `docs/` folder a contributor can read and a 70-file dump of session
artifacts. **The root is closed; a `pre-commit` guard enforces it.**

## The split

| | Product documentation | Agentic coding notes |
| :--- | :--- | :--- |
| **Location** | `docs/` root | `docs/agentic-notes/<release>/` |
| **Audience** | Users and contributors | The next agent session |
| **Answers** | "How does the shipped software behave?" | "What did we measure, try, and refute?" |
| **Lifetime** | Maintained; kept true | Frozen at its commit; never restated |
| **Adding one** | Needs maintainer sign-off | Just write it in the right folder |

**The test:** if it records what a session *measured, tried, or refuted*, it is a note. If
it tells someone *how the shipped software behaves*, it is a product doc.

## Product docs (the complete allowlist)

`README.md`, `OPTIMIZERS.md`, `GENERATION_PIPELINE.md`, `KNOWN_LIMITATIONS.md`.

That is the whole list. It is duplicated in `tools/check_docs_layout.sh`, which the
`pre-commit` hook runs. Adding a fifth means editing that script **in the same commit,
with a reason** — and it means you have convinced the maintainer, not just yourself.

## Writing a note

```
docs/agentic-notes/v0.4.7/MY_MEASUREMENT.md
```

- `<release>` is the release the work is **for**, not the one that was current when you
  started. Targeting v0.4.7? Create `v0.4.7/`. Do **not** park it in `v0.4.6/` because
  that folder happens to exist.
- Raw JSON/CSV backing a report goes beside the report, or in `<release>/data/` if bulky.
  **But data a LATER RELEASE will diff — frozen baselines, gate tallies, transition matrices —
  goes in `measurements/<release>/` instead**, written by `tools/harvest_measurements.py`. The
  split is by audience, same as everything else here: a note's data is read *with the note*; a
  measurement is read *against the next release*. See `measurements/README.md`.
- A research program that genuinely spans releases gets its own topic folder at the
  `agentic-notes/` level — `injectivity/` is the precedent. Use this sparingly; per-release
  is the default.
- Add a row to `docs/agentic-notes/README.md` when you create a new folder.

**Name the commit your numbers were measured at.** A figure without a commit is an order
of magnitude, not a measurement.

## Things that are not notes and not docs

- `spec/process/`, `spec/handoffs/`, `scratchpad/` — gitignored and `pre-commit`-blocked.
  See `.agents/rules/git-workflow.md`.
- `docs/social_media/` — gitignored; owned by the `social-post` skill.
- Sweep output — never under `docs/`, never under `/tmp`. It goes in the dataset
  directory as `results-*/`.

## When a note graduates

Do **not** promote a note by moving it. When a finding becomes something a user needs,
write it into the product doc in the user's language and leave the note where it is as
the evidence trail. Cross-link the note from the product doc if the derivation matters.

## If the guard blocks your commit

```
❌ COMMIT BLOCKED — new file at the docs/ root: docs/FOO.md
```

You wrote a session note to the root. Move it:

```bash
git restore --staged docs/FOO.md
mkdir -p docs/agentic-notes/v0.4.7
git mv docs/FOO.md docs/agentic-notes/v0.4.7/   # or plain mv if never committed
```

Do **not** reach for `--no-verify` — it is banned for routine commits and skips the
`commit-msg` trailer rewrite as well. See `.agents/rules/git-workflow.md`.
