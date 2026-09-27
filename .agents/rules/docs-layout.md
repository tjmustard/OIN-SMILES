# Documentation Layout (OIN-SMILES)

`docs/` holds two kinds of file with two different audiences. Keeping them apart is the
difference between a `docs/` folder a contributor can read and a 70-file dump of session
artifacts. **The root is closed; a `pre-commit` guard enforces it.**

## The split

| | Product documentation | Agentic coding notes |
| :--- | :--- | :--- |
| **Location** | `docs/` root | `spec/process/` + `spec/handoffs/<release>/` (gitignored) |
| **Audience** | Users and contributors | The next agent session |
| **Answers** | "How does the shipped software behave?" | "What did we measure, try, and refute?" |
| **Lifetime** | Maintained; kept true | Local to the machine; never restated |
| **Adding one** | Needs maintainer sign-off | Just write it in the right folder |

**The test:** if it records what a session *measured, tried, or refuted*, it is a note. If
it tells someone *how the shipped software behaves*, it is a product doc.

## Product docs (the complete allowlist)

`README.md`, `OPTIMIZERS.md`, `GENERATION_PIPELINE.md`, `KNOWN_LIMITATIONS.md`.

That is the whole list. It is duplicated in `tools/check_docs_layout.sh`, which the
`pre-commit` hook runs. Adding a fifth means editing that script **in the same commit,
with a reason** — and it means you have convinced the maintainer, not just yourself.

## Writing a note (owner decision, 2026-09-26)

```
spec/process/process_2026MMDD_<slug>.md     # the session narrative (/hyper-process-document)
spec/handoffs/v0.4.NN/NEXT.md               # next steps for the following session
```

- Both directories are **gitignored** and `pre-commit`-blocked: a note stays on the machine that
  wrote it. That is deliberate. Anything a later release or another machine needs goes in the
  tracked record instead:
  - **data a LATER RELEASE will diff** — frozen baselines, gate tallies, transition matrices —
    in `measurements/<release>/`, written by `tools/harvest_measurements.py` (see
    `measurements/README.md`);
  - **the decision and its evidence** — in `CHANGELOG.md` and, for a lever, its evidence block
    in `src/oinsmiles/oin/levers.py`;
  - **shipped behaviour** — in a product doc.
- `<release>` is the release the work is **for**, not the one that was current when you started.
- `docs/agentic-notes/` holds the notes written before 2026-09-26 and stays as the historical
  evidence trail (with its `README.md` index). **Add nothing new there.**

**Name the commit your numbers were measured at.** A figure without a commit is an order
of magnitude, not a measurement.

## Things that are not notes and not docs

- `scratchpad/` — gitignored and `pre-commit`-blocked. (`spec/process/` and `spec/handoffs/`
  are too, but they ARE the notes home now — see above.) See `.agents/rules/git-workflow.md`.
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

You wrote a session note to the root. Move it to the notes home:

```bash
git restore --staged docs/FOO.md
mv docs/FOO.md spec/process/process_2026MMDD_foo.md   # gitignored: it leaves the commit
```

Do **not** reach for `--no-verify` — it is banned for routine commits and skips the
`commit-msg` trailer rewrite as well. See `.agents/rules/git-workflow.md`.
