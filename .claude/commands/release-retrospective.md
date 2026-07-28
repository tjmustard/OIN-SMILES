---
description: "Refresh the v0.4.4→ release retrospective (markdown + GitHub Pages + Artifact) for a new release"
---

Update the OIN-SMILES release retrospective. Self-contained — do not read any other skill or
handoff to run this.

## Argument

`$ARGUMENTS` is either a version (`v0.4.13`) or `--full`. If empty, ask which release to append.

- **`v0.4.13` — append mode, the common case.** Read *only* the new release's `CHANGELOG` stanza,
  its `docs/agentic-notes/<release>/` folder, and any new `tmCAT-tmPHOTO_xyz_dataset/results-*`
  directory. Add one per-release section. Update the two ELI20 headline tables, both verdict
  blocks, the ladder rail, the state strip, the "What is not known" list and the gap table.
  Leave every earlier per-release section byte-identical unless a source now contradicts it.
- **`--full` — regenerate.** Re-derive every figure from the sources below. Use this when an
  earlier number has been corrected, which in this project happens.

## Where things live — the file-placement contract

**This is settled. Do not re-litigate it, and do not invent a side branch.**

| branch | contents |
|---|---|
| **`main`** | **all source** — the report markdown, the page body, the build script, this command |
| **`gh-pages`** | **only the generated `index.html` + `.nojekyll`.** Nothing else, ever. |

Commit the source **directly to `main`**. Three reasons it cannot live anywhere else:

1. **A slash command only exists on the branch you have checked out.** On a side branch or on
   `gh-pages`, `/release-retrospective` is simply absent when anyone works on `main`.
2. **The report cites lane docs, `levers.py` and `CHANGELOG.md`** — all on `main`. Splitting them
   leaves the next updater without one side.
3. **`gh-pages` is an orphan branch and that is the point.** No shared history with `main` means a
   session rewriting `main` cannot disturb the published site. **Never merge `main` into it**, and
   never put source on it — everything there is publicly served.

A `docs/release-retrospective` branch was tried once and was wrong: it conflated *isolation from
another session's checkout* (which a **worktree** solves) with *needing separate history* (never
true). It was merged to `main` and deleted.

## Step 0 — isolate with a WORKTREE, not a branch (do not skip)

**The primary checkout is frequently on another session's swimlane, and `main` gets rewritten
underneath you** — both happened while this report was first written. Check before touching
anything:

```bash
git -C <repo> rev-parse --abbrev-ref HEAD     # may be a sibling's swimlane, NOT main
git -C <repo> worktree list
```

If the primary checkout is already on a clean `main`, work there. Otherwise add a worktree
**checked out on `main`** — no `-b`, no new branch:

```bash
git worktree add ../oin-retrospective main
```

After **every** commit, re-read `git log --oneline -1` and confirm your commit is still there; if
it vanished, recover it with `git reflog` + `git cherry-pick`.

## Step 1 — read the sources (do not guess a number)

| what | where |
|---|---|
| Per-release narrative and figures | `CHANGELOG.md` |
| Honest baseline + transition matrix | `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/bucket_report_both.md` |
| Cohort bucket reports | `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.{4-sl4,4-regression,5-rebaseline,6-sweep}/` |
| Lane docs for the new release | `docs/agentic-notes/<release>/` |
| **What actually SHIPS** | `src/oinsmiles/oin/levers.py` — `_DEFAULT_ON` and `_HELD_OFF` |
| Gap decomposition and ladder | `docs/agentic-notes/ROADMAP_100_100.md` |
| Suite count | measure it (`discover tests/unit`) or read it from the release's own lane doc — **never copy the previous release's number** |

## Step 2 — write

Two files, kept in sync:

- `docs/agentic-notes/v0.4.12/RETROSPECTIVE_v0.4.4_to_v0.4.12.md` — the full record
- `docs/agentic-notes/v0.4.12/retrospective.page.html` — the page **body only**, no
  `<!doctype>` / `<html>` / `<head>` / `<body>` tags (the build script and the Artifact
  publisher both supply the document shell; the script refuses the file if it finds one)

Per-release section skeleton — match the existing sections exactly:

```
### <version> — <one-line characterisation>
**What shipped.**   levers PROMOTED vs added-but-OFF, stated separately
**Accuracy:**       open with "Real gain" / "Zero, by design" / "Zero, gated" and the number
**Speed:**          default-path only; note separately anything measured but shipped OFF
**Refuted.**        equal weight to wins — most releases here produce more of these
**What it cost.**   regressions, CPU-h, known limitations
Suite: NNN OK.
```

## Step 3 — standing caveats, never drop them

These are why the report is trustworthy. Carry every one forward.

1. **The bucket reports come from four different cohorts** (6719 / 3917 / 936 / 5000). They are
   **not** a time series. Every accuracy figure carries its N and cohort. Never draw a line
   between two cohorts — the figure deliberately does not connect them.
2. **`metrics.elapsed_s` is nested** (a top-level read silently yields `0`) **and is a SUM** over
   up to three SIGKILLed attempts. Never quote `max(elapsed_s)` as a single-run duration.
3. **The corpus speed figure** (`994/5000 = 19.88%` over 30 s, median 7.19 s) is from the v0.4.6
   sweep. Mark it stale until a sweep says otherwise.
4. **Suite count is a rigour proxy, not accuracy.**
5. **A lever measured but shipped default-OFF changes nothing for a user.** Say so explicitly
   every time — that distinction is the spine of this report.

## Step 4 — commit the source to `main`

```bash
git add docs/agentic-notes/v0.4.12/RETROSPECTIVE_v0.4.4_to_v0.4.12.md \
        docs/agentic-notes/v0.4.12/retrospective.page.html
git commit                       # normally — --no-verify is banned
git log --oneline -1             # re-read: confirm your commit survived
```

Do **not** push `main` — see Guardrails.

## Step 5 — rebuild and push the GitHub Pages branch

Only the generated file crosses over. Nothing else belongs on `gh-pages`.

```bash
git worktree add ../oin-ghpages gh-pages
python3 tools/build_retrospective_page.py \
    --source docs/agentic-notes/v0.4.12/retrospective.page.html \
    --out ../oin-ghpages/index.html
git -C ../oin-ghpages commit -am "site: retrospective through <version>"
git -C ../oin-ghpages push origin gh-pages
```

Pushing `gh-pages` **is authorised** — the maintainer asked for GitHub Pages hosting and the site
is already live at <https://tjmustard.github.io/OIN-SMILES/>. That authorisation covers this
branch only; `main` is still do-not-push.

Confirm the rebuild landed rather than assuming it:

```bash
gh api repos/tjmustard/OIN-SMILES/pages/builds/latest --jq '{status, error: .error.message}'
curl -s https://tjmustard.github.io/OIN-SMILES/ | grep -c '<some string you just added>'
```

## Step 6 — republish the Artifact to the SAME url

Read the `<!-- artifact-url: ... -->` comment at the top of the markdown report and pass it as
the `Artifact` tool's `url` parameter. Without it a new URL is minted and the old link goes
stale. If the comment is missing or reads `(not yet published)`, find the URL with
`Artifact action:"list"` before publishing.

- Publish `docs/agentic-notes/v0.4.12/retrospective.page.html`.
- Keep `favicon` (`⚗️📉`) and the `<title>` **stable** — users find the tab by its icon.
- Load the `artifact-design` skill before editing the page, and `dataviz` before touching any
  chart code. If you add a chart colour, run
  `node scripts/validate_palette.js "<hex>" --mode light --surface "#FCFCFB"` and the `--mode
  dark --surface "#151B1D"` counterpart. The current data colour is `#008C9E` light /
  `#00A6BC` dark — both pass the lightness band, chroma floor and contrast checks.
- Write the returned URL back into **both** files if it changed.

## Guardrails

- **Read-only on the corpus.** Do **not** run a sweep, benchmark or A/B unless explicitly asked:
  a full corpus sweep is ~55 CPU-h, the frozen 328-molecule runtime benchmark ~1–2 CPU-h.
- **Pushing is split, and the split is deliberate.** `gh-pages` **may** be pushed — that is how the
  site updates, and the maintainer authorised it. **`main` may not**: report the push command and
  let them run it. A sibling session has been observed pushing `main` against instruction, so do
  not treat a moved remote `main` as permission.
- **Commit normally — `--no-verify` is banned.** `pre-commit` runs the docs-layout guard and
  `ruff`; `commit-msg` rewrites the trailer.
- **Never create a side branch for this work.** Source goes on `main`, generated HTML on
  `gh-pages`, nothing else on either. See the file-placement contract at the top.
- The `docs/` root is closed to four product docs. Everything here goes in a subdirectory.
- If a figure and a source disagree, **re-measure or drop the figure** — do not pick one.
