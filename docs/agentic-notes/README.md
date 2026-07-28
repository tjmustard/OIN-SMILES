# Agentic coding notes

The working record of what agentic coding sessions measured, tried, and refuted while
building OIN-SMILES. These are **not** product documentation. They are dated, they
contradict each other across releases, and a large share of them exist specifically to
record that an approach **did not work** — which is what stops the next session spending
a week rediscovering it.

For how the shipped software behaves, see [`../README.md`](../README.md) and the four
product docs beside it.

---

## Where does my file go?

**The test:** if it records what a session *measured, tried, or refuted*, it is a note.
If it tells a user or contributor *how the shipped software behaves*, it is a product doc.

| You are writing | It goes in |
| :--- | :--- |
| A measurement report, sweep result, or A/B | `agentic-notes/<release>/` |
| A lane / swimlane write-up | `agentic-notes/<release>/` |
| A refuted hypothesis or negative result | `agentic-notes/<release>/` |
| A status snapshot or session handoff record | `agentic-notes/<release>/` |
| Raw JSON backing a report | beside that report, or `<release>/data/` |
| A number a LATER RELEASE will diff (baseline, gate tally, transition matrix) | `measurements/<release>/` — via `tools/harvest_measurements.py` |
| A multi-release research program | its own topic folder (see `injectivity/`) |
| A behaviour a user would hit | `docs/KNOWN_LIMITATIONS.md` |
| How a shipped subsystem works | `docs/GENERATION_PIPELINE.md` or a new product doc |

`<release>` is the release the work is **for**, not the release that was current when you
started. Work targeting v0.4.7 goes in `v0.4.7/` — create the folder; do not park it in
`v0.4.6/` because that folder already exists.

**The `docs/` root is closed.** A `pre-commit` guard (`tools/check_docs_layout.sh`)
rejects new files there. If something genuinely belongs at the root as product
documentation, get maintainer sign-off and add it to the allowlist in that script — in
the same commit, with a reason.

### When a note graduates

Notes do not get promoted by moving them. When a finding becomes something a user needs,
**write it into the product doc in the user's language** and leave the note where it is as
the evidence trail. Cross-link the note from the product doc if the derivation matters.

---

## Index

### Per release

| Folder | Contents |
| :--- | :--- |
| [`v0.3.7/`](v0.3.7) | Round-trip residual tail — R4 triage & fix. |
| [`v0.4.1/`](v0.4.1) | Full-corpus failure-mode distribution (a screening floor, not a headline). |
| [`v0.4.2/`](v0.4.2) | Per-molecule before→after for the v0.4.2 accuracy wave. |
| [`v0.4.3/`](v0.4.3) | Structure-quality wave: the elimination study, the distortion research, the vdW acceptance gate. Raw telemetry in [`v0.4.3/data/`](v0.4.3/data). |
| [`v0.4.4/`](v0.4.4) | Accuracy release: the gap decomposition, OIN-direct assembly validation, encoder robustness (SL5), reliability (SL4). |
| [`v0.4.5/`](v0.4.5) | The canonical-OIN release — 18 lane reports: canonical body/slots, encode-fail histogram, boron cage, atom-count hydrogen, valence order/search, renumbering instability, perf, the promotion gate. |
| [`v0.4.6/`](v0.4.6) | The 5,000-molecule sweep, the accuracy **instrument** (both error directions measured), the eta accept-scored A/B, and the H-faithful negative result. JSON artifacts alongside. |
| [`v0.4.7/`](v0.4.7) | Five lanes, four of them negative: the `OIN_ACCEPT_SCORED` **do-not-promote** verdict and the attachment check that is its missing safety condition, the frozen slow cohort + two-arm byte-identity gate, the encode floor's three cost regimes, and the boron generation ceiling (whose own fast-fail proposal this release refutes). |

### Cross-cutting

| Folder / file | Contents |
| :--- | :--- |
| [`injectivity/`](injectivity) | The Y1/Y2/Y3 encoder blind-spot program. Spans v0.4.5–v0.4.6; kept together because the argument runs across releases. |
| [`v0.4.5-retrospective/`](v0.4.5-retrospective/README.md) | Master retrospective for v0.4.5 + v0.4.6 — five WAVE documents and eighteen LANE documents. Start at its own README. |
| `ROADMAP_100_100.md` | The release ladder toward `byte_exact = 100%` and `max(elapsed_s) < 30 s`. Git-durable twin of the gitignored `spec/handoffs/roadmap-100-100/`. |

---

## Reading these honestly

Three cautions that apply to the whole tree:

1. **Numbers are measured at a commit.** Most reports name the commit; where one does not,
   treat the number as an order of magnitude and re-measure.
2. **A later note may refute an earlier one** without the earlier one being edited. The
   retrospective and `ROADMAP_100_100.md` carry the current reconciliation.
3. **`v0.4.6/METRIC_FALSE_POSITIVES.md` is load-bearing for every accuracy figure written
   before it.** The reported accuracy was over-stated by roughly 5.7 points; reports
   predating that measurement have not been restated.
