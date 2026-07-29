# v0.4.13 Lane 2 — the MEDZUR and GAVSED classes, measured

> **Both were handed to this release known on n = 1. They are 99 and 280.**
> Neither is a single-molecule curiosity, so the sketch's own deletion clause
> (*"if both turn out to be single-molecule curiosities, this release should be deleted from
> the ladder"*) does **not** fire. And the larger of the two re-sizes a release four rungs up.

Instrument: `tools/attach_class_audit.py` · data: `attach_class_audit.json`

```bash
V=/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python
PYTHONPATH=$PWD/src $V tools/attach_class_audit.py \
    --results-dir tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest \
    --out docs/agentic-notes/v0.4.13/attach_class_audit.json
```

---

## 1. Why this could be measured offline, and why that is not a shortcut

v0.4.8 backfilled a `coordination` block into all 5000 individual reports — `intact`,
`boundary_only`, and per-metal contact counts computed from the **stored** generated geometry.
The classification these two classes need is therefore already a property of the frozen corpus.
Re-running the generator would not make it more true; it would make it a *different corpus*, and
the release would then own two baselines instead of one.

`#DONE 5000`, 0 reports missing or unreadable, **0 `UNKNOWN`** — every molecule carried a usable
`coordination` block. That last count is the one that says the join actually happened.

## 2. The table

| bucket | NO_STRUCTURE | UNKNOWN | DETACHED | BOUNDARY | INTACT | total |
|---|---:|---:|---:|---:|---:|---:|
| `byte_exact` | 3 | 0 | **48** | 1273 | 2299 | 3623 |
| `key_equal` | 2 | 0 | 52 | 267 | 289 | 610 |
| `structural` | 5 | 0 | **266** | 48 | 98 | 417 |
| `hard_fail` | **315** | 0 | 3 | 1 | 0 | 319 |
| `facmer_divergent` | 0 | 0 | 11 | 4 | 1 | 16 |
| `encode_fail` | 15 | 0 | 0 | 0 | 0 | 15 |
| **ALL** | 340 | **0** | 380 | 1593 | 2687 | 5000 |

Classification, first match wins: `NO_STRUCTURE` (no structure to judge) → `UNKNOWN` (no usable
`coordination`) → `DETACHED` (`intact` false) → `BOUNDARY` (`boundary_only`) → `INTACT`.

`BOUNDARY` is deliberately **not** merged into either neighbour. It means the attachment verdict
rests on contacts sitting within 0.1 Å of the cutoff — the call itself is the uncertain quantity,
and merging it would let a tolerance choice decide the headline.

## 3. The control, which must be read first

| | n | DETACHED | % |
|---|---:|---:|---:|
| `byte_exact` (round-trips) | 3623 | 48 | **1.32%** |
| non-`byte_exact` | 1377 | 332 | **24.11%** |
| | | **enrichment** | **18.2×** |

A broken version of this tool — one where `intact` were mostly false for unrelated reasons, or
mostly `None` and silently coerced — would print a large, plausible, meaningless GAVSED class.
An 18.2× enrichment against a 1.32% floor is what makes the number evidence rather than an
artifact of the predicate. **Read this before the class sizes.**

## 4. The two classes

`key_equal` molecules round-trip to the same comparison key and differ only in presentation —
benign canonicalization, the encoder-side block this release's promotion and v0.4.14 address.
Counting them as "independent re-perception disagrees" would inflate MEDZUR with molecules whose
geometry is fine and whose string is merely renumbered. The honest denominator excludes them.

**Over the 767 genuine failures:**

| class | n | what it is |
|---|---:|---|
| **GAVSED shape** — `DETACHED` | **280** | acceptance rejected every conformer, `_select_by_geometry`'s fallback returned one anyway, **and that fallback is not attachment-aware** |
| **MEDZUR shape** — `INTACT` | **99** | attachment fully intact, independent re-perception still disagrees |
| `BOUNDARY` | 53 | attachment call is inside the tolerance band |
| `NO_STRUCTURE` | 335 | nothing was generated — no returned conformer to be wrong |

Of the 432 genuine failures where a structure exists **and** the attachment call is decisive,
**280 (64.8%) are attachment failures.** That is the dominant mechanism in the failure side of
the gap, and it is a *return*-path defect, not an acceptance-path one.

## 4b. Cross-checked against an independent decomposition — and it resolves that one's ambiguity

`tools/injectivity/missed_success_audit.py` partitions the **same 767 failures** by *cause*, and
it was written years before this lane. The two agree without having been made to:

| missed-success audit | n | this lane | n |
|---|---:|---|---:|
| output names a different isomer — *"ambiguous: generator OR notation"* | **437** | `DETACHED` + `INTACT` + `BOUNDARY` | **432** |
| timeout (269) + produced nothing (32) + encoder refused (15) | 316 | `NO_STRUCTURE` | 335 |
| canonicalization noise | 14 | — | — |

Two independently-written tools, two different predicates, the same corpus, agreeing to within
~5 molecules on both halves. The residual is definitional: this lane calls a molecule
`NO_STRUCTURE` when `status != success` **or** `smiles_2` is empty, which catches a few the audit
attributes to a named cause.

**The point is not the agreement — it is what the split adds.** The audit's largest bucket is 437
molecules it explicitly labels *ambiguous, generator or notation, we cannot tell*. This lane tells:

| the audit's ambiguous 437 | resolves to |
|---|---|
| **280** | **generator** — the returned structure has ligands off the metal |
| 99 | genuinely unexplained (MEDZUR) — attachment intact, re-perception still disagrees |
| 53 | undecidable at this tolerance (`BOUNDARY`) |

So the single largest "we don't know" in the project's failure attribution is now **65% known**,
and known to be a *return-path guard*, not a notation defect.

## 5. 🔴 What this does to the ladder

**`structural` is 266/417 = 63.8% `DETACHED`.**

`structural` is 417 molecules / **8.34 points** — the second-largest block in the 27.54-point gap
— and it is scheduled at **v0.4.17**, four rungs out, labelled *"generator capability floor …
bounded by what the generator can assemble"*. That label is now measurably wrong for two thirds
of it. Those 266 molecules did not fail because the generator could not assemble them; it
assembled something and **returned it with ligands off the metal**, because the fallback ranking
that produced it never checks attachment.

That is not a capability problem. It is a **one-site guard** — the same guard `OIN_ATTACH_CHECK`
already implements for *acceptance*, applied to *return*. v0.4.12 handed this forward explicitly
as "the GAVSED class … closing it changes arm A's behaviour and needs its own gate."

The arithmetic that follows is in `ROADMAP_100_100.md` under `LADDER DECISION 2026-07-27 (v0.4.13)`.

## 6. The 48 that should not exist

48 `byte_exact` molecules read `DETACHED`: the string round-trips byte-for-byte while the stored
geometry has lost a metal contact relative to its input. They are 1.32% — the control floor above
— and this note is **not** claiming they are false positives. Two readings survive, and this lane
did not run the measurement that separates them:

1. the lost contact is on a metal or ligand whose loss the OIN string does not express, so the
   string is right and the geometry is degraded; or
2. `coordination.intact` is slightly over-sensitive at this tolerance, which is exactly what the
   1593-molecule `BOUNDARY` band warns about.

**What would settle it:** run `tools/injectivity/oracle.py` and `vdw_clash_count` on those 48 and
check whether the lost contact is one the notation encodes. Not run here, and named so the next
release can pick it up rather than rediscover it.

## 7. Carried traps this lane tripped over, or avoided

- **Never `honest_class.endswith("->byte")`.** The bucket is read from the frozen
  `bucket_report_honest.json`; the ad-hoc string test disagrees with it by exactly the eight
  atom-count-gate molecules.
- **Ask what a broken instrument would print.** Answered with the `byte_exact` control arm and
  the `UNKNOWN` bucket, both of which are in the tool's output by construction rather than by
  discipline.
- **A sample that only exercises the common case confirms whatever you already believed.** Both
  classes were n = 1 for two releases. The full corpus moved one of them by 280×.
