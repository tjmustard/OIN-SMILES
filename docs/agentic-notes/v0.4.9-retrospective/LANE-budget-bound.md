# Lane 1 retrospective — make the budget a bound

Findings: `docs/agentic-notes/v0.4.9/BUDGET_BOUND_v0.4.9.md`.
Premise refutation: `docs/agentic-notes/v0.4.9/ELAPSED_S_IS_A_SUM_v0.4.9.md`.
This records the route, including the design that was measured and thrown away.

## The lane refuted its own charter twice

**First: the justification.** v0.4.9 was chartered on *"759.9 s against a 300 s budget. That
single number is this release's justification: the budget is not a budget."* It is arithmetic on
a sum — the harness runs up to three separately SIGKILLed attempts and adds their wall-clock into
one field. Split by `tier_passed`, all **4658** single-attempt rows finish within **0.2 s** of
their 300 s cap. The advisory-timeout defect is real (it is in the code, and two direct-call
probes measure it), but the corpus number was never evidence for it.

**Second: the design.** The charter said *"establish where the time is actually spent before
choosing a mechanism — this is the lane's first step, not an implementation detail."* The lane
did that, on two molecules, got a clear and consistent answer (`embed.get_embedding`, 75–77% of
generation on both), threaded the deadline into it, shipped, and measured:

> **ε = +48.4 s on a 30 s budget. Max ratio 2.65× OFF → 2.61× ON. The bound changed almost
> nothing.**

Two molecules agreeing is not a population. `VAFMIA_comp_0` has a third cost profile —
**`chirality._reparse_cip_label_once`, 77.8 s of a 78.6 s generation, 99%** — split roughly
evenly between the `accept_fn` re-encode *inside* the attempt loop and `_select_by_geometry`,
which runs in the adapter **after `generate_3d_structures` has returned** and was therefore
structurally unreachable from the deadline.

**The lesson, stated so the next lane does not repeat it:** *a bound threaded into whichever
function profiled expensive last is not a bound.* Enforcement belongs at the level the budget is
promised at — the whole generation — with checks wherever control returns, not in the one
function two samples happened to indict. This is the charter's own "a sample that only exercises
the common case confirms whatever you already believed", one level up: **the sample was of
molecules, and the belief was about which function to fix.**

## What was NOT needed, against expectation

- **No `fork` + `RLIMIT_CPU`.** The charter anticipated that `EmbedMolecule` might force it.
  `get_embedding` turned out to be a nested *Python* loop, and the genuinely atomic unit is one
  `accept_fn` re-encode — a Python call the deadline can decline to *start*, though not
  interrupt.
- **No CBC `timeLimit`.** The charter's prime suspect is **2.1%** on FOSNEI and **0.35 s** on
  VAFMIA. The topology memo added earlier already collapsed it. Adding a `timeLimit` would have
  measured as no change while looking like a fix.

## What ε actually is, and why it is not zero

ε is **one in-flight `accept_fn` re-encode** — up to ~24 s on `VAFMIA`, dominated by
`_reparse_cip_label_once` at ~2.4 s a call. The bound can decline to start the next one; it
cannot interrupt the one running. Bounding *inside* the CIP/perception layer is a real change to
perception behaviour and is not something to bolt on at the end of a release whose stated job is
"change when work stops, not how long it takes".

**Stated honestly rather than papered over**, from the 32-molecule A/B at a 30 s budget:

| | OFF | ON (design 1) | ON (design 2, shipped) |
|---|---:|---:|---:|
| max ratio | 2.63× | 2.61× | **2.09×** |
| ε | — | +48.4 s | **+32.8 s** |
| molecules over budget | 11 | 11 | **11** |

**The bound compresses the tail; it does not remove it.** The same 11 molecules exceed the
budget in every arm — the bound reduces how far over they go, not how many go over. Byte-identity
held on all 28 that finished in both arms, and no late success was converted to a failure.
**Any claim that v0.4.9 delivers `max(elapsed_s) < 30 s` would be false.**

## Scope held

- **The 22% win was found and not taken.** `get_embedding`'s outer loop calls
  `alternative_ace_mol_list.index(...)` and discards the result — 3711 calls, 198
  eigendecompositions, 15.99 s of 72.02 s on `CAHQEJ`. It belongs to v0.4.10, with the
  measurement attached, because landing an optimization inside the release that claims to
  optimize nothing would muddy its own A/B.
- **The lever is default OFF and was not promoted.** Promotion is an accuracy decision: a 30 s
  bound recovers 37.8 CPU-h per 5k sweep and costs **251 passes = 5.02 points of `byte_exact`**,
  against a headline goal of 100%.
- **15 molecules are out of reach entirely.** They are hard-killed *inside the encoder*
  (`exceeded 300s while encoding`) and never reach the generator. Found by Lane 2. Any claim
  that a generator-side bound delivers Goal B has to exclude them.
