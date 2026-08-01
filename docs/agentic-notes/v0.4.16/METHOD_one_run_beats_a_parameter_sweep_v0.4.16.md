# The knee curve is arithmetic, not a parameter sweep — and the instrument was wrong twice

**v0.4.16's two transferable results are about method, not about the lever.** Both are the same
shape: an instrument that would have printed a plausible, precise, self-consistent number, caught
only because something else was checked against it.

---

## 1. One instrumented run replaces a run per candidate bound

The charter budgeted **~1–2 h per point** to plot recovered-molecules against the search bound.
That is unnecessary, and the reason is structural rather than clever.

In `generator3d/__init__.py`:

- `incumbent_hit` is set to the **first** conformer whose `accept_fn` verdict is
  `ACCEPT_INCUMBENT`, and it is returned as the **sole pool member** — *regardless of how much
  longer the pool fills*.
- `early_hit` breaks both fill loops the instant a string-exact conformer appears.

⇒ Bounding the search at *N* changes the answer for **exactly one class of molecule**: those whose
string-exact hit lies beyond *N*. Every molecule that never hits is answer-neutral by construction,
and every molecule whose hit is at or below *N* is recovered identically.

So recording, per molecule, **the smallest bound that still recovers it** (`min_bound`) plus an
elapsed stamp at each post-incumbent evaluation makes **both** curves — recovered(*N*) and
runtime(*N*) — arithmetic over a **single** unbounded run.

**Three design decisions that make the derived numbers trustworthy rather than merely cheap:**

1. **Record the ANSWER, not a counter.** The telemetry emits `min_bound` — the bound that
   recovers this molecule — not the raw `_since_incumbent`. A future reader of a frozen JSON
   cannot see the loop, so any off-by-one convention held only in a comment is a bug waiting for a
   second reader. (The first draft of that loop *did* advance the counter on the evaluation that
   recorded the incumbent, shifting every point of the curve by one.)
2. **`has_incumbent` is the honesty field.** A molecule that never records an incumbent falls
   through to the energy-sorted pool, where truncation is **not** answer-neutral. Those are
   **excluded and counted**, and the summary states the exclusion on its own line — "those are not
   zeros; they are unmeasurable by this method."
3. **Bound 0 is a wiring gate that a broken bound cannot pass.** It returns the incumbent the
   instant it is recorded, which is byte-identically what the lever-OFF arm returns. A broken bound
   and a working one print the same recovered count at large *N*; **they differ at 0.** It is
   pinned as a unit test *and* run as the first, blocking arm of the live confirmation.

⚠ **A derived curve is a PREDICTION.** Reading one as an end-to-end result is exactly the hole
v0.4.13 and v0.4.14 fell into with the offline re-score. The derivation says *which* bound to test;
a live arm says whether it holds.

### The run conditions were verified comparable, not assumed

`OIN3DGenerator(timeout=)` is advisory and the embed loop checks its deadline *between* attempts,
so CPU starvation shrinks the pool — which biases **accuracy**, not merely timing. Against the
frozen v0.4.15 arm over the same molecules:

```
v0.4.15 ON   total 2854 s   median 18.88 s
v0.4.16 run  total 2828 s   median 17.90 s
ratio        0.99x total, 0.95x median
```

---

## 2. 🔴 A normalizer that was wrong 57 times in 109 — and printed a clean table

Lane 2 classifies why a generated structure's re-perception disagrees, and the load-bearing
question is *did the heavy-atom graph survive?* — because that separates **perception** (bond
orders, aromaticity, H, charge) from **construction** (a bond actually broke).

The first implementation normalized **strings**: strip bracket atoms to their element, drop bond
symbols, uppercase. It produced a complete, precise, plausible table.

Checked against an RDKit **canonical** heavy-graph comparison, it disagrees on **57 of 109**
molecules. A coin flip. SMILES ring-closure digits and atom ordering are arbitrary labels, so an
identical graph written two ways reads as different.

| | broken | corrected |
|---|---:|---:|
| `SKELETON` (→ construction) | **74 / 172 (43%)** | **4 / 172 (2.3%)** |

**A factor of 18, and both look like a finished measurement.** The broken version would have
confirmed the roadmap's existing assumption — that `structural` is construction work — which is
precisely why it would not have been questioned.

Two more of the same shape in the same instrument, each found the same way:

- **23 bodies carry a `RAW:` sentinel** the key builder prefixes when canonicalization failed. It
  is not part of the SMILES. Parsed verbatim they read as *"RDKit cannot parse this molecule"* —
  13.4% of the population reported as UNCLASSIFIED, when the truth was a missing four-character
  strip.
- **Chiral tags survive `MolToSmiles`**, so an **inverted stereocentre read as a broken graph**.
  `DIVZOY_comp_0` differs only as `[P@@]` vs `[P@]`. Comparing with *and* without stereo separates
  them — and `STEREO_INVERSION` turns out to be a distinct mechanism in the same family as the
  enantiomer class v0.4.17 owns.

### The generalisable rule, restated

This project already says *"ask what a BROKEN version of your instrument would print; if it is the
same thing, you have measured nothing."* v0.4.16 adds the harder case:

> **When a normalizer decides your headline, validate it against an independent canonical
> comparison — not against your reading of a handful of examples.**

Reading examples is what caught it *after* the fact: `BEDLII_comp_0` had been read during planning
as aromatic-perception drift, and the tool called it construction. **The disagreement between a
read example and the instrument is the signal.** Neither alone was enough — the eyeball read was
too shallow on a 300-character macrocycle, and the instrument was confidently wrong.
