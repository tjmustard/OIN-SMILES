# RedTeam Report — Stereo Roadmap Phase 1: "Preserve the Signal" (Winding Plumbing)

- **Target**: `spec/active/Draft_PRD.md` (v0.1.0 Draft)
- **Reviewer**: Red Team Agent
- **Date**: 2026-07-03
- **Verdict**: Draft is unusually disciplined and mostly correct. The §5.1 consumer
  map was verified complete against a repo-wide grep. **One material correctness
  gap** (winding preservation to `ParsedOIN` is conditionally silent — CQ/WI below)
  and **one semantic ambiguity** (`>`/`<`/`^` are the *heading* marker, not a pure
  "winding" concept) should be resolved before compilation. Two of the three
  stress items in §9 are cleared by evidence and downgraded.

**Evidence gathered (grounding for every finding below):**
- `parse_inline_string` callers: `generation/oin_parser.py:477`, `oin/parser.py:25`,
  `tests/unit/test_inline.py:12,26`. No others. §5.1 map is **complete**.
- `SLOT_REGEX` sites: definition `oin/inline.py:44`; `finditer` `:350`; `sub("",…)` `:356`.
  A second, independent slot-suffix stripper exists at `oin/inline.py:89`
  (generate-side V2.4→inline), **not** covered by §5.1 — see §5 analysis.
- ProcessPoolExecutor boundary (`molassembler_adapter.py:2220-2233`) pickles a **plain
  dict of primitives** — `ParsedOIN`/`OINVector`/`SlotAssignment` never cross it.
  Template path (`:2198`) runs **in-process**. Execution-Note (b) is a non-issue.
- Ferrocene fixture (`tests/candidate_outputs/ferrocene_oin.txt`):
  `[Fe_LIN].[cH]{0>}1[cH]{0}…` → `geo_code="LIN"`, which **is** in `TEMPLATES`
  (`oin_parser.py:410`). US-002's example path is live — but only because LIN has a template.

---

## §1 Introduction & Goals — Analysis

### Clarifying Questions
- **Root-cause completeness.** §1.1 says the signal is "destroyed at the first parse
  step" because `SLOT_REGEX` doesn't capture the suffix. Confirmed for the *parse-back*
  path. But there is a **second** suffix-stripping site — `oin/inline.py:89`:
  `int(slot_str.replace('^','').replace('>','').replace('<',''))` in the V2.4→inline
  *generate* path. Is the Draft's "first parse step" claim scoped to OIN→XYZ only, and
  is the generate-side stripper intentionally left alone (it is out of the vector_data
  contract, but it *is* a winding-erasure site)?
- **Diagnostic invariance.** US-004.2 requires `test_haptic_face_winding` to stay
  `expectedFailure`. That test compares **re-encoded 3D output**, not `ParsedOIN`.
  Confirm the team accepts that Phase 1 can make `ParsedOIN.vectors[].winding` correct
  while the diagnostic still (correctly) fails — i.e. the two are decoupled and no one
  will "notice winding is now captured" and prematurely flip the test.

### What-If Scenarios
- **`^` legacy heading marker.** The generate side (`oin/inline.py:245`) *normalizes*
  `^` → `>`, proving `^` is a legal historical suffix. The PRD regex `([><])?` (and
  today's `[><]?`) does **not** match `^`. A hand-authored or legacy `{0^}` marker would
  therefore be neither captured **nor stripped**, leaking a literal `^` into the output
  SMILES and yielding `winding is None`. Pre-existing, but the Draft's regex change is
  the moment to decide whether `^` is in-alphabet.

### Points for Improvement
- Add one sentence to §1.1 acknowledging the `oin/inline.py:89` generate-side stripper
  as a *known, out-of-scope* second erasure site, so a future reader doesn't "discover"
  it and think the map is wrong.
- State the alphabet of the winding suffix explicitly: `{>, <}` captured, `^` **rejected
  by regex** (unchanged). If `^` should map to `>`, widen to `([><^])?` and normalize.

---

## §2 Confidence Mandate — Analysis

### Clarifying Questions
- 9/10 is defensible for the enumerated blast radius. But the confidence rests on
  "winding is provably preserved" — which, as §5 below shows, is **conditionally**
  preserved at the `ParsedOIN` layer. Should the confidence be re-anchored to US-001
  (unconditionally solid) vs US-002 (gated on `tmpl_vectors is not None`)?

### What-If Scenarios
- The residual "1 point" is assigned entirely to RISK-2 (the `oin/parser.py:34` broken
  unpack). The genuinely riskier item — silent winding drop for template-less geometries
  — is unaccounted for in the score.

### Points for Improvement
- Split the confidence: "US-001 capture: 10/10; US-002 threading to `ParsedOIN`: 8/10,
  gated by the template-presence branch (see NFR below)."

---

## §3 Scope — Analysis

### Clarifying Questions
- §3.1.4 threads winding "in the inline branch (`parse`, ~line 485)". Verified: the copy
  site is the `for … in vector_data` loop at `oin_parser.py:485`, **inside**
  `if tmpl_vectors is not None:` (`:484`). Is it acceptable that winding is copied onto
  `OINVector` **only when the geometry code has a `TEMPLATES` entry**? For any
  `geo_code` not in `TEMPLATES` (notably `"NON"`, and any eta geometry lacking a
  template), `vectors = []` and **winding is silently dropped at the `ParsedOIN`
  boundary** — the exact layer §1.2 promises "provably preserved."

### What-If Scenarios
- **Template-less haptic ligand (the motivating class).** Phase 1 exists to unblock
  Phase 3 *haptic-face control*. Haptic/eta ligands are the ones most likely to route
  through `geo_code="NON"` or an eta path with no `TEMPLATES` match. In that case Phase 1
  delivers `winding` on the `SlotAssignment` (US-001 ✓) but **not** on any `OINVector`
  (US-002 ✗) — the guarantee evaporates precisely for the ligand family Phase 3 targets.
  Ferrocene survives only because `LIN` happens to be a template key.
- **Slot-index overflow drop.** Even inside the template branch, `oin_parser.py:485`
  guards `if slot_idx < len(tmpl_vectors)`. A `slot_idx` ≥ template length (geometry/OIN
  mismatch, or a higher-coordination fragment) silently skips the `OINVector` **and its
  winding**. A second conditional erasure the Draft doesn't call out.

### Points for Improvement
- Promote to an explicit **in-scope acceptance criterion**: "winding survives to
  `ParsedOIN.vectors[].winding` for **all** geometries that produce vectors, including
  the `tmpl_vectors is None` / `NON` path." Either (a) build `OINVector`s (with winding,
  minimal `vector`/placeholder) even when `tmpl_vectors is None`, or (b) explicitly
  document winding-at-`ParsedOIN` as **template-geometry-only** in §3.2 Out-of-Scope so
  Phase 3 knows not to rely on it for `NON`/eta.
- Add a negative test: a `NON`-geometry (or template-less eta) OIN with `{0>}` — assert
  the *documented* behavior (winding present, or explicitly `[]`), so the gate is pinned.

---

## §4 User Stories — Analysis

### Clarifying Questions
- US-002 AC-1 uses ferrocene and will pass (LIN is templated). Does the team want a
  **second** US-002 fixture on a template-less geometry to prevent the §3 gap from
  hiding behind a lucky fixture choice?

### What-If Scenarios
- **US-003 positional back-compat vs the actual asserts.** `test_inline.py:21,30` assert
  `assertEqual(vectors[0], (1, 0, 0))`. `SlotAssignment(1,0,0,None) == (1,0,0)` is
  **False** (length 4 ≠ 3), so these asserts *must* change (Draft correctly notes this in
  RISK-4). But US-003 AC-2 claims "`sa[0..2]` still returns lig_rank/atom_idx/slot
  (positional back-compat)" as a *safety net for missed consumers*. That net only helps
  **slice/index** readers; every **unpack** consumer (`a, b, c = sa`) breaks loudly on
  arity — which is actually desirable, but the Draft frames arity-safety as graceful
  degradation. It is graceful for `sa[i]`, fatal for `a,b,c = sa`. Confirm the framing.

### Points for Improvement
- Reword US-003/RISK-1: "positional **index** reads stay correct; positional **unpack**
  reads fail fast (intended)." Precision here prevents a reviewer from assuming an
  unpack-site is safe.

---

## §5 Technical Specifications — Analysis

### Clarifying Questions
- **C2 (`oin/parser.py:34`) truly dead?** Verified: `parse()` inline branch reaches
  `for rank, slot in vector_data:` where `vector_data` holds **3-tuples** → `ValueError`
  on the first non-empty ligand. This branch *is* reachable via `OINParser().parse()` on
  any inline string with ligands, so it is not "unreachable code" — it is **reachable but
  always-raising**. RISK-2 calls it "cannot run today," which is accurate in outcome
  (raises) but not in reachability. Is there any caller that wraps this in try/except and
  depends on the raise? (grep shows the legacy `SMILESToXYZ`/TD-003 path; confirm no test
  asserts the `ValueError`.)
- **`SlotAssignment` export surface.** §5.4 defines it in `oin/inline.py`. `oin/parser.py`
  and `generation/oin_parser.py` both import from `oin.inline`. Confirm no import cycle:
  `oin/parser.py` already does `from .inline import OINInlineHandler` (`:16`), so adding
  `SlotAssignment` to that import is safe.

### What-If Scenarios
- **Pickle / ProcessPoolExecutor (Execution-Note b) — cleared.** The DG worker receives
  `args` = a dict of `str/int/list[tuple[float,…]]` only (`molassembler_adapter.py:2220`).
  Neither `OINVector` nor `SlotAssignment` is submitted across the boundary, and the
  template path runs in-process (`:2198`). A module-level `NamedTuple` would pickle fine
  regardless. **No pickle risk exists** — recommend deleting this stress item from §9.
- **Silent field-shift on a missed *keyword* constructor.** The Draft's safety argument
  (NamedTuple keeps `[0..2]` valid) protects readers. It does **not** protect a *writer*
  that constructs `SlotAssignment(lig_rank, atom_idx, slot)` positionally omitting winding
  — which is fine (defaults `None`) — but a writer that positionally passes a 4th value
  meaning something else would silently populate `winding`. Only the one producer at
  `inline.py:353` constructs it, so blast radius is 1; low risk, worth a one-line comment.
- **`test_oin_generation.py` collateral (unlisted).** `tests/unit/test_oin_generation.py:15-17`
  reads `parsed.vectors[0].atom_in_fragment_idx` etc. via **named** attributes. Adding a
  defaulted `winding` field to the `OINVector` dataclass does **not** break it (verified),
  but §5.3's blast-radius list omits it. It's a verify-unchanged site, not an edit — add
  it to the "verify unchanged" column for completeness.

### Points for Improvement
- In the C2 row of §5.1, change "already broken (…) → cannot run today" to "reachable but
  always-raising (`ValueError`)"; then the bugfix rationale is exact.
- Add `_build_connected_smiles` (`molassembler_adapter.py:1683`) **and**
  `test_oin_generation.py` to an explicit "verify-unchanged" checklist so the audit phase
  has a closed set.
- Note the second suffix stripper (`oin/inline.py:89`) in §5.1 as an out-of-scope
  winding-erasure site to forestall "map is incomplete" objections in `/hyper-resolve`.

---

## §6 Negative Constraints — Analysis

### Clarifying Questions
- "DO NOT normalize winding to bool/enum — store literal `'>'`/`'<'`/`None`." Consistent
  and good. But given `^` normalizes to `>` on the **generate** side, is the parse side
  expected to *also* accept `^` and store it verbatim (violating a clean two-symbol
  alphabet), or reject it? The constraint is silent on `^`.

### What-If Scenarios
- **"Winding available, not consumed" is only true where vectors exist.** §6 constraint
  "Winding is *available* to `_build_connected_smiles`/`_template_generate`, not
  *consumed*" is correct for the template path — but for `tmpl_vectors is None`, winding
  is not even *available* (no vectors). The negative constraint is technically satisfied
  (nothing consumes it) yet the positive guarantee (§1.2) is not met. These two framings
  should be reconciled.

### Points for Improvement
- Add: "DO NOT allow the `tmpl_vectors is None` branch to silently discard winding without
  it being an explicit, documented behavior."

---

## §7 Risks & Mitigation — Analysis

### Clarifying Questions
- The risk register omits the two conditional-drop paths (`tmpl_vectors is None`;
  `slot_idx >= len(tmpl_vectors)`). Should these be **RISK-5**?

### What-If Scenarios
- **RISK-1 "degrades gracefully" overstated.** True for index reads, false for unpack
  reads (see §4). The only real unpack consumers (C1 `:485`, C2 `:34`, C3 asserts) are
  all in the edit set, so nothing is silently *missed* — but the mitigation's stated
  mechanism (graceful positional degradation) is not what actually saves us; **exhaustive
  enumeration** is. Credit the correct mitigation.

### Points for Improvement
- Add **RISK-5 — winding silently dropped at `ParsedOIN` for template-less geometries**;
  mitigation = the new negative test + documented behavior from §3 improvements.
- Reword RISK-1 mitigation to lead with "the consumer set is closed and fully enumerated
  (§5.1)"; positional-index safety is a secondary net, not the primary guarantee.

---

## §8 Success Metrics — Analysis

### Clarifying Questions
- "Zero diff in generated XYZ for any existing generation test." Is there a *committed*
  golden-XYZ corpus this can be diffed against, or is behavioral inertness only asserted
  transitively via "tests stay green"? A byte-diff harness would make the inertness claim
  provable rather than assumed.

### What-If Scenarios
- A green suite does **not** prove byte-identical XYZ if no test asserts on exact
  coordinates. The `NamedTuple` change alters `ParsedOIN` construction ordering nowhere,
  so inertness is very likely — but "OK" from `unittest discover` is a weaker signal than
  the §8 wording implies.

### Points for Improvement
- Add a concrete inertness check: run the existing generation fixtures pre/post and
  `diff` the emitted XYZ (or hash it). Cheap, and it converts §8's strongest claim from
  assertion to evidence.
- Add the US-002-on-`NON`-geometry assertion (from §3) to the metrics list.

---

## §9 Execution Note — Analysis (meta)

- **(a) §5.1 completeness** → **CONFIRMED complete** for the `vector_data`/`parse_inline_string`
  contract (grep-verified: 3 non-test consumers + 1 test). One *adjacent* site
  (`oin/inline.py:89`, generate-side) and one *collateral* test
  (`test_oin_generation.py`) should be named for closure, but neither is a missed edit.
- **(b) NamedTuple pickling / ProcessPoolExecutor** → **CLEARED / non-issue.** The worker
  boundary carries a primitives-only dict; `OINVector`/`SlotAssignment` never pickle.
  Recommend removing this from the stress list.
- **(c) No winding-driven behavior leaks into placement** → **VALID, keep.** Reinforce
  with the byte-diff inertness harness (§8).

**Net new issue the Execution Note did not anticipate:** the winding guarantee at the
`ParsedOIN` layer is **conditional on `tmpl_vectors is not None`** — the single most
important thing to resolve before this Draft is compiled, because it silently fails for
the eta/haptic ligand family Phase 3 depends on.

---

## Top Findings (ranked for `/hyper-resolve`)

1. **[Correctness gap] Winding is silently dropped at `ParsedOIN` when `geo_code` has no
   `TEMPLATES` entry** (`oin_parser.py:484` gate) and when `slot_idx >= len(tmpl_vectors)`
   (`:485`). Ferrocene passes only because `LIN` is templated. Resolve scope + add a
   template-less negative test. *(→ new RISK-5, US-002 second fixture, §3 acceptance.)*
2. **[Semantic ambiguity] The suffix is the *heading* marker, and its alphabet includes
   `^`** (generate side normalizes `^`→`>`, `oin/inline.py:245`). PRD regex `([><])?`
   silently ignores `^`. Decide: two-symbol alphabet (reject `^`) or three (`([><^])?`
   + normalize). Name the second stripper at `oin/inline.py:89`.
3. **[Precision] RISK-1 / US-003 "graceful positional degradation"** protects index reads,
   not unpack reads; the real safety is the closed, enumerated consumer set. Reword.
4. **[De-risk] Execution-Note (b) pickle concern is a non-issue** (primitives-only worker
   boundary). Remove it and spend the attention on #1.
5. **[Evidence] §8 inertness** should be proven with a pre/post XYZ byte-diff, not inferred
   from a green suite.

---

**Final Action:** Run `/hyper-resolve` to triage these findings — priority to Finding #1
(conditional winding drop) and #2 (`^`/heading alphabet), which change the acceptance
criteria; #3–#5 are hardening and wording.
