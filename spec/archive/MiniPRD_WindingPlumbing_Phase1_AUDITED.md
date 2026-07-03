# MiniPRD: Winding Plumbing — Stereo Roadmap Phase 1

**Hypergraph Node ID:** `module_winding_plumbing_phase1`
**Parent Node:** `system_oin_generation`
**Parent SuperPRD:** `spec/compiled/SuperPRD_StereoPhase1_Winding.md`
**Roadmap Reference:** `spec/worklog/ROADMAP-stereo.md` § Phase 1
**Depends On:** TASK-10 (Phase 0 diagnostics, DONE 2026-07-03)
**Execution Order:** 1 of 1 (single-MiniPRD effort)
**Estimated Effort:** 3–5 hours
**Priority:** P1

---

## 1. The Confidence Mandate

**Confidence Score:** 10/10

**Rationale:**
- ✅ Consumer set is a **closed, grep-verified** list (1 producer, 3 consumers, 2 tests) —
  Red Team CONFIRMED completeness against a repo-wide grep.
- ✅ The one design fork (tuple representation) is resolved: `NamedTuple SlotAssignment`.
- ✅ The Red Team's material gap (winding dropped for template-less geometries) is resolved
  by a `ParsedOIN.winding_by_slot` dict populated outside the `tmpl_vectors` gate — this
  keeps `NON`/eta covered **and** keeps the adapter's `vectors` iteration untouched
  (structural inertness, not merely tested).
- ✅ The `^` alphabet ambiguity is resolved: regex `([><^])?` + normalize `^`→`>`.
- ✅ Inertness is provable via a pre/post XYZ byte-diff harness, not inferred.

---

## 2. Atomic User Stories

- **US-001:** As the parse layer, I capture the winding suffix (`>`/`<`/`^`) in
  `SLOT_REGEX` so it survives the first parse step, normalizing `^`→`>`.
- **US-002:** As a downstream generator, I read winding off `ParsedOIN` for **every**
  geometry — `OINVector.winding` on template paths, `ParsedOIN.winding_by_slot` on
  template-less (`NON`/eta) paths — without re-parsing the OIN string.
- **US-003:** As a maintainer, I get a self-documenting `SlotAssignment` NamedTuple whose
  index reads stay valid and whose unpack reads fail fast on arity.
- **US-004:** As the test suite, I stay green, add coverage for winding capture + both
  threading legs, keep `test_haptic_face_winding` `expectedFailure`, and prove generation
  output is byte-identical.

---

## 3. Implementation Plan (Task List)

- [ ] **Task 1 (regex capture + normalize):** In `src/oinsmiles/oin/inline.py:44`, change
  `SLOT_REGEX` from `\{(\d+)[><]?\}` to `\{(\d+)([><^])?\}`. Where the slot is parsed
  (`:350` finditer region), read `group(2)`; if it is `'^'`, set it to `'>'`; else keep
  `'>'`/`'<'`/`None`.
- [ ] **Task 2 (SlotAssignment type):** Define
  `class SlotAssignment(NamedTuple): lig_rank: int; atom_idx: int; slot: int; winding: Optional[str] = None`
  at module scope in `oin/inline.py` (importable). Add a one-line comment: only the
  producer at `:353` constructs it; a positional 4th arg means winding.
- [ ] **Task 3 (producer):** At `oin/inline.py:353`, replace the 3-tuple append with
  `vector_data.append(SlotAssignment(lig_rank, atom_idx, slot, winding))`, sourcing
  `winding` from Task 1. Update `parse_inline_string`'s return annotation to
  `Tuple[str, str, List[SlotAssignment]]`.
- [ ] **Task 4 (OINVector field):** In `generation/oin_parser.py`, add
  `winding: Optional[str] = None` (defaulted, last) to the `OINVector` dataclass.
- [ ] **Task 5 (ParsedOIN universal channel):** Add
  `winding_by_slot: Dict[int, Optional[str]] = field(default_factory=dict)` to `ParsedOIN`.
- [ ] **Task 6 (C1 — threading, `oin_parser.py:485`):** Rewrite the loop to
  `for sa in vector_data:` using named attrs. **Outside** the `if tmpl_vectors is not None:`
  gate, populate `winding_by_slot[sa.slot] = sa.winding` for every `sa`. **Inside** the
  gate (unchanged control flow), set `OINVector(..., winding=sa.winding)`. Do **not** emit
  `OINVector`s on the `else`/`None` branch.
- [ ] **Task 7 (C2 bugfix — `oin/parser.py:34`):** Replace the broken
  `for rank, slot in vector_data:` with `for sa in vector_data:` using `sa.lig_rank`,
  `sa.slot`. Add `from .inline import SlotAssignment` if a type reference is needed
  (import cycle already safe — `oin/parser.py:16` imports from `.inline`). Do **not**
  propagate winding here (legacy w-tag path, out of scope).
- [ ] **Task 8 (C3 test literals — `tests/unit/test_inline.py:21,30`):** Update
  `assertEqual(vectors[0], (1, 0, 0))` → `assertEqual(vectors[0], SlotAssignment(1, 0, 0, None))`.
- [ ] **Task 9 (new capture tests):** Add unit tests asserting `parse_inline_string`
  yields `winding == '>'` for `{0>}`, `'<'` for `{0<}`, `'>'` for `{0^}` (normalized),
  and `None` for plain `{0}`.
- [ ] **Task 10 (new threading tests):** (a) ferrocene/`LIN` → `ParsedOIN.vectors[]`
  heading `OINVector.winding == '>'`, non-heading `None`. (b) a **template-less `NON`**
  OIN with `{0>}` → `ParsedOIN.winding_by_slot[0] == '>'` while `ParsedOIN.vectors == []`.
- [ ] **Task 11 (inertness harness):** Add a test/harness that runs the existing
  generation fixtures pre/post and asserts byte-identical emitted XYZ (hash or `diff`).
- [ ] **Task 12 (green + red gates):** `uv run python -m unittest discover tests/unit`
  and `discover tests` both **OK**; `test_haptic_face_winding` still `expectedFailure`.

---

## 4. The Negative Space (Constraints)

- **DO NOT** change any 3D coordinate, RMSD, or placement outcome (byte-diff gate, Task 11).
- **DO NOT** emit `OINVector`s on the `tmpl_vectors is None` path — winding rides
  `winding_by_slot` there so the adapter's `vectors` iteration is provably untouched.
- **DO NOT** un-mark `test_haptic_face_winding` (`expectedFailure` is Phase 3's gate).
- **DO NOT** normalize winding to bool/enum; store literal `'>'`/`'<'`/`None` (`'^'`→`'>'`
  on capture only).
- **DO NOT** propagate winding through the legacy V2.4 `w:` path or the generate-side
  stripper `oin/inline.py:89` (both out of scope).
- **DO NOT** modify `_build_connected_smiles` / `_template_generate` /
  `_stitch_eta_fragment` — verify-unchanged sites only.
- **DO NOT** expand C2's behavior beyond the unpack bugfix.

---

## 5. Integration Tests & Verification

- **Test 1 (Deterministic — capture):** `parse_inline_string("[Fe_LIN].[cH]{0>}1...")`
  → the `SlotAssignment` for the heading atom has `winding == '>'`; a plain `{0}` element
  has `winding is None`; `{0<}` → `'<'`; `{0^}` → `'>'`.
- **Test 2 (Deterministic — template threading):** `OINParser().parse(ferrocene_oin)`
  → `ParsedOIN.vectors` contains an `OINVector` with `winding == '>'` for the heading
  atom, `None` for non-heading vectors.
- **Test 3 (Deterministic — NON threading, the Red Team gap):** parse a `NON`-geometry
  OIN carrying `{0>}` → `ParsedOIN.winding_by_slot[0] == '>'` **and**
  `ParsedOIN.vectors == []`. Pins that the universal channel covers the template-less path.
- **Test 4 (Inertness):** emitted XYZ for every existing generation fixture is
  byte-identical pre/post this change.
- **Test 5 (Suite invariants):** `discover tests/unit` and `discover tests` both OK;
  `expected failures=1` (`test_haptic_face_winding`) preserved; test count increased.
