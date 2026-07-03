# SuperPRD — Stereo Roadmap Phase 1: "Preserve the Signal" (Winding Plumbing)

## Metadata
- **Project Name**: OIN-SMILES — Generation-side stereochemistry, Phase 1
- **Version**: 1.0.0 (Compiled — post `/hyper-resolve`)
- **Status**: Ready for `/hyper-execute`
- **Owner**: Architect Agent (Fable) / Thomas Mustard
- **Roadmap ref**: `spec/worklog/ROADMAP-stereo.md` § Phase 1
- **Depends on**: TASK-10 (Phase 0 diagnostics, DONE 2026-07-03)
- **Provenance**: Draft_PRD.md (v0.1.0) → RedTeam_Report.md (2026-07-03) → this compilation
- **Child MiniPRD**: `spec/compiled/MiniPRD_WindingPlumbing_Phase1.md`

---

## 1. Introduction & Goals

### 1.1 Problem Statement
The OIN→XYZ generation direction silently discards eta-ligand **winding direction**
(`{n>}` = CW heading, `{n<}` = CCW heading). Phase 0 proved this empirically
(`tests/unit/test_stereo_roundtrip_diagnostics.py::test_haptic_face_winding`):
flipping ferrocene ring-0 from `{0>}` to `{0<}` produces **byte-identical** 3D output,
because the parse layer never sees the direction suffix.

Root cause is narrow and mechanical: `SLOT_REGEX` (`src/oinsmiles/oin/inline.py:44`,
`\{(\d+)[><]?\}`) *matches* the `>`/`<` suffix but does not **capture** it, so
`parse_inline_string` strips it during the same pass that strips the slot marker. The
signal is destroyed at the first parse step; every downstream layer is blameless — it
was never given the data.

**Known, out-of-scope second erasure site (Red Team #2):** the generate-side
V2.4→inline stripper at `oin/inline.py:89`
(`int(slot_str.replace('^','').replace('>','').replace('<',''))`) also erases winding.
It is **not** part of the OIN→XYZ `vector_data` contract this phase repairs and is left
untouched — recorded here so a future reader does not mistake the §5 map for incomplete.

### 1.2 Solution Overview
Phase 1 is **plumbing only — no geometry/placement behavior changes.** Capture the
winding suffix in the regex, carry it through `parse_inline_string`'s returned data on a
new `winding` field, and surface it on `ParsedOIN` so the adapter *could* read it in
Phase 3. The value is made available and **provably preserved for all geometries**;
nothing yet consumes it to alter output.

**The §1.2 preservation guarantee is unconditional** (Red Team Finding #1 resolution):
winding survives to `ParsedOIN` for **every** geometry that carries a heading marker,
**including template-less `NON`/eta ligands** — the exact ligand family Phase 3 targets.
This is delivered via a dedicated `ParsedOIN.winding_by_slot` dict populated on all
parse paths, so the guarantee does not depend on the geometry having a `TEMPLATES` entry.

### 1.3 Target Audience
Internal: the OIN→XYZ generation pipeline and its maintainers. Phase 1 unblocks Phase 3
(haptic-face control), the first phase permitted to *use* winding to change geometry.

---

## 2. Confidence Mandate
**Confidence Score**: 10/10 (post-resolution)

- **US-001 (capture at `parse_inline_string`)**: 10/10 — mechanical regex change,
  fully enumerated.
- **US-002 (threading to `ParsedOIN`)**: 10/10 — the universal `winding_by_slot`
  channel removes the previous `tmpl_vectors is not None` conditionality that the Red
  Team flagged (Finding #1). Inertness on the `NON` path is structural (vectors
  emission unchanged), not merely tested.
- Residual risk formerly assigned to the `oin/parser.py:34` broken unpack (RISK-2) and
  the conditional winding drop (Red Team #1) are both resolved below.

All Red Team findings (#1–#5) carry a documented decision — see §7 and §9.

---

## 3. Scope

### 3.1 In-Scope
1. **Capture winding in `SLOT_REGEX`** (`oin/inline.py:44`): `\{(\d+)[><]?\}` →
   `\{(\d+)([><^])?\}`. Group 2 ∈ {`'>'`, `'<'`, `'^'`, `None`}; on capture, normalize
   `'^'` → `'>'` to match the generate side (`oin/inline.py:245`). Stored value is
   therefore always `'>'`, `'<'`, or `None`.
2. **New element type `SlotAssignment`** — a `typing.NamedTuple` co-located with its
   producer in `oin/inline.py`, exported for consumers:
   `SlotAssignment(lig_rank: int, atom_idx: int, slot: int, winding: Optional[str] = None)`.
3. **`parse_inline_string` returns `List[SlotAssignment]`** for `vector_data`,
   populating `winding` from the (normalized) `SLOT_REGEX` group 2.
4. **Universal winding channel on `ParsedOIN`** (Red Team #1 resolution):
   `ParsedOIN` gains `winding_by_slot: Dict[int, Optional[str]]`, populated on **every**
   parse path — template *and* template-less (`tmpl_vectors is None` / `NON` / eta).
   For the template path, `OINVector` additionally gains a defaulted
   `winding: Optional[str] = None` copied from `sa.winding`.
5. **Migrate all consumers** of the `vector_data` element shape (see §5.1) mechanically.
6. **Tests** — update the exact-equality asserts in `tests/unit/test_inline.py`; add unit
   tests proving winding capture (`>`, `<`, `^`→`>`, `None`), arrival on
   `ParsedOIN.vectors[].winding` for a template geometry (ferrocene/LIN), **and** arrival
   on `ParsedOIN.winding_by_slot` for a template-less (`NON`) geometry.
7. **Byte-diff inertness harness** (Red Team #5): run the existing generation fixtures
   pre/post and assert byte-identical emitted XYZ, as a hard gate on §6 inertness.

### 3.2 Out-of-Scope (explicit)
- **Any placement/geometry/RMSD change.** No builder reads `winding` this phase.
- **Flipping `test_haptic_face_winding` to passing.** It is the **Phase 3** target and
  MUST remain `@unittest.expectedFailure` at the end of Phase 1. It compares re-encoded
  3D output, not `ParsedOIN`, so it correctly stays red even though `ParsedOIN` winding
  is now captured — the two are decoupled and no one should prematurely flip it.
- **The legacy V2.4 sidecar `w:` tag parse path** (`oin_parser.py:524+`, `oin/parser.py`
  w-tag reconstruction). Winding stays `None` on that path.
- **The generate-side stripper `oin/inline.py:89`** — a known winding-erasure site left
  untouched (§1.1).
- **P/N `@/@@` CIP enforcement, haptic-face computation, Zone-A stereo** — Phases 2/3/4.
- **`molassembler_adapter._build_connected_smiles` / `_template_generate` /
  `_stitch_eta_fragment` logic** — verified *unchanged*, not modified.

---

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
| :-- | :-- | :-- | :-- |
| US-001 | As the parse layer, I capture the winding suffix so it survives the first parse step. | 1. `parse_inline_string("[Fe_LIN].[cH]{0>}1...")` yields a `SlotAssignment` with `winding == '>'`.<br>2. Plain `{0}` → `winding is None`.<br>3. `{0<}` → `winding == '<'`.<br>4. `{0^}` → `winding == '>'` (normalized). | High |
| US-002 | As a downstream generator, I can read winding off `ParsedOIN` without re-parsing — for **every** geometry. | 1. `OINParser().parse(ferrocene_oin)` → `ParsedOIN` whose `.vectors` includes an `OINVector` with `winding == '>'` for the heading atom; non-heading vectors `winding is None`.<br>2. **Template-less `NON` geometry** with `{0>}` → `ParsedOIN.winding_by_slot[0] == '>'` even though `.vectors == []`.<br>3. `winding_by_slot` is populated on all paths. | High |
| US-003 | As a maintainer, I get a self-documenting element type, not a bare tuple whose 4th slot is a mystery. | 1. `SlotAssignment` is a `NamedTuple` with named fields incl. `winding`.<br>2. Positional **index** reads `sa[0..2]` stay correct; positional **unpack** (`a,b,c = sa`) fails fast on arity (intended — see RISK-1).<br>3. `sa.winding` defaults to `None` when omitted at construction. | Medium |
| US-004 | As the test suite, I stay green and the Phase-3 gap stays visibly red. | 1. `discover tests/unit` → OK, count up by the new plumbing tests.<br>2. `test_haptic_face_winding` remains `expectedFailure` (unchanged).<br>3. `discover tests` (root) → OK.<br>4. Byte-diff harness reports zero XYZ delta across generation fixtures. | High |

---

## 5. Technical Specifications

### 5.1 Consumer Map (the blast radius that matters — Red Team CONFIRMED complete)

`vector_data` element shape is the load-bearing contract. Every site that constructs,
unpacks, or asserts on it (grep-verified closed set: 1 producer + 3 consumers + tests):

| # | Site | Today | After |
| :-- | :-- | :-- | :-- |
| P | `oin/inline.py:353` (producer) | `vector_data.append((lig_rank, atom_idx, slot))` | `append(SlotAssignment(lig_rank, atom_idx, slot, winding))` |
| C1 | `generation/oin_parser.py:485` | `for lig_rank, atom_in_fragment_idx, slot_idx in vector_data:` (inside `if tmpl_vectors is not None:`) | `for sa in vector_data:` → use `sa.lig_rank/.atom_idx/.slot/.winding`; set `OINVector(..., winding=sa.winding)`. **Additionally**, populate `winding_by_slot[sa.slot] = sa.winding` for **all** `sa` **outside** the `tmpl_vectors` gate so `NON`/eta paths surface winding. |
| C2 | `oin/parser.py:34` | `for rank, slot in vector_data:` — **reachable but always-raising** (2-way unpack of 3-tuples → `ValueError`); TD-003 `SMILESToXYZ` path | `for sa in vector_data:` → use `sa.lig_rank`, `sa.slot`. **Bugfix only**; winding NOT propagated (legacy w-tag reconstruction, out of scope). |
| C3 | `tests/unit/test_inline.py:21,30` | `assertEqual(vectors[0], (1, 0, 0))` — **False** vs a 4-field NamedTuple, so update is mandatory | `assertEqual(vectors[0], SlotAssignment(1, 0, 0, None))` (== `(1,0,0,None)`) |

**Verify-unchanged sites (closed set for the audit phase; NOT edits):**
- `molassembler_adapter._build_connected_smiles` (`:1683`) — reads
  `OINVector.atom_in_fragment_idx` only; a defaulted `OINVector.winding` does not touch
  it. `winding_by_slot` is a new field it never iterates.
- `tests/unit/test_oin_generation.py:15-17` — reads `parsed.vectors[0]` named attrs;
  adding a defaulted field does not break it.
- `oin/inline.py:89` (generate-side stripper) — out of scope, unchanged.

**Why a `NamedTuple` and not a plain 4-tuple or sidecar dict** (resolved fork): a plain
4-tuple forces positional unpack with an undocumented 4th slot; a sidecar dict for the
*element* adds a second per-element channel to keep in sync. The `NamedTuple` is the
minimal self-documenting change that keeps index reads valid. (Note: the `NON`-path
`winding_by_slot` dict is a deliberate **`ParsedOIN`-level** channel, chosen precisely
because emitting placeholder `OINVector`s on that path would perturb the adapter — see
§7 RISK-5.)

### 5.2 Data flow (unchanged topology, new fields riding along)
```
OIN string
  → OINInlineHandler.parse_inline_string      [SLOT_REGEX captures + normalizes winding]
      returns (smiles, geometry, List[SlotAssignment])   # winding on each element
  → generation.OINParser.parse
      template path:  builds ParsedOIN.vectors: List[OINVector]  # winding on each
      ALL paths:      populates ParsedOIN.winding_by_slot         # universal channel
  → molassembler_adapter._build_connected_smiles / _template_generate
      reads atom_in_fragment_idx only          # winding present but UNUSED (Phase 3)
```

### 5.3 API Contracts / Schema
```python
# src/oinsmiles/oin/inline.py
class SlotAssignment(NamedTuple):
    lig_rank: int
    atom_idx: int
    slot: int
    winding: Optional[str] = None   # '>' (CW) | '<' (CCW) | None  ('^' normalized to '>')

# SLOT_REGEX: r"\{(\d+)([><^])?\}"   # group(1)=slot, group(2)=winding-or-None (pre-normalize)

@staticmethod
def parse_inline_string(inline_string: str
    ) -> Tuple[str, str, List["SlotAssignment"]]: ...

# src/oinsmiles/generation/oin_parser.py
@dataclass
class OINVector:
    atom_idx: int
    vector: Tuple[float, float, float]
    fragment_idx: int
    atom_in_fragment_idx: int
    winding: Optional[str] = None            # NEW, defaulted; populated on template path

@dataclass
class ParsedOIN:
    # ... existing fields ...
    winding_by_slot: Dict[int, Optional[str]] = field(default_factory=dict)  # NEW, ALL paths
```

### 5.4 Dependencies
No new dependencies. Standard library `typing.NamedTuple` / `dataclasses`; existing
RDKit / molassembler untouched.

---

## 6. Negative Constraints (The "Do NOTs")
- **DO NOT** change any 3D coordinate, RMSD, or placement outcome. Phase 1 is
  observationally inert on generation output (enforced by the §7/§8 byte-diff harness).
- **DO NOT** emit `OINVector`s on the `tmpl_vectors is None` path — surface winding there
  via `winding_by_slot` only, so the adapter's `vectors` iteration is provably untouched.
- **DO NOT** un-mark `test_haptic_face_winding`; it stays `expectedFailure`.
- **DO NOT** normalize winding to a bool/enum. Store literal `'>'`/`'<'`/`None`
  (`'^'` collapses to `'>'` on capture only).
- **DO NOT** propagate winding through the legacy V2.4 sidecar `w:` path or the
  generate-side stripper (`oin/inline.py:89`). Leave `None` / untouched.
- **DO NOT** allow any `winding`-carrying path to silently discard the value without it
  being explicit and covered by a test (Red Team #1 / #6).
- **DO NOT** modify `_build_connected_smiles` / `_template_generate` /
  `_stitch_eta_fragment` logic. Winding is *available* to them, not *consumed*.
- **DO NOT** touch `METAL_REGEX` or the `@desc` back-compat tolerance (D-3, prior work).

---

## 7. Risks & Mitigation
- **RISK-1 — a `vector_data` consumer is missed.** → **Mitigation:** the complete
  consumer set is a **closed, grep-enumerated** list (§5.1) — that is the primary
  guarantee. Positional-**index** safety (`sa[0..2]` stays correct on a NamedTuple) is a
  *secondary* net for slice/index readers only; positional-**unpack** readers
  (`a,b,c = sa`) fail fast on arity, which is desirable, not "graceful degradation."
- **RISK-2 — C2 (`oin/parser.py:34`) is reachable but always-raising.** → **Mitigation:**
  it is the incomplete `SMILESToXYZ` (TD-003) path; a 2-way unpack of 3-tuples raises
  `ValueError` on any inline string with ligands. No test asserts that raise. Fixing to
  named iteration is strictly an improvement; do not expand its behavior.
- **RISK-3 — accidental Phase-3 behavior leak.** → **Mitigation:** §3.2 + §6 make the
  `expectedFailure` invariant an explicit acceptance criterion (US-004.2), and the
  byte-diff harness (§8) would catch any placement change.
- **RISK-4 — `test_inline.py` equality asserts** now compare against a `NamedTuple`.
  → **Mitigation:** `SlotAssignment(1,0,0,None) == (1,0,0,None)` holds; `== (1,0,0)` is
  False (length differs), so updating to 4-length literals is mandatory.
- **RISK-5 — winding silently dropped at `ParsedOIN` for template-less geometries**
  (the Red Team's headline gap). → **Mitigation:** `winding_by_slot` is populated on
  **all** paths, *outside* the `tmpl_vectors is not None` gate, so `NON`/eta ligands
  surface winding without emitting placeholder `OINVector`s (which would perturb the
  adapter). A dedicated negative/positive test on a `NON` geometry pins the behavior.
  Also covers the `slot_idx >= len(tmpl_vectors)` overflow drop, since `winding_by_slot`
  is keyed by slot independent of template length.

---

## 8. Success Metrics
- `uv run python -m unittest discover tests/unit` → **OK**, `expected failures=1`
  preserved (`test_haptic_face_winding`), count up by the new plumbing tests.
- `uv run python -m unittest discover tests` (root) → **OK**.
- New unit assertions: `winding == '>'` / `'<'` / `'^'→'>'` / `None` captured by
  `parse_inline_string`; the same value on `ParsedOIN.vectors[].winding` (template path,
  US-002.1) **and** on `ParsedOIN.winding_by_slot` for a `NON` geometry (US-002.2).
- **Byte-diff inertness gate:** run existing generation fixtures pre/post and `diff`
  (or hash) the emitted XYZ — **zero delta**. Converts the inertness claim from assertion
  to evidence (Red Team #5).

---

## 9. Red Team Disposition Log (all findings resolved)
| RT # | Finding | Decision |
| :-- | :-- | :-- |
| 1 | Winding silently dropped at `ParsedOIN` when geometry has no `TEMPLATES` entry (`oin_parser.py:484` gate) and on `slot_idx` overflow. | **Resolved — universal channel.** Add `ParsedOIN.winding_by_slot` populated on all paths, outside the template gate; keyed by slot (overflow-immune). NON/eta now covered. New test pins it. |
| 2 | Suffix is the *heading* marker; alphabet includes `^` (generate side normalizes `^`→`>`). Regex `([><])?` ignores `^`. | **Resolved — three-symbol.** Regex `([><^])?`; normalize `^`→`>` on capture. Named the second stripper `oin/inline.py:89` as out-of-scope. |
| 3 | RISK-1 / US-003 "graceful positional degradation" overstated (protects index reads, not unpack). | **Resolved — reworded.** Primary safety = closed enumerated consumer set; index-safe / unpack-fail-fast framing (US-003.2, RISK-1). |
| 4 | Execution-Note (b) pickle / ProcessPoolExecutor concern is a non-issue (worker boundary carries a primitives-only dict; `OINVector`/`SlotAssignment` never pickle). | **Resolved — deleted** as a stress item; noted as CLEARED. |
| 5 | §8 inertness inferred from a green suite, not proven. | **Resolved — byte-diff harness** added as a hard gate (§7/§8). |
