# Draft PRD — Stereo Roadmap Phase 1: "Preserve the Signal" (Winding Plumbing)

## Metadata
- **Project Name**: OIN-SMILES — Generation-side stereochemistry, Phase 1
- **Version**: 0.1.0 (Draft)
- **Status**: Draft
- **Owner**: Architect Agent (Fable) / Thomas Mustard
- **Roadmap ref**: `spec/worklog/ROADMAP-stereo.md` § Phase 1
- **Depends on**: TASK-10 (Phase 0 diagnostics, DONE 2026-07-03)

## 1. Introduction & Goals

### 1.1 Problem Statement
The OIN→XYZ generation direction silently discards eta-ligand **winding direction**
(`{n>}` = CW heading, `{n<}` = CCW heading). Phase 0 proved this empirically
(`tests/unit/test_stereo_roundtrip_diagnostics.py::test_haptic_face_winding`):
flipping ferrocene ring-0 from `{0>}` to `{0<}` produces **byte-identical** 3D
output, because the parse layer never sees the direction suffix.

Root cause is narrow and mechanical: `SLOT_REGEX`
(`src/oinsmiles/oin/inline.py:44`, `\{(\d+)[><]?\}`) *matches* the `>`/`<`
suffix but does not **capture** it, so `parse_inline_string` strips it during
the same pass that strips the slot marker. The signal is destroyed at the first
parse step; every downstream layer (`OINParser.parse` → `ParsedOIN` → the
molassembler adapter) is therefore blameless — it was never given the data.

### 1.2 Solution Overview
Phase 1 is **plumbing only — no geometry/placement behavior changes.** Capture
the winding suffix in the regex, carry it through `parse_inline_string`'s
returned data on a new `winding` field, and thread it into `ParsedOIN.vectors`
so the adapter *could* read it in Phase 3. The value is made available and
provably preserved; nothing yet consumes it to alter output.

The design fork — how to change the `vector_data` element type without
scattering breakage — is resolved: `vector_data` becomes a list of a new
`typing.NamedTuple` (`SlotAssignment`) with `winding: Optional[str] = None`.
This gives named `.winding` access, keeps positional/index reads working
(a `NamedTuple` *is* a tuple), and makes every consumer's migration a
mechanical, reviewable edit.

### 1.3 Target Audience
Internal: the OIN→XYZ generation pipeline and its maintainers. Phase 1 unblocks
Phase 3 (haptic-face control), which is the first phase permitted to *use*
winding to change geometry.

## 2. Confidence Mandate
**Confidence Score**: 9/10

The scope is tightly bounded and the full consumer set was enumerated from the
codebase (see §5). The one genuine design decision (tuple representation) is
resolved. The residual 1 point is the pre-existing broken 2-tuple unpack in
`oin/parser.py:34`, which the contract change forces us to touch — see RISK-2.

**Clarifying Questions** (all resolved during interview):
- [x] Tuple contract change strategy → **`NamedTuple` with `winding=None` default** (user-selected).
- [x] Winding value representation → literal suffix char `'>'` / `'<'` / `None`
  (architect decision; faithful to source, matches the generate-side and the
  Phase-0 diagnostic strings — no premature normalization to enum/bool).
- [x] Scope of threading → the **inline** parse path only (`SLOT_REGEX`). The
  legacy V2.4 sidecar `w:` tag path (`oin_parser.py:524+`) is out of scope.

## 3. Scope

### 3.1 In-Scope
1. **Capture winding in `SLOT_REGEX`** — `src/oinsmiles/oin/inline.py:44`:
   `\{(\d+)[><]?\}` → `\{(\d+)([><])?\}`. Group 2 = `'>'`, `'<'`, or `None`.
2. **New element type `SlotAssignment`** — a `typing.NamedTuple` co-located with
   its producer in `oin/inline.py`, exported for consumers:
   `SlotAssignment(lig_rank: int, atom_idx: int, slot: int, winding: Optional[str] = None)`.
3. **`parse_inline_string` returns `List[SlotAssignment]`** for `vector_data`,
   populating `winding` from `SLOT_REGEX` group 2 (default `None` when the
   marker is a plain `{n}`).
4. **Thread winding into `ParsedOIN`** — `generation/oin_parser.py`:
   `OINVector` gains `winding: Optional[str] = None`; the inline branch
   (`parse`, ~line 485) copies `sa.winding` into each `OINVector`.
5. **Migrate all consumers** of the element shape (see §5.1) to the
   `SlotAssignment` form, mechanically.
6. **Tests** — update the two exact-equality asserts in
   `tests/unit/test_inline.py`; add new unit tests proving winding capture
   (`>`, `<`, `None`) and its arrival on `ParsedOIN.vectors[].winding`.

### 3.2 Out-of-Scope (explicit)
- **Any placement/geometry/RMSD change.** No builder reads `winding` this phase.
- **Flipping `test_haptic_face_winding` to passing.** That test is the **Phase 3**
  target and MUST remain `@unittest.expectedFailure` at the end of Phase 1.
- **The legacy V2.4 sidecar `w:` tag parse path** (`oin_parser.py:524+`,
  `oin/parser.py`'s tag reconstruction). Winding stays `None` on that path.
- **P/N `@/@@` CIP enforcement, haptic-face computation, Zone-A stereo** —
  Phases 2/3/4 respectively.
- **`molassembler_adapter._build_connected_smiles` logic** — it reads
  `OINVector.atom_in_fragment_idx`, never `winding`; it is verified *unchanged*,
  not modified.

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
| :-- | :-- | :-- | :-- |
| US-001 | As the parse layer, I capture the winding suffix so it survives the first parse step. | 1. `parse_inline_string("[Fe_LIN].[cH]{0>}1...")` yields a `SlotAssignment` with `winding == '>'` for the heading marker.<br>2. A plain `{0}` marker yields `winding is None`.<br>3. `{0<}` yields `winding == '<'`. | High |
| US-002 | As a downstream generator, I can read winding off `ParsedOIN` without re-parsing the OIN string. | 1. `OINParser().parse(ferrocene_oin)` returns `ParsedOIN` whose `.vectors` includes an `OINVector` with `winding == '>'` for the heading atom.<br>2. Non-heading vectors carry `winding is None`. | High |
| US-003 | As a maintainer, I get a self-documenting element type, not a bare tuple whose 4th slot is a mystery. | 1. `SlotAssignment` is a `NamedTuple` with named fields incl. `winding`.<br>2. `sa[0..2]` still returns `lig_rank/atom_idx/slot` (positional back-compat).<br>3. `sa.winding` defaults to `None` when omitted at construction. | Medium |
| US-004 | As the test suite, I stay green and the Phase-3 gap stays visibly red. | 1. `discover tests/unit` → OK, and count increases by the new plumbing tests.<br>2. `test_haptic_face_winding` remains `expectedFailure` (unchanged).<br>3. `discover tests` (root) → OK. | High |

## 5. Technical Specifications (The Blueprint)

### 5.1 Architecture & Consumer Map (the blast radius that matters)

`vector_data` element shape is the load-bearing contract. Every site that
constructs, unpacks, or asserts on it:

| # | Site | Today | After |
| :-- | :-- | :-- | :-- |
| P | `oin/inline.py:353` (producer) | `vector_data.append((lig_rank, atom_idx, slot))` | `append(SlotAssignment(lig_rank, atom_idx, slot, winding))` |
| C1 | `generation/oin_parser.py:485` | `for lig_rank, atom_in_fragment_idx, slot_idx in vector_data:` | `for sa in vector_data:` → use `sa.lig_rank/.atom_idx/.slot/.winding`; set `OINVector(..., winding=sa.winding)` |
| C2 | `oin/parser.py:34` | `for rank, slot in vector_data:` — **already broken** (2-way unpack of 3-tuples → `ValueError`); dead TD-003 `SMILESToXYZ` path | `for sa in vector_data:` → use `sa.lig_rank`, `sa.slot`. **Bugfix only**; winding NOT propagated (legacy w-tag reconstruction, out of scope). |
| C3 | `tests/unit/test_inline.py:21,30` | `assertEqual(vectors[0], (1, 0, 0))` | `assertEqual(vectors[0], SlotAssignment(1, 0, 0, None))` (or equivalently the 4-tuple `(1, 0, 0, None)`) |

**Why a `NamedTuple` and not a plain 4-tuple or a sidecar dict** (resolved fork):
a plain 4-tuple forces every reader onto positional unpack with an undocumented
4th slot; a sidecar dict adds a second data channel callers must keep in sync
(easy to desync). The `NamedTuple` is the minimal change that is also
self-documenting and keeps index-based reads valid, so any consumer missed in
review still gets correct `[0..2]` values rather than a silent shift.

**`_build_connected_smiles` (`molassembler_adapter.py:1683`)** groups `ParsedOIN.vectors`
by `fragment_idx` and reads only `atom_in_fragment_idx`. Adding `OINVector.winding`
(a defaulted field) does not touch it. It is an explicit **no-change / verify-unchanged**
site, not an edit site.

### 5.2 Data flow (unchanged topology, new field riding along)
```
OIN string
  → OINInlineHandler.parse_inline_string        [SLOT_REGEX now captures winding]
      returns (smiles, geometry, List[SlotAssignment])   # winding on each element
  → generation.OINParser.parse (inline branch)
      builds ParsedOIN.vectors: List[OINVector]  # winding copied onto each OINVector
  → molassembler_adapter._build_connected_smiles / _template_generate
      reads atom_in_fragment_idx only            # winding present but UNUSED (Phase 3)
```

### 5.3 System Graph Blast Radius
`spec/compiled/architecture.yml` tracks the HACF framework, not this project
(see project memory), so the authoritative blast radius is the source set:
- `src/oinsmiles/oin/inline.py` — `SLOT_REGEX`, `SlotAssignment` (new), `parse_inline_string` (signature/return element type).
- `src/oinsmiles/generation/oin_parser.py` — `OINVector` (new field), `OINParser.parse` inline branch.
- `src/oinsmiles/oin/parser.py` — C2 bugfix (named-attr iteration).
- `tests/unit/test_inline.py` — updated + new assertions.
- `tests/unit/` — new plumbing test(s) for the `ParsedOIN`-threading path (US-002).

### 5.4 API Contracts / Schema
```python
# src/oinsmiles/oin/inline.py
class SlotAssignment(NamedTuple):
    lig_rank: int
    atom_idx: int
    slot: int
    winding: Optional[str] = None   # '>' (CW) | '<' (CCW) | None

# SLOT_REGEX: r"\{(\d+)([><])?\}"   # group(1)=slot, group(2)=winding-or-None

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
    winding: Optional[str] = None   # NEW, defaulted for back-compat
```

### 5.5 Dependencies
No new dependencies. Standard library `typing.NamedTuple`; existing RDKit /
molassembler untouched.

## 6. Negative Constraints (The "Do NOTs")
- **DO NOT** change any 3D coordinate, RMSD, or placement outcome. Phase 1 is
  observationally inert on generation output.
- **DO NOT** un-mark `test_haptic_face_winding`; it must stay `expectedFailure`
  (it is Phase 3's acceptance gate, not Phase 1's).
- **DO NOT** normalize winding to a bool/enum. Store the literal `'>'`/`'<'`/`None`.
- **DO NOT** propagate winding through the legacy V2.4 sidecar `w:` path
  (`oin_parser.py:524+`, `oin/parser.py` w-tag reconstruction). Leave `None`.
- **DO NOT** modify `_build_connected_smiles` / `_template_generate` /
  `_stitch_eta_fragment` logic. Winding is *available* to them, not *consumed*.
- **DO NOT** touch `METAL_REGEX` or the `@desc` back-compat tolerance (D-3).

## 7. Risks & Mitigation
- **RISK-1 — a `vector_data` consumer is missed and silently misbehaves.**
  → **Mitigation:** the `NamedTuple` keeps `[0..2]` positional reads correct, so
  a missed positional reader degrades gracefully rather than shifting fields.
  §5.1 enumerates the complete set from a repo-wide grep; new tests exercise
  both the `inline.py` and `oin_parser.py` legs.
- **RISK-2 — C2 (`oin/parser.py:34`) is already broken and may be "load-bearing"
  in a way grep can't see.** → **Mitigation:** it's the incomplete `SMILESToXYZ`
  (TD-003) path and *cannot* run today (2-way unpack of 3-tuples raises).
  Fixing it to named iteration is strictly an improvement; do not expand its
  behavior (no winding, no w-tag semantics change) to avoid scope creep.
- **RISK-3 — accidental Phase-3 behavior leak** (e.g. wiring winding into
  placement to "make the test pass"). → **Mitigation:** §3.2 + §6 make the
  `expectedFailure` invariant an explicit acceptance criterion (US-004.2).
- **RISK-4 — `test_inline.py` equality asserts** now compare against a
  `NamedTuple`. → **Mitigation:** `SlotAssignment(1,0,0,None) == (1,0,0,None)`
  holds; update the literals to 4-length. (Note `== (1,0,0)` is False — length
  differs — so the update is mandatory, not optional.)

## 8. Success Metrics
- `uv run python -m unittest discover tests/unit` → **OK**, `expected failures=1`
  preserved (`test_haptic_face_winding`), test count up by the new plumbing tests.
- `uv run python -m unittest discover tests` (root) → **OK**.
- New unit assertions: `winding == '>'` / `'<'` / `None` captured by
  `parse_inline_string`, and the same value present on the corresponding
  `ParsedOIN.vectors[].winding` (US-001, US-002).
- Zero diff in generated XYZ for any existing generation test (behavioral
  inertness of Phase 1).

## 9. Execution Note (for /hyper-redteam → /hyper-resolve)
This Draft is a single-MiniPRD-sized effort. The Red Team should stress:
(a) completeness of the §5.1 consumer map, (b) whether `NamedTuple` equality/
pickling interacts with any test or the molassembler `ProcessPoolExecutor`
worker boundary, and (c) that no winding-driven behavior sneaks into placement.
