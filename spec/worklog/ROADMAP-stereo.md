# ROADMAP: Generation-side stereochemistry (OIN→XYZ)

Status: Phase 0 pending (TASK-10). Phases 1–4 each require a HACF MiniPRD
(`/hyper-architect` → `/hyper-redteam` → `/hyper-resolve`) before
implementation — this roadmap is direction, not spec.

## The problem

The XYZ→OIN direction faithfully encodes stereo (P/N `@/@@` via
`ChiralityRecoveryUtility`, eta winding `{n>}`/`{n<}`, isomer via slot
ordering). The OIN→XYZ direction silently drops most of it:

| Signal | Where lost | Consequence |
|---|---|---|
| Winding `{n>}`/`{n<}` | `SLOT_REGEX` (`src/oinsmiles/oin/inline.py:44`) matches but doesn't capture | haptic face / ring orientation not constrainable |
| P/N CIP enforcement | nothing on the generation path reads CIP or verifies `@/@@` post-embed | chiral phosphines/amines may invert |
| Haptic face | `_stitch_eta_fragment` only flips ring normal toward metal (`molassembler_adapter.py:627`) | prochiral face of substituted Cp/arene arbitrary |
| Metal isomer | OK — carried by slot ordering into template placement / permutation search | (not a gap) |

Architecture reality check (established 2026-07-02): molassembler is only the
FALLBACK builder. The primary path is RDKit ETKDG per-fragment + Kabsch
alignment onto geometry templates (`_template_generate`,
`molassembler_adapter.py:908`). Ligand `@/@@` tags DO ride along in fragment
SMILES (isomericSmiles=True) but nothing verifies they survive
sanitize/AddHs/embed.

**Standing decision (NOTES.md D-4):** do NOT evaluate alternative structure
builders until Phases 1–3 produce data. Neither current builder has ever been
GIVEN the stereo signals; judging them (or shopping for replacements) before
fixing our own information flow would be premature.

## Phases

### Phase 0 — Diagnostics (TASK-10, lightweight task, no MiniPRD)
Three `expectedFailure` round-trip tests (chiral-P BDPP, chiral-N BDNN,
ferrocene winding). Output: concrete mismatch data per pathway. Each later
phase flips its test to a hard assert.

### Phase 1 — Preserve the signal (plumbing, MiniPRD)
- `SLOT_REGEX` → `\{(\d+)([><])?\}` (capture winding).
- `parse_inline_string` vector_data tuples gain a winding field (default
  `None` — back-compat for all existing callers).
- Thread winding through `generation/oin_parser.py::OINParser.parse` →
  `ParsedOIN` so the adapter can see it.
- No behavior change in placement yet; acceptance = existing suite green +
  new unit tests for the parse plumbing.

### Phase 2 — Verify/enforce ligand `@/@@` through ETKDG (MiniPRD)
- Experiment first: generate from BDNN OIN and from its `@↔@@`-flipped twin
  via `_template_generate`; assign CIP from the resulting 3D
  (`AssignStereochemistryFrom3D`) and check the codes are opposite.
- If tags are dropped: fix at the fragment-SMILES→mol boundary (where
  `_stitch_fragment` builds its embed mol); if ETKDG ignores them: post-embed
  check-and-reflect of the offending stereocenter.
- Flips TASK-10 tests (a)/(b) for the template path.

### Phase 3 — Haptic face control (MiniPRD)
- In `_stitch_eta_fragment` (`molassembler_adapter.py:509`): after Kabsch
  placement, compute the signed circulation of the ring atoms (in SMILES
  order) about the metal→centroid axis; if it disagrees with the Phase-1
  winding marker, mirror the fragment across the ring plane before final
  placement. Same for `_stitch_multi_eta_fragment` (`:110`).
- Flips TASK-10 test (c).
- Open design question for the MiniPRD: does winding as defined in V3.6
  (`_determine_winding`, `utils/oin_aligner.py:589`) actually pin the
  prochiral face for substituted rings, or only ring direction? May need a
  substituted-Cp fixture (e.g. methylcyclopentadienyl) to make face identity
  observable.

### Phase 4 — Zone-A P/N stereo + builder decision (MiniPRD)
- Zone-A = P/N atoms bonded directly to the metal; XYZ→OIN currently CLEARS
  their tags (`core/chirality.py:154-157`, per ChiralPNStereocenters spec)
  because fragment-local CIP is ill-defined without the metal. Generation
  therefore has no signal for metal-bound stereocenters; decide encoding
  (extend OIN? derive from slot geometry?) and enforcement (post-embed
  check-and-reflect with metal present, or molassembler atom
  stereopermutators in `_molassembler_worker`, `molassembler_adapter.py:1418`).
- THEN, with Phases 1–3 data in hand: keep ETKDG+templates, extend
  molassembler usage, or adopt an alternative builder. This is the only point
  where "alternative structural building tool" gets decided.

## Fixtures & oracles

- `tests/fixtures/PdCl2-RR-BDPP.xyz` (P-chiral), `PdCl2-RR-BDNN.xyz`
  (N-chiral), `PdCl2-R-BINAP.xyz` (axial; end-state marker is the skipped
  `tests/unit/test_axial_chiral.py` test), `ferrocene.xyz`, `fac/mer_irppy3.xyz`.
- Golden strings: `tests/candidate_outputs/*.txt|.smi` (v3.7 style after
  TASK-04).
- RDKit CIP assignment from 3D is the stereo oracle (per ChiralPNStereocenters
  resolution H-1).
