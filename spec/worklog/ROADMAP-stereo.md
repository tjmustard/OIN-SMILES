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

### Phase 0 — Diagnostics (TASK-10) — DONE 2026-07-03
Three round-trip tests in `tests/unit/test_stereo_roundtrip_diagnostics.py`.
Results, which **revise the plan below**:
- **Ligand `@/@@` on the primary (template) path is NOT lost** — chiral-P
  (BDPP) and chiral-N (BDNN) round-trips passed byte-for-byte, so those two
  tests are now plain passing tests, not `expectedFailure`.
- **BUT the BDPP/BDNN fixtures don't test what their names imply.** Their
  stereocenters are the **backbone carbons** (`C[C@@H](C[C@H](C)...`); the P
  and N atoms each carry two identical phenyl groups and are NOT CIP centers.
  So "chiral phosphine/amine generation" is still **unverified** — we proved
  carbon `@/@@` survives, not P/N. **A dedicated fixture where P (or N) is
  itself the stereocenter (three distinct substituents) is required** before
  Phases 2/4 can be validated. This is the single most important Phase-0
  finding.
- **Winding IS lost (haptic test failed as expected).** Flipping ferrocene
  ring-0 `{0>}`→`{0<}` produced byte-identical 3D output either way:
  generation ignores input winding and re-derives it from geometry, because
  `SLOT_REGEX` (`oin/inline.py:44`) drops the `>`/`<` suffix. This is the real
  live gap; that test stays `expectedFailure` until Phase 3.

### Phase 1 — Preserve the signal (plumbing, MiniPRD)
- `SLOT_REGEX` → `\{(\d+)([><])?\}` (capture winding).
- `parse_inline_string` vector_data tuples gain a winding field (default
  `None` — back-compat for all existing callers).
- Thread winding through `generation/oin_parser.py::OINParser.parse` →
  `ParsedOIN` so the adapter can see it.
- No behavior change in placement yet; acceptance = existing suite green +
  new unit tests for the parse plumbing.

### Phase 2 — Verify/enforce ligand `@/@@` through ETKDG (MiniPRD) — SUPERSEDED 2026-07-03, folds into Phase 4
- **TASK-20 finding (2026-07-03):** the Phase-2 prerequisite fixture (a P/N
  atom with three distinct substituents *and* the metal bond — see below) is
  DEFINITIONALLY a Zone-A atom (bonded directly to the metal). Built and
  tested against `Rh-RR-DIPAMP-Cl2.xyz` ((R,R)-DIPAMP RhCl2, both P atoms
  genuine CIP stereocentres): `XYZToSMILES().convert()` produced
  `[Rh_SPL].Cc1ccccc1P{0}(CCP{1}(c1ccccc1)c1ccccc1C)c1ccccc1.[Cl]{2}.[Cl]{3}`
  — **no `@`/`@@` on either P atom.** Root cause:
  `ChiralityRecoveryUtility.recover()` (`core/chirality.py:155-158`)
  unconditionally clears the chiral tag on any P/N atom with
  `total_degree < 4` in the post-fragmentation ligand mol, and a
  metal-binding P/N atom always has exactly 3 fragment-local neighbours — so
  it is *always* Zone A, for *any* fixture. `PseudoAtomStrategy`
  (`chirality.py:22`), which exists to backfill the missing 4th substituent
  for this exact case, is never invoked anywhere (dead code).
- **Consequence: Phase 2's ETKDG verify/enforce experiment is inapplicable
  as scoped.** There is no `@`/`@@` on a metal-bound P/N atom in ANY input
  OIN for the generation side to preserve or lose — the signal never reaches
  generation because it is stripped one stage earlier, at encoding. Phase 2
  and Phase 4 are therefore **the same decision point**: Phase 4 must decide
  how (or whether) Zone-A stereocenters get encoded at all before there is
  anything for a Phase-2-style ETKDG experiment to verify.
- Original Phase-2 plan (kept for reference, now superseded): generate from a
  fixture's OIN and its `@↔@@`-flipped twin via `_template_generate`; assign
  CIP from the resulting 3D (`AssignStereochemistryFrom3D`) and check the
  codes are opposite; if tags are dropped, fix at the fragment-SMILES→mol
  boundary or post-embed check-and-reflect. This remains valid **once Phase 4
  produces an encoding that actually carries a P/N tag past the fragment
  boundary** — at that point re-run TASK-20-style tests as the acceptance
  check for the ETKDG side.

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

- `tests/fixtures/PdCl2-RR-BDPP.xyz`, `PdCl2-RR-BDNN.xyz` — **misleadingly
  named**: the "RR" stereocenters are the backbone CARBONS, not the P/N atoms
  (which carry two identical phenyls → not CIP centers). Useful as
  carbon-chirality pass-through fixtures only. **Do not treat as P/N-center
  coverage** (Phase-0 finding, 2026-07-03).
- `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz` (also under `tests/integration/`) —
  **fixture built 2026-07-03 (TASK-20):** RhCl2(DIPAMP), P atoms are genuine
  CIP stereocentres (metal + 3 distinct groups), built independently in
  Avogadro. Confirms the P/N-center gap is real and encoding-side (see Phase
  2/4). Golden: `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt`.
- `PdCl2-R-BINAP.xyz` (axial; end-state marker is the skipped
  `tests/unit/test_axial_chiral.py` test), `ferrocene.xyz`, `fac/mer_irppy3.xyz`.
- Golden strings: `tests/candidate_outputs/*.txt|.smi` (v3.7 style after
  TASK-04).
- RDKit CIP assignment from 3D is the stereo oracle (per ChiralPNStereocenters
  resolution H-1).
