
# Active Context

## Purpose
Captures the current state of OIN-SMILES development. Updated after significant task completions. Read first to understand where to pick up.

## Current State (as of 2026-06-07)

### Release Status
- **v0.2.0** — Released 2026-03-07: Molassembler backend, P/N stereocenter encoding, CLI, OIN v3.6
- **v0.2.1** — In `[Unreleased]`: CLI fix, P/N fixtures (BINAP, BDNN, BDPP), TiCat1/3/4 ETKDG fix, Direct Parser audit completions
- **v0.2.2** — Planned: Direct Parser bugfixes (5 P0/P1 blockers, see audit doc below)

### No Active Sprint
No MiniPRD is currently in progress. Next work requires creating `MiniPRD_DirectParser_Bugfixes_v0.2.2.md` via `/hyper-architect`.

## Direct Parser — Deferred to v0.2.2
**Audit doc**: `spec/audit/DirectParser_IntegrationAudit_20260506.md`

Integration is blocked by 5 issues in `src/oinsmiles/generation/oin_parser.py`:
1. **(P0 Blocker)** Fragment rank ↔ atom index mapping missing → "bond to self" errors
2. **(P1)** Polydentate ligand connectivity not handled
3. **(P0 Blocker)** No permutation/isomerism selection (cis vs. trans, fac vs. mer)
4. **(P1)** Eta bond translation to atom indices broken
5. **(P2)** Missing test coverage for direct parser

**Production pipeline (current)**: Uses legacy `OINParser.parse()` + `MolassemblerAdapter.generate()` — all integration tests pass.

## Recent Completions

### v0.2.1 Work (in [Unreleased])
- **Direct Parser MiniPRD audits (4 of 5)**: RegexPreprocessor, ASTTokenization, MolassemblerInstantiation, FragmentMapping — all audited, archived. Integration MiniPRD deferred.
- **CLI fix**: `oin2xyz` command updated to access `.xyz` from `GeneratedStructure` return type
- **P/N stereocenters fixtures**: PdCl2-R-BINAP, PdCl2-RR-BDNN, PdCl2-RR-BDPP — all pass round-trip RMSD < 1.0 Å
- **TiCat1/3/4 3D generation fix**: ETKDG + de-aromatization strategy for aromatic η5 ligands. See `docs/ETKDG_AROMATIC_FIX.md`.
- **`_extract_oin_constraints()` rename**: 31 call sites updated, `fragment_to_atom_mapping` added to return tuple

### Toolchain (HACF)
- HACF updated to v0.5.1 (post-v0.5.0 skills: `hyper-contextualize`, `hyper-handoff`, `hyper-grill-docs`)
- Agent instruction files (AGENTS.md, CLAUDE.md, GEMINI.md) framing-banner aligned

## Known Limitations (v0.2.1)
- TiCat1/3/4 round-trip bonding inference fails (SINGLE bonds instead of AROMATIC — de-aromatization trade-off)
- `SMILESToXYZ` in `translator.py` is incomplete (dummy atoms, not real SMILES parsing) — TD-003
- `XYZToSMILES.convert()` defined twice — second shadows first — TD-001

## Key Files for Next Session
- `src/oinsmiles/generation/oin_parser.py` — Direct Parser implementation (blockers above)
- `spec/audit/DirectParser_IntegrationAudit_20260506.md` — full analysis of 5 blockers
- `spec/compiled/SuperPRD.md` — system of record (v1.1.0)
- `CHANGELOG.md` `[Unreleased]` — pending v0.2.1 release entries
