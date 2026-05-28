# Active Context: OIN-SMILES Project

## Session State
- Last updated: 2026-05-27
- Previous clear: 2026-05-05T19:15:00Z (via /hyper-clear)

## Purpose
Captures current work state, recent completions, and next steps for the OIN-SMILES molecular chemistry library.

## Project Status: v0.2.0 Released (March 7, 2026) — v0.2.1 Fixes Complete

Core v0.2.0 feature set:
- ✅ XYZ ↔ OIN lossless round-trip conversion (OIN v3.6)
- ✅ SCINE Molassembler integration (template-based 3D generation + DG fallback)
- ✅ P/N stereocenter encoding (CIPAssigner, ChiralityRecoveryUtility, PseudoAtomStrategy)
- ✅ CLI (`oin-smiles xyz2oin`, `oin-smiles oin2xyz`)
- ✅ Return type change: `OIN3DGenerator.generate()` → `GeneratedStructure(xyz, mol)`

Post-release v0.2.1 bug fixes (all committed):
- ✅ **2026-05-05**: CLI `oin2xyz` fixed to access `.xyz` from `GeneratedStructure` (commit `74bab57`)
- ✅ **2026-05-05**: Disable `_stitch_multi_eta_fragment`; document TiCat1/3/4 limitations
- ✅ **2026-05-05**: Fix Ir(ppy)3 regression; improve ETKDG fallback sanitization
- ✅ **2026-04-24**: Add multi-eta ligand support for ansa-metallocenes
- ✅ **2026-04-03**: Fix missing bonds in eta-ligand generated mol files (ferrocene, TiCp2Me2)
- ✅ **2026-03-28**: TiCat1/3/4 aromatic η-ligand support via de-aromatization + ETKDG strategy

## Current Focus: v0.2.2 Direct Parser (In Progress)

The Direct Parser replaces the legacy `OINParser.parse()` + `MolassemblerAdapter.generate()` pipeline with a single integrated path. Integration is DEFERRED due to 5 blockers identified in audit (see `spec/audit/DirectParser_IntegrationAudit_20260506.md`).

### MiniPRD Audit Status
**Completed (audited & archived):**
- ✅ MiniPRD_DirectParser_RegexPreprocessor_AUDITED.md (2026-05-06) — `_extract_oin_constraints()` 3-tuple return, 14 tests passing
- ✅ MiniPRD_DirectParser_ASTTokenization_AUDITED.md (2026-05-06) — `tokenize_unsanitized_smiles()`, 16 tests passing
- ✅ MiniPRD_DirectParser_MolassemblerInstantiation_AUDITED.md (2026-05-06) — `construct_molassembler_mol()`, 20 tests passing
- ✅ MiniPRD_DirectParser_FragmentMapping_v0.2.2_AUDITED.md (2026-05-10) — `_extract_oin_constraints()` fragment mapping, 55 tests passing

**Pending execution (next up):**
- ⏳ MiniPRD_DirectParser_Permutation_v0.2.2.md — Permutation selection (cis/trans, fac/mer)
- ⏳ MiniPRD_DirectParser_EtaBonds_v0.2.2.md — Eta-bond translation from OIN vertex indices
- ⏳ MiniPRD_DirectParser_Polydentate_v0.2.2.md — Polydentate ligand connectivity
- ⏳ MiniPRD_DirectParser_Tests_v0.2.2.md — Integration test suite
- ⏳ MiniPRD_DirectParser_Verification.md — RMSD and round-trip fidelity validation

## Known Limitations & Technical Debt
- TD-001: `XYZToSMILES.convert()` defined twice — second shadows first (refactor needed)
- TD-002: `SMILESToXYZ` incomplete in translator.py (dummy implementation)
- TD-003: `OINInlineHandler.generate_inline_string()` has `pass` stub
- TD-005: `TEMPLATES`/`TEMPLATE_SPECS` duplicated (consolidate to shared constants)
- **pyproject.toml version mismatch**: shows `0.1.0` but CHANGELOG has `[0.2.0] - 2026-03-07`

## Test Commands
```bash
uv sync                                                    # Install dependencies
uv run python -m unittest discover tests                  # All unit tests
uv run python tests/integration/verify_roundtrip.py      # Round-trip validation (XYZ→OIN→XYZ)
uv run python tests/integration/verify_roundtrip.py --include-tmqm  # Include tmQM dataset
```

## Audited & Archived MiniPRDs

**v0.2.0 Release (All 5 Complete, 2026-05-05):**
- ✅ MiniPRD_MolassemblerSpike_AUDITED.md — Molassembler import, picklability, ProcessPoolExecutor
- ✅ MiniPRD_MolassemblerAdapter_AUDITED.md — Template placement, DG fallback, timeout
- ✅ MiniPRD_ChiralEncoding_AUDITED.md — CIPAssigner, ChiralityRecoveryUtility, PseudoAtomStrategy
- ✅ MiniPRD_ChiralTests_AUDITED.md — Unit/integration tests for P/N stereocenter encoding
- ✅ MiniPRD_CLI_AUDITED.md — `oin-smiles` CLI with `xyz2oin` and `oin2xyz` (archived 2026-05-05)

## Files Recently Modified
- `README.md` — Redesigned with logo, badges, ToC, emoji section headers (2026-05-27)
- `CHANGELOG.md` — Added CLI fix entry and Direct Parser audit completions
- `media/OIN-SMILES-logo-dark.webp` — Added WebP logo for README
- `media/OIN-SMILES-logo-light.webp` — Added WebP logo (light variant)
- `src/oinsmiles/cli.py` — Fixed `_cmd_oin2xyz` to access `.xyz` attribute (commit `74bab57`)
- `src/oinsmiles/generation/molassembler_adapter.py` — ETKDG fix for TiCat1/3/4
- `src/oinsmiles/generation/oin_parser.py` — Direct Parser components (regex, AST, Molassembler, fragment mapping)
- `tests/integration/verify_roundtrip.py` — Enhanced diagnostics, MOL/SDF output, DG strategy benchmarking
