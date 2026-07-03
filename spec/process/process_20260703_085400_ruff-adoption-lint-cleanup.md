# Process Document: Ruff Adoption & Lint Cleanup

**Generated:** 2026-07-03T08:54:00-07:00
**Session Focus:** Add Ruff configuration to the project and bring the codebase to zero `ruff check` errors under Google docstring conventions, without changing runtime behavior.

## Problem Statement

The OIN-SMILES codebase had no linter or formatter configured. The goal was to adopt Ruff (lint + format) with a specific rule set — pycodestyle (E/W), pyflakes (F), isort (I), and pydocstyle (D) with the Google docstring convention — and then actually make the code comply: a clean `ruff check` and consistent `ruff format` output, while keeping all tests green.

## Starting State

- Git HEAD: `7edae02` ("Add Phase 0 stereo diagnostics; correct roadmap fixture assumptions").
- Working tree already had uncommitted, in-progress work across several files (notably `tests/unit/test_inline.py` and `tests/unit/test_oin_generation.py`, plus `src/oinsmiles/generation/oin_parser.py`, `src/oinsmiles/oin/inline.py`, `src/oinsmiles/oin/parser.py`).
- `pyproject.toml` had no `[tool.ruff]` section.
- Test baseline: **55 unit tests passing** (`uv run python -m unittest discover tests`). The suite emits expected SMILES-parse-error log lines from negative-path tests; these are not failures.
- A first Ruff run (after config was added) reported **1095 errors** across `src/` + `tests/` (639 in `src/`, 456 in `tests/`), ~708 auto-fixable.

## Approach & Methodology

Spec-light, iterative, and safety-first. The work was planned in Plan Mode and executed in three phases with a test run as a gate after each phase: (1) a mechanical auto-fix pass, (2) hand-review of behavior-sensitive issues, (3) docstring authoring and line-length wrapping. The guiding principle was to separate the high-volume, zero-risk mechanical churn (whitespace, quotes, import sorting) from the small set of edits that could change behavior, and to prefer scoped config (`per-file-ignores`) over rewriting vendored or diagnostic code that legitimately deviates from library style.

## Steps Taken

1. **Added the Ruff config** to `pyproject.toml` per the user's spec (line-length 100; `select = ["E","F","W","I","D"]`; `ignore = ["D100","D104"]`; Google pydocstyle; `quote-style = "double"`). The one substitution was `known-first-party = ["oinsmiles"]`, resolved from `[tool.hatch.build.targets.wheel]` (`packages = ["src/oinsmiles"]`). Verified the config parsed with `uvx ruff check src/ --statistics`.

2. **Established a test baseline** with `uv run python -m unittest discover tests` → 55 pass. Re-run after every subsequent phase to catch regressions immediately.

3. **Scoped-out the noise (Phase 1 config).** Added `[tool.ruff.lint.per-file-ignores]`: `"tests/**" = ["D", "E402", "E501"]` (test/diagnostic scripts don't need library docstrings and legitimately do `sys.path` setup before imports / carry long data lines), and `"src/oinsmiles/utils/xyz2mol_local.py" = ["D"]` after confirming from its header that it is vendored (Jensen et al. xyz2mol, "Modified by Maria Harris Rasmussen 2024").

4. **Ran the auto-fixers (Phase 1).** `uvx ruff check src/ tests/ --fix` applied 672 safe fixes (import sorting, safe docstring-format rules); `uvx ruff format src/ tests/` reformatted 48 files. Re-ran tests → still 55 green. Remaining: 162 errors.

5. **Hand-fixed behavior-sensitive issues (Phase 2)**, each reviewed individually rather than auto-deleted:
   - Located the `F811`/`E402` sources with targeted `ruff --select` runs before editing.
   - `core/translator.py`: removed the dead duplicate `convert()` stub (the known TD-001 double-definition where the second shadowed the first).
   - `utils/xyz2mol.py`: consolidated a duplicate mid-file `from .oin_aligner import OINSanitizer` into the existing top-level import — after grepping `oin_aligner.py` to confirm no circular import back to `xyz2mol`.
   - `utils/xyz2mol.py`: renamed the ambiguous variable `I` (E741) to `inertia` in `_align_to_pai`, verifying all `\bI\b` references in the function were the inertia tensor before a scoped replace.
   - Converted 4 bare `except:` → `except Exception:` (`oin_aligner.py`, `tests/integration/rmsd_utils.py`).
   - Removed genuinely-dead assignments (`me_atom_start_idx`; a leftover `expected` in `test_writer.py`; an unused `fw`/`_patch_fw` binding), renamed a retained-for-intent computation to `_zone_a_indices`, and marked two intentional cases with `# noqa` (a `has_multi_eta` state flag; an RDKit-availability import probe in `test_chiral_p.py`).
   - After this, `tests/` reached **zero** Ruff errors; `src/` had 114 remaining (all docstrings + line length).

6. **Authored docstrings and wrapped long lines (Phase 3)** for the project's own modules only (vendored file already scope-ignored):
   - Wrote ~40 Google-style docstrings for undocumented public classes / methods / functions / `__init__`s (`core/graph.py`, `core/translator.py`, `oin/writer.py`, `oin/parser.py`, `oin/inline.py`, `generation/engine.py`, `generation/oin_parser.py`, `generation/molassembler_adapter.py`, `utils/oin_aligner.py`, `utils/xyz2mol.py`, `cli.py`).
   - Reformatted multi-line docstring summaries to satisfy D205/D415/D200/D212.
   - Wrapped ~36 over-length comments and strings that `ruff format` will not reflow: split the 220-char `MetalNon_Hg` SMARTS constant via implicit string concatenation (byte-identical), split debug `f"[DEBUG] …"` log strings the same way, and broke long comments across lines.

7. **Final verification.** `uvx ruff format src/ tests/` (1 residual file reformatted), `uvx ruff check src/ tests/` → **All checks passed!**, `uvx ruff format --check` → clean, `uv run python -m unittest discover tests` → **55 pass**, plus an import + `OINParser.parse()` smoke test across every heavily-edited module.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| `known-first-party = ["oinsmiles"]` | Leave placeholder | It is the real import name (source at `src/oinsmiles`) |
| Relax `D`/`E402`/`E501` for `tests/**` via `per-file-ignores` | Rewrite every test/diagnostic script | Diagnostic scripts legitimately use `sys.path` setup and long data lines; docstrings on test methods add no value |
| Scope-ignore `D` for `xyz2mol_local.py` | Rewrite its docstrings | It is vendored third-party code; preserving upstream style is correct |
| Split mechanical vs. behavior-sensitive into phases with test gates | One big `--fix --unsafe-fixes` run | Keeps the risky edits small, reviewed, and easy to isolate if a regression appears |
| Consolidate duplicate import to top of file | Add `# noqa: E402` in place | Confirmed no circular import, so the real fix is cleaner than suppression |
| `# noqa` for `has_multi_eta` and the RDKit import probe | Delete them | Both are intentional (a documented state flag; an availability check) |
| Wrap `MetalNon_Hg` SMARTS via implicit string concat | `# noqa: E501` | Keeps the constant under the length limit without breaking the SMARTS semantics |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| Ruff config | `pyproject.toml` | updated (`[tool.ruff]`, `[tool.ruff.lint]`, `per-file-ignores`, pydocstyle, isort, format) |
| Metal-center graph model | `src/oinsmiles/core/graph.py` | docstrings added |
| Translator | `src/oinsmiles/core/translator.py` | removed dead `convert()` stub; docstrings; comment wrap |
| XYZ→OIN pipeline | `src/oinsmiles/utils/xyz2mol.py` | import consolidation; `I`→`inertia`; docstrings; SMARTS/comment wraps |
| OIN aligner | `src/oinsmiles/utils/oin_aligner.py` | bare-except fix; docstrings; comment wraps |
| Molassembler adapter | `src/oinsmiles/generation/molassembler_adapter.py` | dead-var removal; `# noqa` flag; docstrings; debug-string wraps |
| Generation OIN parser | `src/oinsmiles/generation/oin_parser.py` | dataclass/class/method docstrings; docstring wraps |
| 3D generator | `src/oinsmiles/generation/engine.py` | `__init__` docstring; comment wrap |
| Inline / parser / writer | `src/oinsmiles/oin/{inline,parser,writer}.py` | docstrings; comment wraps |
| CLI | `src/oinsmiles/cli.py` | `main()` docstring; description-string wrap |
| Test suite (broad) | `tests/**` | `ruff format` reformat; 5 hand-fixes (F401/F841/E722) |

_No git commit was made this session — all changes are in the working tree._

## Results & Outcomes

- `uvx ruff check src/ tests/` → **All checks passed!** (down from 1095 errors).
- `uvx ruff format --check src/ tests/` → **53 files already formatted** (clean).
- `uv run python -m unittest discover tests` → **55 tests pass**, unchanged from baseline.
- Import + `OINParser.parse('[Pt_SPL].N{0}.[Cl]{1}')` smoke test succeeds across all edited modules, confirming the import-consolidation and inertia-rename edits did not break the runtime path.

## How to Reproduce

Prerequisites: `uv` installed, repo on a branch at (or based on) commit `7edae02`, no other uncommitted work you care about mixing in.

1. Add the Ruff sections to `pyproject.toml` (line-length 100; `select = ["E","F","W","I","D"]`; `ignore = ["D100","D104"]`; `known-first-party = ["oinsmiles"]`; Google pydocstyle; `quote-style = "double"`), plus `per-file-ignores`: `"tests/**" = ["D","E402","E501"]` and `"src/oinsmiles/utils/xyz2mol_local.py" = ["D"]`.
2. Record the baseline: `uv run python -m unittest discover tests` (expect 55 OK) and `uvx ruff check src/ tests/ --statistics`.
3. Run `uvx ruff check src/ tests/ --fix` then `uvx ruff format src/ tests/`. Re-run the tests.
4. Fix the remaining behavior-sensitive items by hand (use `uvx ruff check <path> --select F811,E402,F841,E722,E741 --output-format=concise` to locate them). Re-run the tests after each cluster.
5. Author docstrings and wrap remaining `E501` lines for the `src/` project modules (skip the vendored file). Use `uvx ruff check src/ --output-format=concise` to drive the list to zero.
6. Verify: `uvx ruff check src/ tests/` → passes; `uvx ruff format --check src/ tests/` → clean; `uv run python -m unittest discover tests` → 55 pass.

Gotchas / order-dependencies:
- Run the fixers and format **before** hand-editing; line numbers shift after reformatting, so re-query Ruff for current locations rather than trusting an earlier report.
- `--select <RULE>` on the CLI can re-surface config-ignored rules (e.g. `D100`); trust a plain `ruff check` for the real error set.
- When splitting a string across lines, use adjacent string literals (implicit concatenation) and verify the pieces concatenate byte-for-byte — critical for the SMARTS constant.

## Patterns & Lessons

- **Phase the cleanup and gate on tests.** Separating mechanical auto-fixes from the ~30 behavior-sensitive edits made the risky part small and auditable, and the per-phase test runs pinpoint any regression immediately.
- **Prefer scoped config over rewriting code that legitimately differs.** `per-file-ignores` for vendored (`xyz2mol_local.py`) and diagnostic (`tests/**`) code is the honest fix; it avoids churning external code and low-value test docstrings.
- **A `# noqa` is the right tool for intentional "violations"** (state flags, import-availability probes) — but only after confirming the code really is intentional.
- **Recommendation for committing:** land the mechanical `ruff format` reformat as its own dedicated commit (separate from feature work) to keep `git blame` clean. Note that `tests/unit/test_inline.py` and `tests/unit/test_oin_generation.py` had pre-existing uncommitted edits that the formatter has now intermingled — review those two before committing.
