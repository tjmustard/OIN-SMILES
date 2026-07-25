# Test Suite

This directory holds the OIN-SMILES test suite.

| Path | Contents |
|---|---|
| `unit/` | Fast, isolated unit tests (`python -m unittest discover tests/unit`). |
| `integration/` | End-to-end round-trip and pipeline tests. |
| `fixtures/` | Static `.xyz` inputs and conformer sets used by tests. |
| `candidate_outputs/` | Human-review artifacts emitted by tests (e.g. axial-chirality encodings). |
| `run_unit.sh`, `run_verification*.sh`, `run_integration.sh` | Convenience runners. |

## Running

```bash
bash tests/run_unit.sh                       # all unit tests (via uv)
uv run python -m unittest discover tests/unit # equivalent
# In a worktree without `uv sync`, use the main checkout's interpreter:
PYTHONPATH=$PWD/src .venv/bin/python -m unittest discover tests/unit
```

## Skipped tests

The suite intentionally skips a small number of tests. Skips are **expected**, not
failures. There are two kinds: one permanent skip, and several that are gated on the
environment (a fixture file, an optional dependency, or an OS capability) and only fire
when that prerequisite is absent.

A typical CI box **without** `mace-torch` and `xtb` reports **3 skips** (`~551 OK / 3
skip`): the permanent axial-chirality skip plus the two optional-tool skips below.

### Permanent skip (always skipped)

| Test | Guard | Why |
|---|---|---|
| `test_axial_chiral.py::TestAxialChiral::test_axial_chiral_descriptor_present` | `@unittest.skip(...)` | Axial chirality is not yet encoded by the pipeline. **Un-skip only after** the pipeline is updated **and** a human reviews `tests/candidate_outputs/axial_chiral_encoded.smi`. This is a known-limitation guard, not a broken test. |

### Optional-dependency skips (skip when the tool is missing)

| Test | Guard | Skips when |
|---|---|---|
| `test_generator3d_units.py::TestASEOptimizer::test_mace_constructs_or_skips` | `self.skipTest(...)` | `mace-torch` (optional extra) is not importable. |
| `test_generator3d_units.py::TestASEOptimizer::test_xtb_optimize_returns_energy_without_calculator` | `self.skipTest(...)` | the `xtb` binary is not on `PATH`. |

Install the optional extras / `xtb` to exercise these.

### Fixture- and capability-gated skips (run when the prerequisite is present)

These carry `@unittest.skipUnless(...)` guards and run normally as long as their
prerequisite is satisfied. They are listed so an unexpected skip can be traced to a
missing fixture or an absent capability.

| Test(s) | Guard | Runs when |
|---|---|---|
| `test_axial_chiral.py`, `test_chiral_p.py`, `test_chiral_n.py` (rdkit-gated cases) | `skipUnless(_RDKIT_AVAILABLE, ...)` | `rdkit` is importable. |
| `test_encoder_robustness.py::test_get_tmc_mol_raises_typed_error`, `::test_convert_propagates_typed_error` | `skipUnless(os.path.exists(FIXTURE), ...)` | `tests/fixtures/RAWJEG_comp_0.xyz` exists. |
| `test_encoder_robustness.py::test_benvog_recovers_via_cpu_budget_fallback` | `skipUnless(fixture exists and hasattr(os, "fork"), ...)` | `tests/fixtures/BENVOG_comp_0.xyz` exists **and** `os.fork` is available (POSIX). |
| `test_xyz2mol_errors.py::test_get_tmc_mol_raises_valueerror_on_unbuildable_ligand` | `skipUnless(_FIXTURE.exists(), ...)` | `tests/fixtures/ticat3_generated_broken.xyz` exists. |

Since these fixtures are checked into `tests/fixtures/`, these tests should run in any
normal checkout; a skip here signals a missing fixture (or a non-POSIX host for the
`fork`-gated case).
