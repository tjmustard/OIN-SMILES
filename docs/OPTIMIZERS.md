# 3D-Generation Optimizers

When converting an OIN-SMILES string to a 3D structure, the default MetalloGen
engine (`OIN3DGenerator`) first builds a force-field-relaxed conformer pool, then
optionally refines the selected conformer with a higher-level optimizer. The
optimizer is chosen with the `optimizer` argument (Python) or `--optimizer` (CLI).

## Summary

| `optimizer` value | Backend | Extra install | Notes |
| :--- | :--- | :--- | :--- |
| `"xtb"` **(default)** | g-xTB via the `xtb` binary (subprocess) | `tools/install_gxtb.sh` | Semi-empirical refinement; **falls back to FF** with a warning if `xtb` is not on `PATH`. |
| `"ff"` / `"none"` / `None` | Force field only (UFF/MMFF) | none | Fastest; no `torch`, no external binary. |
| `"mace-omol-0-extra-large-1024"` | MACE MLIP (ASE + LBFGS) | `uv sync --extra mace` + weights | Most accurate; **fails loudly** if `mace-torch` or the weights are missing. |
| `"mace-omol25"` | MACE MLIP (OMol25 checkpoint) | `uv sync --extra mace` + HF-gated weights | As above, using the Hugging-Face OMol25 checkpoint. |

## Default: g-xTB

```bash
bash tools/install_gxtb.sh          # puts the `xtb` (g-xTB) binary on PATH
oin-smiles oin2xyz "<OIN string>"   # optimizer defaults to "xtb"
```

The g-xTB path shells out to `xtb <struct>.xyz --gxtb --opt`. If the binary is not
found, or the run fails, generation degrades gracefully to the FF geometry (a
warning is printed) — it never hard-fails on a missing optimizer.

## FF-only (no torch, no binary)

```python
from oinsmiles.generation.engine import OIN3DGenerator
gen = OIN3DGenerator(optimizer="ff")   # or optimizer=None
result = gen.generate("<OIN string>")
```

This is the lightest path and matches the default `uv sync` install (no `torch`,
no `mace-torch`).

## Opt-in: MACE (machine-learning interatomic potential)

MACE gives the most accurate geometries but pulls a pinned CUDA-11.8 `torch`.

```bash
uv sync --extra mace                 # installs mace-torch + torch
bash tools/install_mace_weights.sh   # downloads the extra-large weights + sets .env
uv run --extra mace oin-smiles oin2xyz "<OIN string>" --optimizer mace-omol-0-extra-large-1024
```

> **Note:** pass `--extra mace` to `uv run` as well, not just `uv sync`. `uv run`
> re-syncs the environment to the default (light) dependency set on each call, so a
> bare `uv run oin-smiles ... --optimizer mace-...` would report `mace-torch is not
> installed`. Alternatively, activate the venv (`source .venv/bin/activate`) after
> `uv sync --extra mace` and call `oin-smiles` directly.

### Weights & environment variables

MACE weights are **not** bundled (they are large binaries, git-ignored under
`models/mace/`). The optimizer locates them via environment variables, which the
installer writes into a repo-root `.env` (auto-loaded via `python-dotenv`):

| `optimizer` | Environment variable |
| :--- | :--- |
| `"mace-omol-0-extra-large-1024"` | `MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH` |
| `"mace-omol25"` | `MACE_OMOL25_MODEL_PATH` |

See [`models/mace/README.md`](../models/mace/README.md) for the download details
(the extra-large checkpoint is a public GitHub release; OMol25 is Hugging-Face-gated).
Unlike the g-xTB path, MACE **fails loudly** when its dependency or weights are
missing, so an accuracy request is never silently downgraded.
