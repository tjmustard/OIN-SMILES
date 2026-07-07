# MACE Models Directory

This directory holds the pre-trained MACE foundation-model weights used by the
**opt-in** MACE optimizer of the MetalloGen 3D backend
(`OIN3DGenerator(optimizer="mace-omol-0-extra-large-1024")` or `"mace-omol25"`).

All `.model`, `.pt`, and `.pth` files here are git-ignored so large binaries never
land in the repository.

> MACE is optional. The default engine uses the fast FF + g-xTB path and needs
> neither `torch` nor these weights. Install the ML stack only if you want MACE:
> `uv sync --extra mace`.

---

## Recommended: MACE-omol-0 (extra large) — automated download

This checkpoint is a **public** asset on the ACEsuit `mace-foundations` GitHub
release, so it can be fetched with no login. From the repository root:

```bash
bash tools/install_mace_weights.sh
```

The script downloads `MACE-omol-0-extra-large-1024.model` (~400 MB) into this
directory and sets `MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH` in a repo-root `.env`
file (loaded automatically at runtime via `python-dotenv`). It is idempotent —
re-running skips the download if the file is already present.

Manual equivalent:

```bash
curl -fL -o models/mace/MACE-omol-0-extra-large-1024.model \
  https://github.com/ACEsuit/mace-foundations/releases/download/mace_omol_0/MACE-omol-0-extra-large-1024.model
export MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH="$(pwd)/models/mace/MACE-omol-0-extra-large-1024.model"
```

Then:

```bash
uv sync --extra mace
uv run --extra mace oin-smiles oin2xyz "<OIN string>" --optimizer mace-omol-0-extra-large-1024
```

(`uv run` re-syncs to the default light environment each call, so pass `--extra
mace` to `uv run` too — or activate the venv and call `oin-smiles` directly.)

---

## Alternative: MACE-OMol25 (Hugging-Face-gated)

The OMol25 checkpoints are hosted on Hugging Face under the FAIR Chemistry
License and **cannot** be downloaded anonymously — they require an approved
access request.

1. **Request access** at the [facebook/OMol25](https://huggingface.co/facebook/OMol25) repository.
2. **Download** the `.model` checkpoint once approved.
3. **Place** it in this directory (`models/mace/`).
4. **Register** the path:
   ```bash
   export MACE_OMOL25_MODEL_PATH="$(pwd)/models/mace/MACE-OMol25.model"
   ```
   then use `optimizer="mace-omol25"`.
