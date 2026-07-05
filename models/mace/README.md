# MACE Models Directory

This directory is intended for storing pre-trained MACE model weights (e.g., MACE-OMol25) required by the `OIN3DGenerator(optimizer="mace-omol25")` backend. 

All `.model`, `.pt`, and `.pth` files placed in this directory are intentionally ignored by git to prevent committing large binaries to the repository.

## How to download MACE-OMol25

The OMol25 checkpoints are hosted on Hugging Face under the FAIR Chemistry License. They cannot be downloaded anonymously via standard scripts.

1. **Request Access**: Go to the [facebook/OMol25 Hugging Face repository](https://huggingface.co/facebook/OMol25) and request access.
2. **Download Model**: Once approved, download the specific `.model` checkpoint you require (e.g., `MACE-OMol25.model` or `MACE-omol-0-extra-large-1024.model`).
3. **Place in this directory**: Move the downloaded file into this folder (`models/mace/`).
4. **Set Environment Variable**: Export the path to the model so the generator can find it:
   
   For standard OMol25:
   ```bash
   export MACE_OMOL25_MODEL_PATH="$(pwd)/models/mace/MACE-OMol25.model"
   ```
   
   For OMol-0 Extra Large:
   ```bash
   export MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH="$(pwd)/models/mace/MACE-omol-0-extra-large-1024.model"
   ```

You are now ready to run the OIN-SMILES generation tests with the MACE backend!
To use the extra large model, simply pass `optimizer="mace-omol-0-extra-large-1024"` to the `OIN3DGenerator`.
