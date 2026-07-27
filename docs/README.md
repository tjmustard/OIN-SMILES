# OIN-SMILES Documentation

Project documentation for OIN-SMILES — lossless conversion between 3D
transition-metal-complex structures and 1D Open Isomer Notation (OIN) SMILES.

Start with the top-level [`README.md`](../README.md) for installation, usage, and
the OIN format overview. This folder holds deeper technical references:

| Document | What it covers |
| :--- | :--- |
| [OPTIMIZERS.md](OPTIMIZERS.md) | The 3D-generation optimizers: FF, the default g-xTB path, and the opt-in MACE MLIP — how to select, install, and configure them. |
| [GENERATION_PIPELINE.md](GENERATION_PIPELINE.md) | The full default OIN → 3D pipeline, stage by stage. |
| [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) | Failure modes that are understood, reproducible, and not bugs in the layer where they surface. |

See also [`CHANGELOG.md`](../CHANGELOG.md) for the release history and
[`CONTRIBUTING.md`](../CONTRIBUTING.md) for development setup.

## Session notes vs. documentation

Everything above describes **how the shipped software behaves** — it is written for
someone using or contributing to OIN-SMILES.

Everything an agentic coding session *produced while getting there* — measurement
reports, per-lane write-ups, A/B results, refuted hypotheses, status snapshots — lives
under [`agentic-notes/`](agentic-notes/README.md), organised by release. Those files are
the working record, not the manual: they are dated, they contradict each other across
releases, and several of them exist specifically to document an approach that **did not
work**.

**Adding a file here?** Read [`agentic-notes/README.md`](agentic-notes/README.md) first.
The `docs/` root is closed to new files; a `pre-commit` guard enforces it.
