# Product Context
## Purpose
To define the "Why" and "Who" of OIN-SMILES.

## Project Goals
- **Universal Format**: Create a standard 1D representation for 3D organometallic structures.
- **Reversibility**: Guarantee 1:1 mapping between `XYZ` and `OIN` string.
- **Automation**: Enable high-throughput pipeline deduplication and generative workflows.

## User Personas
- **Pipeline Engineer**: Needs fast, deterministic canonicalization for millions of structures.
- **Generative Chemist**: Needs an intuitive string format to convey complex 3D constraints (hapticity, chelates) to generation engines (Architector).

## Success Metrics
- **Round-Trip Fidelity**: 100% reconstruction of stereochemistry.
- **Uniqueness**: 0% collision for distinct isomers.
- **Drift Prevention**: No "SMILES Drift" (e.g., `[cH]` -> `c`).
