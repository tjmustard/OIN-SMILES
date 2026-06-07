
# Product Context

## Purpose
High-level "Why" of the OIN-SMILES project. Describes what the product does, who it serves, and the core value proposition.

## What It Is
OIN-SMILES is a Python library for **lossless conversion between 3D XYZ molecular structures and 1D SMILES strings** for Transition Metal Complexes (TMCs). It implements **Open Isomer Notation (OIN) v3.6** — an extended SMILES format that encodes coordination geometry, slot assignments, hapticity, winding direction, and P/N stereochemistry so that exact 3D reconstruction is possible from the string alone.

## The Problem It Solves
Standard SMILES is lossy for TMCs: coordination geometry (cis vs. trans), isomer identity, and P/N stereochemistry are destroyed when a 3D XYZ structure is flattened to 1D. OIN-SMILES fixes this by encoding the missing information in a compact inline string annotation, enabling exact round-trips.

## User Personas
- **Computational chemists** who need to store or exchange TMC structures without losing isomer identity
- **Cheminformatics / ML researchers** building TMC datasets where isomer labels must survive round-trips
- **Database curators** who need a canonical, human-readable string representation that is also machine-parseable

## Core Value Proposition
Given an XYZ file for a TMC, `XYZToSMILES().convert()` returns an OIN v3.6 string. Given that string, `OIN3DGenerator().generate()` reconstructs an XYZ block with the same coordination geometry, isomer, and P/N stereochemistry — verified by RMSD < 1.0 Å on all curated examples.

## Success Metrics
- Round-trip RMSD < 1.0 Å for all curated examples (currently 25+ complexes)
- OIN string identity preserved across XYZ → OIN → XYZ → OIN
- P/N stereocenter CIP codes survive round-trip
- `oin-smiles xyz2oin` / `oin-smiles oin2xyz` CLI works end-to-end

## Scope (as of v0.2.x)
- **In scope**: Square planar, octahedral, tetrahedral, TBP TMCs; P/N stereocenters; η-ligands (Cp, indenyl, ansa-metallocenes); CLI; Python API
- **Out of scope**: Arbitrary organics, proteins, periodic systems, GPU acceleration
