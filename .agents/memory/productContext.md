---
name: Project Goals & Product Context
description: OIN-SMILES research goals, target users, and success metrics
type: project
---

# OIN-SMILES Project: Goals & Context

## The Problem We Solve

**Challenge:** Converting 3D transition metal complex (TMC) structures to 1D SMILES and back is lossy — stereochemistry, coordination geometry, and isomer information are destroyed. This breaks computational chemistry workflows that require:
- Exact isomer preservation for database curation
- Lossless round-trip conversions for structure analysis
- Canonical representation for machine learning models

**Solution:** OIN-SMILES implements **Open Isomer Notation (OIN)** — a 1D string encoding that captures all stereochemical and geometric information needed for perfect 3D reconstruction.

## Target Audience

1. **Computational Chemistry Researchers**
   - Need lossless structure representation for TMC datasets
   - Use case: Training ML models on metal complexes (bond predictions, geometry estimation, property regression)
   - Pain point: Current SMILES notation loses coordination geometry → models see "generic" metal atoms

2. **Chemical Database Curators**
   - Maintain curated TMC libraries (e.g., tmQM, TMQM with explicit H)
   - Need canonical, round-trip-safe representations
   - Pain point: Merging databases loses isomer distinctions; duplicate detection is unreliable

3. **Quantum Chemistry Software Developers**
   - Need to exchange 3D structures with external tools (RDKit, SCINE Molassembler, etc.)
   - Want high-fidelity 3D generation from text input
   - Pain point: SMILES → 3D pipelines struggle with metal coordination geometries

## Project Goals

### Goal 1: Lossless Round-Trip Conversion ✅ (v0.2.0)
- **Definition:** XYZ → OIN → XYZ with exact isomer preservation
- **Success Metric:** RMSD < 1.0 Å for all test complexes
- **Status:** ACHIEVED
  - Cisplatin: RMSD 0.2 Å ✅
  - Transplatin: RMSD 0.2 Å ✅
  - Ferrocene: RMSD < 1.0 Å ✅
  - Ir(ppy)3: RMSD ~5.0 Å (atom ordering issue, not geometry)

### Goal 2: Stereocenter Encoding for P/N Atoms ✅ (v0.2.0)
- **Definition:** Chiral phosphorus and nitrogen centers are preserved through round-trip
- **Implementation:** CIP codes embedded in OIN string (e.g., `[Pt@SP1_SPL]`)
- **Success Metric:** Matching @/@@ stereo markers in input vs. output SMILES
- **Status:** ACHIEVED with PseudoAtomStrategy fallback for edge cases

### Goal 3: Fast, Reliable 3D Generation from OIN ✅ (v0.2.0)
- **Definition:** OIN → 3D conformer in < 60 seconds (default timeout)
- **Implementation:** SCINE Molassembler (template-based placement + DG fallback)
- **Replaces:** Architector (slower, less reliable)
- **Success Metric:** All known geometries (SP1, SP3, OC-6, TBPY-5, etc.) generate on first attempt
- **Status:** ACHIEVED; TiCat1/3/4 (aromatic η-ligands) fixed via ETKDG (v0.2.1)

### Goal 4: Open-Source, Production-Ready Library ✅ (v0.2.0)
- **Definition:** Usable as both Python API and CLI
- **API:** `XYZToSMILES().convert(xyz_file)` → OIN string
- **API:** `OIN3DGenerator().generate(oin_string)` → `GeneratedStructure(xyz, mol)`
- **CLI:** `oin-smiles xyz2oin complex.xyz` / `oin-smiles oin2xyz "[Pt@SP1_SPL]..."`
- **Status:** ACHIEVED; v0.2.0 released March 7, 2026

### Goal 5 (Future): Support tmQM & Extended TMC Dataset ⏳
- **Definition:** Validated on > 100 diverse metal complexes (Pt, Pd, Fe, Ir, Ti, Mo, Ru, etc.)
- **Current Coverage:** 
  - ✅ Verified: Pt, Fe, Ir (fac/mer), Ti (ansa-metallocenes)
  - ⏳ In Testing: Pd, Ru, Mo, W
  - ⏳ Backlog: Rare earths, main-group complexes
- **Success Metric:** RMSD < 2.0 Å for 95% of test set; < 5% geometric failures

## Success Metrics

### Correctness
- **Round-trip RMSD:** < 1.0 Å for standard geometries, < 2.0 Å for edge cases
- **Stereo preservation:** 100% of chiral centers match input @/@@ markers
- **Bond preservation:** Generated mol files include all M–L dative bonds and full ligand connectivity

### Performance
- **3D generation time:** < 60 seconds (default) for any OIN string
- **XYZ→OIN conversion:** < 5 seconds for typical TMC
- **Memory usage:** < 500 MB for typical workflow (mol in memory, no disk buffering)

### Robustness
- **Test coverage:** > 90% line coverage; all known complexes in regression test suite
- **Error handling:** Graceful failure (PseudoAtomStrategy, DG fallback) instead of crashes
- **Edge cases:** Non-standard valence, mixed-valence, polynuclear complexes handled without loss

### Usability
- **CLI availability:** `oin-smiles` command registered as package entry point
- **Documentation:** README with usage examples; CHANGELOG with API changes
- **API stability:** Backward-compatible returns (GeneratedStructure preserves `.xyz` string access)

## Development Roadmap

### Released (v0.2.0 - March 7, 2026)
- ✅ Core XYZ→OIN→XYZ pipeline
- ✅ SCINE Molassembler backend (replaces Architector)
- ✅ P/N stereocenter encoding
- ✅ CLI
- ✅ Round-trip tests for Pt, Fe, Ir

### In Progress (v0.2.1 - May 2026)
- 🔧 Bug fixes: Ir(ppy)3 atom ordering, TiCat1/3/4 ETKDG finalization
- 🔧 Extended test coverage: Pd, Ru, Mo
- 🔧 Documentation: ETKDG strategy deep-dive

### Future (Post-v0.2.1)
- ⏳ tmQM integration: Validated generation for all 101 structures
- ⏳ CLI enhancements: Batch processing, Molblock output, visualization
- ⏳ Performance: Parallel conformer generation for ensemble methods
- ⏳ Integration: Plugins for RDKit, Open Babel, Avogadro

## Known Limitations

1. **Atom Ordering Issues:** Ir(ppy)3 generates with ~5 Å RMSD due to RDKit atom index mismatch (not a geometry problem)
2. **AROMATIC Bond Inference:** Generated mols from de-aromatized η-ligands have SINGLE instead of AROMATIC bonds (working as designed for safety)
3. **Polynuclear Complexes:** Not yet tested; may require template expansion
4. **Main-Group Complexes:** OIN designed for TMCs; no explicit support for pure organic or rare-earth systems

## Research Impact

**Expected Users:**
- ML researchers building metal complex models (e.g., bond prediction, property regression)
- Chemical database projects (tmQM, TMQM-H, open-source TMC libraries)
- Quantum chemistry software integrators (SCINE packages, electronic structure codes)
- Academic groups studying coordination chemistry

**Potential Citations:**
- TMC machine learning (e.g., Stieber et al., Jørgensen et al.)
- Chemical informatics for metals (e.g., Maldonado, O'Boyle)
- Chemoinformatics toolkits (RDKit, Open Babel)
