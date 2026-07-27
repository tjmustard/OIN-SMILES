# RedTeam_Report: OIN to Molassembler Direct Parser

**Generated:** 2026-05-05  
**Reviewer Role:** Red Team (Adversarial Architecture Review)  
**Status:** READY FOR TRIAGE via `/hyper-resolve`

---

## Executive Summary

The Draft PRD proposes replacing RDKit-based ETKDG 3D generation with a direct Molassembler parser. The core architecture is sound, but **7 critical ambiguities and 3 architectural blindspots** require resolution before implementation. Most severe: 
- **Determinism guarantee is unproven** (SCINE DG seeding not specified)
- **Timeout handling is TBD** (blocks deployment)
- **Regex parser has silent failure modes** (malformed OIN crashes in extraction)
- **Test coverage is fixture-limited** (TiCat1/3/4 only; unknown shape generalization)

---

## 1. Introduction & Goals Analysis

### Clarifying Questions

1. **Aromaticity vs. Haptic Bonding Root Cause:** The problem statement claims RDKit "cannot kekulize extracted Cp/indenyl rings" due to Hückel enforcement. But kekulization is a *sanitization step*, not a bonding step. Is the real issue that:
   - RDKit rejects Cp rings as non-aromatic after kekulization (breaking the graph)?
   - De-aromatization corrupts the SMILES representation?
   - Or the eta bonds themselves are incompatible with RDKit's bond model?
   
   **Need:** Precise root cause diagnosis. If it's purely SMILES → Mol conversion, why not use `sanitize=False` without de-aromatization?

2. **"Preserves chemical graph integrity" vs. Unsanitized AST:** The solution claims unsanitized AST "preserves" the graph, but doesn't explain *how*. RDKit's unsanitized parsing may produce:
   - Implicit hydrogens not counted
   - Valence errors not reported
   - Aromatic flags not set
   
   Will downstream SMILES generation (via `.GetSmiles()` on the mol after Molassembler construction) still produce valid OIN-SMILES?

3. **Audience Scope Creep:** "Internal use only" contradicts the fact that `OIN3DGenerator.generate()` is part of the public API (`SMILESToXYZ`). How does "internal" reconcile with customer calls to `SMILESToXYZ.convert(oin_smiles)`?

### What-If Scenarios

1. **Malformed OIN Injection:** A user typo in OIN-SMILES (e.g., `[Pd_SQP].[Cl]{99}` where vertex index 99 is out of bounds) hits the regex extractor. Current spec says regex "assumes valid format" — what happens?
   - **Scenario A (Silent Failure):** Regex extracts `{99}`, constraint dict stores it, Molassembler validation fails deep in shape injection with an opaque error
   - **Scenario B (Crash):** Index 99 escapes bounds check, Molassembler segfaults on vertex assignment
   - **Risk:** User gets non-actionable error; no path to debug "what OIN was invalid?"

2. **Floating-Point Non-Determinism:** SCINE's DG engine produces coordinates. Across different hardware (GPU vs. CPU), different BLAS libraries (MKL vs. OpenBLAS), or different DG iterations, are coordinates bit-identical?
   - **Scenario:** Test passes on CI (Intel Xeon), fails in user's environment (ARM Mac)
   - **RMSD fidelity:** If conformer changes by 0.2Å due to hardware difference, does RMSD still meet ≤ 1.0Å threshold? Tolerance is only 5x the noise.

3. **Molassembler Version Incompatibility:** Pinned to `< 5.0.0` but no upper bound on 3.x or 4.x API drift. 
   - **Scenario:** User installs `scine_molassembler==4.5.0` (hypothetically), API for `set_shape()` changes signature silently
   - **Failure Mode:** Code compiles, but shape assignment produces wrong geometry

4. **Atom Map Index Collision:** When replacing shape codes and vertex indices with RDKit atom maps (e.g., `[Pd:1]`), how are atom map numbers assigned?
   - **Scenario:** OIN-SMILES already contains atom maps: `[Pd:5].[Cl:6]` + new extraction adds `[Pd:1].[Cl:2]`
   - **Result:** Atom map dict has duplicates or off-by-one indices; RDKit parses incorrectly

5. **Partial Molecule Construction Failure:** Sequential `mol.push_back(atom)` succeeds for 90 atoms, then fails on atom 91 with memory error or invalid atomic number.
   - **Current state:** Mol object is half-constructed, inconsistent
   - **Error handling:** Who reverts? Molassembler transaction semantics not specified

### Points for Improvement

1. **Explicit Error Boundary for Invalid OIN:** Instead of "assumes valid format," implement a pre-validation regex that:
   - Ensures all brace-enclosed indices are integers in range [0, max_vertices_for_shape)
   - Rejects malformed patterns (unmatched braces, non-numeric indices, etc.)
   - Raises `ValueError` with human-readable message before Molassembler touch
   
   **Example:**
   ```python
   def validate_oin_format(oin_smiles: str, shape_vertices: dict) -> None:
       """Fail fast on invalid OIN syntax; precondition error not runtime error."""
       # 1. Extract shape codes
       # 2. Look up max_vertices per shape
       # 3. Validate all eta indices against max
       # Raise on first violation with shape + max + found_index
   ```

2. **Document Determinism Guarantee:** Add explicit section to spec:
   - Is DG seeded? If yes, document seed value and how it's preserved across versions
   - Are coordinates deterministic across hardware/BLAS? How were they verified?
   - Test protocol: "Run on N different machines, compute coordinate variance, confirm < X Å"

3. **Clarify Unsanitized Molecule Guarantees:** Define what properties are preserved:
   - Implicit hydrogens: counted or not?
   - Aromatic flags: set by RDKit or deferred to Molassembler?
   - Chiral tags (@/@@ from original SMILES): preserved through `sanitize=False`?

---

## 2. Confidence Mandate Analysis

### Clarifying Questions

1. **Why 8/10 and not 9/10?** The spec says "No critical" open questions, yet lists multiple architectural unknowns:
   - Determinism not proven
   - Timeout handling TBD
   - Test coverage limited
   - Regex error modes not specified
   
   Is 8/10 accounting for these, or is confidence actually lower?

2. **"Error boundaries are explicit" — which boundaries?** The spec mentions vertex index validation but not:
   - Regex extraction failures
   - RDKit parsing failures (e.g., invalid valence in unsanitized mol)
   - Molassembler instantiation failures (e.g., shape type not in library)
   - DG timeout or convergence failures

3. **Verification Protocol Determinism:** "Verification protocol is deterministic" — but RMSD is a floating-point computation. Is deterministic defined as "produces same RMSD value to bit precision" or "within tolerance"?

### What-If Scenarios

1. **Confidence Decay Under Novel Ligand Types:** Spec is validated on TiCat1/3/4 (Ti with eta-Cp ligands). What happens with:
   - **Ir(ppy)₃:** Ir with phenylpyridine (no eta bonds, but the memory notes this has atom ordering issues already)
   - **Re(CO)₅:** Rhenium with CO ligands (no eta bonds, but linear CO geometry)
   - **Ni-Phosphine:** Nickel with bidentate phosphine (requires dihedral sampling, not single DG)
   
   **Risk:** Confidence drops to 5/10 on novel geometry types; success metrics claim 100% coverage but only 3 shapes tested.

2. **Timeout as a Blocker:** DG timeout is "TBD" — this is a known open issue. If timeout is not specified, the code cannot be deployed:
   - Implement without timeout → infinite waits on complex molecules
   - Implement with arbitrary timeout → unpredictable failures
   - **Scenario:** Timeout=30s works fine until user passes a 500-atom ligand; timeout triggers; user gets `TimeoutError` instead of 3D structure
   - **Confidence impact:** This pushes confidence to 6/10 (known blocker).

### Points for Improvement

1. **Disaggregate Confidence Score by Component:** Current 8/10 is a blunt average. Break down:
   - Regex extraction: 9/10 (well-defined patterns, but silent failures not handled)
   - RDKit unsanitized parsing: 7/10 (unknown edge cases with valence, implicit H)
   - Molassembler instantiation: 8/10 (API is proven, but version stability unclear)
   - DG determinism: 5/10 (not measured, hardware-dependent)
   - Test coverage: 6/10 (only 3 fixtures, shape generalization unknown)
   - **Weighted confidence: ~7/10** (higher than stated; likely underestimated risk surface)

2. **Establish Determinism Baseline Before Implementation:**
   - [ ] Run DG 10x on each TiCat fixture; compute coordinate variance
   - [ ] Run on 3 different hardware platforms (CI Linux, local Mac, cloud instance)
   - [ ] Document variance and confirm < 0.1Å tolerance
   - If variance is > 0.1Å, rethink RMSD success metric (1.0Å is too tight)

3. **Resolve Timeout Handling as Precondition:**
   - Document exact timeout values (soft/hard limits)
   - Implement timeout in *this* PRD, not in follow-up
   - Specify retry/fallback behavior: fail fast or try alternative?

---

## 3. Scope Analysis

### Clarifying Questions

1. **"Assumes valid OIN format" — what validates it?** Scope says "Validation of OIN-SMILES syntax [is] pre-parser." But who calls the validator?
   - If `OIN3DGenerator.generate(oin_smiles)` calls it internally, then validation is in-scope
   - If callers are responsible, then callers will fail silently and blame the parser
   
   **Recommendation:** Validation is in-scope; move it from out-of-scope to in-scope.

2. **"Candidate artifact routing [is] out-of-scope"** — but the existing molassembler_adapter.py does template-based placement + DG fallback for multi-conformer routing. Is this logic being **completely** removed?
   - If yes, conformers are purely DG-based (deterministic, one output only)
   - If no, the existing template logic survives somewhere else
   
   **Clarification needed:** Where does template-based placement live after this refactor?

3. **"Support for non-standard polyhedral shapes beyond SCINE's shape library [is] out-of-scope"** — but what happens if an OIN-SMILES references an unknown shape?
   - **Scenario:** User passes `[Pd_UNKNOWN_SHAPE]` → regex extracts `UNKNOWN_SHAPE` → Molassembler's shape lookup fails
   - **Current error handling:** None specified; would be caught only when calling `mol.set_shape()`
   - **Scope creep risk:** Should this fail in regex validation (fast) or Molassembler instantiation (slow, opaque error)?

### What-If Scenarios

1. **Scope Boundary Violation: Template Logic Migration:** If template-based placement is removed but downstream code still calls it, blast radius is massive:
   - `OIN3DGenerator.generate()` currently may return template-based 3D structures
   - Callers may expect template geometry for specific complexes (e.g., Cisplatin uses known template)
   - **After refactor:** All complexes go through DG; if DG produces different geometry than template, RMSD changes
   - **User impact:** Existing integrations break silently if they rely on template positions

2. **Out-of-Scope Shape Support = Runtime Crashes:** Unknown shape codes are not pre-validated:
   - User typo: `[Pd_SQP]` vs. `[Pd_SQP_WEIRD]`
   - Regex extracts both fine
   - Molassembler lookup fails with opaque error: "Shape not in library" (vs. helpful "SQP_WEIRD not recognized")

### Points for Improvement

1. **Move Validation into In-Scope Explicitly:**
   ```python
   def validate_oin_syntax(oin_smiles: str) -> tuple[str, dict]:
       """Pre-validate OIN format before parser; fail fast with human errors."""
       # 1. Extract and validate shape codes against SCINE_SHAPES
       # 2. Extract and validate vertex indices against max per shape
       # 3. Ensure atom maps don't collide with new indices
       # Raises ValueError(f"Invalid shape {shape}: choose from {available}")
   ```

2. **Clarify Template Logic Fate:**
   - If removed entirely: Document breaking change for callers relying on template geometry
   - If moved to separate code path: Draw edge in architecture.yml from `atom_oin3d_generator` to both `atom_direct_parser` and `atom_template_engine`

3. **Pre-Validate Shape Library Membership:**
   - Build shape validation into regex extractor or separate validator
   - Fail with: "Unknown shape SQP_WEIRD. Available: [SQP, OC, TBP, ...]"

---

## 4. User Stories Analysis

### Clarifying Questions

1. **US-001: "Regex extracts all ... correctly"** — what does "correctly" mean? Spec has three regex patterns; each has potential ambiguities:
   - Shape pattern `_([A-Z0-9]+)` — no termination anchor. Does `_SQP_weird` match as `SQP_weird`? Or both `SQP` and `weird`?
   - Vertex pattern `\{([0-9,\s]+)\}` — matches `{1,2,3}` or `{1 , 2}` (spaces OK?). What about `{1,}` (trailing comma)?
   - Eta pattern `\|eta:([0-9,\s]+)\|` — what if multiple eta blocks exist? `|eta:1,2||eta:3,4|`? Does regex match first only?

2. **US-002: "Sequential mol.push_back() completes without errors"** — what is "without errors"?
   - RDKit errors (invalid atomic number)?
   - Molassembler errors (atom limit exceeded)?
   - Silent failures (Molassembler silently ignores bad atoms)?
   
   Spec should specify error type and recovery strategy.

3. **US-003: "Each eta-listed atom receives eta bond from metal center"** — but what if the metal center is not yet in the molecule?
   - **Scenario:** Parser loops atoms first (adds 10 atoms), then bonds. Eta block references atom at index 5. Atom 5 exists, but metal center is at index 0, which was already added.
   - **Question:** Are eta bonds added in order, or out-of-order? If out-of-order, what if metal center isn't added yet?

4. **US-004: "RMSD ≤ 1.0 Å for all tested complexes"** — but baseline was never measured. Previous memory says "TBD (measure v0.2.0 baseline)."
   - Is this acceptance criterion even testable without a baseline?
   - **New question:** What is the 1.0Å derived from? Literature? Empirical tolerance?

5. **US-005: "Reject invalid vertex indices with clear fatal errors"** — but what about partial molecule state?
   - **Scenario:** Parser adds 5 atoms successfully, then tries to add invalid vertex index on atom 6
   - **Current state:** Molassembler mol has 5 atoms, is inconsistent with intent
   - **Error message:** "Vertex index 99 invalid for SquarePlanar (max 3)" — clear, but what happened to the 5 atoms already added?

### What-If Scenarios

1. **Regex Pattern Ambiguity = Silent Corruption:** The patterns lack anchors and termination checks:
   - Input: `[Pd_SQP_WEIRD].[Cl]{1}|eta:1,2`
   - Shape pattern `_([A-Z0-9]+)` matches `_SQP` (greedy, first match)
   - Result: `_WEIRD` is not extracted, but not flagged as error
   - Constraint dict is incomplete; Molassembler instantiation uses wrong shape

2. **Multiple Eta Blocks Edge Case:** OIN format may allow `|eta:1||eta:2|` (two separate eta blocks):
   - Regex `\|eta:([0-9,\s]+)\|` matches first only
   - Second eta block is ignored silently
   - Molassembler lacks bonds for atoms in second block

3. **Atom Map Collision = Index Shift:** Original SMILES has atom maps; extraction adds new atom maps:
   - Input: `[Pd:5].[Cl:6]|eta:1,2`
   - Extraction replaces with: `[Pd:1].[Cl:2]` (new maps)
   - RDKit parses both: atom maps 1,2,5,6 are all present
   - Constraint dict keyed by extracted indices (1,2), but actual atoms are at indices 5,6
   - Molassembler shape injection targets wrong atom

4. **Unsanitized Molecule = Implicit H Bombs:** RDKit unsanitized parsing does not set explicit hydrogens:
   - Input: `[C]` (carbon, no explicit H)
   - RDKit parses as C(implicit: 4), Molassembler counts explicit bonds only
   - If Molassembler adds bonds based on atom, it may add extra bonds (valence error)

### Points for Improvement

1. **Regex Patterns Need Anchors and Validation:**
   ```python
   # Bad (ambiguous):
   shape_pattern = r'_([A-Z0-9]+)'
   
   # Good (anchored, greedy only to terminator):
   shape_pattern = r'_([A-Z0-9]+?)(?=[.\[\]])'  # Stop at next delimiter
   # Or validate shape membership after extraction
   ```

2. **Multi-Eta-Block Support:**
   - Test case: `[Pd].[Cl]{1}|eta:1||eta:2|`
   - Use `re.findall()` instead of `re.search()` to catch all blocks
   - Merge all eta indices into single set

3. **Atom Map Collision Resolution:**
   - Scan original SMILES for existing atom maps before extraction
   - Allocate new maps in disjoint range (e.g., if original has maps 1-10, new maps start at 11)
   - Document mapping in constraint dict: `{ 'old_to_new_atom_map': { 5: 11, 6: 12 } }`

4. **Transaction Semantics for Partial Molecules:**
   - Specify: "If any atom/bond operation fails, entire Molassembler molecule is discarded"
   - Implement via try/finally or copy-construct molecule
   - Document in API: "GeneratedStructure is all-or-nothing; no partial results"

---

## 5. Technical Specifications Analysis

### Clarifying Questions

1. **Regex Pre-processor: What Happens to Extracted Segments?** The spec says "Replace extracted segments with standard RDKit atom maps." But:
   - Where are extracted segments stored? Constraint dict is mentioned, but is it `dict[int, Any]` or `dict[str, str]`?
   - How are shape codes and vertex indices mapped back after RDKit parsing?
   - **Example:** OIN `[Pd_SQP].[Cl]{1}` → stripped SMILES `[Pd:X].[Cl:Y]` where X,Y are what values?

2. **Unsanitized RDKit Parsing: Implicit Hydrogens?** Calling `Chem.MolFromSmiles(smiles, sanitize=False)`:
   - RDKit still counts implicit hydrogens even without sanitization
   - But valence is not enforced
   - **Question:** If unsanitized SMILES has invalid valence (e.g., `[C-2]`), does RDKit parse it or reject it?

3. **Molassembler Atom/Bond Addition Order:** The spec iterates atoms first, then bonds. But:
   - Are bonds added in the same order as RDKit's mol.GetBonds()?
   - What if a bond references an atom index that hasn't been added yet?
   - **Scenario:** Bonds are in order [0-1, 1-2, 2-0]. Atoms are added in order [0, 1, 2]. Bond [0-1] is added after atom 1, so both endpoints exist. OK. But what if bonds are [0-1, 2-1, 1-0]? The second bond [2-1] is added before atom 2 is added.

4. **Shape Assignment Timing:** The spec assigns shape *after* atoms and bonds are added. But:
   - Does Molassembler allow shape assignment on an atom that already has bonds?
   - What if the number of bonds doesn't match the shape (e.g., square planar requires 4 bonds, but atom has 3)?
   - **Error handling:** Who validates shape compatibility with existing bonds?

5. **Eta Bond Addition: Indices Stability?** After atoms and bonds are added, eta bonds are added. But:
   - Are the indices in the eta block *original* indices (from OIN-SMILES) or *Molassembler* indices?
   - If RDKit stripped some atoms (implicit H?), do indices shift?

### What-If Scenarios

1. **Greedy Regex Matching = Over-Extraction:** The shape pattern `_([A-Z0-9]+)` is greedy:
   - Input: `[Pd_SQP_][Cl]` (shape code `SQP`, underscore, then `Cl`)
   - Pattern matches `_SQP_` (greedy to end)
   - Extracted shape: `SQP_` (includes trailing underscore)
   - Shape lookup fails: no shape named `SQP_`

2. **Implicit H Valence Bomb:** Unsanitized parsing allows invalid valence:
   - Input: `[C:1].[C:2]` (two carbons, no bonds)
   - RDKit parses each as C(implicit H: 4) 
   - When added to Molassembler, atom1 and atom2 each expect 4 explicit bonds
   - If only 1 bond [1-2] is added, atom1 and atom2 each have 3 missing bonds
   - Molassembler shape assignment fails: "Atom 1 has 1 bond but shape requires 4"

3. **Bond Index Out-of-Bounds:** If RDKit's bond order doesn't match atom addition order:
   - Atoms added: [Pd, Cl, Cl, N, N] (5 atoms, indices 0-4)
   - Bonds from RDKit: [(0,1), (0,2), (0,3), (0,4), (1,2)] — *last bond is between two ligands, not metal*
   - When adding last bond (1,2) to Molassembler, both endpoints exist, so no crash
   - But bond (1,2) is not an eta bond, so Molassembler may reject it (ligand-ligand bond without metal mediation)

4. **Shape Assignment After Bonds = Geometry Lock:** If shape is assigned after bonds:
   - Metal atom starts with 4 bonds (square planar)
   - `set_shape()` is called to enforce SquarePlanar geometry
   - But geometry may be incompatible with the bonds already present (e.g., bonds form tetrahedral, not square planar)
   - **Error:** Molassembler cannot deform bonds, only constrains shape; if shape & bonds are incompatible, who adjusts?

5. **Conformer Generation After All Constraints:** The spec calls DG *after* all atoms, bonds, shapes, and eta bonds are added:
   - Molassembler has now a fully constrained molecule
   - DG engine must satisfy all constraints (eta bonds, shape, standard bonds)
   - **What if constraints are over-constrained?** (e.g., shape + bonds + eta bonds form impossible geometry)
   - **Error:** DG may fail or produce geometry far from constraints

### Points for Improvement

1. **Detailed Atom/Bond Construction Protocol:**
   ```python
   def construct_molassembler_mol(atoms: list, bonds: list, shape_constraints: dict) -> masm.Molecule:
       """
       PRECONDITION: atoms are in RDKit order; bonds reference RDKit indices
       STEP 1: Create empty mol
       STEP 2: Add atoms (push_back preserves index order)
       STEP 3: Validate all bond endpoints exist
       STEP 4: Add standard bonds (single, double, triple, aromatic)
       STEP 5: Assign shape to metal atom (if shape_constraints has entry)
       STEP 6: Add eta bonds (now atom indices are stable)
       STEP 7: Call DG engine (all constraints set)
       POSTCONDITION: mol is ready for coordinate generation
       ERROR: If any step fails, raise with context (step #, atom index, bond index, shape)
       """
   ```

2. **Explicit Bond Validation:**
   ```python
   def validate_bonds_against_atoms(bonds: list, atom_count: int):
       """Fail fast on out-of-bounds bond indices."""
       for i, j, bond_type in bonds:
           if i >= atom_count or j >= atom_count or i < 0 or j < 0:
               raise ValueError(f"Bond ({i},{j}) has invalid endpoints; atom count is {atom_count}")
   ```

3. **Shape-Bond Compatibility Check:**
   ```python
   def validate_shape_compatibility(shape: str, bond_count_on_metal: int):
       """Ensure metal atom has expected number of bonds for shape."""
       expected_bonds = SHAPE_COORDINATION_NUMBER[shape]
       if bond_count_on_metal != expected_bonds:
           raise ValueError(f"{shape} expects {expected_bonds} bonds; metal atom has {bond_count_on_metal}")
   ```

4. **Deterministic DG Seeding:**
   ```python
   # Document seed strategy
   DG_SEED = 42  # Fixed seed for deterministic conformers
   # Verify: Run 10x on same molecule, compare coordinates
   # Document variance
   ```

---

## 6. Negative Constraints Analysis

### Clarifying Questions

1. **"DO NOT sanitize the RDKit molecule after AST extraction"** — but what if downstream code expects sanitized input?
   - The spec claims unsanitized preserves haptic topology
   - But SMILES generation via `.GetSmiles()` may behave differently on sanitized vs. unsanitized
   - **Question:** What SMILES is returned by unsanitized mol vs. sanitized mol for `[Pd].[Cl]{1}`?

2. **"DO NOT assume external callers"** — but `SMILESToXYZ` is a public API class. How is this "internal"?
   - `src/oinsmiles/__init__.py` exports `SMILESToXYZ`
   - Users call `SMILESToXYZ().convert(oin_smiles)`
   - This *is* external
   - **Contradiction:** Either design for internal use only (hide from API) or design robustly for external use

3. **"DO NOT perform candidate artifact routing"** — but what if DG fails or produces bad geometry?
   - Current molassembler_adapter has template + DG fallback
   - New parser removes fallback
   - **What happens if DG fails?** Raise exception? Return None? Return best-effort structure?

4. **"DO NOT modify OIN format specification"** — but regex changes how OIN is *parsed*. Does this change the spec?
   - Spec says regex patterns are `_([A-Z0-9]+)`, `\{([0-9,\s]+)\}`, etc.
   - If spec is silent on *how* shape codes are formatted, then spec is unchanged
   - But if a downstream tool relies on different parsing logic, it breaks

### What-If Scenarios

1. **Sanitization Dependency Hidden in SMILES Generation:** Unsanitized mol → RDKit's `.GetSmiles()`:
   - Unsanitized: `[Pd].[Cl]` might generate SMILES without aromatic flags
   - Sanitized: `[Pd].[Cl]` generates proper SMILES with valence checks
   - If SMILES is fed back into another tool (e.g., OIN validation), it may fail

2. **External Caller Violates "Internal Use Only":** 
   - User imports `SMILESToXYZ` from public API
   - Calls `convert(oin_smiles)` with complex novel ligand
   - Molassembler fails with opaque error: "Shape not in library" or "DG timeout"
   - **User impact:** API is public but error handling is internal-only; user has no recourse

3. **DG Failure Without Fallback:**
   - Input: `[Pd].[Cl]{1}|eta:1,2,3` (but only 3 atoms total, vertices for SquarePlanar require 4)
   - Molassembler shape assignment succeeds (no validation yet)
   - DG engine starts, realizes constraint is over-determined
   - **Without fallback:** DG fails; return GeneratedStructure with None for mol; user gets empty result

4. **Negative Constraint Contradiction = Silent Failure:**
   - "DO NOT sanitize" + "DO NOT assume external callers" = implies unsanitized is OK externally
   - But unsanitized may break downstream SMILES-based tools
   - **Result:** User's tool chain breaks silently; blames parser, not understands constraint violation

### Points for Improvement

1. **Reconcile "Internal Use Only" with Public API:**
   - Option A: Move `OIN3DGenerator.generate()` to internal module, hide from `__init__.py`
   - Option B: Design for robust external use; remove "internal" from constraints
   - **Recommendation:** Option B — public API must be defensive
   - Add constraint: "DO NOT assume OIN format is valid; validate and fail fast"

2. **Document Sanitization Semantics:**
   ```python
   # In spec: when is unsanitized mol acceptable?
   # Example:
   # - perception_tmc.py: sanitizes with Chem.SanitizeMol() — expects valid valence
   # - OINInlineHandler: calls Chem.MolToSmiles() — expects sanitized for canonical SMILES
   # - ChiralityRecoveryUtility: reads _OIN_CIPCode property — expects property set by CIPAssigner
   #
   # This parser: produces unsanitized mol
   # Downstream impact: any code expecting sanitized input must be adapted
   ```

3. **Fallback Strategy for DG Failure:**
   ```python
   # Define behavior:
   # 1. If DG succeeds: return GeneratedStructure(xyz, mol)
   # 2. If DG times out: raise TimeoutError (let user handle timeout)
   # 3. If DG fails (over-constrained, etc.): raise ConstraintError (list violations)
   # 4. NO fallback to de-aromatization or heuristics
   ```

4. **Explicitly Exempt Negative Constraints from External Callers:**
   - Add new constraint: "Negative constraints apply to implementation only; callers may assume sanitized output via defensive wrapping"
   - Or: guarantee unsanitized mol is wrapped in sanitization layer before returning to API

---

## 7. Risks & Mitigation Analysis

### Clarifying Questions

1. **Molassembler Version Compatibility: Upper Bound but No Strategy?**
   - Spec pins `>= 2.0.0, < 5.0.0` (wide range, covers 2.x, 3.x, 4.x)
   - Mitigation says "add CI test for version compatibility"
   - **But what test?** Does it:
     - Run existing tests against all versions?
     - Just check that package imports?
     - Actually call `set_shape()` and `add_bond()` on all versions?

2. **"RMSD Exceeds 1.0 Å: Document Scope"** — this is not a mitigation, it's a deferral.
   - What happens on deployment if a novel complex produces RMSD = 1.5Å?
   - Is 1.0Å hard requirement or soft guideline?
   - **Missing:** Escalation path. Who decides if 1.5Å is acceptable? What's the SLA?

3. **"Conformer Generation Timeout: Policy TBD"** — this is a blocking blocker.
   - Timeout policy is mentioned twice as TBD (in risks section and success metrics)
   - **Before implementation:** Must resolve timeout values (soft/hard limits), retry strategy
   - Spec cannot be "ready for implementation" with TBD timeout

4. **"Invalid OIN Format: Parser Assumes Pre-Validated"** — but who validates?
   - Risk: "Regex extraction fails on malformed OIN-SMILES"
   - Mitigation: "Document that parser assumes pre-validated OIN format"
   - **This is not a mitigation, it's passing the buck.** If callers don't validate, failures are silent.

### What-If Scenarios

1. **Version Compatibility: Silent API Drift:**
   - Spec is written for Molassembler 2.0.0
   - User installs 4.8.0 (within range)
   - API for `set_shape(metal_idx, SCINE_SHAPE)` changed to `set_shape(metal_idx, SCINE_SHAPE, options={})`
   - Code compiles (optional parameter), but silently ignores options
   - **Result:** Geometry changes; RMSD increases; user is confused

2. **RMSD Variance = Tolerance Creep:**
   - Test fixture TiCat1 has RMSD = 0.8Å (pass)
   - Deployment: user's complex produces RMSD = 1.3Å (fail? or pass with warning?)
   - Without clear threshold, team argues: "1.0Å is soft limit" vs. "hard requirement"
   - **Result:** Inconsistent behavior; no SLA

3. **Timeout Threshold Too Tight:**
   - Spec says "< 10 seconds per complex"
   - A 100-atom ligand takes 15 seconds
   - Timeout triggers at 10s; user gets `TimeoutError`
   - **Escalation:** Team adjusts timeout to 30s, discovers new regressions on other machines
   - **Result:** No deterministic timeout; per-hardware tuning required

4. **Malformed OIN = Silent Corruption:**
   - User typo: `[Pd_SQP][Cl]{99}` (index 99 out of bounds for SquarePlanar, max 3)
   - Regex extracts shape `SQP`, indices `[99]`
   - No validation happens (precondition error in spec)
   - Molassembler shape assignment called with metal_idx=0, shape=SQP (OK)
   - Eta bond assignment called: `mol.add_bond(0, 99, BondType.eta)` — **out of bounds crash or silent fail**
   - **Result:** User blames parser, doesn't know input was malformed

### Points for Improvement

1. **Version Compatibility CI:**
   ```python
   # In CI, test against multiple versions:
   for version in ['2.0.0', '3.1.0', '4.8.0']:
       # Install version
       # Run: OIN3DGenerator.generate(oin_smiles) on each test fixture
       # Check: RMSD <= 1.0 Å
       # Log: any API deprecation warnings
   ```

2. **RMSD Threshold: Hard SLA, Not Soft Guideline:**
   - Define: "Success means RMSD ≤ 1.0Å OR documented exception"
   - For exceptions, add test case to fixtures and re-baseline
   - Document variance range: "RMSD typically 0.6–0.9Å; outliers documented"

3. **Timeout Mitigation: Resolve Before Implementation:**
   ```python
   # Spec must include:
   TIMEOUT_SOFT_LIMIT = 30  # seconds; log warning if exceeded
   TIMEOUT_HARD_LIMIT = 60  # seconds; raise TimeoutError
   
   # Test on actual hardware:
   # 1. Baseline DG time for TiCat1/3/4 on CI (Linux)
   # 2. Baseline on user environments (estimate 2x-10x variance)
   # 3. Set soft limit = 2x max baseline; hard limit = 10x baseline
   # 4. Document: "If timeout triggers, complex geometry is not supported"
   ```

4. **Input Validation = Not Optional:**
   ```python
   # Move validation from out-of-scope to in-scope:
   def validate_oin_smiles(oin_smiles: str) -> None:
       """Fail fast on malformed OIN before Molassembler touch."""
       # 1. Extract shape codes; validate against SHAPE_LIBRARY
       # 2. Extract vertex indices; validate against shape's max vertices
       # 3. Reject with human-readable error (not silent)
       if not is_valid_oin(oin_smiles):
           raise ValueError(f"Invalid OIN format: {oin_smiles}. Valid shapes: {SHAPE_LIBRARY.keys()}")
   ```

---

## 8. Success Metrics Analysis

### Clarifying Questions

1. **"100% of valid ... parse without error"** — Tautology?
   - If "valid" is pre-defined as "passes validation," then 100% is guaranteed by construction
   - **Real metric:** "99% of OIN-SMILES produced by `XYZToSMILES` parse without error"
   - But success metrics don't reference `XYZToSMILES`, only "OIN-SMILES strings"

2. **"RMSD ≤ 1.0 Å for all tested complexes"** — but baseline was never measured.
   - Spec says "measure v0.2.0 baseline" (TBD from memory)
   - **How can RMSD be < 1.0Å if baseline is unknown?**
   - Is 1.0Å derived from literature, empirical tolerance, or arbitrary?

3. **"OIN→XYZ→OIN preserves chemical graph"** — but this is not a new success metric, it's an existing requirement.
   - Current v0.2.0 has round-trip validation (memory says "integration tests for round-trip")
   - Spec doesn't add new round-trip tests; just claims "preserves" without proof
   - **Metric is untestable without actual round-trip execution**

4. **"Conformer generation completes in < 10 seconds per complex"** — but timeout is TBD, and baseline is not measured.
   - How can "< 10 seconds" be a metric if no baseline exists?
   - Is 10 seconds a hard SLA or a guess?

5. **"Internal OIN3DGenerator interface remains unchanged"** — but return type *did* change from `str` to `GeneratedStructure`.
   - Memory says: "return type changed from `str` to `GeneratedStructure`"
   - Is this API change "internal" (allowed) or external (breaking)?
   - **Metric is misleading**

### What-If Scenarios

1. **Metric 1 = Circular Logic:**
   - Test: "Run parser on 100 valid OIN-SMILES; 100% parse without error"
   - Result: Pass (tautology)
   - **What about invalid OIN?** No metric covers error handling
   - **Failure:** Metric doesn't catch regression (malformed OIN → silent corruption)

2. **Metric 2 = Baseline Assumption Without Proof:**
   - Test: "RMSD ≤ 1.0Å for TiCat1, TiCat3, TiCat4"
   - But: no ground-truth geometry is documented
   - **What if DG produces different geometry than old method (ETKDG)?** New RMSD might be 0.9Å (pass) but geometry changed unexpectedly
   - **Failure:** Metric doesn't catch silent behavioral change

3. **Metric 3 = Dependent on External Roundtrip Code:**
   - Test: "OIN→XYZ→OIN round-trip" passes
   - But: round-trip depends on `xyz2mol` (unchanged), `OINInlineHandler` (unchanged), etc.
   - If this spec's parser doesn't break those, round-trip still works
   - **Metric is redundant; doesn't test *this* parser**

4. **Metric 4 = Missing Baseline = Meaningless:**
   - Test: "Conformer generation completes in < 10 seconds"
   - Baseline: unknown (timeout policy TBD)
   - **If timeout is 5 seconds, metric fails. If 30 seconds, metric is loose.**
   - **Failure:** Metric doesn't establish performance boundary; just claims < 10s without justification

5. **Metric 5 = Contradiction:**
   - Spec says: "Internal OIN3DGenerator interface remains unchanged"
   - Memory says: "return type changed from `str` to `GeneratedStructure`"
   - **Metric fails by contradiction**

### Points for Improvement

1. **Replace Tautological Metric with Coverage Metric:**
   ```python
   # Old: "100% of valid OIN-SMILES parse without error"
   # New: "100% of OIN-SMILES produced by XYZToSMILES (v0.2.0) parse without error"
   # Test: for each complex in fixtures (Cisplatin, Ferrocene, fac-Ir(ppy)₃, TiCat1/3/4, PdCl2-R-BINAP, etc.):
   #   1. xyz2mol → OIN string
   #   2. OINParser.generate → parsed_oin
   #   3. Assert no exception raised
   ```

2. **Baseline RMSD Measurement as Precondition:**
   ```python
   # Before finalizing metric, measure:
   # 1. Ground-truth geometry (literature, prior experiments)
   # 2. DG-generated geometry for each fixture
   # 3. Compute RMSD between ground-truth and DG
   # 4. Document variance range (e.g., 0.5–0.9Å typically)
   # 5. Set threshold = max(observed) + 20% tolerance = 1.1Å
   # 6. Success metric: "RMSD ≤ 1.1Å for all fixtures; variance < 0.3Å"
   ```

3. **Explicit Round-Trip Metric:**
   ```python
   # For each fixture:
   # 1. xyz2mol → OIN1
   # 2. OINParser + Molassembler → XYZ2
   # 3. xyz2mol(XYZ2) → OIN2
   # 4. Assert: OIN1 == OIN2 (or OIN1.canonical() == OIN2.canonical())
   # Metric: "100% of fixtures pass round-trip within tolerance"
   ```

4. **Performance Metric Tied to Baseline:**
   ```python
   # Measure on CI:
   # 1. Baseline DG time for each fixture (TiCat1/3/4)
   # 2. Compute mean: e.g., 3.5 seconds
   # 3. Set threshold = mean * 3 = 10.5 seconds (conservative)
   # Metric: "DG completes within 3x baseline time (measured on CI); no timeouts on test fixtures"
   ```

5. **Fix API Change Contradiction:**
   ```python
   # Metric 5 revised:
   # "OIN3DGenerator.generate() signature unchanged; return type GeneratedStructure = v0.2.0 contract"
   # (Acknowledge that internal implementation changes; API contract is stable)
   ```

---

## Final Triage Summary

### 🔴 Critical Blockers (Resolve Before Implementation)

| Item | Blocker | Action |
|------|---------|--------|
| **Timeout Handling** | TBD in spec; cannot deploy without timeout policy | Specify soft/hard limits (e.g., 30s/60s); measure baseline on CI |
| **Determinism Guarantee** | SCINE DG seeding not documented; hardware variability unknown | Measure coordinate variance across 3 machines; document seed/variance < 0.1Å |
| **RMSD Baseline** | Success metric assumes 1.0Å threshold without baseline; baseline measurement is TBD | Measure ground-truth vs. DG on TiCat1/3/4; set threshold = max(observed) + tolerance |
| **Validation Error Handling** | Regex assumes valid OIN; malformed input crashes silently | Add `validate_oin_syntax()` function; fail fast with human-readable errors |

### 🟡 High-Priority Improvements (Defer to MiniPRD/Review)

| Item | Risk | Action |
|------|------|--------|
| **Regex Pattern Ambiguity** | Over-extraction or silent misses on edge cases | Anchor patterns; test against corpus of real OIN-SMILES |
| **Molassembler Version Compatibility** | API drift between 2.x, 3.x, 4.x | Add CI matrix test for versions 2.0, 3.1, 4.8 |
| **Unsanitized Molecule Semantics** | Downstream code may expect sanitized input | Document guaranteed properties (valence, implicit H, aromaticity); test integration points |
| **Atom Map Collision** | Index shifts if original SMILES has atom maps | Implement collision detection; allocate new maps in disjoint range |
| **Partial Molecule State** | If construction fails partway, mol is inconsistent | Specify transaction semantics (all-or-nothing); implement try/finally |
| **Public API vs. "Internal Use"** | Contradiction — SMILESToXYZ is exported, not internal | Reconcile: either hide from API or design for robust external use |
| **Shape Compatibility Validation** | Shape + existing bonds may be incompatible | Add pre-DG validation: shape coordination number vs. actual bond count |
| **DG Over-Constraint Handling** | If constraints are impossible, DG fails with opaque error | Define fallback (fail fast? return best-effort?) or prove constraints are always compatible |

### 🟢 Ready-for-MiniPRD (Clarifications Noted)

- Regex pre-processor (with anchors, collision detection, validation)
- RDKit unsanitized parsing (document implicit H semantics)
- Molassembler instantiation (with transaction semantics)
- Integration into `OIN3DGenerator.generate()` (no breaking changes)
- Verification against TiCat1/3/4 fixtures (baseline measured first)

---

## Next Steps

**User Action:** Review this RedTeam_Report, prioritize blockers, and invoke `/hyper-resolve` to triage findings into final SuperPRD.

**Red Team Stance:** This proposal is architecturally sound but operationally incomplete. The core mutation (replace ETKDG with direct Molassembler parser) is well-motivated. However, 7 critical unknowns must be resolved:
1. ✓ Timeout policy
2. ✓ Determinism proof
3. ✓ RMSD baseline
4. ✓ Input validation
5. ✓ Molassembler version compatibility
6. ✓ Unsanitized molecule guarantees
7. ✓ Public API semantics

**Confidence Reassessment:** 8/10 → **6/10 pending blockers** (drops when concrete unknowns are exposed). Once blockers are resolved, confidence rises to **8/10+**.

---

**Report Generated By:** Red Team Agent  
**Status:** Ready for `/hyper-resolve` triage
