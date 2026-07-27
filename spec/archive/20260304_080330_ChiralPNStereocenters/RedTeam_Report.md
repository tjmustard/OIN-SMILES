# Red Team Report: Chiral P/N Stereocenter Support — Pass 2
**Analyzed:** `spec/active/Draft_PRD.md`
**Blast Radius Source:** `spec/compiled/architecture.yml` (36 nodes)
**Date:** 2026-03-04
**Pass:** 2 of 2 (architecture.yml now available; PRD updated to make Molassembler sole backend)

---

## Executive Summary

The PRD incorporates the major architectural decision from Pass 1 (Molassembler as sole backend, Architector removed, `has_chiral_centers()` routing eliminated). This resolves C-3 from Pass 1. However, **3 new critical blockers** have been introduced or clarified by the update. The most severe is a direct internal contradiction across four sections describing the Zone A encoding strategy — the data flow diagrams still describe the rejected pseudo-atom approach, including an invalid element symbol that will cause a runtime failure. The second critical issue is that `OINInlineHandler.parse_inline_string()` will silently corrupt `@`/`@@` via RDKit canonicalization, and no fix is specified. A third critical issue is that the sanitization prerequisite for `CIPAssigner.assign_all()` is missing, meaning CIP assignment will silently return no results on the current pipeline's output.

**Severity distribution:** 3 Critical · 4 High · 5 Medium · 4 Low

**Progress from Pass 1:**
- ✅ C-3 (routing detection broken) — Resolved by Molassembler-for-all decision
- ✅ H-1 (Molassembler output interface incompatibility) — Partially resolved by `GeneratedStructure` typed wrapper
- ✅ H-3 (blast radius missing `translator.py` / `inline.py`) — Now included in §5.2
- ⚠️ C-1 (Oganesson wrong R/S) — Partially addressed in Resolved Trade-offs (§5.1) and §6, but **data flow diagrams and API contracts still contradict the fix**
- ⚠️ C-2 (@/@@ destroyed before transplant) — Partially addressed in §5.1 Resolved Trade-offs, but **no fix specified for the OINInlineHandler RDKit round-trip corruption**
- ❌ M-4 (SCINE version pin / performance budget) — Still missing from §5.5 and §7

---

## Finding Registry

| ID | Severity | Section | Summary |
|---|---|---|---|
| **C-1** | **Critical** | §1.2A / §5.1 / §5.4 / §6 | Four sections describe mutually incompatible Zone A encoding strategies. Data flow diagrams use `[Zz]` pseudo-atom. Resolved Trade-offs prohibit it. §6 explicitly bans Oganesson. API contracts define `PseudoAtomStrategy(PSEUDO_ATOMIC_NUM=118)`. Implementers have no authoritative source. |
| **C-2** | **Critical** | §5.1 / §7 | `OINInlineHandler.parse_inline_string()` round-trips SMILES through `MolFromSmiles → MolToSmiles(canonical=True)`. This silently corrupts `@`/`@@` via atom reordering. Risk acknowledged in §7 but no fix is specified — only "verify empirically." |
| **C-3** | **Critical** | §5.1 / §5.4 | `CIPAssigner.assign_all()` will produce silent no-op results because the current pipeline returns a mol with `with_stereo=False`. RDKit's `AssignStereochemistry()` requires a sanitized mol. This is not addressed anywhere in the PRD. |
| **H-1** | High | §1.1 / §5.3 | Rh-BINAP test fixture is chemically invalid for Zone A chiral P validation. BINAP's chirality is axial (out-of-scope per §3.2), not P-centered. BINAP P atoms have no independent CIP assignment from RDKit. This test will silently pass the "best-effort" acceptance criterion while validating nothing. |
| **H-2** | High | §5.4 | Timeout implementation for Molassembler is unspecified. If Molassembler's C++ bindings hold the GIL during distance geometry, `threading.Thread.join(timeout)` cannot interrupt the call. `concurrent.futures.ProcessPoolExecutor` is required but not mandated. |
| **H-3** | High | §5.1 | xTB invocation mechanism is undefined in `MolassemblerAdapter`. xTB is "optional" in §5.5 but appears as a mandatory step in the §5.1 decode flow. Behavior when xTB is absent is unspecified. |
| **H-4** | High | §4 | US-001 acceptance criterion uses RDKit's re-derived CIP from the same AC2mol output as the circular reference. If `chiral_stereo_check()` sets wrong stereo flags, RDKit CIP on the resulting Mol is also wrong, and the test passes with a wrong answer. No external ground truth is specified. |
| **M-1** | Medium | §5.2 | Blast radius table assigns "@/@@ survival verification" to `generation/oin_parser.py`. The `@`/`@@` corruption actually occurs in `oin/inline.py` (`OINInlineHandler`). The downstream dependency `oin/parser.py` (`OINParser_oin`) is also not listed, despite depending on `OINInlineHandler`. |
| **M-2** | Medium | §5.5 | `scine-molassembler` PyPI package name is unconfirmed. Actual package may be `scine_molassembler` (underscore) or require the full SCINE meta-package. Spike must confirm before `pyproject.toml` is modified. |
| **M-3** | Medium | §5.5 | xTB (LGPL-3.0) may violate §6's "BSD/MIT/Apache only" constraint. LGPL is not in the permitted list. Either the constraint must be updated or xTB's licensing status clarified. |
| **M-4** | Medium | §8 | Success Metric 6 (120s SLA) has no regression baseline against retired Architector. A 3× slowdown passes the metric. No performance regression criterion defined. |
| **M-5** | Medium | §7 | `chiral_stereo_check()` in `perception_core.py` sets stereo flags before `CIPAssigner.assign_all()` runs. Both modify stereo on the same Mol. Precedence rule undefined; absent from blast radius. |
| **L-1** | Low | §5.4 | `GeneratedStructure.get_xyz_string()` has no format spec. Callers cannot depend on this without knowing: standard XYZ? Angstroms? atom count header? |
| **L-2** | Low | §5.4 | CLI missing `--seed`, `--no-xtb`, `--version`, and `--charge` for the `generate` subcommand (charge needed for xTB; seed for reproducibility). |
| **L-3** | Low | §4 | US-005 does not specify stdin support. `oin-smiles convert -` reading from stdin is standard UNIX CLI practice for computational chemistry pipelines. |
| **L-4** | Low | §8 | Metric 4 (≥1 Zone A chiral phosphine) is too low to distinguish "reliably works" from "got lucky." Should be ≥2 with independently verified CIP assignments. |

---

## § 1. Introduction & Goals — Analysis

### Clarifying Questions

1. **[C-1 origin]** §1.2A states the Zone A encoding mechanism as:
   > "replace the metal neighbor with a pseudo-atom `[Zz]` (heavy dummy)"

   Yet §5.1 Resolved Trade-offs states: *"Option B (pseudo-atom) is NOT used as the primary mechanism"* and §6 states: *"DO NOT use Oganesson (Z=118) or any other high-Z pseudo-atom."* And §5.4 defines `PseudoAtomStrategy.PSEUDO_ATOMIC_NUM: int = 118 # Oganesson`. Which of these four sections is authoritative?

2. §1.1 lists "Chelating chiral diphosphines (e.g., BINAP-Rh, chiral P bonded to metal)" as an in-scope example. §3.2 explicitly out-scopes "Axial chirality (atropisomers, e.g., BINAP's biaryl axial chirality)." BINAP's P atoms have no independent R/S assignment — their CIP designation derives entirely from the axial chirality. Will `CIPAssigner.assign_all()` return a `_CIPCode` for BINAP's phosphorus atoms? If not, the primary motivating example cannot be demonstrated.

3. §1.2 claims this feature provides "lossless" encoding for chiral complexes. But OIN strings generated by the *old* system for chiral complexes contain no `@`/`@@`. Users with legacy OIN databases will get different strings when re-encoding the same XYZ files. Is this backward-compatibility break documented anywhere? Is there a format version bump?

### What-If Scenarios

- **W1.1 — Amine inversion lability:** Chiral nitrogen stereocenters in non-bridgehead amines (e.g., tertiary amines in a ligand side chain) undergo rapid pyramidal inversion at room temperature. RDKit will assign `_CIPCode` to these atoms from 3D coordinates, but the designation is physically meaningless. The PRD provides no filtering for configurally labile N atoms.

- **W1.2 — Legacy OIN string database:** A user has 10,000 OIN strings for chiral TMCs generated by the old system (no `@`/`@@`). They run `generate()` on all of them. All are now processed by Molassembler (Molassembler-for-all). Molassembler ignores the absent `@`/`@@` and produces an arbitrary enantiomer. The user has no way to know their database is chemically incorrect. No migration path, no warning.

- **W1.3 — AC2mol misassigns bond order at chiral P:** If `AC2mol` assigns incorrect bond order to a P=O bond (treating phosphine oxide as phosphine), RDKit will not perceive P as tetrahedral. `CIPAssigner.assign_all()` returns no `_CIPCode`. Chirality is silently dropped. The failure is upstream of CIPAssigner, so the PRD's warning-on-failure logic never fires.

### Points for Improvement

- **P1.1 (Critical):** Rewrite §1.2A Zone A encoding description to match Option C (pre-fragmentation CIP + `recover_chirality_tag()`). Remove all references to `[Zz]` and pseudo-atom substitution.
- **P1.2 (High):** Replace BINAP-Rh examples throughout with a genuine tetrahedral Zone A chiral phosphine, e.g., a Pd or Rh complex with a monodentate phosphine carrying three distinguishable substituents (not BINAP whose chirality is axial).
- **P1.3 (Medium):** Add a filtering rule: N atoms are only encoded with `@`/`@@` if they are configurationally stable (ring-constrained, bridgehead, or explicit inclusion list). Document this criterion.

---

## § 2. Confidence Mandate — Analysis

### Clarifying Questions

1. Open question: "What seed API does Molassembler expose?" — yet §5.4 already specifies `seed: Optional[int] = 42`. If the Spike determines no seed API exists, is the `seed` parameter silently ignored, or does it raise `NotImplementedError`? The behavior in the absence of a seed API must be defined now, not deferred to the Spike.

2. Open question: "Is `scine-molassembler` pip-installable?" — yet §5.1 Resolved Trade-offs and §6 mandate it as the only backend. If the Spike answer is "no," the entire architecture is blocked. The PRD needs a fallback decision rule for a "no" answer, not just "TBD."

3. A 7th open question is missing: "Does `ParsedOIN` carry enough information for `MolassemblerAdapter` to identify the local index of chiral P/N atoms in each ligand fragment?" `ParsedOIN.vectors` provides `atom_in_fragment_idx` for binding atoms only. The chiral P/N atom's local index in the fragment SMILES is not separately encoded. For a chelating phosphine with two binding atoms and one chiral P, MolassemblerAdapter has no explicit field identifying the chiral atom's local index.

### What-If Scenarios

- **W2.1 — Spike determines Molassembler does not support LIN:** LIN geometry is used for Ferrocene (bis-η⁵-Cp). If Molassembler cannot generate haptic η⁵-ring coordination, Ferrocene fails. The PRD says "scope must be adjusted" but provides no decision tree. This leaves the project with no generation path for one of its primary verified test cases and no documented resolution.

### Points for Improvement

- **P2.1 (High):** Define the Spike's explicit go/no-go decision tree: (a) ≥7/9 geometries → proceed; (b) 5–6/9 → add geometry-mapping MiniPRD; (c) <5/9 → halt and re-evaluate. Encode this in the Spike MiniPRD definition.
- **P2.2 (Medium):** Specify `seed` behavior when API absent: "emit `UserWarning('Molassembler seed not supported; output may be non-deterministic')`. Do not raise."
- **P2.3 (Medium):** Add open question 7 about `ParsedOIN` carrying chiral atom indices.

---

## § 3. Scope — Analysis

### Clarifying Questions

1. §3.1 requires Molassembler to support all 9 OIN geometry codes, including LIN. LIN in OIN represents bis-haptic η⁵ coordination (Ferrocene). Does Molassembler's concept of LIN match this? Molassembler's LIN geometry typically describes a 2-coordinate linear metal with two distinct monodentate ligands — not two delocalized ring centroids. This needs to be verified in the Spike explicitly.

2. §3.2 out-scopes "Metal-centered chirality (Λ/Δ octahedral, handled by OIN geometry tags)." The V3.6 winding direction tag (`>` / `<`) encodes Λ/Δ. When `MolassemblerAdapter` generates an OCT structure from an OIN string with winding, must it respect the winding? If not, all octahedral complexes with Λ/Δ assignments will produce incorrect enantiomers.

### What-If Scenarios

- **W3.1 — Geometry mismatch between OIN tag and Molassembler perception:** OINDiscreteAligner assigns geometry via RMSD fit (competitive detection). Molassembler assigns geometry via its own stereopermutator algorithm from connectivity. For a distorted octahedral complex that OINDiscreteAligner classifies as OCT, Molassembler may perceive it differently. No reconciliation mechanism is specified.

- **W3.2 — Shell metacharacters in OIN string via CLI:** OIN strings contain `{`, `}`, `[`, `]`, `@`, `;` characters. When passed as a shell argument to `oin-smiles generate "[Fe_LIN].c{0>}..."`, these characters may be interpreted by the shell before reaching the Python process. No quoting documentation or input sanitization is specified.

### Points for Improvement

- **P3.1 (High):** Explicitly define in §3.1 what "LIN geometry support" means for Molassembler: must it handle delocalized η⁵ ring centroids, or only strict 2-coordinate linear metals? Add this to the Spike checklist.
- **P3.2 (Medium):** Add to §3.1: "Verify Molassembler respects V3.6 winding direction (`>` / `<`) for OCT complexes." Λ/Δ fidelity is part of the existing OIN specification.
- **P3.3 (Low):** Add CLI input sanitization to §3.1 scope: reject `<xyz_file>` paths containing `..`; treat `oin_string` as a raw string literal (no shell interpolation).

---

## § 4. User Stories — Analysis

### Clarifying Questions

1. **US-001 circular validation:** Acceptance criterion — "CIP code matches RDKit assignment on input XYZ." But the RDKit Mol for the input XYZ is built by `get_tmc_mol()` → `chiral_stereo_check()`. If `chiral_stereo_check()` incorrectly assigns stereocenters (known weakness of distance-based graph stereo detection), RDKit's CIP assignment on that Mol is also wrong. The test system-under-test and the reference oracle share the same potentially-flawed code path. What external ground truth validates CIP correctness?

2. **US-004 atom-ordering assumption:** `OIN(XYZ₁) == OIN(Molassembler(OIN(XYZ₁)))` assumes `OINDiscreteAligner` produces the same canonical fragment SMILES regardless of input atom ordering. The existing Ir(ppy)₃ failure (~5Å RMSD from atom reordering) shows Molassembler-like backends can reorder atoms. Will `RDKit.MolToSmiles()` on a Molassembler-reordered structure produce the same canonical SMILES? This is an unproven assumption.

3. **US-002 systematic failure invisible:** Acceptance criterion includes "if assignment fails, warning logged and no @/@@ emitted." This means if `recover_chirality_tag()` systematically fails for all Zone A chiral P (due to a bug or RDKit API mismatch), every US-002 test passes — they just all emit warnings. There is no failure criterion for systematic Zone A CIP failure.

### What-If Scenarios

- **W4.1 — Meso ligand with two chiral P:** A bidentate phosphine has chiral P at both binding atoms with opposite configurations (meso compound). `CIPAssigner.assign_all()` assigns R to P₁ and S to P₂. `recover_chirality_tag()` must apply `@` and `@@` to different atoms in the same ligand SMILES. The PRD discusses a single chiral center per ligand in all examples. Polychiral ligands are not addressed.

- **W4.2 — US-004 passes with wrong enantiomer:** If `recover_chirality_tag()` applies `@@` when it should apply `@` (off-by-one in the CIP → SMILES mapping), the OIN string is stable (same wrong string each re-encoding) but chemically incorrect. US-004 passes without detecting the error.

### Points for Improvement

- **P4.1 (High):** US-001 must specify an external ground truth: CSD entry with published R/S, or CCDC-validated assignment. The reference cannot be derived from the same adjacency-matrix-based RDKit path.
- **P4.2 (High):** US-002 must add a negative acceptance criterion: "If `recover_chirality_tag()` produces warnings on ≥3 of the Zone A test cases, the test suite reports systemic failure — not a series of acceptable best-effort degradations."
- **P4.3 (Medium):** US-004 round-trip must be split: (a) OIN string stability (string-identity of `@`/`@@`), tested independently of (b) geometric fidelity (RMSD ≤ threshold vs. input XYZ).
- **P4.4 (Low):** US-005 must specify stdin support: `oin-smiles convert -` reads from stdin.

---

## § 5.1 Architecture & Resolved Trade-offs — Analysis

---

### [CRITICAL C-1] Four-Way Internal Contradiction on Zone A Encoding Strategy

The PRD contains mutually incompatible specifications of the Zone A chirality encoding mechanism across four separate sections:

| Location | Strategy Described |
|---|---|
| §1.2A data flow text | `"replace metal neighbor with pseudo-atom [Zz] (heavy dummy)"` |
| §5.1 encoding diagram | `"→ replace metal neighbor with pseudo-atom [Zz]"` |
| §5.1 Resolved Trade-offs | `"Option B (pseudo-atom) is NOT used as the primary mechanism"` |
| §5.4 API contract | `PseudoAtomStrategy.PSEUDO_ATOMIC_NUM: int = 118  # Oganesson` |
| §6 Negative Constraint | `"DO NOT use Oganesson (Z=118) or any other high-Z pseudo-atom"` |

An implementer reading §1.2A and §5.1 diagrams will implement pseudo-atom substitution. An implementer reading §5.1 Resolved Trade-offs and §6 will implement Option C (pre-fragmentation CIP + `recover_chirality_tag()`). An implementer reading §5.4 will implement `PseudoAtomStrategy` with Z=118. These are three different implementations of the same feature.

**Additionally:** `[Zz]` is not a recognized element in RDKit's periodic table. `Chem.MolFromSmiles("[Zz]")` returns `None`. Any code path using this symbol will silently fail at runtime.

**Mandatory resolution:** §1.2A encoding flow and §5.1 encoding diagram must be rewritten to describe Option C exclusively. The `PseudoAtomStrategy` class in §5.4 must either be (a) deleted, or (b) renamed to `ChiralityRecoveryStrategy` and reimplemented without pseudo-atom substitution.

---

### [CRITICAL C-2] `OINInlineHandler.parse_inline_string()` Corrupts `@`/`@@` via RDKit Canonicalization

The inline handler converts `N{0}` → `[N:1000]`, calls `Chem.MolFromSmiles(processed_frag)`, then generates clean SMILES via `Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)`.

This is a **known chirality corruption vector** in RDKit:

1. RDKit's canonical SMILES algorithm can reorder atoms relative to the input.
2. `@`/`@@` in SMILES are defined relative to the atom traversal order at write time. If RDKit reorders atoms during canonicalization, a `@` may become `@@`, or be omitted if perceived as unspecified.
3. Atom map numbers (used for slot injection) do not preserve chirality — they are a separate property.

This is noted in §7 as a risk: *"@/@@ markers destroyed during `OINInlineHandler.generate_inline_string()` injection step."* But §7 provides no fix — only: *"Verify empirically during Spike."*

Empirical verification without a specified fix is not a mitigation. Three possible fixes exist; exactly one must be selected and specified before implementation:

- **Option A:** Extract `@`/`@@` from the original fragment SMILES before slot injection; re-apply via `RDKit.SetStereoFrom3D()` after canonicalization (requires 3D coords to be passed through).
- **Option B:** Use `Chem.MolFromSmiles(smiles, sanitize=False)` and avoid re-canonicalization entirely; preserve atom order.
- **Option C:** Parse the original SMILES and canonical SMILES in parallel; map chiral atoms by element+neighbor topology and transfer `@`/`@@` explicitly.

**Mandatory resolution:** Select a fix option and document it. Add a unit test: *"apply `parse_inline_string()` to `[P@@H]{0}c1ccccc1`; verify output SMILES contains `@@`."*

---

### [CRITICAL C-3] `CIPAssigner.assign_all()` Will Produce Silent No-Op on Current Pipeline Output

The PRD specifies `CIPAssigner.assign_all()` is called after `get_tmc_mol()`. The current pipeline calls:

```python
tmc_mol, xyz_coords = get_tmc_mol(path, charge, with_stereo=False)
```

The `with_stereo=False` parameter means the returned Mol is NOT sanitized for stereochemistry. RDKit's `AssignStereochemistry()` (the basis of `CIPAssigner.assign_all()`) requires:

1. The mol must have aromatic perception completed (`SANITIZE_SETAROMATICITY`).
2. Ring information must be computed (`SANITIZE_FINDRADICALS`).
3. Valence must be assigned (`SANITIZE_PROPERTIES`).

Without these, `AssignStereochemistry()` silently returns without setting `_CIPCode` on any atom. No error is raised. The CIPAssigner call succeeds but does nothing. All `has_chiral_pn()` calls return False. All `@`/`@@` encoding is silently skipped.

**Mandatory resolution:** `CIPAssigner.assign_all()` must call `Chem.SanitizeMol(mol)` (or the minimum required flags: `SANITIZE_SETAROMATICITY | SANITIZE_SETHYBRIDIZATION | SANITIZE_SYMMRINGS`) before `AssignStereochemistry()`. This must be documented in the API contract.

---

### Remaining § 5.1 Clarifying Questions

4. The decoding flow shows: *"xTB optimization → refined structure"* and *"Post-optimization CIP check: if inverted → skip xTB, return pre-optimization structure."* But `MolassemblerAdapter.convert()` must capture the pre-xTB structure before calling xTB. `GeneratedStructure` has no `_pre_xtb_snapshot` attribute. How is the revert implemented?
5. The `timeout_seconds` parameter: the docstring applies it to "DG generation." But xTB is the computationally expensive step for large complexes. Does `timeout_seconds` wrap the entire `convert()` call (DG + xTB) or only DG?

### What-If Scenarios

- **W5.1.1 — Pre-fragmentation CIP on Mol with incorrect bond orders:** For the Ferrocene case, the full complex Mol may have bond orders that RDKit cannot sanitize correctly (metal-cyclopentadienyl bonds). If `Chem.SanitizeMol()` raises an exception or silently skips problematic atoms, `CIPAssigner.assign_all()` returns a partially-assigned Mol. Some `_CIPCode` values may be set; others absent. The pipeline proceeds without warning.

- **W5.1.2 — xTB inverts chirality, revert not possible:** If xTB is called synchronously inside `convert()` and the pre-xTB `GeneratedStructure` object is not retained in a local variable before xTB runs, the revert cannot occur. The API contract provides no mechanism to capture the pre-optimization state.

### Points for Improvement

- **P5.1.1 (Critical):** Rewrite §1.2A and §5.1 encoding diagrams to describe Option C only. Remove `[Zz]` and all pseudo-atom references from all sections.
- **P5.1.2 (Critical):** Specify the @/@@ fix for `OINInlineHandler`. Select one of the three options above and add it to the blast radius with a required unit test.
- **P5.1.3 (Critical):** Add mol sanitization to `CIPAssigner.assign_all()` API contract.
- **P5.1.4 (High):** Specify the `timeout_seconds` scope: it must wrap both DG and xTB. Specify the implementation mechanism: `concurrent.futures.ProcessPoolExecutor` is required if Molassembler holds the GIL.
- **P5.1.5 (High):** Add `_pre_xtb_atoms: Optional[List[Tuple[str,float,float,float]]]` attribute to `GeneratedStructure` to support xTB revert.

---

## § 5.2 System Graph Blast Radius — Analysis

### Clarifying Questions

1. The blast radius table lists `generation/oin_parser.py` as "Modified — verify @/@@ markers survive `parse_inline_string()`." However, `parse_inline_string()` lives in `oin/inline.py` (`atom_OINInlineHandler`), not `generation/oin_parser.py`. The modification is assigned to the wrong file.

2. `pyproject.toml` change: "remove `architector` if it is a direct dep." This is a conditional. The Spike must confirm whether `architector` is listed in `pyproject.toml`; the removal should be **mandatory if present**, not conditional.

3. `atom_chiral_stereo_check` in `perception_core.py` sets stereo flags during `AC2mol` → before `CIPAssigner.assign_all()` runs. If both set conflicting stereo information on the same Mol, the second caller wins — but which caller is authoritative? This interaction is not in the blast radius.

### What-If Scenarios

- **W5.2.1 — `chiral_stereo_check` conflicts with `CIPAssigner`:** The pipeline is: `xyz2AC → AC2mol → chiral_stereo_check → get_tmc_mol() returns Mol → CIPAssigner.assign_all()`. Both modify stereo properties on the same object. If `chiral_stereo_check` sets `@` and `CIPAssigner` computes `@@` for the same atom (because one uses 3D distance geometry and the other uses CIP algorithm on the full complex), the result is non-deterministic depending on call order.

- **W5.2.2 — `oin/parser.py` uncovered:** `atom_OINParser_oin` depends on `atom_OINInlineHandler`. Any change to `OINInlineHandler.parse_inline_string()` (required by C-2 fix) cascades into `atom_OINParser_oin`. This downstream impact is absent from the blast radius.

### Points for Improvement

- **P5.2.1 (Medium):** Correct the blast radius table: move the "@/@@ survival" modification entry from `generation/oin_parser.py` to `oin/inline.py`.
- **P5.2.2 (Medium):** Add `src/oinsmiles/utils/perception_core.py` to blast radius (read impact: define precedence between `chiral_stereo_check` and `CIPAssigner.assign_all()`).
- **P5.2.3 (Medium):** Add `src/oinsmiles/oin/parser.py` to blast radius as a downstream dependency of `OINInlineHandler`.
- **P5.2.4 (Low):** Change `pyproject.toml` "if it is a direct dep" to "confirm in Spike; removal is mandatory if present."

---

## § 5.3 Execution Checklist (MiniPRDs) — Analysis

### Clarifying Questions

1. The Spike has no defined STOP condition. If Molassembler supports only 4/9 OIN geometry codes, what happens? §7 says "scope must be adjusted or a geometry-mapping layer added" — but adding a geometry-mapping MiniPRD is not listed in §5.3, and there is no documented process for halting the PRD and revisiting the architecture.

2. The MiniPRD ordering is Spike → ChiralEncoding → MolassemblerAdapter → CLI → ChiralTests. ChiralEncoding modifies `OINSanitizer` (encoding path), which must be validated before MolassemblerAdapter (decoding path) is built on top of it. But the integration test between encoding and decoding only appears in MiniPRD 5 (ChiralTests). A bug in ChiralEncoding discovered during ChiralTests requires unwinding MiniPRDs 3 and 4.

### What-If Scenarios

- **W5.3.1 — Rollback scenario:** MiniPRD 3 removes `ArchitectorAdapter` and `ArchitectorWrapper`. Three weeks later, `scine-molassembler` releases a breaking API change. `MolassemblerAdapter` no longer compiles. The project has no generation backend. Architector has been deleted. No rollback plan is documented.

- **W5.3.2 — Spike succeeds on 7/9 geometries, PBP fails:** PBP (pentagonal bipyramidal, 7-coordinate) is unsupported by Molassembler. The PRD creates a 7-coordinate OIN→XYZ call that raises `MolassemblerGeometryError`. This is a correct, explicit failure — but it is still a scope reduction relative to the current Architector-based system. No user-facing documentation change is required by §5.3.

### Points for Improvement

- **P5.3.1 (High):** Define the Spike go/no-go gate (see P2.1). Without this, the Spike has no exit condition.
- **P5.3.2 (Medium):** Insert an encoding integration gate: after MiniPRD 2 (ChiralEncoding) and before MiniPRD 3 (MolassemblerAdapter), require that round-trip encoding tests pass with existing achiral fixtures. This catches ChiralEncoding regressions before the decoding backend is built.
- **P5.3.3 (Medium):** Define a rollback plan: Architector deprecation PR should not be merged until Molassembler passes all 5 existing achiral round-trip tests. Pin the deprecation to a green CI gate.

---

## § 5.4 API Contracts / Schema — Analysis

### Clarifying Questions

1. **`PseudoAtomStrategy.PSEUDO_ATOMIC_NUM: int = 118`** is defined in the API contract. §6 explicitly prohibits this. (See C-1.)

2. **`recover_chirality_tag()` placement:** It is a static method of `PseudoAtomStrategy` (pseudo-atom strategy class), but the Resolved Trade-offs describe it as part of Option C (not the pseudo-atom approach). Why is the function for the accepted Option C nested inside the class for the rejected Option B?

3. **`MolassemblerAdapter.convert()` Step 4 — "Run xTB optimization":** xTB is listed as "optional" in §5.5. What is the behavior when xTB is not installed? Return pre-xTB structure? Raise? Silently skip? Unspecified.

4. **`GeneratedStructure.get_xyz_string()`:** What format? Standard XYZ (N-line, comment-line, symbol x y z in Ångstroms)? No format contract means callers cannot depend on this method.

5. **`MolassemblerAdapter(seed=42)`:** If Molassembler exposes no seed API, what does `seed=42` do? Silent no-op? `UserWarning`? `NotImplementedError`?

### What-If Scenarios

- **W5.4.1 — GIL-holding Molassembler blocks timeout:** If Molassembler's C++ extensions hold the Python GIL during distance geometry, `threading.Thread.join(timeout=120)` cannot interrupt the blocked thread. The process hangs for minutes or until OOM. `concurrent.futures.ProcessPoolExecutor` spawns a subprocess (solves the GIL problem) but cannot share in-process RDKit Mol objects without serialization. The design of the timeout mechanism fundamentally constrains what data can be passed to/from Molassembler.

- **W5.4.2 — `write_file()` on bad path:** `GeneratedStructure.write_file(path)` raises `OSError` on a read-only path or full disk. The CLI calls this for `--output FILE`. No exception contract is defined. The CLI must catch `OSError` explicitly (violation of §7 low-risk: "CLI error handling: unformatted Python traceback").

- **W5.4.3 — `has_chiral_pn()` called after `OINSanitizer`:** `OINSanitizer.generate_robust_smiles()` calls `UpdatePropertyCache()` which may clear RDKit atom properties including `_CIPCode`. If `has_chiral_pn()` is called after sanitization, it returns False even for confirmed chiral atoms.

### Points for Improvement

- **P5.4.1 (Critical):** Remove `PseudoAtomStrategy` class entirely. Replace with `ChiralityRecoveryUtility` (standalone module-level functions, not a class) that implements Option C: `assign_from_full_complex(full_mol, fragment_mol, chiral_atom_full_idx) -> Tuple[str, str]` returns `(smiles_with_chirality, cip_code)`.
- **P5.4.2 (High):** Specify xTB-missing behavior in `MolassemblerAdapter.convert()` docstring: "If xTB is not installed, step 4 is skipped; pre-xTB Molassembler structure is returned with `UserWarning`."
- **P5.4.3 (High):** Specify timeout implementation in `MolassemblerAdapter.__init__()`: `concurrent.futures.ProcessPoolExecutor` with `future.result(timeout=timeout_seconds)`.
- **P5.4.4 (Medium):** Add format spec to `get_xyz_string()`: "Standard XYZ format: line 1 = atom count (int), line 2 = comment string, lines 3..N+2 = `<symbol> <x> <y> <z>` in Ångstroms."
- **P5.4.5 (Medium):** Specify seed-API-absent behavior: "emit `UserWarning('Molassembler seed not supported; output may be non-deterministic')`. Do not raise."
- **P5.4.6 (Low):** Add `--seed INT`, `--no-xtb`, `--version`, and `--charge INT` (for `generate` subcommand) to CLI options spec.
- **P5.4.7 (Low):** Add exception contract to `write_file()`: "Raises `OSError` on write failure. CLI callers must catch and format as user-facing error message."

---

## § 5.5 Dependencies — Analysis

### Clarifying Questions

1. Is the actual PyPI package name `scine-molassembler` (hyphen) or `scine_molassembler` (underscore)? PyPI package names are case-insensitive but underscore/hyphen matters for some pip versions. The Spike must confirm: `pip install scine-molassembler` actually installs the correct package.

2. `xtb ≥ 22.1` is listed as "optional (already in pyproject.toml)." xTB version 22.1 refers to the standalone `xtb` binary, not the `xtb-python` Python package. Which interface does `MolassemblerAdapter` use: subprocess call to the `xtb` binary on PATH, or the `xtb` Python package? This determines how the "optional" behavior is detected.

3. After `ArchitectorAdapter` is removed, does `mendeleev` (used only in `ArchitectorAdapter` for covalent radii lookup) become an orphaned dependency? The dependency list in §5.5 omits this audit.

### What-If Scenarios

- **W5.5.1 — LGPL conflict with §6 constraint:** The `xtb` binary is licensed LGPL-3.0. §6 states "Only BSD/MIT/Apache dependencies are permitted." LGPL is not in this list. If `xtb` is already in `pyproject.toml` (§5.5 says so), the license constraint in §6 is already violated by an existing dependency. Either the constraint must be updated to permit LGPL, or `xtb` must be optional with a documented fallback path.

- **W5.5.2 — ARM64 / Windows build failure:** SCINE's C++ extensions require compilation with Boost, Eigen, and potentially OpenMP. If a pre-built wheel is unavailable for the target platform, `pip install scine-molassembler` will attempt a source build and likely fail. This is a hard blocker since Architector has been removed and there is no fallback.

### Points for Improvement

- **P5.5.1 (Medium):** Add to Spike: "Confirm exact `pip install` command for Molassembler on Ubuntu 22.04 x86_64 in the project `.venv`. Record the installed version number for pinning."
- **P5.5.2 (Medium):** Clarify xTB interface (subprocess vs. Python package) and update "optional" handling specification.
- **P5.5.3 (Medium):** Audit `pyproject.toml` after Architector removal; confirm `mendeleev`, `architector`, and `xTB` are correctly updated.
- **P5.5.4 (Medium):** Resolve §6 license constraint vs. LGPL-licensed xTB: either update the constraint to "BSD/MIT/Apache/LGPL" or mark xTB as optional with skip behavior documented.

---

## § 6. Negative Constraints — Analysis

### Clarifying Questions

1. "DO NOT modify the OIN string syntax" — `@`/`@@` are standard SMILES. Adding them to OIN strings changes the encoding output for chiral complexes. For a user with a legacy OIN database, re-encoding the same XYZ files produces different strings. Is this constraint about syntax (no new OIN-specific tags) or output stability (same XYZ always → same OIN)? These are different constraints with different implications.

2. "DO NOT hardcode chiral P/N ligand lists" — sound in principle. But what IS the algorithmic detection criterion? Options: (a) any atom where `AssignStereochemistry()` sets `_CIPCode`, (b) any P or N with exactly 4 non-hydrogen substituents, (c) `RDKit.Chem.FindMolChiralCenters()`. These produce different sets. The implementation team needs the exact predicate.

3. The `PseudoAtomStrategy` class in §5.4 directly violates "DO NOT use Oganesson (Z=118)." No note in §6 clarifies that §5.4 should be updated. An implementer reading §6 after §5.4 will be confused about which is correct.

### What-If Scenarios

- **W6.1 — DO NOT retain Architector + Spike geometry failure:** If the Spike determines Molassembler cannot generate LIN (Ferrocene) or PBP (7-coordinate), the "DO NOT retain Architector" constraint leaves the project with no generation path for those geometries. The only permitted behavior is `MolassemblerGeometryError`, which means the system regresses for those cases vs. the current Architector baseline. No escape valve is defined.

### Points for Improvement

- **P6.1 (Medium):** Clarify the OIN syntax constraint: "DO NOT add new OIN-specific tags. `@`/`@@` are standard SMILES and are permitted. Existing achiral OIN strings are NOT guaranteed backward-compatible encoding with the new system."
- **P6.2 (Medium):** Add a positive directive to complement the DO NOTs: "DO implement Zone A chirality recovery exclusively via Option C (pre-fragmentation `CIPAssigner.assign_all()` + post-sanitization `recover_chirality_tag()`). No pseudo-atom substitution."
- **P6.3 (Low):** Specify the exact algorithmic chirality detection predicate: "`has_chiral_pn()` returns True if and only if `atom.GetAtomicNum() in (7, 15)` AND `atom.GetPropsAsDict().get('_CIPCode')` is not None, after `CIPAssigner.assign_all()` has been called on the full sanitized complex mol."
- **P6.4 (Low):** Add a cross-reference note in §6: "Note: §5.4 `PseudoAtomStrategy` class violates this constraint and must be removed (see C-1)."

---

## § 7. Risks & Mitigation — Analysis

### Clarifying Questions

1. The "atom ordering" risk is rated "Medium." The Ir(ppy)₃ Architector case already fails (~5Å RMSD) due to atom reordering — this is a documented, confirmed failure in the existing system. Molassembler will have the same or worse behavior because it generates structures from scratch. Why is this "Medium" rather than "High"?

2. "Verify canonical OIN string is identical regardless of Molassembler atom ordering" — what is the verification mechanism? If `OINDiscreteAligner` depends on `RDKit.MolToSmiles()` for canonical fragment SMILES, and Molassembler changes atom ordering within a fragment, the canonical SMILES changes, and the OIN string changes. The risk has no mitigation strategy.

### Missing Risks

| Risk | Severity | Gap |
|---|---|---|
| Thread safety of Molassembler C++ backend | Medium | Concurrent `generate()` calls (Jupyter parallel cells, pytest-xdist) may corrupt Molassembler's internal state. Not addressed. |
| SCINE project maintenance | Medium | `scine-molassembler` is academic software (ETH Zurich). Infrequent releases. If abandoned (no Python 3.13 wheel), project loses sole generation backend with no fallback since Architector is removed. |
| RDKit CIP algorithm version sensitivity | Low | CIP assignment has changed between RDKit minor versions. Test fixture `@`/`@@` ground truth must record which RDKit version produced the expected values. |
| `chiral_stereo_check` + `CIPAssigner` conflict | Medium | Both modify stereo on the same Mol. Precedence undefined. Not in risk table. |

### Points for Improvement

- **P7.1 (Medium):** Escalate "atom ordering" from Medium to High. Add mitigation: "During Spike, test whether `OINDiscreteAligner` produces identical OIN strings for the same structure with permuted atom ordering. If not, add atom-ordering normalization before `OINDiscreteAligner`."
- **P7.2 (Medium):** Add thread safety risk: "MolassemblerAdapter should document single-threaded use. Concurrent calls must be serialized by caller or protected with `threading.Lock` if the C++ backend is not thread-safe."
- **P7.3 (Medium):** Add SCINE maintenance risk with mitigation: "Pin to `scine-molassembler>=X.Y,<X+1`. If package shows no Python-version-compatible release within 18 months, initiate backend migration."
- **P7.4 (Low):** Add RDKit version sensitivity: "Fixture expected `@`/`@@` must record generating RDKit version. CI enforces `rdkit >= 2025.9.3` before running chiral tests."

---

## § 8. Success Metrics — Analysis

### Clarifying Questions

1. Metric 3: "≥ 95% of well-defined test cases" — the fixture set is 2–3 complexes. 95% of 2 = 1.9 ≈ 2. This is effectively "all cases must pass." The 95% qualifier is meaningless at this sample size. What constitutes "well-defined"?

2. Metric 4: "At least 1 Zone A chiral phosphine test case" — one case cannot distinguish a correct implementation from a lucky coincidence on a specific geometry.

3. Metric 2 (achiral regression): "Geometrically correct structures" — what is the quantitative definition of "correct"? The existing Architector baselines (Cisplatin 0.471Å, Ferrocene 0.586Å, etc.) are the implicit reference, but neither these values nor an RMSD threshold appear in §8.

4. Metric 6: "All `generate()` calls complete within 120 seconds" — with no baseline, a system that takes 119 seconds (vs. Architector's 30 seconds) passes this metric. No regression floor is defined.

### What-If Scenarios

- **W8.1 — Systematic Zone A failure passes Metric 4:** If `recover_chirality_tag()` works correctly for exactly 1 Zone A test fixture (the one in the test suite) but fails silently on all others, Metric 4 passes. The metric does not distinguish reliable from incidental correctness.

- **W8.2 — Molassembler 3× slower than Architector:** All generate() calls finish in <120s. Metric 6 passes. But the 120s SLA is 4× worse than Architector's typical runtime for Cisplatin (~30s). No user-facing SLA regression is caught.

### Points for Improvement

- **P8.1 (High):** Metric 2 must add an RMSD threshold: "Molassembler output for Cisplatin, Transplatin, Cis-PtCl₂(en), Ferrocene must achieve RMSD ≤ 2.0 Å vs. input XYZ after atom-order alignment."
- **P8.2 (High):** Metric 4 must require ≥ 2 Zone A chiral P test cases: one monodentate, one bidentate (chelating). Both must have independently verified CIP assignments (not self-validated by the same RDKit path).
- **P8.3 (Medium):** Metric 6 must add a regression baseline: "Molassembler `generate()` time ≤ 2× the Architector baseline measured on all 5 existing fixtures before Architector removal."
- **P8.4 (Medium):** Metric 3: define "well-defined test case" — "a complex where `RDKit.AssignStereochemistry()` on the manually-verified fragment returns a non-None `_CIPCode`, cross-referenced against CSD or published assignment."
- **P8.5 (Low):** Add Metric 8: "Thread safety smoke test — concurrent `generate()` from 4 threads on the Cisplatin fixture; verify no memory corruption or inconsistent output."

---

## Blast Radius from architecture.yml

**Architecture graph: 36 nodes (1 System · 5 Module · 30 Atomic)**

Nodes transitioning `clean → dirty` due to this PRD:

| Node ID | Current Status | File | Reason |
|---|---|---|---|
| `atom_OINSanitizer` | clean | `oin_aligner.py` | Zone A CIP recovery logic added; `recover_chirality_tag()` called post-sanitization |
| `atom_OINDiscreteAligner` | clean | `oin_aligner.py` | Must pass `@`/`@@` through to serialization without corruption |
| `atom_get_tmc_mol` | clean | `perception_tmc.py` | `CIPAssigner.assign_all()` inserted; mol sanitization required pre-call |
| `atom_get_oin_string` | clean | `perception_tmc.py` | Must preserve `@`/`@@` through inline generation path |
| `atom_OIN3DGenerator` | clean | `generation/engine.py` | Rewired exclusively to `MolassemblerAdapter` |
| `atom_OINInlineHandler` | dirty | `oin/inline.py` | C-2 fix required: `@`/`@@` preservation through RDKit round-trip |
| `atom_XYZToSMILES` | dirty | `core/translator.py` | Wired to `CIPAssigner` (already dirty TD-001) |
| `atom_SMILESToXYZ` | dirty | `core/translator.py` | Wired to `MolassemblerAdapter` (already dirty TD-003) |

Blast radius table in §5.2 is missing these nodes (gap M-1):
- `atom_OINParser_oin` (`oin/parser.py`) — downstream of `OINInlineHandler`
- `atom_chiral_stereo_check` (`perception_core.py`) — precedence conflict with `CIPAssigner`

Nodes removed (must be deleted from architecture graph):
- `atom_ArchitectorAdapter` (`generation/architector_adapter.py`)
- `atom_ArchitectorWrapper` (`generation/wrapper.py`)

New nodes to add to architecture graph after implementation:
- `atom_CIPAssigner` + `atom_recover_chirality_tag` (`core/chirality.py`)
- `atom_MolassemblerAdapter` + `atom_GeneratedStructure` (`generation/molassembler_adapter.py`)
- `atom_CLI_main` (`cli.py`)

**Impact: 8 existing nodes modified, 2 removed, 5 added = 15 of 36 nodes affected (42% of graph)**

---

## Mandatory Pre-Implementation Actions

These must be resolved **before any MiniPRD coding begins**:

1. **Rewrite §1.2A and §5.1 diagrams** to describe Option C exclusively. Remove all `[Zz]` references.
2. **Remove `PseudoAtomStrategy`** from §5.4. Replace with `ChiralityRecoveryUtility` implementing Option C without pseudo-atom substitution.
3. **Specify the @/@@ fix** for `OINInlineHandler.parse_inline_string()` — select one of the three fix options and document it in §5.1 and §5.2 blast radius.
4. **Add mol sanitization prerequisite** to `CIPAssigner.assign_all()` API contract in §5.4.
5. **Replace Rh-BINAP** with a genuine Zone A chiral phosphine example throughout (§1.1, §5.3, §8).
6. **Define the Spike go/no-go gate** (≥7/9 geometries = proceed; define fallback for each failure case).
7. **Confirm `scine-molassembler` PyPI name** via pip lookup before adding to `pyproject.toml`.

---

Red Team analysis complete. **Start a new conversation and run `/resolve` to begin triaging the vulnerabilities.**
