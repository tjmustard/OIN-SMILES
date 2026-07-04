# Design brief — canonicalize eta-ring fragment order + heading atom (XYZ→OIN)

Status: DESIGN BRIEF — seeds a `/hyper-architect` MiniPRD (NOT a quick patch).
The diagnostic is already done (measured below); architect can go straight to
options. Red-team focus is spelled out because this touches load-bearing
canonical-form code.

## The failing test

`tests/unit/test_stereo_roundtrip_diagnostics.py::test_haptic_face_golden_match`
(the sole remaining `@unittest.expectedFailure` in the suite). It asserts that
`Ferrocene-halide-face.xyz` → OIN(1) → generate 3D → re-encode reproduces the
pinned golden string byte-for-byte. It does not — but the mismatch is a
canonicalization relabeling, NOT a stereo/generation bug (the physical winding
is preserved; proven by the content-anchored per-ring tests in the same file,
which are hard passes).

## Measured diff (2026-07-03, Fable — real data, do not re-derive)

```
GOLDEN  : [Fe_LIN].Oc{0<}1[cH]{0}c{0}(Cl)c{0}(Br)c{0}1I . Oc{1}1[cH]{1}c{1}(I)c{1<}(Br)c{1}1Cl
REENCODE: [Fe_LIN].Oc{0}1[cH]{0<}c{0}(I)c{0}(Br)c{0}1Cl . Oc{1}1[cH]{1}c{1}(Cl)c{1}(Br)c{1<}1I
```

The complex is a ferrocene with two DIFFERENT pentahalo-substituted Cp rings
(one reads O,H,Cl,Br,I; the other O,H,I,Br,Cl). Matching rings by CONTENT (not
by position) shows exactly two relabelings, neither reversing winding:

1. **Fragment order swap.** The ring with substituent sequence `O,H,Cl,Br,I` is
   fragment 0 in the golden but fragment 1 in the re-encode (the other ring vice
   versa). Same physical structure, different ligand ordering.
2. **Heading-atom drift.** Within a ring of fixed content, the winding marker
   `<` sits on a different ring atom between the two encodings (e.g. the
   O-carbon in the golden, the I/H-carbon in the re-encode). Traversal direction
   and marker CHARACTER are identical → winding sense unchanged; only WHICH atom
   is the visible heading moved.

Both encodings are individually deterministic (XYZ→OIN on the same file is
stable; all 4 generate→re-encode runs produced the identical string). The
non-determinism is between two physically-equivalent structures (hand-built
fixture vs. generated), which the encoder's canonical form fails to collapse to
one string.

## Root causes (exact file:line — both in `src/oinsmiles/utils/oin_aligner.py`,
`OINDiscreteAligner`)

**RC1 — fragment order is arrival-order, not content-canonical.**
`:250` `base_sort_key = (i, first_binding_atom_mass, lig["smiles"])`. Leading
with `i` (the fragment's enumerate index = xyz2mol emission order) makes the
homogeneous sort essentially preserve arrival order. Two eta rings with equal
binding-atom mass (both C) tie on mass; `lig["smiles"]` would distinguish them
but is itself heading/traversal-dependent (see RC2), and `i` dominates anyway.
So the two rings order by whichever xyz2mol emitted first — not stable across
fixture-vs-generated.

**RC2 — heading atom is geometry-dependent for substituted eta rings.**
`:540-560` picks `best_idx` = the ring atom whose centroid→atom vector, rotated
into the template frame, maximally aligns with the slot's `ref_vec`. This is a
GEOMETRIC choice keyed on the ring's absolute 3D orientation, which the
generator (ETKDG placement) sets differently from the hand-built fixture. The
`SYMMETRIC_LIGANDS` forced-heading override at `:573` (which forces heading =
lowest constituent index, giving determinism for unsubstituted Cp/ferrocene)
does NOT fire for these substituted rings, so they fall through to the unstable
geometric pick.

## What to decide (the MiniPRD's design space)

**For RC1 (fragment order):** give homogeneous eta groups a CONTENT-canonical
sort identity that is independent of heading choice and arrival order — e.g. a
canonical ring signature (canonical SMILES of the ring computed independent of
the OIN heading; or the sorted substituent multiset combined with the winding
sense). Must re-order ONLY genuinely ambiguous same-mass eta fragments; must
leave every other fragment's order byte-identical.

**For RC2 (heading atom):** for substituted (asymmetric) eta rings, replace the
geometry-dependent `best_idx` with a content-canonical heading choice —
analogous to the existing `SYMMETRIC_LIGANDS` forced-heading, but content-based
(e.g. highest-CIP-priority substituent carbon, or the atom bearing the
canonically-first substituent), not "lowest local index." Whatever rule is
chosen must be computable identically from the fixture and the generated
structure (i.e. from ring TOPOLOGY/substituents, not 3D orientation).

## The correctness trap (RED-TEAM FOCUS — do not skip)

The canonicalization must be a PURE RELABELING that preserves winding. It must
NOT normalize away a genuinely reflected ring face. Concretely: after choosing a
canonical heading atom and canonical fragment order, the winding character
(`>`/`<`) must STILL be computed from geometry via `_determine_winding` /
`signed_circulation` — the fix changes WHICH atom is heading and WHICH ring is
listed first, never the winding computation. A test that would catch a
regression here already exists but is deliberately kept as a `@unittest.skip`:
`test_haptic_face_r2_geometric_fallback_never_auto_substituted` (its docstring
explains WHY it must never silently replace the exact-match assertion). The
MiniPRD must keep that safety property: a real reflected face must still produce
a different, winding-flipped string, not get canonicalized into a match.

## Constraints / non-regression

- Every existing golden must stay byte-identical: `test_regression_stability`
  (cisplatin, transplatin, cis-PtCl2(en), fac/mer-Ir(ppy)3), the ferrocene /
  ansa-metallocene haptic goldens, and all `tests/candidate_outputs/*`. The
  symmetric-eta forced-heading path and all non-eta fragments MUST be untouched
  by construction (gate any new logic on "substituted/asymmetric eta group").
- Oracle: the pinned golden `tests/candidate_outputs/Ferrocene-halide-face_oin.txt`.
  Acceptance = un-`expectedFailure` `test_haptic_face_golden_match` and it PASSES
  (generate→re-encode reproduces the golden byte-for-byte), with the per-ring
  content-anchored tests still passing and the R2 skip still meaningful.
- Do not weaken or delete the R2 geometric-fallback skip.

## Deliverable

MiniPRD via `/hyper-architect` (node `atom_oin_aligner`), then red-team →
resolve → execute (Sonnet) → audit. Acceptance test:
`test_haptic_face_golden_match` flips from xfail to a real pass; full suite green
(`discover tests/unit` skipped=3, **expected failures=0**); no golden shifts.
