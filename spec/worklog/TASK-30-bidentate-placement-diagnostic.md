# TASK-30: Bidentate/polydentate placement fidelity — diagnostic + decision

Status: DONE (2026-07-03, Sonnet) — **MIX, resolved further than "likely":**
root cause is **B** (genuine conformation mismatch — empirically confirmed,
0/72 swept angles reach a safe H-to-metal distance, and the production
algorithm's chosen angle already IS the global optimum of the full 360°
sweep, so this is NOT an objective/grid bug). BUT the required fix is
**cheap**, not the expensive "constrained re-embed" — forcing the DG
fallback (already-existing code) for DIPAMP produces a clean chelate AND
a byte-identical round trip with zero new engineering. So the actionable
fix is B's other listed option ("make DG fallback the path for
incompatible-bite polydentates"), achieved via a lightweight guard change
(element-aware/H-inclusive `:1578` check, or a tighter reuse of the
org-vs-target bite-distance delta already computed at `:1497`) — Sonnet-tier,
NOT a full HACF MiniPRD. Full measured numbers: see NOTES.md Log,
2026-07-03 "TASK-30 diagnostic" entry.
Depends on: none (independent of the stereo roadmap; unblocks the two xfail'd
round-trips `test_p_stereocenter_roundtrip` + `test_haptic_face_golden_match`)
Suggested model: Sonnet (diagnostic needs measurement care)

## Goal

The full XYZ→OIN→XYZ→OIN round-trip fails for chelating ligands. Root cause
(pinned 2026-07-03, see NOTES.md): `_stitch_fragment` places a bidentate
ligand from only its 2 binding vectors, which underdetermines the 3D rotation;
the ligand lands at a bad angle and collides with the metal. For DIPAMP, 3
ligand **H** atoms end up 1.39–1.65 Å from Rh (perceived as Rh–H hydrides by
the XYZ→OIN re-encoder → wrong topology). **No atoms are lost** — element
census is identical; this is purely a placement-geometry problem in OIN's own
template path (NOT molassembler, NOT the OIN format).

This task MEASURES the failure precisely and DECIDES the fix class — do not
commit to a fix until the diagnostic says which mechanism dominates.

## Context (exact code, no prior repo knowledge needed)

`src/oinsmiles/generation/molassembler_adapter.py`, function `_stitch_fragment`
(def ~:1319):
- `:1508` `Rotation.align_vectors(t_centered, c_centered)` — Kabsch on the
  binding atoms only. scipy warns "poorly defined" for ≤2 vectors (the bite-axis
  rotational DOF).
- `:1511-1571` — EXISTING bite-axis optimizer: sweeps 360° in 5° steps around
  the bite axis, keeping the angle that avoids heavy-atom clashes with
  `forbidden_positions` and maximizes min-distance-from-metal (origin). Note its
  min-dist objective DOES include H atoms (`np.linalg.norm(positions_aligned,
  axis=1).min()`), but its clash check (`_has_clash`) only looks at
  `forbidden_positions` (already-placed OTHER fragments), not the metal.
- `:1573-1584` — rejection guard: if any non-binding **heavy** atom is < 1.7 Å
  from the metal, return `None` → DG fallback. **`symbols[i] != "H"` means H
  collisions never trigger this guard** — the likely reason DIPAMP is accepted
  despite H atoms at 1.4 Å.
- Isolated-ligand conformation: `_stitch_fragment` embeds each ligand in
  isolation (ETKDGv3, fixed seed) BEFORE Kabsch. If the isolated conformation's
  bite distance/backbone dihedral doesn't match the chelate's cis bite, NO rigid
  rotation can fix the collision — that would be a conformation problem, not a
  rotation-DOF problem.

Fixtures that exercise this: `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz` (bidentate
P^P), and per the code comment `:1575` the ppy chelates in fac/mer-Ir(ppy)3.

## Diagnostic steps (measure, don't fix)

Write a throwaway script (put it in the scratchpad, not the repo) that, for
DIPAMP and one Ir(ppy)3 fixture:
1. Runs `XYZToSMILES().convert()` → OIN, then `OIN3DGenerator().generate()`,
   and reports the generated structure's per-element nearest-neighbour
   distances to the metal (confirm/deny the H-collision; quantify).
2. Instruments `_stitch_fragment` (temporary prints or a copy) to report, for
   the chelate fragment: did the bite-axis optimizer run? what `best_angle` did
   it pick? what was the min-H-to-metal distance at that angle vs the best
   achievable over the full 360° sweep? — i.e. **is a good angle available but
   not chosen (objective/grid problem), or is every angle bad (conformation
   problem)?** This is THE decision question.
3. Check what the DG fallback (`_molassembler_worker`) produces for DIPAMP if
   the guard is made to fire (temporarily lower the guard / force fallback):
   does molassembler give a clean chelate, or does it also fail? This tells us
   whether "tighten the guard so it falls back to DG" is even a viable fix.

## Decision → fix class (record in NOTES.md + this file)

- **A — objective/guard bug (cheap, likely partial):** a good bite-axis angle
  exists but isn't chosen, and/or the guard misses H collisions. Fix = include
  H (element-aware radii) in the `:1578` rejection guard, and/or add a
  metal-clash term to the bite-axis objective, and/or finer grid. Lightweight
  task, Haiku/Sonnet. Only sufficient if a good angle actually exists.
- **B — conformation mismatch (deeper):** every rigid angle collides → the
  isolated ligand conformation can't span the bite. Fix = constrained re-embed
  of the chelate to the target bite distance, or make the DG fallback the
  primary path for polydentates. This is a design change → HACF MiniPRD
  (`/hyper-architect`), with DIPAMP round-trip as the acceptance test.
- Likely a mix (guard fix as a safety net + B for real fidelity). The
  diagnostic numbers decide the split.

## Acceptance (of THIS task — measurement only)

- A written decision (A vs B vs mix) in `spec/worklog/NOTES.md` Log + this
  file's Status, backed by the measured numbers (best-achievable vs chosen
  H-to-metal distance; DG-fallback quality).
- NO `src/` change in this task. NO new committed test yet (the fix task adds
  the round-trip assertion; the two existing xfails already track the symptom).

## Constraints / DO NOT

- Do NOT modify `src/` — measurement only. Do NOT flip the existing xfails.
- Do NOT assume the fix is A without the step-2 evidence that a good angle
  exists — that is the whole point of the diagnostic.

## On completion

Set `Status:`, append a dated Log entry to `spec/worklog/NOTES.md` with the
measured numbers and the A/B/mix decision, and (if B) note that the follow-on
is a HACF MiniPRD, not a lightweight task.
