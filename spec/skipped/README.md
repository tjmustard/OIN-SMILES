# spec/skipped/

Specs that were fully compiled (Draft → Red Team → Resolve) and marked
"Ready for Implementation," but were never executed and are not currently
planned to be. Distinct from `spec/archive/`, which holds specs for work
that **was** implemented and audited.

A file lands here when either:
- It was superseded by a later spec before any implementation commit
  followed it, or
- Its blocking dependency was itself deferred/archived and the spec was
  never updated to reflect that, leaving it orphaned.

Do not treat anything in this folder as pending or "next up." If this work
is picked up again, it needs a fresh `/hyper-architect` pass against the
current codebase, not a resurrection of these files as-is.

## Contents (moved 2026-07-03)

The "Direct Parser v0.2.2" remediation chain — compiled 2026-05-07 to fix
5 blockers identified in `spec/audit/DirectParser_IntegrationAudit_20260506.md`,
but never implemented. `parse_oin_direct()` (`src/oinsmiles/generation/engine.py`)
still exists but is uncalled in production; `OIN3DGenerator.generate()` runs
only the legacy `OINParser.parse()` + `MolassemblerAdapter.generate()` path.
Only 1 of 5 blockers (FragmentMapping) was ever executed and audited — see
`spec/archive/MiniPRD_DirectParser_FragmentMapping_v0.2.2_AUDITED.md`.

- `SuperPRD_DirectParser.md` (v0.2.1, design-only) — superseded by the v0.2.2
  SuperPRD below (see its `Supersedes:` field)
- `SuperPRD_DirectParser_v0.2.2.md` — the 5-blocker remediation plan
- `MiniPRD_DirectParser_Polydentate_v0.2.2.md`
- `MiniPRD_DirectParser_EtaBonds_v0.2.2.md`
- `MiniPRD_DirectParser_Tests_v0.2.2.md`
- `MiniPRD_DirectParser_Permutation_v0.2.2.md`
- `MiniPRD_DirectParser_Verification.md` — orphaned: still labeled v0.2.1
  and blocked by the already-deferred/archived Integration MiniPRD; never
  updated to the v0.2.2 chain
