"""One place where every v0.4.5 encoder lever's default lives.

WHY THIS EXISTS
===============
By the end of v0.4.5 the encoder had nine ``OIN_*`` levers read from ``os.environ`` at nine
scattered call sites, each spelling its own default. Two spellings were already in use across
the codebase and they behave differently:

    os.environ.get("OIN_EMIT_AXIAL")                  # truthy -> "0" ENABLES it
    os.environ.get("OIN_EARLY_EXIT", "1") != "0"      # "0" disables, anything else enables

The first form is a trap: ``OIN_EMIT_AXIAL=0`` turns the lever **on**, because ``"0"`` is a
non-empty string. Anyone opting out the obvious way gets the opposite of what they asked for.

Promoting six levers to default-ON meant touching nine sites and getting the sense right at
every one. Centralizing it makes a promotion a one-line change to ``_DEFAULT_ON`` and makes
the shipped configuration readable in a single place — which matters because "which levers
are on?" is the first question anyone debugging a string difference will ask.

SEMANTICS
=========
``"0"``, ``"false"``, ``"no"``, ``"off"`` and the empty string disable; anything else enables.
So ``OIN_EMIT_AXIAL=0`` now does what it looks like. Callers may also pass an explicit
override (typically from ``ff_params``) which wins over the environment, using a **membership**
test at the call site so an explicit ``False`` can opt out — the pattern
``metallogen_adapter.py``'s ``OIN_EARLY_EXIT`` promotion established in v0.4.4.
"""

import os

#: Levers that ship ENABLED. Everything not listed here defaults to disabled.
#:
#: The six below were promoted in v0.4.5 on the evidence in
#: ``docs/PROMOTION_GATE_v0.4.5.md``: on a 300-molecule seed-42 sample, all six together took
#: byte-stability under rotation/renumbering from 58.1% to 69.6% (+35 molecules) and cut
#: comparison-key instability from 60 molecules to 16 — 1-in-5 to roughly 1-in-19 — with every
#: veto passing: fac/mer and cis/trans still distinct raw AND at key level, goldens
#: byte-identical on the default path, the mirror guard green, and ``geometry_tag_shift``
#: showing 0/298 ``[M_XXX]`` changes.
#:
#: What they have in common, and why it made them safe to promote: each one **repairs a
#: renumbered presentation without rewriting the canonical answer.** That is why the corpus
#: shows no churn. Levers that ADD information to the string (``OIN_EMIT_AXIAL``,
#: ``OIN_EMIT_LOCKED_DONOR``) are a different kind of change — the generator must then be able
#: to reproduce what they emit — and stay opt-in.
#: OIN_BORON_CAGE joined this set in v0.4.6. It is the one lever in the release with a directly
#: measured accuracy gain: on the 36 `XYZToSMILES failed` molecules of the 936-molecule
#: re-baseline, 34 are `electron-deficient boron cluster` and the lever takes them from **0/36
#: encoding to 34/36**, at 0.2-4.2s each. The boron lane separately measured 48/48 round-tripping
#: (docs/BORON_CAGE_v0.4.5.md).
#:
#: It is NOT a free win and was held back through v0.4.5 for a real reason: it moves 14 molecules
#: that were SCORED AS PASSING to failing. That is correct -- those 14 passed while describing the
#: WRONG GRAPH (VEJXOZ invents a C=B double bond) -- but it trades 14 silent false positives for
#: 14 loud honest failures, so a headline pass rate can move either way. Promoted here because a
#: notation that emits nothing for 34 molecules is worse than one that fails audibly for 14, and
#: because correctness is the point of a lossless notation.
#:
#: ⚠ A THIRD cost, measured 2026-07-26 and missing from the original pricing: the 34 molecules do
#: not GENERATE. Sample of 10 (docs/BORON_CAGE_v0.4.5.md §10, tools/boron_gen_times.jsonl): 0/10
#: produce a 3D structure, 3 fail instantly on an unsupported geometry code, and 6 burn the entire
#: generation cap producing nothing. XIQKOY_comp_0 is the two-point proof -- lever OFF it fails in
#: 0.87s with UncoordinatedFragmentError on a disconnected cage fragment; lever ON it encodes a
#: correct coordinated B10 cage and then runs past 340s. So the promotion moves this class from
#: failing INSTANTLY to failing SLOWLY: ~34 molecules x up to the 300s budget, roughly 2.8 CPU-hours
#: added per full sweep for zero additional passes. Still the right call -- the 34/34 encoder table
#: is real and a right-graph loud failure beats a wrong-graph silent pass -- but anyone reading
#: "34 now encode and round-trip" should know that "round-trip" there is NOTATION-level, and that
#: assembling a polyhedral borane cage is an open generator3d problem.
_DEFAULT_ON = frozenset(
    {
        "OIN_BORON_CAGE",
        "OIN_CANONICAL_BODY",
        "OIN_CANONICAL_PERCEPTION",
        "OIN_CANONICAL_SLOTS",
        "OIN_CANONICAL_ETA_WINDING",
        "OIN_STABLE_METAL_AC",
        "OIN_STABLE_STEREO",
    }
)

_FALSEY = frozenset({"0", "", "false", "no", "off"})

#: Deliberately NOT promoted, with the reason, so nobody has to reconstruct it.
_HELD_OFF = {
    "OIN_EMIT_AXIAL": (
        "emits a new atropisomer token the generator must reproduce; promoting converts a "
        "silent false positive into a loud false negative. Evidence to promote is recorded, "
        "and the key's _AXIAL_TOKEN_RE fold must be removed in the same commit. ALSO INTERACTS "
        "with OIN_CANONICAL_PERCEPTION, now default-ON: the Y2 cohort numbers backing that "
        "evidence (single-axis 22/22, corpus mirror audit 37/37) were measured with perception "
        "OFF. Perception feeds _is_atropisomer_candidate, which gates its steric wall on "
        "`not GetIsAromatic()` -- measured on YESKOZ, hindered axes go 2 -> 1 under canonical "
        "perception because more of the macrocycle reads aromatic. No emitted string moves "
        "today (YESKOZ's axes are non-stereogenic, so its token is empty either way; BINAP is "
        "unaffected at 1 hindered / 1 emitting), but axial.py's safety argument covers only the "
        "GENERATOR reading FEWER aromatic atoms, not the encoder reading more. Re-measure both "
        "cohorts with perception ON before promoting."
    ),
    "OIN_EMIT_METAL_CONFIG": (
        "Y1 P1 metal-centred Delta/Lambda helicity, as a trailing |mc:+|/|mc:-| sidecar. The "
        "descriptor is BUILT and validated (oin/metal_config.py, 16 tests): it detects "
        "Delta/Lambda on ZUMNEC (chiral tris-catecholato -- emits, and INVERTS under "
        "reflection) and correctly emits nothing for JEGKOW (square planar, achiral). Held "
        "opt-in for the standard information-ADDING reason: the generator must reproduce what "
        "is emitted, so promoting converts a silent collapse of Delta/Lambda enantiomers into "
        "a loud round-trip failure. Promote only with generator support plus a corpus "
        "population measurement. Note that compare.py's key does not know the token yet, so "
        "lever-ON round trips will report mismatches until it does."
    ),
    "OIN_EMIT_LOCKED_DONOR": (
        "same trade for metal-locked N/P donor configuration (P3), AND currently INCOMPATIBLE "
        "with OIN_CANONICAL_BODY, which is default-ON: canonical_body_emit reparses the body, "
        "and sanitizing the metal-free fragment clears the [N@] on a 2-degree amine -- RDKit "
        "sees a freely inverting amine, the very behaviour this descriptor works around. So "
        "with both on, the tag is stamped and then discarded. P3 therefore works only with "
        "OIN_CANONICAL_BODY=0 today (which is how tests/unit/test_locked_donor.py runs). "
        "⚠ THE OBVIOUS FIX IS WRONG -- tried in v0.4.6 and MEASURED wrong, do not repeat it. "
        "Copying the chiral tag onto the reparsed donor (the correspondence is available; "
        "_reparse_once's Guard 2 already proves same element and same heavy degree) does make P3 "
        "emit under OIN_CANONICAL_BODY, and POJJOP passes. But setting a tag AFTER the sanitize "
        "introduces a stereocentre the canonical ranker did not account for, which moves the "
        "canonical WRITE ORDER -- and @/@@ is a parity relative to that order. On RIFGUJ_comp_2 "
        "the three ring-CARBON tags then flip between a structure and its mirror, and the "
        "geometry says they must not: those carbons are pseudo-asymmetric (lowercase `s`, a "
        "RELATIVE all-cis descriptor) and read identically for the structure and its reflection. "
        "A correct fix must preserve the tag WITHOUT perturbing the ranking -- keep the donor "
        "bracketed through the sanitize, or re-derive parity from the parent geometry once the "
        "write order is fixed. Guarded by "
        "test_locked_donor.py::TestRifgujRingCarbonsArePseudoAsymmetric."
    ),
    "OIN_ACCEPT_SCORED": (
        "makes pool acceptance use the predicate the SCORE uses -- "
        "`get_oin_string(gen.mol, coords)`, i.e. `_reencode_oin_fast` -- instead of also "
        "requiring the independent `XYZToSMILES().convert` re-perception. This is the mechanism "
        "behind OIN_ETA_EARLY_EXIT's 'fires but ineffective' result: the eta early targets do "
        "fire, and then step 2 rejects the conformer anyway. MEASURED on HIDCIH_comp_1 with a "
        "forced full pool fill (tools/probe_accept_gap.py): of 48 conformers, the scored "
        "predicate matched 46 -- the FIRST at pool index 0, t=1.66s -- while independent "
        "re-perception matched only 2, first at index 25, t=49.4s. The unpatched run spends 96s "
        "reaching that conformer. So 44 conformers were scored-successes that acceptance threw "
        "away. Held OFF for two reasons that are the whole substance of the trade: (1) step 2 is "
        "the ONLY test in the predicate that does not share the generator's own connectivity, so "
        "dropping it makes acceptance circular in the same way the reported metric already is -- "
        "it buys latency, not genuine losslessness, and 2/48 independent re-perception is itself "
        "a finding worth keeping visible; (2) accepting earlier bypasses _select_by_geometry's "
        "clash-first ranking (clash.VDW_ACCEPTANCE_ENABLED is ON), so structure quality must be "
        "an arm of the promotion A/B next to pass-rate and runtime. "
        "MEASURED, 22-molecule stratified cohort, two runs (tools/ab_accept_scored.py):\n"
        "  run 1, no hard cap -- pass 18/22 BOTH arms, ZERO regressions, zero fixes; median "
        "16.01s -> 3.63s; total 1980.7s -> 626.0s; >30s 10 -> 3.\n"
        "  run 2, hard cap + honest clash metric -- median 13.87s -> 5.54s; total 1202.9s -> "
        "351.3s; >30s 8 -> 2; clash 16 over 1/17 mols (7 SEVERE, worst_overlap 0.4344) -> 2 over "
        "2/19 mols (0 severe, worst 0.7283). Its '2 fixes' (GAVSED, QIDKUL) are arm-A cap "
        "truncations, not real fixes -- both passed uncapped in run 1 at 302.5s / 390.3s.\n"
        "So all three gate arms hold, and the structure-quality concern above was REFUTED on this "
        "cohort: arm B is better, not worse (severe clashes 7 -> 0). That is consistent with "
        "_select_by_geometry's own comment that clash-minimising selection can splay donors to the "
        "edge of the gate; arm A picked a 16-clash POVPIA conformer that arm B never reached for. "
        "⚠ THE 'PROMOTE' READING OF THE ABOVE IS SUPERSEDED -- see docs/ACCEPT_SCORED_v0.4.7.md "
        "(v0.4.7 lane L2-promote), which added the three arms this cohort A/B lacked. Its G2 gate "
        "shows the pass-rate arm above is CIRCULAR: `passed` is computed with "
        "`get_oin_string(gen.mol, coords)` -- the same predicate the lever accepts on -- so "
        "'18/22 both arms, zero regressions' cannot detect what dropping step 3 costs, and is not "
        "evidence of no loss. Measured with a genuinely independent arm (full "
        "`XYZToSMILES().convert()` vs the input OIN): indep 15/20 -> 7/20, **8 regressions, 0 "
        "fixes**, one-way. And they are NOT cosmetic -- 6 of the 8 lose haptic coordination "
        "outright with the metal geometry tag degrading in lockstep, 1 reassigns the donor atom, 1 "
        "detaches a hydrogen. Its G3 gate does PASS (emitted `smiles_2` byte-identical 20/20), and "
        "G1 confirms the quality improvement above. So: byte-identical notation, CHANGED GEOMETRY. "
        "Standing recommendation from that lane is PROMOTE-WITH-SCOPE -- opt in for throughput and "
        "metric-fidelity work, leave OFF for correctness/losslessness work, do NOT flip the global "
        "default on this evidence, because the cost is invisible to the very metric that would "
        "police it. Runtime numbers above were taken on a contended box (a 5k sweep alongside), so "
        "read them as ratios within a run, not absolutes."
    ),
    "OIN_ETA_EARLY_EXIT": (
        "lets an ETA molecule short-circuit the conformer pool on the winding criterion that "
        "actually judges it. MEASURED (pool.attempts_spent): Ferrocene "
        "spends 32 attempts / 32 pool "
        "slots and never short-circuits, while non-eta CisPlatin accepts on attempt 0 -- but "
        "Ferrocene is a golden and round-trips via _select_by_geometry(honor_winding=True). So the "
        "cheap early-exit predicate (canonical_roundtrip_key) is STRICTER than the one that "
        "decides success, so eta molecules pay the full widened pool for nothing. Eta is 63.5% "
        "of the >30s runtime tail vs 19.8% of the fast set, so this is the runtime lever. It "
        "cannot accept anything the pipeline would reject: it applies the final selection's "
        "own test earlier. Promotion gate: a corpus A/B answering 'does any molecule that "
        "currently passes stop passing?', plus the attempts-spent distribution. Two fixtures "
        "is the sample size that produced four wrong answers about this tail already."
    ),
    "OIN_H_FAITHFUL": (
        "The OIN_CANONICAL_BODY interaction that blocked this in v0.4.5 is FIXED: both of "
        "canonical_body_emit's MolToSmiles writes now go through h_faithful_smiles, so the "
        "canonical body no longer discards the repair. It stays off for a DIFFERENT and better "
        "reason -- promoting it buys NOTHING measurable. A/B over the 45-molecule `Atom count "
        "mismatch` population: match 8 / mismatch 37 with the lever off, and match 8 / mismatch "
        "37 with it on. Identical. That class is not a write/read-fidelity defect; the divergence "
        "sits between the perceived parent and the emitted string (perceived_H == input_H in "
        "36/45) and is heterogeneous -- 28/45 at dH +1..+3, 4 at dH 0 where hydrogen is not the "
        "issue, and three large losses (-14, -16, -36) with no shared explanation. Two aggregate "
        "hypotheses have already been refuted (see docs/V046_HFAITHFUL_FINDINGS.md); the next "
        "step is per-ATOM provenance, not another aggregate. Promote only with evidence that it "
        "moves a real population."
    ),
    "OIN_RESCUE_STUCK_RING": (
        "its one molecule (ASISAX) encodes but is not renumbering-stable, so promoting moves "
        "it between buckets rather than fixing it."
    ),
}


def lever_enabled(name: str, override=None) -> bool:
    """Is ``name`` enabled?

    ``override`` (typically from ``ff_params``) wins over the environment when it is not
    ``None``; pass it only when the caller genuinely has an explicit setting, using a
    membership test on the params dict so an explicit ``False`` can opt out.

    Unset falls back to :data:`_DEFAULT_ON`. ``"0"``, ``"false"``, ``"no"``, ``"off"`` and the
    empty string disable — so ``OIN_FOO=0`` disables, which the older bare-truthiness reads
    got backwards.
    """
    if override is not None:
        return bool(override)
    raw = os.environ.get(name)
    if raw is None:
        return name in _DEFAULT_ON
    return raw.strip().lower() not in _FALSEY


def default_on() -> frozenset:
    """The levers that ship enabled. Useful for provenance stamping in reports."""
    return _DEFAULT_ON


def held_off() -> dict:
    """Levers deliberately left opt-in, mapped to why. Keep this honest."""
    return dict(_HELD_OFF)
