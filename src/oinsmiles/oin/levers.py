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
#: ``docs/agentic-notes/v0.4.5/PROMOTION_GATE_v0.4.5.md``: on a 300-molecule seed-42
#: sample, all six together took
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
#: (docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md).
#:
#: It is NOT a free win and was held back through v0.4.5 for a real reason: it moves 14 molecules
#: that were SCORED AS PASSING to failing. That is correct -- those 14 passed while describing the
#: WRONG GRAPH (VEJXOZ invents a C=B double bond) -- but it trades 14 silent false positives for
#: 14 loud honest failures, so a headline pass rate can move either way. Promoted here because a
#: notation that emits nothing for 34 molecules is worse than one that fails audibly for 14, and
#: because correctness is the point of a lossless notation.
#:
#: ⚠ A THIRD cost, measured 2026-07-26 and missing from the original pricing: these molecules
#: mostly do not GENERATE. 33 of 34 measured (docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md §10,
#: tools/boron_gen_times.jsonl): only **2/33** produce a 3D structure, 25 burn the entire generation
#: cap producing nothing (~2.1 CPU-h wasted per 5k sweep), 6 fail instantly. An earlier 10-molecule
#: sample read 0/10 and was used to propose a blanket boron fast-fail; the full measurement REFUTES
#: that -- a blanket fast-fail would cost the 2 that work, and no clean discriminator was found.
#: XIQKOY_comp_0 is the two-point proof -- lever OFF it fails in
#: 0.87s with UncoordinatedFragmentError on a disconnected cage fragment; lever ON it encodes a
#: correct coordinated B10 cage and then runs past 340s. So the promotion moves MOST of this class
#: from failing INSTANTLY to failing SLOWLY: 25 molecules x up to the 300s budget, ~2.1 CPU-h
#: per full sweep, for the 2 passes it does buy. Still the right call -- the 34/34 encoder table
#: is real and a right-graph loud failure beats a wrong-graph silent pass -- but anyone reading
#: "34 now encode and round-trip" should know that "round-trip" there is NOTATION-level, and that
#: assembling a polyhedral borane cage is an open generator3d problem.
#:
#: OIN_INDEP_SCORE joined this set in v0.4.8, and it is the odd one out: it changes no encoder
#: or generator output at all. It changes what the HARNESS reports about them. The harness had
#: been scoring a round trip with ``get_oin_string(gen_result.mol, coords)`` -- the GENERATOR's
#: own bond graph -- which is not merely inaccurate but CIRCULAR: ``gen_result.mol`` is exactly
#: the artifact that would have to be wrong for the test to fail. Measured consequences, both
#: directions, 936-molecule cohort (docs/agentic-notes/v0.4.6/METRIC_FALSE_POSITIVES.md): 61
#: false positives (9.6%; 28.1% of HAPTIC inputs) where the graph asserts bonds the coordinates
#: do not support, and 8 false negatives where it drops stereo they do.
#:
#: Both directions close with ONE call -- a full ``XYZToSMILES().convert()`` of the generated
#: XYZ, re-perceiving bonds AND stereo from coordinates alone.
#:
#: ⚠ THE COST ARGUMENT THAT HELD THIS OFF WAS NEVER MEASURED, AND IS WRONG. The old rationale
#: priced the second encode at "0.4-1.5 s/molecule, so a 5k sweep pays 1-2 CPU-hours" and cited
#: that as a reason to wait. Measured over the full corpus (tools/honest_rescore.py, 4688 stored
#: structures): **0.33 s/molecule median**, and the whole corpus re-scores in ~5 minutes. Against
#: a sweep whose own median is 7.19 s and whose p95 is the 300 s timeout -- the ~55 CPU-hour
#: range -- the honest metric costs a low single-digit percentage of the run it corrects. There
#: was never a cost case for scoring dishonestly.
#:
#: SCOPE, deliberately narrow: this promotion changes what is REPORTED, not what is ACCEPTED.
#: The harness's tier ladder and the generator's ``accept_fn`` both keep their existing
#: predicates. Moving those changes runtime and the failure mix, and doing it in the same
#: release that re-baselines the number would make both unmeasurable -- the identical confound
#: that let the ``OIN_ACCEPT_SCORED`` A/B report "zero regressions" while the honest arm read 8.
#:
#: OIN_CANONICAL_DONOR_FOLD and OIN_FOLD_PARITY_VETO joined this set in v0.4.13, together and
#: never separately. They are ONE promotion with two names: the fold is the change that pays,
#: the veto is the condition under which it is safe, and the pair is the smallest unit that is
#: both. Read the coupling invariant below before touching either.
#:
#: WHAT IT BUYS. v0.4.11 built the within-fragment donor fold, measured +7.86 ``byte_exact``
#: points across 393 molecules -- one direction, 0 molecules moving any other way, comparison
#: key unchanged on 992 strings, both gate arms byte-identical -- and then REFUTED it: the fold
#: collapses enantiomers in 221 of those same 393 gains (56.2%), and 19 of a uniform 250-molecule
#: draw (7.6%). v0.4.12 built the reflection-parity veto that separates the two sets. Through the
#: shipped predicate, 171 of the 393 gains survive: **+3.42 points**, 72.46% -> 75.88%.
#:
#: WHY THE NUMBER IS TRUSTED. It was reached twice, by different mechanisms, and the two agree to
#: a single molecule. v0.4.11 bounded the safe set from ABOVE by counting collapses in a mirror
#: audit ("at most ~172 of 393"); v0.4.12 counted SURVIVORS from below through the veto itself
#: (171). An upper bound from damage and a lower bound from survival meeting at 171/172 is
#: stronger evidence than either alone, because a bug that inflated one would have to inflate the
#: other in the same direction by the same amount.
#:
#: THE GATE THAT PROMOTED IT. Mirror audit ``tools/mirror_audit_donor_fold.py`` on a uniform draw
#: reads 19 -> 0 at seed 7 AND at seed 11, with ``achiral_or_preexisting_fold`` unmoved (157 and
#: 141) -- so the veto is not buying its zero by declining everything. The seed-11 BASELINE
#: independently reproduces 19/250 on a disjoint draw, which is what establishes 7.6% as a
#: property of the corpus rather than of seed 7.
#:
#: 🔴 COUPLING INVARIANT -- ``OIN_CANONICAL_DONOR_FOLD`` IS ONLY SAFE WITH ``OIN_FOLD_PARITY_VETO``.
#: The fold ON with the veto OFF is the exact configuration v0.4.11 refuted: it silently collapses
#: enantiomers, and NEITHER ``byte_exact`` NOR the round-trip comparison key can see it, because
#: ``compare._parse_vertex_colors`` folds that axis deliberately. ``test_levers`` pins the pair
#: in this set together; if you demote one, demote both, and if a future release needs them
#: apart, the mirror audit is the arbiter, not the gate arms.
#:
#: SCOPE. There is exactly one production path into the fold: ``perception_tmc.get_oin_string``
#: calls ``fold_parity.resolve``, and ``canonicalize_oin_slots`` has no other caller in ``src/``.
#: The veto declines to fold -- returning the rotation-only labeling, the pre-v0.4.11 answer --
#: whenever it cannot build its own evidence (no conformer, no atom/coordinate pairing, a failed
#: reconstruction, a failed mirror encode). Every failure mode of the instrument therefore lands
#: on the SAFE side, which is why a corpus reading of ``declined_* = 0`` is what proves it alive
#: rather than the separation count.
#:
#: ⚠ WHAT THIS PROMOTION COSTS, stated because it was accepted rather than avoided. It changes
#: the default answer for 171 molecules, which VOIDS the carry-forward licence: v0.4.13 owes a
#: full re-sweep and re-frozen ARM 1 / ARM 2 goldens, and no table measured before it may be
#: quoted beside the new headline. It also lands in the same release as the harness
#: false-negative measurement (Lane 1/Lane 2), so a reader cannot attribute a headline move to
#: the fold rather than to the harness without reading both apart. That confound is real, was
#: raised, and was accepted by the project owner; the release doc states the two separately and
#: the sweep that produces the headline runs with the Lane levers OFF.
#: OIN_RESONANCE_DONOR_FOLD joined this set in v0.4.14. It widens OIN_CANONICAL_DONOR_FOLD's
#: donor-equivalence test to the fragment's constitutional skeleton, so acac written ketone/enol
#: and a sulfonate's oxygens stop being read as two inequivalent donors. It is therefore subject
#: to the SAME coupling invariant as the fold it widens: only safe with OIN_FOLD_PARITY_VETO on,
#: pinned by test_resonance_donor_fold::TestResonanceFoldInheritsTheVetoCoupling.
#:
#: WHAT IT BUYS: **+1.42 points -- 78 gains, 7 losses, net +71 molecules**, byte_exact
#: 75.88% -> 77.30%. MEASURED end-to-end over all 182 molecules the lever can affect.
#:
#: 🔴 THE FIRST NUMBER THIS LEVER SHIPPED WITH (+1.56, bad_direction = 0) WAS WRONG. It came from an
#: offline re-score licensed by fold_key_invariance.py reading 0 keys changed. That verdict bounds
#: the ACCEPTANCE step -- accept_fn decides by key -- but the generator's input is the OIN STRING,
#: so a slot relabeling changes ParsedOIN, the CoordMap, and the POOL ITSELF. An offline re-score
#: holds the generated structure fixed, so for a molecule that already round-trips it can only ever
#: print "still fine": it reports bad_direction = 0 whether or not losses exist.
#: Measured honestly end-to-end (tools/generator_ab_honest.py, seeded generator), over the COMPLETE
#: affected population rather than a sample -- tools/lever_string_movers.py derives it from
#: coordinates (93 molecules move encode(input); union with the generated-side movers is 182), and
#: every other molecule is unchanged BY CONSTRUCTION because it receives an identical input string:
#:     n = 182 of 182   78 REAL GAINS   7 REAL LOSSES   net +71
#: So no sweep was needed: the affected set is derivable, and an A/B over it is exact.
#: The offline sim had the GAINS exactly right (78) and was blind only to the losses -- which is its
#: defining limitation, not bad luck: it re-encodes a FIXED structure, so for a molecule that
#: already round-trips it can only ever print "still fine".
#: ⚠ The 7 losses are NOT the encoder emitting a wrong string. The encoder is canonical in both
#: arms; the GENERATOR, handed an equally valid but differently-LABELLED input, builds a worse
#: structure. The generator's output depends on the slot labeling of its input and it should not.
#: Every future canonicalization lever pays this, which is why it is a v0.4.15 lane.
#: See docs/agentic-notes/v0.4.14/GENERATOR_NEUTRALITY_HAS_A_HOLE_v0.4.14.md.
#:
#: WHY THE GAINS ARE THE RIGHT ONES. All 78 land inside the 103-molecule class the v0.4.13 fold
#: provably cannot reach; 0 come from the 222 the parity veto reverts, 0 from rdkit_canonical,
#: 0 from anywhere else. The lever hits its target population and nothing adjacent.
#:
#: 🔴 THE GATE, AND ITS COVERAGE -- which is the part v0.4.13 got wrong and this states first.
#: A uniform cat/ draw (n=250, seed 7) contains **1 of the 179 movers = 0.4% coverage**, so its
#: identical before/after tally is NOT safety evidence; it says only that the lever does not
#: damage molecules it never touches. (It is still run, and it reproduces v0.4.13's frozen
#: 157/92/1 exactly, which is what proves the instrument and the draw are alive.)
#: The gate that can SEE this change is a mover-enriched cohort at **179/179 = 100% coverage**:
#:     distinct_both_arms 108 · achiral_or_preexisting_fold 71 · REGRESSIONS 0, both arms
#:     per-molecule verdicts differing between arms: 0
#: `achiral_or_preexisting_fold` unmoved at 71 is what rules out buying the zero by declining
#: everything. And the number that separates "safe" from "never fired on anything chiral":
#: **33 of the 78 gains are on molecules the encoder resolves as chiral (distinct_both_arms) and
#: STILL resolves after the widening.** The other 45 are achiral, where folding cannot collapse
#: anything.
#:
#: ⚠ WHAT IT COSTS. It changes the default answer for 78 molecules, so arm2's goldens move:
#: 7 of 325 rows re-frozen in gate_v049_arm2_golden.tsv, 1 of 100 in gate_v047_arm2_golden.tsv.
#: (Coverage of the moved population is the larger, different count -- 12/325 and 2/100 -- because
#: the v0.4.9 golden predates the donor fold and some of its labelings coincide with today's.)
#: ARM 1 is unaffected (0 of 62 fixtures are movers) and stays byte-identical, #DONE 62.
#: ⚠ arm2's field 3 is a FRESH STOCHASTIC GENERATION, so a changed out-hash is not attributable to
#: an encoder lever. Two rows flipped round-trip status in the re-freeze and NEITHER is this
#: lever's doing -- the control on the frozen sweep structures reads HEKFEL True->True and FOJJUM
#: False->False. Recorded as comments inside the goldens so a future diff does not misread them.
_DEFAULT_ON = frozenset(
    {
        "OIN_BORON_CAGE",
        "OIN_CANONICAL_BODY",
        "OIN_CANONICAL_DONOR_FOLD",
        "OIN_CANONICAL_PERCEPTION",
        "OIN_CANONICAL_SLOTS",
        "OIN_CANONICAL_ETA_WINDING",
        "OIN_FOLD_PARITY_VETO",
        "OIN_INDEP_SCORE",
        "OIN_RESONANCE_DONOR_FOLD",
        "OIN_STABLE_METAL_AC",
        "OIN_STABLE_STEREO",
    }
)

_FALSEY = frozenset({"0", "", "false", "no", "off"})

#: Deliberately NOT promoted, with the reason, so nobody has to reconstruct it.
_HELD_OFF = {
    "OIN_PREFILTER_ADVISORY": (
        "v0.4.13 Lane 1. `_reencode_key_matches` is two-stage: step 1 re-serializes the generated "
        "geometry through the GENERATOR'S OWN contract-mol connectivity and rejects on a key "
        "mismatch, on the stated theory that 'a MISMATCH here is a reliable geometry-is-wrong "
        "signal'. This lever makes step 1 ADVISORY -- a cheap mismatch falls through to the "
        "strict independent test instead of returning False.\n"
        "THE THEORY HAS A COUNTER-EXAMPLE. On AROHIA_comp_0 the cheap test matches 0/48 "
        "conformers while the strict independent test matches 16/48. Because the cheap veto "
        "returns False first, those 16 are unreachable IN BOTH ARMS OF EVERY A/B EVER RUN on "
        "this. The defect runs in the PESSIMISTIC direction -- it makes the project look worse "
        "than it is -- which is why nobody chased it.\n"
        "SCOPE: acceptance only. The SCORING half of this defect was closed by OIN_INDEP_SCORE "
        "in v0.4.8 -- on the frozen corpus, cheap-fails-but-independent-passes is 28/5000 and "
        "the honest metric already counts every one correctly. Do not re-open it.\n"
        "🔴 v0.4.14 MEASURED BOTH HALVES OF THAT GATE. The prevalence is LOW and the latency "
        "objection is REFUTED -- so this lever is still held off, but NOT for the reason that "
        "was on record here.\n"
        "PREVALENCE (49 of 50 stratified-sample molecules measured, seed 7, 1 hard timeout; 5 "
        "INSTRUMENT_DEAD excluded so the live denominator is 44): the cheap prefilter vetoed "
        "261 conformers and the strict test disagreed with it on 4. That is 4/261 = 1.5% of "
        "vetoes and 2/44 = 4.5% of molecules. AROHIA_comp_0's 0/48-vs-16/48 is an OUTLIER, not "
        "the typical case -- which is precisely what n=1 could not tell anyone.\n"
        "LATENCY, the claim this entry used to make: the lever is measured FASTER, not slower. "
        "Total -478.8s (-11.6%) over 49 molecules, and it decomposes: the 2 molecules where it "
        "ACTED contribute -498.0s (HUTCOQ_comp_0 alone is -495.0s, an early exit that stopped a "
        "long pool fill), while the 42 where the predicate ran and it did NOTHING contribute "
        "+19.9s total -- about +0.47s/molecule. Anyone citing cost as the reason to keep this "
        "off is citing a refuted number.\n"
        "STILL HELD OFF because the ACCURACY case is not made: the single fail->pass recovery in "
        "the sample is a stochastic A/B on a stochastic generator and cannot carry a corpus "
        "projection. What would settle it is a NON-A/B recovery count -- force the pool full "
        "(tools/probe_accept_gap.py) on the 26 molecules where the prefilter actually vetoes and "
        "count conformers the strict test would have accepted. Bounded work; no sweep needed.\n"
        "⚠ A lever that never fires and a lever that fires and finds nothing BOTH report zero "
        "overrides. Gate on adapter.prefilter_veto_overridden being non-zero on AROHIA_comp_0 "
        "(cheap 0/48, strict 16/48) before quoting any corpus number."
    ),
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
        "⚠ THE 'PROMOTE' READING OF THE ABOVE IS SUPERSEDED -- see "
        "docs/agentic-notes/v0.4.7/ACCEPT_SCORED_v0.4.7.md "
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
    "OIN_ATTACH_CHECK": (
        "OIN_ACCEPT_SCORED's missing safety condition -- 'accept the first conformer the score "
        "credits THAT STILL HAS ITS LIGANDS ATTACHED'. Coordinate-only: every coordination SITE "
        "the generator claims must retain at least one atom inside bonding distance, measured on "
        "_get_tmc_mol_impl's own donor set (xyz2AC_obabel tolerance 0.5, THEN the "
        "aromatic-ring-carbon filter). It NEVER reads a bond object as evidence of attachment -- "
        "a detached ligand keeps its bond, so any GetBonds()-based check certifies exactly what "
        "it exists to catch. Falsified before it was built "
        "(docs/agentic-notes/v0.4.7/ATTACH_CHECK_v0.4.7.md §1) on 21 "
        "molecules / 40 accepted conformers in BOTH arms: separates 7 of the 8 known indep "
        "regressions with 0 false positives over 22 round-tripping conformers. The two predicates "
        "the promote lane proposed BOTH fail -- count-based wrongly rejects 11 of 22, set-based 3 "
        "-- because the raw donor-ATOM count is not conserved against the encoder's perception "
        "even on good structures (MEDZUR claims 10, the coordinate path perceives 7, and it "
        "round-trips). Known residual, ship it honestly: POVPIA is NOT caught (metal sphere "
        "intact; a hydrogen detaches and C-N reads as C=N), so this is 7/8, never 8/8. Cost 7-81 "
        "ms/conformer against the 48-57 s strict test it replaces. GATES RE-RUN with it on "
        "(docs/agentic-notes/v0.4.7/ATTACH_CHECK_v0.4.7.md §3): on the n=100 guard population it "
        "ERASES G1's failure "
        "(clash 77->82 and severe 5->5, against the bare lever's 77->105 and 5->14) and recovers "
        "17 of the 26 G2 regressions, with G3 still 97/97 byte-identical. It does NOT clear the "
        "promotion bar: G2 still regresses 9 molecules (population rate 29.5% -> 11.4%, not 0), "
        "only 43% of the speedup survives (3.33x -> 1.43x vs default), >30s goes 1 -> 14, and "
        "DURPAH_comp_0 stops passing by exceeding the budget. Hence default-OFF. ⚠ BUT the "
        "PAIRING is a separate and clearer call: against the bare lever this check is strictly "
        "one-directional -- 17 indep fixes, ZERO regressions, clash 106->83, severe 14->5 -- so "
        "anyone enabling OIN_ACCEPT_SCORED should enable this too. Residual is three classes, "
        "only one of them predicted: POVPIA (ligand-internal, unreachable by design), MEDZUR "
        "(attachment intact and indep still disagrees -- class size UNMEASURED), and GAVSED "
        "(acceptance rejected everything, and _select_by_geometry's fallback ranking never "
        "consults this check -- the check guards ACCEPTANCE, not RETURN)."
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
        "is the sample size that produced four wrong answers about this tail already.\n"
        "🔴 v0.4.12: THAT PROMOTION GATE IS VOID, NOT UNRUN, AND THIS LEVER IS RUNTIME-INERT "
        "AS SITED. It lives in _select_by_geometry_impl, which runs AFTER "
        "generate_3d_structures has already filled the entire pool, so it cannot reduce embed "
        "count by construction -- its own A/B says exactly that (Ferrocene: fires, attempts "
        "32 -> 32). Spending a corpus A/B to re-measure a documented structural null would "
        "have bought nothing. The criterion was moved to the site that CAN stop pool filling; "
        "see OIN_ETA_ACCEPT_EXIT. Keep this one only as the marker of where the boundary is."
    ),
    "OIN_ETA_ACCEPT_EXIT": (
        "v0.4.12 Lane 2: OIN_ETA_EARLY_EXIT's winding criterion, relocated into accept_fn -- "
        "the only site consulted per conformer DURING pool filling, and therefore the only one "
        "that can shorten it. Eta is 53.1% of the whole >30s tail (528 of 994), so this is the "
        "runtime lever the ladder has been pointing at since v0.4.6.\n"
        "IT IS A CONJUNCTION, DELIBERATELY. Accepting on winding alone would stop the pool "
        "before _select_by_geometry's clash-first ranking ran -- structurally the SAME defect "
        "as OIN_ACCEPT_SCORED, which cost 26 independent-re-perception regressions on n=100. "
        "So a conformer must ALSO classify as the requested geometry AND pass "
        "conformer_ligands_attached. That attachment check is unconditional inside this branch "
        "rather than gated on OIN_ATTACH_CHECK, because v0.4.7's one unambiguous finding was "
        "'never run a scored-acceptance lever without it' -- here it is part of what the "
        "predicate MEANS, not an option.\n"
        "Promotion gate: the four gates re-measured from scratch (the v0.4.7 stored arms no "
        "longer exist, and they predate both the honest metric and v0.4.10's memoization), "
        "PLUS a fifth arm. G1-G4 are structurally blind to metal configuration: G3 compares "
        "smiles_2 bytes and compare.py folds |mc:|, so an arm could return the opposite "
        "enantiomer and every gate would report identical -- the exact shape of the v0.4.11 "
        "refutation. Report token_for_mol per accepted conformer per arm; any divergence "
        "blocks promotion."
    ),
    "OIN_ATTACH_RETURN": (
        "v0.4.15 Lane 1: v0.4.7's attachment predicate on the RETURN path, where it has never "
        "run. OIN_ATTACH_CHECK's own entry above names this gap as its third residual class -- "
        "'GAVSED (acceptance rejected everything, and _select_by_geometry's fallback ranking "
        "never consults this check -- the check guards ACCEPTANCE, not RETURN)'. Deliberately a "
        "SEPARATE lever: OIN_ATTACH_CHECK governs acceptance and carries its own promotion "
        "history and measured costs, and conflating the two would make both unmeasurable and a "
        "future demotion impossible to scope.\n"
        "WHY IT IS A PREFERENCE, NOT A REJECTION. `_attach_rank` is a sort key that demotes a "
        "detached conformer below an attached one at all three exits of "
        "_select_by_geometry_impl (the geometry-classified sort, the eta-winding scan, and the "
        "lowest-energy fallback). When NO conformer in the pool holds its sites, the pre-lever "
        "answer -- lowest energy -- is still returned. So the lever is non-regressive in "
        "accuracy by construction and its whole cost is latency; converting the no-attached "
        "case into a loud failure is OIN_ATTACH_RETURN_STRICT's job, separately.\n"
        "MEASURED BEFORE IT WAS BUILT (tools/attach_return_preflight.py on the v0.4.14 baseline "
        "sweep, measurements/v0.4.15/attach_return_preflight.json). The guard's own predicate -- "
        "a claimed coordination SITE holding no atom, coordinate-only, OIN slot count as the "
        "claim -- fires on structural/DETACHED 289/301 = 96.0% (5.78 pts) against "
        "byte_exact/DETACHED 1/52 = 1.9% and byte_exact/INTACT 0/250. 50x separation, and the "
        "exposure is 52x SMALLER than the charter's upper bound, because coordination_report "
        "('did the donor set change', input vs generated) and ligands_attached ('did a site go "
        "empty', generated alone) are different tests: the 52 lose 1-3 light/ambiguous donors "
        "(H, Si, B, F) with the actual count often unchanged or higher, while the 301 lose whole "
        "multi-carbon haptic groups (13x 'Ni: lost {C:6} (7->1)').\n"
        "COORDINATE-ONLY, ALWAYS. A detached ligand keeps its bond, so any GetBonds()-based "
        "attachment test certifies exactly what this exists to catch. Known residual, ship it "
        "honestly: POVPIA is NOT caught (metal sphere intact; a hydrogen detaches and C-N reads "
        "as C=N) -- 7/8, never 8/8. Cost 7-81 ms/conformer, now also on the return path; the eta "
        "winding scan KEEPS its early exit (a full scan there would pay _reencode_oin at 48-57 s "
        "per haptic conformer, more than the lane can win) and steps over a detached match "
        "while remembering it."
    ),
    "OIN_ATTACH_RETURN_STRICT": (
        "v0.4.15 Lane 1, second stage: when no conformer in the pool holds every coordination "
        "site the OIN claims, raise instead of returning the lowest-energy detached one. "
        "Requires OIN_ATTACH_RETURN (it only reaches the branch that lever creates).\n"
        "SEPARATE FROM OIN_ATTACH_RETURN ON PURPOSE, because the two have opposite risk "
        "profiles. The preference cannot lower accuracy; this CAN, and is meant to: it converts "
        "a silently-wrong structure into an honest hard_fail, which is the trade v0.4.6's boron "
        "promotion and v0.4.8's honest metric both took deliberately. Expect "
        "structural -> hard_fail movement and a LOWER byte_exact. Report the transition matrix, "
        "never a net.\n"
        "It RAISES rather than returning None: both _select_by_geometry call sites do "
        "get_xyz_string(chosen_mol) on the next line, so a None would surface as an "
        "AttributeError from inside the xyz writer -- a real failure reported as a mystery. The "
        "ValueError matches the idiom the OIN-direct path already uses for 'no conformers'."
    ),
    "OIN_ACCEPT_STRING_EXACT": (
        "v0.4.15 Lane 2: acceptance compares the OIN STRING, not only the round-trip key.\n"
        "🔴 THE DEFECT. `compare._parse_vertex_colors` folds reflection DELIBERATELY -- it "
        "colours every donor atom of a ligand with that ligand's WHOLE canonical body, so two "
        "same-coloured vertices are interchangeable and a transposition between them (an ODD "
        "permutation, i.e. a reflection) is invisible to `_polyhedron_signature`. And `accept_fn` "
        "decides by that key. So the generator can build the MIRROR IMAGE, acceptance takes it, "
        "and the harness files the result as a benign same-isomer string difference -- in a "
        "bucket whose name says the difference is benign. Third instance of one rule in this "
        "project, after v0.4.8 (scored vs honest) and v0.4.11 (the donor fold): A LOSSY KEY MUST "
        "NEVER BE REUSED AS AN ACCEPTANCE PREDICATE FOR AN AXIS IT FOLDS.\n"
        "MEASURED, and this is what picked the predicate. Over the 183 known MIRROR_MATCH "
        "molecules, re-read from the v0.4.14 baseline sweep's own smiles_1/smiles_2_indep: "
        "normalized strings DIFFER 183/183, keys AGREE 183/183. The mechanism is visible in the "
        "strings -- AFADOC_comp_0 `O{5}`/`O{4}` swapped, AGAVIQ_comp_0 `P{3}`/`P{4}` swapped. On "
        "the baseline sweep the class re-derives to 201 of 242 (83.1%).\n"
        "WHY NOT metal_config, which the charter proposed. An acceptance test needs a REFERENCE "
        "handedness, and the generator's only input is the OIN string. "
        "`_select_by_geometry_impl` already has a helicity-aware branch and it is dead by "
        "construction: `parse_metal_config_token(parsed.original_oin)` is None for the whole "
        "corpus because OIN_EMIT_METAL_CONFIG is off, and emitting the token is a strictly larger "
        "change (v0.4.16). The string needs no new descriptor and no new emission.\n"
        "WHY NORMALIZED, NOT RAW. Raw carries the metal `@OH`/`@SP` labels the encoder documents "
        "as atom-order-dependent and irreproducible, so raw equality would reject nearly "
        "everything. `normalize_oin_for_comparison` strips exactly those and keeps slots and "
        "winding: strictly between the key (too loose) and raw (too strict).\n"
        "NON-REGRESSIVE BY CONSTRUCTION, via ACCEPT_INCUMBENT. An accepted conformer is returned "
        "as the SOLE pool member (`return [early_hit]`), so plainly REJECTING a key-equal "
        "conformer would let the energy-sorted pool hand back a DIFFERENT mol -- a real "
        "regression. Instead the predicate returns generator3d.ACCEPT_INCUMBENT: the pool keeps "
        "filling, and if nothing string-exact turns up the FIRST incumbent is returned, which is "
        "byte-identically what the pre-lever run returned. The cost is therefore latency, not "
        "accuracy -- expect `>30 s` UP, because this fires on all 365 key_equal molecules and "
        "every one of them currently stops the pool early.\n"
        "SCOPE IS DELIBERATELY WIDER THAN THE ENANTIOMERS (owner-accepted 2026-07-29). It fires "
        "on `slot_renumber` too, which also fails byte_exact and is therefore the same defect; "
        "ceiling 365 molecules / 7.30 pts rather than the enantiomer class's 4.02.\n"
        "⚠ THE COMPARISON KEY CANNOT SCORE THIS LEVER. It folds reflection, so an A/B scored with "
        "it reports 'no change' BY CONSTRUCTION -- a hard blocker, not a caution. `byte_exact` is "
        "reflection-blind for the same reason, so a molecule that stops being MIRROR_MATCH may "
        "not change bucket at all: report the MIRROR_MATCH count "
        "(tools/veto_residue_chirality.py) as the primary number and byte_exact as secondary. "
        "Delta/Lambda is also not the only chirality axis -- atropisomers (OIN_EMIT_AXIAL, Y2) "
        "and bound secondary amines (P3) are NOT covered, so report a count rather than claiming "
        "the class goes to zero."
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
        "hypotheses have already been refuted (see "
        "docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md); the next "
        "step is per-ATOM provenance, not another aggregate. Promote only with evidence that it "
        "moves a real population."
    ),
    "OIN_RESCUE_STUCK_RING": (
        "its one molecule (ASISAX) encodes but is not renumbering-stable, so promoting moves "
        "it between buckets rather than fixing it."
    ),
    "OIN_BORON_GEN_FASTFAIL": (
        "detects a COORDINATED deltahedral boron-cage ligand fragment (the same B-B-B triangle "
        "motif OIN_BORON_CAGE gates on) on a geometry outside the safe set, BEFORE 3D generation "
        "starts, and raises immediately instead of running the embed attempt loop. "
        "⚠ READ THIS FIRST: THE BLANKET FAST-FAIL THIS LEVER WAS BUILT TO ARGUE FOR IS REFUTED, "
        "and the lever survives only as a scoped, default-OFF instrument. Its lane measured a "
        "48-molecule sweep at a 30s soft / 120s hard cap and concluded '47/48 produce nothing, "
        "the predicate is 48/48 with zero mismatches, geometry code is the only signal not "
        "refuted'. The fuller measurement on main (71443eda, tools/boron_gen_times.jsonl, 60s "
        "cap) contradicts all three: 2/33 DO generate -- RAWJEG (LIN, 1.2s) and ULODUU "
        "([Zr_TET], 61.8s). ULODUU is invisible to a 30s cap, so its success is BUDGET-DEPENDENT "
        "rather than either sweep being wrong, and that is itself the finding: the class boundary "
        "moves with the compute you give it. Cross-tabulating geometry against outcome kills the "
        "discriminator outright -- LIN 1 success / 3 failures, TET 1 / 4, everything else 0 / 15 "
        "-- so EVERY geometry with a success also has failures. Geometry separates nothing. "
        "The safe set is therefore {LIN, TET}, which preserves both real passes for 4 of the 25 "
        "recoverable cap-burners (PAQBOZ/PAQCAM/RIWKAK/RONQOD, all TET, all 60-64s). Excluding "
        "LIN costs nothing: its 3 failures (MODZUA/RANCIU/RULBUV) already fail in 0.00-0.01s on "
        "the uncoordinated-fragment path, so there is no budget there to reclaim. What is left "
        "is a predicate that fires where nothing has YET been seen to generate: an absence of "
        "evidence, not a rule. "
        "MEASURED (docs/agentic-notes/v0.4.7/BORON_GEN_CEILING_v0.4.7.md, "
        "tools/boron_gen_sweep34.sh + "
        "tools/boron_gen_sweep_14silent.sh, 48-molecule sample at the 30s cap): 40 molecules burn "
        "the whole embed budget (measured up to 2.3x OVER the requested cap; see the "
        "embed-time-budget finding below) before discovering they cannot assemble. That waste is "
        "real and is ~2.1 CPU-hours per 5,000-molecule sweep at the 300s budget. "
        "Counter-attributed (tools/perf_attribute.py, not clock-only): the cost "
        "is a full-complex PuLP/CBC bond-order-and-charge solve "
        "(chem.Molecule.get_valid_molecule(method='pulp') -> compute_chg_and_bo_pulp.py, called "
        "from embed.get_alternative_molecule, memoized only per dummy-metal `option` -- so "
        "priming each of the 3 default options costs a fresh unbounded CBC solve). A cage vertex "
        "has no 2c-2e Lewis structure, which is presumably why that solve is slow rather than "
        "terminating on a clean optimum. "
        "THE NAIVE VERSION OF THIS PREDICATE (bare motif check, no other gate) WAS MEASURED "
        "WRONG: RAWJEG_comp_0 ([Hg_LIN], a monodentate cage + Cl-, no haptic anything) carries "
        "the motif and STILL produces a structure in 2.5s. Hapticity "
        "and fragment/complex size were ALSO tried as discriminators and ALSO refuted (HAXJOG is "
        "monodentate like RAWJEG but fails; PEKQII/SEMTOV/VEJXOZ are smaller than every failing "
        "molecule in the 34-set but still fail), and geometry has since joined them via ULODUU. "
        "So all four candidate discriminators are now refuted by a counter-example: there is no "
        "known signal that separates the cages that assemble from the cages that do not. "
        "MODZUA/RANCIU carry the motif on an UNCOORDINATED fragment (0 binding "
        "vectors, an outer-sphere counterion) and already hit the existing, differently-labelled "
        "UncoordinatedFragmentError just as fast; the predicate requires >=1 vector so it does "
        "not relabel that already-correct failure. The predicate reads ONLY parsed.fragments / "
        "parsed.vectors / parsed.geo_code (the OIN STRING's own text and structured slot data, "
        "re-derived fresh here) -- never anything the GENERATOR constructs (no metal_complex, no "
        "embedded conformer) -- specifically so it cannot certify a defect the generator's own "
        "graph would paper over (the failure mode a sibling v0.4.7 lane hit with a bond-derived "
        "coordination check). Held OFF, and NOT CLOSE to promotion. The exclusion set rests on "
        "n=1 in each geometry (RAWJEG for LIN, ULODUU for TET), the sample is 48 molecules with "
        "no corpus-wide scan behind it, and -- decisively -- the discriminator the gate is built "
        "on has been refuted, so the next counter-example is expected rather than surprising. "
        "Promotion gate, which nobody should expect to clear as stated: a corpus-wide scan of "
        "every molecule carrying >=3 boron confirming zero produce a structure once flagged, PLUS "
        "a MECHANISM that explains why RAWJEG and ULODUU assemble when same-geometry siblings do "
        "not. Correlation at n=1 per cell is what produced the refuted rule in the first place. "
        "The compute this lever was meant to recover is better taken by v0.4.9's enforced "
        "generation bound, which reclaims it without needing a discriminator that does not exist."
    ),
    "OIN_ENFORCE_BUDGET": (
        "makes OIN3DGenerator(timeout=) a BOUND instead of a hint. Unset, the timeout becomes "
        "embed_time_budget, a deadline checked only at the TOP of the embed attempt loop "
        "(generator3d/__init__.py:350, :457, :529) -- so an in-flight attempt always runs to "
        "completion, and everything wrapping the loop is unbounded. Two direct-call probes "
        "measured the consequence: 60 s asked, 60.7-137.9 s spent on an eta sample (GOHWOQ "
        "2.3x) and 60.0-172.8 s on the boron set (2.9x). With the lever ON the same deadline is "
        "additionally checked inside get_embedding's two nested Python loops, after "
        "clean_geometry, and before the post-loop selection, and budget exhaustion with an "
        "EMPTY pool raises BudgetExhaustedError instead of degrading to a generic "
        "no_conformers. A non-empty pool is still returned -- a bound should stop work, not "
        "throw away a usable answer it already has. "
        "⚠ WHAT THE CHARTER GOT WRONG, and the reason this lever is shaped the way it is. "
        "v0.4.9 was chartered on 'max(elapsed_s) = 759.9 s against a 300 s budget'. That figure "
        "is arithmetic on a SUM: tools/test_dataset_roundtrip.py runs up to three separately "
        "SIGKILLed attempts per molecule and adds their wall-clock into one field. Split by "
        "tier_passed, all 4658 single-attempt rows in the 5k sweep finish at <= 300.2 s against "
        "a 300 s cap -- the harness enforces to eps ~ 0.2 s. The advisory-timeout defect is real "
        "(it is in the code, and the two probes above measure it directly), but 759.9 s is not "
        "evidence for it. See docs/agentic-notes/v0.4.9/ELAPSED_S_IS_A_SUM_v0.4.9.md. "
        "⚠ AND THE CHARTER'S SUSPECTS WERE WRONG TOO. It named the unbounded PuLP/CBC solve and "
        "the 48-57 s accept_fn re-encode as the likely sinks. Profiled on FOSNEI_comp_0 (the "
        "759.9 s worst case) at a 300 s budget: CBC is 1.74 s of 82.44 s (2.1%) -- the topology "
        "memo already collapsed it -- and the accept_fn re-encode is 0.63 s (0.8%). The sink is "
        "embed.get_embedding at 61.5 s of SELF time across 10 calls, plus clean_geometry.ff_clean "
        "at 15.2 s. That matters for the mechanism: get_embedding is a nested PYTHON loop over "
        "AllChem.EmbedMolecule, not one long native call, so a deadline checked inside it is real "
        "enforcement and eps is one EmbedMolecule, not one whole attempt. No fork/RLIMIT_CPU "
        "machinery is needed, and none is added. "
        "MEASURED: docs/agentic-notes/v0.4.9/BUDGET_BOUND_v0.4.9.md. "
        "Default OFF and NOT promoted in v0.4.9. Promotion is an accuracy decision, not a "
        "runtime one: 93.1% of honest passes already finish under 30 s, so a 30 s bound recovers "
        "37.8 CPU-h per 5000-molecule sweep but costs 251 passes = 5.02 points of byte_exact, "
        "against a headline goal of 100%. At 300 s it costs 3 passes for 3.1 CPU-h. The right "
        "default is whatever the release that OWNS the accuracy/compute trade decides, with the "
        "cap-vs-cost curve in front of it -- not this one, which only makes the trade expressible."
    ),
    "OIN_MEMO_CIP_REPARSE": (
        "memoises core.chirality._reparse_cip_label_once on its own arguments. v0.4.10 is a "
        "byte-identical-by-construction release and this is the change that fits its rule most "
        "exactly: the function's three arguments are immutable scalars (smiles, probe, "
        "fill_deficit), its body builds a FRESH rdkit mol from the SMILES on every call, it "
        "touches no module or global state, and it returns a str or None. So a memo returns "
        "precisely what the function would have computed -- only the wall-clock changes. That is "
        "the same sentence compute_chg_and_bo_pulp's topology memo is justified by. "
        "WHY IT HITS RATHER THAN MERELY BEING SAFE: both call sites key on a SMILES derived from "
        "the OIN template (metallogen_adapter._template_sp3_label, reached per conformer through "
        "accept_fn -> _reencode_key_matches -> build_contract_mol) or from a metal-free fragment "
        "(_reparse_aromatic_cip_label). Neither string carries coordinates, so every conformer of "
        "a molecule generates the SAME key and the repeat traffic is structural, not incidental. "
        "MEASURED (docs/agentic-notes/v0.4.9/BUDGET_BOUND_v0.4.9.md §2): on VAFMIA_comp_0 "
        "([Cu_LIN], bis-adamantyl NHC) this function is 77.78 s of SELF time across 32 calls at "
        "2.43 s each -- 99% of that molecule's generation. It is also what sets v0.4.9's epsilon: "
        "OIN_ENFORCE_BUDGET can decline to START the next accept_fn call but cannot interrupt the "
        "one in flight, and that in-flight call is ~24 s of this function, which is why the bound "
        "holds to 2.09x rather than 1.0x. "
        "⚠ THE COST IS BIMODAL BY MOLECULE, which is the trap this release inherited from v0.4.9: "
        "the same accept_fn re-encode is 62% of VAFMIA and 0.8% of FOSNEI_comp_0. A lever aimed "
        "at whichever function profiled expensive last optimises one molecule class and measures "
        "nothing on the corpus -- v0.4.9 shipped that mistake once (epsilon barely moved) and "
        "recorded it. Do not read a VAFMIA number as a corpus number. "
        "Cross-molecule retention is safe by the purity argument above; the cache is bounded at "
        "_CIP_REPARSE_MEMO_MAX so a long single-interpreter sweep cannot grow it without limit, "
        "and _reparse_cip_memo_clear() lets a per-molecule gate guarantee isolation the way "
        "_ac2bo_memo_clear() does for perception."
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
