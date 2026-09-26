#!/usr/bin/env python
"""Freeze a release's comparison numbers into the tracked ``measurements/`` tree.

WHY THIS EXISTS
===============
``tmCAT-tmPHOTO_xyz_dataset/`` is **gitignored in its entirety** (``.gitignore:85``) -- the
2.1 GB dataset and every ``results-*/`` sweep. ``spec/handoffs/`` is gitignored too. So for six
consecutive releases the project wrote excellent prose about numbers that lived on exactly one
disk, and:

    🔴 v0.4.11's mirror-audit JSON -- the 19/250 measurement that chartered the whole of
       v0.4.12 -- is GONE. It was recoverable only because the TOOL was committed and
       deterministic, at ~50 minutes of CPU to re-run.

``CLOSEOUT.md`` §4 mandated git-durable *prose* and said nothing about *data*, which is exactly
how the lapse went unnoticed: v0.4.3 committed 2 data files, v0.4.6 committed 5, and
v0.4.7-v0.4.12 committed **zero**.

This tool is the other half of that ritual step. It is meant to be run automatically from
CLOSEOUT §4b and manually via ``/freeze-measurements``.

⚠ THE TREE IS PUBLIC. ``origin`` is a public GitHub repo and a sibling session pushes ``main``
regularly, so anything harvested here becomes public and permanent. Two consequences: local
absolute paths must never leak (checked below), and size discipline is not cosmetic -- hence
the caps.

WHERE IT WRITES, AND WHY NOT ``$PWD``
=====================================
Always ``<main checkout>/measurements/``, resolved from ``git rev-parse --git-common-dir``,
**never** the current worktree. A worktree's files vanish on ``git worktree remove`` -- which is
how this session nearly lost its own data, having created and then removed
``../oin-v0412-release``. Untracked files also do not survive ``git clean -fd``. The main
checkout plus a commit is the only durable destination.

Same resolution trick v0.4.9 used to stop ``gate_v047.sh`` silently selecting a sibling
project's venv.

SELECTION IS AN ALLOWLIST, NOT A SIZE SWEEP
===========================================
Measured: a naive "every small .json/.md under a cap" rule selects **12,756 files / 19 MB** of
per-molecule probe dumps and case registries. The allowlist below yields **29 files / 838 KB**
across the whole historical tree -- the artifacts a *later release actually diffs*.

Usage
=====
    python tools/harvest_measurements.py --backfill --dry-run
    python tools/harvest_measurements.py --release v0.4.12 --from <dir> [--from <dir2>]
    python tools/harvest_measurements.py --link          # read-access symlink in a worktree
"""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

#: Filename patterns worth preserving forever. Each is a *comparison* artifact: a later release
#: diffs it to learn whether it regressed. Report-backing evidence for a single write-up belongs
#: in ``docs/agentic-notes/<release>/`` instead -- see ``measurements/README.md`` for the test.
ALLOW = [
    "bucket_report*.md",
    "FROZEN.md",
    "SOURCE",
    "MANIFEST*",
    "VALIDATION.md",
    "README.md",
    "triage_hard_fails.md",
    "transitions.json",
    "*_DONE",
    "*TREND.tsv",
    "CASE_REGISTRY.md",
    # v0.4.12 onward: the per-release instruments' own output, harvested from a scratchpad.
    "mirror_audit*.json",
    "transition_*.json",
    "ab_*.json",
    "cohort_*.json",
    "audit_*.json",
    # v0.4.13. Added because this release's instruments are EXACTLY the artifact class whose
    # loss motivated this tool -- the two mirror arms are the same kind of file as v0.4.11's
    # vanished 19/250 audit, and `fold_transition_veto.json` names all 171 molecules behind the
    # +3.42 headline, which no prose in `docs/` reproduces.
    # ⚠ `mirror_*` deliberately, NOT `mirror_arm*`. The narrower pattern was written first and
    # silently dropped `mirror_cat_{promoted,noveto}.json` -- the arms that reproduce v0.4.12's
    # PUBLISHED gate (19 -> 0), i.e. half the evidence the v0.4.13 promotion rests on. Caught only
    # because the dry-run's file list was read line by line against the scratchpad. A too-narrow
    # allowlist fails exactly like a broken instrument: it prints a plausible total.
    "mirror_*.json",
    "fold_*.json",
    "attach_class_audit.json",
    "prefilter_*.json",
    # v0.4.14. The end-to-end A/B outputs are the load-bearing ones: `reso_full_ab.json` names
    # every one of the 78 gains and 7 losses behind the +1.42 headline, and `v0413_atrisk_ab.json`
    # is the measurement that corrected v0.4.13's published +3.42 to ~+2.82. Neither is reproducible
    # from prose. `reso_movers_exact.json` is the coordinate-derived affected population -- the
    # thing that made a ~55 CPU-h sweep unnecessary -- and a later release re-deriving it is exactly
    # the diff this tree exists to support.
    "veto_*.json",
    "resonance_*.json",
    "reso_*.json",
    "*_ab.json",
    "generator_ab*.json",
    # v0.4.14 baseline sweep. The run itself is 261 MB and gitignored; these are the parts a later
    # release actually diffs. `per_molecule_extract.tsv` is the load-bearing one -- it re-derives the
    # bucket table, the key_equal sub-split and the runtime percentiles WITHOUT the raw run, which
    # is the whole point of this tree (v0.4.11's mirror-audit JSON is gone for want of exactly this).
    "per_molecule_extract.tsv",
    # v0.4.14 residue: the golden re-freeze rows (provenance for which arm2 rows moved and why)
    # and the SUPERSEDED scored A/B, kept precisely because it is the artifact behind a documented
    # error -- it was scored with arm2's circular predicate and produced a wrong conclusion.
    "arm2_refreeze_rows.tsv",
    "gain_generator_ab.tsv",
    "RUN.md",
    # v0.4.15. ⚠ `attach_*`, not `attach_class_audit.json` -- the narrow name was already in this
    # list and would have silently dropped `attach_return_preflight.json`, which is the file that
    # sized Lane 1's exposure at 1/52 rather than the charter's 52/52. That is the same too-narrow
    # allowlist failure the `mirror_arm*` note above records, one release later.
    "attach_*.json",
    "string_exact_*.json",
    # The v0.4.15 three-arm A/B output. ⚠ CHECKED BEFORE THE ARMS FINISHED, and it is as well:
    # none of `*_ab.json`, `generator_ab*.json`, `ab_*.json` or `cohort_*.json` matches
    # `lane1_pop_*.json`, so the release's LOAD-BEARING measurement would have been dropped in
    # silence -- the third time this allowlist has been one pattern too narrow (see the
    # `mirror_arm*` and `attach_class_audit.json` notes above). Verified by running the real
    # fnmatch over the real filenames rather than by reading the list.
    "lane1_*.json",
    "lane2_*.json",
    "both_*.json",
    # v0.4.16. `knee_*` is the recovered-vs-bound evidence behind the OIN_STRING_EXACT_BOUND
    # promote/hold call -- the per-molecule `min_bound` ordinals that no prose reproduces, and the
    # only artifact a later release can diff to ask whether the knee moved. `pop_*.txt` is the
    # SAMPLE MEMBERSHIP: a rate without its sample is not reproducible, and v0.4.15's own
    # populations live in measurements/ for exactly that reason.
    #
    # ⚠ VERIFIED BY RUNNING fnmatch OVER THIS RELEASE'S ACTUAL FILENAMES, not by reading the list.
    # Checked before adding: `knee_shard1.json`, `knee_curve_v0416.json` and `pop_all.txt` matched
    # NOTHING here, so a fourth consecutive release would have silently dropped its evidence while
    # the dry-run printed a plausible total. That is the same failure mode as a broken instrument.
    "knee_*.json",
    "pop_*.txt",
    # v0.4.16, second pass. Both of these were written into `measurements/` BY HAND because the
    # harvester would not have taken them, which is the same failure the ALLOW list keeps having --
    # noticed here only because the patterns were re-run through fnmatch rather than eyeballed.
    #
    # `cohort-*_manifest.json` is the SAMPLE DEFINITION: seed, n, dedup priority and all N molecule
    # names. Every accuracy figure this project has ever published is a rate over one of these, and
    # a rate without its sample is not reproducible. `MANIFEST*` did not match (wrong case, wrong
    # position) and neither did `cohort_*.json` (underscore, not hyphen).
    # `sweep_extract_*.jsonl.gz` is a whole sweep's per-molecule rows -- the thing that lets a later
    # release re-derive a bucket table instead of re-quoting prose.
    "cohort-*_manifest.json",
    "sweep_extract_*.jsonl.gz",
    # v0.4.17 census (C1-C4). The per-molecule verdicts of four instruments that never touched the
    # generator, gzipped so each sits under PER_FILE_CAP (raw they are 1-3.5 MB each), plus the
    # 5,000-row attribution table they join into and the small cross-tabs. A later release asks
    # "did the fault partition move?" and can only answer it by diffing THESE, per molecule.
    # ⚠ VERIFIED BY RUNNING fnmatch OVER results-census/ (the `census_*` prefix is written by
    # tools/census/attribution_table.py --gz precisely so ONE pattern covers the set).
    "census_*.jsonl.gz",
    "attribution_table.tsv.gz",
    "attribution_summary.json",
    "collisions.json",
    "string_sufficiency_summary.json",
    "g_crosstab.json",
    # v0.4.17 L1 (exact donor fold) + the v0.4.17 baseline sweep. Everything is staged by
    # tools/v0417/freeze_stage.py under ONE prefix, scrubbed BEFORE it is gzipped -- this tool
    # copies a .gz verbatim, so a path inside one would be published unscrubbed.
    # ⚠ VERIFIED BY RUNNING fnmatch OVER THE STAGED FILENAMES, not by reading this list: before
    # these four lines, 19 of the 20 staged files matched NOTHING (only `bucket_report_honest.md`
    # did), and the dry run would have printed a plausible 15 KB total for a release whose
    # evidence is 1 MB. `v0417_pop_*.txt` also does not match `pop_*.txt` -- a prefix moves it.
    "v0417_*.jsonl.gz",
    "v0417_*.json.gz",
    "v0417_*.json",
    "v0417_*.txt",
    "v0417_*.md",
    "v0417_*.tsv",
    "v0417_*.tsv.gz",
    # v0.4.18: the census re-run on the v0.4.17 sweep, staged by tools/census/stage_reattribution.py
    # under ONE prefix. Its table cannot be called `attribution_table.tsv.gz`: PROVENANCE is keyed
    # on the filename alone and that name already says "the v0.4.14 sweep".
    "v0418_census_*",
    # v0.4.18 L2 (eta construction), staged by tools/v0418/freeze_stage.py under ONE prefix.
    "v0418_l2_*",
    # v0.4.18 promotion: the sweep of record + the /refreeze-goldens run, staged by
    # tools/v0418/freeze_stage_sweep.py into their OWN release (v0.4.18-sweep).
    "v0418_sweep_*",
    "v0418_arm2_*",
    # v0.4.18 RELEASE: the release sweep (both promotions), its /refreeze-goldens run, and the
    # census re-run on it (tools/v0418/freeze_stage_sweep.py --release; stage_reattribution.py
    # --results results-v0.4.18-release-census --prefix v0418_rcensus_).
    "v0418_rsweep_*",
    "v0418_rarm2_*",
    "v0418_rcensus_*",
    # v0.4.18 selection lane (OIN_ETA_RETARGET), staged by tools/v0418/freeze_stage_selection.py.
    "v0418_sel_*",
    # v0.4.19 serializer lane (OIN_H_FAITHFUL + OIN_RC1_PROPAGATE, measured, NOT promoted), staged by
    # tools/v0419/freeze_stage.py.
    "v0419_ser_*",
    # v0.4.19 RELEASE: the candidate sweep (= the release sweep, the three serializer levers ON),
    # its post-processing under the levers and the /refreeze-goldens run
    # (tools/v0419/freeze_stage_release.py, release v0.4.19-release-sweep), and the census on it
    # (stage_reattribution.py --results results-v0.4.19-release-census --prefix v0419_rcensus_).
    "v0419_rsweep_*",
    "v0419_rarm2_*",
    "v0419_rcensus_*",
    # v0.4.19 E2 lane (OIN_N_VALENCE_2 + OIN_CANONICAL_RESONANCE, perception order; measured, NOT
    # promoted), staged by tools/v0419/freeze_stage_e2.py, release v0.4.19-e2.
    "v0419_e2_*",
]

#: Directories that are raw inputs or bulk per-molecule output. Never harvested.
PRUNE_DIRS = {
    "individual_reports",
    "structures",
    "cat",
    "photo",
    "rebaseline_inputs",
    "regression_inputs",
    "reports",
    "__pycache__",
    ".git",
}

#: Extensions/names that are never comparison artifacts regardless of the allowlist.
DENY = ["*.log", "*.xyz", "*.oin", "*.bak", "*.pre-rebuild-bak", "worker_pids.txt", "*.pyc"]

PER_FILE_CAP = 512 * 1024
TOTAL_CAP = 5 * 1024 * 1024

#: Provenance for names this project's tools produce. Anything unmatched is reported
#: ``UNKNOWN`` rather than given a plausible-looking guess -- a figure without its command is an
#: order of magnitude, not a measurement.
PROVENANCE = [
    (
        r"^census_g_verdict(_control-mirror)?\.jsonl\.gz$|^g_crosstab\.json$",
        "tools/census/g_vs_input.py [--control mirror] [--crosstab]"
        "   (census C1: the neutral ruler -- element-labelled distance graph + donor-direction"
        " spheres + signed volumes, NO oinsmiles import -- judges every generated structure"
        " against its input; the mirror control gives every input's ruler chirality)",
    ),
    (
        r"^census_e_selfconsistency\.jsonl\.gz$|^census_veto_probe(_chiral|_renumber)?\.jsonl\.gz$",
        "tools/census/e_selfconsistency.py / tools/census/veto_probe.py"
        "   (census C2: eleven encodes per input -- identity, rewrite, rotation, 3 renumberings,"
        " 3 x 0.02 A noise, mirror; the veto probe re-encodes with OIN_FOLD_PARITY_VETO and the"
        " fold off to name the mechanism behind achiral molecules with two strings)",
    ),
    (
        r"^census_parseback(_gen)?\.jsonl\.gz$|^census_pflags\.jsonl\.gz$"
        r"|^collisions\.json$|^string_sufficiency_summary\.json$",
        "tools/census/string_sufficiency.py {parseback [--side gen] | collide | pflags"
        " --charge-probe | summary}   (census C3: smiles_1 read back to a graph with NO 3D vs"
        " the ruler's input graph; natural-twin collision scan; perception flags with the"
        " stated charge honoured through a patched entry point)",
    ),
    (
        r"^mirror_probe\.json$",
        "tools/census/mirror_probe.py   (census C1: does E(mirror x) differ from E(x) on a"
        " sample of ruler-achiral vs ruler-chiral inputs -- the estimate C2 then measured over"
        " the whole cohort)",
    ),
    (
        r"^v0417_sweep_(RUN\.md|run_config\.json)$",
        "tools/v0417/launch_sweep.sh -> tools/run_sweep.sh <cohort-v0.4.5-5k> <out> 6 300"
        "   (THE v0.4.17 BASELINE SWEEP: commit 814abff3, SHIPPED DEFAULTS -- the lever block in"
        " run_config.json is EMPTY on purpose, that is the thing under test -- 6 shards 1-BASED,"
        " --mol-timeout 300, BLAS=1. RUN.md is hand-written provenance)",
    ),
    (
        r"^v0417_sweep_two_numbers\.txt$",
        "tools/v0417/post_sweep.sh -> tools/v0417/sweep_two_numbers.py --sweep <sweep>"
        "   (THE v0.4.17 HEADLINE: self-consistent 4136/5000 = 82.72%, VERIFIED 3730/5000 = 74.60%."
        " The predicate's control runs first and must reproduce the census on the v0.4.14 sweep of"
        " record, 5000 / 3858 / 3462, or nothing is printed)",
    ),
    (
        r"^v0417_(sweep|ab_(off|on))_g_verdict\.jsonl\.gz$",
        "tools/census/g_vs_input.py --cohort <cohort> --sweep <run> --out <run>"
        "   (the census's neutral ruler on the GENERATED STRUCTURES of {0}: the structure half of"
        " the VERIFIED predicate. No oinsmiles import)",
    ),
    (
        r"^v0417_sweep_parseback\.jsonl\.gz$",
        "tools/census/string_sufficiency.py parseback --sweep <sweep> --out <sweep>"
        "   (smiles_1 of the NEW sweep read back to a graph with no 3D: the string half of the"
        " VERIFIED predicate, RE-DERIVED rather than carried from the census -- 1 of 5000 verdicts"
        " moved under the relabeling)",
    ),
    (
        r"^v0417_autofold_audit\.json\.gz$",
        "tools/v0417/autofold_audit.py   (L1a OFFLINE GATE, exact because the slot post-pass is a"
        " pure string function: bucket fold vs exact fold applied to the census's stored"
        " rotation-only strings. Positive control 5792/5792; 550 of the 568 vetoed achiral pairs"
        " unified; 18 kept split are chiral by wrap; 17 of the 36 protected pairs unified)",
    ),
    (
        r"^v0417_automorphism_extension_check\.json$",
        "tools/v0417/automorphism_extension_check.py   (does pendant stripping invent a symmetry?"
        " 13716 pruned automorphisms over 7968 fragment checks, 0 without a full-graph extension;"
        " negative control refused 7953/7953)",
    ),
    (
        r"^v0417_e_selfconsistency_exact\.jsonl\.gz$",
        "OIN_EXACT_DONOR_FOLD=1 tools/census/e_selfconsistency.py --name e_selfconsistency_exact"
        "   (census C2 re-run unchanged under the lever, same seeds: eleven presentations of every"
        " input. Diff against census_e_selfconsistency.jsonl.gz: achiral-and-differs 670 -> 120,"
        " slot-level renumber drift 494 -> 33, KEY-level 255 -> 255 with membership changed on 0)",
    ),
    (
        r"^v0417_generated_side_movers\.jsonl\.gz$|^v0417_pop_generated_side_movers\.txt$",
        "tools/v0417/generated_side_movers.py   (every stored generated structure of the v0.4.14"
        " sweep encoded lever OFF and ON: 286 of 4743 move. OFF == the sweep of record on"
        " 4743/4743)",
    ),
    (
        r"^v0417_pop_(input_side_movers|ab_union_movers)\.txt$",
        "tools/v0417/exact_fold_live_report.py (input side, 392) / union with the generated side"
        "   (SAMPLE MEMBERSHIP of the generator A/B: 504 = 392 u 286, 174 in both. Every other"
        " molecule gets the same two strings in both arms and is unchanged by construction --"
        " which the sweep then measured: 4494 of 4496)",
    ),
    (
        r"^v0417_generator_ab_report\.txt$|^v0417_ab_(gains|losses)\.txt$"
        r"|^v0417_ab_(off|on)_bucket_report_honest\.json\.gz$",
        "tools/v0417/run_generator_ab.sh <union movers> ; tools/v0417/post_generator_ab.sh"
        "   (L1b: the REAL sweep harness over the 504 movers, OFF and ON arms run SIMULTANEOUSLY,"
        " commit de4a9d2c. Self-consistent 276 gains / 3 losses, VERIFIED 269 / 4; dead-lever"
        " check 392; noise floor 3/504 vs the sweep of record)",
    ),
    (
        r"^v0417_arm2_refreeze_rows\.tsv$",
        "tools/v0417/arm2_refreeze.py splice --write   (THE RE-FREEZE: 51 of 425 ARM 2 golden rows,"
        " one reason each -- LEVER 36 / STALE+LEVER 4 / STALE 11. Columns: old row, lever-off row,"
        " re-frozen fields 2-6; fields 7+ of a golden are preserved, never spliced)",
    ),
    (
        r"^v0417_arm2_field2_audit_v04[79]\.jsonl\.gz$",
        "tools/v0417/arm2_field2_audit.py --golden <golden> --cohort-dir <cohort>   (field 2 of"
        " EVERY golden row, encode-only, fresh process per lever setting: shipped vs"
        " OIN_EXACT_DONOR_FOLD=0. 15 of 425 rows were stale BEFORE v0.4.17; a worker importing a"
        " foreign oinsmiles aborts the run, as does a lever that fires on 0 rows)",
    ),
    (
        r"^v0417_arm2_stale_cause\.(jsonl\.gz|txt)$",
        "tools/v0417/arm2_stale_cause.py   (which lever left the rows stale: the pre-lane encoder,"
        " then each other default-ON lever set to 0 ALONE, 14 molecules x 13 arms."
        " OIN_CANONICAL_DONOR_FOLD=0 restores 14 of 14 -> stale since v0.4.13)",
    ),
    (
        r"^v0417_arm2_(rows_on_full_gate\.tsv\.gz|diff_full_gate_vs_old_goldens\.txt)$",
        "tools/v0417/run_arm2_refreeze.sh on ; tools/v0417/arm2_refreeze.py diff|rows   (the FULL"
        " gate, all 425 rows of both ARM 2 goldens, shipped defaults, 6 shards, commit f5f507ef,"
        " against the goldens as they stood: 38 field-2 + 13 field-3-only mismatches, 2 killed)",
    ),
    (
        r"^v0417_arm2_rows_off_control\.tsv$",
        "tools/v0417/run_arm2_refreeze.sh off   (OIN_EXACT_DONOR_FOLD=0 on the 53 unreproduced"
        " rows, commit 6cf9b185: a healthy old row comes back byte-for-byte -- 36 did)",
    ),
    (
        r"^v0417_arm2_rows_verify_real_gate\.tsv$",
        "tools/v0417/run_arm2_refreeze.sh verify   (tools/gate_v047.sh arm2 ITSELF reading the"
        " re-frozen goldens, commit 1a3a1075: 51 of 51 re-frozen rows pass; the only mismatches are"
        " EQEROI and MUKGUW, SIGKILLed at --hard-timeout in all three runs and never re-frozen)",
    ),
    (
        r"^v0418_census_(attribution_table\.tsv\.gz|attribution_summary\.json|attribution\.txt)$",
        "tools/census/attribution_table.py --sweep <results-v0.4.17-sweep> --out <reattr> + one path"
        " per instrument   (THE CENSUS RE-RUN ON THE v0.4.17 SWEEP: 5000 rows, one fault each,"
        " UNATTRIBUTED 0. FAIL 864 = 17.28 pts: G 9.62 / E2 3.50 / E1 3.38 / P 0.58 / data 0.20."
        " NONE = 3730 = the VERIFIED count sweep_two_numbers.py printed, or the tool aborts)",
    ),
    (
        r"^v0418_census_control_on_the_record_sweep\.txt$",
        "tools/census/attribution_table.py --write-to <scratch>   (CONTROL, run first: the"
        " parameterised tool on the v0.4.14 record sweep reproduces measurements/v0.4.17-census/"
        "attribution_table.tsv.gz BYTE-FOR-BYTE)",
    ),
    (
        r"^v0418_census_reattribution_diff\.(txt|json)$",
        "tools/census/reattribution_diff.py --old <census table> --new <new table> --movers <504>"
        "   (where each fault WENT, split by mover: 4482 of 4496 non-movers keep theirs -- the"
        " noise floor of attribution is 14 molecules / 0.28 pts)",
    ),
    (
        r"^v0418_census_(parseback_gen|pflags)\.jsonl\.gz$|^v0418_census_attach_class_audit\.json\.gz$",
        "tools/census/string_sufficiency.py parseback --side gen | pflags --charge-probe ;"
        " tools/attach_class_audit.py   (the three SWEEP-DEPENDENT instruments, derived on the"
        " v0.4.17 sweep. The other inputs are frozen where they were made: g_verdict + parseback in"
        " v0.4.17-sweep/, e_selfconsistency_exact in v0.4.17/, the mirror ruler + twin clusters in"
        " v0.4.17-census/)",
    ),
    (
        r"^v0418_sel_selection_sim\.(txt|json)$",
        "tools/v0418/selection_sim.py   (OFFLINE BOUND, and it CAN express a loss: the L2 A/B left"
        " three structures per eta molecule with KNOWN verdicts. An input-free 'every binding atom"
        " the string DECLARES is inside the contact cutoff' choice between them = 668 -> 723 verified"
        " (+57 / -2), oracle 740. Adding SURPLUS contacts as a criterion costs 20 losses)",
    ),
    (
        r"^v0418_sel_(retarget_probe\.txt|probe5_classes\.json)$",
        "tools/v0418/eta_path_probe.py --arms dg,retarget   (fresh-process SAMPLE: 63 oracle-headroom"
        " molecules + 40 verified controls. Honest round trip 22 -> 52 (+31 / -1), controls 40/40;"
        " lever OFF reproduces the v0.4.18 sweep on 103/103 structures. A SAMPLE, NOT A VERDICT)",
    ),
    (
        r"^v0418_sel_eta_ab_report_retarget\.(txt|json)$",
        "tools/v0418/run_selection_ab.sh -> post_eta_ab.sh retarget -> eta_ab_report.py --on-arm"
        " retarget --base 4305,3883   (THE HARNESS A/B of OIN_ETA_RETARGET over all 1146 eta-bound"
        " molecules, commit 2f0b68d5: self-consistent +49/-4, VERIFIED +43/-5 -> 87.00% / 78.42%"
        " projected; excluding rows on the 300 s budget boundary +44/-4 and +40/-5. THE OFF ARM WAS"
        " BUILT from the v0.4.18 sweep of record, not run, so the report's noise-floor line is"
        " circular here. NO DEFAULT CHANGED: the lever ships OFF)",
    ),
    (
        r"^v0419_ser_pop1b?_report\.txt$|^v0419_ser_pop1b?\.jsonl\.gz$",
        "tools/v0419/lever_parseback.py --names <pop1_names> --classes <pop1_classes> --arms"
        " shipped,hfaith,rc1off,both | shipped,rc1off,rc1prop,fix2   (THE PARSE-BACK PROBE: 296"
        " molecules from the v0.4.18 release census -- E1_HCOUNT 84, E1_GRAPH/SPHERE_DIFF 76,"
        " E1_GRAPH/LIGAND_DIFF 56, 80 verified controls -- encoded live in a fresh process per arm,"
        " string parsed back at the adapter level with the census's ruler. Control: live shipped =="
        " sweep smiles_1 296/296. OIN_H_FAITHFUL=1 repairs 56/84 E1_HCOUNT; OIN_RC1_PROPAGATE=1"
        " repairs 14/76 SPHERE_DIFF, byte-identical to forcing RC1's fail-safe; controls 0 changed)",
    ),
    (
        r"^v0419_ser_changed_set\.txt$|^v0419_ser_changed_single\.jsonl\.gz$",
        "tools/v0419/serializer_ab_report.py changed --write-set ; tools/v0419/lever_parseback.py"
        " --arms shipped,hfaith,rc1prop,fix2 over that set   (the COMPLETE changed-string set: base"
        " under both levers != release-sweep smiles_1 on 80 of 5,000 -- H_FAITHFUL 64, RC1_PROPAGATE"
        " 15, 1 encode-timeout row -- and 4,909 unchanged incl. every verified pass)",
    ),
    (
        r"^v0419_ser_esc_fix2_summary\.txt$|^v0419_ser_e_selfconsistency_fix2\.jsonl\.gz$",
        "OIN_H_FAITHFUL=1 OIN_RC1_PROPAGATE=1 tools/census/e_selfconsistency.py --sweep"
        " <results-v0.4.18-release-sweep> --cpu 8   (THE CANONICALITY AUDIT under both levers, all"
        " 5,000 x 11 encodes, 13.85 CPU-h: determinism, rewrite, rotation, 3 renumberings, 3 noise,"
        " mirror. Against the shipped run (measurements/v0.4.17-exactfold e_selfconsistency): every"
        " column identical, mirror 2x2 unchanged, 0 became fragile)",
    ),
    (
        r"^v0419_ser_rescore_fix2\.jsonl\.gz$",
        "OIN_H_FAITHFUL=1 OIN_RC1_PROPAGATE=1 tools/honest_rescore.py --results-dir"
        " <results-v0.4.18-release-sweep> --workers 4 --timeout 300   (the OFFLINE re-score: every"
        " stored generated structure re-encoded under the levers; joined with the lever base strings"
        " for the 4,909 unchanged rows: 0 self-consistent losses, 5 false-pass gains, VERIFIED 3920"
        " -> 3920)",
    ),
    (
        r"^v0419_ser_serializer_ab_report\.(txt|json)$|^v0419_ser_ab_commits\.tsv$",
        "tools/v0419/run_serializer_ab.sh <changed_set> -> post_serializer_ab.sh ->"
        " serializer_ab_report.py all   (THE LIVE A/B over the 80 changed-string molecules, ON arm"
        " commit d3adba15, OFF = the release sweep's own rows: self-consistent +29/-1, VERIFIED"
        " +24/-0 with each arm's parse-back under its own reader; projection 86.94 -> 87.60% /"
        " 78.40 -> 78.88%. READER COUPLING: the same ON strings verify on 15 rows under the SHIPPED"
        " reader (OIN_H_FAITHFUL gates metallogen_adapter too). NO DEFAULT CHANGED)",
    ),
    (
        r"^v0419_ser_ab_on_(bucket_report_honest\.json|g_verdict\.jsonl|parseback\.jsonl|parseback_shippedreader\.jsonl)\.gz$",
        "tools/roundtrip_bucket_report.py --score honest ; tools/census/g_vs_input.py ;"
        " tools/census/string_sufficiency.py parseback under the ON levers | under the shipped"
        " reader   (the ON arm's own instruments; the OFF arm's are the release sweep's, frozen in"
        " measurements/v0.4.18-release-sweep)",
    ),
    (
        r"^v0419_ser_pop1c_report\.txt$|^v0419_ser_pop1c\.jsonl\.gz$",
        "tools/v0419/lever_parseback.py --arms shipped,cap,fix3 over the same 296 rows   (the"
        " parse-back probe for OIN_CAP_IGNORES_METAL, hydrogen exempt: LIGAND_DIFF 14/56,"
        " SPHERE_DIFF 7/76, HCOUNT 3/84 repaired; with the two serializer levers 95 of 216 E1"
        " strings describe their input; controls 80/80, 0 errors)",
    ),
    (
        r"^v0419_ser_changed_set_(cap|fix3)\.txt$|^v0419_ser_changed_single_(cap|fix3)\.jsonl\.gz$",
        "tools/v0419/serializer_ab_report.py changed --esc <arm esc> --write-set ;"
        " lever_parseback.py --arms shipped,hfaith,rc1prop,cap over that set   (the changed-string"
        " set of the named arm -- cap = OIN_CAP_IGNORES_METAL alone, fix3 = all three levers -- and"
        " the per-lever attribution of each moved string)",
    ),
    (
        r"^v0419_ser_esc_(cap|fix3)_summary\.txt$|^v0419_ser_e_selfconsistency_(cap|fix3)\.jsonl\.gz$",
        "<arm levers> tools/census/e_selfconsistency.py --sweep <results-v0.4.18-release-sweep>"
        " --cpu 6   (the canonicality audit under the named arm's levers, all 5,000 x 11 encodes)",
    ),
    (
        r"^v0419_ser_rescore_(cap|fix3)\.jsonl\.gz$",
        "<arm levers> tools/honest_rescore.py --results-dir <results-v0.4.18-release-sweep>"
        " --timeout 300   (the offline re-score under the named arm's levers)",
    ),
    (
        r"^v0419_ser_serializer_ab_report_(cap|fix3)\.(txt|json)$|^v0419_ser_ab_commits_(cap|fix3)\.tsv$",
        "tools/v0419/run_serializer_ab.sh <changed_set> <arm> -> post_serializer_ab.sh <arm> ->"
        " serializer_ab_report.py all   (the live A/B of the named arm over its changed-string set;"
        " OFF = the release sweep's own rows; each arm's parse-back under its own reader. NO"
        " DEFAULT CHANGED)",
    ),
    (
        r"^v0419_ser_ab_(cap|fix3)_on_(bucket_report_honest\.json|g_verdict\.jsonl|parseback\.jsonl|parseback_shippedreader\.jsonl)\.gz$",
        "tools/roundtrip_bucket_report.py --score honest ; tools/census/g_vs_input.py ;"
        " tools/census/string_sufficiency.py parseback under the arm's levers | under the shipped"
        " reader   (the named arm's own instruments)",
    ),
    (
        r"^v0418_sel_ab_(off|retarget)_(g_verdict\.jsonl|bucket_report_honest\.json)\.gz$",
        "tools/roundtrip_bucket_report.py --score honest ; tools/census/g_vs_input.py   (per arm: the"
        " honest bucket and the neutral ruler's verdict for every one of the 1146; ab_off is the"
        " v0.4.18 sweep's own rows. With v0.4.18-census/v0418_census_attribution_table.tsv.gz these"
        " re-derive the headline: tools/v0418/freeze_stage_selection.py --verify <this dir>)",
    ),
    (
        r"^v0418_sel_ab_rows\.tsv\.gz$|^v0418_sel_ab_commits\.tsv$",
        "tools/v0418/freeze_stage_selection.py   (MADE at staging time: per molecule and arm the"
        " sha256 of the generated structure -- the lever changed 158 of 1067 -- coordination.intact"
        " and elapsed_s, plus the sweep of record's sha256; and where each arm came from)",
    ),
    (
        r"^v0418_rsweep_(RUN\.md|run_config\.json)$",
        "SWEEP_TAG=v0.4.18-release tools/v0418/launch_sweep.sh -> tools/run_sweep.sh"
        " <cohort-v0.4.5-5k> <out> 6 300   (THE v0.4.18 RELEASE SWEEP -- both eta promotions,"
        " OIN_ETA_COVALENT_TARGET + OIN_VDW_EXEMPT_BINDING and OIN_ETA_RETARGET, as SHIPPED"
        " DEFAULTS: the lever block of run_config.json is EMPTY. Commit cdb76eca, 6 shards 1-BASED,"
        " --mol-timeout 300, BLAS=1. RUN.md is hand-written provenance)",
    ),
    (
        r"^v0418_rsweep_(two_numbers\.txt|bucket_report_honest\.md)$",
        "SWEEP_TAG=v0.4.18-release tools/v0418/post_sweep.sh -> roundtrip_bucket_report.py --score"
        " honest ; tools/v0417/sweep_two_numbers.py --sweep <sweep>   (THE v0.4.18 HEADLINE, see"
        " RUN.md for the two numbers. The '(record ...)' labels INSIDE two_numbers.txt are the"
        " v0.4.14 record's)",
    ),
    (
        r"^v0418_rsweep_vs_ab\.(txt|json)$",
        "tools/v0418/sweep_vs_ab.py --record results-v0.4.18-sweep --ab-on"
        " results-v0.4.18-eta-selection/ab_retarget --predicted 4305+45,3883+38   (CHECKS the two"
        " claims the selection A/B's projection stood on: non-eta structures byte-identical to the"
        " sweep it replaces, eta structures byte-identical to the A/B's retarget arm)",
    ),
    (
        r"^v0418_rsweep_(g_verdict|parseback)\.jsonl\.gz$|^v0418_rsweep_bucket_report_honest\.json\.gz$",
        "tools/census/g_vs_input.py --sweep <sweep> ; tools/census/string_sufficiency.py parseback"
        " --sweep <sweep> ; roundtrip_bucket_report.py --score honest   (the neutral ruler on every"
        " generated structure of the release sweep, its smiles_1 read back with no 3D, and the"
        " per-molecule honest bucket: gunzip the three and sweep_two_numbers.two_numbers() returns"
        " the headline from the frozen tree alone)",
    ),
    (
        r"^v0418_rarm2_field2_audit_v04[79]\.(jsonl\.gz|txt)$",
        "tools/v0417/arm2_field2_audit.py --lever OIN_ETA_RETARGET --generator-side   (field 2 of"
        " EVERY golden row, encode-only, fresh process per lever setting; for a generator-side lever"
        " the dead-lever abort is INVERTED: any moved row would mean the lever reaches the encoder)",
    ),
    (
        r"^v0418_rarm2_(diff_full_gate_vs_goldens\.txt|rows_on_full_gate\.tsv\.gz|rows_off_control\.tsv)$",
        "tools/v0417/run_arm2_refreeze.sh on|off ; tools/v0417/arm2_refreeze.py diff|rows   (the FULL"
        " ARM 2 gate, all 425 rows, shipped defaults, after the OIN_ETA_RETARGET promotion. The"
        " fresh rows carry the new HONEST observation column (field 9): sha256 of an independent"
        " XYZToSMILES().convert of the generated xyz, beside the gated field 3)",
    ),
    (
        r"^v0418_rarm2_(lever_fired_but_unseen\.txt|golden_comment_block\.txt|honest_column\.txt)$",
        "tools/v0418/arm2_lever_fired.py ; tools/v0418/arm2_honest_column.py   (did the lever fire"
        " inside the gate -- structure sha vs the previous full run -- and did the gated field see"
        " it; and what the honest column says across all 425 rows: how many rows the gate passes"
        " that the harness's predicate would fail. The comment block is what was spliced into both"
        " goldens)",
    ),
    (
        r"^v0418_rcensus_(attribution_table\.tsv\.gz|attribution_summary\.json|attribution\.txt)$",
        "tools/v0418/run_census_release.sh -> tools/census/attribution_table.py --sweep"
        " <results-v0.4.18-release-sweep> + one path per instrument   (THE CENSUS RE-RUN ON THE"
        " RELEASE SWEEP: 5000 rows, one fault each, UNATTRIBUTED 0; NONE must equal the VERIFIED"
        " count sweep_two_numbers.py printed, or the tool aborts)",
    ),
    (
        r"^v0418_rcensus_control_on_the_record_sweep\.txt$",
        "tools/census/attribution_table.py --write-to <scratch>   (CONTROL, run first: the tool on the"
        " v0.4.14 record sweep reproduces measurements/v0.4.17-census/attribution_table.tsv.gz"
        " BYTE-FOR-BYTE, checked with cmp)",
    ),
    (
        r"^v0418_rcensus_reattribution_diff\.(txt|json)$",
        "tools/census/reattribution_diff.py --old <v0.4.17 re-attribution table> --new <this table>"
        " --movers <structure differs from the v0.4.17 sweep>   (where each fault WENT across both"
        " v0.4.18 promotions, split by mover; the non-mover change count is the noise floor of"
        " attribution)",
    ),
    (
        r"^v0418_rcensus_(parseback_gen|pflags)\.jsonl\.gz$|^v0418_rcensus_attach_class_audit\.json\.gz$",
        "tools/census/string_sufficiency.py parseback --side gen | pflags --charge-probe ;"
        " tools/attach_class_audit.py   (the three SWEEP-DEPENDENT instruments, derived on the"
        " release sweep. C2 e_selfconsistency is REUSED from v0.4.17/ (the encoder is unchanged:"
        " smiles_1 moved on 0 of 5000 in every sweep since); the mirror ruler + twin clusters from"
        " v0.4.17-census/)",
    ),
    (
        r"^v0418_sweep_(RUN\.md|run_config\.json)$",
        "tools/v0418/launch_sweep.sh -> tools/run_sweep.sh <cohort-v0.4.5-5k> <out> 6 300   (THE"
        " v0.4.18 BASELINE SWEEP: commit 2d2708db, SHIPPED DEFAULTS -- the lever block of"
        " run_config.json is EMPTY on purpose, that is the thing under test -- 6 shards 1-BASED,"
        " --mol-timeout 300, BLAS=1. RUN.md is hand-written provenance)",
    ),
    (
        r"^v0418_sweep_(two_numbers\.txt|bucket_report_honest\.md)$",
        "tools/v0418/post_sweep.sh -> roundtrip_bucket_report.py --score honest ;"
        " tools/v0417/sweep_two_numbers.py --sweep <sweep>   (THE v0.4.18 HEADLINE: self-consistent"
        " 4305/5000 = 86.10%, VERIFIED 3883/5000 = 77.66%. The '(record ...)' labels INSIDE"
        " two_numbers.txt are the v0.4.14 record's -- the previous baseline is 82.72 / 74.60)",
    ),
    (
        r"^v0418_sweep_vs_ab\.(txt|json)$",
        "tools/v0418/sweep_vs_ab.py   (CHECKS the two claims the A/B's projection stood on: 3713/3713"
        " non-eta structures byte-identical to the v0.4.17 sweep, 1068/1068 eta structures"
        " byte-identical to the A/B's ON arm, smiles_1 moved on 0 of 5000. Predicted 4302 / 3880,"
        " measured 4305 / 3883: the +3 is seven rows that were hard_fail at the 300 s budget)",
    ),
    (
        r"^v0418_sweep_(g_verdict|parseback)\.jsonl\.gz$|^v0418_sweep_bucket_report_honest\.json\.gz$",
        "tools/census/g_vs_input.py --sweep <sweep> ; tools/census/string_sufficiency.py parseback"
        " --sweep <sweep> ; roundtrip_bucket_report.py --score honest   (the neutral ruler on every"
        " generated structure of the v0.4.18 sweep, its smiles_1 read back with no 3D, and the"
        " per-molecule honest bucket: gunzip the three and sweep_two_numbers.two_numbers() returns"
        " 5000 / 4305 / 3883 from the frozen tree alone)",
    ),
    (
        r"^v0418_arm2_field2_audit_v04[79]\.(jsonl\.gz|txt)$",
        "tools/v0417/arm2_field2_audit.py --lever <the pair> --generator-side   (field 2 of EVERY"
        " golden row, encode-only, fresh process per lever setting: 421 SAME + 4 sentinels, moved 0,"
        " stale 0. For a generator-side lever the dead-lever abort is INVERTED: any moved row would"
        " mean the lever reaches the encoder)",
    ),
    (
        r"^v0418_arm2_(diff_full_gate_vs_goldens\.txt|rows_on_full_gate\.tsv\.gz|rows_off_control\.tsv)$",
        "tools/v0417/run_arm2_refreeze.sh on|off ; tools/v0417/arm2_refreeze.py diff|rows   (the FULL"
        " ARM 2 gate, all 425 rows, shipped defaults, commit 2921aa67: 100/100 + 272/272 compared"
        " rows reproduce, 0 re-frozen. The off control is the two rows SIGKILLed at the hard timeout"
        " -- killed with the levers at 0 too, so not the levers)",
    ),
    (
        r"^v0418_arm2_(lever_fired_but_unseen\.txt|golden_comment_block\.txt)$",
        "tools/v0418/arm2_lever_fired.py   (WHY THE ZERO IS A FINDING ABOUT THE GATE: against the"
        " v0.4.17 full run on pre-lever code, 138 of 198 eta rows built a DIFFERENT structure and 0"
        " of 202 non-eta rows did, yet the gated field 3 differs on 0 of 400 -- the runner re-encodes"
        " through the generator's own bond graph. The comment block is what was spliced into both"
        " goldens; 0 data rows changed, manifests unchanged)",
    ),
    (
        r"^v0418_l2_eta_distance_(identity|audit)\.txt$|^v0418_l2_eta_distance_audit\.jsonl\.gz$",
        "tools/v0418/eta_distance_audit.py [--identity]   (HOW FAR from the metal the generator puts"
        " each ligand, no atom mapping: ligands matched by formula, k shortest M-X distances."
        " Identity control 0.000 A. On the v0.4.17 sweep: eta groups +0.728 A in DETACHED, +0.170 A"
        " in verified passes; sigma -0.313 A. The .jsonl.gz is the run's rows ROUNDED to 3 dp by"
        " tools/v0418/freeze_stage.py to fit the per-file cap)",
    ),
    (
        r"^v0418_l2_(eta_(path|target|exempt|scoped|rescoped|all3)_probe|noneta_identity_probe)\.txt$"
        r"|^v0418_l2_probe[34]_classes\.json$",
        "tools/v0418/eta_path_probe.py --names <file> --arms <...>   (fresh-process SAMPLES, honest"
        " round trip: path = DG vs the rigid placer (ruled out); target/exempt = M1, M2 and both on"
        " 39 detached + 36 verified controls (4 -> 26, controls 36/36); scoped/rescoped = the"
        " exemption's scope; noneta = 25/25 non-eta structures byte-identical; all3 = the third"
        " lever on the A/B's own 47 losses + 45 gains. SAMPLES, NOT VERDICTS: the 36 controls"
        " could not see the 7.5% loss rate the harness A/B then measured)",
    ),
    (
        r"^v0418_l2_eta_ab_report(_on3)?\.(txt|json)$",
        "tools/v0418/run_eta_ab.sh [on3] -> post_eta_ab.sh [on3] -> eta_ab_report.py [--on-arm on3]"
        "   (THE HARNESS A/B over all 1146 eta-bound molecules. M1+M2: self-consistent +195/-29,"
        " VERIFIED +189/-39 -> 86.04% / 77.60% projected. All three levers (_on3): +203/-51,"
        " +179/-59 -- DOMINATED. Noise floor ZERO: smiles_1 identical, OFF == the sweep of record"
        " on 1063/1063 structures. NO DEFAULT CHANGED: both levers ship OFF)",
    ),
    (
        r"^v0418_l2_ab_(off|on|on3)_(g_verdict\.jsonl|bucket_report_honest\.json)\.gz$",
        "tools/roundtrip_bucket_report.py --score honest ; tools/census/g_vs_input.py   (per arm:"
        " the honest bucket and the neutral ruler's verdict on THAT ARM'S structure for every one"
        " of the 1146. With v0.4.18-census/v0418_census_attribution_table.tsv.gz these re-derive"
        " both headline pairs: tools/v0418/freeze_stage.py --verify <this dir>)",
    ),
    (
        r"^v0418_l2_ab_rows\.tsv\.gz$|^v0418_l2_ab_commits\.tsv$",
        "tools/v0418/freeze_stage.py   (MADE at staging time from 3438 structure files + 3438"
        " harness reports too large to freeze: per molecule and arm the sha256 of the generated"
        " structure (the dead-lever check: 798 / 875 differ from OFF), coordination.intact,"
        " elapsed_s, and the sweep of record's sha256 (the noise floor); and the commit each arm"
        " was launched from -- src/ is identical between them)",
    ),
    (
        r"^attribution_(table\.tsv\.gz|summary\.json)$",
        "tools/census/attribution_table.py --gz   (census C4: the join -- one fault per molecule,"
        " first rule that fires, passes included; bucket column reproduces bucket_report_honest;"
        " UNATTRIBUTED printed as a number)",
    ),
    (
        r"^mirror_audit_seed(\d+)_veto_(on|off)\.json$",
        "tools/mirror_audit_donor_fold.py --dataset <cat> --n 250 --seed {0}"
        "   (OIN_FOLD_PARITY_VETO={1})",
    ),
    (
        r"^transition_(fold|veto)\.json$",
        "tools/fold_transition_sim.py --sweep <frozen sweep> --arm {0}",
    ),
    (
        r"^mirror_(movers|cat)_(control_resoff|resonance)\.json$",
        "tools/mirror_audit_donor_fold.py --dataset <{0} cohort> --n {0}"
        "   (OIN_RESONANCE_DONOR_FOLD={1}; the `movers` draw is the MOVER-ENRICHED cohort at"
        " 179/179 coverage -- the `cat` draw holds only 1 mover in 250 = 0.4% and is a"
        " general-population control, NOT evidence about this lever)",
    ),
    (
        r"^veto_outcomes\.json$",
        "tools/veto_outcome_audit.py --sweep <frozen sweep> --dataset <cat> --dataset <photo>"
        "   (which of fold_parity's FIVE outcomes each reverted molecule got: 222/222"
        " vetoed_collapse, 0 no_evidence, over 393/393 movers)",
    ),
    (
        r"^veto_residue_chirality\.json$",
        "tools/veto_residue_chirality.py --outcomes veto_outcomes.json --sweep <frozen sweep>"
        "   (183/222 MIRROR_MATCH on the v0.4.8 corpus; 201/242 re-derived on the v0.4.14"
        " baseline sweep: the round trip built the ENANTIOMER)",
    ),
    (
        r"^(lane1|lane2|both)_pop_(.+)\.json$",
        "tools/run_v0415_arms.sh {0}   (v0.4.15 three-arm A/B over population {1}, via"
        " generator_ab_honest.py -- REAL generation, scored by re-perceiving the WRITTEN XYZ."
        " `both` holds OIN_ATTACH_RETURN on and varies OIN_ACCEPT_STRING_EXACT, i.e. Lane 2's"
        " MARGINAL effect. ⚠ The first run of these arms was VOID: the tool's own"
        " sys.path.insert overrode PYTHONPATH and every arm imported main's oinsmiles, where"
        " neither lever exists, so both sides ran identical code and returned a flawless null."
        " The runner now invokes the tool out of the arm's own checkout and refuses to start"
        " unless oinsmiles resolves there)",
    ),
    (
        r"^selection_pool_probe.*\.json$|^string_exact_pool_probe.*\.json$",
        "tools/selection_pool_probe.py --lever <NAME> --molecules-file <population>"
        "   (counts, per molecule, how many candidate conformers a lever EXAMINED AND REJECTED."
        " A zero-recovery result is only a finding once the rejection count is NON-ZERO;"
        " otherwise it is a wiring failure. ⚠ run it out of the tree under test -- it has the"
        " same sys.path.insert property)",
    ),
    (
        r"^attach_return_preflight\.json$",
        "tools/attach_return_preflight.py --sweep <frozen sweep>"
        "   (v0.4.15 Lane 1 pre-flight. Does the GUARD's predicate -- a claimed coordination"
        " SITE holding nothing, coordinate-only -- fire where `coordination_report` says"
        " DETACHED? Target structural/DETACHED 289/301 = 96.0%; exposure byte_exact/DETACHED"
        " 1/52 = 1.9%, so the charter's 52-molecule exposure over-stated it 52x; control"
        " byte_exact/INTACT 250/250 SITES_HELD, which is what makes the 1.9% a reading rather"
        " than a broken claim count)",
    ),
    (
        r"^resonance_transition\.json$",
        "tools/resonance_transition_sim.py --sweep <frozen sweep> --baseline-byte-exact 3794"
        "   (OFFLINE re-score. Reported +78/0 losses; SUPERSEDED by reso_full_ab.json --"
        " an offline re-score cannot express a loss)",
    ),
    (
        r"^resonance_key_invariance\.json$",
        "tools/fold_key_invariance.py --sweep <frozen sweep> --lever OIN_RESONANCE_DONOR_FOLD"
        " --holding OIN_CANONICAL_DONOR_FOLD"
        "   (9669 compared, 228 moved, 0 keys changed. ⚠ bounds ACCEPTANCE, not embedding)",
    ),
    (
        r"^arm2_refreeze_rows\.tsv$",
        "tools/gate_arm2_roundtrip_one.py over the 13 v0.4.14 golden movers -- the rows spliced"
        " into gate_v047/v049_arm2_golden.tsv (fields 1-6 only; field 7 is the --band column)",
    ),
    (
        r"^gain_generator_ab\.tsv$",
        "SUPERSEDED. tools/gate_arm2_roundtrip_one.py A/B over 20 v0.4.14 gains -- scored with"
        " arm2's CIRCULAR predicate (get_oin_string(result.mol)) and therefore wrong. Kept as the"
        " artifact behind a documented error; the correct measurement is reso_full_ab.json",
    ),
    (
        r"^per_molecule_extract\.tsv$",
        "derived from results-v0.4.14-sweep/individual_reports (N=5000): molecule, status,"
        " tier_passed, bucket under BOTH scored and honest, subclass, and the NESTED"
        " metrics.elapsed_s. Re-derives the authoritative table and the runtime percentiles"
        " without the 261 MB run",
    ),
    (
        r"^RUN\.md$",
        "hand-written provenance for the v0.4.14 baseline sweep: 6 shards 1-BASED,"
        " --mol-timeout 300, shipped lever defaults, BLAS threads capped to 1, and why"
        " systemd's OOMPolicy could not be used on a --scope",
    ),
    (
        r"^bucket_report_PASS1_authoritative\.md$",
        "tools/roundtrip_bucket_report.py --results-dir <sweep> --score honest, frozen"
        "   (THE v0.4.14 ABSOLUTE BASELINE: byte_exact 3858/5000 = 77.16%, gap 22.84)",
    ),
    (
        r"^reso_sample_definitions\.json$",
        "hand-assembled index: seed + membership for every sample behind a v0.4.14 figure,"
        " plus the class splits. Exists because the scratchpad they were drawn in lives under"
        " /tmp and is DELETED between sessions -- a rate without its sample is not reproducible",
    ),
    (
        r"^v0413_gains_ab\.json$",
        "tools/generator_ab_honest.py --lever OIN_CANONICAL_DONOR_FOLD over 25 of v0.4.13's 171"
        " claimed gains (seed 17), resonance held OFF in both arms"
        "   (22/25 = 88% real => ~150 gains; with ~30 losses, v0.4.13's true net is ~+2.42 pts)",
    ),
    (
        r"^reso_movers_exact\.json$",
        "tools/lever_string_movers.py --lever OIN_RESONANCE_DONOR_FOLD --holding"
        " OIN_CANONICAL_DONOR_FOLD --holding OIN_FOLD_PARITY_VETO"
        "   (93 of 5000 move encode(input) -- the coordinate-derived affected population)",
    ),
    (
        r"^reso_full_ab\.json$",
        "tools/generator_ab_honest.py --lever OIN_RESONANCE_DONOR_FOLD over all 182 affected"
        "   (THE v0.4.14 HEADLINE: 78 gains, 7 losses, net +71 = +1.42 pts, n=182 of 182)",
    ),
    (
        r"^v0413_atrisk_ab\.json$",
        "tools/generator_ab_honest.py --lever OIN_CANONICAL_DONOR_FOLD over 40 of v0.4.13's 197"
        " at-risk molecules (seed 13)"
        "   (6 losses = 15% => ~30 over the population => v0.4.13's true net ~+2.82, NOT +3.42)",
    ),
    (
        r"^(atrisk_generator_ab|hekfel_honest_ab|generator_ab_honest)\.json$",
        "tools/generator_ab_honest.py --lever OIN_RESONANCE_DONOR_FOLD (sampling pass,"
        " superseded by reso_full_ab.json)",
    ),
    (
        r"^ab_.*\.json$",
        "tools/ab_accept_scored.py --cohort <cohort> --lever <lever> --timeout 150 --hard-cap 240",
    ),
    (r"^cohort_.*\.json$", "cohort manifest, built from the frozen sweep"),
    # v0.4.13's instruments.
    (
        r"^mirror_arm([AB])_(\w+)\.json$",
        "tools/mirror_audit_donor_fold.py --dataset <cohort-v0.4.5-5k> --n 250 --seed 7"
        "   (arm {0}: {1}; the noveto arm sets OIN_FOLD_PARITY_VETO=0) -- mixed cat+photo draw,"
        " reads 33 collapses -> 0",
    ),
    (
        r"^mirror_cat_(\w+)\.json$",
        "tools/mirror_audit_donor_fold.py --dataset <cat> --n 250 --seed 7"
        "   ({0}; the noveto arm sets OIN_FOLD_PARITY_VETO=0) -- CAT-ONLY draw, reproduces"
        " v0.4.12's published 19 -> 0 with achiral unmoved at 157",
    ),
    (
        r"^fold_transition_(fold|veto)\.json$",
        "tools/fold_transition_sim.py --sweep <frozen sweep> --arm {0}"
        "   --dataset <cat> --dataset <photo>   (ABSOLUTE roots: the relative default"
        " silently excludes every mover from a worktree)",
    ),
    (
        r"^fold_key_invariance\.json$",
        "tools/fold_key_invariance.py --sweep <frozen sweep>"
        "   (does the fold ever change the round-trip KEY? 0 => generator-neutral"
        " => an offline re-score is exact and no sweep is owed)",
    ),
    (
        r"^attach_class_audit\.json$",
        "tools/attach_class_audit.py --results-dir <frozen sweep>"
        "   (MEDZUR/GAVSED split, with the byte_exact control arm)",
    ),
    (
        r"^prefilter_.*\.json$",
        "tools/prefilter_prevalence.py --xyz <input> | --cohort <dir>"
        "   (OIN_PREFILTER_ADVISORY two-arm; needs a QUIET machine — it reports a latency cost)",
    ),
    (r"^bucket_report.*\.md$", "tools/roundtrip_bucket_report.py --results-dir <dir>"),
    (r"^FROZEN\.md$|^SOURCE$", "hand-written provenance record for a frozen sweep"),
    (r"^triage_hard_fails\.md$", "tools/triage_hard_fails.py"),
    (r"^CASE_REGISTRY\.md$|^.*TREND\.tsv$", "tools/rebuild_summary.py / milestone_report.py"),
    # --- v0.4.19 RELEASE ------------------------------------------------------------------------
    (
        r"^v0419_rsweep_run_config\.json$",
        "tools/v0419/launch_candidate_sweep.sh -> tools/run_sweep.sh <cohort-v0.4.5-5k> <out> 6 300"
        "   (THE v0.4.19 RELEASE SWEEP: the full 5,000 with OIN_H_FAITHFUL + OIN_RC1_PROPAGATE +"
        " OIN_CAP_IGNORES_METAL ON via the unit's environment -- the lever block of run_config.json"
        " is EXACTLY those three -- at commit 1f61681c, promoted to default-ON at 0809c7a9 (a"
        " byte-identical read). 6 shards 1-BASED, --mol-timeout 300, BLAS=1. Launched twice: the"
        " first run was aborted at 528/5000 because src/ was edited in its worktree)",
    ),
    (
        r"^v0419_rsweep_(two_numbers\.txt|bucket_report_honest\.md)$",
        "SWEEP_TAG=v0.4.19-candidate tools/v0419/post_candidate_sweep.sh -> roundtrip_bucket_report.py"
        " --score honest ; tools/v0417/sweep_two_numbers.py --sweep <sweep>   (THE v0.4.19 HEADLINE:"
        " 4,389 / 5,000 = 87.78% self-consistent, 3,956 / 5,000 = 79.12% VERIFIED; predicted 4,378 /"
        " 3,951. The '(record ...)' labels INSIDE two_numbers.txt are the v0.4.14 record's)",
    ),
    (
        r"^v0419_rsweep_vs_projection\.(txt|json)$|^v0419_rsweep_unchanged_string_movers\.json$",
        "tools/v0419/candidate_vs_projection.py --sweep <sweep>   (the sweep against the projection"
        " ROW BY ROW: 40 differences, 5 on the 300 s boundary, every one an unchanged-string row --"
        " OIN_CAP_IGNORES_METAL reaches the GENERATOR through its re-encode of candidate conformers:"
        " 44 of 4,891 unchanged-string rows built a different structure, listed in"
        " unchanged_string_movers.json with both runs' honest buckets)",
    ),
    (
        r"^v0419_rsweep_(g_verdict|parseback|parseback_gen)\.jsonl\.gz$|^v0419_rsweep_bucket_report_honest\.json\.gz$",
        "OIN_H_FAITHFUL=1 OIN_RC1_PROPAGATE=1 OIN_CAP_IGNORES_METAL=1 tools/census/string_sufficiency.py"
        " parseback [--side gen] --sweep <sweep> ; tools/census/g_vs_input.py --sweep <sweep> ;"
        " roundtrip_bucket_report.py --score honest   (the neutral ruler on every generated structure,"
        " the strings read back UNDER THE LEVERS -- the parse-back reader is what OIN_H_FAITHFUL"
        " gates -- and the per-molecule honest bucket: gunzip the three and"
        " sweep_two_numbers.two_numbers() returns the headline from the frozen tree alone)",
    ),
    (
        r"^v0419_rsweep_parseback_shippedreader\.jsonl\.gz$",
        "OIN_H_FAITHFUL=0 OIN_RC1_PROPAGATE=0 OIN_CAP_IGNORES_METAL=0 tools/census/string_sufficiency.py"
        " parseback --sweep <sweep>   (the SAME strings read by the pre-promotion reader: what a"
        " build without OIN_H_FAITHFUL makes of a string written with it -- the writer + reader"
        " coupling, kept beside the promoted reader's verdict, never hidden)",
    ),
    (
        r"^v0419_rarm2_arm1_(on|off_control)\.tsv$",
        "tools/gate_arm1_encode.py --expect-n 65, shipped defaults / the three levers at 0   (ARM 1"
        " at 65 fixtures: the control reproduces the v0.4.17 golden on all 62 rows; the promotion"
        " moves 5, each attributed by one lever at a time in the golden's header)",
    ),
    (
        r"^v0419_rarm2_field2_audit_v04[79]\.(jsonl\.gz|txt)$",
        "tools/v0417/arm2_field2_audit.py --lever OIN_H_FAITHFUL,OIN_RC1_PROPAGATE,OIN_CAP_IGNORES_METAL"
        " --golden tools/gate_v04{7,9}_arm2_golden.tsv --cohort-dir <cohort>   (/refreeze-goldens"
        " step 1, encode-only, every row: v047 100 SAME; v049 295 SAME + 26 MOVED_BY_LEVER + 4"
        " sentinels, 0 STALE. v047's '0 fired' abort is genuine: none of its 100 molecules is in"
        " the 109-row changed-string set)",
    ),
    (
        r"^v0419_rarm2_(diff_full_gate_vs_goldens\.txt|splice_dry_run\.txt|rows_on_full_gate\.tsv\.gz|rows_off_control\.tsv|rows_verify_real_gate\.tsv|golden_comment_block\.txt)$",
        "REFREEZE_LEVER=OIN_H_FAITHFUL,OIN_RC1_PROPAGATE,OIN_CAP_IGNORES_METAL REFREEZE_TAG=v0.4.19"
        " tools/v0417/run_arm2_refreeze.sh on|off|verify ; tools/v0417/arm2_refreeze.py diff|splice"
        " --accept-det MAHTOE_comp_0|rows   (the FULL gate, shipped defaults, 6 shards: v047 100/100"
        " reproduce; v049 owes 26 rows, all field 2, exactly the audit's; the control brings every"
        " one back -> 26 LEVER, 0 STALE; MAHTOE (NO_STRUCTURE_DET, now built and gated) decided by"
        " hand; the real gate reads the re-frozen goldens)",
    ),
    (
        r"^v0419_rcensus_(attribution_table\.tsv\.gz|attribution_summary\.json|attribution\.txt)$",
        "tools/v0419/run_census_release.sh -> tools/census/attribution_table.py --sweep"
        " <results-v0.4.19-candidate-sweep> + one path per instrument   (THE CENSUS ON THE v0.4.19"
        " RELEASE SWEEP; C2 = the lane's fix3 canonicality audit, since the encoder changed)",
    ),
    (
        r"^v0419_rcensus_control_on_the_record_sweep\.txt$",
        "tools/census/attribution_table.py --write-to <control>   (the parameterised tool must"
        " reproduce the frozen v0.4.17 census table byte-for-byte before this release's is believed)",
    ),
    (
        r"^v0419_rcensus_reattribution_diff\.(txt|json)$|^v0419_rcensus_movers_vs_v0418\.txt$",
        "tools/census/reattribution_diff.py --old <v0.4.18 release census> --new <this> --movers"
        " movers_vs_v0418.txt   (fault transitions against the v0.4.18 release census, split by"
        " whether the generated STRUCTURE differs from the v0.4.18 release sweep's)",
    ),
    (
        r"^v0419_rcensus_(parseback_gen|pflags)\.jsonl\.gz$|^v0419_rcensus_attach_class_audit\.json\.gz$",
        "tools/census/string_sufficiency.py parseback --side gen / pflags --charge-probe ;"
        " tools/attach_class_audit.py   (the per-sweep instruments the attribution reads)",
    ),
    # --- v0.4.19 E2 lane (branch research/v0419-e2; every run records its commit in commits.tsv)
    (
        r"^v0419_e2_(fragile105\.txt|esc105_e2b_final\.jsonl\.gz|e2b_final_fragility\.jsonl\.gz)$",
        "tools/census/e_selfconsistency.py --only <the 105 census E2_P_FRAGILE rows> under"
        " OIN_N_VALENCE_2=1 OIN_CANONICAL_RESONANCE=1 (cut 1) -> tools/v0419/e2_fragility_report.py"
        " --names fragile105.txt   (88 stable, renumbering fixed 94/100, 0 newly fragile)",
    ),
    (
        r"^v0419_e2_n2_valence_scan\.json\.gz$",
        "an ad hoc scan during the XIVMEX bisection (scratchpad script, not committed): per molecule,"
        " ligands carrying a two-coordinate nitrogen and the AC2BO valence-product size with and"
        " without the valence-2 option",
    ),
    (
        r"^v0419_e2_(e_selfconsistency_e2b\.jsonl\.gz|esc_e2b_summary\.txt|canonicality_audit_e2b\.(txt|json)"
        r"|e2b_cohort_fragility\.txt|changed_set_e2b\.txt|became_fragile_e2b\.txt|movers_e2b\.txt)$",
        "CUT 1 (valence 2 appended for every two-coordinate N, + RES) on the whole cohort:"
        " tools/census/e_selfconsistency.py --sweep results-v0.4.19-candidate-sweep --cpu 8 ->"
        " OIN_AB_BASELINE=v0.4.19 tools/v0419/serializer_ab_report.py audit/changed ;"
        " tools/v0419/e2_fragility_report.py vs the fix3 run   (343 strings moved, 167 VERIFIED"
        " passes; fragile fixed 215 / broken 53, noise +5/-44 -- REFUTED)",
    ),
    (
        r"^v0419_e2_changed_single_e2b(\.jsonl\.gz|_summary\.txt)$",
        "tools/v0419/lever_parseback.py --names changed_set_e2b.txt --arms shipped,n2,res,e2b"
        "   (single-lever attribution of cut 1's 343: N2 alone 292, RES alone 44)",
    ),
    (
        r"^v0419_e2_(e_selfconsistency_e2d_movers\.jsonl\.gz|esc_e2d_movers_summary\.txt"
        r"|e2d_movers_vs_(shipped|e2b)\.txt)$",
        "CUT 2 (a second pass inside AC2BO, 14825a5f) on the 443 movers: e_selfconsistency.py --only"
        " <movers_e2b.txt> -> e2_fragility_report.py vs fix3 and vs cut 1   (251 still moved --"
        " REFUTED)",
    ),
    (
        r"^v0419_e2_movers_single_e2e(\.jsonl\.gz|_summary\.txt)$",
        "CUT 3 (a ladder-level fallback on a bond-order guess, 2a6a4b6a): lever_parseback.py --names"
        " movers_e2b.txt --arms shipped,n2,e2b   (177 of 443 moved, 59 verified passes -- REFUTED)",
    ),
    (
        r"^v0419_e2_(movers_single_e2f(\.jsonl\.gz|_summary\.txt)|n2_scoped_movers\.json)$",
        "CUT 4 (the option scoped to pyrrolide nitrogens, 7dfaba98): lever_parseback.py --names"
        " movers_e2b.txt --arms shipped,n2,res,e2b   (N2 moves 84, 66 already fragile)",
    ),
    (
        r"^v0419_e2_(e_selfconsistency_(e2f|res)\.jsonl\.gz|esc_(e2f|res)_summary\.txt"
        r"|(e2f|res)_cohort_fragility\.txt|changed_(e2f|res)\.txt|changed_set_(e2f|res)\.txt)$",
        "the whole cohort at 7dfaba98 under cut 4 + RES (e2f) and under RES alone (res):"
        " e_selfconsistency.py --cpu 8 -> e2_fragility_report.py vs fix3 ; serializer_ab_report.py"
        " changed   (e2f: 136 moved, fragile fixed 120 / broken 5; res: 56 moved, 81 / 1)",
    ),
    (
        r"^v0419_e2_rescore_(e2f|res)\.jsonl\.gz$",
        "tools/honest_rescore.py --results-dir results-v0.4.19-candidate-sweep under each arm's"
        " levers   (the offline column: unchanged-string rows, structure held fixed -- BLIND to"
        " generator reach, which a perception lever has)",
    ),
    (
        r"^v0419_e2_(serializer_ab_report_(e2f|res)\.(txt|json)|ab_(e2f|res)_on_(bucket_report_honest\.json"
        r"|g_verdict\.jsonl|parseback\.jsonl|parseback_shippedreader\.jsonl)\.gz|commits\.tsv)$",
        "tools/v0419/run_serializer_ab.sh <changed_set> e2f|res -> post_serializer_ab.sh e2f|res"
        " (OFF = the v0.4.19 release sweep's own rows)   (e2f: live self +69/-7, VERIFIED +66/-7,"
        " offline +17/-2 and +11/-1, projection 89.32% / 80.50%; res: +23/-3, +20/-3, +18/-1,"
        " +12/-0, 88.52% / 79.70%. NO DEFAULT CHANGED: both levers ship OFF)",
    ),
]


def repo_root() -> Path:
    """The MAIN checkout, resolved from anywhere -- including inside a linked worktree.

    ``--git-common-dir`` is the shared ``.git`` directory; from a worktree it is an absolute
    path into the main checkout, and from the main checkout it is the relative ``.git``. Its
    parent is the main working tree in both cases.
    """
    common = subprocess.check_output(["git", "rev-parse", "--git-common-dir"], text=True).strip()
    return Path(common).resolve().parent


def in_worktree() -> bool:
    git_dir = subprocess.check_output(["git", "rev-parse", "--git-dir"], text=True).strip()
    return "worktrees" in Path(git_dir).resolve().parts


def _wanted(name: str) -> bool:
    if any(fnmatch.fnmatch(name, d) for d in DENY):
        return False
    return any(fnmatch.fnmatch(name, a) for a in ALLOW)


def collect(src: Path) -> list[tuple[int, Path]]:
    """Allowlisted files under ``src``, each within the per-file cap."""
    out: list[tuple[int, Path]] = []
    for dirpath, dirnames, files in os.walk(src):
        # ⚠ `INVALID*` is pruned by PREFIX, not by exact name, and it is not cosmetic. v0.4.15 kept
        # its void arms beside the real ones as `INVALID-wrong-tree/`, and those files carry the
        # SAME BASENAMES -- so without this the harvester copies both and one silently overwrites
        # the other, publishing a measurement that was already known to be wrong. Caught in the
        # dry run only because the file list was read line by line and showed
        # `lane1_pop_L1_target_site_lost.json` twice at two different sizes.
        dirnames[:] = [
            d for d in dirnames if d not in PRUNE_DIRS and not d.upper().startswith("INVALID")
        ]
        for f in files:
            if not _wanted(f):
                continue
            p = Path(dirpath) / f
            try:
                size = p.stat().st_size
            except OSError:
                continue
            if size <= PER_FILE_CAP:
                out.append((size, p))
    return sorted(out, key=lambda t: str(t[1]))


def provenance_for(name: str) -> str:
    for pattern, template in PROVENANCE:
        m = re.match(pattern, name)
        if m:
            try:
                return template.format(*m.groups())
            except (IndexError, KeyError):
                return template
    return "UNKNOWN"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def scrub(text: str, root: Path) -> tuple[str, bool]:
    """Rewrite local absolute paths to repo-relative form. Returns ``(text, still_leaks)``.

    ⚠ This REWRITES rather than refuses, and the distinction matters. The tools that produced
    these artifacts wrote absolute paths into their own headers, so a refuse-on-`/home/` rule
    rejects exactly the files worth keeping -- measured: it dropped every `bucket_report_*.md`
    and the v0.4.8 `SOURCE`, i.e. the frozen baseline this whole tree exists to preserve, while
    happily keeping the `ARM_A_DONE` sentinels. Losing the data to protect the path is the wrong
    trade when the path is trivially removable.

    ``<repo>/`` is stripped (the reference stays meaningful, since the tree lives in the repo)
    and any surviving home directory is masked. Anything still matching afterwards is a genuine
    leak and the caller drops the file -- so the guarantee is kept, just not by throwing the
    baby out.
    """
    out = text.replace(str(root) + "/", "").replace(str(root), ".")
    out = re.sub(r"/home/[^/\s]+/", "<HOME>/", out)
    out = re.sub(r"/tmp/claude-\d+/[^\s\"')]*", "<SCRATCH>", out)
    still = "/home/" in out or "/tmp/claude" in out
    return out, still


#: Renames applied on the way in. The v0.4.12 scratchpad called the same kind of run
#: ``audit_baseline_seed7`` and ``audit_base_seed11``, and "baseline" is ambiguous anyway --
#: mirror_audit_donor_fold.py runs BOTH fold arms internally, so what differed between those
#: runs was the ambient veto, not the fold.
RENAME = {
    "audit_baseline_seed7.json": "mirror_audit_seed07_veto_off.json",
    "audit_veto_seed7.json": "mirror_audit_seed07_veto_on.json",
    "audit_base_seed11.json": "mirror_audit_seed11_veto_off.json",
    "audit_veto_seed11.json": "mirror_audit_seed11_veto_on.json",
    "transition_fold_baseline.json": "transition_fold.json",
    "ab_pilot_eta.json": "ab_eta_accept_stale_cohort.json",
    "ab_pilot2_eta.json": "ab_eta_accept_realpop.json",
    "cohort_pilot.json": "cohort_pilot_stale.json",
}

#: Scratch intermediates: real files, but not artifacts anyone will diff.
SKIP_NAMES = {"movers6.json", "movers40.json"}


def write_release(
    dest: Path, release: str, picks: list[tuple[int, Path, str]], dry: bool, root: Path
):
    rel_dir = dest / release
    rows = []
    for size, src, newname in picks:
        rows.append((newname, size, src))
        if dry:
            continue
        rel_dir.mkdir(parents=True, exist_ok=True)
        try:
            text = src.read_text(errors="strict")
        except (OSError, UnicodeDecodeError):
            shutil.copy2(src, rel_dir / newname)  # binary/unreadable: copy verbatim
            continue
        cleaned, _still = scrub(text, root)
        (rel_dir / newname).write_text(cleaned)

    if dry:
        return rows

    lines = [
        f"# `measurements/{release}` — frozen comparison artifacts",
        "",
        "Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.",
        "",
        "| file | bytes | sha256 | produced by |",
        "|---|---:|---|---|",
    ]
    for name, size, src in sorted(rows):
        lines.append(f"| `{name}` | {size} | `{sha256(rel_dir / name)}` | {provenance_for(name)} |")
    # ⚠ The generated index is itself published, so it gets the same scrubbing as the harvested
    # files. Missing this leaked the full scratchpad path (`/tmp/claude-<uid>/-home-<user>-...`)
    # into a public repo on the first run -- the guard has to cover what the guard writes.
    src_lines, _ = scrub("\n".join(f"{n}  <-  {s}" for n, _sz, s in sorted(rows)), root)
    lines += ["", "Source paths at harvest time:", "", "```", src_lines, "```"]
    (rel_dir / "README.md").write_text("\n".join(lines) + "\n")
    return rows


def infer_release(dirname: str) -> str:
    """``results-v0.4.8-honest`` -> ``v0.4.8-honest``; ``results-capstone-v042`` -> ``capstone-v042``.

    The ``results-`` prefix is always stripped -- it is noise inside a tree whose every entry is
    a result. Names that carry no version at all keep their remaining text rather than being
    forced into a version-shaped slot they do not fit.
    """
    return re.sub(r"^results-", "", dirname)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", help="target release folder, e.g. v0.4.12")
    ap.add_argument(
        "--from", dest="sources", action="append", default=[], help="source directory (repeatable)"
    )
    ap.add_argument(
        "--backfill",
        action="store_true",
        help="harvest every tmCAT-tmPHOTO_xyz_dataset/results-* dir",
    )
    ap.add_argument("--dry-run", action="store_true", help="print the selection, write nothing")
    ap.add_argument(
        "--link", action="store_true", help="create .measurements-main in this worktree and exit"
    )
    ap.add_argument("--force", action="store_true", help="proceed past the total-size cap")
    args = ap.parse_args()

    root = repo_root()
    dest = root / "measurements"
    print(f"main checkout : {root}")
    print(f"destination   : {dest}")
    if in_worktree():
        print("⚠ running inside a WORKTREE — writing to the main checkout above, not here.")

    if args.link:
        link = Path.cwd() / ".measurements-main"
        if link.is_symlink() or link.exists():
            print(f"already present: {link}")
            return 0
        link.symlink_to(dest)
        print(f"linked {link} -> {dest}")
        return 0

    jobs: list[tuple[str, Path]] = []
    if args.backfill:
        ds = root / "tmCAT-tmPHOTO_xyz_dataset"
        for d in sorted(ds.glob("results-*")):
            if d.is_dir():
                jobs.append((infer_release(d.name), d))
    if args.sources:
        if not args.release:
            ap.error("--from requires --release")
        for s in args.sources:
            jobs.append((args.release, Path(s).expanduser()))
    if not jobs:
        ap.error("nothing to do: pass --backfill and/or --release with --from")

    grand_total, grand_files, unresolved = 0, 0, []
    by_release: dict[str, list] = {}
    for release, src in jobs:
        if not src.is_dir():
            print(f"  ⚠ missing source, skipped: {src}")
            continue
        picks = []
        for size, p in collect(src):
            if p.name in SKIP_NAMES:
                continue
            try:
                text = p.read_text(errors="strict")
            except (OSError, UnicodeDecodeError):
                picks.append((size, p, RENAME.get(p.name, p.name)))
                continue
            _cleaned, still = scrub(text, root)
            if still:
                # Scrubbing could not make it safe -- that IS a refusal, and a rare one.
                print(f"  🔴 REFUSED (local path survives scrubbing, tree is public): {p}")
                continue
            picks.append((size, p, RENAME.get(p.name, p.name)))
        if not picks:
            continue
        total = sum(s for s, _p, _n in picks)
        grand_total += total
        grand_files += len(picks)
        print(f"\n{release}: {len(picks)} files, {total / 1024:.0f} KB   <- {src}")
        for size, p, newname in picks:
            tag = "" if provenance_for(newname) != "UNKNOWN" else "   [provenance UNKNOWN]"
            rename = f"  (was {p.name})" if newname != p.name else ""
            print(f"    {size:>8}  {newname}{rename}{tag}")
            if provenance_for(newname) == "UNKNOWN":
                unresolved.append(f"{release}/{newname}")
        # ⚠ ACCUMULATE, never write per --from. write_release() rebuilds README.md from the batch
        # it is handed and OVERWRITES, so calling it once per source left the index describing only
        # the LAST source: v0.4.14 harvested 16 files and its README documented 2. The files were
        # all present, which is what makes this shape dangerous -- the data looks complete and the
        # provenance record, the entire point of the tree, silently is not.
        by_release.setdefault(release, []).extend(picks)

    for release, picks in by_release.items():
        write_release(dest, release, picks, args.dry_run, root)

    print(f"\nTOTAL: {grand_files} files, {grand_total / 1024:.0f} KB")
    if unresolved:
        print(
            f"⚠ {len(unresolved)} file(s) with UNKNOWN provenance — recorded as such, not guessed."
        )
    if grand_total > TOTAL_CAP and not args.force:
        print(
            f"🔴 REFUSED: {grand_total / 1024 / 1024:.1f} MB exceeds the {TOTAL_CAP // 1024 // 1024} MB cap."
        )
        print("   measurements/ is a comparison tree in a PUBLIC repo, not a data dump.")
        print("   Narrow the selection, or pass --force if this is genuinely warranted.")
        return 1
    if args.dry_run:
        print("(dry run — nothing written)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
