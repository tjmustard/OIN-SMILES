#!/usr/bin/env python3
"""Group v0.4.0 round-trip failures into a curated v0.4.1 accuracy backlog.

Consumes the ``case_registry.json`` written by ``tools/classify_failures.py`` plus
the raw ``individual_reports/`` in the same output dir, and produces three files in
that dir:

* ``V0.4.1_ACCURACY_BACKLOG.md`` -- failures bucketed into tiers
  (A real accuracy defects / B generator robustness / C wontfix / D noise), each
  class with its molecule list, a representative example, and a copy-paste repro.
* ``V0.4.1_TREND.tsv`` -- one row appended per run so we can watch the class mix
  stabilize over the life of the accumulator.
* ``.v041_state.json`` -- last snapshot, used to print the per-run delta.

Strictly read-only on the round-trip harness's outputs: it reads
``case_registry.json`` and ``individual_reports/`` and writes only its own three
files. Safe to run repeatedly alongside a live ``--continue`` accumulator.

Usage:
    python tools/group_v041_backlog.py --output-dir <results-dir>
"""

import argparse
import glob
import json
import os
from collections import Counter, defaultdict
from datetime import datetime

# --- taxonomy -> tier -----------------------------------------------------
# Tier A: deterministic encode / bond-order / stereo defects. Quick-insensitive,
#         these are the real v0.4.1 accuracy targets.
# Tier B: generator robustness. --quick FF-only inflates these; verify at full
#         fidelity before treating any as an accuracy defect.
# Tier C: notation limits we've decided not to fix (document instead).
# Tier D: harness / provenance noise -- re-run, don't debug.
TIER_A = {
    "donor_H_atom_count",
    "garbled_aromatic",
    "macrocycle_perception",
    "atom_stereo",
    "EZ_bond_stereo",
    "winding_flip",
    "geometry_or_fragment_change",
    "H_on_terminal_oxo_imido",
    "eta_diene_localization",
    "frag_vector_IndexError",
    "geometry_NON",
    "kekulize_encode_crash",
    "xyz2mol_none_crash",
    "encode_crash_other",
    "rmsd_mapping_failed",
    "high_rmsd",
    "string_mismatch_other",
}
TIER_B = {"no_conformers", "timeout", "gen_exception_other"}
TIER_C = {"carborane_unsupported"}
TIER_D = {"interrupted", "junk_row", "stale_harness_TypeError", "rmsd_sentinel"}

TIER_TITLE = {
    "A": "Tier A · Real accuracy defects (v0.4.1 targets)",
    "B": "Tier B · Generator robustness (triage — quick-sensitive)",
    "C": "Tier C · Wontfix / notation limit",
    "D": "Tier D · Harness / provenance noise (re-run, don't debug)",
}
TIER_BLURB = {
    "A": "Deterministic encode / bond-order / stereo bugs. Quick-insensitive — these "
    "are the disjoint fix sessions for v0.4.1.",
    "B": "`--quick` FF-only inflates these; re-run a sample at full fidelity before "
    "treating any as a real accuracy target.",
    "C": "3c2e cages / exotic main-group notation. Document in `docs/KNOWN_LIMITATIONS.md`.",
    "D": "Interrupted rows, junk, legacy sentinels. Not defects — re-run to clear.",
}


def tier_of(cls):
    if cls in TIER_A:
        return "A"
    if cls in TIER_B:
        return "B"
    if cls in TIER_C:
        return "C"
    if cls in TIER_D:
        return "D"
    return "A"  # unknown class -> treat conservatively as an accuracy target


def load_reports(indiv_dir):
    """molecule -> raw report dict (for smiles/error detail)."""
    reps = {}
    for fp in glob.glob(os.path.join(indiv_dir, "*.json")):
        try:
            with open(fp) as f:
                rep = json.load(f)
        except Exception:
            continue
        mol = rep.get("molecule")
        if mol:
            reps[mol] = rep
    return reps


def elapsed_of(mol, reps):
    """metrics.elapsed_s for a molecule (harness wall-clock), or None if absent.

    Added to the report schema alongside the --quick 60s->30s change; lets us tell a
    genuinely pathological embed (tens of seconds up against the timeout) from a fast
    quick-mode artifact in the robustness tier.
    """
    v = (reps.get(mol, {}).get("metrics") or {}).get("elapsed_s")
    return float(v) if isinstance(v, (int, float)) else None


def representative(cls, mols, reps):
    """A short evidence block for the first molecule in a class."""
    mol = sorted(mols)[0]
    rep = reps.get(mol, {})
    s1 = (rep.get("smiles_1") or "").strip()
    s2 = (rep.get("smiles_2") or "").strip()
    err = (rep.get("error") or "").split("\n")[0].strip()
    el = elapsed_of(mol, reps)
    tag = f" ({el:.0f}s)" if el is not None else ""
    lines = [f"  - representative: `{mol}`{tag}"]
    if s1 and s2 and s1 != s2:
        lines.append(f"    - exp: `{s1[:160]}`")
        lines.append(f"    - got: `{s2[:160]}`")
    elif err:
        lines.append(f"    - error: `{err[:180]}`")
    lines.append(
        f"    - repro: `python tools/test_dataset_roundtrip.py --only {mol} "
        f"--quick --output-dir /tmp/rt-{mol}`"
    )
    return lines


def mol_list(cls, mols, reps):
    """Comma list of molecules; Tier B (robustness) annotated with wall-clock."""
    if tier_of(cls) != "B":
        return ", ".join(f"`{m}`" for m in mols)
    items = []
    for m in sorted(mols, key=lambda x: -(elapsed_of(x, reps) or 0)):
        el = elapsed_of(m, reps)
        items.append(f"`{m}`" + (f" ({el:.0f}s)" if el is not None else ""))
    return ", ".join(items)


def build_backlog(rows, reps, output_dir):
    failures = [r for r in rows if r["class"] != "success"]
    n_reports = len(rows)
    n_success = n_reports - len(failures)

    by_class = defaultdict(list)
    for r in failures:
        by_class[r["class"]].append(r["molecule"])

    tiers = defaultdict(list)  # tier -> [(cls, [mols])]
    for cls, mols in by_class.items():
        tiers[tier_of(cls)].append((cls, sorted(mols)))
    for t in tiers:
        tiers[t].sort(key=lambda cm: (-len(cm[1]), cm[0]))

    pass_pct = 100 * n_success / max(n_reports, 1)
    md = [
        "# v0.4.1 Accuracy Backlog — tmCAT/tmPHOTO round-trip (v0.4.0 stream)",
        "",
        f"Generated {datetime.now().isoformat(timespec='seconds')} "
        "by `tools/group_v041_backlog.py` from the live `case_registry.json`.",
        "",
        f"**{n_reports} molecules covered · {n_success} pass / {len(failures)} fail "
        f"· {pass_pct:.1f}% pass.**",
        "",
        "Failure classes bucketed by whether they are real accuracy defects, generator "
        "robustness (quick-sensitive), notation limits, or harness noise. This file is "
        "regenerated every watch tick; when the class mix stabilizes it becomes the "
        "v0.4.1 fix-wave scoping doc.",
        "",
        "| tier | classes | molecules |",
        "|---|---|---|",
    ]
    for t in ("A", "B", "C", "D"):
        cls_list = tiers.get(t, [])
        n_mol = sum(len(m) for _, m in cls_list)
        md.append(f"| {t} | {len(cls_list)} | {n_mol} |")
    md.append("")

    for t in ("A", "B", "C", "D"):
        cls_list = tiers.get(t, [])
        if not cls_list:
            continue
        md += [f"## {TIER_TITLE[t]}", "", f"_{TIER_BLURB[t]}_", ""]
        for cls, mols in cls_list:
            md.append(f"### `{cls}` — {len(mols)} molecule(s)")
            md.append("")
            md += representative(cls, mols, reps)
            md.append(f"  - molecules: {mol_list(cls, mols, reps)}")
            md.append("")

    with open(os.path.join(output_dir, "V0.4.1_ACCURACY_BACKLOG.md"), "w") as f:
        f.write("\n".join(md) + "\n")

    return n_reports, n_success, failures, by_class, tiers


def append_trend(output_dir, n_reports, n_success, failures, by_class, tiers):
    ts = datetime.now().isoformat(timespec="seconds")
    pass_pct = 100 * n_success / max(n_reports, 1)
    tier_counts = {t: sum(len(m) for _, m in tiers.get(t, [])) for t in ("A", "B", "C", "D")}
    top = ",".join(
        f"{cls}:{n}" for cls, n in Counter({c: len(m) for c, m in by_class.items()}).most_common()
    )
    path = os.path.join(output_dir, "V0.4.1_TREND.tsv")
    new = not os.path.exists(path)
    with open(path, "a") as f:
        if new:
            f.write(
                "timestamp\tn_reports\tn_success\tpass_pct\tn_fail\t"
                "tierA\ttierB\ttierC\ttierD\tclasses\n"
            )
        f.write(
            f"{ts}\t{n_reports}\t{n_success}\t{pass_pct:.1f}\t{len(failures)}\t"
            f"{tier_counts['A']}\t{tier_counts['B']}\t{tier_counts['C']}\t"
            f"{tier_counts['D']}\t{top}\n"
        )


def delta_report(output_dir, rows, failures, by_class):
    """Print the change vs the last snapshot, then persist the new snapshot."""
    state_path = os.path.join(output_dir, ".v041_state.json")
    cur_fail = {r["molecule"]: r["class"] for r in failures}
    cur_all = {r["molecule"] for r in rows}
    cur_classes = set(by_class)

    prev = {}
    if os.path.exists(state_path):
        try:
            with open(state_path) as f:
                prev = json.load(f)
        except Exception:
            prev = {}
    prev_fail = prev.get("fail", {})
    prev_all = set(prev.get("all", []))
    prev_classes = set(prev.get("classes", []))

    new_fail = sorted(m for m in cur_fail if m not in prev_fail)
    resolved = sorted(m for m in prev_fail if m not in cur_fail)
    newly_covered = len(cur_all - prev_all)
    new_classes = sorted(cur_classes - prev_classes)

    print(
        f"coverage {len(cur_all)} (+{newly_covered})  "
        f"fails {len(cur_fail)}  classes {len(cur_classes)}"
    )
    if new_fail:
        print("  NEW fails: " + ", ".join(f"{m}[{cur_fail[m]}]" for m in new_fail))
    if resolved:
        print("  cleared:   " + ", ".join(resolved))
    if new_classes:
        print("  NEW class: " + ", ".join(new_classes))
    if not (new_fail or resolved or new_classes or newly_covered):
        print("  (no change)")

    with open(state_path, "w") as f:
        json.dump(
            {"all": sorted(cur_all), "fail": cur_fail, "classes": sorted(cur_classes)},
            f,
        )
    return bool(new_fail or resolved or new_classes)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--output-dir", required=True, help="Results dir with case_registry.json")
    args = ap.parse_args()
    output_dir = os.path.abspath(args.output_dir)

    with open(os.path.join(output_dir, "case_registry.json")) as f:
        rows = json.load(f)
    reps = load_reports(os.path.join(output_dir, "individual_reports"))

    n_reports, n_success, failures, by_class, tiers = build_backlog(rows, reps, output_dir)
    append_trend(output_dir, n_reports, n_success, failures, by_class, tiers)
    delta_report(output_dir, rows, failures, by_class)


if __name__ == "__main__":
    main()
