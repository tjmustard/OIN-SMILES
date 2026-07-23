"""SL4 hard-fail triage: classify the FF-only re-run of the v0.4.4 BASELINE hard_fail
cohort into disjoint populations on the SL0 key, and report how many reclaimed and why.

The v0.4.4 BASELINE froze 306 ``hard_fail`` molecules (status != success). SL4 demotes
RMSD to a diagnostic and adds a gated no-progress cutoff, then re-runs the cohort FF-only.
This tool reads that re-run's ``individual_reports/`` and answers the acceptance question:
how far did the cohort shrink, and *why* -- RMSD-demote reclaim vs passed-outright (budget /
determinism drift) vs still-genuinely-failed -- with the SL2 eta overlap kept separate so it
is not double-counted.

Usage:
    PYTHONPATH=src python tools/triage_hard_fails.py \
        --rerun-dir  <ff_only_rerun_dir_with_individual_reports> \
        --baseline-json <SL0 bucket_report.json> \
        [--output-dir <dir>]

Outputs ``triage_hard_fails.md`` / ``.json`` (default: the re-run dir).
"""

import argparse
import json
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

from oinsmiles.oin.compare import canonical_roundtrip_key  # noqa: E402

# Ordered error-head taxonomy for a still-failed re-run (first hit wins). Mirrors the
# harness error strings and tools/classify_failures.py, collapsed to SL4's disjoint
# reliability populations.
_STILL_FAILED_HEADS = [
    ("no_conformer", ("StructuralAssemblyError", "MetalloGen failed", "no conformer")),
    ("timeout", ("TimeoutException",)),
    ("atom_count", ("Atom count",)),
    ("string_mismatch", ("String mismatch",)),
    ("rmsd_mapping_failed", ("RMSD mapping failed",)),
    ("gen_exception", ("Generation/Verification failed",)),
]


def _still_failed_population(error: str) -> str:
    e = error or ""
    for name, needles in _STILL_FAILED_HEADS:
        if any(n in e for n in needles):
            return name
    return "other"


def _load_reports(rerun_dir):
    """Return {molecule: report} from <rerun_dir>/individual_reports/, skipping any
    report whose JSON is unreadable (a concurrent-write race on a duplicated input)."""
    ind = os.path.join(rerun_dir, "individual_reports")
    out = {}
    corrupt = []
    if not os.path.isdir(ind):
        return out, corrupt
    for fn in os.listdir(ind):
        if not fn.endswith(".json"):
            continue
        try:
            with open(os.path.join(ind, fn)) as f:
                r = json.load(f)
            out[r["molecule"]] = r
        except (json.JSONDecodeError, KeyError, OSError):
            corrupt.append(fn)
    return out, corrupt


def triage(rerun_dir, baseline_json):
    baseline = json.load(open(baseline_json))
    hard_fail = [r for r in baseline if r.get("bucket") == "hard_fail"]
    reports, corrupt = _load_reports(rerun_dir)

    rows = []
    for b in hard_fail:
        mol = b["molecule"]
        eta = bool(b.get("eta"))
        rr = reports.get(mol)
        if rr is None:
            rows.append(
                {"molecule": mol, "eta": eta, "outcome": "not_rerun", "population": "not_rerun"}
            )
            continue
        status = rr.get("status")
        m = rr.get("metrics") or {}
        if status == "success":
            # Reclaimed. Reason: did the RMSD demote carry it (string-exact + over-gate),
            # or did the fresh FF re-run pass the tightness gate outright?
            reason = "rmsd_demote" if m.get("rmsd_over_gate") else "passed_outright"
            # Sanity: a success must key-match (the harness gated on it).
            key_ok = None
            if rr.get("smiles_1") and rr.get("smiles_2"):
                try:
                    key_ok = canonical_roundtrip_key(rr["smiles_1"]) == canonical_roundtrip_key(
                        rr["smiles_2"]
                    )
                except Exception:
                    key_ok = None
            rows.append(
                {
                    "molecule": mol,
                    "eta": eta,
                    "outcome": "reclaimed",
                    "population": f"reclaimed:{reason}",
                    "rmsd": m.get("rmsd"),
                    "clash_count": m.get("clash_count"),
                    "key_ok": key_ok,
                }
            )
        else:
            pop = _still_failed_population(rr.get("error"))
            rows.append(
                {
                    "molecule": mol,
                    "eta": eta,
                    "outcome": "still_failed",
                    "population": pop,
                    "clash_count": m.get("clash_count"),
                }
            )
    return rows, corrupt, len(hard_fail)


def _counts(rows, key, pred=lambda r: True):
    out = {}
    for r in rows:
        if pred(r):
            out[r[key]] = out.get(r[key], 0) + 1
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


def _render_md(rows, corrupt, n_baseline):
    reclaimed = [r for r in rows if r["outcome"] == "reclaimed"]
    still = [r for r in rows if r["outcome"] == "still_failed"]
    not_rerun = [r for r in rows if r["outcome"] == "not_rerun"]
    key_bad = [r for r in reclaimed if r.get("key_ok") is False]

    def pct(n):
        return f"{100 * n / n_baseline:.1f}%" if n_baseline else "-"

    lines = []
    lines.append("# SL4 hard-fail triage — FF-only re-run of the BASELINE cohort\n")
    lines.append(f"Baseline `hard_fail` cohort: **{n_baseline}**  ")
    lines.append(
        f"Re-run reports found: **{len(rows) - len(not_rerun)}**  (not re-run: {len(not_rerun)})  "
    )
    if corrupt:
        lines.append(f"Unreadable reports skipped (write races): {len(corrupt)}  ")
    lines.append("")
    lines.append(
        f"- **Reclaimed (now key-exact success): {len(reclaimed)}  ({pct(len(reclaimed))} of cohort)**"
    )
    lines.append(f"- Still failed: {len(still)}")
    if not_rerun:
        lines.append(f"- Not re-run: {len(not_rerun)}")
    if key_bad:
        lines.append(f"- **WARNING: {len(key_bad)} reclaimed rows do NOT key-match (investigate)**")
    lines.append("")

    lines.append("## Reclaim reason\n")
    lines.append("| reason | count | eta | non-eta |")
    lines.append("|---|---:|---:|---:|")
    for reason in ("reclaimed:rmsd_demote", "reclaimed:passed_outright"):
        sub = [r for r in reclaimed if r["population"] == reason]
        e = sum(1 for r in sub if r["eta"])
        lines.append(f"| {reason.split(':')[1]} | {len(sub)} | {e} | {len(sub) - e} |")
    lines.append("")

    lines.append("## Still-failed populations (disjoint)\n")
    lines.append("| population | count | eta | non-eta |")
    lines.append("|---|---:|---:|---:|")
    for pop, cnt in _counts(still, "population").items():
        sub = [r for r in still if r["population"] == pop]
        e = sum(1 for r in sub if r["eta"])
        lines.append(f"| {pop} | {cnt} | {e} | {len(sub) - e} |")
    lines.append("")

    # clash-free fraction among reclaimed (the "no-clash" half, measured not gated)
    clashvals = [r.get("clash_count") for r in reclaimed if isinstance(r.get("clash_count"), int)]
    if clashvals:
        clashfree = sum(1 for c in clashvals if c == 0)
        lines.append("## Reclaimed structure quality (diagnostic)\n")
        lines.append(
            f"Clash-free (clash_count == 0): **{clashfree}/{len(clashvals)}** "
            f"({100 * clashfree / len(clashvals):.1f}%); "
            f"mean clash_count {sum(clashvals) / len(clashvals):.2f}."
        )
        lines.append("")

    lines.append("## Eta split (SL2 overlap)\n")
    r_eta = sum(1 for r in reclaimed if r["eta"])
    s_eta = sum(1 for r in still if r["eta"])
    lines.append(f"- Reclaimed: {r_eta} eta / {len(reclaimed) - r_eta} non-eta")
    lines.append(
        f"- Still failed: {s_eta} eta / {len(still) - s_eta} non-eta "
        f"(the eta residual is SL2's target; re-measure after SL2 lands)"
    )
    lines.append("")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--rerun-dir", required=True, help="Dir with individual_reports/ from the FF-only re-run."
    )
    ap.add_argument(
        "--baseline-json",
        required=True,
        help="SL0 bucket_report.json (source of the hard_fail set).",
    )
    ap.add_argument("--output-dir", default=None, help="Where to write (default: the re-run dir).")
    args = ap.parse_args()

    rows, corrupt, n_baseline = triage(args.rerun_dir, args.baseline_json)
    out_dir = args.output_dir or args.rerun_dir
    os.makedirs(out_dir, exist_ok=True)
    md = _render_md(rows, corrupt, n_baseline)
    with open(os.path.join(out_dir, "triage_hard_fails.md"), "w") as f:
        f.write(md)
    with open(os.path.join(out_dir, "triage_hard_fails.json"), "w") as f:
        json.dump({"n_baseline": n_baseline, "corrupt": corrupt, "rows": rows}, f, indent=1)
    print(md)
    print(f"\nWrote {out_dir}/triage_hard_fails.md / .json")


if __name__ == "__main__":
    main()
