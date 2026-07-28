"""How big are the MEDZUR and GAVSED classes? (v0.4.13 Lane 2)

Both were handed to this release known on **n = 1 each**, from the v0.4.7 attachment lane:

* **GAVSED shape** -- acceptance rejected every conformer, so ``_select_by_geometry``'s
  geometry-ranked fallback returned one anyway, **and that fallback is not attachment-aware**.
  The returned structure has ligands off the metal. The attachment check guards *acceptance*,
  not *return*, so this is unreachable by tightening the check.
* **MEDZUR shape** -- attachment fully intact and independent re-perception still disagrees.
  Predicted by no model, found by that lane, and left unexplained. The attachment check is
  simply not the binding constraint here.

WHY THIS RUNS OFFLINE, AND WHY THAT IS NOT A SHORTCUT
=====================================================
v0.4.8 backfilled a ``coordination`` block into every individual report -- ``intact``,
``boundary_only``, and per-metal contact counts computed from the STORED generated geometry.
The classification this tool needs is therefore already a property of the frozen corpus, and
re-running the generator would not make it more true; it would make it a different corpus.
So this is a JSON join, it runs in seconds, and it is reproducible by anyone with the frozen
directory.

THE CONTROL ARM IS THE POINT, NOT A COURTESY
============================================
The standing rule this project pays for repeatedly:

    Before quoting an instrument, ask what a BROKEN version of it would print. If that is the
    same thing, you have not measured anything yet.

A broken version of this tool -- one where ``coordination.intact`` were mostly ``False`` for
reasons unrelated to failure, or mostly ``None`` and silently coerced -- would print a large,
plausible, entirely meaningless GAVSED class. Two things make the output falsifiable:

1. **The ``byte_exact`` control.** The same split is computed over molecules that round-trip
   perfectly. If ``DETACHED`` is common *there*, then detachment does not discriminate and the
   GAVSED number means nothing. The release must read the two side by side or not at all.
2. **``UNKNOWN`` is never folded into another class.** A missing or unparseable ``coordination``
   block is its own bucket with its own count. Every table prints its denominator.

Classification, in priority order (first match wins):

    NO_STRUCTURE  the run produced no structure to judge (status != success, or smiles_2 None).
                  Neither class -- there is no returned conformer whose attachment could be wrong.
    UNKNOWN       a structure exists but carries no usable `coordination` block.
    DETACHED      coordination.intact is False              -> the GAVSED shape
    BOUNDARY      intact, but every contact is within the cutoff tolerance (boundary_only)
    INTACT        intact and not boundary                   -> the MEDZUR shape

BOUNDARY is kept separate from INTACT deliberately. `coordination.boundary_only` means the
verdict rests on contacts sitting within 0.1 A of the cutoff, i.e. the attachment call itself is
the uncertain quantity. Merging it into either neighbour would let a tolerance choice decide the
headline.

Usage
-----
    V=$PWD/.venv/bin/python; export PYTHONPATH=$PWD/src
    $V tools/attach_class_audit.py \
        --results-dir tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest \
        --out docs/agentic-notes/v0.4.13/attach_class_audit.json

Reads ``bucket_report_honest.json`` from ``--results-dir`` for the AUTHORITATIVE bucket, never
an ad-hoc ``honest_class.endswith("->byte")`` -- that shortcut disagrees with the bucket report
by exactly the eight atom-count-gate molecules (v0.4.9's standing trap).
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

CLASSES = ("NO_STRUCTURE", "UNKNOWN", "DETACHED", "BOUNDARY", "INTACT")


def classify(report: dict) -> str:
    """Which attachment class does this molecule's RETURNED structure fall in?"""
    if report.get("status") != "success" or not report.get("smiles_2"):
        return "NO_STRUCTURE"
    coord = report.get("coordination")
    if not isinstance(coord, dict) or coord.get("intact") is None:
        return "UNKNOWN"
    if not coord.get("intact"):
        return "DETACHED"
    if coord.get("boundary_only"):
        return "BOUNDARY"
    return "INTACT"


def load_buckets(results_dir: Path) -> dict:
    """molecule -> (bucket, subclass) from the authoritative frozen report."""
    path = results_dir / "bucket_report_honest.json"
    if not path.exists():
        sys.exit(f"FATAL: no authoritative bucket report at {path}")
    rows = json.loads(path.read_text())
    if not isinstance(rows, list) or not rows:
        sys.exit(f"FATAL: {path} is not a non-empty list of per-molecule records")
    return {r["molecule"]: (r.get("bucket"), r.get("subclass")) for r in rows}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    results_dir = args.results_dir.resolve()
    buckets = load_buckets(results_dir)
    reports_dir = results_dir / "individual_reports"
    if not reports_dir.is_dir():
        sys.exit(f"FATAL: no individual_reports/ under {results_dir}")

    # bucket -> class -> [molecules]
    table: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    missing_report = []
    n_seen = 0

    for molecule, (bucket, _subclass) in sorted(buckets.items()):
        path = reports_dir / f"{molecule}.json"
        if not path.exists():
            missing_report.append(molecule)
            continue
        try:
            report = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            missing_report.append(molecule)
            continue
        n_seen += 1
        table[bucket][classify(report)].append(molecule)

    # ---- report -------------------------------------------------------------
    all_buckets = sorted(table)
    print(f"# attachment-class audit -- {results_dir.name}")
    print(f"\nmolecules in bucket report: {len(buckets)}")
    print(f"individual reports read:    {n_seen}")
    print(f"reports missing/unreadable: {len(missing_report)}")
    if missing_report:
        print(f"  (first 10: {missing_report[:10]})")

    header = f"\n| {'bucket':<18} | " + " | ".join(f"{c:>12}" for c in CLASSES) + " |   total |"
    print(header)
    print("|" + "-" * 20 + "|" + "|".join("-" * 14 for _ in CLASSES) + "|---------|")
    totals: Counter = Counter()
    for bucket in all_buckets:
        row = table[bucket]
        n = sum(len(row[c]) for c in CLASSES)
        cells = " | ".join(f"{len(row[c]):>12}" for c in CLASSES)
        print(f"| {bucket:<18} | {cells} | {n:>7} |")
        for c in CLASSES:
            totals[c] += len(row[c])
    cells = " | ".join(f"{totals[c]:>12}" for c in CLASSES)
    print(f"| {'ALL':<18} | {cells} | {sum(totals.values()):>7} |")

    # ---- the two named classes, with their control --------------------------
    failing = [b for b in all_buckets if b != "byte_exact"]
    gavsed = sum(len(table[b]["DETACHED"]) for b in failing)
    medzur = sum(len(table[b]["INTACT"]) for b in failing)
    boundary = sum(len(table[b]["BOUNDARY"]) for b in failing)
    no_struct = sum(len(table[b]["NO_STRUCTURE"]) for b in failing)
    unknown = sum(len(table[b]["UNKNOWN"]) for b in failing)
    n_failing = gavsed + medzur + boundary + no_struct + unknown

    ctrl = table.get("byte_exact", {})
    n_ctrl = sum(len(ctrl.get(c, [])) for c in CLASSES)
    ctrl_detached = len(ctrl.get("DETACHED", []))

    print(f"\n## The two classes, over the {n_failing} NON-byte_exact molecules\n")
    print(f"  GAVSED shape (DETACHED, returned anyway) : {gavsed:>5}")
    print(f"  MEDZUR shape (INTACT, indep disagrees)   : {medzur:>5}")
    print(f"  BOUNDARY (attachment call is uncertain)  : {boundary:>5}")
    print(f"  NO_STRUCTURE (nothing to judge)          : {no_struct:>5}")
    print(f"  UNKNOWN (no coordination block)          : {unknown:>5}")

    # ---- the cut that matters: key_equal is NOT a disagreement ---------------
    # `key_equal` molecules round-trip to the SAME comparison key and differ only in
    # presentation -- benign canonicalization, the encoder-side block that the donor fold and
    # v0.4.14 address. Counting them as "independent re-perception disagrees" inflates the
    # MEDZUR class with molecules whose geometry is fine and whose string is merely renumbered.
    # The MEDZUR shape as v0.4.7 defined it is a GENUINE round-trip failure with attachment
    # intact, so the honest denominator excludes key_equal.
    genuine = [b for b in failing if b != "key_equal"]
    g_detached = sum(len(table[b]["DETACHED"]) for b in genuine)
    g_intact = sum(len(table[b]["INTACT"]) for b in genuine)
    g_boundary = sum(len(table[b]["BOUNDARY"]) for b in genuine)
    g_nostruct = sum(len(table[b]["NO_STRUCTURE"]) for b in genuine)
    g_total = g_detached + g_intact + g_boundary + g_nostruct
    print(f"\n## Excluding key_equal (benign canonicalization) -- {g_total} genuine failures\n")
    print(f"  GAVSED shape (DETACHED)  : {g_detached:>5}")
    print(f"  MEDZUR shape (INTACT)    : {g_intact:>5}")
    print(f"  BOUNDARY                 : {g_boundary:>5}")
    print(f"  NO_STRUCTURE             : {g_nostruct:>5}")

    print("\n## CONTROL -- the same split over byte_exact (molecules that round-trip)\n")
    if n_ctrl:
        pct_ctrl = 100.0 * ctrl_detached / n_ctrl
        pct_fail = 100.0 * gavsed / n_failing if n_failing else 0.0
        print(f"  byte_exact DETACHED : {ctrl_detached:>5} / {n_ctrl} ({pct_ctrl:.2f}%)")
        print(f"  failing    DETACHED : {gavsed:>5} / {n_failing} ({pct_fail:.2f}%)")
        if pct_ctrl > 0:
            print(f"  enrichment          : {pct_fail / pct_ctrl:.1f}x")
        print(
            "\n  Read this BEFORE the class sizes above. If the two percentages are close, "
            "\n  detachment does not discriminate and the GAVSED number is not evidence."
        )
    else:
        print("  🔴 NO byte_exact CONTROL AVAILABLE -- the class sizes above are unfalsifiable.")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "results_dir": str(results_dir),
            "n_bucket_report": len(buckets),
            "n_reports_read": n_seen,
            "n_missing": len(missing_report),
            "missing": missing_report,
            "table": {b: {c: table[b][c] for c in CLASSES} for b in all_buckets},
            "summary": {
                "gavsed_detached": gavsed,
                "medzur_intact": medzur,
                "boundary": boundary,
                "no_structure": no_struct,
                "unknown": unknown,
                "n_failing": n_failing,
                "control_byte_exact_n": n_ctrl,
                "control_byte_exact_detached": ctrl_detached,
                "genuine_failures": {
                    "n": g_total,
                    "gavsed_detached": g_detached,
                    "medzur_intact": g_intact,
                    "boundary": g_boundary,
                    "no_structure": g_nostruct,
                },
            },
        }
        args.out.write_text(json.dumps(payload, indent=2))
        print(f"\nwrote {args.out}")

    print(f"\n#DONE {n_seen}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
