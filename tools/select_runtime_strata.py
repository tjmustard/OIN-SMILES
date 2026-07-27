"""Draw a STRATIFIED runtime benchmark cohort from ONE sweep results dir.

This is the v0.4.9 replacement for the selection half of ``tools/select_slow_byte_exact.py``.
That tool takes **top-N by elapsed_s** under a byte-exact predicate; this one takes a
**stratified sample across runtime bands x eta**. They are different selection algorithms,
so this is a second tool rather than a flag on the first -- the v0.4.7 slow-100 cohort must
stay re-derivable exactly as it was built.

WHY TOP-N IS THE WRONG SHAPE FOR GOAL B
=======================================
Goal B is ``max(elapsed_s) < 30 s`` -- a **p100 over the corpus**. A top-N cohort cut at
93.06 s (what v0.4.7 froze) cannot see the 30-93 s band, which is where most of the 994
over-30 s molecules live, and it has no fast control at all, so a change that speeds up the
tail by slowing everything else measures as a pure win. A p100 target needs a benchmark that
can detect a regression *anywhere*, not only where the current worst cases are.

AND WHY THE v0.4.7 PREDICATE IS THE WRONG PREDICATE
===================================================
The v0.4.7 cohort selected on ``smiles_1 == smiles_2`` -- the **scored** verdict, which
v0.4.8 showed over-states ``byte_exact`` by 10.34 points corpus-wide and by far more in the
slow tail. Cross-tabulated against the honest baseline, that cohort is
**28 byte->byte / 59 byte->FAIL / 4 byte->key / 8 no_structure / 1 key->FAIL**: 59% of it is
a molecule whose round trip only *looked* exact because the verdict read the generator's own
bond graph.

This tool therefore does NOT filter on any pass/fail predicate. A runtime benchmark must
contain the molecules that FAIL, because those are the slow ones -- excluding them is how
the 30-300 s band went missing in the first place. It records ``honest_class`` on every row
instead, so a green gate can never be misread as "the chemistry is right".

SAMPLING: SYSTEMATIC, ON THE ELAPSED-SORTED RANK
=================================================
Within each (band, eta) cell the members are sorted by ``elapsed_s`` and picks are spaced
evenly over the rank, ``round(i * (n-1) / (quota-1))``. Systematic rank sampling, deliberately,
because it is the only one of the three candidates that gets both properties this benchmark
needs:

  * **range coverage** -- picks span the cell's whole elapsed distribution, so a regression
    at either edge of a band is visible (a seeded random draw can miss an edge entirely);
  * **density preservation** -- rank is uniform, so more picks land where the distribution
    is dense (an even *quantile* spread would over-weight the sparse tail of each cell).

And it needs no RNG, so the draw is reproducible from the source dir alone -- no seed to
record, lose, or get wrong.

``metrics.elapsed_s`` IS NESTED, AND IT IS A SUM
=================================================
Two traps, both live:

1. ``elapsed_s`` lives inside ``report["metrics"]``. A top-level read does not raise -- it
   returns your ``.get`` default (typically 0) for every row, and a benchmark built on zeros
   looks perfectly reproducible. This tool indexes ``metrics`` explicitly and EXCLUDES a
   missing/malformed value rather than defaulting it.
2. ``elapsed_s`` is the **sum over up to three supervised attempts** (PASS 1 ``UFF_1``, then
   PASS 2 ``tier1`` and ``tier5``), each under its own SIGKILL -- see
   ``tools/test_dataset_roundtrip.py`` and
   ``docs/agentic-notes/v0.4.9/ELAPSED_S_IS_A_SUM_v0.4.9.md``. So a 759.9 s row is
   300 + 300 + 160, not one 2.5x overrun, and the bands below are bands of **per-molecule
   total** time, which is what Goal B is stated over. ``tier_passed`` is recorded per row so
   a consumer can tell a one-attempt row from a three-attempt one.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/select_runtime_strata.py \\
        --results-dir tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest \\
        --out spec/handoffs/v0.4.9/cohort_v049_strata.json \\
        --names-out spec/handoffs/v0.4.9/cohort_v049_names.txt
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys

HAPTIC = re.compile(r"\{\d+[<>]\}")

#: Runtime bands, on the per-molecule TOTAL ``metrics.elapsed_s`` (see the module docstring:
#: that figure is a sum across tiers). ``<30 s`` is the goal boundary itself; ``30-60 s`` is
#: the band the v0.4.7 cohort omitted entirely; ``>=300 s`` is where the budget is already
#: being violated and where 49.8% of the sweep's compute goes for a 1.0% pass rate.
BANDS: list[tuple[str, float, float]] = [
    ("fast", 0.0, 30.0),
    ("b30_60", 30.0, 60.0),
    ("b60_300", 60.0, 300.0),
    ("b300_plus", 300.0, float("inf")),
]

#: Per-(band, eta) quotas. Deliberately NOT proportional to corpus share: the ``fast`` cell
#: is a **control** -- it exists to catch "the tail got faster and everything else got
#: slower", which nothing in the v0.4.7 cohort could see -- while the three slow cells are
#: what Goal B is actually about and need enough members to resolve a shift. 90 fast
#: molecules detect a systematic slowdown comfortably; see the blind-spot section of
#: docs/agentic-notes/v0.4.9/RUNTIME_BENCHMARK_v0.4.9.md for what it does NOT detect.
QUOTAS: dict[tuple[str, bool], int] = {
    ("fast", True): 30,
    ("fast", False): 60,
    ("b30_60", True): 40,
    ("b30_60", False): 40,
    ("b60_300", True): 45,
    ("b60_300", False): 40,
    ("b300_plus", True): 35,
    ("b300_plus", False): 35,
}

#: Named molecules that must be in the cohort whatever the sampling picks, each with the
#: reason it is pinned. A benchmark that loses its own worst case to a rank stride is not a
#: benchmark. Missing pins are reported, never fatal -- a pin that is absent from the source
#: sweep is a fact about the source, and silently dropping it is what would be wrong.
PINNED: dict[str, str] = {
    "FOSNEI_comp_0": "worst observed elapsed_s in the 5k sweep (759.9 s)",
    "GOHWOQ_comp_0": "the 2.3x advisory-timeout overrun (60 s asked, 137.9 s spent)",
    "RAWJEG_comp_0": "boron cage that DOES assemble (~2.5 s) -- refutes the blanket fast-fail",
    "ULODUU_comp_0": "the second boron success (TET, 4 slots), budget-dependent at ~61.8 s",
    "XIQKOY_comp_0": "boron two-point proof: lever OFF fails in 0.87 s, ON runs past 340 s",
}


def band_for(elapsed_s: float) -> str:
    for name, lo, hi in BANDS:
        if lo <= elapsed_s < hi:
            return name
    raise AssertionError(f"unbanded elapsed_s={elapsed_s!r}")


def load_rows(results_dir: str) -> list[dict]:
    """Load every report that carries a usable ``metrics.elapsed_s``.

    No pass/fail predicate is applied -- see the module docstring. A report whose
    ``elapsed_s`` is missing or malformed is EXCLUDED and counted, never defaulted to 0.
    """
    indiv = os.path.join(results_dir, "individual_reports")
    if not os.path.isdir(indiv):
        sys.exit(f"error: {indiv} not found")

    rows: list[dict] = []
    n_files = 0
    n_no_elapsed = 0
    for fn in sorted(os.listdir(indiv)):
        if not fn.endswith(".json"):
            continue
        n_files += 1
        path = os.path.join(indiv, fn)
        try:
            with open(path) as f:
                r = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            print(f"warning: could not parse {path}: {e}", file=sys.stderr)
            continue

        metrics = r.get("metrics")
        elapsed = metrics.get("elapsed_s") if isinstance(metrics, dict) else None
        if not isinstance(elapsed, (int, float)) or isinstance(elapsed, bool):
            n_no_elapsed += 1
            continue

        s1 = r.get("smiles_1") or ""
        rows.append(
            {
                "name": fn[: -len(".json")],
                "elapsed_s": float(elapsed),
                "eta": bool(HAPTIC.search(s1)),
                "band": band_for(float(elapsed)),
                "status": r.get("status"),
                "tier_passed": r.get("tier_passed"),
                "honest_class": r.get("honest_class"),
            }
        )

    if n_files == 0:
        sys.exit(
            f"error: 0 reports found under {indiv} -- refusing to draw a cohort from an "
            "empty corpus (a serene '0/0' that exits 0 has destroyed a full A/B here before)"
        )
    print(
        f"# loaded {len(rows)}/{n_files} reports ({n_no_elapsed} excluded for a "
        f"missing/malformed metrics.elapsed_s)",
        file=sys.stderr,
    )
    return rows


def systematic(cell: list[dict], quota: int) -> list[dict]:
    """Take exactly ``quota`` members of ``cell`` by even spacing on the elapsed-sorted rank.

    ``cell`` is sorted ascending by ``elapsed_s`` with ``name`` as the tie-break, so the draw
    is a pure function of the source dir -- no seed to record, lose, or get wrong. When the
    cell is at or under quota the whole cell is taken.

    Spacing is ``round(i * (n-1) / (quota-1))`` rather than an integer stride ``cell[::k]``.
    An integer stride under-delivers whenever ``n`` is not a multiple of the quota AND drops
    the cell's upper edge: 151 members at quota 40 gives ``k=4``, 38 picks, and a last pick at
    rank 148 -- losing exactly the three slowest molecules a p100 benchmark most wants. The
    fractional form hits both endpoints and the quota exactly.
    """
    cell = sorted(cell, key=lambda r: (r["elapsed_s"], r["name"]))
    if len(cell) <= quota:
        return list(cell)
    if quota <= 1:
        return cell[-quota:] if quota else []
    n = len(cell)
    idx = sorted({round(i * (n - 1) / (quota - 1)) for i in range(quota)})
    return [cell[j] for j in idx]


def draw(rows: list[dict]) -> tuple[list[dict], dict]:
    by_name = {r["name"]: r for r in rows}
    selected: dict[str, dict] = {}

    # Pins first, so a pinned molecule occupies its cell's quota rather than being added on
    # top of a full cell -- otherwise the strata counts in the manifest would be a lie.
    #
    # A pin ABSENT from the source sweep is still taken. The cohort is a benchmark, not a
    # subsample: GOHWOQ, ULODUU and XIQKOY are named evidence for this release and three of
    # them are simply not in the frozen seed-42 5k draw. Dropping them would quietly remove
    # the release's own fixtures from the release's own benchmark. They carry
    # ``elapsed_s: None`` and ``band: "unbaselined"`` so no consumer can mistake them for
    # rows the source dir has a runtime baseline for -- they have none, and the golden
    # builder must be told to expect that (``--allow-unbaselined``).
    pins_found, pins_unbaselined = [], []
    for name, reason in PINNED.items():
        r = by_name.get(name)
        if r is None:
            selected[name] = {
                "name": name,
                "elapsed_s": None,
                "eta": None,
                "band": "unbaselined",
                "status": None,
                "tier_passed": None,
                "honest_class": None,
                "pinned": reason,
            }
            pins_unbaselined.append(name)
            continue
        selected[name] = dict(r, pinned=reason)
        pins_found.append(name)

    cell_stats = {}
    for band_name, _lo, _hi in BANDS:
        for eta in (True, False):
            quota = QUOTAS[(band_name, eta)]
            cell = [r for r in rows if r["band"] == band_name and r["eta"] is eta]
            already = [r for r in cell if r["name"] in selected]
            remaining_quota = max(0, quota - len(already))
            candidates = [r for r in cell if r["name"] not in selected]
            picked = systematic(candidates, remaining_quota)
            for r in picked:
                selected[r["name"]] = dict(r, pinned=None)
            cell_stats[f"{band_name}/{'eta' if eta else 'non-eta'}"] = {
                "population": len(cell),
                "quota": quota,
                "pinned_already_in_cell": len(already),
                "selected": len(already) + len(picked),
            }

    out = sorted(selected.values(), key=lambda r: r["name"])
    meta = {
        "cells": cell_stats,
        "pins_found": sorted(pins_found),
        "pins_unbaselined": sorted(pins_unbaselined),
    }
    return out, meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--results-dir",
        required=True,
        help="THE single sweep results dir to draw from. Never intersect two -- that "
        "mistake destroyed the v0.4.7 cohort once already (see COHORT_v0.4.7.md).",
    )
    ap.add_argument("--out", default=None, help="write the full JSON manifest here")
    ap.add_argument("--names-out", default=None, help="write one <name>.xyz per line here")
    args = ap.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    rows = load_rows(results_dir)
    selected, meta = draw(rows)

    if not selected:
        sys.exit("error: 0 molecules selected -- refusing to write an empty cohort")

    payload = {
        "source_results_dir": results_dir,
        "selection": "stratified systematic on elapsed-sorted rank; no pass/fail predicate",
        "bands": [
            {"name": n, "lo_s": lo, "hi_s": None if hi == float("inf") else hi}
            for n, lo, hi in BANDS
        ],
        "quotas": {f"{b}/{'eta' if e else 'non-eta'}": q for (b, e), q in QUOTAS.items()},
        "corpus_n": len(rows),
        "n_selected": len(selected),
        "cells": meta["cells"],
        "pins_found": meta["pins_found"],
        "pins_unbaselined": meta["pins_unbaselined"],
        "molecules": selected,
    }

    text = json.dumps(payload, indent=2)
    if args.out:
        with open(os.path.abspath(args.out), "w") as f:
            f.write(text + "\n")
        print(f"# wrote {os.path.abspath(args.out)}", file=sys.stderr)
    else:
        print(text)

    if args.names_out:
        with open(os.path.abspath(args.names_out), "w") as f:
            for r in selected:
                f.write(f"{r['name']}.xyz\n")
        print(f"# wrote {os.path.abspath(args.names_out)}", file=sys.stderr)

    for label, st in meta["cells"].items():
        print(
            f"#   {label:22s} population={st['population']:5d} quota={st['quota']:3d} "
            f"selected={st['selected']:3d}",
            file=sys.stderr,
        )
    if meta["pins_unbaselined"]:
        print(
            f"# NOTE {len(meta['pins_unbaselined'])} pin(s) are not in the source sweep and "
            f"carry NO runtime baseline: {meta['pins_unbaselined']}",
            file=sys.stderr,
        )
    print(f"#DONE {len(selected)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
