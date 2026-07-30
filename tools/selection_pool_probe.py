#!/usr/bin/env python3
"""v0.4.15: is the remaining gap a SELECTION problem or a CONSTRUCTION problem?

WHY THIS EXISTS
===============
v0.4.15 built two selection-side levers and both recovered exactly ZERO over 1107
molecule-pairs, with ``GENERATED output moved = 0`` everywhere:

    OIN_ATTACH_RETURN        prefer an ATTACHED conformer on the return path
    OIN_ACCEPT_STRING_EXACT  accept only a conformer whose OIN STRING matches

🔴 A LEVER THAT WAS NEVER REACHED PRINTS EXACTLY THAT TOO. So a null here is worthless until the
levers are shown to FIRE. This probe supplies that, quantitatively: it counts, per molecule, how
many candidate conformers each lever examined and rejected. A non-zero rejection count with a
zero recovery count is the finding -- *the pool contains candidates and none of them is better* --
and it points at construction, not selection.

WHAT IT READS
=============
Telemetry counters (``OIN_TELEMETRY=1``), which are recorded at the exact decision sites:

  adapter.attach_return_winding_skip          a winding-matching conformer was DETACHED
  adapter.attach_return_winding_all_detached  ...and every one of them was, so the first came back
  adapter.attach_return_fallback_promoted     a LATER conformer was attached and got promoted (a WIN)
  adapter.attach_return_none_attached         no conformer in the pool held its sites
  adapter.attach_return_all_scored_detached   the best geometry-classified conformer was detached
  pool.accept_incumbent_recorded              a key-matching conformer was not string-exact
  pool.accept_incumbent_returned              ...and none ever was, so the incumbent came back
  adapter.string_exact_early_exit_incumbent   same, on the selection re-scan

``*_promoted`` is the only counter that means the lever CHANGED the answer. If it is 0 while the
skip counters are large, selection had candidates and rejected all of them.

Usage:
    python tools/selection_pool_probe.py --molecules-file <list> --lever <NAME> \
        --cohort-dir <dir> [--limit 25] --out-json out.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

WIN_SITES = (
    "adapter.attach_return_fallback_promoted",  # Lane 1 changed the returned conformer
)
FIRE_SITES = (
    "adapter.attach_return_winding_skip",
    "adapter.attach_return_winding_all_detached",
    "adapter.attach_return_none_attached",
    "adapter.attach_return_all_scored_detached",
    "pool.accept_incumbent_recorded",
    "pool.accept_incumbent_returned",
    "adapter.string_exact_incumbent",
    "adapter.string_exact_early_exit_incumbent",
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--molecules-file", required=True)
    ap.add_argument("--cohort-dir", required=True)
    ap.add_argument("--lever", required=True)
    ap.add_argument("--limit", type=int, default=25)
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    os.environ["OIN_TELEMETRY"] = "1"
    os.environ[args.lever] = "1"

    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    from oinsmiles import XYZToSMILES
    from oinsmiles.generation import _telemetry
    from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen

    names = [
        n.strip() for n in open(args.molecules_file) if n.strip() and not n.lstrip().startswith("#")
    ][: args.limit]

    rows = []
    for i, mol in enumerate(names, 1):
        p = os.path.join(args.cohort_dir, f"{mol}.xyz")
        if not os.path.exists(p):
            rows.append({"molecule": mol, "error": "input_missing"})
            continue
        try:
            oin_in = XYZToSMILES().convert(p)
        except Exception as exc:  # noqa: BLE001 -- unencodable input is data
            rows.append({"molecule": mol, "error": f"encode_in:{type(exc).__name__}"})
            continue
        with _telemetry.collecting():
            gen = OIN3DGeneratorMetallogen(
                optimizer=None, ensemble_size=1, timeout=args.timeout, ff_params=None
            )
            err = None
            try:
                gen.generate(oin_in)
            except Exception as exc:  # noqa: BLE001 -- a failed generation is data
                err = f"generate:{type(exc).__name__}"
            counts = _telemetry.counts()
        fired = {k: v for k, v in counts.items() if k in FIRE_SITES}
        won = {k: v for k, v in counts.items() if k in WIN_SITES}
        rows.append(
            {
                "molecule": mol,
                "error": err,
                "fired": fired,
                "won": won,
                "n_rejected": sum(fired.values()),
            }
        )
        print(
            f"  [{i}/{len(names)}] {mol} rejected={sum(fired.values())} won={sum(won.values())}",
            flush=True,
        )

    scored = [r for r in rows if "fired" in r]
    n_fired = sum(1 for r in scored if r["n_rejected"] > 0)
    n_won = sum(1 for r in scored if sum(r["won"].values()) > 0)
    total_rejected = sum(r["n_rejected"] for r in scored)
    out = {
        "lever": args.lever,
        "molecules_file": os.path.basename(args.molecules_file),
        "n_requested": len(names),
        "n_scored": len(scored),
        "n_molecules_where_lever_FIRED": n_fired,
        "n_molecules_where_lever_CHANGED_the_answer": n_won,
        "total_candidate_conformers_rejected": total_rejected,
        "rows": rows,
    }
    with open(args.out_json, "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"\n=== {args.lever} over {os.path.basename(args.molecules_file)} ===")
    print(f"  scored                              {len(scored)}")
    print(f"  lever FIRED on                      {n_fired}")
    print(f"  lever CHANGED the returned answer   {n_won}")
    print(f"  candidate conformers REJECTED       {total_rejected}")
    if n_fired and not n_won:
        print(
            "\n  ⇒ selection had candidates and rejected ALL of them. The pool does not contain a\n"
            "    better conformer, so this is a CONSTRUCTION finding, not a selection bug."
        )
    if not n_fired:
        print("\n  🔴 the lever never fired -- this is a WIRING failure, not a finding.")
    print(f"\nwrote {args.out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
