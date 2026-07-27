#!/usr/bin/env python
"""Measure the harness's HAPTIC FALSE-POSITIVE rate on the default path.

The open question
----------------
`docs/agentic-notes/v0.4.7/ACCEPT_SCORED_v0.4.7.md` §4.7 established that the harness's success metric **cannot see a
detached haptic ligand**: it scores a round trip with `get_oin_string(gen_result.mol, coords)`,
which reuses the *generator's own bond graph*. If the generated geometry has let a Cp ring drift
off the metal, that bond graph still says it is coordinated, so the molecule scores a PASS.

That document names the consequence and explicitly leaves it open: *"Nobody currently knows the
false-positive rate with the lever OFF."* This measures it. It matters more than any lever,
because it sets whether the project's reported accuracy is the accuracy — and chasing 100% on an
instrument with unknown false positives is chasing a number that is not real.

Method (generator-free — reads stored structures, runs no 3D generation)
-----------------------------------------------------------------------
For each report that the harness scored a SUCCESS and for which a `*_generated.xyz` was stored:

1. Re-encode the **input** xyz with the CURRENT encoder -> `oin_in`.
2. Re-encode the **stored generated** xyz with the CURRENT encoder -> `oin_indep`. This is a full
   `XYZToSMILES().convert()`: bonds are perceived from coordinates alone, so it is blind to what
   the generator believed.

Both sides go through the same encoder, so a difference cannot be a version artifact — that is
why the input is re-encoded rather than read from the report's `smiles_1`.

3. Classify the difference. The distinction that matters is NOT "did the string change" (v0.4.5's
   `[[reencode-vs-harness-smiles2]]` already showed this path inflates `structural` ~19x through
   harmless presentation drift). It is **did the coordination survive**:

   `HAPTIC_LOST`   input had haptic slots, the independent encode has strictly fewer
                   -> a ring left the metal. The harness scored this a PASS. FALSE POSITIVE.
   `GEO_DEGRADED`  the `[El_GEO]` metal geometry tag changed -> the coordination polyhedron
                   itself differs. Also a false positive.
   `DENTICITY_LOST`     fewer distinct slots overall -> a donor detached.
   `KEY_DIFF_COORD_OK`  key differs but coordination is intact -> presentation drift, not a
                        false positive. Counted separately and NOT charged to the metric.
   `KEY_MATCH`          independent re-perception agrees. A genuine pass.

Usage:
    python tools/haptic_false_positive.py --results-dir DIR [--shard I:N] [--limit N]
                                          [--json out.json] [--haptic-only]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.oin.compare import canonical_roundtrip_key  # noqa: E402

#: A haptic slot HEAD -- `{n>}` / `{n<}`. The winding char marks the heading atom of an eta ring,
#: so one of these per haptic ligand. `oin/inline.py`'s detector uses the same shape.
_HAPTIC_RE = re.compile(r"\{(\d+)[<>]\}")
#: Any slot marker, haptic or not.
_SLOT_RE = re.compile(r"\{(\d+)[<>]?\}")
#: The metal's geometry tag, e.g. `[Fe_SAND]` -> `SAND`.
_GEO_RE = re.compile(r"\[[A-Z][a-z]?_([A-Z]{2,5})\]")


def haptic_slots(oin: str) -> set:
    return set(_HAPTIC_RE.findall(oin or ""))


def all_slots(oin: str) -> set:
    return set(_SLOT_RE.findall(oin or ""))


def geo_tag(oin: str):
    m = _GEO_RE.search(oin or "")
    return m.group(1) if m else None


def hapticity(oin: str) -> dict:
    """Slot number -> how many atoms carry it, i.e. the eta order of each haptic ligand.

    Needed because `haptic_slots` only sees the winding HEAD (`{n>}`), which survives a partial
    slip. An eta5 Cp that has slipped to eta2 keeps its head and its slot number, so full
    detachment and hapticity REDUCTION are different failures and only this distinguishes them.
    """
    counts: dict = {}
    for n in _SLOT_RE.findall(oin or ""):
        counts[n] = counts.get(n, 0) + 1
    return counts


def classify(oin_in: str, oin_indep: str) -> str:
    """Name the difference between the input encode and the independent re-encode."""
    if canonical_roundtrip_key(oin_in) == canonical_roundtrip_key(oin_indep):
        return "KEY_MATCH"
    g_in, g_out = geo_tag(oin_in), geo_tag(oin_indep)
    if g_in != g_out:
        return "GEO_DEGRADED"
    h_in, h_out = haptic_slots(oin_in), haptic_slots(oin_indep)
    if h_in and not h_out.issuperset(h_in):
        return "HAPTIC_LOST"
    s_in, s_out = all_slots(oin_in), all_slots(oin_indep)
    if len(s_out) < len(s_in):
        return "DENTICITY_LOST"
    # A haptic slot that survives but binds FEWER atoms is a ring slip (eta5 -> eta2), which the
    # slot-set checks above cannot see -- the head and the slot number both persist.
    cin, cout = hapticity(oin_in), hapticity(oin_indep)
    for n in h_in:
        if cout.get(n, 0) < cin.get(n, 0):
            return "HAPTICITY_REDUCED"
    return "KEY_DIFF_COORD_OK"


#: Classes that mean the harness credited a pass the geometry does not support.
FALSE_POSITIVE = {"HAPTIC_LOST", "GEO_DEGRADED", "DENTICITY_LOST", "HAPTICITY_REDUCED"}


def find_generated(results_dir: str, mol: str):
    for pat in (f"{mol}_generated.xyz", f"{mol}*_generated.xyz"):
        hits = sorted(glob.glob(os.path.join(results_dir, "structures", pat)))
        if hits:
            return hits[0]
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--shard", help="I:N")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--json")
    ap.add_argument(
        "--status",
        default="success",
        help="which reported status to audit. 'success' measures FALSE POSITIVES (passes whose "
        "coordination did not survive); 'failed' measures FALSE NEGATIVES (reported failures that "
        "independent re-perception says round-trip fine). Both halves are needed before any "
        "reported accuracy figure can be trusted.",
    )
    ap.add_argument(
        "--haptic-only",
        action="store_true",
        help="restrict to molecules whose input encode carries a haptic slot",
    )
    args = ap.parse_args()

    reports = sorted(glob.glob(os.path.join(args.results_dir, "individual_reports", "*.json")))
    if args.shard:
        i, n = (int(x) for x in args.shard.split(":"))
        reports = [r for k, r in enumerate(reports) if k % n == (i - 1) % n]
    if args.limit:
        reports = reports[: args.limit]

    conv = XYZToSMILES()
    rows, skipped = [], {"not_success": 0, "no_structure": 0, "no_input": 0, "encode_fail": 0}

    for path in reports:
        try:
            rep = json.load(open(path))
        except Exception:
            continue
        mol = rep.get("molecule") or os.path.basename(path)[:-5]
        if rep.get("status") != args.status:
            skipped["not_success"] += 1
            continue
        inp = rep.get("input_xyz")
        if not inp or not os.path.exists(inp):
            skipped["no_input"] += 1
            continue
        gen = find_generated(args.results_dir, mol)
        if not gen:
            skipped["no_structure"] += 1
            continue
        try:
            oin_in = conv.convert(inp)
            oin_indep = XYZToSMILES().convert(gen)
        except Exception as e:
            skipped["encode_fail"] += 1
            rows.append(
                {"molecule": mol, "case": "ENCODE_FAIL", "error": f"{type(e).__name__}: {e}"}
            )
            continue

        if args.haptic_only and not haptic_slots(oin_in):
            continue
        case = classify(oin_in, oin_indep)
        rows.append(
            {
                "molecule": mol,
                "case": case,
                "haptic_in": sorted(haptic_slots(oin_in)),
                "haptic_indep": sorted(haptic_slots(oin_indep)),
                "geo_in": geo_tag(oin_in),
                "geo_indep": geo_tag(oin_indep),
                "n_slots_in": len(all_slots(oin_in)),
                "n_slots_indep": len(all_slots(oin_indep)),
                "hapticity_in": hapticity(oin_in),
                "hapticity_indep": hapticity(oin_indep),
                "oin_in": oin_in,
                "oin_indep": oin_indep,
            }
        )
        print(f"  {mol:22s} {case}", flush=True)

    counts: dict = {}
    for r in rows:
        counts[r["case"]] = counts.get(r["case"], 0) + 1
    scored = [r for r in rows if r["case"] != "ENCODE_FAIL"]
    fp = [r for r in scored if r["case"] in FALSE_POSITIVE]
    haptic = [r for r in scored if r.get("haptic_in")]
    haptic_fp = [r for r in haptic if r["case"] in FALSE_POSITIVE]

    if args.status == "failed":
        # Inverted reading: the harness said FAIL. A KEY_MATCH here means independent
        # re-perception disagrees -- the molecule round-trips and was reported as a failure.
        fn = [r for r in scored if r["case"] == "KEY_MATCH"]
        print("\n================ FALSE-NEGATIVE RATE (reported failures) ================")
        print(f"  harness-scored FAILED molecules measured : {len(scored)}")
        for c in sorted(counts):
            print(f"    {c:20s} {counts[c]:5d}")
        if scored:
            print(
                f"\n  FALSE NEGATIVES (reported FAIL, independently round-trips): "
                f"{len(fn)}/{len(scored)} = {100 * len(fn) / len(scored):.1f}%"
            )
        print(f"  skipped: {skipped}")
        if args.json:
            with open(args.json, "w") as fh:
                json.dump({"rows": rows, "counts": counts, "skipped": skipped}, fh, indent=1)
            print(f"  wrote {args.json}")
        return 0

    print("\n================ HAPTIC FALSE-POSITIVE RATE (default path) ================")
    print(f"  harness-scored SUCCESS molecules measured : {len(scored)}")
    for c in sorted(counts):
        print(f"    {c:20s} {counts[c]:5d}")
    if scored:
        print(
            f"\n  FALSE POSITIVES (coordination not supported by the geometry): "
            f"{len(fp)}/{len(scored)} = {100 * len(fp) / len(scored):.1f}%"
        )
    if haptic:
        print(
            f"  ...restricted to HAPTIC inputs: {len(haptic_fp)}/{len(haptic)} = "
            f"{100 * len(haptic_fp) / len(haptic):.1f}%"
        )
    print(f"  skipped: {skipped}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"rows": rows, "counts": counts, "skipped": skipped}, fh, indent=1)
        print(f"  wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
