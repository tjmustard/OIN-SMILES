"""The two numbers a sweep must be read with: SELF-CONSISTENT and VERIFIED. Reads only.

    self-consistent   ``byte_exact`` under the honest score:  E(x) == E(G(E(x)))
    VERIFIED          the census's fault ``NONE``: a pass whose STRING describes the input
                      (parse-back clean) AND whose GENERATED STRUCTURE the neutral ruler finds
                      graph-isomorphic to the input with no stereo element mirrored or changed

The census measured 77.16% / 69.24% on the v0.4.14 sweep of record. A canonicality lever can raise
the first number by teaching the second encode to agree with the first, whatever was built; only
the second number says the structure is right. Every sweep from v0.4.17 on reports both.

THE CONTROL IS NOT OPTIONAL. ``--control`` runs this file's own predicate over the sweep of record
with the census's own ruler verdicts and REQUIRES 3,858 and 3,462. This script re-implements
``tools/census/attribution_table.py``'s rule order in six lines; if those six lines are wrong, the
control is what says so -- before a new sweep's number is believed, not after.

A sweep directory needs (``tools/v0417/post_sweep.sh`` makes all three):

    bucket_report_honest.json    roundtrip_bucket_report.py --score honest
    g_verdict.jsonl              tools/census/g_vs_input.py --sweep <dir> --out <dir>
    parseback.jsonl              tools/census/string_sufficiency.py parseback --sweep <dir> --out <dir>

``parseback.jsonl`` is re-derived from the NEW ``smiles_1`` rather than carried from the census: a
slot relabeling should not move a parse-back verdict, and "should not" is measured here.

    <main>/.venv/bin/python tools/v0417/sweep_two_numbers.py --control
    <main>/.venv/bin/python tools/v0417/sweep_two_numbers.py --sweep <results-v0.4.17-sweep>
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
RECORD, CENSUS = MAIN / "results-v0.4.14-sweep", MAIN / "results-census"
ISO = ("ISO", "ISO_MARGINAL", "ISO_CLASH")
BAD_STEREO = ("MIRROR", "MIRROR_PARTIAL", "DIFFERENT")
N, RECORD_SELF, RECORD_VERIFIED = 5000, 3858, 3462


def _jsonl(path):
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not rows or not rows[-1].startswith("#DONE"):
        sys.exit(f"ABORT: {path} has no #DONE trailer -- refusing a partial file")
    return {r["molecule"]: r for r in map(json.loads, rows[:-1])}


def _string_fault(pb):
    """``attribution_table.py`` rule 2, verbatim in effect: does the string describe the input?"""
    if pb is None:
        return None
    if int(pb.get("n_unslotted") or 0):
        return "DATA_MULTI" if pb.get("parser_graph") in ISO else "P_DETACHED"
    adapter = pb.get("adapter_graph")
    if adapter is not None and adapter not in ISO:
        return "E1_GRAPH"
    return "E1_HCOUNT" if pb.get("adapter_h_decoration") == "DIFF" else None


def two_numbers(sweep: Path, verdicts: Path, parseback: Path):
    B = {
        r["molecule"].removesuffix(".xyz"): r
        for r in json.loads((sweep / "bucket_report_honest.json").read_text())
    }
    G, P = _jsonl(verdicts), _jsonl(parseback)
    why = Counter()
    n_self = n_ver = 0
    for m, row in B.items():
        if row["bucket"] != "byte_exact":
            continue
        n_self += 1
        sf, g = _string_fault(P.get(m)), G.get(m)
        if sf:
            why[sf] += 1
        elif not g or g.get("graph") not in ISO:
            why["structure: graph differs (or no ruler verdict)"] += 1
        elif g.get("stereo") in BAD_STEREO:
            why[f"structure: stereo {g['stereo']}"] += 1
        else:
            n_ver += 1
    return len(B), n_self, n_ver, why, Counter(r["bucket"] for r in B.values())


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--sweep", type=Path)
    ap.add_argument("--control", action="store_true")
    args = ap.parse_args()

    n, s, v, _why, _b = two_numbers(RECORD, CENSUS / "g_verdict.jsonl", CENSUS / "parseback.jsonl")
    ok = (n, s, v) == (N, RECORD_SELF, RECORD_VERIFIED)
    print(
        f"CONTROL on the sweep of record: n={n} self={s} verified={v}   "
        f"expected {N} / {RECORD_SELF} / {RECORD_VERIFIED}   {'OK' if ok else 'MISMATCH'}"
    )
    if not ok:
        sys.exit("ABORT: this predicate does not reproduce the census. Nothing below is printed.")
    if args.control or not args.sweep:
        return

    n, s, v, why, buckets = two_numbers(
        args.sweep, args.sweep / "g_verdict.jsonl", args.sweep / "parseback.jsonl"
    )
    print(f"\n{args.sweep.name}   DENOMINATOR n={n}  (must be {N})")
    print("   buckets:", dict(buckets.most_common()))
    print(
        f"   self-consistent  {s:5d} / {n} = {100 * s / n:.2f}%   (record {100 * RECORD_SELF / N:.2f}%)"
    )
    print(
        f"   VERIFIED         {v:5d} / {n} = {100 * v / n:.2f}%   (record {100 * RECORD_VERIFIED / N:.2f}%)"
    )
    print(f"   passes the metric cannot verify: {s - v}")
    for k, c in why.most_common():
        print(f"      {c:5d}  {k}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
