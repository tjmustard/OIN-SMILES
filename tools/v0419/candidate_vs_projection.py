"""v0.4.19 -- the candidate sweep against the projection that predicted it, row by row. Reads only.

The projection (docs/agentic-notes/v0.4.19/SERIALIZER_LANE.md section 6b) was assembled per row:

    unchanged-string rows  the release sweep's own outcome, corrected by the OFFLINE re-score
                           (tools/honest_rescore.py under the levers): serializer_ab_report_fix3.json
                           "offline" gains / losses / verified_gains / verified_losses
    changed-string rows    the live changed-set arm: "ab" self_gains / self_losses /
                           verified_gains / verified_losses  (109 rows)

so every molecule has a PROJECTED self-consistent and VERIFIED outcome. This script computes the
MEASURED outcome of the candidate sweep with sweep_two_numbers' own predicate (honest byte_exact;
parse-back clean under the levers; ruler ISO with no stereo element mirrored or changed), then
prints the two 2x2s and every disagreeing row with its elapsed_s in both runs -- the owner's
condition is "lands on its prediction", and the only rows allowed to differ are the ones the
300 s budget can flip on load.

    <main>/.venv/bin/python tools/v0419/candidate_vs_projection.py --sweep <results-v0.4.19-candidate-sweep>
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "v0417"))
from sweep_two_numbers import BAD_STEREO, ISO, _jsonl, _string_fault  # noqa: E402

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
RECORD = MAIN / "results-v0.4.18-release-sweep"
CENSUS = MAIN / "results-v0.4.18-release-census"
LANE = MAIN / "results-v0.4.19-serializer"
REPORT = LANE / "serializer_ab_report_fix3.json"
AB_ON = LANE / "ab_fix3" / "ab_on"
BUDGET_S = 280.0  # a row that spent this long in either run is on the 300 s boundary


def buckets(sweep):
    return {
        r["molecule"].removesuffix(".xyz"): r
        for r in json.loads((sweep / "bucket_report_honest.json").read_text())
    }


def elapsed(d, m):
    p = d / "individual_reports" / f"{m}.json"
    try:
        r = json.loads(p.read_text())
        v = (r.get("metrics") or {}).get("elapsed_s", r.get("elapsed_s"))  # NESTED, and a SUM
        return float(v) if v is not None else None
    except Exception:  # noqa: BLE001
        return None


def measured(sweep):
    B = buckets(sweep)
    G, P = _jsonl(sweep / "g_verdict.jsonl"), _jsonl(sweep / "parseback.jsonl")
    self_set, ver_set, why = set(), set(), {}
    for m, row in B.items():
        if row["bucket"] != "byte_exact":
            why[m] = f"bucket {row['bucket']}"
            continue
        self_set.add(m)
        sf, g = _string_fault(P.get(m)), G.get(m)
        if sf:
            why[m] = sf
        elif not g or g.get("graph") not in ISO:
            why[m] = "structure: graph differs"
        elif g.get("stereo") in BAD_STEREO:
            why[m] = f"structure: stereo {g['stereo']}"
        else:
            ver_set.add(m)
            why[m] = "NONE"
    return self_set, ver_set, why, B


def projected():
    rep = json.loads(REPORT.read_text())
    B = buckets(RECORD)
    self_set = {m for m, r in B.items() if r["bucket"] == "byte_exact"}
    ver_set = set()
    with open(CENSUS / "attribution_table.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["fault"] == "NONE":
                ver_set.add(row["molecule"])
    off, ab = rep["offline"], rep["ab"]
    self_set |= set(off["gains"]) | set(ab["self_gains"])
    self_set -= set(off["losses"]) | set(ab["self_losses"])
    ver_set |= {x[0] if isinstance(x, list) else x for x in off["verified_gains"]}
    ver_set |= set(ab["verified_gains"])
    ver_set -= {x[0] if isinstance(x, list) else x for x in off["verified_losses"]}
    ver_set -= set(ab["verified_losses"])
    changed = set(rep["changed"]["molecules"])
    assert (len(self_set), len(ver_set)) == (
        rep["projection"]["self"],
        rep["projection"]["verified"],
    ), (
        len(self_set),
        len(ver_set),
        rep["projection"],
    )
    return self_set, ver_set, changed


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--sweep", type=Path, required=True)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    p_self, p_ver, changed = projected()
    m_self, m_ver, why, B = measured(args.sweep)
    n = len(B)
    print(f"candidate sweep {args.sweep.name}: n={n}")
    print(
        f"  self-consistent  projected {len(p_self)}  measured {len(m_self)}  ({100 * len(m_self) / n:.2f}%)"
    )
    print(
        f"  VERIFIED         projected {len(p_ver)}  measured {len(m_ver)}  ({100 * len(m_ver) / n:.2f}%)"
    )

    def other_run(m):
        return AB_ON if m in changed else RECORD

    out = {
        "projected": {"self": len(p_self), "verified": len(p_ver)},
        "measured": {"self": len(m_self), "verified": len(m_ver)},
        "rows": [],
    }
    for name, P, M in (("self-consistent", p_self, m_self), ("VERIFIED", p_ver, m_ver)):
        only_p, only_m = sorted(P - M), sorted(M - P)
        print(
            f"\n== {name}: projected-not-measured {len(only_p)}   measured-not-projected {len(only_m)}"
        )
        tally = Counter()
        for tag, rows in (("LOST vs projection", only_p), ("GAINED vs projection", only_m)):
            for m in rows:
                e_c, e_o = elapsed(args.sweep, m), elapsed(other_run(m), m)
                boundary = any(
                    e is not None and e >= BUDGET_S for e in (e_c, e_o)
                ) or "timeout" in why.get(m, "")
                tally[(tag, "budget-boundary" if boundary else "NOT boundary")] += 1
                o_bucket = (buckets(other_run(m)).get(m) or {}).get("bucket")
                print(
                    f"   {tag:20s} {m:16s} {'changed-set' if m in changed else 'unchanged':11s} "
                    f"cand {why.get(m, '?'):32s} other bucket {o_bucket!s:14s} "
                    f"elapsed cand {e_c!s:>8} / other {e_o!s:>8}"
                    f"{'   BOUNDARY' if boundary else '   *** not a boundary row'}"
                )
                out["rows"].append(
                    {
                        "axis": name,
                        "tag": tag,
                        "molecule": m,
                        "changed_set": m in changed,
                        "candidate_why": why.get(m),
                        "elapsed_candidate": e_c,
                        "elapsed_other": e_o,
                        "boundary": boundary,
                    }
                )
        print("   ", dict(tally))
    if args.out:
        args.out.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
