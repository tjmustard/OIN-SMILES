"""GENERATOR REACH of a perception lever: rows whose input string did NOT move between two sweeps but
whose generated structure did. Reads only.

A perception lever changes more than ``smiles_1``: the generator re-encodes its candidate conformers
to accept one, so the same string can build a different structure (v0.4.19: OIN_CAP_IGNORES_METAL
alone rebuilt 44 of 4,891 unchanged-string rows). The offline re-score holds the structure fixed and
cannot see this -- this table is what it misses, measured.

For every molecule with the same ``smiles_1`` in both sweeps, the sha256 of
``structures/<m>_generated.xyz`` is compared; a structure built in ONE sweep only counts as a mover
too (kind ``built_one_side``: the same string reached a different outcome), as the ad hoc v0.4.19
table did. Each mover is listed with its bucket and its self-consistent / VERIFIED outcome in both
sweeps (``candidate_vs_projection.measured``, i.e. sweep_two_numbers' predicate).

    <main>/.venv/bin/python tools/v0419/generator_reach.py \\
        --record <results-v0.4.19-candidate-sweep> --candidate <results-v0.4.19-e2-candidate-sweep> \\
        --out <candidate>/unchanged_string_structure_movers.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from candidate_vs_projection import measured  # noqa: E402


def _sha(p: Path):
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--record", type=Path, required=True)
    ap.add_argument("--candidate", type=Path, required=True)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    r_self, r_ver, _, rB = measured(args.record)
    c_self, c_ver, _, cB = measured(args.candidate)
    same_string = [
        m
        for m in sorted(rB)
        if m in cB and rB[m].get("smiles_1") and rB[m]["smiles_1"] == cB[m].get("smiles_1")
    ]
    both_built = identical = neither = 0
    movers = []
    for m in same_string:
        a = _sha(args.record / "structures" / f"{m}_generated.xyz")
        b = _sha(args.candidate / "structures" / f"{m}_generated.xyz")
        if a is None and b is None:
            neither += 1
            continue
        if a is not None and b is not None:
            both_built += 1
            if a == b:
                identical += 1
                continue
        movers.append(
            {
                "molecule": m,
                "kind": "different_structure" if a and b else "built_one_side",
                "bucket_record": rB[m]["bucket"],
                "bucket_candidate": cB[m]["bucket"],
                "self": f"{'pass' if m in r_self else 'fail'}->{'pass' if m in c_self else 'fail'}",
                "verified": f"{'pass' if m in r_ver else 'fail'}->{'pass' if m in c_ver else 'fail'}",
            }
        )
    kinds = Counter(x["kind"] for x in movers)
    print(
        f"unchanged-string rows {len(same_string)}; built in neither {neither}; built in both"
        f" {both_built} (identical {identical}); MOVERS {len(movers)} {dict(kinds)}"
    )
    for axis in ("self", "verified"):
        t = Counter(x[axis] for x in movers)
        print(f"  {axis:9s} {dict(t)}   net {t['fail->pass'] - t['pass->fail']:+d}")
    for x in movers:
        print(
            f"   {x['molecule']:16s} {x['bucket_record']:>16s} -> {x['bucket_candidate']:<16s}"
            f" self {x['self']:10s} VERIFIED {x['verified']}"
        )
    if args.out:
        args.out.write_text(
            json.dumps(
                {
                    "record": args.record.name,
                    "candidate": args.candidate.name,
                    "unchanged_string": len(same_string),
                    "built_in_both": both_built,
                    "built_in_neither": neither,
                    "identical": identical,
                    "movers": movers,
                },
                indent=1,
            )
        )


if __name__ == "__main__":
    main()
