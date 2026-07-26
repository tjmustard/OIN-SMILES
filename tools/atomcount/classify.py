"""Phase-1 classifier for the v0.4.5 `atom_count` hard-fail class.

Pure file arithmetic: compares the input XYZ against the stored generated XYZ by
element histogram, and against the recorded OIN string.  No generation, no RDKit
embedding -- runs in seconds over all 74 molecules.
"""

from __future__ import annotations

import argparse
import collections
import json
import re
from pathlib import Path

RESULTS = Path(
    "/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042"
)
WORKLIST = Path(
    "/home/tjmustard/Documents/GitHub/OIN-SMILES/spec/handoffs/v0.4.5/hard_fail_worklists.json"
)


def read_xyz(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    n = int(lines[0].split()[0])
    out = []
    for ln in lines[2 : 2 + n]:
        parts = ln.split()
        if parts:
            out.append(parts[0])
    return out


def hist(elems: list[str]) -> collections.Counter:
    return collections.Counter(elems)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    work = json.loads(WORKLIST.read_text())["atom_count"]
    rows = []
    for entry in work:
        mol = entry["molecule"]
        rep = json.loads((RESULTS / "individual_reports" / f"{mol}.json").read_text())
        inp = Path(rep["input_xyz"])
        gen = RESULTS / "structures" / f"{mol}_generated.xyz"
        m = re.search(r"Input (\d+) != Gen (\d+)", entry["error"])
        row = {
            "molecule": mol,
            "reported_input": int(m.group(1)),
            "reported_gen": int(m.group(2)),
            "stage": (re.search(r"mismatch at (\S+?)\.", entry["error"]) or [None, "?"])[1],
            "elapsed_s": entry["elapsed_s"],
            "oin_1": rep.get("smiles_1"),
            "oin_2": rep.get("smiles_2"),
            "oin_stable": rep.get("smiles_1") == rep.get("smiles_2"),
            "input_xyz": str(inp),
            "gen_exists": gen.exists(),
        }
        if inp.exists():
            ih = hist(read_xyz(inp))
            row["input_n"] = sum(ih.values())
            row["input_hist"] = dict(ih)
        if gen.exists():
            gh = hist(read_xyz(gen))
            row["gen_n"] = sum(gh.values())
            row["gen_hist"] = dict(gh)
        if inp.exists() and gen.exists():
            diff = {}
            for el in set(ih) | set(gh):
                d = gh.get(el, 0) - ih.get(el, 0)
                if d:
                    diff[el] = d
            row["diff"] = diff
            row["h_only"] = set(diff) <= {"H"}
            row["heavy_delta"] = sum(v for k, v in diff.items() if k != "H")
        rows.append(row)

    # ---- summary -------------------------------------------------------
    print(f"n = {len(rows)}")
    have = [r for r in rows if "diff" in r]
    print(f"stored generated structure present: {len(have)}/{len(rows)}")

    print("\n-- stored-gen vs reported-gen consistency --")
    same = sum(1 for r in have if r["gen_n"] == r["reported_gen"])
    print(f"  stored gen atom count == reported gen: {same}/{len(have)}")
    samein = sum(1 for r in have if r["input_n"] == r["reported_input"])
    print(f"  input xyz atom count == reported input: {samein}/{len(have)}")

    print("\n-- differing elements (multiset of element symbols that differ) --")
    c = collections.Counter(tuple(sorted(r["diff"])) for r in have)
    for k, v in c.most_common():
        print(f"  {v:3d}  {', '.join(k) if k else '(none)'}")

    print("\n-- H-only vs heavy-atom involvement --")
    honly = [r for r in have if r["h_only"]]
    heavy = [r for r in have if not r["h_only"]]
    print(f"  H-only differences : {len(honly)}")
    print(f"  involves heavy atom: {len(heavy)}")

    print("\n-- per-element net delta summed over the class --")
    tot = collections.Counter()
    for r in have:
        tot.update(r["diff"])
    for el, d in sorted(tot.items(), key=lambda kv: -abs(kv[1])):
        print(f"  {el:3s} {d:+5d}")

    print("\n-- OIN re-encode stability (smiles_1 == smiles_2) --")
    print(f"  stable: {sum(1 for r in rows if r['oin_stable'])}/{len(rows)}")

    print("\n-- detail --")
    for r in sorted(have, key=lambda r: (not r["h_only"], r["gen_n"] - r["input_n"])):
        tag = "H-only" if r["h_only"] else "HEAVY "
        d = ", ".join(f"{k}{v:+d}" for k, v in sorted(r["diff"].items()))
        print(
            f"  {r['molecule']:22s} in={r['input_n']:4d} gen={r['gen_n']:4d} "
            f"d={r['gen_n'] - r['input_n']:+4d} {tag} [{d}]"
            f"{'  oin-stable' if r['oin_stable'] else ''}"
        )

    if args.out:
        args.out.write_text(json.dumps(rows, indent=1, default=str))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
