"""Build the ARM 2 (round-trip) frozen golden manifest for `tools/gate_v047.sh`.

The v0.4.7 slow-100 cohort (see `docs/COHORT_v0.4.7.md`) was selected specifically
because ``smiles_1``/``smiles_2`` are ALREADY known and byte-exact, independently, in
two separate v0.4.5 sweep runs (`tools/select_slow_byte_exact.py`'s selection
predicate). So the golden values this tool freezes come from those EXISTING reports
-- no new generation is run here. Re-running the full round trip for all 62
molecules today would cost the exact wall-clock this cohort was chosen to be
expensive at, which is the quiet-phase sweep's job, not this lane's.

For each cohort molecule this tool:
  1. Loads its report from BOTH source dirs (``--source-a``, ``--source-b``).
  2. Asserts ``smiles_1``/``smiles_2`` agree byte-for-byte between the two sources
     -- a re-check, not a re-trust, of the selection step's own predicate.
  3. Emits ``name<TAB>sha256(smiles_1)<TAB>sha256(smiles_2)<TAB>len1<TAB>len2<TAB>eta``.

Sorted by name, followed by ``#DONE <n>``. Same TSV shape as ARM 1's manifest so both
arms are diffed the same way by ``gate_v047.sh``.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/gate_v047_build_arm2_golden.py \
        --names-file cohort_names.txt \
        --source-a <DS>/results-v0.4.5-sweep-partial-2697mols \
        --source-b <DS>/results-v0.4.5-rebaseline \
        --out tools/gate_v047_arm2_golden.tsv
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys

HAPTIC = re.compile(r"\{\d+[<>]\}")


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def load_report(results_dir: str, name: str) -> dict | None:
    path = os.path.join(results_dir, "individual_reports", f"{name}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--names-file", required=True, help="one basename (no .xyz) or name.xyz per line"
    )
    ap.add_argument("--source-a", required=True)
    ap.add_argument("--source-b", required=True)
    ap.add_argument("--out", default=None, help="also write the TSV to this path")
    args = ap.parse_args()

    names = []
    with open(args.names_file) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            names.append(line[: -len(".xyz")] if line.endswith(".xyz") else line)
    names = sorted(set(names))
    if not names:
        sys.exit(f"error: 0 names read from {args.names_file}")

    lines = []
    mismatches = []
    missing = []
    for name in names:
        ra = load_report(args.source_a, name)
        rb = load_report(args.source_b, name)
        if ra is None or rb is None:
            missing.append((name, ra is None, rb is None))
            continue
        s1a, s2a = ra.get("smiles_1"), ra.get("smiles_2")
        s1b, s2b = rb.get("smiles_1"), rb.get("smiles_2")
        if not (s1a and s2a and s1b and s2b):
            missing.append((name, "empty smiles"))
            continue
        if s1a.strip() != s1b.strip() or s2a.strip() != s2b.strip():
            mismatches.append(name)
            continue
        oin1, oin2 = s1a.strip(), s2a.strip()
        eta = "eta" if HAPTIC.search(oin1) else "-"
        lines.append(f"{name}\t{sha(oin1)}\t{sha(oin2)}\t{len(oin1)}\t{len(oin2)}\t{eta}")

    if missing:
        print(
            f"warning: {len(missing)} names missing/incomplete in one source: {missing}",
            file=sys.stderr,
        )
    if mismatches:
        sys.exit(
            f"error: {len(mismatches)} names disagree between {args.source_a} and "
            f"{args.source_b} (should be impossible -- these were selected FOR "
            f"byte-exact agreement): {mismatches}"
        )
    if len(lines) == 0:
        sys.exit("error: 0 golden rows produced -- refusing to write an empty golden manifest")

    manifest = "\n".join(lines)
    out_text = manifest + f"\n# MANIFEST_SHA256={sha(manifest)}\n#DONE {len(lines)}"
    print(out_text)
    if args.out:
        with open(args.out, "w") as f:
            f.write(out_text + "\n")
        print(f"# wrote {os.path.abspath(args.out)}", file=sys.stderr)


if __name__ == "__main__":
    main()
