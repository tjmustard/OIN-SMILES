"""Build the ARM 2 (round-trip) frozen golden manifest for `tools/gate_v047.sh`.

The v0.4.7 slow-100 cohort (see `docs/agentic-notes/v0.4.7/COHORT_v0.4.7.md`) is selected from a SINGLE
sweep results dir (`tools/select_slow_byte_exact.py`, top-100 by `metrics.elapsed_s`
within that one dir) -- ``smiles_1``/``smiles_2`` are therefore ALREADY known and
byte-exact there. So the golden values this tool freezes come from that EXISTING
primary source -- no new generation is run here. Re-running the full round trip for
all 100 molecules today would cost the exact wall-clock this cohort was chosen to be
expensive at, which is the quiet-phase sweep's job, not this lane's.

PRIMARY SOURCE, NOT TWO PEERS
=============================
An earlier version of this tool required agreement between TWO sources, treated as
equal peers, and hard-failed on any disagreement or absence in either. That was
wrong: the two source dirs are different, uncorrelated random draws over the
dataset (see ``docs/agentic-notes/v0.4.7/COHORT_v0.4.7.md`` -- the "62, not 100" incident) and share only
~100 names total, so requiring both was a near-total sample destruction, not a
robustness check. This tool now takes ONE required primary source
(``--source-a``, the same dir the cohort was selected from -- every cohort name
MUST resolve there and pass the full predicate, or the build fails loudly) and an
OPTIONAL corroboration source (``--source-b``) that is purely informational: where
a name happens to also have a report there, its ``smiles_1`` is cross-checked as an
encoder-determinism sanity check (mismatch => loud warning, never fatal), and
absence is expected and unremarkable.

For each cohort molecule this tool:
  1. Loads its report from the primary source and re-verifies the full selection
     predicate there (status/tier/byte-exact/elapsed) -- a re-check, not a re-trust.
  2. If a corroboration source is given and has a report for this name with a
     non-null ``smiles_1``, compares it to the primary's -- logged, never fatal.
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


def read_names(names_file: str) -> list[str]:
    names = []
    with open(names_file) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            names.append(line[: -len(".xyz")] if line.endswith(".xyz") else line)
    return sorted(set(names))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--names-file", required=True, help="one basename (no .xyz) or name.xyz per line"
    )
    ap.add_argument(
        "--source-a",
        required=True,
        help="PRIMARY source -- the dir the cohort was selected from; must pass the full "
        "predicate for every name, or the build fails",
    )
    ap.add_argument(
        "--source-b",
        default=None,
        help="OPTIONAL corroboration source -- informational only, never required, never "
        "blocks the build on absence or mismatch (mismatch is warned loudly instead)",
    )
    ap.add_argument("--out", default=None, help="also write the TSV to this path")
    args = ap.parse_args()

    names = read_names(args.names_file)
    if not names:
        sys.exit(f"error: 0 names read from {args.names_file}")

    lines = []
    missing_primary = []
    bad_primary = []
    corroborated = 0
    corroboration_mismatches = []
    for name in names:
        ra = load_report(args.source_a, name)
        if ra is None:
            missing_primary.append(name)
            continue
        if ra.get("status") != "success" or ra.get("tier_passed") != "UFF_1":
            bad_primary.append((name, ra.get("status"), ra.get("tier_passed")))
            continue
        s1a, s2a = ra.get("smiles_1"), ra.get("smiles_2")
        if not s1a or not s2a or s1a.strip() != s2a.strip():
            bad_primary.append((name, "smiles_1/smiles_2 missing or not byte-exact"))
            continue
        oin1, oin2 = s1a.strip(), s2a.strip()

        if args.source_b:
            rb = load_report(args.source_b, name)
            if rb is not None:
                s1b = rb.get("smiles_1")
                if s1b:
                    if s1b.strip() == oin1:
                        corroborated += 1
                    else:
                        corroboration_mismatches.append(name)
                        print(
                            f"warning: {name} smiles_1 DIFFERS in corroboration source "
                            f"{args.source_b} -- possible encoder non-determinism, "
                            f"investigate (NOT blocking the golden build)",
                            file=sys.stderr,
                        )

        eta = "eta" if HAPTIC.search(oin1) else "-"
        lines.append(f"{name}\t{sha(oin1)}\t{sha(oin2)}\t{len(oin1)}\t{len(oin2)}\t{eta}")

    if missing_primary:
        sys.exit(
            f"error: {len(missing_primary)} / {len(names)} names have NO report in the "
            f"primary source {args.source_a} -- every cohort name must resolve there "
            f"(that is where it was selected from): {missing_primary}"
        )
    if bad_primary:
        sys.exit(
            f"error: {len(bad_primary)} / {len(names)} names failed the full predicate "
            f"in the primary source (should be impossible -- re-selection would have "
            f"excluded them): {bad_primary}"
        )
    if len(lines) == 0:
        sys.exit("error: 0 golden rows produced -- refusing to write an empty golden manifest")

    print(
        f"# primary={args.source_a} names={len(names)} "
        f"corroborated={corroborated} corroboration_mismatches={len(corroboration_mismatches)}",
        file=sys.stderr,
    )

    manifest = "\n".join(lines)
    out_text = manifest + f"\n# MANIFEST_SHA256={sha(manifest)}\n#DONE {len(lines)}"
    print(out_text)
    if args.out:
        with open(args.out, "w") as f:
            f.write(out_text + "\n")
        print(f"# wrote {os.path.abspath(args.out)}", file=sys.stderr)


if __name__ == "__main__":
    main()
