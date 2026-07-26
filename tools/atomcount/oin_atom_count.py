"""How many atoms does an OIN string *imply*?  (generator-free)

The `atom_count` hard-fail class is the harness's LAST gate, reached only after
``canonical_roundtrip_key`` already declared the round trip equal -- and that key
deliberately folds implicit-vs-explicit H.  So the class is exactly the set where
the key says "same" and the raw H count says "different".

This probe adjudicates between the two candidate culprits without running the
generator:

  * if the OIN-implied atom count matches the *generated* structure, the string is
    what the generator faithfully built, and the input XYZ is the outlier
    (crystallographic H, or an H the encoder dropped);
  * if it matches the *input*, the generator is at fault.
"""

from __future__ import annotations

import argparse
import collections
import json
import re
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

SLOT_RE = re.compile(r"\{\d+[><^]?\}")
METAL_RE = re.compile(r"\[([A-Z][a-z]?)(?:@[A-Z0-9]+)?_([A-Z]{3})\]")
AXIAL_RE = re.compile(r"\s*\|ax:[+\-]*\|")


def oin_implied_counts(oin: str) -> tuple[collections.Counter, list[str]]:
    """Element histogram implied by an OIN string, H's included.

    Returns (histogram, list of fragments that failed to parse).
    """
    s = AXIAL_RE.sub("", oin.strip())
    s = METAL_RE.sub(r"[\1]", s)
    s = SLOT_RE.sub("", s)
    hist: collections.Counter = collections.Counter()
    bad: list[str] = []
    for frag in s.split("."):
        frag = frag.strip()
        if not frag:
            continue
        mol = Chem.MolFromSmiles(frag)
        if mol is None:
            mol = Chem.MolFromSmiles(frag, sanitize=False)
            if mol is None:
                bad.append(frag)
                continue
            try:
                Chem.SanitizeMol(
                    mol,
                    sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
                )
            except Exception:
                bad.append(frag)
                continue
        try:
            molh = Chem.AddHs(mol)
        except Exception:
            bad.append(frag)
            continue
        for a in molh.GetAtoms():
            hist[a.GetSymbol()] += 1
    return hist, bad


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--classify-json", type=Path, required=True)
    args = ap.parse_args()
    rows = json.loads(args.classify_json.read_text())

    verdict: collections.Counter = collections.Counter()
    print(
        f"{'molecule':22s} {'inXYZ':>6s} {'oin':>6s} {'gen':>6s} "
        f"{'oin-in':>7s} {'oin-gen':>8s}  verdict"
    )
    out = []
    for r in rows:
        hist, bad = oin_implied_counts(r["oin_1"])
        n_oin = sum(hist.values())
        n_in = r.get("input_n")
        n_gen = r.get("gen_n")
        r["oin_n"] = n_oin
        r["oin_hist"] = dict(hist)
        r["oin_unparsed"] = bad
        if bad:
            v = "UNPARSED-FRAG"
        elif n_gen is not None and n_oin == n_gen and n_oin != n_in:
            v = "string==gen (encoder/input)"
        elif n_oin == n_in and n_gen is not None and n_oin != n_gen:
            v = "string==input (generator)"
        elif n_gen is None:
            v = "no-stored-gen: string==input" if n_oin == n_in else "no-stored-gen: string!=input"
        else:
            v = "three-way-different"
        verdict[v] += 1
        r["verdict"] = v
        print(
            f"{r['molecule']:22s} {str(n_in):>6s} {n_oin:6d} {str(n_gen):>6s} "
            f"{n_oin - (n_in or 0):+7d} "
            f"{(n_oin - n_gen) if n_gen is not None else 0:+8d}  {v}"
        )
        out.append(r)

    print("\n-- verdicts --")
    for k, v in verdict.most_common():
        print(f"  {v:3d}  {k}")

    args.classify_json.with_suffix(".oin.json").write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
