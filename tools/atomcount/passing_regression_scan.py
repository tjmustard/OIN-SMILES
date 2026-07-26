"""Would the atom_count fixes break a molecule that currently PASSES?

A passing round trip means the generated structure had exactly as many atoms as the
input.  So for any passing molecule, the adapter-implied atom count computed with the
CURRENT code must still equal its input atom count -- if it does not, the fix has moved
that molecule from pass to fail.  That makes this a real before/after without needing a
second checkout: the "before" arm is the frozen sweep's own success verdict.

Targets the population Fix B actually touches (fragments carrying a 5-membered
heteroaromatic ring -- thiophene, pyrrole, furan, pyrazole, imidazole), since those are
the rings it stops charging, plus a random control sample.
"""

from __future__ import annotations

import argparse
import collections
import json
import random
import re
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

from oinsmiles.generation import metallogen_adapter as MA  # noqa: E402
from oinsmiles.generation.oin_parser import OINParser  # noqa: E402

RESULTS = Path(
    "/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042"
)
SLOT_RE = re.compile(r"\{\d+[><^]?\}")
#: a 5-membered aromatic ring containing a heteroatom -- exactly what Fix B stops charging
HETERO5 = Chem.MolFromSmarts("[a;r5]1[a;r5][a;r5][a;r5][a;r5]1")


def reparse(s):
    m = Chem.MolFromSmiles(s)
    if m is not None:
        return m
    m = Chem.MolFromSmiles(s, sanitize=False)
    if m is None:
        return None
    try:
        Chem.SanitizeMol(
            m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        )
        return m
    except Exception:
        pass
    m = Chem.MolFromSmiles(s, sanitize=False)
    try:
        m.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(m)
        return m
    except Exception:
        return None


def count(m):
    if m is None:
        return None
    return sum(1 + (0 if a.GetAtomicNum() == 1 else a.GetTotalNumHs()) for a in m.GetAtoms())


def adapter_total(oin: str):
    parsed = OINParser().parse(oin)
    metal_frag, specs, _geo = MA._prepare_ligand_fragments(parsed)
    total = count(reparse(re.sub(r"_[A-Z0-9]+", "", metal_frag)))
    if total is None:
        return None
    for smi, _w in specs:
        t = count(reparse(smi))
        if t is None:
            return None
        total += t
    return total


def has_hetero5(oin: str) -> bool:
    for frag in oin.split("."):
        m = reparse(SLOT_RE.sub("", re.sub(r"_[A-Z0-9]+", "", frag)))
        if m is None:
            continue
        for match in m.GetSubstructMatches(HETERO5):
            if any(m.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in match):
                return True
    return False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit-hetero", type=int, default=60)
    ap.add_argument("--limit-control", type=int, default=40)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    reports = sorted((RESULTS / "individual_reports").glob("*.json"))
    random.seed(args.seed)
    random.shuffle(reports)

    hetero, control = [], []
    for rp in reports:
        if len(hetero) >= args.limit_hetero and len(control) >= args.limit_control:
            break
        try:
            r = json.loads(rp.read_text())
        except Exception:
            continue
        if r.get("status") != "success" or not r.get("smiles_1"):
            continue
        try:
            if has_hetero5(r["smiles_1"]):
                if len(hetero) < args.limit_hetero:
                    hetero.append(r)
            elif len(control) < args.limit_control:
                control.append(r)
        except Exception:
            continue

    verdict: collections.Counter = collections.Counter()
    for label, group in (("hetero5", hetero), ("control", control)):
        print(f"\n=== {label}: {len(group)} currently-PASSING molecules ===")
        for r in group:
            inp = Path(r["input_xyz"])
            if not inp.exists():
                continue
            n_in = int(inp.read_text().splitlines()[0].split()[0])
            try:
                tot = adapter_total(r["smiles_1"])
            except Exception as e:  # noqa: BLE001
                verdict[f"{label}:probe-error"] += 1
                print(f"  {r['molecule']:22s} probe error {type(e).__name__}")
                continue
            if tot is None:
                verdict[f"{label}:unparseable"] += 1
                continue
            if tot == n_in:
                verdict[f"{label}:still-ok"] += 1
            else:
                verdict[f"{label}:REGRESSED"] += 1
                print(
                    f"  {r['molecule']:22s} in={n_in} adapter={tot} d={tot - n_in:+d}  <== REGRESSION"
                )

    print("\n-- verdicts --")
    for k, v in sorted(verdict.items()):
        print(f"  {v:4d}  {k}")
    bad = sum(v for k, v in verdict.items() if "REGRESSED" in k)
    print(f"\nregressions among currently-passing molecules: {bad}")


if __name__ == "__main__":
    main()
