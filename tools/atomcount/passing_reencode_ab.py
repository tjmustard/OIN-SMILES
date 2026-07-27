"""Fix A's encoder-side blast radius on molecules that currently PASS.

`passing_regression_scan.py` reuses the OIN strings the frozen sweep recorded, so it
only exercises the adapter half of `OIN_H_FAITHFUL`. This re-encodes each molecule in
both arms, which is the only way to see whether the lever changes a passing molecule's
notation -- and, if it does, whether the adapter still builds the input's atom count.

A changed string is not automatically a regression: the lever only rewrites a fragment
whose serialization was already losing hydrogen. What must not change is the built atom
count for a molecule that already round-tripped.
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import random
import re
import signal
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

from oinsmiles.generation import metallogen_adapter as MA  # noqa: E402
from oinsmiles.generation.oin_parser import OINParser  # noqa: E402
from oinsmiles.oin import hydrogen as HY  # noqa: E402
from oinsmiles.utils.perception_tmc import get_oin_string, get_tmc_mol  # noqa: E402

RESULTS = Path(
    "/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042"
)


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


class TO(Exception):
    pass


def _a(_s, _f):
    raise TO()


def arm(path: str, lever: bool):
    if lever:
        os.environ["OIN_H_FAITHFUL"] = "1"
    else:
        os.environ.pop("OIN_H_FAITHFUL", None)
    assert HY.hydrogen_faithfulness_enabled() is lever
    tmc, coords = get_tmc_mol(path, 0, with_stereo=True)
    oin = get_oin_string(tmc, coords)
    return oin, adapter_total(oin)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=25)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--timeout", type=int, default=150)
    args = ap.parse_args()

    reports = sorted((RESULTS / "individual_reports").glob("*.json"))
    random.seed(args.seed)
    random.shuffle(reports)
    signal.signal(signal.SIGALRM, _a)

    v: collections.Counter = collections.Counter()
    done = 0
    for rp in reports:
        if done >= args.limit:
            break
        try:
            r = json.loads(rp.read_text())
        except Exception:
            continue
        if r.get("status") != "success":
            continue
        inp = Path(r["input_xyz"])
        if not inp.exists():
            continue
        n_in = int(inp.read_text().splitlines()[0].split()[0])
        signal.alarm(args.timeout)
        try:
            oin_off, tot_off = arm(str(inp), False)
            oin_on, tot_on = arm(str(inp), True)
        except TO:
            v["encode-timeout"] += 1
            signal.alarm(0)
            continue
        except Exception as e:  # noqa: BLE001
            v[f"encode-fail:{type(e).__name__}"] += 1
            signal.alarm(0)
            continue
        finally:
            signal.alarm(0)
        done += 1
        changed = oin_off != oin_on
        v["string-changed" if changed else "string-identical"] += 1
        if tot_on is None:
            v["unparseable-on"] += 1
            continue
        if tot_on == n_in:
            v["built-count-still-correct"] += 1
        else:
            v["BUILT-COUNT-REGRESSED"] += 1
            print(f"  {r['molecule']:22s} in={n_in} off={tot_off} on={tot_on}  <== REGRESSION")
            print(f"      off: {oin_off}")
            print(f"      on : {oin_on}")
        if changed:
            print(f"  {r['molecule']:22s} string changed (built count {tot_on} == input {n_in})")
    print(f"\n-- {done} currently-passing molecules re-encoded in both arms --")
    for k, val in sorted(v.items()):
        print(f"  {val:4d}  {k}")


if __name__ == "__main__":
    main()
