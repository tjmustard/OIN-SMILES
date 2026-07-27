"""End-to-end atom-count audit without any 3D work: encode, then ask the adapter.

Encodes the input XYZ with the current code and lever settings, then runs the
resulting OIN string through the generator's front half and counts the atoms
MetalloGen would be asked to build.  Because the adapter-implied count was measured
to equal the harness's reported generated count in 59/59 auditable molecules
(docs/agentic-notes/v0.4.5/ATOM_COUNT_v0.4.5.md Sec 1), this predicts the atom-count gate without paying
for embedding -- which makes it usable as an A/B arm.
"""

from __future__ import annotations

import argparse
import collections
import json
import re
import signal
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

from oinsmiles.generation import metallogen_adapter as MA  # noqa: E402
from oinsmiles.generation.oin_parser import OINParser  # noqa: E402
from oinsmiles.utils.xyz2mol import get_oin_string, get_tmc_mol  # noqa: E402


def reparse(s):
    m = Chem.MolFromSmiles(s)
    if m is None:
        m = Chem.MolFromSmiles(s, sanitize=False)
        if m is None:
            return None
        try:
            Chem.SanitizeMol(
                m,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
            )
        except Exception:
            return None
    return m


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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worklist", required=True)
    ap.add_argument("--key", default="atom_count")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--timeout", type=int, default=150)
    ap.add_argument(
        "--dataset",
        default="/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset",
    )
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    wl = json.loads(Path(args.worklist).read_text())[args.key]
    if args.limit:
        wl = wl[: args.limit]
    root = Path(args.dataset)
    signal.signal(signal.SIGALRM, _a)

    rows = []
    verdict: collections.Counter = collections.Counter()
    for e in wl:
        name = e["molecule"]
        hits = list(root.glob(f"*/{name}.xyz"))
        if not hits:
            continue
        p = hits[0]
        n_in = int(p.read_text().splitlines()[0].split()[0])
        row = {"molecule": name, "input_n": n_in}
        signal.alarm(args.timeout)
        try:
            tmc, coords = get_tmc_mol(str(p), 0, with_stereo=True)
            oin = get_oin_string(tmc, coords)
            row["oin"] = oin
            tot = adapter_total(oin)
            if tot is None:
                row["status"] = "unparseable"
            else:
                row["status"] = "ok"
                row["adapter_n"] = tot
                row["delta"] = tot - n_in
        except TO:
            row["status"] = "timeout"
        except Exception as ex:  # noqa: BLE001
            row["status"] = f"fail: {type(ex).__name__}: {str(ex)[:60]}"
        finally:
            signal.alarm(0)
        v = (
            "MATCH"
            if row.get("delta") == 0
            else (row["status"] if row["status"] != "ok" else f"delta {row['delta']:+d}")
        )
        verdict[v if v.startswith(("MATCH", "timeout", "unparseable")) else "MISMATCH"] += 1
        rows.append(row)
        print(
            f"  {name:22s} in={n_in:4d} {row['status'][:22]:22s} delta={row.get('delta', 'NA')}",
            flush=True,
        )

    print("\n-- summary --")
    for k, v in verdict.most_common():
        print(f"  {v:4d}  {k}")
    if args.out:
        args.out.write_text(json.dumps(rows, indent=1, default=str))


if __name__ == "__main__":
    main()
