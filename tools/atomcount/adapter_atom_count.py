"""How many atoms does the *adapter* build from an OIN string?

Runs the generator's own front half -- ``OINParser`` then
``metallogen_adapter._prepare_ligand_fragments`` -- and counts the atoms of the
per-fragment SMILES it hands MetalloGen, hydrogens included.  That number is what
the generated XYZ will contain, so comparing it against the input XYZ separates an
encoder defect (the string already implies the wrong count) from a generator
defect (the string is right and construction lost atoms).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

from oinsmiles.generation.metallogen_adapter import _prepare_ligand_fragments  # noqa: E402
from oinsmiles.generation.oin_parser import OINParser  # noqa: E402


def count_with_hs(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is None:
            return None, None
        try:
            Chem.SanitizeMol(
                m,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
            )
        except Exception:
            return None, None
    try:
        mh = Chem.AddHs(m)
    except Exception:
        return None, None
    heavy = sum(1 for a in mh.GetAtoms() if a.GetAtomicNum() > 1)
    hs = sum(1 for a in mh.GetAtoms() if a.GetAtomicNum() == 1)
    return heavy, hs


def implied(oin: str):
    parsed = OINParser().parse(oin)
    metal_frag, specs, geo = _prepare_ligand_fragments(parsed)
    total = 0
    detail = []
    mh, mhs = count_with_hs(metal_frag)
    if mh is None:
        return None, [("METAL-UNPARSED", metal_frag)]
    total += mh + mhs
    detail.append((metal_frag, mh, mhs))
    for smi, _w in specs:
        h, hs = count_with_hs(smi)
        if h is None:
            return None, detail + [("UNPARSED", smi)]
        total += h + hs
        detail.append((smi, h, hs))
    return total, detail


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("molecules", nargs="+")
    ap.add_argument(
        "--results",
        default="/home/tjmustard/Documents/GitHub/OIN-SMILES/"
        "tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042",
    )
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    R = Path(args.results)

    print(f"{'molecule':22s} {'in':>4s} {'adapter':>8s} {'genXYZ':>7s} {'a-in':>5s}  verdict")
    for name in args.molecules:
        rep = json.loads((R / "individual_reports" / f"{name}.json").read_text())
        oin = rep["smiles_1"]
        inp = Path(rep["input_xyz"])
        n_in = int(inp.read_text().splitlines()[0].split()[0])
        g = R / "structures" / f"{name}_generated.xyz"
        n_gen = int(g.read_text().splitlines()[0].split()[0]) if g.exists() else None
        try:
            n_ad, detail = implied(oin)
        except Exception as e:  # noqa: BLE001
            print(
                f"{name:22s} {n_in:4d} {'ERR':>8s} {str(n_gen):>7s}        {type(e).__name__}: {e}"
            )
            continue
        if n_ad is None:
            print(f"{name:22s} {n_in:4d} {'UNPARSED':>8s} {str(n_gen):>7s}")
            if args.verbose:
                for d in detail:
                    print("      ", d)
            continue
        if n_ad != n_in and n_gen is not None and n_ad == n_gen:
            v = "ENCODER: string already implies the wrong count"
        elif n_ad == n_in and n_gen is not None and n_gen != n_in:
            v = "GENERATOR: string is right, construction lost/gained"
        elif n_ad != n_in and n_gen is None:
            v = "ENCODER (no stored gen to confirm)"
        elif n_ad == n_in and n_gen is None:
            v = "string is right (no stored gen)"
        else:
            v = "three-way different"
        print(f"{name:22s} {n_in:4d} {n_ad:8d} {str(n_gen):>7s} {n_ad - n_in:+5d}  {v}")
        if args.verbose:
            for d in detail:
                print("      ", d)


if __name__ == "__main__":
    main()
