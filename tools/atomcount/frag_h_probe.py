"""Dump, per ligand fragment, the H bookkeeping at the moment the encoder freezes it.

`OINSanitizer.generate_robust_smiles` step 1 (`utils/oin_aligner.py:283-289`) freezes
`GetTotalNumHs()` -- explicit + RDKit's implicit-valence *guess* -- as an explicit
bracket count.  This probe shows both halves of that sum separately, which is how you
tell an input-derived H from a phantom one.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

from oinsmiles.utils import oin_aligner  # noqa: E402
from oinsmiles.utils.xyz2mol import get_oin_string, get_tmc_mol  # noqa: E402

_orig = oin_aligner.OINSanitizer.generate_robust_smiles
LOG: list = []


def patched(ligand_mol, binding_indices_in_ligand):
    pre = []
    m = Chem.Mol(ligand_mol)
    try:
        m.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    for a in m.GetAtoms():
        pre.append(
            {
                "idx": a.GetIdx(),
                "sym": a.GetSymbol(),
                "arom": a.GetIsAromatic(),
                "exp": a.GetNumExplicitHs(),
                "imp": a.GetNumImplicitHs(),
                "tot": a.GetTotalNumHs(),
                "val": a.GetTotalValence(),
                "noimp": a.GetNoImplicit(),
                "chg": a.GetFormalCharge(),
                "rad": a.GetNumRadicalElectrons(),
                "deg": a.GetDegree(),
                "binder": a.GetIdx() in set(binding_indices_in_ligand),
            }
        )
    out = _orig(ligand_mol, binding_indices_in_ligand)
    LOG.append({"pre": pre, "smiles": out[0] if isinstance(out, tuple) else out})
    return out


oin_aligner.OINSanitizer.generate_robust_smiles = staticmethod(patched)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("molecules", nargs="+")
    ap.add_argument(
        "--dataset",
        default="/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset",
    )
    ap.add_argument("--only-phantom", action="store_true", help="print only imp>0 atoms")
    args = ap.parse_args()

    for name in args.molecules:
        hits = list(Path(args.dataset).glob(f"*/{name}.xyz"))
        if not hits:
            print(f"{name}: NOT FOUND")
            continue
        LOG.clear()
        tmc, coords = get_tmc_mol(str(hits[0]), 0, with_stereo=True)
        oin = get_oin_string(tmc, coords)
        print(f"\n=== {name} ===\n    oin: {oin}")
        for frag in LOG:
            phantom = [a for a in frag["pre"] if a["imp"] > 0]
            if args.only_phantom and not phantom:
                continue
            print(f"  frag -> {frag['smiles']}")
            for a in frag["pre"]:
                if args.only_phantom and a["imp"] == 0:
                    continue
                flag = ""
                if a["imp"] > 0:
                    flag = (
                        "  <== PHANTOM implicit H frozen by step 1"
                        if a["binder"]
                        else "  <== implicit H (non-binder)"
                    )
                print(
                    f"    {a['idx']:3d} {a['sym']:>2s} arom={int(a['arom'])} deg={a['deg']} "
                    f"exp={a['exp']} imp={a['imp']} tot={a['tot']} val={a['val']} "
                    f"noimp={int(a['noimp'])} chg={a['chg']} rad={a['rad']} "
                    f"binder={int(a['binder'])}{flag}"
                )


if __name__ == "__main__":
    main()
