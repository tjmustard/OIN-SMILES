"""Trace the hydrogen count of one molecule through every encoder/generator stage.

Stages, in order:
  1. input XYZ                       -- ground truth
  2. encoder intent (`kmol`)         -- what OINSanitizer believes, per fragment
  3. emitted fragment SMILES         -- re-parsed, to expose write/read asymmetry
  4. the OIN string's fragments      -- after oin/inline.py's de-bracketing
  5. the adapter's prepared fragment -- what MetalloGen is actually asked to build
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

from oinsmiles.generation.metallogen_adapter import _prepare_ligand_fragments  # noqa: E402
from oinsmiles.generation.oin_parser import OINParser  # noqa: E402
from oinsmiles.utils import oin_aligner  # noqa: E402
from oinsmiles.utils.perception_tmc import get_oin_string, get_tmc_mol  # noqa: E402

_orig = oin_aligner.OINSanitizer.generate_robust_smiles
CALLS: list = []


def _reparse(s):
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


def _totals(m):
    if m is None:
        return None, None
    heavy = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1)
    hs = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 1)
    hs += sum(a.GetTotalNumHs() for a in m.GetAtoms() if a.GetAtomicNum() > 1)
    return heavy, hs


def patched(ligand_mol, binding_indices_in_ligand):
    out = _orig(ligand_mol, binding_indices_in_ligand)
    smiles, kmol = out if isinstance(out, tuple) else (out, None)
    CALLS.append((smiles, kmol, sorted(binding_indices_in_ligand)))
    return out


oin_aligner.OINSanitizer.generate_robust_smiles = staticmethod(patched)
SLOT_RE = re.compile(r"\{\d+[><^]?\}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("molecules", nargs="+")
    ap.add_argument(
        "--results",
        default="/home/tjmustard/Documents/GitHub/OIN-SMILES/"
        "tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042",
    )
    args = ap.parse_args()
    R = Path(args.results)

    for name in args.molecules:
        rep = json.loads((R / "individual_reports" / f"{name}.json").read_text())
        inp = Path(rep["input_xyz"])
        n_in = int(inp.read_text().splitlines()[0].split()[0])
        CALLS.clear()
        tmc, coords = get_tmc_mol(str(inp), 0, with_stereo=True)
        oin = get_oin_string(tmc, coords)
        print(f"\n=== {name}: input {n_in} atoms | reported {rep['error']}")
        print(f"    oin: {oin}")

        seen = {}
        print("\n  [stage 2/3] encoder fragments (intent vs re-parse of its own string)")
        for smiles, kmol, binders in CALLS:
            if smiles in seen:
                continue
            seen[smiles] = True
            ih, ihs = _totals(kmol)
            rb = _reparse(smiles)
            rh, rhs = _totals(rb)
            mark = "" if (ihs == rhs) else f"   <== WRITE/READ H DRIFT {ihs} -> {rhs}"
            print(f"    {smiles:52s} intentH={ihs} reparseH={rhs} heavy={ih}{mark}")

        print("\n  [stage 4/5] OIN fragments -> adapter-prepared fragments")
        parsed = OINParser().parse(oin)
        try:
            metal_frag, specs, geo = _prepare_ligand_fragments(parsed)
        except Exception as e:  # noqa: BLE001
            print(f"    adapter failed: {type(e).__name__}: {e}")
            continue
        tot = 0
        mh, mhs = _totals(_reparse(metal_frag))
        tot += (mh or 0) + (mhs or 0)
        print(f"    metal {metal_frag:46s} heavy={mh} H={mhs}")
        for i, frag in enumerate(parsed.fragments):
            if i == parsed.metal_fragment_idx:
                continue
            bare = SLOT_RE.sub("", frag)
            b = _reparse(bare)
            bh, bhs = _totals(b)
            print(f"    oin-frag  {frag:48s} (slots stripped) heavy={bh} H={bhs}")
        for smi, _w in specs:
            m = _reparse(smi)
            h, hs = _totals(m)
            tot += (h or 0) + (hs or 0)
            print(f"    adapter   {smi:48s} heavy={h} H={hs}")
        print(f"\n    adapter total = {tot}   input = {n_in}   delta = {tot - n_in:+d}")


if __name__ == "__main__":
    main()
