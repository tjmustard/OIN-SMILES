"""Is the axial token invariant to how the molecule's bond orders were PERCEIVED?

The encoder and the generator build their mols by different routes -- ``xyz2mol`` bond-order
perception from interatomic distances on one side, ``build_contract_mol``'s per-fragment
transfer from the OIN fragment SMILES on the other -- and for a metalloporphyrin the two
disagree about the macrocycle. Comparing a requested token against a generated one is only
meaningful if the descriptor does not depend on that disagreement.

This sweep applies a *worst-case* perception perturbation to the same coordinates -- every
non-metal bond flattened to a plain single bond, aromatic flags and formal charges cleared --
and asserts the token is unchanged. Coordinates and connectivity are untouched, so the
molecule's handedness is untouched; only the perception is. Any structure that reports a
MISMATCH has a token the round trip cannot verify.

Run:  PYTHONPATH=$PWD/src python -m tools.injectivity.axial_perception_sweep \
          [--fixtures] [--dataset DIR --n 400]
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import random
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
for _p in (str(REPO), str(REPO / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from rdkit import Chem  # noqa: E402

from oinsmiles.core.constants import TRANSITION_METALS_NUM  # noqa: E402
from oinsmiles.oin.axial import axial_token, detect_axial_axes  # noqa: E402

DEFAULT_DATASET = Path("/home/tjmustard/Documents/GitHub/tmCat-tmPhoto/tmCAT-tmPHOTO_xyz_dataset")
OUT_DIR = REPO / "results-injectivity-y2"
SAMPLE_SEED = 42


@contextlib.contextmanager
def _silence():
    with open(os.devnull, "w") as devnull:
        old = os.dup(2)
        os.dup2(devnull.fileno(), 2)
        try:
            yield
        finally:
            os.dup2(old, 2)
            os.close(old)


def delocalize(mol: Chem.Mol) -> Chem.Mol:
    """Same atoms, same coordinates, same connectivity -- bond orders and aromaticity erased.

    Metal bonds are left alone: they carry the coordination topology, and flattening a dative
    M-donor bond to a covalent single bond would close chelate rings that do not exist,
    changing ring membership rather than merely the perception.
    """
    p = Chem.RWMol(mol)
    for b in p.GetBonds():
        if (
            b.GetBeginAtom().GetAtomicNum() in TRANSITION_METALS_NUM
            or b.GetEndAtom().GetAtomicNum() in TRANSITION_METALS_NUM
        ):
            continue
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in p.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
    out = p.GetMol()
    out.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(out)
    return out


def check(path: Path) -> dict:
    from oinsmiles.utils.xyz2mol import get_tmc_mol

    with _silence():
        mol, _ = get_tmc_mol(path, 0, with_stereo=False)
    base = axial_token(mol)
    pert = axial_token(delocalize(mol))
    axes = detect_axial_axes(mol)
    return {
        "name": path.stem,
        "path": str(path),
        "token": base,
        "perturbed": pert,
        "invariant": base == pert,
        "n_axes": len(axes),
        "n_emitting": sum(1 for a in axes if a.emits),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fixtures", action="store_true", help="sweep tests/fixtures/*.xyz")
    ap.add_argument(
        "--skip",
        nargs="*",
        default=["BENVOG_comp_0"],
        help="fixture stems to skip. BENVOG's macrocycle burns unbounded CPU in "
        "ResonanceMolSupplier (v0.4.4 SL5 gave the ENCODER a forked budget for it; "
        "get_tmc_mol here has no such guard), so it never returns.",
    )
    ap.add_argument("--dataset", type=Path, default=None)
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--out", default="axial_perception_sweep")
    args = ap.parse_args(argv)

    files: list[Path] = []
    if args.fixtures or args.dataset is None:
        files += sorted((REPO / "tests" / "fixtures").glob("*.xyz"))
    if args.dataset is not None:
        ds = args.dataset if args.dataset != Path("default") else DEFAULT_DATASET
        pool = []
        for sub in ("cat", "photo"):
            d = ds / sub
            if d.is_dir():
                pool.extend(sorted(d.glob("*.xyz")))
        rng = random.Random(SAMPLE_SEED)
        files += sorted(rng.sample(pool, min(args.n, len(pool))))
    skip = set(args.skip or ())
    files = [f for f in files if f.stem not in skip]

    rows, failed = [], 0
    for i, p in enumerate(files, 1):
        try:
            r = check(p)
        except Exception as e:
            failed += 1
            print(f"  SKIP {p.stem}: {type(e).__name__}: {e}", flush=True)
            continue
        rows.append(r)
        if not r["invariant"] or r["token"]:
            print(
                f"  {'OK  ' if r['invariant'] else 'MISM'} {r['name']:22} "
                f"token={r['token']!r:10} perturbed={r['perturbed']!r:10} "
                f"axes={r['n_axes']} emitting={r['n_emitting']}",
                flush=True,
            )
        if i % 50 == 0:
            print(f"  ... {i}/{len(files)}", flush=True)

    mism = [r for r in rows if not r["invariant"]]
    emitting = [r for r in rows if r["token"]]
    summary = {
        "n_checked": len(rows),
        "n_load_failures": failed,
        "n_emitting": len(emitting),
        "n_mismatch": len(mism),
        "mismatches": mism,
        "rows": rows,
    }
    OUT_DIR.mkdir(exist_ok=True)
    (OUT_DIR / f"{args.out}.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"\nperception-invariance: {len(rows) - len(mism)}/{len(rows)} invariant "
        f"({len(emitting)} emitting a token); {len(mism)} MISMATCH; {failed} unloadable"
    )
    return 0 if not mism else 1


if __name__ == "__main__":
    raise SystemExit(main())
