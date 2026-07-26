"""Measure serialization H-faithfulness of the encoder, per molecule.

`OINSanitizer.generate_robust_smiles` writes each ligand fragment with
`MolToSmiles`.  For an atom in the SMILES organic subset, RDKit's writer emits a
BARE symbol whenever it judges brackets unnecessary -- and a bare symbol is read
back as "fill to the next allowed valence with hydrogen".  `SetNoImplicit(True)`
does not force a bracket.  So an atom the encoder believes has 0 H can serialize
bare and re-parse with 1 H: a phantom hydrogen, baked into the OIN string, which
the generator then faithfully builds.

This probe re-parses every fragment the encoder emits and reports the per-atom
H delta, so the predicted atom-count delta can be compared against the delta the
round-trip harness actually reported.
"""

from __future__ import annotations

import argparse
import collections
import json
import signal
import sys
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

from oinsmiles.utils import oin_aligner  # noqa: E402
from oinsmiles.utils.xyz2mol import get_oin_string, get_tmc_mol  # noqa: E402

_orig = oin_aligner.OINSanitizer.generate_robust_smiles
FRAGS: list = []


def _reparse(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    if m is not None:
        return m
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    if m is None:
        return None
    try:
        Chem.SanitizeMol(
            m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        )
    except Exception:
        return None
    return m


def patched(ligand_mol, binding_indices_in_ligand):
    out = _orig(ligand_mol, binding_indices_in_ligand)
    smiles, kmol = out if isinstance(out, tuple) else (out, None)
    rec = {"smiles": smiles, "sites": [], "err": None}
    if kmol is not None:
        try:
            order = list(
                kmol.GetPropsAsDict(includePrivate=True, includeComputed=True)[
                    "_smilesAtomOutputOrder"
                ]
            )
        except Exception:
            order = None
        back = _reparse(smiles)
        if back is None:
            rec["err"] = "reparse-failed"
        elif order is None:
            rec["err"] = "no-output-order"
        elif back.GetNumAtoms() != kmol.GetNumAtoms():
            rec["err"] = f"atom-count {kmol.GetNumAtoms()} -> {back.GetNumAtoms()}"
        else:
            binders = set(binding_indices_in_ligand)
            for pos, orig_idx in enumerate(order):
                ia = kmol.GetAtomWithIdx(int(orig_idx))
                ba = back.GetAtomWithIdx(pos)
                want = ia.GetTotalNumHs()
                got = ba.GetTotalNumHs()
                if want != got:
                    rec["sites"].append(
                        {
                            "sym": ia.GetSymbol(),
                            "want": want,
                            "got": got,
                            "delta": got - want,
                            "val": ia.GetTotalValence(),
                            "arom": bool(ia.GetIsAromatic()),
                            "deg": ia.GetDegree(),
                            "chg": ia.GetFormalCharge(),
                            "rad": ia.GetNumRadicalElectrons(),
                            "binder": int(orig_idx) in binders,
                            "ring": ia.IsInRing(),
                        }
                    )
    FRAGS.append(rec)
    return out


oin_aligner.OINSanitizer.generate_robust_smiles = staticmethod(patched)


class TimeoutErr(Exception):
    pass


def _alarm(_s, _f):
    raise TimeoutErr()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worklist", default=None)
    ap.add_argument("--molecules", nargs="*", default=[])
    ap.add_argument(
        "--dataset",
        default="/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset",
    )
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    names = list(args.molecules)
    reported = {}
    if args.worklist:
        wl = json.loads(Path(args.worklist).read_text())["atom_count"]
        import re as _re

        for e in wl:
            names.append(e["molecule"])
            m = _re.search(r"Input (\d+) != Gen (\d+)", e["error"])
            reported[e["molecule"]] = (int(m.group(1)), int(m.group(2)))

    root = Path(args.dataset)
    signal.signal(signal.SIGALRM, _alarm)
    results = []
    for name in names:
        hits = list(root.glob(f"*/{name}.xyz"))
        if not hits:
            results.append({"molecule": name, "status": "no-input"})
            continue
        p = hits[0]
        n_in = int(p.read_text().splitlines()[0].split()[0])
        FRAGS.clear()
        rec = {"molecule": name, "input_n": n_in}
        if name in reported:
            rec["reported"] = reported[name]
            rec["reported_delta"] = reported[name][1] - reported[name][0]
        signal.alarm(args.timeout)
        try:
            tmc, coords = get_tmc_mol(str(p), 0, with_stereo=True)
            oin = get_oin_string(tmc, coords)
            rec["status"] = "ok"
            rec["oin"] = oin
        except TimeoutErr:
            rec["status"] = "encode-timeout"
        except Exception as e:  # noqa: BLE001
            rec["status"] = f"encode-fail: {type(e).__name__}: {e}"
        finally:
            signal.alarm(0)
        sites = [s for f in FRAGS for s in f["sites"]]
        rec["frag_errs"] = [f["err"] for f in FRAGS if f["err"]]
        rec["sites"] = sites
        rec["predicted_delta"] = sum(s["delta"] for s in sites)
        rec["n_phantom"] = sum(s["delta"] for s in sites if s["delta"] > 0)
        rec["n_dropped"] = -sum(s["delta"] for s in sites if s["delta"] < 0)
        results.append(rec)
        print(
            f"{name:22s} {rec['status'][:18]:18s} in={n_in:4d} "
            f"pred_delta={rec['predicted_delta']:+3d} "
            f"reported_delta={rec.get('reported_delta', 'NA'):>4} "
            f"sites={len(sites)} fragerr={len(rec['frag_errs'])}",
            flush=True,
        )

    ok = [r for r in results if r["status"] == "ok"]
    print(f"\n=== {len(ok)}/{len(results)} encoded ===")
    exact = [r for r in ok if r.get("reported_delta") == r["predicted_delta"]]
    print(f"predicted delta EXACTLY equals reported delta: {len(exact)}/{len(ok)}")
    signm = [
        r
        for r in ok
        if r.get("reported_delta") is not None
        and r["predicted_delta"] != 0
        and (r["predicted_delta"] > 0) == (r["reported_delta"] > 0)
    ]
    print(f"same sign, nonzero:                            {len(signm)}/{len(ok)}")
    print(
        f"predicted delta == 0 (unexplained):            {sum(1 for r in ok if r['predicted_delta'] == 0)}/{len(ok)}"
    )

    print("\n-- phantom-H sites by (element, valence, aromatic, binder) --")
    c = collections.Counter()
    for r in ok:
        for s in r["sites"]:
            if s["delta"] > 0:
                c[(s["sym"], s["val"], s["arom"], s["binder"], s["delta"])] += 1
    for k, v in c.most_common(30):
        print(f"  {v:4d}  sym={k[0]:2s} val={k[1]} arom={int(k[2])} binder={int(k[3])} +{k[4]}H")

    print("\n-- dropped-H sites --")
    c2 = collections.Counter()
    for r in ok:
        for s in r["sites"]:
            if s["delta"] < 0:
                c2[(s["sym"], s["val"], s["arom"], s["binder"], s["delta"])] += 1
    for k, v in c2.most_common(30):
        print(f"  {v:4d}  sym={k[0]:2s} val={k[1]} arom={int(k[2])} binder={int(k[3])} {k[4]}H")

    if args.out:
        args.out.write_text(json.dumps(results, indent=1, default=str))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    sys.exit(main())
