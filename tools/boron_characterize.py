#!/usr/bin/env python
"""Boron-cage spike: characterize the 34 `boron_cluster` encode_fail molecules.

Question 1 of the spike: are these all genuinely 3c-2e cages, and *where exactly*
does each one fail (obabel AC / AC2BO / SanitizeMol / MolToSmiles)?

Every molecule runs in an isolated subprocess, because ``AC2BO`` calls
``sys.exit()`` when an atom's AC degree exceeds every entry in its
``atomic_valence`` list -- which is precisely what a cage boron vertex does.
A ``sys.exit`` cannot be caught by ``except Exception``, so an in-process loop
would die on the first molecule.

Usage:
    PYTHONPATH=src python tools/boron_characterize.py --dataset-dir <abs path>
    PYTHONPATH=src python tools/boron_characterize.py --worker MOL --dataset-dir <...>
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import traceback

BORON_COHORT = [
    "AVOFIB",
    "BEKLUA",
    "BEKMIP",
    "CAKBEW",
    "CAKBOG",
    "COZCEZ",
    "GANYEZ",
    "GOHWOQ",
    "HAXJAS",
    "HAXJOG",
    "ICEZIC",
    "JABGAX",
    "JAFMIP",
    "JAFTAO",
    "JAFTES",
    "MAFSIY",
    "MODZUA",
    "OZAREO",
    "PAQBOZ",
    "PAQCAM",
    "PAYTUH",
    "RANCIU",
    "RANMUR",
    "RAWJEG",
    "RIWKAK",
    "RIWKEO",
    "RONPES",
    "RONQET",
    "RONQOD",
    "RULBUV",
    "ULODUU",
    "XUKRIF",
    "YIBZIV",
    "YIVLAQ",
]

PER_MOL_TIMEOUT_S = 180


def find_xyz(mol: str, dataset_dir: str) -> str | None:
    for sub in ("cat", "photo"):
        p = os.path.join(dataset_dir, sub, f"{mol}_comp_0.xyz")
        if os.path.exists(p):
            return p
    return None


def characterize(mol_name: str, dataset_dir: str) -> dict:
    """Walk the encoder's own perception stages and record the first failure site."""
    import warnings

    warnings.filterwarnings("ignore")
    from collections import Counter

    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdmolops
    from rdkit.Chem.MolStandardize import rdMolStandardize

    RDLogger.DisableLog("rdApp.*")

    from oinsmiles.utils.xyz2mol import (
        TRANSITION_METALS_NUM,
        MetalNon_Hg,
        get_basic_mol,
        get_lig_mol,
        get_proposed_ligand_charge,
        params,
    )

    out: dict = {"mol": mol_name, "stage": None}
    path = find_xyz(mol_name, dataset_dir)
    if path is None:
        out["stage"] = "NO_FILE"
        return out
    out["path"] = path

    # --- raw composition -------------------------------------------------------
    with open(path) as fh:
        lines = fh.read().splitlines()
    n_declared = int(lines[0].split()[0])
    syms = [ln.split()[0] for ln in lines[2 : 2 + n_declared] if ln.split()]
    out["natoms"] = len(syms)
    out["formula"] = "".join(f"{el}{n}" for el, n in sorted(Counter(syms).items()))
    out["nB_raw"] = sum(1 for s in syms if s == "B")

    # --- stage 1: obabel adjacency + basic mol ---------------------------------
    try:
        basic, _xyz = get_basic_mol(path, 0)
    except BaseException as e:  # noqa: BLE001
        out["stage"] = "xyz2AC_obabel"
        out["err"] = f"{type(e).__name__}: {e}"
        return out
    out["stage_obabel"] = "OK"

    for a in basic.GetAtoms():
        a.SetIntProp("__origIdx", a.GetIdx())
    AC_full = Chem.rdmolops.GetAdjacencyMatrix(basic)

    metal_idx = [a.GetIdx() for a in basic.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM]
    out["metals"] = [basic.GetAtomWithIdx(i).GetSymbol() for i in metal_idx]

    # --- stage 1b: boron connectivity in the FULL basic mol, pre-disconnection.
    # Distinguishes "obabel's AC never saw the cage edges" from "the metal
    # disconnector fragmented an intact cage".
    b_idx = [a.GetIdx() for a in basic.GetAtoms() if a.GetAtomicNum() == 5]
    full_bdeg = []
    full_bb = 0
    for bi in b_idx:
        a = basic.GetAtomWithIdx(bi)
        full_bdeg.append(a.GetDegree())
    for bnd in basic.GetBonds():
        if bnd.GetBeginAtom().GetAtomicNum() == 5 and bnd.GetEndAtom().GetAtomicNum() == 5:
            full_bb += 1
    out["full_B_degrees"] = sorted(full_bdeg)
    out["full_nBB"] = full_bb
    out["full_maxBdeg"] = max(full_bdeg) if full_bdeg else 0
    # how many B carry an H in the raw file (cage vertices are B-H)
    out["full_B_H"] = sum(
        1
        for bi in b_idx
        if any(n.GetAtomicNum() == 1 for n in basic.GetAtomWithIdx(bi).GetNeighbors())
    )
    # boron-metal edges
    out["full_B_metal_edges"] = sum(1 for bi in b_idx for mi in metal_idx if AC_full[bi, mi])
    # connected components restricted to the B/C heavy-atom cage skeleton
    out["n_frags_basic"] = len(rdmolops.GetMolFrags(basic))

    # --- stage 2: metal disconnect + fragment inventory ------------------------
    try:
        mdis = rdMolStandardize.MetalDisconnector(params)
        mdis.SetMetalNon(Chem.MolFromSmarts(MetalNon_Hg))
        frags = mdis.Disconnect(basic)
        frag_mols = rdmolops.GetMolFrags(frags, asMols=True)
    except BaseException as e:  # noqa: BLE001
        out["stage"] = "MetalDisconnector"
        out["err"] = f"{type(e).__name__}: {e}"
        return out

    frag_recs = []
    for i, f in enumerate(frag_mols):
        m = Chem.Mol(f)
        borons = [a for a in m.GetAtoms() if a.GetAtomicNum() == 5]
        if not borons:
            continue
        nbb = 0
        for bnd in m.GetBonds():
            if bnd.GetBeginAtom().GetAtomicNum() == 5 and bnd.GetEndAtom().GetAtomicNum() == 5:
                nbb += 1
        # cage-vertex degrees: heavy-atom degree of each boron in the FRAGMENT
        bdeg = sorted(a.GetDegree() for a in borons)
        # carbon atoms bonded to >=2 boron (carborane cage carbons)
        cage_C = [
            a.GetIdx()
            for a in m.GetAtoms()
            if a.GetAtomicNum() == 6
            and sum(1 for n in a.GetNeighbors() if n.GetAtomicNum() == 5) >= 2
        ]
        # is any boron bonded to the metal (metallaborane)?
        b_orig = {a.GetIntProp("__origIdx") for a in borons}
        b_metal = any(AC_full[bi, mi] for bi in b_orig for mi in metal_idx)
        frag_rec = {
            "frag": i,
            "natoms": m.GetNumAtoms(),
            "nB": len(borons),
            "nBB": nbb,
            "B_degrees": bdeg,
            "maxBdeg": max(bdeg),
            "n_cage_C": len(cage_C),
            "B_bound_to_metal": bool(b_metal),
        }
        try:
            frag_rec["smiles_nosan"] = Chem.MolToSmiles(m)
        except BaseException as e:  # noqa: BLE001
            frag_rec["smiles_nosan"] = f"<{type(e).__name__}>"

        # --- stage 3: get_lig_mol on this boron fragment ----------------------
        lig_coord = [
            a.GetIdx()
            for a in m.GetAtoms()
            if any(AC_full[a.GetIntProp("__origIdx"), mi] for mi in metal_idx)
        ]
        frag_rec["n_donors"] = len(lig_coord)
        try:
            proposed = get_proposed_ligand_charge(f)
        except BaseException as e:  # noqa: BLE001
            proposed = None
            frag_rec["huckel_err"] = f"{type(e).__name__}: {e}"
        frag_rec["proposed_charge"] = proposed
        try:
            lm, lc = get_lig_mol(Chem.Mol(m), proposed if proposed is not None else 0, lig_coord)
            frag_rec["get_lig_mol"] = "None" if lm is None else "OK"
            frag_rec["final_charge"] = lc
        except SystemExit as e:
            frag_rec["get_lig_mol"] = f"SYS_EXIT({e.code})"
        except BaseException as e:  # noqa: BLE001
            frag_rec["get_lig_mol"] = f"{type(e).__name__}: {str(e)[:120]}"
        frag_recs.append(frag_rec)

    out["boron_frags"] = frag_recs
    out["n_boron_frags"] = len(frag_recs)

    # --- stage 4: the real encoder, to record the honest terminus -------------
    try:
        from oinsmiles import XYZToSMILES

        oin = XYZToSMILES().convert(path)
        out["stage"] = "ENCODES"
        out["oin"] = oin
    except SystemExit as e:
        out["stage"] = "SYS_EXIT"
        out["err"] = f"sys.exit({e.code})"
    except BaseException as e:  # noqa: BLE001
        tb = traceback.extract_tb(sys.exc_info()[2])
        last = None
        for fr in tb:
            if "oinsmiles" in fr.filename:
                last = fr
        out["stage"] = "CONVERT_FAIL"
        out["etype"] = type(e).__name__
        out["loc"] = f"{os.path.basename(last.filename)}:{last.lineno}:{last.name}" if last else "?"
        out["err"] = str(e).replace("\n", " ")[:300]
    return out


def drive(dataset_dir: str, out_path: str, only: list[str] | None) -> None:
    todo = only or BORON_COHORT
    results = []
    for mol in todo:
        cmd = [
            sys.executable,
            os.path.abspath(__file__),
            "--worker",
            mol,
            "--dataset-dir",
            dataset_dir,
        ]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=PER_MOL_TIMEOUT_S)
            line = next((ln for ln in proc.stdout.splitlines() if ln.startswith("{")), None)
            if line:
                rec = json.loads(line)
            else:
                rec = {
                    "mol": mol,
                    "stage": "WORKER_DIED",
                    "rc": proc.returncode,
                    "err": (proc.stderr or "")[-300:],
                }
        except subprocess.TimeoutExpired:
            rec = {"mol": mol, "stage": "TIMEOUT"}
        results.append(rec)
        bf = rec.get("boron_frags") or []
        summary = ";".join(
            f"nB={f['nB']},nBB={f['nBB']},maxdeg={f['maxBdeg']},lig={f.get('get_lig_mol')}"
            for f in bf
        )
        print(
            f"{mol:10s} {rec.get('formula', '?'):28s} {rec.get('stage'):14s} {summary}",
            flush=True,
        )
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nwrote {out_path}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", default=None)
    ap.add_argument("--dataset-dir", required=True)
    ap.add_argument("--out", default="tools/boron_characterize.json")
    ap.add_argument("--only", default=None, help="comma-separated subset")
    args = ap.parse_args()
    if args.worker:
        print(json.dumps(characterize(args.worker, args.dataset_dir)), flush=True)
    else:
        drive(
            args.dataset_dir,
            args.out,
            args.only.split(",") if args.only else None,
        )


if __name__ == "__main__":
    main()
