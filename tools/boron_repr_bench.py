#!/usr/bin/env python
"""Can an intact boron cage be represented, serialized and RE-PARSED at all?

Builds the cage graph straight from the *unpruned* covalent-radius adjacency
matrix (which `tools/boron_ac_probe.py` showed recovers textbook-exact
topologies), then benches several representations:

  A  plain SINGLE bonds, full sanitize                     (expected: valence error)
  B  plain SINGLE bonds, SANITIZE_ALL ^ SANITIZE_PROPERTIES
  C  plain SINGLE bonds, SANITIZE_ALL ^ PROPERTIES ^ KEKULIZE
  D  cage B-B bonds as DATIVE, full sanitize
  E  cage B-B bonds as UNSPECIFIED (zero-order), full sanitize
  F  plain SINGLE + noImplicit/numExplicitHs pinned, then full sanitize

For each: does it sanitize, does it serialize, and does the SMILES **re-parse to
the same graph** (atom count, element multiset, bond multiset)?  Round-trippability
is the requirement, not chemical elegance.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np


def find_xyz(mol: str, dataset_dir: str) -> str | None:
    for sub in ("cat", "photo"):
        p = os.path.join(dataset_dir, sub, f"{mol}_comp_0.xyz")
        if os.path.exists(p):
            return p
    return None


def raw_ac(atoms, xyz, tolerance=0.5):
    from rdkit import Chem

    pt = Chem.GetPeriodicTable()
    n = len(atoms)
    coords = np.asarray(xyz, dtype=float)
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    rcov = np.array([pt.GetRcovalent(int(z)) for z in atoms])
    thr = rcov[:, None] + rcov[None, :] + tolerance
    return ((d <= thr) & ~np.eye(n, dtype=bool)).astype(int)


def graph_fingerprint(m):
    """Order-independent graph identity: element multiset + bonded-element-pair multiset."""
    from collections import Counter

    els = Counter(a.GetSymbol() for a in m.GetAtoms())
    bonds = Counter()
    for b in m.GetBonds():
        pair = tuple(sorted((b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol())))
        bonds[(pair, str(b.GetBondType()))] += 1
    return (
        m.GetNumAtoms(),
        tuple(sorted(els.items())),
        tuple(sorted((str(k), v) for k, v in bonds.items())),
    )


def build_cage_frag(mol_name, dataset_dir, tolerance=0.5):
    """Extract the largest boron-cage connected component from the UNPRUNED AC.

    Metal atoms are excluded, so what comes back is the cage ligand as the
    geometry actually presents it.
    """
    from rdkit import Chem

    from oinsmiles.utils.perception_core import read_xyz_file
    from oinsmiles.utils.perception_tmc import TRANSITION_METALS_NUM

    path = find_xyz(mol_name, dataset_dir)
    atoms, _c, xyz = read_xyz_file(path)
    AC = raw_ac(atoms, xyz, tolerance)
    n = len(atoms)
    metal = {i for i, z in enumerate(atoms) if z in TRANSITION_METALS_NUM}

    # connected components over non-metal atoms
    seen, comps = set(), []
    for s in range(n):
        if s in seen or s in metal:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in range(n):
                if AC[u, v] and v not in seen and v not in metal:
                    seen.add(v)
                    stack.append(v)
        comps.append(comp)
    # pick the component with the most borons
    comps.sort(key=lambda c: -sum(1 for i in c if atoms[i] == 5))
    comp = sorted(comps[0])
    remap = {g: k for k, g in enumerate(comp)}

    rw = Chem.RWMol()
    for g in comp:
        a = Chem.Atom(int(atoms[g]))
        a.SetNoImplicit(True)
        rw.AddAtom(a)
    edges = []
    for ii, g in enumerate(comp):
        for h in comp[ii + 1 :]:
            if AC[g, h]:
                edges.append((remap[g], remap[h]))
    return rw, edges, [int(atoms[g]) for g in comp]


BB = 5  # boron atomic number


def make_mol(rw_proto, edges, elems, bond_for_bb, pin_hs=False):
    from rdkit import Chem

    rw = Chem.RWMol()
    for z in elems:
        a = Chem.Atom(z)
        a.SetNoImplicit(True)
        rw.AddAtom(a)
    for i, j in edges:
        bt = Chem.BondType.SINGLE
        if elems[i] == BB and elems[j] == BB:
            bt = bond_for_bb
        rw.AddBond(i, j, bt)
    m = rw.GetMol()
    if pin_hs:
        for a in m.GetAtoms():
            a.SetNumExplicitHs(0)
            a.SetNoImplicit(True)
    return m


def try_repr(label, m, sanitize_ops, dataset_note=""):
    from rdkit import Chem

    rec = {"repr": label, "note": dataset_note}
    m = Chem.Mol(m)
    rec["n_atoms"] = m.GetNumAtoms()
    rec["n_bonds"] = m.GetNumBonds()
    try:
        if sanitize_ops is None:
            Chem.SanitizeMol(m)
        else:
            Chem.SanitizeMol(m, sanitizeOps=sanitize_ops)
        rec["sanitize"] = "OK"
    except Exception as e:  # noqa: BLE001
        rec["sanitize"] = f"{type(e).__name__}: {str(e)[:110]}"
        return rec
    try:
        smi = Chem.MolToSmiles(m)
        rec["smiles"] = smi
        rec["smiles_len"] = len(smi)
    except Exception as e:  # noqa: BLE001
        rec["serialize"] = f"{type(e).__name__}: {str(e)[:110]}"
        return rec
    rec["serialize"] = "OK"

    fp_before = graph_fingerprint(m)
    # re-parse, unsanitized (the only honest way to read back a graph the
    # valence rule rejects)
    m2 = Chem.MolFromSmiles(smi, sanitize=False)
    if m2 is None:
        rec["reparse"] = "MolFromSmiles returned None"
        return rec
    try:
        if sanitize_ops is None:
            Chem.SanitizeMol(m2)
        else:
            Chem.SanitizeMol(m2, sanitizeOps=sanitize_ops)
        rec["reparse"] = "OK"
    except Exception as e:  # noqa: BLE001
        rec["reparse"] = f"resanitize {type(e).__name__}: {str(e)[:80]}"
    fp_after = graph_fingerprint(m2)
    rec["roundtrip_graph_equal"] = fp_before == fp_after
    rec["atoms_before_after"] = [fp_before[0], fp_after[0]]
    if fp_before != fp_after:
        rec["fp_before"] = str(fp_before)[:300]
        rec["fp_after"] = str(fp_after)[:300]
    # canonical SMILES idempotence: re-serialize the re-parsed mol
    try:
        smi2 = Chem.MolToSmiles(m2)
        rec["smiles_idempotent"] = smi2 == smi
        if smi2 != smi:
            rec["smiles2"] = smi2
    except Exception as e:  # noqa: BLE001
        rec["smiles_idempotent"] = f"{type(e).__name__}"
    return rec


def bench(mol_name, dataset_dir, tolerance=0.5):
    import warnings

    warnings.filterwarnings("ignore")
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")

    rw, edges, elems = build_cage_frag(mol_name, dataset_dir, tolerance)
    nB = sum(1 for z in elems if z == BB)
    nBB = sum(1 for i, j in edges if elems[i] == BB and elems[j] == BB)
    degs = {}
    for i, j in edges:
        degs[i] = degs.get(i, 0) + 1
        degs[j] = degs.get(j, 0) + 1
    maxBdeg = max((degs.get(i, 0) for i, z in enumerate(elems) if z == BB), default=0)
    head = {
        "mol": mol_name,
        "frag_atoms": len(elems),
        "frag_nB": nB,
        "frag_nBB": nBB,
        "frag_maxBdeg": maxBdeg,
    }

    ALL = Chem.SanitizeFlags.SANITIZE_ALL
    PROPS = Chem.SanitizeFlags.SANITIZE_PROPERTIES
    KEK = Chem.SanitizeFlags.SANITIZE_KEKULIZE

    cases = [
        ("A single / full sanitize", make_mol(rw, edges, elems, Chem.BondType.SINGLE), None),
        (
            "B single / ^PROPERTIES",
            make_mol(rw, edges, elems, Chem.BondType.SINGLE),
            ALL ^ PROPS,
        ),
        (
            "C single / ^PROPERTIES ^KEKULIZE",
            make_mol(rw, edges, elems, Chem.BondType.SINGLE),
            ALL ^ PROPS ^ KEK,
        ),
        ("D B-B DATIVE / full", make_mol(rw, edges, elems, Chem.BondType.DATIVE), None),
        (
            "E B-B UNSPECIFIED / full",
            make_mol(rw, edges, elems, Chem.BondType.UNSPECIFIED),
            None,
        ),
        (
            "F single+pinned H / full",
            make_mol(rw, edges, elems, Chem.BondType.SINGLE, pin_hs=True),
            None,
        ),
        (
            "G B-B ZERO / full",
            make_mol(rw, edges, elems, Chem.BondType.ZERO),
            None,
        ),
        (
            "H B-B ZERO / ^PROPERTIES",
            make_mol(rw, edges, elems, Chem.BondType.ZERO),
            ALL ^ PROPS,
        ),
    ]
    out = [try_repr(lbl, m, ops) for lbl, m, ops in cases]
    return head, out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset-dir", required=True)
    ap.add_argument("--mols", default="AVOFIB,OZAREO,CAKBEW,GANYEZ,MODZUA")
    ap.add_argument("--out", default="tools/boron_repr_bench.json")
    args = ap.parse_args()
    allres = []
    for m in args.mols.split(","):
        head, recs = bench(m, args.dataset_dir)
        print(
            f"\n=== {head['mol']}  cage frag: {head['frag_atoms']} atoms, "
            f"nB={head['frag_nB']}, B-B={head['frag_nBB']}, maxBdeg={head['frag_maxBdeg']}"
        )
        for r in recs:
            print(
                f"  {r['repr']:34s} san={r.get('sanitize', '-')[:38]:38s} "
                f"ser={r.get('serialize', '-')[:20]:20s} "
                f"rt={r.get('reparse', '-')[:26]:26s} "
                f"graph_eq={r.get('roundtrip_graph_equal')} "
                f"idem={r.get('smiles_idempotent')}"
            )
        allres.append({"head": head, "results": recs})
    with open(args.out, "w") as fh:
        json.dump(allres, fh, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main()
