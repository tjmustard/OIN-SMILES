"""Conformation-robust chirality oracle for eta (metallocene-like) complexes.

WHY THIS EXISTS
===============
``tools/injectivity/oracle.py::is_distinct_enantiomer`` mirrors the whole structure and
asks for the best proper-rotation superposition. On a *rigid* fixture that is exactly
right, and it separated the Y1 fixtures cleanly (achiral ~0.05 A, chiral 3-4 A). On a
**flexible** molecule it is useless: an ansa-metallocene's pendant aryls, tBu rotors and
methoxy groups sit at whatever torsions the crystal chose, so even a genuinely *achiral*
meso complex fails to superimpose on its own mirror image. Run on the Lane 3 cases it
reports "chiral, RMSD 2.6-4.4 A" for every single one, cap-hitting the automorphism
enumerator -- i.e. it measures **conformational** chirality and cannot see the
**configurational** question we actually need answered.

WHAT THIS DOES INSTEAD
======================
The stereochemistry an eta winding marker encodes lives entirely in a *rigid* core:
the metal, the coordinated ring systems, the ansa bridge, and the first atom of each
ring substituent. Substituents beyond that first shell are rotors and carry no
configurational information -- but they cannot simply be deleted either, because the
first shell is what breaks a ring's internal reversal symmetry (strip an
indenyl's 2-methyl and 4-aryl and the bare indenyl regains the in-plane C2 that makes
its two faces equivalent).

So: build the core = metal + every fused ring system containing a metal-coordinated ring
atom + the atoms bridging two such systems + one shell of substituent atoms; then run the
mirror / proper-rotation / graph-automorphism test on that. First-shell atoms on a planar
aromatic ring are themselves in-plane, so the core is rigid and the RMSD is meaningful.

Verdict: ``core_rmsd`` near zero => the configuration is **achiral** (its mirror is the
same compound, so the two mirror-related OIN spellings must be canonicalized to one).
Large => **chiral** (the two spellings are genuine enantiomers and folding them would
destroy stereochemistry).

Usage:
    PYTHONPATH=src python tools/eta_core_chirality.py --xyz A.xyz [--xyz B.xyz]
    PYTHONPATH=src python tools/eta_core_chirality.py --results-dir DIR --only MOL1,MOL2
"""

import argparse
import contextlib
import glob
import json
import os
import sys

import numpy as np
from rdkit import Chem

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from oinsmiles.core.constants import TRANSITION_METALS_NUM  # noqa: E402
from tools.injectivity.oracle import _kabsch_proper_rmsd  # noqa: E402

#: A ring atom is "coordinated" if it lies within this distance of the metal (A).
ETA_CUTOFF = 2.95

#: Below this core mirror-RMSD the configuration is called achiral.
CORE_ACHIRAL_RMSD = 0.35


@contextlib.contextmanager
def _silence_fds():
    with open(os.devnull, "w") as devnull:
        old_out, old_err = os.dup(1), os.dup(2)
        try:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            os.dup2(old_out, 1)
            os.dup2(old_err, 2)
            os.close(old_out)
            os.close(old_err)


def _metal_free_view(mol, metal_idx):
    """A copy of ``mol`` with the metal's bonds cut, for through-LIGAND path searches.

    ``Chem.GetShortestPath`` returns only THE shortest path, and between two eta rings
    that is almost always the 2-bond hop via the metal -- which would make every
    ansa-metallocene look unbridged. Cutting the metal's bonds forces the search through
    the covalent bridge. Atom indices are preserved.
    """
    rw = Chem.RWMol(mol)
    for nb in [a.GetIdx() for a in mol.GetAtomWithIdx(int(metal_idx)).GetNeighbors()]:
        if rw.GetBondBetweenAtoms(int(metal_idx), int(nb)) is not None:
            rw.RemoveBond(int(metal_idx), int(nb))
    out = rw.GetMol()
    try:
        out.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(out)
    except Exception:
        pass
    return out


def _fused_system(mol, seed_atoms):
    """Every atom reachable from ``seed_atoms`` through ring bonds (the fused system)."""
    ri = mol.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings()]
    out = set()
    frontier = set(seed_atoms)
    while frontier:
        a = frontier.pop()
        if a in out:
            continue
        out.add(a)
        for r in rings:
            if a in r:
                frontier |= r - out
    return out


def build_core(mol, coords):
    """Atom indices of the rigid configurational core, or None if there is no eta ring.

    core = metal + fused ring systems bearing a coordinated ring atom + bridge atoms
    joining two such systems + one shell of substituent atoms.
    """
    metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM]
    if not metals:
        return None
    m = metals[0]
    ri = mol.GetRingInfo()

    eta_seeds = [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetIdx() != m
        and ri.NumAtomRings(a.GetIdx()) > 0
        and float(np.linalg.norm(coords[a.GetIdx()] - coords[m])) <= ETA_CUTOFF
    ]
    if not eta_seeds:
        return None

    # One fused system per connected group of seeds.
    systems = []
    unassigned = set(eta_seeds)
    while unassigned:
        seed = next(iter(unassigned))
        sysatoms = _fused_system(mol, [seed])
        systems.append(sysatoms)
        unassigned -= sysatoms

    core = set(systems[0]).union(*systems[1:]) if systems else set()
    core.add(m)

    # Bridge: the shortest through-bond path between two eta systems, metal excluded.
    if len(systems) >= 2:
        for i in range(len(systems)):
            for j in range(i + 1, len(systems)):
                best = None
                for a in systems[i]:
                    for b in systems[j]:
                        try:
                            path = Chem.GetShortestPath(mol, int(a), int(b))
                        except Exception:
                            continue
                        if not path or m in path:
                            continue
                        if best is None or len(path) < len(best):
                            best = path
                if best:
                    core |= set(best)

    # One shell of substituent atoms: enough to break a ring's reversal symmetry,
    # not enough to admit a rotor.
    shell = set()
    for idx in list(core):
        for nb in mol.GetAtomWithIdx(int(idx)).GetNeighbors():
            if nb.GetAtomicNum() != 1:
                shell.add(nb.GetIdx())
    core |= shell
    return sorted(core), m


def core_chirality(mol, coords, core_idx, metal_idx, *, max_autos=50000):
    """(mirror RMSD over proper rotations + core automorphisms, n_autos, cap_hit)."""
    amap = {int(a): i for i, a in enumerate(core_idx)}
    sub = Chem.RWMol()
    conf_pos = []
    for a in core_idx:
        at = mol.GetAtomWithIdx(int(a))
        new = Chem.Atom(at.GetAtomicNum())
        # Formal charge is deliberately NOT copied, for the same reason bond order is
        # flattened below: perception localizes a Cp anion's -1 onto ONE arbitrary ring
        # carbon, which collapses ferrocene's 200 core automorphisms to 8 and makes the
        # mirror test report the achiral molecule as chiral at 0.81 A. Charge placement
        # on a delocalized anion is a Kekule-style artifact, not a structural feature.
        new.SetNoImplicit(True)
        sub.AddAtom(new)
        conf_pos.append(coords[int(a)])
    seen_bonds = set()
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in amap and j in amap:
            # Bond ORDER is deliberately flattened to SINGLE: an arbitrary Kekule
            # choice must not make two automorphic rings look inequivalent.
            sub.AddBond(amap[i], amap[j], Chem.BondType.SINGLE)
            seen_bonds.add(frozenset((amap[i], amap[j])))

    # Metal->eta bonds are usually NOT explicit in the perceived mol (haptic
    # coordination has no single sigma bond). Without them the core graph is
    # disconnected and its automorphism group is badly under-enumerated -- which
    # inflates the mirror RMSD and would call ferrocene chiral. Add one bond from the
    # metal to every coordinated ring atom so the automorphism group is the
    # coordination-aware one.
    for a in core_idx:
        if a == metal_idx:
            continue
        if float(np.linalg.norm(coords[int(a)] - coords[metal_idx])) <= ETA_CUTOFF:
            key = frozenset((amap[metal_idx], amap[int(a)]))
            if key not in seen_bonds:
                sub.AddBond(amap[metal_idx], amap[int(a)], Chem.BondType.SINGLE)
                seen_bonds.add(key)
    core = sub.GetMol()
    try:
        core.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    P = np.asarray(conf_pos, dtype=float)
    mirror = P.copy()
    mirror[:, 2] *= -1.0
    n = core.GetNumAtoms()
    matches = core.GetSubstructMatches(
        core, uniquify=False, useChirality=False, maxMatches=max_autos
    )
    best = np.inf
    for mt in matches:
        if len(mt) != n:
            continue
        perm = np.fromiter(mt, dtype=int, count=n)
        best = min(best, _kabsch_proper_rmsd(P, mirror[perm]))
    if best is np.inf:
        best = _kabsch_proper_rmsd(P, mirror)
    return float(best), len(matches), len(matches) >= max_autos


def core_is_rigid(mol, coords, core_idx, metal_idx):
    """(is_rigid, reason). A non-rigid core makes the mirror RMSD meaningless.

    Two ways the core moves even though its *configuration* is fixed:

    * **free ring rotation** -- an eta ring not tied to a second one by a covalent bridge
      spins about the metal->centroid axis, so a mirror image differs by an azimuth that
      no graph automorphism can absorb (a substituted ring's rotations are not graph
      automorphisms);
    * **a flexible bridge** -- an sp3 chain joining two eta systems has real torsions.

    Both are *conformational*, not configurational, so the honest answer in either case is
    UNDECIDABLE rather than a large RMSD dressed up as "chiral".
    """
    ri = mol.GetRingInfo()
    eta_atoms = [
        int(a)
        for a in core_idx
        if a != metal_idx
        and ri.NumAtomRings(int(a)) > 0
        and float(np.linalg.norm(coords[int(a)] - coords[metal_idx])) <= ETA_CUTOFF
    ]
    systems = []
    unassigned = set(eta_atoms)
    while unassigned:
        s = _fused_system(mol, [next(iter(unassigned))])
        systems.append(s)
        unassigned -= s
    if len(systems) < 2:
        return True, ""

    # A ring whose azimuthal rotation IS a graph automorphism (an unsubstituted Cp)
    # costs nothing -- the automorphism enumeration already absorbs it. Only a
    # SUBSTITUTED ring's rotation is a genuine unquenched degree of freedom.
    def _substituted(system):
        return any(
            nb.GetIdx() not in system and nb.GetAtomicNum() not in (1, *TRANSITION_METALS_NUM)
            for a in system
            for nb in mol.GetAtomWithIdx(int(a)).GetNeighbors()
        )

    free = _metal_free_view(mol, metal_idx)
    core_set = {int(a) for a in core_idx}
    for i in range(len(systems)):
        for j in range(i + 1, len(systems)):
            path = None
            for a in systems[i]:
                for b in systems[j]:
                    try:
                        p = Chem.GetShortestPath(free, int(a), int(b))
                    except Exception:
                        continue
                    if p and (path is None or len(p) < len(path)):
                        path = p
            if path is None:
                if _substituted(systems[i]) or _substituted(systems[j]):
                    return False, "unbridged substituted eta rings (free ring rotation)"
                continue
            rotatable = sum(
                1
                for k in range(len(path) - 1)
                if not mol.GetBondBetweenAtoms(int(path[k]), int(path[k + 1])).IsInRing()
                and path[k] in core_set
                and path[k + 1] in core_set
                and mol.GetAtomWithIdx(int(path[k])).GetDegree() > 1
                and mol.GetAtomWithIdx(int(path[k + 1])).GetDegree() > 1
            )
            # A one-atom bridge (Me2Si, CMe2) has no torsional freedom between the rings;
            # two or more rotatable links do.
            if rotatable > 2:
                return False, f"flexible bridge ({rotatable} rotatable links)"
    return True, ""


def analyze(xyz_path):
    from pathlib import Path

    from oinsmiles.utils.xyz2mol import get_tmc_mol

    with _silence_fds():
        mol, _ = get_tmc_mol(Path(xyz_path), 0, with_stereo=False)
    coords = mol.GetConformer().GetPositions()
    built = build_core(mol, coords)
    if built is None:
        return {"error": "no eta ring found"}
    core_idx, metal_idx = built
    rmsd, n_autos, cap = core_chirality(mol, coords, core_idx, metal_idx)
    rigid, reason = core_is_rigid(mol, coords, core_idx, metal_idx)
    if not rigid:
        verdict = "UNDECIDABLE"
    else:
        verdict = "CHIRAL" if rmsd > CORE_ACHIRAL_RMSD else "ACHIRAL"
    return {
        "n_atoms": mol.GetNumAtoms(),
        "n_core": len(core_idx),
        "core_rmsd": round(rmsd, 3),
        "core_automorphisms": n_autos,
        "cap_hit": cap,
        "rigid": rigid,
        "note": reason,
        "verdict": verdict,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--xyz", action="append", default=[])
    ap.add_argument("--results-dir")
    ap.add_argument("--only")
    ap.add_argument("--json-out")
    args = ap.parse_args()

    jobs = [(os.path.basename(p), {"single": p}) for p in args.xyz]
    if args.results_dir:
        src = os.path.abspath(args.results_dir)
        wanted = {m.strip() for m in (args.only or "").split(",") if m.strip()}
        for rp in sorted(glob.glob(os.path.join(src, "individual_reports", "*.json"))):
            with open(rp) as f:
                rep = json.load(f)
            mol = rep.get("molecule")
            if wanted and mol not in wanted:
                continue
            paths = {}
            if rep.get("input_xyz") and os.path.exists(rep["input_xyz"]):
                paths["input"] = rep["input_xyz"]
            gen = os.path.join(src, "structures", f"{mol}_generated.xyz")
            if os.path.exists(gen):
                paths["generated"] = gen
            if paths:
                jobs.append((mol, paths))

    out = {}
    print(f"{'molecule':20s} {'which':10s} {'verdict':12s} {'rmsd':>7s} {'core':>5s} {'autos':>7s}")
    for name, paths in jobs:
        out[name] = {}
        for label, p in paths.items():
            try:
                r = analyze(p)
            except Exception as e:  # noqa: BLE001
                r = {"error": f"{type(e).__name__}: {e}"}
            out[name][label] = r
            if "error" in r:
                print(f"{name:20s} {label:10s} ERROR {r['error'][:50]}")
            else:
                print(
                    f"{name:20s} {label:10s} {r['verdict']:12s} {r['core_rmsd']:7.3f} "
                    f"{r['n_core']:5d} {r['core_automorphisms']:7d}"
                    f"{'  CAP' if r['cap_hit'] else ''}"
                    f"{'  ' + r['note'] if r['note'] else ''}"
                )
    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump(out, f, indent=2)
        print(f"Wrote {args.json_out}")


if __name__ == "__main__":
    main()
