"""Per-atom provenance for the ``Atom count mismatch`` class.

WHY PER-ATOM, AND NOT ANOTHER POPULATION COUNT
----------------------------------------------
Two aggregate hypotheses have already died on this class
(``docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md``). ``OIN_H_FAITHFUL`` was built for
it and bought exactly nothing: A/B over the 45-molecule population read match 8 / mismatch 37
with the lever off and **the identical 8 / 37** with it on. What that A/B did establish is
that the divergence sits between the *perceived parent* and the *emitted string* --
``perceived_H == input_H`` in 36/45 -- so perception is not the suspect. The recorded next
step was per-atom provenance. This is it.

WHAT IT MEASURES
----------------
Three things per molecule, none of which consults a bond graph the generator produced:

1. **Which atoms gained hydrogen**, by diffing a purely geometric signature -- for every heavy
   atom, ``(element, #H attached, sorted heavy-neighbour elements)``. Multiset diffs need no
   atom correspondence between the two files, which matters because the generator does not
   preserve heavy-atom order.
2. **Whether the metal lost contacts**, via ``oin.coordination.coordination_report`` -- the
   coordinate-only instrument, ~2.2 ms, that reads distances and neither bond graph.
3. **Whether the ENCODER can tell the two structures apart at all**: re-encode the generated
   XYZ independently and compare to ``smiles_1``. This is the decisive one. If a structure with
   two extra hydrogens encodes to the same OIN as the input, then no string comparison --
   scored, honest, or key-based -- can ever catch this class, and the atom-count gate is the
   only instrument that can.

⚠ The comparison in (3) must use an INDEPENDENT re-encode, never the sweep's ``smiles_2``.
``smiles_2`` comes from ``get_oin_string(gen_result.mol, coords)`` -- the generator's own bond
graph -- so "the strings match" there is circular and proves nothing about the geometry.

Usage:
    PYTHONPATH=src python tools/atom_count_provenance.py \\
        --results-dir tmCAT-tmPHOTO_xyz_dataset/results-v0.4.6-sweep
"""

import argparse
import collections
import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from oinsmiles.oin.coordination import coordination_report, parse_xyz  # noqa: E402

#: Cordero covalent radii (A) for the elements this corpus actually contains. Used ONLY to
#: decide heavy-atom adjacency for the signature; metal coordination is left to
#: ``coordination.py``, which has a calibrated slack for exactly that job.
_R = {
    "H": 0.31,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Ge": 1.20,
    "As": 1.19,
    "Se": 1.20,
    "Br": 1.20,
    "Sn": 1.39,
    "Sb": 1.39,
    "Te": 1.38,
    "I": 1.39,
    "Al": 1.21,
    "Ga": 1.22,
    "In": 1.42,
    "Tl": 1.45,
    "Pb": 1.46,
    "Bi": 1.48,
    "Li": 1.28,
    "Na": 1.66,
    "K": 2.03,
    "Mg": 1.41,
    "Ca": 1.76,
    "Be": 0.96,
    "Ti": 1.60,
    "V": 1.53,
    "Cr": 1.39,
    "Mn": 1.39,
    "Fe": 1.32,
    "Co": 1.26,
    "Ni": 1.24,
    "Cu": 1.32,
    "Zn": 1.22,
    "Y": 1.90,
    "Zr": 1.75,
    "Nb": 1.64,
    "Mo": 1.54,
    "Ru": 1.46,
    "Rh": 1.42,
    "Pd": 1.39,
    "Ag": 1.45,
    "Cd": 1.44,
    "Hf": 1.75,
    "Ta": 1.70,
    "W": 1.62,
    "Re": 1.51,
    "Os": 1.44,
    "Ir": 1.41,
    "Pt": 1.36,
    "Au": 1.36,
    "Hg": 1.32,
}
_BOND_SLACK = 1.15


def _signature(symbols, coords):
    """``{(element, nH, (neighbour elements...)): count}`` over heavy atoms, from geometry only.

    Each hydrogen is assigned to its NEAREST heavy atom rather than to everything within a
    cutoff: a bridging or agostic H would otherwise be counted twice and manufacture a
    difference that is not there.
    """
    n = len(symbols)
    heavy = [i for i in range(n) if symbols[i] != "H"]
    hydro = [i for i in range(n) if symbols[i] == "H"]

    def dist(i, j):
        return math.dist(coords[i], coords[j])

    n_h = collections.Counter()
    for h in hydro:
        if heavy:
            n_h[min(heavy, key=lambda j: dist(h, j))] += 1

    sig = collections.Counter()
    for i in heavy:
        nbrs = sorted(
            symbols[j]
            for j in heavy
            if j != i
            and dist(i, j) < (_R.get(symbols[i], 1.4) + _R.get(symbols[j], 1.4)) * _BOND_SLACK
        )
        sig[(symbols[i], n_h[i], tuple(nbrs))] += 1
    return sig


def _elements(symbols):
    return collections.Counter(symbols)


def analyse(rep, results_dir, reencode=True):
    """Provenance for one report. Returns a dict, or None when nothing can be checked."""
    mol = rep["molecule"]
    gen_path = os.path.join(results_dir, "structures", f"{mol}_generated.xyz")
    if not os.path.exists(gen_path) or not os.path.exists(rep.get("input_xyz", "")):
        return None

    in_text = open(rep["input_xyz"]).read()
    gen_text = open(gen_path).read()
    in_syms, in_xyz = parse_xyz(in_text)
    gen_syms, gen_xyz = parse_xyz(gen_text)

    delta_el = _elements(gen_syms)
    delta_el.subtract(_elements(in_syms))
    delta_el = {k: v for k, v in delta_el.items() if v}

    sig_in, sig_gen = _signature(in_syms, in_xyz), _signature(gen_syms, gen_xyz)
    gained = sig_gen - sig_in
    lost = sig_in - sig_gen

    # A transition is "protonation in place" when the SAME element with the SAME heavy
    # neighbours appears with one more hydrogen. Anything else means the heavy-atom skeleton
    # itself moved, which is a different (and worse) finding than an added H.
    protonated_in_place, skeleton_moved = [], []
    remaining = collections.Counter(gained)
    for (el, nh, nbrs), cnt in list(lost.items()):
        key = (el, nh + 1, nbrs)
        take = min(cnt, remaining.get(key, 0))
        if take:
            protonated_in_place.append(
                {"element": el, "from_H": nh, "neighbours": list(nbrs), "n": take}
            )
            remaining[key] -= take
    for (el, nh, nbrs), cnt in remaining.items():
        if cnt > 0:
            skeleton_moved.append({"element": el, "H": nh, "neighbours": list(nbrs), "n": cnt})

    coord = coordination_report(in_text, gen_text)

    indep, indep_matches = None, None
    if reencode:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
        from oinsmiles import XYZToSMILES

        try:
            indep = XYZToSMILES().convert(gen_path)
            indep_matches = indep == rep.get("smiles_1")
        except Exception as e:  # noqa: BLE001 - an encode failure is a datum
            indep = f"ENCODE_FAILED: {type(e).__name__}"
            indep_matches = None

    return {
        "molecule": mol,
        "n_in": len(in_syms),
        "n_gen": len(gen_syms),
        "delta_elements": delta_el,
        "protonated_in_place": protonated_in_place,
        "skeleton_moved": skeleton_moved,
        "coordination_intact": coord.get("intact"),
        "coordination_reason": coord.get("reason"),
        "smiles_1": rep.get("smiles_1"),
        "smiles_2_indep": indep,
        "indep_equals_input_oin": indep_matches,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--only", default=None, help="comma-separated molecule ids")
    ap.add_argument("--json-out", default=None)
    ap.add_argument("--no-reencode", action="store_true", help="skip the independent re-encode")
    args = ap.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    only = set(args.only.split(",")) if args.only else None

    rows = []
    for f in sorted(glob.glob(os.path.join(results_dir, "individual_reports", "*.json"))):
        rep = json.load(open(f))
        if only is not None:
            if rep.get("molecule") in only:
                rows.append(rep)
        elif "Atom count mismatch" in (rep.get("error") or ""):
            rows.append(rep)

    out = []
    for rep in sorted(rows, key=lambda r: r["molecule"]):
        res = analyse(rep, results_dir, reencode=not args.no_reencode)
        if res is None:
            print(f"{rep['molecule']:22s}  SKIPPED (no stored structure)")
            continue
        out.append(res)
        print(
            f"\n=== {res['molecule']}   {res['n_in']} -> {res['n_gen']} atoms   {res['delta_elements']}"
        )
        for p in res["protonated_in_place"]:
            print(
                f"    +H in place : {p['n']}x  {p['element']}({p['from_H']}H -> {p['from_H'] + 1}H)  nbrs={p['neighbours']}"
            )
        for s in res["skeleton_moved"]:
            print(f"    skeleton    : {s['n']}x  {s['element']} {s['H']}H  nbrs={s['neighbours']}")
        print(
            f"    coordination: intact={res['coordination_intact']}  {res['coordination_reason'] or ''}"
        )
        if not args.no_reencode:
            print(f"    independent re-encode == input OIN?  {res['indep_equals_input_oin']}")

    if out:
        n = len(out)
        n_flag = sum(1 for r in out if r["coordination_intact"] is False)
        n_same = sum(1 for r in out if r["indep_equals_input_oin"] is True)
        n_prot = sum(1 for r in out if r["protonated_in_place"])
        h_only = sum(1 for r in out if set(r["delta_elements"]) == {"H"})
        print(f"\n{'=' * 72}")
        print(f"population                                        : {n}")
        print(f"element delta is HYDROGEN ONLY                    : {h_only}/{n}")
        print(f"at least one atom protonated in place             : {n_prot}/{n}")
        print(f"coordination flags a lost metal contact           : {n_flag}/{n}")
        if not args.no_reencode:
            print(f"INDEPENDENT re-encode still equals the input OIN  : {n_same}/{n}")
            print(
                "\nThat last line is the verdict. Where it is TRUE, the encoder maps a structure\n"
                "with extra hydrogens onto the SAME OIN as the input -- so no string comparison,\n"
                "scored or honest, can detect it, and the atom-count gate is the only instrument\n"
                "that can."
            )

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nwrote {args.json_out}")


if __name__ == "__main__":
    main()
