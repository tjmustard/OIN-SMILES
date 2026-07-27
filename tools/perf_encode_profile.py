"""cProfile + deterministic-counter attribution for the ENCODE step (XYZ -> OIN).

Why this exists: every perf wave in this project targeted *generation*. Nobody had
profiled the encoder, yet a single eta molecule's bare ``XYZToSMILES().convert()``
was measured at 46-71s, and a round trip runs the encode twice. See
``docs/agentic-notes/v0.4.5/ENCODER_PERF_v0.4.5.md``.

Two modes, both on one molecule:

  ``--mode profile``   cProfile over ``XYZToSMILES().convert()``. Under host
                       contention the *absolute* seconds are noisy, but cumulative
                       time is still valid for **relative** attribution within one
                       run (which function dominates). Reported as such.

  ``--mode counters``  monkeypatched call counters only -- exact and
                       contention-robust. Counts the encoder's plausible hot spots:
                       RDKit ``SanitizeMol`` / ``MolToSmiles`` / ring perception /
                       ``GetSubstructMatches``, ``ResonanceMolSupplier``,
                       xyz2mol_local's ``AC2BO`` / ``get_UA_pairs`` /
                       ``nx.max_weight_matching``, the aligner's ``_map_to_template``
                       / ``_brute_force_symmetries``, and the eta canonicalizers.

Usage:
    PYTHONPATH=src .venv/bin/python tools/perf_encode_profile.py \
        --dataset tmCAT-tmPHOTO_xyz_dataset --molecule QIDKUL_comp_0 --mode both
"""

from __future__ import annotations

import argparse
import cProfile
import json
import os
import pstats
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


def find_xyz(dataset_dir: str, molecule: str) -> str:
    if os.path.exists(molecule):
        return molecule
    for sub in ("cat", "photo", "regression_inputs", "cohort-v0.4.5-5k", "."):
        p = os.path.join(dataset_dir, sub, f"{molecule}.xyz")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{molecule} not found under {dataset_dir}")


def atom_count(xyz_path: str) -> int:
    with open(xyz_path) as f:
        return int(f.readline().strip())


COUNTER_NAMES = [
    "SanitizeMol",
    "MolToSmiles",
    "GetSymmSSSR",
    "GetSSSR",
    "FastFindRings",
    "GetSubstructMatches",
    "ResonanceMolSupplier",
    "AssignStereochemistry",
    "AssignAtomChiralTagsFromStructure",
    "rdCIPLabeler",
    "CanonicalRankAtoms",
    "AC2BO",
    "get_UA_pairs",
    "max_weight_matching",
    "xyz2AC_obabel",
    "_map_to_template",
    "_brute_force_symmetries",
    "canonical_eta_set_representative",
    "_orientation_symmetry_graph",
    "get_oin_string",
    "get_tmc_mol",
    "BO2mol",
    "AC2mol",
    "valences_not_too_large",
    "chiral_stereo_check",
]


class Counters(dict):
    def bump(self, name):
        self[name] = self.get(name, 0) + 1


def _wrap(counters, name, obj, attr):
    orig = getattr(obj, attr, None)
    if orig is None:
        return None

    def wrapper(*a, **kw):
        counters.bump(name)
        return orig(*a, **kw)

    wrapper.__name__ = getattr(orig, "__name__", attr)
    setattr(obj, attr, wrapper)
    return (obj, attr, orig)


def _wrap_timed(counters, name, obj, attr, timings):
    """Like _wrap but also accumulates inclusive wall time under `name`."""
    orig = getattr(obj, attr, None)
    if orig is None:
        return None

    def wrapper(*a, **kw):
        counters.bump(name)
        t0 = time.perf_counter()
        try:
            return orig(*a, **kw)
        finally:
            timings[name] = timings.get(name, 0.0) + (time.perf_counter() - t0)

    wrapper.__name__ = getattr(orig, "__name__", attr)
    setattr(obj, attr, wrapper)
    return (obj, attr, orig)


def instrument(counters, timings):
    """Monkeypatch encoder hot-spot candidates to count. Returns restore list."""
    import networkx as nx
    from rdkit import Chem
    from rdkit.Chem import rdCIPLabeler, rdmolops

    from oinsmiles.utils import oin_aligner as al_mod
    from oinsmiles.utils import xyz2mol as x2m_mod
    from oinsmiles.utils import xyz2mol_local as loc_mod

    restores = []

    def add(r):
        if r:
            restores.append(r)

    # --- RDKit level (module attribute patch; call sites that did
    #     `from rdkit import Chem` then `Chem.SanitizeMol(...)` are covered) ---
    add(_wrap_timed(counters, "SanitizeMol", Chem, "SanitizeMol", timings))
    add(_wrap_timed(counters, "MolToSmiles", Chem, "MolToSmiles", timings))
    add(_wrap_timed(counters, "GetSymmSSSR", Chem, "GetSymmSSSR", timings))
    add(_wrap_timed(counters, "GetSSSR", Chem, "GetSSSR", timings))
    add(_wrap(counters, "FastFindRings", Chem, "FastFindRings"))
    add(_wrap_timed(counters, "ResonanceMolSupplier", Chem, "ResonanceMolSupplier", timings))
    add(_wrap_timed(counters, "AssignStereochemistry", Chem, "AssignStereochemistry", timings))
    add(
        _wrap_timed(
            counters,
            "AssignAtomChiralTagsFromStructure",
            Chem,
            "AssignAtomChiralTagsFromStructure",
            timings,
        )
    )
    add(_wrap_timed(counters, "rdCIPLabeler", rdCIPLabeler, "AssignCIPLabels", timings))
    add(_wrap_timed(counters, "CanonicalRankAtoms", Chem, "CanonicalRankAtoms", timings))
    add(_wrap_timed(counters, "CanonicalRankAtoms", rdmolops, "CanonicalRankAtoms", timings))

    # GetSubstructMatches is a bound C++ method on Mol -- can't patch per-instance
    # cheaply, so count via the free function form used in this codebase where
    # possible. Mol.GetSubstructMatches is patchable on the class in Python RDKit.
    try:
        add(_wrap_timed(counters, "GetSubstructMatches", Chem.Mol, "GetSubstructMatches", timings))
    except (AttributeError, TypeError):
        pass

    # --- xyz2mol_local (bond-order perception) ---
    for fn in (
        "AC2BO",
        "get_UA_pairs",
        "xyz2AC_obabel",
        "BO2mol",
        "AC2mol",
        "valences_not_too_large",
        "chiral_stereo_check",
        "xyz2AC",
        "get_BO",
        "BO_is_OK",
        "charge_is_OK",
        "get_UA",
        "get_bonds",
        "_ordered_valences",
        "get_atomic_charge",
    ):
        add(_wrap_timed(counters, fn, loc_mod, fn, timings))
    add(_wrap_timed(counters, "max_weight_matching", nx, "max_weight_matching", timings))
    # nx.max_weight_matching may be imported by name inside xyz2mol_local
    add(_wrap_timed(counters, "max_weight_matching", loc_mod, "max_weight_matching", timings))

    # --- xyz2mol (encoder driver) ---
    for fn in ("get_oin_string", "get_tmc_mol"):
        add(_wrap_timed(counters, fn, x2m_mod, fn, timings))

    # --- aligner ---
    for fn in (
        "_map_to_template",
        "_brute_force_symmetries",
        "canonical_eta_set_representative",
        "_orientation_symmetry_graph",
    ):
        add(_wrap_timed(counters, fn, al_mod, fn, timings))
        cls = getattr(al_mod, "OINDiscreteAligner", None)
        if cls is not None:
            add(_wrap_timed(counters, fn, cls, fn, timings))

    return restores


def restore(restores):
    for obj, attr, orig in reversed(restores):
        setattr(obj, attr, orig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="tmCAT-tmPHOTO_xyz_dataset")
    ap.add_argument("--molecule", required=True)
    ap.add_argument("--mode", default="both", choices=["profile", "counters", "both"])
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args()

    xyz = find_xyz(args.dataset, args.molecule)
    n = atom_count(xyz)
    load = os.getloadavg()
    print(f"molecule={args.molecule} atoms={n} path={xyz}")
    print(f"loadavg_at_start={load[0]:.2f},{load[1]:.2f},{load[2]:.2f}")

    from oinsmiles import XYZToSMILES

    result = {"molecule": args.molecule, "atoms": n, "loadavg_start": load}

    if args.mode in ("counters", "both"):
        counters = Counters()
        timings: dict[str, float] = {}
        restores = instrument(counters, timings)
        try:
            t0 = time.perf_counter()
            oin = XYZToSMILES().convert(xyz)
            dt = time.perf_counter() - t0
        finally:
            restore(restores)
        print(f"\n=== COUNTERS (exact, contention-robust) ===  wall={dt:.2f}s")
        print(f"oin_len={len(oin)}")
        for k in sorted(counters, key=lambda k: -counters[k]):
            t = timings.get(k)
            tstr = f"  incl_wall={t:8.2f}s" if t is not None else ""
            print(f"  {k:38s} {counters[k]:8d}{tstr}")
        result["counters"] = dict(counters)
        result["inclusive_wall"] = timings
        result["wall_counters_mode"] = dt
        result["oin"] = oin

    if args.mode in ("profile", "both"):
        pr = cProfile.Profile()
        t0 = time.perf_counter()
        pr.enable()
        oin = XYZToSMILES().convert(xyz)
        pr.disable()
        dt = time.perf_counter() - t0
        print(
            f"\n=== cPROFILE (relative attribution only; host may be contended) === wall={dt:.2f}s"
        )
        st = pstats.Stats(pr)
        st.sort_stats("cumulative")
        st.print_stats(args.top)
        print("--- by tottime ---")
        st.sort_stats("tottime")
        st.print_stats(args.top)
        result["wall_profile_mode"] = dt
        result["oin"] = oin

    result["loadavg_end"] = os.getloadavg()
    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump(result, f, indent=2)
        print(f"\nwrote {args.json_out}")


if __name__ == "__main__":
    main()
