#!/usr/bin/env python
"""§6.5 FALSIFICATION: does a coordinate-only donor-set predicate separate arm A's accepted
conformer from arm B's on the known `OIN_ACCEPT_SCORED` regressions?

THE TRAP THIS AVOIDS (docs/agentic-notes/v0.4.7/ACCEPT_SCORED_v0.4.7.md §6.1)
--------------------------------------------------------
`metallogen_adapter._coordination_vectors` derives donors from
`metal.GetBonds()` -- the GENERATOR'S OWN GRAPH. A ligand that has physically left the
coordination sphere KEEPS ITS BOND OBJECT, so any attachment check built on that path is
blind by construction and would certify exactly the structures it exists to catch.

Here `GetBonds()` supplies only the REFERENCE (which atoms the generator CLAIMS are donors).
The MEASUREMENT is entirely coordinate-derived:

    AC, _ = xyz2AC_obabel(atomic_numbers, coords, tolerance=0.5)   # "modified to capture
    actual = nonzero(AC[metal_idx])                                #  haptic bonds"

That is definitionally the computation `XYZToSMILES().convert()` uses to decide coordination
(`get_basic_mol` -> `_get_tmc_mol_impl`: `np.nonzero(GetAdjacencyMatrix(mol)[tmc_idx, :])`),
so it is the computation whose change produced the six haptic `indep` failures. A detached
ligand stays in `claim` and vanishes from `actual`; the mismatch is the signal.

CANDIDATE PREDICATES, all scored side by side so the choice is made on data
--------------------------------------------------------------------------
    count   |actual| == |claimed|                 -- §6.4 predicts 6/8
    subset  claimed subseteq actual               -- set-based, §6.4 predicts 7/8
    setequal claimed == actual                    -- strictest; expected to over-reject,
                                                     since `_get_tmc_mol_impl` FILTERS
                                                     aromatic ring carbons adjacent to a
                                                     coordinating heteroatom out of the donor
                                                     set and this raw AC row does not.

Usage:
    python tools/attach_probe.py --dump spec/handoffs/v0.4.7/dump21 --json out.json
    python tools/attach_probe.py --inputs spec/handoffs/v0.4.7/cohort_attach21.json
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import GetPeriodicTable  # noqa: E402


def read_xyz(path: str):
    """(atomic_numbers, coords) from an XYZ file. Order is preserved and load-bearing:
    it must match the generator mol's atom order, which is how `res.xyz` is written."""
    pt = GetPeriodicTable()
    lines = open(path).read().splitlines()
    n = int(lines[0].strip())
    z, xyz = [], []
    for i in range(n):
        f = lines[2 + i].split()
        z.append(pt.GetAtomicNumber(f[0]))
        xyz.append([float(v) for v in f[1:4]])
    return z, np.array(xyz)


def actual_donor_set(znums, coords, tolerance: float = 0.5):
    """The metal's donor set AS THE INDEPENDENT ENCODER COMPUTES IT -- from coordinates only.

    Returns (metal_idx, set_of_donor_indices, distances). No bond object is consulted.
    """
    from oinsmiles.utils.perception_core import xyz2AC_obabel
    from oinsmiles.utils.perception_tmc import TRANSITION_METALS_NUM

    metal_idx = next((i for i, z in enumerate(znums) if z in TRANSITION_METALS_NUM), None)
    if metal_idx is None:
        return None, set(), {}
    AC, _ = xyz2AC_obabel(list(znums), coords, tolerance=tolerance)
    donors = set(int(j) for j in np.nonzero(AC[metal_idx])[0])
    d = {
        int(j): round(float(np.linalg.norm(coords[j] - coords[metal_idx])), 3)
        for j in sorted(donors)
    }
    return metal_idx, donors, d


def fast_donor_set(znums, coords, tolerance: float = 0.5):
    """The metal's donor row by the DISTANCE CRITERION ALONE -- `xyz2AC_obabel`'s first pass,
    without building the full N x N matrix and without the valence-cap pruning loop.

    `xyz2AC_obabel` prunes an atom's longest bonds when its degree exceeds max(atomic_valence),
    and that loop CAN touch the metal row (the code's own DUDREA_comp_0 example drops a bridging
    Y-H). So this is an approximation of `actual_donor_set` and the two are measured against
    each other on every conformer rather than assumed equal -- `ac_row_divergence` below.
    """
    from oinsmiles.utils.perception_tmc import TRANSITION_METALS_NUM

    pt = GetPeriodicTable()
    metal_idx = next((i for i, z in enumerate(znums) if z in TRANSITION_METALS_NUM), None)
    if metal_idx is None:
        return None, set()
    rad = np.array([pt.GetRcovalent(int(z)) for z in znums])
    d = np.linalg.norm(coords - coords[metal_idx], axis=1)
    hit = (d <= rad[metal_idx] + rad + tolerance) & (np.arange(len(znums)) != metal_idx)
    return metal_idx, set(int(j) for j in np.nonzero(hit)[0])


_HETEROATOM_DONORS = {7, 8, 15, 16}  # N, O, P, S -- _get_tmc_mol_impl's own set


def encoder_donor_set(znums, coords, tolerance: float = 0.5):
    """`_get_tmc_mol_impl`'s donor set EXACTLY: the AC metal row, then its aromatic-ring-carbon
    filter (a ring carbon whose neighbour is a coordinating N/O/P/S is dropped -- the
    heteroatom is the real donor, e.g. the phenylene bridge of a bidentate phosphine).

    `actual_donor_set` omits that filter and therefore over-counts donors on chelates. Whether
    that over-count matters is measured, not assumed: both sets are scored side by side.
    """
    from oinsmiles.utils.perception_core import xyz2AC_obabel
    from oinsmiles.utils.perception_tmc import TRANSITION_METALS_NUM

    metal_idx = next((i for i, z in enumerate(znums) if z in TRANSITION_METALS_NUM), None)
    if metal_idx is None:
        return None, set()
    AC, proto = xyz2AC_obabel(list(znums), coords, tolerance=tolerance)
    raw = set(int(j) for j in np.nonzero(AC[metal_idx])[0])
    rw = Chem.RWMol(proto)
    n = len(AC)
    for i in range(n):
        for j in range(i + 1, n):
            if AC[i, j]:
                rw.AddBond(i, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    try:
        Chem.GetSymmSSSR(mol)
    except Exception:
        return metal_idx, raw
    keep = set()
    for idx in raw:
        a = mol.GetAtomWithIdx(idx)
        if a.GetAtomicNum() == 6 and a.IsInRing():
            if any(
                nb.GetAtomicNum() in _HETEROATOM_DONORS and nb.GetIdx() in raw
                for nb in a.GetNeighbors()
            ):
                continue
        keep.add(idx)
    return metal_idx, keep


_HAPTIC_GROUP_CUTOFF = 1.6  # A -- metallogen_adapter._HAPTIC_GROUP_CUTOFF, kept in step


def _groups(indices, coords):
    """Transitive single-linkage grouping at 1.6 A -- the encoder's hapticity reduction
    (`_reduce_haptic_positions` / `oin_aligner._reduce_hapticity`). A whole Cp ring is one
    group, i.e. ONE coordination site, which is the unit the OIN's slot numbering uses."""
    idx = list(indices)
    seen, out = set(), []
    for a in range(len(idx)):
        if a in seen:
            continue
        stack, comp = [a], []
        while stack:
            c = stack.pop()
            if c in seen:
                continue
            seen.add(c)
            comp.append(idx[c])
            for k in range(len(idx)):
                if (
                    k not in seen
                    and float(np.linalg.norm(coords[idx[c]] - coords[idx[k]]))
                    < _HAPTIC_GROUP_CUTOFF
                ):
                    stack.append(k)
        out.append(sorted(comp))
    return out


def _oin_slots(oin: str | None) -> int:
    """How many distinct coordination SLOTS the requested OIN states. An eta ring writes the
    same slot number on every ring atom, so this counts coordination SITES, not donor atoms."""
    import re

    if not oin:
        return 0
    return len(set(int(x) for x in re.findall(r"\{(\d+)[<>]?\}", oin)))


def score_one(xyz_path: str, claim_path: str, tolerance: float = 0.5) -> dict:
    claim = json.load(open(claim_path))
    znums, coords = read_xyz(xyz_path)
    t0 = time.monotonic()
    midx, actual, dists = actual_donor_set(znums, coords, tolerance)
    ms = round((time.monotonic() - t0) * 1000, 1)
    t1 = time.monotonic()
    _, fast_actual = fast_donor_set(znums, coords, tolerance)
    fast_ms = round((time.monotonic() - t1) * 1000, 3)
    t2 = time.monotonic()
    _, filt = encoder_donor_set(znums, coords, tolerance)
    filt_ms = round((time.monotonic() - t2) * 1000, 1)

    claimed = set(claim.get("claimed_donors") or [])
    pt = GetPeriodicTable()
    sym = {i: pt.GetElementSymbol(znums[i]) for i in range(len(znums))}

    # ⚠ ALIGNMENT GUARD, and it is load-bearing. `res.xyz` is written from the RAW MetalloGen
    # mol while `res.mol` is the CONTRACT mol, so the two share an atom ordering only by
    # convention (the harness's own `get_oin_string(gen.mol, coords)` relies on it). If that
    # convention ever broke, `claimed` would index into the wrong atoms and this probe would
    # manufacture a separation out of a numbering bug. Compare the element symbols the two
    # sources give for the same indices; a mismatch invalidates the row rather than scoring it.
    aligned = (
        claim.get("metal_idx") == midx
        and [sym.get(i) for i in sorted(claimed)] == list(claim.get("claimed_elements") or [])
        and claim.get("natoms") == len(znums)
    )
    missing = sorted(claimed - actual)  # generator says bonded, geometry says NOT
    extra = sorted(actual - claimed)  # geometry says bonded, generator does not claim it

    n_slots = _oin_slots(claim.get("oin_in"))
    actual_groups = _groups(sorted(actual), coords)
    claim_groups = _groups(sorted(claimed), coords)
    # A claimed coordination SITE is lost when NOT ONE of its atoms is still within bonding
    # distance of the metal. Ring slip (5 of 5 -> 3 of 5) leaves the site present; a ring that
    # has drifted out of the sphere leaves it empty.
    lost_groups = [[sym.get(i) for i in g] for g in claim_groups if not (set(g) & actual)]

    filt_groups = _groups(sorted(filt), coords)
    lost_groups_f = [[sym.get(i) for i in g] for g in claim_groups if not (set(g) & filt)]
    return {
        "molecule": claim["molecule"],
        "natoms": len(znums),
        "metal": sym.get(midx),
        "metal_idx_xyz": midx,
        "metal_idx_claim": claim.get("metal_idx"),
        "n_claimed": len(claimed),
        "n_actual": len(actual),
        "claimed": sorted(claimed),
        "actual": sorted(actual),
        "missing": missing,
        "missing_elements": [sym.get(i) for i in missing],
        "extra": extra,
        "extra_elements": [sym.get(i) for i in extra],
        "missing_dists": {
            str(i): dists.get(i, round(float(np.linalg.norm(coords[i] - coords[midx])), 3))
            for i in missing
        },
        "donor_dist_range": ([min(dists.values()), max(dists.values())] if dists else None),
        "aligned": aligned,
        "n_slots_oin": n_slots,
        "n_sites_actual": len(actual_groups),
        "n_sites_claimed": len(claim_groups),
        "sites_lost": lost_groups,
        "p_count": len(actual) == len(claimed),
        "p_subset": claimed <= actual,
        "p_setequal": claimed == actual,
        # P3 -- the coordinate-derived donors group into exactly the number of coordination
        # SITES the requested OIN states. This is the predicate §6.2 actually argues for: the
        # six haptic failures are all "metal geometry tag degrades by exactly one donor", i.e.
        # a change in this number.
        "p_sites": len(actual_groups) == n_slots and n_slots > 0,
        # P4 -- every coordination site the generator CLAIMS still retains at least one atom
        # inside bonding distance. Tolerant of eta5 -> eta3 ring slip (which the round-trip key
        # forgives) while still catching a ring that has left the metal entirely.
        "p_sitecov": not lost_groups,
        # ...and the same two predicates computed on the FILTERED donor set, i.e. exactly the
        # set `_get_tmc_mol_impl` hands downstream.
        "n_actual_filt": len(filt),
        "n_sites_filt": len(filt_groups),
        "filt_ms": filt_ms,
        "p_sites_f": len(filt_groups) == n_slots and n_slots > 0,
        "p_sitecov_f": not lost_groups_f,
        "check_ms": ms,
        "fast_ms": fast_ms,
        "n_actual_fast": len(fast_actual),
        # Does skipping the valence-cap pruning loop change the METAL's row? Measured, not
        # assumed -- the cheap form is only usable if this is empty everywhere.
        "ac_row_divergence": sorted(fast_actual ^ actual),
        "indep_passed": claim.get("indep_passed"),
        "passed": claim.get("passed"),
    }


def score_dump(dump_root: str, tolerance: float) -> dict:
    arms = sorted(d for d in os.listdir(dump_root) if os.path.isdir(os.path.join(dump_root, d)))
    out: dict = {"tolerance": tolerance, "arms": {}}
    for arm in arms:
        rows = []
        for cp in sorted(glob.glob(os.path.join(dump_root, arm, "*.claim.json"))):
            xp = cp[: -len(".claim.json")] + ".xyz"
            if not os.path.exists(xp):
                continue
            rows.append(score_one(xp, cp, tolerance))
        out["arms"][arm] = rows
    return out


def report(res: dict) -> None:
    arms = list(res["arms"])
    if len(arms) != 2:
        for arm in arms:
            print(f"\n=== {arm} ===")
            for r in res["arms"][arm]:
                print(
                    f"  {r['molecule']:16s} {r['metal']:>3s} claimed={r['n_claimed']:2d} "
                    f"actual={r['n_actual']:2d} sites={r['n_sites_actual']}/{r['n_slots_oin']} "
                    f"lost={r['sites_lost']} indep={r['indep_passed']}"
                )
        return

    a, b = arms
    A = {r["molecule"]: r for r in res["arms"][a]}
    B = {r["molecule"]: r for r in res["arms"][b]}
    mols = sorted(set(A) & set(B))

    bad_align = [
        f"{m}[{arm}]" for m in mols for arm, r in ((a, A[m]), (b, B[m])) if not r.get("aligned")
    ]
    print(
        "\n  ALIGNMENT GUARD (res.xyz order == res.mol order): "
        + (
            f"OK on all {2 * len(mols)}"
            if not bad_align
            else f"*** VIOLATED: {bad_align} -- rows below are NOT interpretable ***"
        )
    )

    print(f"\n=== §6.5 FALSIFICATION  ({a} vs {b}, tolerance {res['tolerance']}) ===")
    print(
        f"{'molecule':16s} {'M':>3s} {'slot':>4s} {'clm':>4s} "
        f"{'Aact':>5s}{'Bact':>5s} {'Astf':>6s}{'Bstf':>6s}  "
        f"{'A_ind':>6s}{'B_ind':>6s}  B_sites_lost"
    )
    for m in mols:
        ra, rb = A[m], B[m]
        print(
            f"{m:16s} {ra['metal']:>3s} {ra['n_slots_oin']:>4d} {ra['n_claimed']:>4d} "
            f"{ra['n_actual']:>5d}{rb['n_actual']:>5d} "
            f"{ra['n_sites_filt']:>6d}{rb['n_sites_filt']:>6d}  "
            f"{str(ra['indep_passed']):>6s}{str(rb['indep_passed']):>6s}  "
            f"{rb['sites_lost'] if rb['sites_lost'] else ''}"
        )

    CANDIDATES = (
        ("P1 count", "p_count"),
        ("P2 subset", "p_subset"),
        ("P3 sites==slots", "p_sites"),
        ("P4 site coverage", "p_sitecov"),
        ("P3f sites==slots FILT", "p_sites_f"),
        ("P4f site cov FILT", "p_sitecov_f"),
        ("P3f and P4f", "__and_f"),
    )

    def val(r, key):
        if key == "__and_f":
            return r["p_sites_f"] and r["p_sitecov_f"]
        return r[key]

    reg = [m for m in mols if A[m]["indep_passed"] and not B[m]["indep_passed"]]
    print(f"\n---- separation on the {len(reg)} INDEP REGRESSIONS ----")
    for name, key in CANDIDATES:
        sep = [m for m in reg if val(A[m], key) and not val(B[m], key)]
        wrong_a = [m for m in reg if not val(A[m], key)]
        miss_b = [m for m in reg if val(B[m], key)]
        print(
            f"  {name:18s} separates {len(sep)}/{len(reg)}  "
            f"| rejects arm A too (FALSE POSITIVE): {wrong_a or 'none'} "
            f"| accepts arm B (MISSED): {miss_b or 'none'}"
        )

    print("\n---- FALSE-POSITIVE control: every indep-PASSING conformer must be ACCEPTED ----")
    n_pass = sum(1 for m in mols for r in (A[m], B[m]) if r["indep_passed"])
    for name, key in CANDIDATES:
        bad = [
            f"{m}[{arm}]"
            for m in mols
            for arm, r in ((a, A[m]), (b, B[m]))
            if r["indep_passed"] and not val(r, key)
        ]
        print(f"  {name:18s} wrongly rejects {len(bad)}/{n_pass}: {bad or 'none'}")

    print("\n---- BONUS: does it also flag conformers that ALREADY fail indep in arm A? ----")
    fails = [
        (m, arm) for m in mols for arm, r in ((a, A[m]), (b, B[m])) if r["indep_passed"] is False
    ]
    for name, key in CANDIDATES:
        caught = [f"{m}[{arm}]" for m, arm in fails if not val(A[m] if arm == a else B[m], key)]
        print(f"  {name:18s} flags {len(caught)}/{len(fails)} indep-FAILING conformers")

    rows = [r for arm in arms for r in res["arms"][arm]]
    ms = [r["check_ms"] for r in rows]
    fms = [r["fast_ms"] for r in rows if "fast_ms" in r]
    div = [r["molecule"] for r in rows if r.get("ac_row_divergence")]
    if ms:
        print(
            f"\n  full  xyz2AC_obabel check: {min(ms):.1f}-{max(ms):.1f} ms "
            f"(median {sorted(ms)[len(ms) // 2]:.1f} ms) over {len(ms)} evaluations"
        )
    if fms:
        print(
            f"  fast  metal-row-only check: {min(fms):.3f}-{max(fms):.3f} ms "
            f"(median {sorted(fms)[len(fms) // 2]:.3f} ms)"
        )
        print(
            f"  metal row identical without the valence-cap pruning loop: "
            f"{len(rows) - len(div)}/{len(rows)}" + (f"  DIVERGES on {div}" if div else "")
        )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", help="directory written by ab_accept_scored.py --dump-xyz")
    ap.add_argument(
        "--inputs",
        help="cohort JSON: run the coordinate-only donor detection on the REAL CRYSTAL "
        "inputs, the 'does it reject real structures?' control",
    )
    ap.add_argument("--tolerance", type=float, default=0.5)
    ap.add_argument("--json")
    args = ap.parse_args()

    if args.inputs:
        raw = json.load(open(args.inputs))
        recs = list(raw.get("eta", [])) + list(raw.get("control", []))
        pt = GetPeriodicTable()
        rows = []
        print(
            f"{'molecule':16s} {'M':>3s} {'atoms':>6s} {'donors':>7s} {'d range (A)':>14s} {'ms':>7s}"
        )
        for rec in recs:
            if not os.path.exists(rec["xyz"]):
                continue
            z, c = read_xyz(rec["xyz"])
            t0 = time.monotonic()
            midx, donors, d = actual_donor_set(z, c, args.tolerance)
            ms = (time.monotonic() - t0) * 1000
            dr = f"{min(d.values()):.2f}-{max(d.values()):.2f}" if d else "-"
            rows.append(
                {
                    "molecule": rec["mol"],
                    "natoms": len(z),
                    "metal": pt.GetElementSymbol(z[midx]) if midx is not None else None,
                    "n_donors": len(donors),
                    "d_min": min(d.values()) if d else None,
                    "d_max": max(d.values()) if d else None,
                    "ms": round(ms, 1),
                }
            )
            print(
                f"{rec['mol']:16s} {rows[-1]['metal'] or '?':>3s} {len(z):>6d} "
                f"{len(donors):>7d} {dr:>14s} {ms:>7.1f}"
            )
        if args.json:
            json.dump({"inputs": rows}, open(args.json, "w"), indent=1)
        return 0

    if not args.dump:
        ap.error("need --dump or --inputs")
    res = score_dump(args.dump, args.tolerance)
    report(res)
    if args.json:
        json.dump(res, open(args.json, "w"), indent=1)
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
