"""Is the donor fold's swap a SYMMETRY of the ligand, or only a swap of two look-alike atoms?

WHY THIS EXISTS
---------------
The census (C2) blamed ``fold_parity.resolve``'s left conjunct and prescribed a geometric
achirality test. Reading ``canonical_slots._donor_swap_permutations`` suggests the defect sits one
level down. The fold exchanges donors **bucket by bucket**: every symmetry class of a fragment is
permuted independently of every other class. A ligand automorphism does not work that way -- the
C2 axis of a linear tetradentate ``t1-i1-i2-t2`` swaps the terminal pair AND the inner pair in one
move. Swapping only the terminals is not an automorphism of the ligand; the labeling it produces
describes a DIFFERENT arrangement (the mirror image for cis-alpha, an impossible isomer for a
planar macrocycle). That is why the fold needed a veto at all.

THE CLAIM UNDER TEST (it is a theorem if the premises hold, so test the premises)
-----------------------------------------------------------------------------
Let ``G`` = the polyhedron's proper rotations, ``A`` = slot permutations induced by TRUE
automorphisms of each fragment (all classes moved together, chirality respected). Every candidate
in ``G x A`` describes the same molecule as the input, so:

* an **exact fold** over ``G x A`` cannot collapse an enantiomer pair -- no veto needed;
* if x is achiral (some improper ``q`` with ``q.L = L.a``, ``a`` in ``A``), the mirror's labeling
  lies in x's own candidate set -- ``S_auto(x) == S_auto(mirror x)`` with no achirality test;
* the ruler's sphere test matches donors by CLASS and ligand id, not by automorphism, so it shares
  the bucket fold's blind spot: a cis-alpha tetradentate with planar donors reads "achiral".

Everything here is OFFLINE and EXACT: the slot post-pass is a pure function of the inline string,
and its candidate set is an orbit, so applying a fold to the stored rotation-only string
(``fold_off``) gives byte-for-byte what the encoder would emit. The positive control proves that
before any other number is read: the SHIPPED bucket fold applied to ``fold_off`` must reproduce
the stored ``veto_off`` string on every record.

    PYTHONPATH=$PWD/src <main>/.venv/bin/python tools/v0417/autofold_audit.py
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter
from pathlib import Path

DEF_OUT = Path(
    "/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/results-census"
)
EXACT = "OIN_EXACT_DONOR_FOLD"


def _load(path):
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not rows or not rows[-1].startswith("#DONE"):
        raise SystemExit(f"{path}: no #DONE trailer -- refusing to read a partial file")
    return [json.loads(ln) for ln in rows[:-1]]


def _canon(s, exact):
    """The encoder's own slot post-pass, with the exact-fold lever forced for this one call.

    Deliberately NOT a re-implementation: the first draft of this tool carried its own copy of
    the automorphism search, which is a second thing that can drift. The lever is written as
    "1"/"0", never unset -- ``lever_enabled`` is the only reader.
    """
    from oinsmiles.oin.canonical_slots import canonicalize_oin_slots

    prior = os.environ.get(EXACT)
    os.environ[EXACT] = "1" if exact else "0"
    try:
        return canonicalize_oin_slots(s)
    finally:
        if prior is None:
            os.environ.pop(EXACT, None)  # restore: `prior` was unset
        else:
            os.environ[EXACT] = prior


def s_bucket(s):
    return _canon(s, exact=False)


def s_auto(s):
    return _canon(s, exact=True)


def audit(recs, title):
    print(f"\n== {title}   n={len(recs)}")
    usable = [r for r in recs if r["fold_off"]["x"] and r["fold_off"]["m"] and r["veto_off"]["x"]]
    print(f"   usable (x and mirror both encoded, all arms): {len(usable)}/{len(recs)}")

    # POSITIVE CONTROL: the offline path must reproduce the encoder, byte for byte.
    ctl = Counter()
    for r in usable:
        for hand in ("x", "m"):
            ctl[s_bucket(r["fold_off"][hand]) == r["veto_off"][hand]] += 1
    print(f"   CONTROL  bucket_fold(fold_off) == stored veto_off : {ctl[True]}/{sum(ctl.values())}")

    tab = Counter()
    rows = []
    for r in usable:
        x, m = r["fold_off"]["x"], r["fold_off"]["m"]
        ax, am = s_auto(x), s_auto(m)
        bx, bm = r["veto_off"]["x"], r["veto_off"]["m"]
        dx, dm = r["default"]["x"], r["default"]["m"]
        vet = "vetoed_collapse" in (r["default"]["outcome_x"], r["default"]["outcome_m"])
        # Does the labeling the fold EMITS describe the input? It does iff it lies in the input's
        # own G x A orbit, i.e. iff the exact fold sends both to one string.
        row = {
            "molecule": r["molecule"],
            "rot_same": x == m,
            "bucket_same": bx == bm,
            "auto_same": ax == am,
            "default_same": dx == dm,
            "vetoed": vet,
            "bucket_x_in_orbit": s_auto(bx) == ax,
            "bucket_m_in_orbit": s_auto(bm) == am,
            "default_x_in_orbit": s_auto(dx) == ax,
            "auto_x": ax,
            "auto_m": am,
        }
        rows.append(row)
        tab[(row["bucket_same"], row["auto_same"])] += 1
    print("   E(mirror)==E(x)  under the BUCKET fold x under the EXACT fold:")
    for b in (True, False):
        for a in (True, False):
            print(f"      bucket_same={b!s:5}  auto_same={a!s:5}  {tab[(b, a)]:5d}")
    n_auto = sum(1 for w in rows if w["auto_same"])
    print(f"   EXACT fold unifies the pair           : {n_auto}/{len(rows)}")
    print(
        "   ... of the pairs a hand was VETOED on : "
        f"{sum(1 for w in rows if w['vetoed'] and w['auto_same'])}"
        f"/{sum(1 for w in rows if w['vetoed'])}"
    )
    off = [w for w in rows if not (w["bucket_x_in_orbit"] and w["bucket_m_in_orbit"])]
    print(
        f"   bucket fold emits a labeling OUTSIDE the input's orbit (x or mirror): {len(off)}/{len(rows)}"
    )
    offd = [w for w in rows if not w["default_x_in_orbit"]]
    print(f"   SHIPPED E(x) is outside x's own orbit : {len(offd)}/{len(rows)}")
    return rows


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--results", type=Path, default=DEF_OUT)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    import oinsmiles

    print("oinsmiles from:", oinsmiles.__file__)
    out = {}
    probe = _load(args.results / "veto_probe.jsonl")
    for g in sorted({r["group"] for r in probe}):
        out[g] = audit([r for r in probe if r["group"] == g], g)
    chiral = _load(args.results / "veto_probe_chiral.jsonl")
    out["chiral & differs (ALL 1,623)"] = audit(chiral, "chiral & differs (ALL 1,623)")
    prot = [
        r
        for r in chiral
        if "vetoed_collapse" in (r["default"]["outcome_x"], r["default"]["outcome_m"])
    ]
    out["the 36 PROTECTED pairs"] = audit(prot, "the PROTECTED pairs (ruler-chiral, a hand vetoed)")
    dest = args.out or (args.results.parent / "results-v0.4.17-exactfold" / "autofold_audit.json")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))
    print("wrote", dest)
    sys.stdout.flush()


if __name__ == "__main__":
    main()
