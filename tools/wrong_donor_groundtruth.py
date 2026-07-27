"""Is a slot-drift pair a SOUNDNESS defect or only a determinism one? (v0.4.5 Lane 9)

``tools/slot_drift_mechanism.py`` splits the residual ``slot_renumber`` pairs at atom level
into ``automorphism`` (benign) and ``DISTINCT_donors`` ("one of the two strings is WRONG").
That second verdict is a claim about *which physical atom is bound where*, and **it cannot be
settled from the strings alone** -- which is what that tool has. This tool settles it from the
3D coordinates.

THE ARGUMENT
============
The slot integer a donor receives comes from a Kabsch fit of the donor direction vectors onto
an idealized template (``OINDiscreteAligner._map_to_template``). That fit is **degenerate by
construction over the polyhedron's proper-rotation group**: if ``g`` is a rotation permuting
the template vertices, then ``template[p o g] == R_g . template[p]``, so aligning the donors
to one is aligning them to the other and the optimal residual is *identical*. The fit
therefore determines the donor -> vertex map only **up to that group**; which member of the
tied coset is emitted is decided by ``_map_to_template``'s strict ``<``, which keeps the
lex-first candidate -- lex over ``itertools.permutations`` of the donor enumeration order,
i.e. over **atom order**.

Consequence: two labelings related by a proper rotation are *equally faithful to the
geometry*. Neither is "wrong about which atom is bound where"; a proper rotation is a change
of reference frame, and it preserves every isomer-level relation (cis/trans, fac/mer) and
every reflection-odd descriptor (metal Delta/Lambda, eta winding, axial sign). Only a
labeling NOT so related can be a soundness defect.

WHAT THIS MEASURES, per molecule
================================
1. Encode the original XYZ, then renumbered variants (renumbering leaves the molecule
   unchanged, so both strings must describe the same physical structure).
2. Recover the fit the encoder actually SELECTED and, for each donor, its ORIGINAL xyz atom
   indices -- read off ``__origIdx``, the only place an atom's un-permuted identity survives.
3. Re-score EVERY permutation with the same ``scipy`` Kabsch the encoder uses (bypassing
   ``_candidate_permutations``, so the answer cannot be a prefilter artifact) and report how
   many tie the winner and how far the nearest non-tied one is.
4. Compute ``pi`` = (base vertex -> variant vertex) for each PHYSICAL donor and test whether
   it extends to a proper rotation of the polyhedron.

  ``pi`` in the group          -> determinism defect. Both strings faithful.
  ``pi`` == identity           -> the geometric fit did not move at all; the drift is
                                  introduced downstream of it.
  ``pi`` NOT in the group      -> a genuine soundness defect: investigate.

The verdict is invariant under the two *further* group relabelings the encoder applies after
the fit (``_permute_and_serialize``'s lex-max, and Lane 2's ``canonicalize_oin_slots``
post-pass), because conjugating a group element by group elements stays in the group. Pass
``GT_RAW=1`` to report the raw fit vertices, which is what "which donor sits at which
template vertex" actually means; without it the Lane 2 post-pass map is composed in.

USAGE
=====
    PYTHONPATH=src OIN_CANONICAL_BODY=1 OIN_CANONICAL_PERCEPTION=1 \
    OIN_CANONICAL_SLOTS=1 OIN_STABLE_METAL_AC=1 GT_RAW=1 \
        .venv/bin/python tools/wrong_donor_groundtruth.py [all|<molecule>] [n_renumberings]

Nothing here changes encoder behaviour: every patch wraps the original and returns its result.
"""

import itertools
import json
import os
import random
import re
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from rdkit import RDLogger  # noqa: E402
from scipy.spatial.transform import Rotation  # noqa: E402

RDLogger.DisableLog("rdApp.*")

from canonicality_probe import read_xyz, write_xyz  # noqa: E402

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.oin import canonical_slots as CS  # noqa: E402
from oinsmiles.utils import oin_aligner as OA  # noqa: E402
from oinsmiles.utils import perception_tmc as X2M  # noqa: E402

#: The 7 molecules Lane 2's classifier flagged ``DISTINCT_donors`` in its stacked arm.
DEFAULT_MOLS = [
    ("IMACOO_comp_0", "cat"),
    ("RUBTIS_comp_0", "cat"),
    ("VEXHIR_comp_0", "cat"),
    ("ZACFER_comp_0", "cat"),
    ("HAVGIW_comp_0", "photo"),
    ("KISQAG_comp_0", "photo"),
    ("ZOSNUS_comp_0", "photo"),
]

DEFAULT_DATASET = os.environ.get("OIN_DATASET", "tmCAT-tmPHOTO_xyz_dataset")

CAP: dict = {}
_o_pai = X2M._align_to_pai
_o_red = OA.OINDiscreteAligner._reduce_hapticity
_o_map = OA.OINDiscreteAligner._map_to_template
_o_can = CS.canonicalize_oin_slots


def _pai(tmc_mol, xyz_coords, metal_idx):
    """Capture tmc_mol index -> ORIGINAL xyz index before the principal-axis rotation.

    ``_align_to_pai`` is the last point that sees ``tmc_mol``, and ``__origIdx`` is the only
    surviving link back to the input file: ``get_tmc_mol`` rebuilds the molecule metal
    fragment first, so a mol index is NOT an xyz line number.
    """
    CAP["origidx"] = {
        a.GetIdx(): (a.GetIntProp("__origIdx") if a.HasProp("__origIdx") else a.GetIdx())
        for a in tmc_mol.GetAtoms()
    }
    return _o_pai(tmc_mol, xyz_coords, metal_idx)


def _red(self):
    """Attach each haptic-reduced virtual donor's global mol indices."""
    vas = _o_red(self)
    for i, lig in enumerate(self.ligands):
        if i == self.metal_idx or not lig.get("binding_atoms"):
            continue
        loc2glob = {ba[3]: ba[0] for ba in lig["binding_atoms"]}
        for va in vas:
            if va["rank"] == i:
                va["_glob"] = tuple(
                    sorted(loc2glob[c] for c in va["constituent_indices"] if c in loc2glob)
                )
    return vas


def _map(self, virtual_atoms, template_vectors):
    res = _o_map(self, virtual_atoms, template_vectors)
    CAP.setdefault("fits", []).append(
        {
            "mapping": res[0],
            "rmsd": res[1],
            "vecs": np.array([a["coords"] for a in virtual_atoms]),
            "template": np.array(template_vectors),
        }
    )
    return res


def _can(s):
    out = _o_can(s)
    CAP.setdefault("postpass", []).append((s, out))
    return out


X2M._align_to_pai = _pai
OA.OINDiscreteAligner._reduce_hapticity = _red
OA.OINDiscreteAligner._map_to_template = _map
CS.canonicalize_oin_slots = _can

TEMPLATE_ARRAYS = {k: np.array(v) for k, v in OA.TEMPLATES.items()}


def encode(conv, path):
    CAP.clear()
    s = conv.convert(path)
    return s, {k: (list(v) if isinstance(v, list) else dict(v)) for k, v in CAP.items()}


def selected_fit(cap, geo):
    """The last fit whose template is the SELECTED geometry's.

    ``_match_geometry_candidates`` evaluates every candidate template of the same
    coordination number (SPL/TET/TPY all have 4 vertices), so the template array itself --
    not its length -- is what identifies the winner.
    """
    t = TEMPLATE_ARRAYS.get(geo)
    if t is None:
        return None
    for f in reversed(cap.get("fits", [])):
        if f["template"].shape == t.shape and np.allclose(f["template"], t):
            return f
    return None


def exact_degeneracy(f, tol=1e-9):
    """(best rssd, set of permutations tied to ``tol``, nearest non-tied rssd).

    Scores every permutation with ``scipy Rotation.align_vectors``, the same call the
    encoder makes, so the result is independent of ``_candidate_permutations``' prefilter.
    """
    v, t = f["vecs"], f["template"]
    norms = v / (np.linalg.norm(v, axis=1)[:, None] + 1e-9)
    rows = []
    for p in itertools.permutations(range(len(t)), len(v)):
        try:
            _R, r = Rotation.align_vectors(t[list(p)], norms)
        except Exception:  # noqa: BLE001
            continue
        rows.append((p, float(r)))
    best = min(r for _p, r in rows)
    tied = {p for p, r in rows if r - best <= tol}
    nxt = min([r for p, r in rows if p not in tied], default=None)
    return best, tied, nxt


def donor_vertices(cap, geo, xyz_order=None):
    """{tuple of ORIGINAL xyz atom indices : template vertex}, plus the fit.

    ``xyz_order[i]`` is the original index of the variant file's atom ``i``, so passing it
    translates a renumbered encoding back onto the original atom numbering.
    """
    f = selected_fit(cap, geo)
    if f is None:
        return None, None
    perm = None if os.environ.get("GT_RAW") else _postpass_map(cap)
    oi = cap["origidx"]
    out = {}
    for vertex, va in enumerate(f["mapping"]):
        if va is None:
            continue
        orig = [oi.get(g, g) for g in (va.get("_glob") or ())]
        if xyz_order is not None:
            orig = [xyz_order[o] for o in orig]
        out[tuple(sorted(orig))] = perm.get(vertex, vertex) if perm else vertex
    return out, f


def _postpass_map(cap):
    pp = cap.get("postpass") or []
    return CS.canonical_slot_map(pp[-1][0]) if pp else None


def report(mol, path, ntrials=12, seed=1234):
    conv = XYZToSMILES()
    syms, coords, comment = read_xyz(path)
    base, cap = encode(conv, path)
    m = re.search(r"\[[A-Za-z]+_([A-Z]{2,3})[\]_]", base or "")
    geo = m.group(1) if m else None
    print(f"\n{'=' * 100}\n### {mol}   natoms={len(syms)}   geometry={geo}")
    print(f"base: {base}")
    bmap, bf = donor_vertices(cap, geo)
    if bmap is None:
        print("!! selected fit not recovered -- cannot judge")
        return {"molecule": mol, "error": "fit_not_recovered"}
    grp = [tuple(x) for x in (CS.geometry_rotation_group(geo) or [])]
    best, tied, nxt = exact_degeneracy(bf)
    print(
        f"FIT  best rssd={best:.6e}   perms tied to 1e-9: {len(tied)}   |rot group|={len(grp)}"
        + (f"   nearest NON-tied={nxt:.6e} (gap {nxt - best:.3e})" if nxt is not None else "")
    )
    pp = cap.get("postpass") or []
    if pp:
        print(f"pre-post-pass: {pp[-1][0]}")
    print("BASE  physical donor (original xyz idx) -> template vertex:")
    for k, v in sorted(bmap.items(), key=lambda kv: kv[1]):
        print(f"      vertex {v}: atoms {k}  [{','.join(syms[i] for i in k)}]")

    rng = random.Random(seed)
    tmp = tempfile.mkdtemp()
    seen, verdicts = {}, []
    for t in range(ntrials):
        order = list(range(len(syms)))
        rng.shuffle(order)
        p = os.path.join(tmp, f"r{t}.xyz")
        write_xyz(p, [syms[i] for i in order], coords[order], comment)
        try:
            got, cap2 = encode(conv, p)
        except Exception:  # noqa: BLE001
            continue
        if not got or got == base or got in seen:
            continue
        seen[got] = True
        print(f"--- variant {len(seen)}  (renumbering {t})")
        print(f"got : {got}")
        gmap, gf = donor_vertices(cap2, geo, xyz_order=order)
        if gmap is None:
            print("      geometry tag changed -- not a pure relabeling")
            verdicts.append("geometry_changed")
            continue
        if set(gmap) != set(bmap):
            print("      !! donor GROUPING differs -- perception moved which atoms are donors")
            verdicts.append("donor_grouping_changed")
            continue
        pi = {bmap[k]: gmap[k] for k in bmap}
        in_group = any(all(g[k] == v for k, v in pi.items()) for g in grp)
        b2, _t2, _n2 = exact_degeneracy(gf)
        print(f"      pi (base vertex -> variant vertex, per PHYSICAL donor) = {pi}")
        print(f"      identity: {all(k == v for k, v in pi.items())}")
        print(f"      extends to a PROPER ROTATION of {geo}: {in_group}")
        print(f"      variant best rssd={b2:.6e}   |delta vs base|={abs(b2 - best):.3e}")
        if not in_group:
            print("      !! SOUNDNESS DEFECT -- this relabeling is not a rotation")
        verdicts.append(
            "identity"
            if all(k == v for k, v in pi.items())
            else ("in_rotation_group" if in_group else "NOT_A_ROTATION")
        )
    if not seen:
        print(f"      no drift reproduced in {ntrials} renumberings")
    print(f"  VERDICTS {mol}: {verdicts}")
    return {
        "molecule": mol,
        "geometry": geo,
        "base": base,
        "best_rssd": best,
        "n_tied": len(tied),
        "group_order": len(grp),
        "nearest_non_tied_gap": (nxt - best) if nxt is not None else None,
        "verdicts": verdicts,
        "variants": list(seen),
    }


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    out = []
    for mol, sub in DEFAULT_MOLS:
        if which not in ("all", mol):
            continue
        path = os.path.join(DEFAULT_DATASET, sub, f"{mol}.xyz")
        if not os.path.isfile(path):
            print(f"### {mol}: MISSING {path} -- the dataset is gitignored; set OIN_DATASET")
            out.append({"molecule": mol, "error": "missing_input"})
            continue
        try:
            out.append(report(mol, path, nt))
        except Exception as e:  # noqa: BLE001
            print(f"### {mol}: FAILED {type(e).__name__}: {e}")
            out.append({"molecule": mol, "error": f"{type(e).__name__}: {e}"})
    tally: dict = {}
    for r in out:
        for v in r.get("verdicts", []) or ([r["error"]] if "error" in r else []):
            tally[v] = tally.get(v, 0) + 1
    print(f"\n{'=' * 100}\nVERDICT TALLY over {len(out)} molecules: {tally}")
    unsound = [r["molecule"] for r in out if "NOT_A_ROTATION" in (r.get("verdicts") or [])]
    print(f"SOUNDNESS DEFECTS (pi not a proper rotation): {unsound or 'NONE'}")
    if os.environ.get("GT_OUT"):
        with open(os.environ["GT_OUT"], "w") as fh:
            json.dump(out, fh, indent=2, default=str)


if __name__ == "__main__":
    main()
