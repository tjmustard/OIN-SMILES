"""Why does slot drift SURVIVE the canonical-slot post-pass? (v0.4.5 Lane 2)

``canonicality_probe.py`` reports *that* a molecule's ``{n}`` slot integers moved under a
re-presentation. This tool says *why*, which is the difference between "the post-pass has a
bug" and "the post-pass cannot reach this by construction". Run it on a probe output
directory.

MECHANISM CLASSES (first stage)
===============================
The post-pass derives its relabeling from ``compare._parse_vertex_colors`` alone, so:

``diff_geometry``
    different geometry tag, hence a different rotation group. Upstream (the 3D fit).
``diff_occupancy``
    a different SET of template vertices is occupied -- the arrangement moved, not just the
    labels. No rotation relates them, so lex-min legitimately differs. Upstream.
``diff_colors``
    the vertex colours differ, so the lex-min target differs. Downstream of the canonical
    ligand body (Lane 1), not of the slot lever.
``same_vcolor_identical``
    **the colored-vertex maps are IDENTICAL.** The post-pass therefore computes the identical
    permutation for both strings, and applying one permutation to two differently-labeled
    strings preserves the difference. This class is outside the post-pass's reach *by
    design*: ``_parse_vertex_colors`` colours every donor of a ligand with that ligand's
    whole body, so which donor of a chelate holds which slot is deliberately invisible to it
    (that blindness is what lets true conformers collapse in the comparison key while fac/mer
    stay distinct).
``postpass_BUG_diverges``
    same occupancy and same colour multiset, different assignment, and re-running the
    post-pass on both does not converge. This one WOULD be a defect. Measured 0.

ATOM-LEVEL VERDICT (second stage, for ``same_vcolor_identical`` only)
=====================================================================
Which donor *atom* of the fragment moved, and are the old and new atoms interchangeable?

``automorphism``
    the two atoms are in the same ``CanonicalRankAtoms(breakTies=False)`` symmetry class of
    their own fragment, so both labelings denote the SAME molecule and the encoder merely
    needs a deterministic choice. (``breakTies=False`` is the right instrument here:
    ``breakTies=True`` settles ties between symmetry-equivalent atoms on the *input* index,
    which is the very dependence being measured. ``includeChirality`` defaults to True, so
    two constitutionally-equivalent branches with different configurations land in
    *different* classes and are correctly NOT called interchangeable.)
``distinct_donors_LOCAL``
    the slots land on donor atoms that this *per-fragment, string-only* test cannot show to
    be interchangeable. **It does NOT mean one of the two strings is wrong** -- see the
    warning on ``atom_verdict``. v0.4.5 Lane 9 settled all 7 of these from the 3D coordinates
    (``tools/wrong_donor_groundtruth.py``) and found **0 soundness defects**: the slot
    labeling is only ever determined up to the polyhedron's proper-rotation group, and a
    rotation acts on every fragment at once, which a per-fragment test cannot see.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/slot_drift_mechanism.py <probe-out-dir> [-v]
"""

import argparse
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from rdkit import Chem, RDLogger  # noqa: E402

from oinsmiles.oin.canonical_slots import canonicalize_oin_slots  # noqa: E402
from oinsmiles.oin.compare import (  # noqa: E402
    _parse_fragment,
    _parse_vertex_colors,
    normalize_oin_for_comparison,
)
from oinsmiles.oin.inline import (  # noqa: E402
    OINInlineHandler,
    _count_smiles_atoms_before,
)

RDLogger.DisableLog("rdApp.*")


def mechanism(base, got):
    """First-stage class for one (base, got) pair."""
    _mb, geo_b, col_b = _parse_vertex_colors(normalize_oin_for_comparison(base))
    _mg, geo_g, col_g = _parse_vertex_colors(normalize_oin_for_comparison(got))
    if geo_b != geo_g:
        return "diff_geometry"
    if set(col_b) != set(col_g):
        return "diff_occupancy"
    if sorted(col_b.values()) != sorted(col_g.values()):
        return "diff_colors"
    if col_b == col_g:
        return "same_vcolor_identical"
    if canonicalize_oin_slots(base) == canonicalize_oin_slots(got):
        return "postpass_converges_but_encoder_did_not"
    return "postpass_BUG_diverges"


def _fragments(oin):
    return [
        f
        for f in OINInlineHandler.METAL_REGEX.sub("", normalize_oin_for_comparison(oin)).split(".")
        if f.strip()
    ]


def _slot_to_atoms(frag):
    """{slot: frozenset of atom indices carrying it} -- ALL occurrences, not just the first.

    A haptic donor stamps its slot on every ring atom (``c{0}1c{0}(C)...``). Keying on the
    *first* occurrence compares whichever ring position happened to be written first, so two
    genuinely interchangeable eta rings can present a methyl-bearing carbon against a
    silyl-bearing one and land in different symmetry classes for no chemical reason. That was
    a measured false positive (ZACFER_comp_0, v0.4.5 Lane 9).
    """
    out: dict[int, set] = {}
    for m in OINInlineHandler.SLOT_REGEX.finditer(frag):
        prefix = OINInlineHandler.SLOT_REGEX.sub("", frag[: m.start()])
        out.setdefault(int(m.group(1)), set()).add(_count_smiles_atoms_before(prefix, len(prefix)))
    return {k: frozenset(v) for k, v in out.items()}


def _pair_fragments(frags_b, frags_g):
    """Pair base fragments with got fragments by BODY TEXT, not by position.

    ``zip(frags_b, frags_g)`` is wrong: the canonical-slot post-pass re-derives fragment
    order from the slot integers, so a molecule whose slots moved generally has its
    fragments in a different order too -- and a complex with two copies of one ligand (two
    dppe, two CO) then gets ligand A compared against ligand B. Measured false positive:
    ZOSNUS_comp_0 (v0.4.5 Lane 9).
    """
    remaining = list(range(len(frags_g)))
    pairs = []
    for a in frags_b:
        body_a = OINInlineHandler.SLOT_REGEX.sub("", a)
        hit = next(
            (j for j in remaining if OINInlineHandler.SLOT_REGEX.sub("", frags_g[j]) == body_a),
            None,
        )
        if hit is None:
            return None  # bodies differ -- caller's `diff_colors` stage should have caught it
        remaining.remove(hit)
        pairs.append((a, frags_g[hit]))
    return pairs


def atom_verdict(base, got):
    """Is the moved donor pair interchangeable *within its own fragment*, or not?

    .. warning::
       **A ``distinct_donors_LOCAL`` verdict does NOT mean one of the two strings is wrong.**
       This is a per-fragment, string-only heuristic, and the question "which physical atom
       sits at which template vertex" is a property of the 3D coordinates that no amount of
       string comparison can settle. Two facts make the local answer unreliable:

       * the emitted slot labeling is only ever determined **up to the coordination
         polyhedron's proper-rotation group** -- the Kabsch fit is exactly degenerate over it
         -- so a relabeling can be a benign change of reference frame that this per-fragment
         test cannot see;
       * a rotation acts on **all** fragments at once. ``RUBTIS_comp_0``'s chelate donors are
         genuinely inequivalent (imine N vs pyridyl N), yet the full relabeling is
         ``(0 1)(2 3)`` in D4 because the COD's two equivalent alkene arms swapped as well.

       v0.4.5 Lane 9 settled all 7 of this tool's ``DISTINCT_donors`` molecules from the 3D
       coordinates with ``tools/wrong_donor_groundtruth.py``: **0 of 7 were soundness
       defects** (4 had a bit-identical donor->vertex map; the other 3 differed by a proper
       rotation, ``|delta rssd| <= 1.2e-14``). Treat this verdict as "needs the 3D
       instrument", never as "the encoder emitted a wrong string".
    """
    frags_b, frags_g = _fragments(base), _fragments(got)
    pairs = _pair_fragments(frags_b, frags_g)
    if pairs is None:
        return "fragment_bodies_differ"
    verdicts = []
    for a, b in pairs:
        sa, sb = _slot_to_atoms(a), _slot_to_atoms(b)
        if sa == sb:
            continue
        mol = _parse_fragment(OINInlineHandler.SLOT_REGEX.sub("", a))
        if mol is None:
            verdicts.append("unparsable")
            continue
        try:
            ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
        except Exception:  # noqa: BLE001
            verdicts.append("unparsable")
            continue

        def _classes(atoms, ranks=ranks):
            if any(i >= len(ranks) for i in atoms):
                return None
            return sorted(ranks[i] for i in atoms)

        ok = True
        for slot, atoms in sa.items():
            other = sb.get(slot)
            ca, cb = _classes(atoms), (_classes(other) if other is not None else None)
            if cb is None or ca is None or ca != cb:
                ok = False
                break
        verdicts.append("automorphism" if ok else "distinct_donors_LOCAL")
    if not verdicts:
        return "no_atom_level_move"
    if "distinct_donors_LOCAL" in verdicts:
        return "distinct_donors_LOCAL"
    if all(v == "automorphism" for v in verdicts):
        return "automorphism"
    return "+".join(sorted(set(verdicts)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("probe_dir", help="a canonicality_probe.py --out directory")
    ap.add_argument("--subclass", default="slot_renumber", help="probe subclass to explain")
    ap.add_argument("-v", "--verbose", action="store_true", help="print each pair")
    args = ap.parse_args()

    path = os.path.join(args.probe_dir, "canonicality_probe_all.json")
    with open(path) as fh:
        data = json.load(fh)
    recs = data["records"] if isinstance(data, dict) and "records" in data else data

    mech, verd, flagged = Counter(), Counter(), []
    for r in recs:
        if not r.get("base") or r.get("stable"):
            continue
        for v in r["variants"]:
            if v.get("subclass") != args.subclass or not v.get("got"):
                continue
            m = mechanism(r["base"], v["got"])
            mech[m] += 1
            if m == "same_vcolor_identical":
                a = atom_verdict(r["base"], v["got"])
                verd[a] += 1
                if a == "distinct_donors_LOCAL":
                    flagged.append(r["molecule"])
            if args.verbose:
                print(f"--- {r['molecule']} [{v['mode']}] {m}")
                print(f"  base: {r['base']}")
                print(f"  got : {v['got']}")
            break

    print(f"\n{args.subclass}: mechanism")
    for k, n in mech.most_common():
        print(f"  {k:45} {n}")
    if verd:
        print("\nof the same_vcolor_identical, atom-level verdict")
        for k, n in verd.most_common():
            print(f"  {k:45} {n}")
    if flagged:
        print(f"\ndistinct_donors_LOCAL -- NOT interchangeable per-fragment ({len(flagged)}):")
        print("  " + ", ".join(flagged))
        print("  This is NOT a soundness verdict. Settle it with:")
        print("    PYTHONPATH=src GT_RAW=1 python tools/wrong_donor_groundtruth.py")
    if mech.get("postpass_BUG_diverges"):
        print("\n!! postpass_BUG_diverges > 0 -- the post-pass itself is at fault. Investigate.")


if __name__ == "__main__":
    main()
