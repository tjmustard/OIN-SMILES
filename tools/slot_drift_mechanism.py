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
``DISTINCT_donors``
    the slots land on inequivalent donor atoms, so one of the two strings is WRONG. A
    soundness defect, not a canonicality one -- and folding these would be over-folding.

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


def _slot_to_atom(frag):
    """{slot: atom index in the slot-stripped fragment SMILES}, first occurrence wins."""
    out = {}
    for m in OINInlineHandler.SLOT_REGEX.finditer(frag):
        prefix = OINInlineHandler.SLOT_REGEX.sub("", frag[: m.start()])
        out.setdefault(int(m.group(1)), _count_smiles_atoms_before(prefix, len(prefix)))
    return out


def atom_verdict(base, got):
    """Second-stage verdict: is the moved donor pair interchangeable, or genuinely distinct?"""
    frags_b = [
        f
        for f in OINInlineHandler.METAL_REGEX.sub("", normalize_oin_for_comparison(base)).split(".")
        if f.strip()
    ]
    frags_g = [
        f
        for f in OINInlineHandler.METAL_REGEX.sub("", normalize_oin_for_comparison(got)).split(".")
        if f.strip()
    ]
    verdicts = []
    for a, b in zip(frags_b, frags_g):
        sa, sb = _slot_to_atom(a), _slot_to_atom(b)
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
        ok = True
        for slot, atom in sa.items():
            other = sb.get(slot)
            if other is None or max(atom, other) >= len(ranks) or ranks[atom] != ranks[other]:
                ok = False
                break
        verdicts.append("automorphism" if ok else "DISTINCT_donors")
    if not verdicts:
        return "no_atom_level_move"
    if "DISTINCT_donors" in verdicts:
        return "DISTINCT_donors"
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
                if a == "DISTINCT_donors":
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
        print(f"\nDISTINCT_donors -- one of the two strings is WRONG ({len(flagged)}):")
        print("  " + ", ".join(flagged))
    if mech.get("postpass_BUG_diverges"):
        print("\n!! postpass_BUG_diverges > 0 -- the post-pass itself is at fault. Investigate.")


if __name__ == "__main__":
    main()
