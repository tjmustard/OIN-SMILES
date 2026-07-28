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

TWO POPULATIONS -- and they are not the same question
=====================================================
``--probe-dir`` reads a ``canonicality_probe.py --out`` directory: pairs built by
*re-presenting one molecule* (renumber, rotate). That is the population every figure in
v0.4.5 was measured on.

``--roundtrip`` reads a sweep's ``individual_reports/*.json`` and pairs ``smiles_1`` (the
input encode) against ``smiles_2_indep`` (the honest re-encode of the *generated* XYZ). This
is the population the round-trip bucket report counts, and until v0.4.11 the taxonomy had
**never been run on it.** The distinction is not cosmetic: a generated conformer can place a
donor on a genuinely different vertex (``diff_occupancy``), which is upstream and unreachable
by any fold, whereas a re-presentation cannot. Measured 2026-07-27 over all 496:
``diff_occupancy`` is **0** and ``same_vcolor_identical`` is **496/496**, so the two
populations agree in shape after all -- but that was a measurement, not a foregone conclusion.

Selection uses the shipped predicate, ``roundtrip_bucket_report.classify(rep, score=...)``,
rather than a local re-implementation of ``blank_num``: two classifiers that must agree are
two classifiers that will drift.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/slot_drift_mechanism.py <probe-out-dir> [-v]
    PYTHONPATH=src .venv/bin/python tools/slot_drift_mechanism.py \\
        --roundtrip <sweep-results-dir> --expect 496 [-v]
"""

import argparse
import glob
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
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


def default_ranks(mol):
    """``breakTies=False`` symmetry classes of the fragment exactly as it was written."""
    try:
        return list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    except Exception:  # noqa: BLE001
        return None


def flattened_ranks(mol):
    """Symmetry classes of a **resonance-insensitive** copy: every bond single, charges zeroed.

    Answers a question the written form cannot: *are these two donors distinct only because
    perception froze one localized resonance form?* A beta-diketonate is the archetype -- acac
    binds through two chemically equivalent oxygens, but ``CC(=O)C=C(C)O`` writes one as a
    ketone and one as an enol, so ``CanonicalRankAtoms`` on the written body puts them in
    different classes and ``atom_verdict`` reports ``distinct_donors_LOCAL``.

    ⚠ **This is a diagnostic, not a folding criterion.** Flattening also erases distinctions
    that are real (an amide N against an amine N), so a fold keyed on it would merge donors
    that genuinely differ -- exactly what the v0.4.11 Rule forbids. Its use is to say *where*
    the residual belongs: a pair that is equivalent only after flattening is a **ligand-body
    canonicalization** problem (``rdkit_canonical``, v0.4.14), not a slot-fold problem.
    """
    em = Chem.RWMol(mol)
    for b in em.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in em.GetAtoms():
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
        a.SetIsAromatic(False)
    m = em.GetMol()
    try:
        Chem.SanitizeMol(
            m, Chem.SanitizeFlags.SANITIZE_SYMMRINGS | Chem.SanitizeFlags.SANITIZE_ADJUSTHS
        )
        return list(Chem.CanonicalRankAtoms(m, breakTies=False))
    except Exception:  # noqa: BLE001
        return None


def atom_verdict(base, got, ranker=default_ranks):
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
        ranks = ranker(mol)
        if ranks is None:
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


def iter_probe_pairs(probe_dir, subclass):
    """``(molecule, mode, base, got)`` for each re-presentation pair labelled ``subclass``."""
    path = os.path.join(probe_dir, "canonicality_probe_all.json")
    with open(path) as fh:
        data = json.load(fh)
    recs = data["records"] if isinstance(data, dict) and "records" in data else data
    for r in recs:
        if not r.get("base") or r.get("stable"):
            continue
        for v in r["variants"]:
            if v.get("subclass") != subclass or not v.get("got"):
                continue
            yield r["molecule"], v.get("mode", ""), r["base"], v["got"]
            break  # one pair per molecule, matching the v0.4.5 counts


def iter_roundtrip_pairs(results_dir, subclass, score="honest"):
    """``(molecule, mode, base, got)`` for each round-trip pair the bucket report labels ``subclass``.

    ``base`` is ``smiles_1`` (the input encode) and ``got`` is ``smiles_2_indep`` under the
    honest score -- a full re-perception of the *generated* XYZ, not the generator's own bond
    graph. Selection delegates to ``roundtrip_bucket_report.classify`` so this tool and the
    bucket report can never disagree about what is in the bucket.
    """
    import roundtrip_bucket_report as rbr

    key = "smiles_2_indep" if score == "honest" else "smiles_2"
    for fp in sorted(glob.glob(os.path.join(results_dir, "individual_reports", "*.json"))):
        try:
            with open(fp) as fh:
                rep = json.load(fh)
        except Exception:  # noqa: BLE001  -- a corrupt report is not this tool's business
            continue
        bucket, sub = rbr.classify(rep, score=score)
        if bucket != "key_equal" or sub != subclass:
            continue
        yield os.path.basename(fp)[:-5], score, rep.get("smiles_1"), rep.get(key)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "probe_dir",
        nargs="?",
        help="a canonicality_probe.py --out directory (re-presentation pairs)",
    )
    ap.add_argument(
        "--roundtrip",
        metavar="RESULTS_DIR",
        help="a sweep results dir; pairs smiles_1 vs smiles_2_indep from individual_reports/",
    )
    ap.add_argument(
        "--score",
        default="honest",
        choices=("honest", "scored"),
        help="--roundtrip only: which round-trip string to read (default honest)",
    )
    ap.add_argument("--subclass", default="slot_renumber", help="probe subclass to explain")
    ap.add_argument(
        "--expect",
        type=int,
        help="abort unless exactly N pairs are selected -- pin the population you meant to measure",
    )
    ap.add_argument("--out-json", help="write the per-molecule classification here")
    ap.add_argument(
        "--explain-distinct",
        action="store_true",
        help="re-test each distinct_donors_LOCAL pair resonance-insensitively, to separate "
        "a frozen resonance form (a ligand-BODY problem) from a genuine donor inequivalence",
    )
    ap.add_argument("-v", "--verbose", action="store_true", help="print each pair")
    args = ap.parse_args()

    if bool(args.probe_dir) == bool(args.roundtrip):
        ap.error("give exactly one of <probe-dir> or --roundtrip RESULTS_DIR")

    if args.roundtrip:
        pairs = list(iter_roundtrip_pairs(args.roundtrip, args.subclass, args.score))
        source = f"{args.roundtrip} (round-trip, --score {args.score})"
    else:
        pairs = list(iter_probe_pairs(args.probe_dir, args.subclass))
        source = f"{args.probe_dir} (re-presentation)"

    print(f"source: {source}")
    print(f"selected: {len(pairs)} pairs labelled {args.subclass!r}")
    if args.expect is not None and len(pairs) != args.expect:
        sys.exit(
            f"ABORT: expected {args.expect} pairs, selected {len(pairs)}. "
            "The population is not the one this run was meant to measure."
        )

    mech, verd, flagged, rows = Counter(), Counter(), [], []
    resonance = Counter()
    for molecule, mode, base, got in pairs:
        m = mechanism(base, got)
        mech[m] += 1
        a = flat = None
        if m == "same_vcolor_identical":
            a = atom_verdict(base, got)
            verd[a] += 1
            if a == "distinct_donors_LOCAL":
                flagged.append(molecule)
                if args.explain_distinct:
                    flat = atom_verdict(base, got, ranker=flattened_ranks)
                    resonance[
                        "resonance_artifact_body_problem"
                        if flat == "automorphism"
                        else "genuinely_distinct_donors"
                    ] += 1
        rows.append(
            {
                "molecule": molecule,
                "mode": mode,
                "mechanism": m,
                "atom_verdict": a,
                "flattened_verdict": flat,
                "base": base,
                "got": got,
            }
        )
        if args.verbose:
            print(f"--- {molecule} [{mode}] {m}" + (f" / {a}" if a else ""))
            print(f"  base: {base}")
            print(f"  got : {got}")

    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump({"source": source, "subclass": args.subclass, "records": rows}, fh, indent=1)
        print(f"\nwrote {args.out_json}")

    total = len(pairs) or 1
    print(f"\n{args.subclass}: mechanism")
    for k, n in mech.most_common():
        print(f"  {k:45} {n:5}  {100 * n / total:5.1f}%")
    if verd:
        n_svi = sum(verd.values())
        print(f"\nof the {n_svi} same_vcolor_identical, atom-level verdict")
        for k, n in verd.most_common():
            print(
                f"  {k:45} {n:5}  {100 * n / n_svi:5.1f}% of svi   {100 * n / total:5.1f}% of all"
            )
        reach = verd.get("automorphism", 0)
        print(
            f"\n  REACHABLE BY A WITHIN-FRAGMENT FOLD: {reach} of {total} "
            f"({100 * reach / total:.1f}%, {reach / 50.0:.2f} points at 5000 molecules)"
        )
    if resonance:
        n_dd = sum(resonance.values())
        print(f"\nof the {n_dd} distinct_donors_LOCAL, re-tested resonance-insensitively")
        for k, n in resonance.most_common():
            print(f"  {k:45} {n:5}  {100 * n / n_dd:5.1f}% of dd   {n / 50.0:.2f} points")
        print(
            "  ^ 'resonance_artifact_body_problem' pairs are equivalent once bond orders are\n"
            "    flattened, so they are a LIGAND-BODY canonicalization gap (rdkit_canonical,\n"
            "    v0.4.14) -- NOT reachable, and not a defect, at the slot-fold seam."
        )
    if flagged:
        print(f"\ndistinct_donors_LOCAL -- NOT interchangeable per-fragment ({len(flagged)}):")
        print("  " + ", ".join(flagged))
        print("  This is NOT a soundness verdict. Settle it with:")
        print("    PYTHONPATH=src GT_RAW=1 python tools/wrong_donor_groundtruth.py")
    if mech.get("postpass_BUG_diverges"):
        print("\n!! postpass_BUG_diverges > 0 -- the post-pass itself is at fault. Investigate.")


if __name__ == "__main__":
    main()
