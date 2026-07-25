"""Attribute every ``rdkit_canonical`` round-trip drift to exactly one root cause.

WHY THIS EXISTS (read before trusting the 500)
==============================================
``tools/roundtrip_bucket_report.py::_key_equal_subclass`` is a **cascade whose last
branch is a fallthrough**: a byte-different / key-equal pair is labelled
``rdkit_canonical`` only because it was not ``fragment_reorder``, not ``slot_renumber``
and not ``winding_star_drift``. Those three tests each require the drift to be
*exclusively* of one kind, so any pair that mixes (say) a renumbered slot with a
reordered fragment falls through to ``rdkit_canonical`` even though **no ligand body
changed at all**. The label is therefore an upper bound on ligand-body drift, not a
measurement of it.

This tool measures it. For each pair it strips the ``{n}`` slot markers, takes the
multiset of ligand *bodies* on each side, and cancels the common ones. What is left is
the drift, and it is attributed to exactly one cause:

``slot_or_order``
    The body multisets are **identical**. Nothing about the ligand bodies drifted; the
    string differs only in slot numbers and/or fragment order. **Lane 2**, not Lane 1.
``reparse_fixable``
    ``_canonical_fragment_smiles`` (``MolFromSmiles`` -> ``MolToSmiles``, with the
    ``_NO_KEKULIZE`` retry) collapses them: aromatic-vs-Kekule, implicit-vs-explicit H,
    bare ring carbon vs ``[CH2]``. **Lane 1, closed by promoting the reparse.**
``ez_chelate_locked``
    Not reparse-fixable, but ``_chelate_locked_fragment_key`` collapses them: a ``/``
    or ``\\`` on a double bond held rigid by a ring closing through the metal.
    **Lane 1, closed by the E/Z clear.**
``resonance_charge``
    Neither collapses it, but the fragment InChIs agree -- same molecule, different
    resonance form or charge site. A reparse **cannot** fix this: ``MolToSmiles`` is
    faithful to the resonance form it is handed. Needs upstream work in ``AC2BO``.
``formal_charge_placement``
    The formal-charge multiset differs (and InChI did not agree/resolve).
``connectivity``
    The heavy-atom graph genuinely differs between the two geometries. This is 3D
    *perception* drift near the covalent-radius cutoff, **not** a serialization problem.
    Explicitly OUT OF SCOPE for Lane 1 -- reported, never folded.
``unattributed``
    Anything left. Should be ~0; if it is not, the taxonomy is incomplete.

Generator-free and reads only stored strings, so it costs seconds and is immune to the
timeout confound in ``spec/handoffs/v0.4.5/BASELINE.md`` §6.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/diagnose_body_drift.py \\
        --bucket-report tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042/bucket_report.json
"""

import argparse
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from rdkit import Chem, RDLogger  # noqa: E402

from oinsmiles.oin.compare import (  # noqa: E402
    _SLOT_RE,
    _canonical_fragment_smiles,
    _chelate_locked_fragment_key,
    _parse_fragment,
    normalize_oin_for_comparison,
)

RDLogger.DisableLog("rdApp.*")

CAUSES = [
    "slot_or_order",
    "reparse_fixable",
    "ez_chelate_locked",
    "resonance_charge",
    "formal_charge_placement",
    "connectivity",
    "unattributed",
]


def _fragments(oin):
    return [f for f in normalize_oin_for_comparison((oin or "").strip()).split(".") if f]


def _body(frag):
    return _SLOT_RE.sub("", frag).strip()


def _multiset_diff(a, b):
    """(a-only, b-only) as sorted lists, common elements cancelled."""
    ca, cb = Counter(a), Counter(b)
    return sorted((ca - cb).elements()), sorted((cb - ca).elements())


def _inchi(body):
    mol = _parse_fragment(body)
    if mol is None:
        return None
    try:
        key = Chem.MolToInchiKey(mol)
    except Exception:
        return None
    return key or None


def _charge_multiset(body):
    mol = _parse_fragment(body)
    if mol is None:
        return None
    return tuple(sorted(a.GetFormalCharge() for a in mol.GetAtoms()))


def _skeleton(body):
    """Canonical SMILES of the heavy-atom graph with bond orders, charges and H's erased.

    Two bodies with the same skeleton have the same connectivity, so any remaining
    difference is a bond-order / charge / H-count question, never a perception question.
    """
    mol = _parse_fragment(body)
    if mol is None:
        return None
    try:
        rw = Chem.RWMol(mol)
        for bond in rw.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
        for atom in rw.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
            atom.SetNumRadicalElectrons(0)
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        out = rw.GetMol()
        out.UpdatePropertyCache(strict=False)
        return Chem.MolToSmiles(out, canonical=True)
    except Exception:
        return None


def _agree(fn, xs, ys):
    """True when fn maps both lists onto the same multiset, with no None results."""
    a = [fn(x) for x in xs]
    b = [fn(y) for y in ys]
    if any(v is None for v in a + b):
        return False
    return sorted(a) == sorted(b)


def _residual_after_body_fix(f1, f2):
    """Which drift survives canonicalizing every ligand body in place.

    Each fragment is reduced to ``(canonical body, its multiset of slot markers)``.
    If the two SEQUENCES agree, Lane 1 alone turns the row byte-exact: both sides now
    hold the same canonical graph, so the same donor atom lands at the same SMILES
    position and the ``{n}`` marker with it. If they agree only as multisets, fragment
    ORDER also drifted; anything else is a slot renumber. Both residuals are Lane 2.

    Caveat: this checks the marker *set* per fragment, not the marker's character
    offset inside the body. That is sound only because the reparse makes both sides the
    same graph -- it is an estimate, not a proof, and the flag-ON re-encode arm is the
    real measurement.
    """

    def sig(frags):
        return [
            (
                _canonical_fragment_smiles(_body(fr)),
                tuple(sorted(m.group(0) for m in _SLOT_RE.finditer(fr))),
            )
            for fr in frags
        ]

    a, b = sig(f1), sig(f2)
    if a == b:
        return "none"
    if sorted(a) == sorted(b):
        return "fragment_order"
    return "slot_renumber"


def attribute(s1, s2):
    """Return (cause, detail) for one byte-different / key-equal pair."""
    f1, f2 = _fragments(s1), _fragments(s2)
    b1, b2 = [_body(f) for f in f1], [_body(f) for f in f2]

    only1, only2 = _multiset_diff(b1, b2)
    detail = {
        "bodies_1": only1,
        "bodies_2": only2,
        "residual_after_body_fix": _residual_after_body_fix(f1, f2),
    }

    if not only1 and not only2:
        return "slot_or_order", detail

    # Lane 1 tier 1: the reparse. Compare the FULL canonical multiset, not just the
    # residual -- a reparse can merge two bodies that cancelled differently above.
    if _agree(_canonical_fragment_smiles, b1, b2):
        return "reparse_fixable", detail

    # Lane 1 tier 2: chelate-locked E/Z. Operates on fragments WITH their {n} markers.
    if _agree(_chelate_locked_fragment_key, f1, f2):
        return "ez_chelate_locked", detail

    # Same molecule, different resonance form / charge site: a reparse cannot fix it.
    if _agree(_inchi, b1, b2):
        return "resonance_charge", detail

    c1 = [_charge_multiset(x) for x in b1]
    c2 = [_charge_multiset(x) for x in b2]
    if not any(v is None for v in c1 + c2) and sorted(c1) != sorted(c2):
        # Charges differ AND InChI did not reconcile them. If the skeleton still
        # agrees this is a placement problem; if not, it is really connectivity.
        if _agree(_skeleton, b1, b2):
            return "formal_charge_placement", detail
        return "connectivity", detail

    if not _agree(_skeleton, b1, b2):
        return "connectivity", detail

    return "unattributed", detail


def main():
    default_report = os.path.join(
        "tmCAT-tmPHOTO_xyz_dataset", "results-capstone-v042", "bucket_report.json"
    )
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bucket-report", default=default_report, help="bucket_report.json to read")
    ap.add_argument(
        "--subclass",
        default="rdkit_canonical",
        help="key_equal subclass to attribute ('all' for every key_equal row)",
    )
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default=None, help="Write per-molecule attribution JSON here")
    ap.add_argument(
        "--examples", type=int, default=3, help="Example pairs to print per cause (0 = none)"
    )
    args = ap.parse_args()

    with open(args.bucket_report) as f:
        rows = json.load(f)

    target = [
        r
        for r in rows
        if r.get("bucket") == "key_equal"
        and (args.subclass == "all" or r.get("subclass") == args.subclass)
    ]
    if args.limit:
        target = target[: args.limit]

    counts = Counter()
    per_mol = []
    examples = {}
    for r in target:
        cause, detail = attribute(r.get("smiles_1"), r.get("smiles_2"))
        counts[cause] += 1
        per_mol.append({"molecule": r.get("molecule"), "cause": cause, **detail})
        examples.setdefault(cause, []).append(r)

    total = max(len(target), 1)
    print(f"{len(target)} key_equal/{args.subclass} rows from {args.bucket_report}\n")
    print(f"{'cause':26s} {'count':>6s} {'share':>8s}")
    print("-" * 42)
    for cause in CAUSES:
        n = counts.get(cause, 0)
        print(f"{cause:26s} {n:6d} {100 * n / total:7.2f}%")
    for cause in sorted(set(counts) - set(CAUSES)):
        print(f"{cause:26s} {counts[cause]:6d} {100 * counts[cause] / total:7.2f}%")
    print("-" * 42)
    print(f"{'total':26s} {len(target):6d}")

    lane1 = counts.get("reparse_fixable", 0) + counts.get("ez_chelate_locked", 0)
    print(f"\nLane 1 addressable (reparse + chelate E/Z): {lane1} ({100 * lane1 / total:.2f}%)")
    print(f"Lane 2 (slot/order only, no body drift):   {counts.get('slot_or_order', 0)}")
    oos = counts.get("connectivity", 0) + counts.get("unattributed", 0)
    print(f"Out of scope (connectivity + unattributed): {oos}")

    print("\nResidual drift after canonicalizing bodies (cross-tab; Lane 2 owns non-'none'):")
    print(f"{'cause':26s} {'none':>6s} {'frag_order':>11s} {'slot_renum':>11s}")
    print("-" * 58)
    for cause in CAUSES:
        rs = [m for m in per_mol if m["cause"] == cause]
        if not rs:
            continue
        res = Counter(m["residual_after_body_fix"] for m in rs)
        print(
            f"{cause:26s} {res.get('none', 0):6d} "
            f"{res.get('fragment_order', 0):11d} {res.get('slot_renumber', 0):11d}"
        )
    closed = sum(
        1
        for m in per_mol
        if m["cause"] in ("reparse_fixable", "ez_chelate_locked")
        and m["residual_after_body_fix"] == "none"
    )
    print(f"\nRows Lane 1 ALONE should turn byte_exact: {closed} ({100 * closed / total:.2f}%)")

    if args.examples:
        for cause in CAUSES:
            rs = examples.get(cause, [])[: args.examples]
            if not rs:
                continue
            print(f"\n### {cause} -- {counts[cause]} rows, {len(rs)} shown")
            for r in rs:
                print(f"  {r['molecule']}")
                print(f"    1: {r['smiles_1']}")
                print(f"    2: {r['smiles_2']}")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as f:
            json.dump({"counts": dict(counts), "molecules": per_mol}, f, indent=1)
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
