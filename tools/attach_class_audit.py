"""How big are the MEDZUR and GAVSED classes? (v0.4.13 Lane 2)

Both were handed to this release known on **n = 1 each**, from the v0.4.7 attachment lane:

* **GAVSED shape** -- acceptance rejected every conformer, so ``_select_by_geometry``'s
  geometry-ranked fallback returned one anyway, **and that fallback is not attachment-aware**.
  The returned structure has ligands off the metal. The attachment check guards *acceptance*,
  not *return*, so this is unreachable by tightening the check.
* **MEDZUR shape** -- attachment fully intact and independent re-perception still disagrees.
  Predicted by no model, found by that lane, and left unexplained. The attachment check is
  simply not the binding constraint here.

WHY THIS RUNS OFFLINE, AND WHY THAT IS NOT A SHORTCUT
=====================================================
v0.4.8 backfilled a ``coordination`` block into every individual report -- ``intact``,
``boundary_only``, and per-metal contact counts computed from the STORED generated geometry.
The classification this tool needs is therefore already a property of the frozen corpus, and
re-running the generator would not make it more true; it would make it a different corpus.
So this is a JSON join, it runs in seconds, and it is reproducible by anyone with the frozen
directory.

THE CONTROL ARM IS THE POINT, NOT A COURTESY
============================================
The standing rule this project pays for repeatedly:

    Before quoting an instrument, ask what a BROKEN version of it would print. If that is the
    same thing, you have not measured anything yet.

A broken version of this tool -- one where ``coordination.intact`` were mostly ``False`` for
reasons unrelated to failure, or mostly ``None`` and silently coerced -- would print a large,
plausible, entirely meaningless GAVSED class. Two things make the output falsifiable:

1. **The ``byte_exact`` control.** The same split is computed over molecules that round-trip
   perfectly. If ``DETACHED`` is common *there*, then detachment does not discriminate and the
   GAVSED number means nothing. The release must read the two side by side or not at all.
2. **``UNKNOWN`` is never folded into another class.** A missing or unparseable ``coordination``
   block is its own bucket with its own count. Every table prints its denominator.

Classification, in priority order (first match wins):

    NO_STRUCTURE  the run produced no structure to judge (status != success, or smiles_2 None).
                  Neither class -- there is no returned conformer whose attachment could be wrong.
    UNKNOWN       a structure exists but carries no usable `coordination` block.
    DETACHED      coordination.intact is False              -> the GAVSED shape
    BOUNDARY      intact, but every contact is within the cutoff tolerance (boundary_only)
    INTACT        intact and not boundary                   -> the MEDZUR shape

BOUNDARY is kept separate from INTACT deliberately. `coordination.boundary_only` means the
verdict rests on contacts sitting within 0.1 A of the cutoff, i.e. the attachment call itself is
the uncertain quantity. Merging it into either neighbour would let a tolerance choice decide the
headline.

Usage
-----
    V=$PWD/.venv/bin/python; export PYTHONPATH=$PWD/src
    $V tools/attach_class_audit.py \
        --results-dir tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest \
        --out docs/agentic-notes/v0.4.13/attach_class_audit.json

Reads ``bucket_report_honest.json`` from ``--results-dir`` for the AUTHORITATIVE bucket, never
an ad-hoc ``honest_class.endswith("->byte")`` -- that shortcut disagrees with the bucket report
by exactly the eight atom-count-gate molecules (v0.4.9's standing trap).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

CLASSES = ("NO_STRUCTURE", "UNKNOWN", "DETACHED", "BOUNDARY", "INTACT")

# --------------------------------------------------------------------------------------------
# v0.4.16 Lane 2: WHY does the independent re-perception disagree?
#
# The attachment split above says whether the ligands are ON the metal. For 187 molecules the
# answer is yes (or boundary) and the round trip still fails, and the roadmap has carried them for
# three releases as "still unexplained" while sizing two later releases against guesses about
# them. This section replaces the guess with a classification of WHICH COMPONENT of the
# round-trip key disagrees, which is a property of the frozen corpus and needs no generation.
#
# 🔴 READ THE CONTROL FIRST, AND THE CONTROL REFUTES THE FRAMING THIS SECTION INHERITED. The
# charter describes `structural`->BOUNDARY as "the attachment call itself is inside the tolerance
# band", i.e. treats BOUNDARY as a cause. Over the byte_exact molecules -- the ones that round-trip
# PERFECTLY -- BOUNDARY is 1367 of 3858 (35.4%) and INTACT is 2431 (63.0%). Neither class
# discriminates; BOUNDARY is the modal state of a passing molecule. A bucket name that asserts a
# cause is a hypothesis, not a measurement.
# --------------------------------------------------------------------------------------------


def _skeleton(body: str, strip_stereo: bool = False):
    """A body reduced to its HEAVY-ATOM GRAPH: no bond orders, charges, H counts or aromaticity.

    The point is to separate two mechanisms that look identical in a raw string diff:

    * same graph, different decoration -> **perception**. The heavy atoms are connected exactly as
      the input had them and the re-perception assigned different bond orders / aromaticity / H
      counts. `BOVCUM_comp_0` reads `CC(=O)c1sc2ccccc2c1O` in and `CC(O)=C1Sc2ccccc2C1=O` out --
      a ketone/enol tautomer shift. No string canonicalization reaches that.
    * different graph -> **construction**. The generated geometry genuinely has different
      heavy-atom connectivity, e.g. `XOYCOE_comp_0` whose chloride fuses into a ligand.

    🔴 THIS GOES THROUGH RDKit CANONICALIZATION, AND THE FIRST VERSION DID NOT. A string-level
    normalizer (strip brackets, drop bond symbols, uppercase) is wrong **57 of 109 times** on this
    very population -- a coin flip -- because SMILES ring-closure digits and atom ordering are
    arbitrary labels, so an identical graph written two ways reads as different. It reported
    `BEDLII_comp_0` as CONSTRUCTION when its heavy graph is unchanged and the real difference is
    aromatic perception. That is a plausible, precise, entirely wrong number of exactly the kind
    this project keeps paying for: measured against a canonical comparison rather than trusted.

    Returns ``None`` when the body cannot be parsed at all -- reported as its own class, never
    silently merged into "different".
    """
    from rdkit import Chem

    body = body or ""
    # ⚠ `RAW:` is a SENTINEL the key builder prefixes when canonical_body_emit could not
    # canonicalize that body -- it is not part of the SMILES. Parsing it verbatim fails on every
    # such body, which is 23 of the 172 here and would have been reported as "RDKit cannot parse
    # this molecule" when the truth is "this tool forgot to strip a marker". Same class of error
    # as the string-level skeleton above: a confident, wrong, self-consistent number.
    if body.startswith("RAW:"):
        body = body[4:]
    mol = Chem.MolFromSmiles(body, sanitize=False)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
        atom.SetIsAromatic(False)
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    if strip_stereo:
        # 🔴 CHIRALITY IS NOT CONNECTIVITY, and conflating them costs a whole mechanism.
        # `DIVZOY_comp_0` differs only as `C[P@@](...)` vs `C[P@](...)` -- one inverted phosphorus
        # stereocentre. With tags left on, `MolToSmiles` writes them and the molecule reads as a
        # different GRAPH, i.e. as construction. It is a STEREO INVERSION, which is the same
        # family as the enantiomer class v0.4.17 owns and wants counting separately, not a broken
        # bond. Comparing both ways is what separates the two.
        for atom in mol.GetAtoms():
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        for bond in mol.GetBonds():
            bond.SetStereo(Chem.BondStereo.STEREONONE)
    try:
        return Chem.MolToSmiles(mol)
    except Exception:  # noqa: BLE001 -- an unwritable graph is data
        return None


def _key_diff_class(rep: dict, key_fn) -> tuple[str, str]:
    """``(class, detail)`` -- which component of the round-trip key disagrees, and how.

    Returns ``UNCLASSIFIED`` rather than guessing when the key cannot be computed or the
    components all agree. That bucket is printed with its count: the charter requires UNKNOWN to
    be reported rather than folded, and an unclassified molecule is exactly an admission that this
    instrument did not reach it.
    """
    s1, s2 = rep.get("smiles_1"), rep.get("smiles_2_indep")
    k1, k2 = key_fn(s1), key_fn(s2)
    if k1 is None or k2 is None:
        return "KEY_UNCOMPUTABLE", ""
    (el1, geo1), _sig1, b1 = (k1[0][0], k1[0][1]), k1[1], k1[2]
    (el2, geo2), _sig2, b2 = (k2[0][0], k2[0][1]), k2[1], k2[2]

    if el1 != el2:
        # Never yet observed. Kept because a silent absence and a measured zero are different
        # claims, and only one of them is evidence.
        return "METAL_ELEMENT", f"{el1}->{el2}"

    sb1, sb2 = sorted(b1), sorted(b2)
    if len(sb1) != len(sb2):
        extra = sorted(set(sb2) - set(sb1)) if len(sb2) > len(sb1) else sorted(set(sb1) - set(sb2))
        bare_h = [e for e in extra if e in ("[H]", "[HH]", "H", "[H][H]")]
        if bare_h and len(bare_h) == len(extra):
            # Hydrogens came off their heavy atom: DOLWEG_comp_0 emits [CH3]{1} and re-perceives
            # as [CH]{1} plus free [H]. A bond-length problem, so construction.
            cls = "H_DETACHED"
        else:
            cls = "LIGAND_SPLIT" if len(sb2) > len(sb1) else "LIGAND_MERGE"
        return cls, f"{len(sb1)}->{len(sb2)}"

    if sb1 != sb2:
        # Two comparisons, because three mechanisms hide behind one string difference. sorted()
        # on both sides: the bodies are a multiset and canonicalization can reorder them
        # relative to the raw-string sort used to pair them up.
        g1, g2 = [_skeleton(b) for b in sb1], [_skeleton(b) for b in sb2]
        f1, f2 = (
            [_skeleton(b, strip_stereo=True) for b in sb1],
            [_skeleton(b, strip_stereo=True) for b in sb2],
        )
        if None in g1 or None in g2 or None in f1 or None in f2:
            # A body RDKit cannot parse at all. Its own class: an unparseable graph is not
            # evidence of a DIFFERENT graph, and folding it into SKELETON would inflate the
            # CONSTRUCTION share with molecules this instrument simply did not reach.
            return "GRAPH_UNPARSEABLE", ""
        if sorted(f1) != sorted(f2):
            return "SKELETON", ""
        # Connectivity agrees once stereo is stripped. Either the chiral tags moved...
        if sorted(g1) != sorted(g2):
            return "STEREO_INVERSION", ""
        # ...or only bond orders / aromaticity / charges / H counts differ.
        return "PERCEPTION", ""

    if geo1 != geo2:
        # Bodies and element agree, only the perceived polyhedron moved. Reported on its own
        # rather than merged into a body class, because it is a DIFFERENT mechanism: the
        # coordination number / geometry classifier read the generated structure differently.
        return "GEOM_CODE", f"{geo1}->{geo2}"

    # Element, bodies and geometry all agree -> the signature is what differs, which is what
    # `facmer_divergent` means. Not a defect of this classifier.
    return "ARRANGEMENT_ONLY", ""


#: What could REACH each class. 🔴 Reachability is a property of a MECHANISM, not of a block --
#: "not reachable by canonicalization" was once recorded and then read as "not reachable", and a
#: 2.26-point block sat out a whole release. So each entry names the mechanism, every time.
_MECHANISM = {
    "PERCEPTION": "PERCEPTION -- heavy-atom skeleton agrees; bond orders/aromaticity/H differ",
    "SKELETON": "CONSTRUCTION -- heavy-atom GRAPH differs even ignoring stereo (RDKit canonical)",
    "STEREO_INVERSION": "CONSTRUCTION/SELECTION -- same connectivity, INVERTED stereocentre",
    "GRAPH_UNPARSEABLE": "UNCLASSIFIED -- a ligand body RDKit could not parse; NOT evidence of a diff",
    "H_DETACHED": "CONSTRUCTION -- hydrogens left their heavy atom (bond length)",
    "LIGAND_SPLIT": "CONSTRUCTION -- a bond broke, one ligand became several",
    "LIGAND_MERGE": "CONSTRUCTION -- fragments fused that the input kept apart",
    "GEOM_CODE": "PERCEPTION/SELECTION -- the coordination polyhedron was read differently",
    "METAL_ELEMENT": "PERCEPTION -- the metal itself was re-perceived as another element",
    "ARRANGEMENT_ONLY": "SELECTION -- same parts, different arrangement (the facmer axis)",
    "KEY_UNCOMPUTABLE": "UNCLASSIFIED -- the key could not be computed for one side",
}


def classify(report: dict) -> str:
    """Which attachment class does this molecule's RETURNED structure fall in?"""
    if report.get("status") != "success" or not report.get("smiles_2"):
        return "NO_STRUCTURE"
    coord = report.get("coordination")
    if not isinstance(coord, dict) or coord.get("intact") is None:
        return "UNKNOWN"
    if not coord.get("intact"):
        return "DETACHED"
    if coord.get("boundary_only"):
        return "BOUNDARY"
    return "INTACT"


def load_buckets(results_dir: Path) -> dict:
    """molecule -> (bucket, subclass) from the authoritative frozen report."""
    path = results_dir / "bucket_report_honest.json"
    if not path.exists():
        sys.exit(f"FATAL: no authoritative bucket report at {path}")
    rows = json.loads(path.read_text())
    if not isinstance(rows, list) or not rows:
        sys.exit(f"FATAL: {path} is not a non-empty list of per-molecule records")
    return {r["molecule"]: (r.get("bucket"), r.get("subclass")) for r in rows}


#: The blocks v0.4.15's reachability map left UNMEASURED -- 3.74 points, and the largest remaining
#: unknown on the board. (bucket, attach_class) -> the roadmap's label for it.
_TARGETS = (
    ("structural", "INTACT", "the MEDZUR class"),
    ("structural", "BOUNDARY", "structural / BOUNDARY"),
    ("facmer_divergent", None, "facmer_divergent (all attach classes)"),
)


def _characterise(results_dir: Path, reports_dir: Path, table) -> dict:
    """Classify every molecule in the unmeasured blocks by WHICH key component disagrees.

    Offline over the frozen sweep, for the same reason the attachment split is: the classification
    is a property of the stored strings, and re-running the generator would not make it more true.
    """
    sys.path.insert(
        0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
    )
    from oinsmiles.oin.compare import canonical_roundtrip_key

    def key_fn(s):
        try:
            return canonical_roundtrip_key(s) if s else None
        except Exception:  # noqa: BLE001 -- an uncomputable key is data, not a crash
            return None

    out = {"blocks": {}, "control_byte_exact": {}}
    print("\n\n# v0.4.16 Lane 2 -- WHY does the independent re-perception disagree?\n")

    for bucket, attach_class, label in _TARGETS:
        classes = [attach_class] if attach_class else list(CLASSES)
        molecules = [m for c in classes for m in table.get(bucket, {}).get(c, [])]
        if not molecules:
            continue
        counts: Counter = Counter()
        details: dict[str, Counter] = defaultdict(Counter)
        example: dict[str, str] = {}
        members: dict[str, list[str]] = defaultdict(list)
        for molecule in sorted(molecules):
            path = reports_dir / f"{molecule}.json"
            if not path.exists():
                counts["KEY_UNCOMPUTABLE"] += 1
                members["KEY_UNCOMPUTABLE"].append(molecule)
                continue
            cls, detail = _key_diff_class(json.loads(path.read_text()), key_fn)
            counts[cls] += 1
            members[cls].append(molecule)
            if detail:
                details[cls][detail] += 1
            example.setdefault(cls, molecule)

        total = sum(counts.values())
        print(f"## {label} -- n={total} ({bucket}/{attach_class or 'ALL'})\n")
        print(f"  {'class':<18} {'n':>5}  {'%':>6}  mechanism")
        for cls, n in counts.most_common():
            pct = f"{100.0 * n / total:.1f}%"
            print(f"  {cls:<18} {n:>5}  {pct:>6}  {_MECHANISM.get(cls, 'UNCLASSIFIED')}")
            print(f"  {'':<18} {'':>5}  {'':>6}  e.g. {example[cls]}")
            if details[cls]:
                top = ", ".join(f"{d} x{c}" for d, c in details[cls].most_common(4))
                print(f"  {'':<18} {'':>5}  {'':>6}  {top}")
        # The denominator, printed every time. A rate without its sample is not reproducible.
        assert total == len(molecules), "every member must land in exactly one class"
        print(f"\n  DENOMINATOR {total} classified, 0 unaccounted\n")
        out["blocks"][f"{bucket}/{attach_class or 'ALL'}"] = {
            "label": label,
            "n": total,
            "counts": dict(counts),
            "examples": example,
            "members": {k: sorted(v) for k, v in members.items()},
        }

    # 🔴 THE CONTROL. Without it the attachment classes above read as causes.
    ctrl = table.get("byte_exact", {})
    n_ctrl = sum(len(ctrl.get(c, [])) for c in CLASSES)
    out["control_byte_exact"] = {c: len(ctrl.get(c, [])) for c in CLASSES}
    print(
        "## CONTROL -- the same attachment classes over the molecules that round-trip PERFECTLY\n"
    )
    for c in CLASSES:
        n = len(ctrl.get(c, []))
        print(
            f"  byte_exact {c:<14} {n:>5}  ({100.0 * n / n_ctrl:.1f}%)" if n_ctrl else f"  {c}: n/a"
        )
    print(
        "\n  🔴 Read this BEFORE any class above is called a CAUSE. BOUNDARY and INTACT are the\n"
        "     MODAL states of a molecule that round-trips perfectly, so neither discriminates.\n"
        "     The charter's 'BOUNDARY = the attachment call is uncertain' is a hypothesis this\n"
        "     control refutes as an explanation for failure."
    )
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument(
        "--characterise",
        action="store_true",
        help="v0.4.16 Lane 2: also classify WHICH key component disagrees, per bucket/class",
    )
    args = ap.parse_args()

    results_dir = args.results_dir.resolve()
    buckets = load_buckets(results_dir)
    reports_dir = results_dir / "individual_reports"
    if not reports_dir.is_dir():
        sys.exit(f"FATAL: no individual_reports/ under {results_dir}")

    # bucket -> class -> [molecules]
    table: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    missing_report = []
    n_seen = 0

    for molecule, (bucket, _subclass) in sorted(buckets.items()):
        path = reports_dir / f"{molecule}.json"
        if not path.exists():
            missing_report.append(molecule)
            continue
        try:
            report = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            missing_report.append(molecule)
            continue
        n_seen += 1
        table[bucket][classify(report)].append(molecule)

    # ---- report -------------------------------------------------------------
    all_buckets = sorted(table)
    print(f"# attachment-class audit -- {results_dir.name}")
    print(f"\nmolecules in bucket report: {len(buckets)}")
    print(f"individual reports read:    {n_seen}")
    print(f"reports missing/unreadable: {len(missing_report)}")
    if missing_report:
        print(f"  (first 10: {missing_report[:10]})")

    header = f"\n| {'bucket':<18} | " + " | ".join(f"{c:>12}" for c in CLASSES) + " |   total |"
    print(header)
    print("|" + "-" * 20 + "|" + "|".join("-" * 14 for _ in CLASSES) + "|---------|")
    totals: Counter = Counter()
    for bucket in all_buckets:
        row = table[bucket]
        n = sum(len(row[c]) for c in CLASSES)
        cells = " | ".join(f"{len(row[c]):>12}" for c in CLASSES)
        print(f"| {bucket:<18} | {cells} | {n:>7} |")
        for c in CLASSES:
            totals[c] += len(row[c])
    cells = " | ".join(f"{totals[c]:>12}" for c in CLASSES)
    print(f"| {'ALL':<18} | {cells} | {sum(totals.values()):>7} |")

    # ---- the two named classes, with their control --------------------------
    failing = [b for b in all_buckets if b != "byte_exact"]
    gavsed = sum(len(table[b]["DETACHED"]) for b in failing)
    medzur = sum(len(table[b]["INTACT"]) for b in failing)
    boundary = sum(len(table[b]["BOUNDARY"]) for b in failing)
    no_struct = sum(len(table[b]["NO_STRUCTURE"]) for b in failing)
    unknown = sum(len(table[b]["UNKNOWN"]) for b in failing)
    n_failing = gavsed + medzur + boundary + no_struct + unknown

    ctrl = table.get("byte_exact", {})
    n_ctrl = sum(len(ctrl.get(c, [])) for c in CLASSES)
    ctrl_detached = len(ctrl.get("DETACHED", []))

    print(f"\n## The two classes, over the {n_failing} NON-byte_exact molecules\n")
    print(f"  GAVSED shape (DETACHED, returned anyway) : {gavsed:>5}")
    print(f"  MEDZUR shape (INTACT, indep disagrees)   : {medzur:>5}")
    print(f"  BOUNDARY (attachment call is uncertain)  : {boundary:>5}")
    print(f"  NO_STRUCTURE (nothing to judge)          : {no_struct:>5}")
    print(f"  UNKNOWN (no coordination block)          : {unknown:>5}")

    # ---- the cut that matters: key_equal is NOT a disagreement ---------------
    # `key_equal` molecules round-trip to the SAME comparison key and differ only in
    # presentation -- benign canonicalization, the encoder-side block that the donor fold and
    # v0.4.14 address. Counting them as "independent re-perception disagrees" inflates the
    # MEDZUR class with molecules whose geometry is fine and whose string is merely renumbered.
    # The MEDZUR shape as v0.4.7 defined it is a GENUINE round-trip failure with attachment
    # intact, so the honest denominator excludes key_equal.
    genuine = [b for b in failing if b != "key_equal"]
    g_detached = sum(len(table[b]["DETACHED"]) for b in genuine)
    g_intact = sum(len(table[b]["INTACT"]) for b in genuine)
    g_boundary = sum(len(table[b]["BOUNDARY"]) for b in genuine)
    g_nostruct = sum(len(table[b]["NO_STRUCTURE"]) for b in genuine)
    g_total = g_detached + g_intact + g_boundary + g_nostruct
    print(f"\n## Excluding key_equal (benign canonicalization) -- {g_total} genuine failures\n")
    print(f"  GAVSED shape (DETACHED)  : {g_detached:>5}")
    print(f"  MEDZUR shape (INTACT)    : {g_intact:>5}")
    print(f"  BOUNDARY                 : {g_boundary:>5}")
    print(f"  NO_STRUCTURE             : {g_nostruct:>5}")

    print("\n## CONTROL -- the same split over byte_exact (molecules that round-trip)\n")
    if n_ctrl:
        pct_ctrl = 100.0 * ctrl_detached / n_ctrl
        pct_fail = 100.0 * gavsed / n_failing if n_failing else 0.0
        print(f"  byte_exact DETACHED : {ctrl_detached:>5} / {n_ctrl} ({pct_ctrl:.2f}%)")
        print(f"  failing    DETACHED : {gavsed:>5} / {n_failing} ({pct_fail:.2f}%)")
        if pct_ctrl > 0:
            print(f"  enrichment          : {pct_fail / pct_ctrl:.1f}x")
        print(
            "\n  Read this BEFORE the class sizes above. If the two percentages are close, "
            "\n  detachment does not discriminate and the GAVSED number is not evidence."
        )
    else:
        print("  🔴 NO byte_exact CONTROL AVAILABLE -- the class sizes above are unfalsifiable.")

    characterisation = None
    if args.characterise:
        characterisation = _characterise(results_dir, reports_dir, table)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "results_dir": str(results_dir),
            "n_bucket_report": len(buckets),
            "n_reports_read": n_seen,
            "n_missing": len(missing_report),
            "missing": missing_report,
            "table": {b: {c: table[b][c] for c in CLASSES} for b in all_buckets},
            "summary": {
                "gavsed_detached": gavsed,
                "medzur_intact": medzur,
                "boundary": boundary,
                "no_structure": no_struct,
                "unknown": unknown,
                "n_failing": n_failing,
                "control_byte_exact_n": n_ctrl,
                "control_byte_exact_detached": ctrl_detached,
                "genuine_failures": {
                    "n": g_total,
                    "gavsed_detached": g_detached,
                    "medzur_intact": g_intact,
                    "boundary": g_boundary,
                    "no_structure": g_nostruct,
                },
            },
        }
        if characterisation is not None:
            payload["characterisation"] = characterisation
        args.out.write_text(json.dumps(payload, indent=2))
        print(f"\nwrote {args.out}")

    print(f"\n#DONE {n_seen}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
