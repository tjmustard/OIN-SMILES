"""Canonical comparison of OIN strings for round-trip equivalence.

These helpers decide whether two OIN strings describe the *same* coordination
compound, collapsing chemically-meaningless notation drift (implicit-vs-explicit
donor H, Kekulé-vs-aromatic rings, NHC carbene bare-``C`` vs ``[CH2]``, which
symmetric carboxylate O carries the binding slot, fragment order, and the
eta-winding labeling ambiguity for symmetric rings) while still distinguishing a
genuinely different structure, metal, geometry, or eta winding multiset.

The module depends only on ``re``, RDKit, and ``OINInlineHandler`` so it can be
imported without pulling in the heavy 3D-generation stack (scine_molassembler,
MACE, ...). The round-trip harness and the interactive verifier both use it.
"""

import re

from rdkit import Chem

from .inline import OINInlineHandler

_METAL_STEREO_RE = re.compile(r"\[([A-Z][a-z]?)@[A-Z0-9]+_([A-Z]{3})\]")


def normalize_oin_for_comparison(oin_string: str) -> str:
    """Normalize an OIN string for round-trip comparison.

    1. Strip atom-ordering-dependent @SP/@OH/@TB stereo descriptors from the
       metal fragment — the slot assignments already encode the isomer geometry;
       the @XY## label depends on XYZ atom ordering and is not reproducible.
    2. Remove empty fragments (consecutive/trailing dots) caused by ligands that
       are present in the XYZ but uncoordinated in the OIN (e.g. H2 in FeH2(CO)4).
    3. Normalize water notation: [OH2] and O are chemically equivalent as bound
       water ligands. The XYZ→OIN pipeline may write O while generated structures
       re-analyzed after H addition write [OH2].
    4. Winding direction markers ({n>} / {n<}) are KEPT and compared verbatim.
       (Historically they were stripped, on the assumption that an eta ligand's
       ring rotation/face could not be reproduced from the OIN alone.) The
       encoder now emits a winding marker per eta ring -- per haptic slot, using
       each ring's actual metal->centroid axis -- so the OIN string losslessly
       encodes eta stereochemistry (an ansa-metallocene's rac/meso, a ring's
       coordinated face). Comparing winding is exactly what lets the round trip
       catch a generated wrong-face / wrong-diastereomer eta ligand, which the
       coordination-sphere RMSD (eta rings reduced to a centroid) cannot see.
    5. Canonicalize slot numbering: for OCT and other symmetric geometries where
       different rotations yield equivalent but numerically different slot assignments,
       renumber slots in order of first appearance. This makes equivalently-rotated
       structures map to the same OIN after normalization (e.g., OCT with N atoms
       at slots {3,5} vs {5,3} both normalize to {0,1} for the N atoms).
    """
    s = _METAL_STEREO_RE.sub(r"[\1_\2]", oin_string)
    # Normalize [OH2] → O (bound water notation equivalence)
    s = s.replace("[OH2]", "O")
    # Winding markers ({n>} / {n<}) are intentionally NOT stripped -- they carry
    # eta-ligand stereochemistry that the round trip must verify (see docstring).
    # Collapse multiple consecutive dots and strip trailing dots
    while ".." in s:
        s = s.replace("..", ".")
    s = s.rstrip(".")

    # Canonicalize slot numbering: renumber slots in order of first appearance
    slot_map = {}
    next_slot = 0

    def replace_slot(match):
        nonlocal next_slot
        old_slot = int(match.group(1))
        winding = match.group(2) or ""  # preserve the {n>}/{n<} marker, if any
        if old_slot not in slot_map:
            slot_map[old_slot] = next_slot
            next_slot += 1
        return "{" + str(slot_map[old_slot]) + winding + "}"

    s = re.sub(r"\{(\d+)([><^]?)\}", replace_slot, s)
    return s


def winding_canonical_key(normalized_oin: str):
    """Canonical comparison key that treats eta-ring winding as a MULTISET.

    Returns ``(winding_stripped_string, sorted_winding_multiset)``.

    Two OIN strings that describe the same molecule but differ only in which of
    two EQUIVALENT eta rings is labeled the lower slot must compare equal. An
    achiral *meso* ansa-metallocene can be written ``{0<}{1>}`` or ``{0>}{1<}``
    -- the two rings are interchangeable, so both are the same structure; only
    which one the encoder happened to call slot 0 differs (a canonicalization
    ambiguity for symmetric rings). Comparing the winding as an order-independent
    multiset makes those equal, while still catching a real error:
      * a diastereomer flip: ``['<','>']`` (meso) vs ``['>','>']`` (rac)
      * an enantiomer flip:  ``['>','>']``       vs ``['<','<']``
    both change the multiset and still fail. The winding-stripped remainder must
    still match exactly, so every non-winding difference is caught as before.
    """
    windings = sorted(re.findall(r"\{\d+([<>])\}", normalized_oin))
    stripped = re.sub(r"\{(\d+)[<>]\}", r"{\1}", normalized_oin)
    return stripped, windings


def _canonical_fragment_smiles(smiles: str) -> str:
    """RDKit-canonical SMILES for one slot-stripped ligand fragment.

    Sanitizing normalizes the notation differences that make two chemically
    identical fragments serialize differently: explicit-vs-implicit donor H
    (``[NH]``/``N``, ``[OH]``/``O``), Kekulé-vs-aromatic rings, and a bare ring
    carbon vs ``[CH2]`` (the NHC-carbene-carbon case). Falls back to a stable
    ``RAW:`` token if the fragment cannot be parsed, so exotic ligands (e.g.
    borane clusters) still compare by string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return "RAW:" + smiles
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            return "RAW:" + smiles
    return Chem.MolToSmiles(mol)


def canonical_roundtrip_key(oin_string: str):
    """Structure-level canonical key for round-trip equivalence.

    Two OIN strings that describe the same coordination compound produce the
    same key: same metal + geometry, the same *multiset* of ligand fragments
    (each canonicalized through RDKit, so implicit-H / Kekulé / carbene notation
    and fragment ordering do not matter), and the same eta *winding multiset*.

    It is strictly stronger than string equality — it collapses chemically
    meaningless notation drift (the dominant round-trip failure mode: explicit
    vs implicit donor H, carbene bare-``C`` vs ``[CH2]``, which symmetric
    carboxylate O carries the slot, fragment order) — while still distinguishing
    a genuinely different structure (allyl double-bond loss, dearomatization,
    P-stereocenter loss), a different metal/geometry (e.g. SPY vs TBP), or a
    different eta winding multiset (a rac/meso or face flip).

    Trade-off: binding-atom identity *within* a fragment is not part of the key,
    so a pure linkage isomer (same fragment bound through a different atom with
    no other geometric change) would compare equal here. The RMSD and atom-count
    checks in the round trip guard that case separately.
    """
    normalized = normalize_oin_for_comparison(oin_string.strip())

    metal_match = OINInlineHandler.METAL_REGEX.search(normalized)
    metal = (metal_match.group(1), metal_match.group(2)) if metal_match else ("", "")
    body = OINInlineHandler.METAL_REGEX.sub("", normalized)

    windings = tuple(
        sorted(">" if w == "^" else w for w in re.findall(r"\{\d+([<>^])\}", normalized))
    )

    frag_keys = []
    for frag in body.split("."):
        clean = OINInlineHandler.SLOT_REGEX.sub("", frag).strip()
        if not clean:
            continue
        frag_keys.append(_canonical_fragment_smiles(clean))

    return (metal, tuple(sorted(frag_keys)), windings)
