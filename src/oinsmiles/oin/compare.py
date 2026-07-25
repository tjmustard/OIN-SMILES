"""Canonical comparison of OIN strings for round-trip equivalence.

These helpers decide whether two OIN strings describe the *same* coordination
compound, collapsing chemically-meaningless notation drift (implicit-vs-explicit
donor H, Kekulé-vs-aromatic rings, NHC carbene bare-``C`` vs ``[CH2]``, which
symmetric carboxylate O carries the binding slot, fragment order, and the benign
slot-relabeling drift between true conformers of a symmetric coordination sphere)
while still distinguishing a genuinely different structure, metal, geometry, eta
winding, or *positional* isomer (fac vs mer, cis vs trans).

The module depends only on ``re``, RDKit, ``OINInlineHandler``, and the numpy-only
``canonical_slots`` so it can be imported without pulling in the heavy 3D-generation
stack (MetalloGen, MACE, ...). The round-trip harness and the interactive verifier
both use it.
"""

import re

from rdkit import Chem

from .canonical_slots import (
    GEOMETRY_VERTICES,
    VERTEX_SENTINEL,
    derive_rotation_group,
    geometry_rotation_group,
    geometry_vertex_count,
    lexmin_vertex_signature,
)
from .inline import OINInlineHandler

# Sanitize everything EXCEPT kekulization. A slot-stripped chelate fragment can be a
# neutral all-carbon eta-Cp/Cp* ring or a bare-``n`` 5-membered azole/pyridyl donor
# ring whose donor N just lost its metal bond -- both raise ``KekulizeException`` on a
# full sanitize, yet RDKit can ring-perceive and canonicalize them fine without
# kekulizing. Skipping only that step lets the chelate-lock E/Z clearing run on
# exactly the fragments that carry the ring-locked slash (see _parse_fragment).
_NO_KEKULIZE = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE


def _parse_fragment(smiles: str):
    """Parse one slot-stripped ligand fragment to an RDKit mol, or ``None``.

    Full sanitize first (the fast, common path). On failure -- almost always a
    ``KekulizeException`` on an aromatic ring that lost its metal-donor context --
    retry with a partial sanitize that skips kekulization. Returns ``None`` only for
    a genuinely unparseable fragment (borane cluster, over-valent ``C#O``), so the
    ``RAW:`` fallback still guards those.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol, sanitizeOps=_NO_KEKULIZE)
    except Exception:
        return None
    return mol


_METAL_STEREO_RE = re.compile(r"\[([A-Z][a-z]?)@[A-Z0-9]+_([A-Z]{3})\]")
#: opt-in axial / atropisomer suffix ` |ax:+-|` (OIN_EMIT_AXIAL); folded out of the key.
_AXIAL_TOKEN_RE = re.compile(r"\s*\|ax:[+\-]*\|")


def normalize_oin_for_comparison(oin_string: str) -> str:
    """Normalize an OIN string for round-trip comparison.

    1. Strip atom-ordering-dependent @SP/@OH/@TB stereo descriptors from the
       metal fragment — the slot assignments already encode the isomer geometry;
       the @XY## label depends on XYZ atom ordering and is not reproducible.
       (Keeping metal @-chirality in the key is DEFERRED — it needs a reproducible
       encoder-side metal stereo descriptor first; see
       ``spec/handoffs/v0.4.4/geometry-canonical-slot-key.md`` §4.)
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

    Slots are NOT renumbered here. The former first-appearance renumber folded
    conformer slot-relabeling drift, but it also folded fac↔mer (which differ only
    in absolute slot numbers), making the key blind to a real positional isomer.
    The conformer drift is folded downstream instead, by the symmetry-aware vertex
    signature (see ``winding_canonical_key`` / ``_polyhedron_signature``), which
    keeps fac and mer distinct. The absolute slots survive here so that signature
    has them to work with; the only other consumers of this string are the
    human-readable ``Exp:/Got:`` round-trip diagnostics.
    """
    # Fold the opt-in axial / atropisomer token (` |ax:+-|`, OIN_EMIT_AXIAL). Like the
    # metal @-stereo above, it is a distinguisher the round-trip key does not yet gate on
    # (the generator cannot reproduce the axis), so the batch key stays blind to it and the
    # harness is unaffected whether or not the emit flag is set. Y2 P2; see
    # docs/INJECTIVITY_Y2_FEASIBILITY.md.
    s = _AXIAL_TOKEN_RE.sub("", oin_string)
    s = _METAL_STEREO_RE.sub(r"[\1_\2]", s)
    # Normalize [OH2] → O (bound water notation equivalence)
    s = s.replace("[OH2]", "O")
    # Winding markers ({n>} / {n<}) are intentionally NOT stripped -- they carry
    # eta-ligand stereochemistry that the round trip must verify (see docstring).
    # Collapse multiple consecutive dots and strip trailing dots
    while ".." in s:
        s = s.replace("..", ".")
    s = s.rstrip(".")
    return s


_SLOT_RE = re.compile(r"\{\d+[><^]?\}")


def _canonical_fragment_smiles(smiles: str) -> str:
    """RDKit-canonical SMILES for one slot-stripped ligand fragment.

    Sanitizing normalizes the notation differences that make two chemically
    identical fragments serialize differently: explicit-vs-implicit donor H
    (``[NH]``/``N``, ``[OH]``/``O``), Kekulé-vs-aromatic rings, and a bare ring
    carbon vs ``[CH2]`` (the NHC-carbene-carbon case). Falls back to a stable
    ``RAW:`` token if the fragment cannot be parsed, so exotic ligands (e.g.
    borane clusters) still compare by string.
    """
    mol = _parse_fragment(smiles)
    if mol is None:
        return "RAW:" + smiles
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return "RAW:" + smiles


def _chelate_locked_fragment_key(frag: str) -> str:
    """Canonical fragment key that clears E/Z on metal-chelate-locked bonds.

    ``frag`` still carries its ``{n}`` binding-slot markers. A double bond held
    rigid by a ring that closes *through the metal* (salicylaldiminate,
    beta-diketiminate, bis-imine, an eta-bound alkene) has no free E/Z, but once
    the metal is stripped the ring opens and RDKit hands the bond a directional
    marker whose sign depends on SMILES traversal. The encoder
    (``_clear_chelate_locked_bond_stereo``) drops that marker on the input, but a
    generated structure whose donor is bonded differently can keep it, so the
    round trip spuriously fails on a bond that was never freely E/Z.

    This reconstructs the chelate rings (a dummy metal bonded to every slot atom)
    and clears E/Z on the double bonds those rings lock, on BOTH the input and
    the generated fragment, so they compare equal. A double bond in *no* metal
    ring -- a pendant, genuinely flippable alkene/imine -- is untouched and still
    distinguishes a real diastereomer. Any failure falls back to the plain
    slot-stripped canonicalization, so behaviour is unchanged for every fragment
    without a metal-locked double bond.
    """
    clean = _SLOT_RE.sub("", frag).strip()
    slots = list(_SLOT_RE.finditer(frag))
    if not clean or not slots or len(slots) > 10:
        return _canonical_fragment_smiles(clean)

    real = _parse_fragment(clean)
    if real is None:
        return _canonical_fragment_smiles(clean)

    try:
        # Build a probe: [Fe] bonded to each slot atom via a ring-closure bond,
        # so the chelate rings that close through the metal become perceivable.
        labels = [f"%9{i}" for i in range(len(slots))]
        out, last = [], 0
        for lab, m in zip(labels, slots):
            out.append(frag[last : m.start()])
            out.append(lab)
            last = m.end()
        out.append(frag[last:])
        probe_smiles = "[Fe]" + "".join(labels) + "." + "".join(out)

        probe = Chem.MolFromSmiles(probe_smiles, sanitize=False)
        if probe is None:
            return Chem.MolToSmiles(real)
        Chem.FastFindRings(probe)
        rings = Chem.GetSymmSSSR(probe)

        # The dummy metal is atom 0; ligand atoms follow in clean's order, so a
        # probe atom index i maps to real atom i-1.
        locked: set[int] = set()
        for ring in rings:
            ring_set = set(ring)
            if 0 in ring_set:
                locked |= ring_set

        for bond in probe.GetBonds():
            if bond.GetBondType() != Chem.BondType.DOUBLE:
                continue
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a == 0 or b == 0 or a not in locked or b not in locked:
                continue
            rbond = real.GetBondBetweenAtoms(a - 1, b - 1)
            if rbond is None or rbond.GetBondType() != Chem.BondType.DOUBLE:
                continue
            rbond.SetStereo(Chem.BondStereo.STEREONONE)
            for end in (rbond.GetBeginAtom(), rbond.GetEndAtom()):
                for nb in end.GetBonds():
                    if nb.GetBondDir() != Chem.BondDir.NONE:
                        nb.SetBondDir(Chem.BondDir.NONE)

        return Chem.MolToSmiles(real)
    except Exception:
        return _canonical_fragment_smiles(clean)


# --- Coordination-polyhedron symmetry, for the fac/mer-aware canonical key -----------
#
# The canonical key must fold *conformer* slot-relabeling drift (a rigid rotation of the
# idealized polyhedron relabels its vertices) while still separating genuine positional
# isomers (fac vs mer). We do that with a symmetry-aware vertex signature: color each
# occupied polyhedron vertex, then take the lexicographically-minimal image over the
# geometry's proper-rotation group. First-appearance slot renumbering (the old approach)
# folded fac/mer too, which is the blindness this replaces.
#
# The machinery itself lives in ``oin.canonical_slots`` (numpy-only, so importing it does
# not compromise this module's light import graph). It moved there in v0.4.5 because the
# ENCODER needs the same lex-min -- plus the permutation that achieves it, which is the
# canonical slot relabeling -- so the geometry table now has a single home (TD-005).
# These aliases keep the key's internals unchanged.
_GEOMETRY_VERTICES = GEOMETRY_VERTICES
_derive_rotation_group = derive_rotation_group
_geometry_rotation_group = geometry_rotation_group
_geometry_vertex_count = geometry_vertex_count
_VERTEX_SENTINEL = VERTEX_SENTINEL


def _donor_element(frag: str, marker_start: int) -> str:
    """Element symbol of the atom immediately preceding a ``{slot}`` marker.

    Handles bracket atoms (``[Br]``, ``[NH]``, ``[cH]``, ``[P@]``), two-letter organic
    ``Cl``/``Br``, and the organic subset (``c``, ``n``, ``O``, ...). Aromatic lower-case
    is folded to the element (``c`` -> ``C``), so the color is at element granularity.
    """
    j = marker_start - 1
    if j >= 0 and frag[j] == "]":
        inner = frag[frag.rfind("[", 0, j) + 1 : j]
        m = re.match(r"[A-Za-z][a-z]?", inner)
        sym = m.group(0) if m else inner
    elif j - 1 >= 0 and frag[j - 1 : j + 1] in ("Cl", "Br"):
        sym = frag[j - 1 : j + 1]
    elif j >= 0:
        sym = frag[j]
    else:
        sym = ""
    return sym[:1].upper() + sym[1:] if sym else sym


def _parse_vertex_colors(normalized: str):
    """Map each occupied polyhedron vertex to a color ``(fragment_body, element, winding)``.

    Returns ``(metal, geometry, {slot: color})``. No chelate grouping: every donor atom of
    a ligand is colored with that ligand's *whole* RDKit-canonical body, so a swap between
    two same-colored donors (a symmetric ligand's interchangeable donors, or the
    guanidinate slot drift in the CETDAI conformer fixture) is invisible -- which is exactly
    what lets true conformers collapse while fac/mer, whose donor *arrangement* differs,
    stay distinct.
    """
    metal_match = OINInlineHandler.METAL_REGEX.search(normalized)
    metal = (metal_match.group(1), metal_match.group(2)) if metal_match else ("", "")
    geo = metal_match.group(2) if metal_match else ""

    vcolor: dict[int, tuple[str, str, str]] = {}
    for frag in OINInlineHandler.METAL_REGEX.sub("", normalized).split("."):
        if not frag.strip():
            continue
        body = _chelate_locked_fragment_key(frag)
        for m in OINInlineHandler.SLOT_REGEX.finditer(frag):
            slot = int(m.group(1))
            winding = m.group(2) or ""
            if winding == "^":
                winding = ">"
            element = _donor_element(frag, m.start())
            if slot not in vcolor:
                vcolor[slot] = (body, element, winding)
            elif winding and not vcolor[slot][2]:
                # An eta ring shares one slot across several atoms; the winding marker may
                # sit on any of them -- keep it wherever it appears.
                prev = vcolor[slot]
                vcolor[slot] = (prev[0], prev[1], winding)
    return metal, geo, vcolor


def _polyhedron_signature(geo: str, vcolor: dict) -> tuple:
    """Lexicographically-minimal colored-vertex tuple over the geometry's rotation group.

    Two slot labelings related by a proper rotation of the idealized polyhedron (conformer
    drift) canonicalize to the same tuple; labelings not so related (fac vs mer) do not.
    Proper rotations preserve winding sense, so a winding char is never flipped.

    Thin wrapper over ``canonical_slots.lexmin_vertex_signature``, which also returns the
    achieving permutation -- that permutation is the canonical slot relabeling the encoder
    applies. The key only needs the signature, so it drops the permutation here and its
    behavior is unchanged.
    """
    signature, _perm = lexmin_vertex_signature(geo, vcolor)
    return signature


def _canonical_key(normalized: str):
    """Fac/mer-aware canonical key from an already-normalized OIN string.

    ``(metal, vertex_signature, fragment_body_multiset)`` -- hashable. The signature carries
    the donor arrangement over the coordination polyhedron (fac vs mer, cis vs trans) and
    the eta winding; the body multiset carries every ligand's RDKit-canonical structure
    (so implicit-H / Kekulé / carbene / E-Z / fragment-order drift collapse, while a real
    structural change or an uncoordinated-fragment difference does not).
    """
    metal, geo, vcolor = _parse_vertex_colors(normalized)
    signature = _polyhedron_signature(geo, vcolor)

    frag_keys = []
    for frag in OINInlineHandler.METAL_REGEX.sub("", normalized).split("."):
        if OINInlineHandler.SLOT_REGEX.sub("", frag).strip():
            frag_keys.append(_chelate_locked_fragment_key(frag))

    return (metal, signature, tuple(sorted(frag_keys)))


def winding_canonical_key(normalized_oin: str):
    """Fac/mer-aware canonical comparison key over an already-normalized OIN string.

    Returns ``(metal, vertex_signature, fragment_body_multiset)`` -- a hashable tuple
    (an improvement over the pre-v0.4.4 unhashable ``(str, list)``). It folds the benign
    slot-relabeling drift between true conformers of a symmetric sphere (via the
    symmetry-aware vertex signature) while distinguishing genuine positional isomers
    (fac vs mer, cis vs trans), geometry, and eta-ring winding (rac/meso, coordinated
    face). Metal ``@SPn`` chirality is intentionally NOT in the key (deferred -- the
    encoder has no reproducible metal stereo descriptor yet; see
    ``spec/handoffs/v0.4.4/geometry-canonical-slot-key.md`` §4).

    ``winding_canonical_key(normalize_oin_for_comparison(oin))`` is the round-trip
    equivalence key used by the convergence/divergence/smoke tests and the dataset harness.
    """
    return _canonical_key(normalized_oin)


def canonical_roundtrip_key(oin_string: str):
    """Structure-level, fac/mer-aware canonical key for round-trip equivalence.

    Two OIN strings that describe the same coordination compound produce the same key:
    same metal + geometry, the same symmetry-canonical vertex signature (donor arrangement
    over the coordination polyhedron, including eta winding -- so fac vs mer, cis vs trans,
    and a winding/face flip all differ), and the same *multiset* of RDKit-canonical ligand
    fragment bodies (so implicit-H / Kekulé / carbene notation and fragment order do not
    matter, and an uncoordinated fragment is still counted).

    It is strictly stronger than string equality and strictly finer than the pre-v0.4.4 key
    (which renumbered slots by first appearance and so was blind to fac/mer). Metal ``@SPn``
    chirality is deferred (see ``winding_canonical_key``). Binding-atom identity *within* a
    fragment is still not in the key, so a pure linkage isomer (same fragment bound through
    a different atom with no other geometric change) is caught by the RMSD and atom-count
    checks in the round trip, not here.
    """
    return _canonical_key(normalize_oin_for_comparison(oin_string.strip()))
