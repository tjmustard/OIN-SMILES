"""Axial / atropisomer descriptor for the OIN encoder (Y2 P2 fix candidate).

The encoder is totally blind to biaryl atropisomerism: R-BINAP and S-BINAP encode to
byte-identical OIN strings
(``docs/agentic-notes/injectivity/INJECTIVITY_Y1_P2_AXIAL.md``). The configuration is,
however, recoverable straight from the 3D geometry as the **signed biaryl dihedral**
(``docs/agentic-notes/injectivity/INJECTIVITY_Y2_FEASIBILITY.md``). RDKit does NOT
perceive it from pure 3D -- per
the RDKit Book an atropisomer bond is only marked when a neighbour bond is *wedged*, and
the configuration is only then read from the coordinates -- so we detect the axis and sign
it ourselves.

This module is the single source of truth for that detection, shared by the encoder's
**opt-in** emit (``OIN_EMIT_AXIAL``) and the Y2 configurational oracle
(``tools/injectivity/config_oracle.py``).

The token is a **canonical** chirality descriptor, not merely a distinguisher: it depends
only on the molecular graph plus the handedness of the geometry, so it is invariant under
input atom renumbering and under any proper rotation, and flips only under reflection
(guards in ``tests/unit/test_axial_emit.py``). That is what lets the generator compare a
freshly embedded conformer's token against the one carried by the requested OIN.

Emitting the token trades a silent round-trip FALSE POSITIVE (a lossy OIN that passes) for
a generator-caused FALSE NEGATIVE, so the emit stays gated OFF by default until the
generator reliably reproduces the axis (see ``_select_by_geometry``'s axial-aware pass).
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from ..core.constants import TRANSITION_METALS_NUM

_logger = logging.getLogger(__name__)

#: the emitted suffix, e.g. ``" |ax:-|"``. Kept here so the encoder, the comparison key and
#: the generator all agree on one spelling.
AXIAL_TOKEN_RE = re.compile(r"\|ax:([+\-]*)\|")

#: an atropisomer sits twisted well away from co-planarity; a planar (conjugation-locked or
#: freely-rotating) biaryl reads |dihedral| near 0 or 180. Rejects planar inter-ring bonds
#: (e.g. the chelated phenyl-pyridyl of ppy at ~0 degrees).
_ATROP_TWIST_LO = 20.0
_ATROP_TWIST_HI = 160.0


def _skeleton(mol: Chem.Mol) -> Chem.Mol:
    """A bond-order-, charge- and aromaticity-erased copy: pure element connectivity.

    **Why the ranks must not come from the molecule as perceived.** The descriptor is only
    useful if the token an encoder computes from the deposited geometry can be compared
    against the token a *generated* conformer yields. Those two mols are built by different
    routes -- ``xyz2mol`` bond-order perception from distances on one side,
    ``build_contract_mol``'s per-fragment transfer on the other -- and for a metalloporphyrin
    they disagree about the macrocycle: the encoder reads an aromatic pyrrolide core on
    Zn(II), the generated mol a neutral localized tautomer with dative bonds pointing the
    other way. Canonical ranks computed on either mol therefore differ, which would pick a
    *different* reference ortho neighbour on each side -- and the reference neighbour sets the
    dihedral SIGN. Comparing two tokens defined against different references is not a
    comparison at all: it can report a match for the mirror-image structure.

    Erasing bond orders, charges and aromatic flags leaves the heavy-atom + explicit-H
    connectivity, which the two routes do agree on. Ranks over that skeleton are therefore
    identical for both perceptions of one structure, while remaining graph-derived -- so the
    token stays invariant under input atom renumbering and under proper rotation, and still
    flips only under reflection. ``tests/unit/test_axial_emit.py`` pins this as an explicit
    perception-invariance guard.

    **Metal bonds are dropped, not down-graded.** Turning a dative M-donor bond into a plain
    single bond closes every chelate ring, and a bidentate biaryl (BINAP) then has its own
    axis bond sitting *inside* the P-M-P ring -- ``IsInRing()`` becomes true and the axis
    disappears. Dropping the bonds also sidesteps a second inconsistency: ``DATIVE``
    direction is begin-to-end, and the two routes write it opposite ways round (the encoder
    N->M, ``build_contract_mol`` M->N), so a rank that saw the bond at all would not be
    route-independent anyway.
    """
    probe = Chem.RWMol(mol)
    metal_bonds = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        for b in probe.GetBonds()
        if b.GetBeginAtom().GetAtomicNum() in TRANSITION_METALS_NUM
        or b.GetEndAtom().GetAtomicNum() in TRANSITION_METALS_NUM
    ]
    for i, j in metal_bonds:
        probe.RemoveBond(i, j)
    for b in probe.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in probe.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
    out = probe.GetMol()
    out.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(out)
    return out


def _ranks(mol: Chem.Mol, *, break_ties: bool = False) -> list[int]:
    return list(Chem.CanonicalRankAtoms(mol, breakTies=break_ties))


def _axis_cut_ranks(skel: Chem.Mol, a1: int, a2: int) -> list[int]:
    """Skeleton symmetry ranks with the axis bond ``a1-a2`` cut.

    Stereogenicity and the reference-neighbour choice are properties of **one half** of the
    molecule -- an end is non-stereogenic when rotating that half 180° about the axis
    reproduces the structure -- so they are judged with the axis bond cut. Two reasons:

    * **It removes a silent coin toss.** The reference ortho neighbour is ``max`` by rank, so
      when the two ortho neighbours *tie* the winner is whichever the neighbour list happens
      to yield first. The two candidates sit ~180° apart, so that choice sets the SIGN. Cut
      ranks are strictly finer than intact ranks (cutting can only distinguish atoms the
      intact graph merged), which resolves ties that carry real asymmetry.
    * **It scopes the question correctly.** Whether *this* end has a local C2 must not depend
      on what the far half looks like.

    It is a refinement, not a repair: a graph automorphism is free to fix a pendant group
    pointwise while permuting the ring it hangs off, so a cut rank can still tie where a
    *geometric* rotation would disturb something remote. A 5,15-diarylporphyrin is exactly
    that case -- the two pyrrole alpha carbons flanking a meso carbon tie either way, so the
    meso-aryl axes are reported non-stereogenic. That verdict happens to be right (both the
    syn and anti configurations are achiral, see
    ``tests/unit/test_axial_emit.py::TestPorphyrinMesoAxesAreNotPerAxisStereogenic``), but it
    is right for a weaker reason than a full automorphism analysis would give. Deciding
    stereogenicity for axes *coupled* through a symmetric core needs that analysis; until it
    exists the gate stays conservative and simply does not emit. See
    ``docs/KNOWN_LIMITATIONS.md``.
    """
    cut = Chem.RWMol(skel)
    cut.RemoveBond(a1, a2)
    out = cut.GetMol()
    out.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(out)
    return _ranks(out)


def _is_trigonal_ring_atom(a: Chem.Atom) -> bool:
    """Is ``a`` a ring atom that can carry an inter-ring axis, judged from connectivity only?

    Replaces an ``atom.GetIsAromatic()`` test that was **not** invariant across the two
    perception routes (see :func:`_skeleton`): a porphyrin meso carbon reads aromatic when
    the encoder perceives the deposited crystal geometry and *aliphatic* on the generated
    side, so the meso-aryl axes of a tetra/di-arylporphyrin were detectable on the input and
    invisible on every generated conformer -- the real cause of the multi-axis round-trip
    failure (the generator does hold both twists; nothing could see them).

    A ring atom bearing exactly **three heavy neighbours and no hydrogen** is trigonal:
    with four valences to spend and only three sigma bonds it must hold a pi bond, so it is
    sp2 whatever bond-order model perceived it. That is a strict *superset* of the previous
    aromatic test -- an aromatic axis end always has two ring neighbours plus the axis
    partner and no H, and an aromatic atom already carrying three ring bonds (a fusion atom)
    has no valence left to bear an axis partner at all -- so no axis that used to be found is
    lost. Metals are excluded so a dative M-donor bond is never mistaken for an axis.

    Hydrogens are counted with ``includeNeighbors=True`` so the test reads the same whether
    the caller's mol carries explicit H atoms (both pipeline routes do) or implicit counts.
    """
    if a.GetAtomicNum() in TRANSITION_METALS_NUM or not a.IsInRing():
        return False
    if a.GetTotalNumHs(includeNeighbors=True) != 0:
        return False
    heavy = sum(1 for n in a.GetNeighbors() if n.GetAtomicNum() > 1)
    return heavy == 3


def _ring_neighbors(a: Chem.Atom, across: int) -> list[Chem.Atom]:
    """Ring neighbours of ``a`` other than the axis partner ``across`` (metals excluded).

    The dihedral reference is picked from this set. For an aromatic axis end this is exactly
    the set the previous ``GetIsAromatic()`` filter produced -- an axis end is never a ring
    fusion atom (no spare valence), so its two ring neighbours lie in the one ring and share
    its aromaticity -- but it survives a perception route that reads the ring as localized.
    """
    return [
        n
        for n in a.GetNeighbors()
        if n.GetIdx() != across and n.IsInRing() and n.GetAtomicNum() not in TRANSITION_METALS_NUM
    ]


@dataclass(frozen=True)
class AxialAxis:
    """A biaryl single-bond axis and its signed dihedral (atropisomer configuration)."""

    a1: int
    a2: int
    sign: int  # +1 / -1: the recoverable configuration
    dihedral_deg: float
    hindered: bool  # twist + ortho-wall heuristic: is this a real atropisomer candidate?
    stereogenic: bool = True  # both ring ends symmetry-asymmetric -> the sign is well-defined

    @property
    def emits(self) -> bool:
        """Whether this axis contributes to the canonical token."""
        return self.hindered and self.stereogenic


def _is_atropisomer_candidate(a1: Chem.Atom, a2: Chem.Atom, dih: float) -> bool:
    """Cheap atropisomer-candidate heuristic (axis SELECTION, not the sign).

    The signed dihedral is the CONFIGURATION side; this decides whether an axis is a real
    atropisomer candidate. Two conditions:

    * **twisted** -- |dihedral| in a hindered range, not near-planar; and
    * **ortho-walled** -- each ring end carries a heavy, non-aromatic, *non-metal* exocyclic
      ortho substituent (the steric wall). Metals are excluded so a chelated biaryl (ppy) is
      not mistaken for a hindered one on the strength of the M-C bond.

    **This is the one part of the descriptor that still reads an aromatic flag**, and
    deliberately so. Relaxing it is not safe in either direction: dropping the aromatic term
    makes every twisted biphenyl "walled" (the ortho ring carbon itself becomes a wall), which
    would emit conformational, non-configurational signs at scale; while tightening it to
    "exocyclic substituent" alone would *stop* a meso-arylporphyrin qualifying, since its
    porphyrin end is walled by the fused pyrrole rings rather than by a substituent.

    The residual asymmetry it leaves is benign in the direction that matters. A generated mol
    reads *fewer* atoms as aromatic, so it is *more* likely to call an axis hindered, never
    less: the generator can therefore never quietly drop an axis the encoder asked for. It
    could in principle report an extra one, which lengthens the token, fails the match, and
    falls through to the unfiltered pool -- loud and non-regressive, never a silent wrong
    answer. See ``docs/KNOWN_LIMITATIONS.md``.
    """
    if not (_ATROP_TWIST_LO <= abs(dih) <= _ATROP_TWIST_HI):
        return False

    def ortho_walled(a: Chem.Atom, across: int) -> bool:
        for nbr in a.GetNeighbors():
            if nbr.GetIdx() == across:
                continue
            for oo in nbr.GetNeighbors():
                if oo.GetIdx() in (a.GetIdx(), across):
                    continue
                if (
                    oo.GetAtomicNum() > 1
                    and not oo.GetIsAromatic()
                    and oo.GetAtomicNum() not in TRANSITION_METALS_NUM
                ):
                    return True
        return False

    return ortho_walled(a1, a2.GetIdx()) and ortho_walled(a2, a1.GetIdx())


def detect_axial_axes(mol: Chem.Mol) -> list[AxialAxis]:
    """Every inter-ring conjugated single bond, with its signed dihedral configuration.

    RDKit does not nominate these from 3D (see module docstring), so we enumerate them
    ourselves. Returns ``[]`` when the mol has no conformer (nothing to sign).

    An axis end must be a **trigonal ring atom** (:func:`_is_trigonal_ring_atom`) rather than
    a *flagged-aromatic* one, because aromatic flags are not stable across the encoder's and
    the generator's perception routes.

    **Canonicality.** The reference neighbour on each ring end is the ortho ring neighbour
    with the highest *symmetry* rank over the connectivity skeleton
    (:func:`_skeleton_ranks`). Those ranks depend only on the molecular graph, so the choice
    -- and hence the dihedral sign -- is invariant under input atom renumbering, under any
    proper rotation of the coordinates, and under which bond-order model perceived the
    molecule; it flips only under reflection. Tie-broken ranks must NOT be used here: RDKit
    breaks ties arbitrarily between symmetry-equivalent atoms, which would pick a reference
    ~180° away on one end of a symmetric biaryl and silently flip the sign.

    An end whose two ortho neighbours are symmetry-*equivalent* has a local C2 through the
    axis: rotating that ring 180° reproduces the molecule, so the axis is **not
    stereogenic** and its sign is meaningless. Such axes are returned with
    ``stereogenic=False`` (kept for diagnostics) and never reach the token.
    """
    if mol.GetNumConformers() == 0:
        return []
    conf = mol.GetConformer()
    skel = _skeleton(mol)
    out: list[AxialAxis] = []
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE or b.IsInRing():
            continue
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if not (_is_trigonal_ring_atom(a1) and _is_trigonal_ring_atom(a2)):
            continue
        n1 = _ring_neighbors(a1, a2.GetIdx())
        n2 = _ring_neighbors(a2, a1.GetIdx())
        if not n1 or not n2:
            continue
        sym = _axis_cut_ranks(skel, a1.GetIdx(), a2.GetIdx())
        # stereogenic only when neither ring end is symmetry-degenerate about the axis
        stereogenic = len({sym[n.GetIdx()] for n in n1}) == len(n1) and len(
            {sym[n.GetIdx()] for n in n2}
        ) == len(n2)
        r1 = max(n1, key=lambda n: sym[n.GetIdx()])
        r2 = max(n2, key=lambda n: sym[n.GetIdx()])
        dih = rdMolTransforms.GetDihedralDeg(
            conf, r1.GetIdx(), a1.GetIdx(), a2.GetIdx(), r2.GetIdx()
        )
        out.append(
            AxialAxis(
                a1=a1.GetIdx(),
                a2=a2.GetIdx(),
                sign=1 if dih > 0 else -1,
                dihedral_deg=round(float(dih), 2),
                hindered=_is_atropisomer_candidate(a1, a2, float(dih)),
                stereogenic=stereogenic,
            )
        )
    return out


def axial_token(mol: Chem.Mol) -> str:
    """Opt-in atropisomer token for the raw OIN string; ``""`` when no hindered axis.

    The signs of every hindered, stereogenic axis, ordered by the skeleton *symmetry* ranks
    of the axis atoms so the token depends only on the molecular graph plus the geometry --
    not on input atom order, and not on which bond-order model perceived the molecule
    (:func:`_skeleton_ranks`). One hindered axis (BINAP) -> ``"-"`` (R) / ``"+"`` (S). No
    qualifying axis -> ``""`` (so ordinary molecules are completely unaffected by the flag).

    **The sort must never depend on the sign.** Symmetry-equivalent axes tie on symmetry
    rank, and an earlier version broke that tie with ``ax.sign`` to keep the string
    order-independent. That silently destroyed the chirality: sorting by sign forces the
    signs into ascending order, so a molecule with two equivalent axes carrying ``-+``
    renders identically to its mirror image carrying ``+-``. The corpus sign-convention
    audit caught exactly this on three multi-axis structures (OJELAQ, YESKOZ, EBUHAN), all
    of which the independent geometric oracle calls chiral. Ties are broken with the
    tie-broken canonical rank instead: also graph-derived, hence renumbering-invariant, but
    it keeps each sign attached to its own axis so reflection flips the token.
    """
    axes = [ax for ax in detect_axial_axes(mol) if ax.emits]
    if not axes:
        return ""
    skel = _skeleton(mol)
    sym, tie = _ranks(skel), _ranks(skel, break_ties=True)
    axes.sort(
        key=lambda ax: (
            tuple(sorted((sym[ax.a1], sym[ax.a2]))),
            tuple(sorted((tie[ax.a1], tie[ax.a2]))),
        )
    )
    return "".join("+" if ax.sign > 0 else "-" for ax in axes)


def parse_axial_token(oin_string: str | None) -> str | None:
    """The axial token carried by an OIN string, or ``None`` when it carries none.

    Absent token -> ``None`` means "no axial constraint requested", which is what makes the
    generator's axial-aware selection self-gating: an OIN encoded without
    ``OIN_EMIT_AXIAL`` simply has nothing to honour, so generation is unchanged.
    """
    if not oin_string:
        return None
    m = AXIAL_TOKEN_RE.search(oin_string)
    return m.group(1) if m else None


def mol_axial_token(mol: Chem.Mol) -> str | None:
    """``axial_token`` for a generated conformer; ``None`` if perception fails.

    Safe to call on generator output: any RDKit failure degrades to ``None`` (treated as
    "cannot tell") rather than raising into the generation path.

    Raw MetalloGen pool conformers are typically UNsanitized, so ring/aromaticity info and
    implicit valences are missing and ``CanonicalRankAtoms`` raises. Retry on a sanitized
    copy, then on a partial sanitize that skips kekulization (the same degradation
    ``oin/compare.py`` uses for metal-bearing fragments).
    """
    last_error: Exception | None = None
    for prepare in (
        lambda m: m,
        lambda m: _sanitized_copy(m),
        lambda m: _sanitized_copy(m, skip_kekulize=True),
    ):
        try:
            probe = prepare(mol)
            if probe is not None:
                return axial_token(probe)
        except Exception as e:  # noqa: PERF203 - each strategy may fail differently
            last_error = e
            continue
    # Returning None here is indistinguishable from "no axis" to a careless caller, and that
    # is precisely how the generator's axial pass was once defeated in silence. Leave a trace
    # so the failure is at least discoverable; callers that care must treat None as
    # "could not tell", never as "no axial token".
    _logger.debug("axial perception failed on all strategies: %r", last_error)
    return None


def _sanitized_copy(mol: Chem.Mol, *, skip_kekulize: bool = False) -> Chem.Mol | None:
    probe = Chem.Mol(mol)
    ops = Chem.SanitizeFlags.SANITIZE_ALL
    if skip_kekulize:
        ops = ops ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
    Chem.SanitizeMol(probe, sanitizeOps=ops)
    return probe
