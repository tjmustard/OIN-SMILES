"""Axial / atropisomer descriptor for the OIN encoder (Y2 P2 fix candidate).

The encoder is totally blind to biaryl atropisomerism: R-BINAP and S-BINAP encode to
byte-identical OIN strings (``docs/INJECTIVITY_Y1_P2_AXIAL.md``). The configuration is,
however, recoverable straight from the 3D geometry as the **signed biaryl dihedral**
(``docs/INJECTIVITY_Y2_FEASIBILITY.md``). RDKit does NOT perceive it from pure 3D -- per
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
    """Every inter-ring aromatic single bond, with its signed dihedral configuration.

    RDKit does not nominate these from 3D (see module docstring), so we enumerate them
    ourselves. Returns ``[]`` when the mol has no conformer (nothing to sign).

    **Canonicality.** The reference neighbour on each ring end is the ortho neighbour with
    the highest *symmetry* rank (``breakTies=False``). Symmetry ranks depend only on the
    molecular graph, so the choice -- and hence the dihedral sign -- is invariant under
    input atom renumbering and under any proper rotation of the coordinates, and flips
    under reflection. Tie-broken ranks must NOT be used here: RDKit breaks ties
    arbitrarily between symmetry-equivalent atoms, which would pick a reference ~180° away
    on one end of a symmetric biaryl and silently flip the sign.

    An end whose two ortho neighbours are symmetry-*equivalent* has a local C2 through the
    axis: rotating that ring 180° reproduces the molecule, so the axis is **not
    stereogenic** and its sign is meaningless. Such axes are returned with
    ``stereogenic=False`` (kept for diagnostics) and never reach the token.
    """
    if mol.GetNumConformers() == 0:
        return []
    conf = mol.GetConformer()
    sym = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    out: list[AxialAxis] = []
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE or b.IsInRing():
            continue
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if not (a1.GetIsAromatic() and a2.GetIsAromatic()):
            continue
        n1 = [n for n in a1.GetNeighbors() if n.GetIdx() != a2.GetIdx() and n.GetIsAromatic()]
        n2 = [n for n in a2.GetNeighbors() if n.GetIdx() != a1.GetIdx() and n.GetIsAromatic()]
        if not n1 or not n2:
            continue
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

    The signs of every hindered, stereogenic axis, ordered by the *symmetry* ranks of the
    axis atoms so the token depends only on the molecular graph plus the geometry -- not on
    input atom order. One hindered axis (BINAP) -> ``"-"`` (R) / ``"+"`` (S). No qualifying
    axis -> ``""`` (so ordinary molecules are completely unaffected by the flag).

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
    sym = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    tie = list(Chem.CanonicalRankAtoms(mol, breakTies=True))
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
