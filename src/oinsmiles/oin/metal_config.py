"""Metal-centred configuration (Y1 blind spot P1: Δ/Λ helicity) — v0.4.6 Lane 5.

WHY
===
The Y1 injectivity audit proved the encoder collapses metal-centred enantiomers: Δ- and
Λ-tris(bidentate) complexes encode byte-identically. Lane 2 measured **0 of 150** molecules
emitting a metal ``@`` tag, so this is not a matter of un-folding something the encoder already
writes — the descriptor does not exist and has to be built.

WHY NOT ``@OHn`` / ``@SPn``
==========================
OpenSMILES' non-tetrahedral tags are defined against the **atom order in the SMILES string**.
That is precisely what makes RDKit's raw ``_chiralPermutation`` non-reproducible across a
re-parse, and it is the defect Lane 8 had to fix for tetrahedral tags. A slot-sequence
descriptor is self-describing instead: it is defined against Lane 2's *canonical* slot
labelling, which is already invariant under rotation and renumbering by construction.

THE DESCRIPTOR
==============
Order the donors by canonical slot index, take the first four, and compute the sign of the
signed volume (scalar triple product) of the three edge vectors from the first:

    sign det[ p1-p0, p2-p0, p3-p0 ]

* **invariant under proper rotation** — ``det(RA) = det(R)det(A) = +det(A)`` for ``det R = +1``;
* **inverts under reflection** — an improper operation has ``det R = -1``;
* **invariant under atom renumbering** — the ordering is by canonical slot, never by atom index;
* **empty for a planar coordination sphere** — the volume degenerates to zero, so square-planar
  complexes emit nothing. That is chemically right, not a gap: a square-planar complex's
  coordination plane IS a mirror plane, so four different donors give diastereomers rather than
  enantiomers, and the correct distinctness operator there is a donor swap, not reflection.
  ``tests/unit/test_metal_stereo_fixtures.py`` records that asymmetry for ``JEGKOW``.

The token is a trailing sidecar, ``|mc:+|`` / ``|mc:-|``, following the landed ``|ax:±|``
precedent: it survives the parser (``generation/oin_parser.py`` already strips sidecars), it
leaves ``[El_GEO]`` untouched so ``METAL_REGEX`` and ``_METAL_STEREO_RE`` are unaffected, and it
cannot collide with a ligand body.

STATUS — DESCRIPTOR AND VALIDATION ONLY, NOT WIRED TO EMIT
==========================================================
This module computes and validates the descriptor. It is deliberately **not** yet called from
``xyz2mol.py``'s emit path, and no lever turns it on, because emitting it is the half that
requires the generator to reproduce what it emits — the trade recorded for every
information-ADDING lever in ``levers.py::_HELD_OFF``. Wiring is the next increment; the
three-property proof has to come first, because the Y2 wave shipped an axial descriptor that
was accidentally reflection-invariant and every single-fixture guard passed.
"""

from __future__ import annotations

import numpy as np

#: Planarity threshold on the **dimensionless** triple product
#: ``det[e1,e2,e3] / (|e1||e2||e3|)``, which lies in [-1, 1] and is therefore comparable across
#: metals and bond lengths.
#:
#: An ABSOLUTE volume threshold was tried first and is wrong in kind, not merely in value: the
#: signed volume scales as (bond length)^3, so a threshold that suits a 2.0 Å Rh-N sphere is
#: meaningless for a 2.7 Å Rh-I one, and there is no single number that works for both.
#:
#: Measured on the two Lane-5 fixtures: ZUMNEC (chiral tris-bidentate) reads -0.862, JEGKOW
#: (square planar, slight crystallographic pucker) reads +0.026. A factor of ~30 separates them,
#: so this threshold is a wide margin rather than a tuned constant.
_PLANARITY_EPS = 0.15

__all__ = ["metal_config_sign", "metal_config_token"]


def metal_config_sign(donor_positions) -> int:
    """``+1`` / ``-1`` for the handedness of a donor set, or ``0`` when planar/degenerate.

    Args:
        donor_positions: donor coordinates **already ordered by canonical slot index**. The
            caller owns that ordering — this function deliberately does no sorting of its own,
            so the invariance property is inherited from Lane 2's canonicalization rather than
            re-derived here (a second derivation is a second thing to drift).

    Fewer than four donors cannot define a chirality, and four coplanar donors define an
    achiral sphere; both return ``0``.
    """
    pts = np.asarray(donor_positions, dtype=float)
    if pts.ndim != 2 or pts.shape[0] < 4 or pts.shape[1] != 3:
        return 0
    p0, p1, p2, p3 = pts[0], pts[1], pts[2], pts[3]
    e1, e2, e3 = p1 - p0, p2 - p0, p3 - p0
    scale = float(np.linalg.norm(e1) * np.linalg.norm(e2) * np.linalg.norm(e3))
    if scale <= 0.0:
        return 0  # coincident donors
    normalized = float(np.dot(np.cross(e1, e2), e3)) / scale
    if abs(normalized) < _PLANARITY_EPS:
        return 0
    return 1 if normalized > 0 else -1


def metal_config_token(donor_positions) -> str:
    """``"|mc:+|"``, ``"|mc:-|"`` or ``""`` — the sidecar for a metal-centred configuration.

    Empty when the sphere is planar or under-determined, which is the same "emit nothing rather
    than emit a meaningless sign" rule the axial lane had to learn: a descriptor that fires on an
    achiral centre is over-sensitive, and over-sensitivity is indistinguishable from a bug once
    it reaches a corpus.
    """
    sign = metal_config_sign(donor_positions)
    if sign == 0:
        return ""
    return "|mc:+|" if sign > 0 else "|mc:-|"
