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

THE DESCRIPTOR — a permutation-invariant pseudoscalar, NOT an ordered signed volume
==================================================================================
The first attempt ordered donors by canonical slot and took the signed volume of the first four.
It failed, measurably, on a HOMOLEPTIC complex: ZUMNEC's six O donors are symmetry-equivalent, so
every scalar ordering key ties, the tie falls to input order, and some resolutions differ by an
ODD permutation — which inverts a signed volume. The sign flipped under pure atom renumbering.

The fix was not a better ordering. It was to **stop needing one**: sum a pseudoscalar over EVERY
ordered 4-tuple (:func:`chirality_index`), which is permutation-invariant by construction.

* **invariant under proper rotation** — the triple product carries ``det R = +1``;
* **inverts under reflection** — an improper operation flips it, dot products unchanged;
* **invariant under any donor relabelling** — measured identical over 6 random permutations;
* **exactly 0 for an achiral arrangement** — a perfect square and an ideal octahedron both return
  ``+0.000e+00``. Achirality falls OUT of the index instead of needing a planarity test beside it,
  which is why square-planar ``JEGKOW`` emits nothing without a special case. That is chemically
  right: its coordination plane is a mirror plane, so four different donors give diastereomers,
  and the correct distinctness operator there is a donor swap, not reflection.

A caution the index also exposed: a signed volume of four LABELLED points is non-zero for any
non-coplanar set, so the old form reported handedness for a *regular* tetrahedron — which is
achiral (Td contains improper operations). "The labelling has an orientation" is not "the shape is
chiral", and only the invariant form can tell them apart.

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

#: Threshold on the normalized permutation-invariant chirality index (see
#: :func:`chirality_index`).
#:
#: Deliberately tiny, and the first value chosen (1e-3) was WRONG -- larger than the real signal,
#: so it would have called ZUMNEC achiral and silently emitted nothing for the very fixture the
#: lane exists for. Measured: an achiral point set cancels to **exactly 0.0** (a perfect square and
#: an ideal octahedron both return +0.000e+00, not rounding noise), while chiral ZUMNEC reads
#: 4.81e-4. The separation is therefore between exact zero and the signal, so the threshold only
#: needs to clear floating-point residue -- not to guess a magnitude.
_CHIRALITY_EPS = 1e-9

__all__ = ["chirality_index", "metal_config_sign", "metal_config_token"]


def chirality_index(donor_positions) -> float:
    """A **permutation-invariant** pseudoscalar for a donor point set.

    This exists because the ordering-based descriptor below cannot work for a HOMOLEPTIC complex.
    Measured on ZUMNEC (tris-catecholato Mo): its six O donors are symmetry-equivalent, so every
    scalar ordering key ties, the tie falls to input order, and some resolutions differ by an
    ODD permutation -- which inverts a signed volume. The sign then flipped under pure atom
    renumbering, 1 -> -1 on 2 of 4 shuffles.

    The framing "find a canonical ordering up to proper rotation only" was the wrong problem.
    **No ordering is needed at all.** Summing a pseudoscalar over EVERY ordered 4-tuple is
    permutation-invariant by construction: relabelling the donors permutes the terms of the sum
    without changing the total.

    The summand is the Osipov-Pickup-Dunmur chirality index form,

        (r_ij x r_kl) . r_il  *  (r_ij . r_jk)  *  (r_jk . r_kl)
        --------------------------------------------------------
                     (|r_ij| |r_jk| |r_kl|)^2

    which is a genuine pseudoscalar: the triple product changes sign under an improper operation
    while every dot product and magnitude is invariant, so the whole index **negates under
    reflection** and is **unchanged under proper rotation**. For an achiral point set the terms
    cancel exactly and it is 0 -- which is the property that makes the square-planar case fall out
    rather than needing a separate planarity test.

    Normalized by the term count so the magnitude is comparable across coordination numbers.
    Cost is O(n^4) in donors; with n <= 8 that is at most a few thousand terms.
    """
    pts = np.asarray(donor_positions, dtype=float)
    if pts.ndim != 2 or pts.shape[0] < 4 or pts.shape[1] != 3:
        return 0.0
    n = pts.shape[0]
    total = 0.0
    terms = 0
    for i in range(n):
        for j in range(n):
            if j == i:
                continue
            r_ij = pts[j] - pts[i]
            n_ij = float(np.linalg.norm(r_ij))
            if n_ij == 0.0:
                continue
            for k in range(n):
                if k in (i, j):
                    continue
                r_jk = pts[k] - pts[j]
                n_jk = float(np.linalg.norm(r_jk))
                if n_jk == 0.0:
                    continue
                for m in range(n):
                    if m in (i, j, k):
                        continue
                    r_kl = pts[m] - pts[k]
                    n_kl = float(np.linalg.norm(r_kl))
                    if n_kl == 0.0:
                        continue
                    r_il = pts[m] - pts[i]
                    pseudo = float(np.dot(np.cross(r_ij, r_kl), r_il))
                    weight = float(np.dot(r_ij, r_jk)) * float(np.dot(r_jk, r_kl))
                    total += pseudo * weight / (n_ij * n_jk * n_kl) ** 2
                    terms += 1
    return total / terms if terms else 0.0


def metal_config_sign(donor_positions) -> int:
    """``+1`` / ``-1`` for the handedness of a donor set, or ``0`` when planar/degenerate.

    Delegates to :func:`chirality_index`, so it needs **no canonical ordering at all** -- the
    donor positions may arrive in any order. That replaced an earlier signed-volume-of-the-first-
    four implementation which required a canonical slot ordering and was measured to flip sign
    under pure atom renumbering on a homoleptic complex (see ``chirality_index``).

    Fewer than four donors cannot define a chirality, and an achiral arrangement (square planar,
    ideal octahedron) cancels to exactly zero; both return ``0``. No separate planarity test is
    needed -- achirality falls out of the index rather than being detected alongside it.
    """
    index = chirality_index(donor_positions)
    if abs(index) < _CHIRALITY_EPS:
        return 0
    return 1 if index > 0 else -1


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
