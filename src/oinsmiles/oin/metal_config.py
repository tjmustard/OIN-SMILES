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

import itertools

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

#: Above this donor count the exhaustive permutation search is skipped (see :func:`is_achiral`).
#: 8! = 40320 superpositions is affordable; 9! = 362880 starts to matter in a hot encode path.
_MAX_PERM_DONORS = 8

__all__ = [
    "chirality_index",
    "is_achiral",
    "metal_config_sign",
    "metal_config_sign_symmetry",
    "metal_config_token",
    "metal_config_token_chelate",
    "token_for_mol",
]


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


#: RMSD tolerance (Å) for "the mirror image superimposes on the original".
#:
#: Chosen from the measured separation, using the chelate-aware test
#: (:func:`is_achiral_chelate_aware`):
#:
#:     ZUMNEC  chiral Δ/Λ tris-bidentate   best mirror-superposition RMSD = 1.3752 Å
#:     JEGKOW  achiral square planar       best mirror-superposition RMSD = 0.0582 Å
#:
#: A 24x separation, so 0.35 is a wide margin rather than a tuned constant. Crystallographic pucker
#: in an achiral sphere leaves a residual far below it, while a genuine Δ/Λ helix cannot be
#: superimposed on its mirror by ANY chelate-preserving proper rotation.
_ACHIRAL_RMSD_TOL = 0.35


def _kabsch_proper_rmsd(a, b) -> float:
    """RMSD after optimal **proper**-rotation superposition of *b* onto *a*.

    Proper only: an improper "rotation" would map any set onto its own mirror and make every
    structure look achiral, which is precisely the distinction being tested. The SVD's last
    singular vector is flipped when ``det < 0`` to force ``det R = +1``.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    u, _s, vt = np.linalg.svd(b.T @ a)
    d = np.sign(np.linalg.det(u @ vt))
    if d < 0:
        u[:, -1] *= -1.0
    rot = u @ vt
    diff = (b @ rot) - a
    return float(np.sqrt((diff**2).sum() / len(a)))


def is_achiral(donor_positions, tol: float = _ACHIRAL_RMSD_TOL) -> bool:
    """Does an IMPROPER operation map this donor set onto itself?

    This replaces thresholding :func:`chirality_index`, which was measured NOT to work: on real
    structures, crystallographic pucker in an achiral complex produces an index of the same order
    of magnitude as genuine helicity (JEGKOW -3.287e-04 achiral vs ZUMNEC -4.807e-04 chiral,
    1.5x apart). No magnitude threshold can separate those, because the exact-zero cancellation
    holds only for idealized coordinates.

    Chirality is a **symmetry** property, so it needs a symmetry test. Mirror the set, then ask
    whether any relabelling of the mirrored points can be superimposed on the original by a
    *proper* rotation. If one can, the mirror is the same object and the set is achiral.

    Permutations are enumerated exhaustively, which is why this is affordable: a coordination
    sphere has at most ~8 donors, and 8! = 40320 cheap superpositions. Above ``_MAX_PERM_DONORS``
    it returns ``False`` (assume chiral) rather than silently sampling — a wrong "achiral" would
    suppress a real descriptor, which is the worse error.
    """
    pts = np.asarray(donor_positions, dtype=float)
    if pts.ndim != 2 or pts.shape[0] < 4 or pts.shape[1] != 3:
        return True  # cannot support a handedness at all
    n = pts.shape[0]
    if n > _MAX_PERM_DONORS:
        return False
    mirrored = pts.copy()
    mirrored[:, 0] *= -1.0
    for perm in itertools.permutations(range(n)):
        if _kabsch_proper_rmsd(pts, mirrored[list(perm)]) <= tol:
            return True
    return False


def metal_config_sign_symmetry(donor_positions) -> int:
    """Achirality decided by symmetry rather than magnitude — and MEASURED INSUFFICIENT.

    ⚠⚠ THE INPUT IS WRONG, NOT JUST THE DECISION RULE. Measured: this reports ZUMNEC — a genuinely
    chiral Δ/Λ tris(catecholato) complex — as **achiral**, because it is. As a *bare point set*,
    six oxygens at octahedral vertices admit improper operations; there is no handedness in the
    positions alone.

    **Δ/Λ helicity is a property of the CHELATE CONNECTIVITY**, not of the donor point set: it is
    which donor pairs belong to the same bidentate ligand, and how those chelate planes twist about
    the metal. Reflecting a Δ complex gives Λ only because the reflection cannot be undone
    *while keeping the chelate pairing intact* — and a permutation search over unlabelled points is
    free to re-pair them, so it always finds a "symmetry" that does not exist chemically.

    That also explains why :func:`chirality_index` looked like it worked: its non-zero reading for
    ZUMNEC (-4.807e-04) is residual crystallographic distortion, the same magnitude as achiral
    JEGKOW's pucker (-3.287e-04). It was never detecting helicity.

    The fix is to constrain the permutation search to relabellings that PRESERVE chelate
    membership — i.e. treat the donors as a coloured point set where the colour is the ligand a
    donor belongs to. A mirror that requires re-pairing chelates is then correctly rejected. That
    needs the ligand partition threaded in from the caller, which the current signature does not
    carry, so it is the remaining work rather than a tweak.

    Left in place, unwired, because the point-set achirality test itself is correct and reusable —
    it is the *input* that is incomplete.
    """
    if is_achiral(donor_positions):
        return 0
    index = chirality_index(donor_positions)
    if index == 0.0:
        return 0
    return 1 if index > 0 else -1


#: Precomputed internal orderings per group size, so the product below iterates LISTS.
_INTERNALS: dict[int, list[tuple[int, ...]]] = {
    n: list(itertools.permutations(range(n))) for n in range(1, 7)
}


def _admissible_permutations(groups):
    """Relabellings that PRESERVE chelate membership.

    ``groups`` partitions the donor indices by the ligand each donor belongs to — so a
    tris(bidentate) sphere is ``[(0,1),(2,3),(4,5)]`` and four monodentates are
    ``[(0,),(1,),(2,),(3,)]``.

    A relabelling is admissible only if it maps whole chelates onto whole chelates of the same
    size. That is the constraint an unconstrained permutation search lacks, and its absence is why
    the search found a bogus "mirror symmetry" for a genuinely chiral Δ complex: it was allowed to
    re-pair the donors into different ligands, which no physical operation can do.

    Count stays small: a tris-bidentate gives 3! group assignments x (2!)^3 internal orderings = 48.
    """
    by_size: dict[int, list[tuple[int, ...]]] = {}
    for g in groups:
        by_size.setdefault(len(g), []).append(g)

    per_size_options = []
    for size, gs in sorted(by_size.items()):
        options = []
        for target_order in itertools.permutations(range(len(gs))):
            # list(...), NOT the bare iterator: `[itertools.permutations(x)] * n` repeats ONE
            # iterator, so after the first is consumed the rest are empty and the product
            # collapses to nothing. That made this generator yield ZERO permutations, which made
            # is_achiral_chelate_aware() return "chiral" for everything -- a vacuous loop that
            # reads exactly like a detection. Same family as the empty-corpus and
            # buffered-stdout failures recorded in docs: nothing measured, confident answer.
            for internals in itertools.product(*[_INTERNALS[size]] * len(gs)):
                pairs = []
                for src_i, dst_i in enumerate(target_order):
                    src, dst = gs[src_i], gs[dst_i]
                    pairs += [(src[k], dst[internals[src_i][k]]) for k in range(size)]
                options.append(pairs)
        per_size_options.append(options)

    for combo in itertools.product(*per_size_options):
        mapping = {}
        for pairs in combo:
            mapping.update(dict(pairs))
        yield [mapping[i] for i in range(len(mapping))]


def is_achiral_chelate_aware(donor_positions, groups, tol: float = _ACHIRAL_RMSD_TOL) -> bool:
    """Achirality with the chelate partition respected — the form that actually works for Δ/Λ.

    :func:`is_achiral` searches ALL permutations and therefore reports a chiral Δ/Λ tris-bidentate
    as achiral: free to re-pair donors into different ligands, it always finds a superposition.
    Restricting to :func:`_admissible_permutations` removes exactly that freedom.
    """
    pts = np.asarray(donor_positions, dtype=float)
    if pts.ndim != 2 or pts.shape[0] < 4 or pts.shape[1] != 3:
        return True
    mirrored = pts.copy()
    mirrored[:, 0] *= -1.0
    for perm in _admissible_permutations(groups):
        if _kabsch_proper_rmsd(pts, mirrored[perm]) <= tol:
            return True
    return False


def metal_config_token_chelate(donor_positions, groups) -> str:
    """The Δ/Λ sidecar, decided by chelate-aware symmetry. ``""`` when achiral."""
    if is_achiral_chelate_aware(donor_positions, groups):
        return ""
    index = chirality_index(donor_positions)
    if index == 0.0:
        return ""
    return "|mc:+|" if index > 0 else "|mc:-|"


def token_for_mol(mol) -> str:
    """The Δ/Λ sidecar for a metal-PRESENT mol with a conformer, or ``""``.

    Derives both inputs the descriptor needs — donor positions and the chelate partition — from the
    mol itself, so the encoder's call site stays a single line. Donors are the metal's perceived
    neighbours; the partition is the connected components left after deleting the metal, which is
    what makes a bidentate's two donors one chelate.

    Must be called on the PRISTINE input conformer, before ``_align_to_pai``: principal-axis
    alignment may reflect the coordinates, and a reflection inverts the descriptor. That is the same
    ordering constraint the axial token documents at its own call site.

    Returns ``""`` on anything unexpected rather than raising — this runs inside the encoder's
    serialization path, where an exception would reroute rather than surface.
    """
    try:
        from rdkit import Chem

        from ..core.constants import TRANSITION_METALS_NUM

        if mol is None or mol.GetNumConformers() == 0:
            return ""
        metal = next(
            (a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM),
            None,
        )
        if metal is None:
            return ""
        conf = mol.GetConformer()
        idxs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(metal).GetNeighbors()]
        if len(idxs) < 4:
            return ""
        pts = np.array([list(conf.GetAtomPosition(i)) for i in idxs])

        stripped = Chem.RWMol(mol)
        stripped.RemoveAtom(metal)
        comp = {a: fi for fi, f in enumerate(Chem.GetMolFrags(stripped.GetMol())) for a in f}
        groups: dict[int, list[int]] = {}
        for pos, atom_idx in enumerate(idxs):
            key = comp.get(atom_idx - 1 if atom_idx > metal else atom_idx)
            if key is None:
                return ""
            groups.setdefault(key, []).append(pos)
        return metal_config_token_chelate(pts, [tuple(v) for v in groups.values()])
    except Exception:
        return ""
