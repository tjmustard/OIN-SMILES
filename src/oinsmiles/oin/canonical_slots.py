"""Canonical coordination-slot labeling: the proper-rotation group and the lex-min relabeling.

WHAT THIS IS FOR
================
A donor's ``{n}`` slot number must be a **graph invariant**, not a property of how the input
XYZ happened to be oriented or numbered. Today it is neither: slots come from a geometric
Kabsch fit of donor direction-vectors to an idealized template, which is distortable and
conformer-dependent. Two conformers of one complex therefore get two different (but
rotation-related) labelings, and the emitted OIN strings differ even though the isomer does
not -- 315 molecules of the v0.4.4 capstone sweep (``slot_renumber``).

THE RULE
========
Fix the coordination geometry from the 3D fit -- that determines the **isomer**: which
vertices are occupied and how they sit relative to each other (adjacent vs opposite =
cis/trans; facial vs meridional = fac/mer). Then choose, **among only the vertex labelings
the geometry's own symmetry permits**, the lexicographically minimal one.

The nuance that makes this correct rather than merely intuitive: **priority cannot be
applied independently of geometry.** Sorting ligands by priority and stamping them onto
slots 0,1,2,... would erase cis/trans and fac/mer, which are *defined* by which vertices
equal-priority ligands occupy. The relabeling must run inside the freedom the geometry
allows -- the **proper-rotation group** (order 24 for an octahedron). Proper rotations
preserve the isomer and the winding sense; they only permute symmetry-equivalent vertices,
which is exactly the freedom that may be canonicalized.

WHY PROPER ROTATIONS ONLY
=========================
An *improper* operation (a reflection) maps a structure to its mirror image. Folding over
reflections would collapse enantiomers -- destroying exactly the stereochemistry the
v0.4.5 injectivity work exists to capture -- and would flip the eta winding sense. So
``derive_rotation_group`` filters to ``det > 0`` for spanning vertex sets. See the
docstring there for why planar/linear sets are handled differently.

RELATIONSHIP TO ``compare.py``
==============================
``oin/compare.py`` has computed the lex-min colored-vertex *signature* since v0.4.4 (that
is what makes the comparison KEY fac/mer-aware). It threw away the permutation that
achieved the minimum. That permutation **is** the canonical slot relabeling, so this module
owns the machinery, returns the permutation as well as the signature, and ``compare.py``
imports from here. Moving it also gives the geometry vertex table a single home
(open debt TD-005): ``utils.oin_aligner.TEMPLATE_SPECS`` keeps its ``ref`` vectors, which
only winding needs, and ``tests/unit/test_canonical_slots.py`` cross-checks the ``pos``
directions against ``GEOMETRY_VERTICES``.

Import graph is deliberately light -- **numpy only**, no RDKit, no aligner -- so both
``compare.py`` (which advertises a light graph) and the encoder can import it freely.
"""

import numpy as np

# Vertex direction vectors per geometry tag; slot i == vertex i. Mirrors the ``pos``
# entries of ``utils.oin_aligner.TEMPLATE_SPECS``, which is the encoder-side table that
# additionally carries the ``ref`` vectors used for winding.
GEOMETRY_VERTICES: dict[str, list[list[float]]] = {
    "LIN": [[0, 0, 1], [0, 0, -1]],
    "TPL": [[0, 1, 0], [0.8660254, -0.5, 0], [-0.8660254, -0.5, 0]],
    "SPL": [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
    "TET": [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]],
    "TPY": [[0, 0, 1], [0, 1, 0], [0.8660254, -0.5, 0], [-0.8660254, -0.5, 0]],
    "TBP": [[0, 0, 1], [0, 0, -1], [0, 1, 0], [0.8660254, -0.5, 0], [-0.8660254, -0.5, 0]],
    "SPY": [
        [0, 0, 1],
        [0.9659258, 0, -0.2588190],
        [-0.9659258, 0, -0.2588190],
        [0, 0.9659258, -0.2588190],
        [0, -0.9659258, -0.2588190],
    ],
    "OCT": [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]],
    "PBP": [
        [0, 0, 1],
        [0, 0, -1],
        [1, 0, 0],
        [0.309017, 0.951057, 0],
        [-0.809017, 0.587785, 0],
        [-0.809017, -0.587785, 0],
        [0.309017, -0.951057, 0],
    ],
    "SQA": [
        [-0.5773503, 0.5773503, 0.5773503],
        [0.5773503, 0.5773503, 0.5773503],
        [-0.5773503, -0.5773503, 0.5773503],
        [0.5773503, -0.5773503, 0.5773503],
        [-0.8141210, 0.0, -0.5773503],
        [0.0, -0.8141210, -0.5773503],
        [0.8141210, 0.0, -0.5773503],
        [0.0, 0.8141210, -0.5773503],
    ],
    "TCT": [
        [0.7555736, 0.0, 0.6550638],
        [-0.3780710, 0.6546229, 0.6546229],
        [-0.3780710, -0.6546229, 0.6546229],
        [0.7555736, 0.0, -0.6550638],
        [-0.3780710, 0.6546229, -0.6546229],
        [-0.3780710, -0.6546229, -0.6546229],
        [-1.0, 0.0, 0.0],
        [0.5821873, -0.8130547, 0.0],
        [0.5821873, 0.8130547, 0.0],
    ],
}

# Sorts above every ASCII letter and bracket, so an occupied vertex always beats an empty
# one in the lex-min and the signature is dominated by where the ligands actually are.
VERTEX_SENTINEL = ("~", "~", "~")

_GROUP_CACHE: dict[str, list | None] = {}


def derive_rotation_group(vertices, tol: float = 1e-3):
    """Proper-rotation vertex permutations of one idealized coordination polyhedron.

    ``vertices`` is the ordered list of vertex direction vectors (slot i == vertex i).
    A vertex permutation that preserves the full Gram (pairwise dot-product) matrix is
    realized by some orthogonal map. For a rank-3 (spanning) vertex set that map is
    unique, so keep the permutation iff it is a proper rotation (``det == +1``). For a
    planar/linear set (SPL, TPL, LIN) the out-of-plane direction is free, so every
    Gram-preserving permutation extends to a proper 3D rotation -- keep them all.

    Returns a sorted list of permutation tuples ``perm`` where ``perm[v]`` is the image
    index of vertex ``v``. Sorted order matters: it makes the lex-min argmin in
    ``lexmin_vertex_signature`` deterministic (see that docstring).
    """
    V = np.asarray(vertices, dtype=float)
    V = V / np.linalg.norm(V, axis=1, keepdims=True)
    n = len(V)
    gram = V @ V.T
    rank = int(np.linalg.matrix_rank(V, tol=1e-6))

    # Prune candidate images by a rotation-invariant per-vertex fingerprint.
    fp = [tuple(np.round(np.sort(gram[i]), 3)) for i in range(n)]
    candidates = [[j for j in range(n) if fp[j] == fp[i]] for i in range(n)]

    # An independent basis pins the (for rank 3, unique) realizing linear map.
    basis: list[int] = []
    for i in range(n):
        if np.linalg.matrix_rank(V[basis + [i]], tol=1e-6) == len(basis) + 1:
            basis.append(i)
        if len(basis) == rank:
            break
    binv = np.linalg.inv(V[basis].T) if rank == 3 else None

    perms: list[tuple[int, ...]] = []
    assign = [-1] * n

    def _backtrack(i: int) -> None:
        if i == n:
            perm = tuple(assign)
            if rank < 3:
                perms.append(perm)
                return
            rot = V[[perm[b] for b in basis]].T @ binv
            if np.linalg.det(rot) > 0 and np.allclose(V @ rot.T, V[list(perm)], atol=tol):
                perms.append(perm)
            return
        used = assign[:i]
        for j in candidates[i]:
            if j in used:
                continue
            if all(abs(gram[i][k] - gram[j][assign[k]]) < tol for k in range(i)):
                assign[i] = j
                _backtrack(i + 1)
                assign[i] = -1

    _backtrack(0)
    return sorted(perms)


def geometry_rotation_group(geo: str):
    """Cached proper-rotation permutation group for a geometry tag, or ``None`` if unknown.

    ``None`` -> the caller keeps absolute slots (identity-only folding). That is
    conservative: an unknown geometry can only *over*-split (miss a benign rotation), never
    wrongly collapse two isomers, so it can never reintroduce fac/mer blindness.
    """
    if geo not in _GROUP_CACHE:
        verts = GEOMETRY_VERTICES.get(geo)
        _GROUP_CACHE[geo] = derive_rotation_group(verts) if verts else None
    return _GROUP_CACHE[geo]


def geometry_vertex_count(geo: str) -> int:
    """Number of vertices the idealized polyhedron has, or 0 for an unknown geometry."""
    verts = GEOMETRY_VERTICES.get(geo)
    return len(verts) if verts else 0


def lexmin_vertex_signature(geo: str, vcolor: dict) -> tuple[tuple, tuple[int, ...]]:
    """Lex-min colored-vertex signature **and** the permutation that achieves it.

    ``vcolor`` maps an occupied slot index to that vertex's color (any comparable value;
    the key uses ``(fragment_body, donor_element, winding)``).

    Two slot labelings related by a proper rotation of the idealized polyhedron (conformer
    drift) produce the same signature; labelings not so related (fac vs mer) do not.

    Returns ``(signature, perm)`` where ``perm[old_slot] == new_slot``.

    **Determinism of the tie-break.** Several permutations can achieve the same minimal
    signature -- that happens exactly when the colored polyhedron has a symmetry, i.e. when
    the choice between them is a genuine automorphism and every choice yields the same
    string. The group arrives sorted and the comparison below is strict ``<``, so the
    FIRST (lexicographically smallest) achieving permutation wins. This must stay a
    property of the *permutation*, never of the color content: in the Y2 wave a tie broken
    on a stereochemical sign made a token reflection-invariant and silently destroyed the
    chirality it encoded. Proper rotations preserve winding sense, so relabeling by one of
    these permutations cannot flip a stereo descriptor.
    """
    group = geometry_rotation_group(geo)
    nverts = geometry_vertex_count(geo)
    if vcolor:
        nverts = max(nverts, max(vcolor) + 1)
    if not group:
        group = [tuple(range(nverts))]

    best = None
    best_perm: tuple[int, ...] = tuple(range(nverts))
    for perm in group:
        arr = [VERTEX_SENTINEL] * nverts
        for slot, color in vcolor.items():
            dest = perm[slot] if slot < len(perm) else slot
            arr[dest] = color
        candidate = tuple(arr)
        if best is None or candidate < best:
            best = candidate
            best_perm = perm
    # `group` is non-empty (guarded above), so the loop always assigns `best`.
    assert best is not None
    return best, best_perm


def canonical_slot_permutation(geo: str, vcolor: dict) -> dict[int, int]:
    """Canonical relabeling ``{old_slot: new_slot}`` for every slot present in ``vcolor``.

    This is the encoder-facing entry point, and the one Lanes 5 and 6 should use to answer
    "what is this donor's canonical slot index?" -- ``canonical_slot_permutation(...)[slot]``.
    Do not re-derive it: a second derivation is a second thing that can drift.

    A dict rather than the bare tuple because the group's permutations only cover the
    geometry's own vertices, while a real structure can carry a slot index beyond that
    count (a coordination number the template does not model). Such slots map to
    themselves; returning a dict makes that explicit instead of leaving callers to
    rediscover the rule.

    An unknown geometry yields the identity mapping, so the caller's slots are unchanged --
    the same conservative degradation as ``geometry_rotation_group``.
    """
    _sig, perm = lexmin_vertex_signature(geo, vcolor)
    return {slot: (perm[slot] if slot < len(perm) else slot) for slot in vcolor}
