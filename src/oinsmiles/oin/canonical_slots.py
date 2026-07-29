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

Import graph is deliberately light -- the module header is **numpy only**, no RDKit, no
aligner -- so both ``compare.py`` (which advertises a light graph) and the encoder can
import it freely. The string-level half of the module (from ``_relabel_slots`` down) needs
``compare`` and ``inline``, and imports them *inside* the functions: that keeps the header
light and keeps ``compare.py`` -> ``canonical_slots.py`` acyclic, at the cost of nothing,
since any caller of those functions has already loaded both.
"""

import re

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
    """Canonical relabeling ``{old_slot: new_slot}`` from a colored-vertex map alone.

    A dict rather than the bare tuple because the group's permutations only cover the
    geometry's own vertices, while a real structure can carry a slot index beyond that
    count (a coordination number the template does not model). Such slots map to
    themselves; returning a dict makes that explicit instead of leaving callers to
    rediscover the rule.

    An unknown geometry yields the identity mapping, so the caller's slots are unchanged --
    the same conservative degradation as ``geometry_rotation_group``.

    .. warning::
       **This is not the function to ask "what is donor *d*'s canonical slot?"** -- use
       :func:`canonical_slot_map`. The lex-min *signature* is invariant, but the
       permutation achieving it need not be unique: whenever the colored polyhedron has a
       nontrivial color-preserving rotation stabilizer, several permutations tie, and the
       tie is broken here on the permutation tuple, which is a property of the *incoming*
       labeling. Two presentations of one molecule arrive with different incoming
       labelings, so they can pick different members of the tie -- and members of the tie
       differ exactly by swapping same-colored donors, which for two chemically distinct
       donors of ONE ligand is a genuinely different emitted string.
       :func:`canonical_slot_map` refines the tie on the rendered output, which is
       invariant by construction. See its docstring for the proof.
    """
    _sig, perm = lexmin_vertex_signature(geo, vcolor)
    return {slot: (perm[slot] if slot < len(perm) else slot) for slot in vcolor}


# --------------------------------------------------------------------------------------
# String level: the encoder post-pass, and the helper Lanes 5/6 consume.
# --------------------------------------------------------------------------------------
#
# Everything below works on a finished **inline** OIN string. That is deliberate: it is
# exactly the representation ``compare._parse_vertex_colors`` already reads, so the
# encoder's canonicalization and the comparison key's canonicalization consume the same
# bytes through the same function and cannot drift apart. RDKit and ``OINInlineHandler``
# are imported lazily inside the functions so the module header stays numpy-only (see the
# module docstring) and so ``compare.py`` -> ``canonical_slots.py`` stays acyclic.

#: The opt-in axial/atropisomer suffix (``OIN_EMIT_AXIAL``) is appended after the inline
#: string is built. It carries no slot, so it is set aside and re-appended verbatim.
_AXIAL_SUFFIX_RE = re.compile(r"\s*\|ax:[+\-]*\|\s*$")

#: Sorts after every real slot index, so an uncoordinated fragment lands at the end.
_NO_SLOT = float("inf")


def _relabel_slots(frag: str, mapping: dict) -> str:
    """Rewrite ``{n}`` / ``{n>}`` / ``{n<}`` / ``{n^}`` to ``{mapping[n]}``, winding intact.

    The winding character is preserved **verbatim**. ``_parse_vertex_colors`` folds ``^``
    to ``>`` for coloring purposes only; the emitted string must keep whichever character
    the aligner computed. That is safe because winding is measured against the ring's own
    metal->centroid axis (``oin_aligner._determine_winding``), never against the slot's
    template direction, so a relabeling cannot invalidate it -- and in any case the group
    contains only proper rotations, which preserve circulation sense.
    """
    from .inline import OINInlineHandler

    def _sub(m):
        slot = int(m.group(1))
        return "{" + str(mapping.get(slot, slot)) + (m.group(2) or "") + "}"

    return OINInlineHandler.SLOT_REGEX.sub(_sub, frag)


def _min_slot(frag: str):
    """Lowest slot index a fragment carries, or ``inf`` if it coordinates nothing."""
    from .inline import OINInlineHandler

    slots = [int(m.group(1)) for m in OINInlineHandler.SLOT_REGEX.finditer(frag)]
    return min(slots) if slots else _NO_SLOT


def _skeleton_ranks(mol):
    """``CanonicalRankAtoms`` classes of a fragment's **constitutional skeleton**, or ``None``.

    This is the v0.4.14 widening (``OIN_RESONANCE_DONOR_FOLD``). It exists because the strict
    ranking below reads a **frozen resonance form** as two inequivalent donors:

        acac  ``CC(=O{0})C=C(C)O{1}``          strict ranks 2 / 3   -> never exchanged
        sulfonate ``O{0}S(=O)(=O{2})...``      strict ranks 2 / 0   -> never exchanged

    Both write one donor as a ketone and its partner as an enol/anion, which is a property of
    the Kekulé structure the perceiver happened to emit, not of the ligand -- the real ligand is
    delocalized and its two donors are the same atom. Measured on the v0.4.8 corpus: **101 of
    the 103** ``key_equal/slot_renumber`` molecules the v0.4.13 fold cannot reach fail on
    exactly this, ``same_colour_DIFFERENT_rank``.

    The skeleton erases only that bookkeeping -- bond orders, aromatic flags, formal charges and
    hydrogen counts -- and keeps **connectivity, element and chiral tag**. So it merges what
    resonance makes equivalent and still refuses what constitution makes different:

        acac O / O            merge      (resonance pair)
        carboxylate O / O     merge      (resonance pair)
        ester  -O- / =O       REFUSE     (one is 2-connected, one terminal)
        ether O / ketone O    REFUSE     (different neighbourhoods)
        amide  N / O          REFUSE     (different elements)

    ⚠ **Chirality is deliberately retained** -- tags survive the flattening and
    ``CanonicalRankAtoms``' ``includeChirality`` default of ``True`` consumes them, so this
    widening does not DISCARD stereochemical information the v0.4.11 strict ranking already
    used. That is the whole of the claim, and the stronger-sounding version is false: measured
    on a diol pair, the C2-symmetric ``(R,R)`` arms do **not** merge (over-conservative, a missed
    fold -- the safe direction) while the meso ``(R,S)`` arms **do**, because they are
    enantiotopic. Folding an enantiotopic pair is a reflection, and the guard against it is
    ``fold_parity``'s per-molecule veto, not this ranking. Zeroing charges and hydrogens is safe
    here; clearing chiral tags would additionally throw away what the strict ranking had.

    Returns ``None`` on any RDKit failure, which the caller treats as "do not widen" -- the
    same conservative degradation as an unknown geometry.
    """
    from rdkit import Chem

    try:
        rw = Chem.RWMol(mol)
        for b in rw.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
        for a in rw.GetAtoms():
            a.SetFormalCharge(0)
            a.SetNoImplicit(True)
            a.SetNumExplicitHs(0)
            a.SetIsAromatic(False)
        skel = rw.GetMol()
        # Ring perception only. A full sanitize would run the valence check, and the flattened
        # skeleton is deliberately valence-invalid (a pyridine N with three single bonds and no
        # charge); rejecting it there would silently disable the widening on every aromatic
        # ligand -- i.e. on most of the population this exists for.
        Chem.SanitizeMol(skel, Chem.SANITIZE_SYMMRINGS | Chem.SANITIZE_ADJUSTHS)
        return list(Chem.CanonicalRankAtoms(skel, breakTies=False))
    except Exception:  # noqa: BLE001 -- an unrankable skeleton simply does not widen the fold
        return None


def _merge_classes(slots, key_of, alt_key_of):
    """Group ``slots`` by ``key_of``, then merge groups that also agree under ``alt_key_of``.

    Union-find rather than a single composite key, because the widening must be a **coarsening**
    of the strict grouping and nothing else. Keying on the pair would let a slot whose alternate
    key differs land in a *different* bucket than it does today, which could make the candidate
    set smaller rather than larger and lose a labeling the shipped encoder can reach. Two slots
    end up together iff they agree under ``key_of`` **or** under ``alt_key_of`` (transitively),
    so the partition can only get coarser.
    """
    parent = {s: s for s in slots}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    for keyfn in (key_of, alt_key_of):
        seen: dict = {}
        for s in slots:
            k = keyfn(s)
            if k is None:
                continue
            if k in seen:
                union(seen[k], s)
            else:
                seen[k] = s

    out: dict = {}
    for s in slots:
        out.setdefault(find(s), []).append(s)
    return list(out.values())


def _donor_swap_permutations(frags: list[str], vcolor: dict, cap: int = 4096) -> list[dict]:
    """Slot permutations exchanging **interchangeable donors within one fragment**.

    This is the v0.4.11 widening (``OIN_CANONICAL_DONOR_FOLD``), and it is the only place in
    this module that folds past the geometry's own proper-rotation group. Two slots of *one*
    fragment may be exchanged when both hold:

    (a) their donor atoms lie in the same ``CanonicalRankAtoms(breakTies=False)`` symmetry
        class of that fragment -- or, with ``OIN_RESONANCE_DONOR_FOLD`` on, in the same class
        of its constitutional skeleton (see :func:`_skeleton_ranks`), and
    (b) they carry the same vertex colour.

    Always includes the identity, so the caller's candidate set is a superset of the
    rotation-only one and the fold can never *lose* a labeling the old code could reach.

    **Why this is not over-folding.** Condition (a) is a graph automorphism of the fragment
    computed with ``includeChirality`` at its default ``True``, so two constitutionally
    equivalent branches with *different* configurations land in different classes and are
    never exchanged. Condition (b) keeps the exchange inside one colour, so the occupied
    vertex set and the colored-vertex signature are both untouched -- the fold acts entirely
    within the signature's kernel, which is precisely why the rotation-only post-pass could
    not reach this class (``same_vcolor_identical``, 496/496 of the v0.4.11 population).
    Nothing that distinguishes an isomer, an enantiomer or a winding sense is a function of
    which of two automorphic donors carries which integer.

    **Why the emitted string stays presentation-invariant.** Let ``G`` be the rotation group
    and ``D`` the set returned here. A second presentation arrives as ``L' = g.L``. The
    buckets above are keyed on *atom* symmetry classes, which the relabeling does not touch,
    so ``D' = g D g^-1``. The candidate set for ``L'`` is
    ``{p.d'.L' : p in G, d' in D'} = {p.g.d.L} = {(p.g).d.L}``, and ``p -> p.g`` is a
    bijection of ``G``, so it equals ``{p.d.L}`` -- the candidate set for ``L``. The minimum,
    hence the emitted string, is identical. This is the same argument
    :func:`canonical_slot_relabeling` makes for the rotation group, extended by one factor.

    ``cap`` bounds the combinatorial product (a fragment with ``k`` interchangeable donors
    contributes ``k!``, and the fragments multiply). Over the cap the fold degrades to the
    identity rather than spending unbounded time -- conservative in the same direction as an
    unknown geometry.
    """
    from itertools import permutations, product

    from rdkit import Chem

    from .compare import _parse_fragment
    from .inline import OINInlineHandler, _count_smiles_atoms_before
    from .levers import lever_enabled

    resonance = lever_enabled("OIN_RESONANCE_DONOR_FOLD")

    buckets_all: list[tuple[int, ...]] = []
    for frag in frags:
        if OINInlineHandler.METAL_REGEX.search(frag):
            continue
        slot_atoms: dict[int, set] = {}
        for m in OINInlineHandler.SLOT_REGEX.finditer(frag):
            prefix = OINInlineHandler.SLOT_REGEX.sub("", frag[: m.start()])
            slot_atoms.setdefault(int(m.group(1)), set()).add(
                _count_smiles_atoms_before(prefix, len(prefix))
            )
        if len(slot_atoms) < 2:
            continue  # one donor cannot be exchanged with anything
        mol = _parse_fragment(OINInlineHandler.SLOT_REGEX.sub("", frag))
        if mol is None:
            continue
        try:
            ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
        except Exception:  # noqa: BLE001  -- an unrankable fragment simply does not fold
            continue
        usable = [
            slot
            for slot, atoms in sorted(slot_atoms.items())
            if not any(i >= len(ranks) for i in atoms)
        ]

        def _strict_key(slot, _ranks=ranks, _atoms=slot_atoms):
            return (tuple(sorted(_ranks[i] for i in _atoms[slot])), vcolor.get(slot))

        skel = _skeleton_ranks(mol) if resonance else None

        def _skel_key(slot, _skel=skel, _atoms=slot_atoms):
            # ``None`` means "this key has no opinion", and `_merge_classes` skips it. That is
            # what makes the lever-off path byte-identical: with no skeleton there is no second
            # pass, so the partition is exactly the strict grouping it was before.
            if _skel is None or any(i >= len(_skel) for i in _atoms[slot]):
                return None
            return (tuple(sorted(_skel[i] for i in _atoms[slot])), vcolor.get(slot))

        for group in _merge_classes(usable, _strict_key, _skel_key):
            if len(group) > 1:
                buckets_all.append(tuple(sorted(group)))

    if not buckets_all:
        return [{}]

    total = 1
    for b in buckets_all:
        for k in range(2, len(b) + 1):
            total *= k
        if total > cap:
            return [{}]

    buckets_all.sort()
    out = []
    for combo in product(*(permutations(b) for b in buckets_all)):
        mapping = {}
        for bucket, perm in zip(buckets_all, combo):
            mapping.update(dict(zip(bucket, perm)))
        out.append(mapping)
    return out


def _render(frags: list[str], metal_pos, mapping: dict):
    """Apply ``mapping`` and re-sort the fragments. Returns ``(string, sort_key_tuple)``.

    Fragment order is ``(minimum canonical slot, fragment text)`` with the metal fragment
    pinned first -- it is ``fragments[0]`` and that is a load-bearing project invariant
    (the generator, the inline parser and the comparison key all assume it).

    ``minimum canonical slot`` is by itself a **total** order on the coordinated fragments,
    because a slot belongs to exactly one fragment; the text term only ever decides between
    two *uncoordinated* fragments. So this sort replaces the old input-atom-index tie-break
    (``perception_tmc.get_input_order_key``) with a property of the molecule rather than of the
    file, and it does so without needing a ligand-body parse.
    """
    relabeled = [_relabel_slots(f, mapping) for f in frags]
    metal_frag = relabeled[metal_pos] if metal_pos is not None else None
    rest = [f for i, f in enumerate(relabeled) if i != metal_pos]
    keyed = sorted((_min_slot(f), f) for f in rest)
    ordered = ([metal_frag] if metal_frag is not None else []) + [f for _k, f in keyed]
    return ".".join(ordered), tuple(keyed)


def canonical_slot_relabeling(oin_string: str) -> tuple[dict[int, int], str]:
    """``({old_slot: new_slot}, canonical_string)`` for one emitted inline OIN string.

    The single implementation behind :func:`canonical_slot_map` (which wants the map) and
    :func:`canonicalize_oin_slots` (which wants the string). Returns the input unchanged
    with an empty map when there is nothing to canonicalize (no metal tag, no slots).

    **Why the minimum is taken over the rendered output and not just the signature.**
    Let ``G`` be the geometry's proper-rotation group and ``L`` the incoming slot labeling.
    A second presentation of the same complex arrives as ``L' = g . L`` for some ``g`` in
    ``G``. The candidate set here is ``{(sig(p.L), render(p.L)) : p in G}``; for ``L'`` it
    is ``{(sig(p.g.L), render(p.g.L)) : p in G}``, and since ``p -> p.g`` is a bijection of
    ``G`` the two sets are **equal**. So the minimum -- hence the emitted string -- is
    identical. Minimizing on the signature alone is not enough: several permutations can
    achieve the minimal signature (see :func:`canonical_slot_permutation`'s warning), and
    picking between them on the incoming labeling is exactly the input dependence this lane
    exists to remove.

    The signature stays the **primary** sort term even though the rendered term alone would
    already be invariant. That is what keeps the encoder and ``compare.py`` in lockstep: the
    emitted labeling achieves ``_polyhedron_signature``'s lex-min, so re-running the key on
    the emitted string finds the identity is already optimal.

    Folding is over **proper rotations only**, so this cannot merge two enantiomers, cannot
    flip a winding sense, and cannot collapse fac into mer (their labelings lie in different
    orbits, so their candidate sets are disjoint).
    """
    from .compare import _parse_vertex_colors, normalize_oin_for_comparison
    from .inline import OINInlineHandler

    m = _AXIAL_SUFFIX_RE.search(oin_string)
    body, suffix = (oin_string[: m.start()], oin_string[m.start() :]) if m else (oin_string, "")

    # Empty fragments (consecutive/trailing dots from an uncoordinated ligand that
    # serialized to nothing) carry no information and cannot be kept "in place" through a
    # re-sort; the comparison key drops them too (normalize_oin_for_comparison step 2).
    frags = [f for f in body.split(".") if f]
    if not frags:
        return {}, oin_string

    metal_pos = next(
        (i for i, f in enumerate(frags) if OINInlineHandler.METAL_REGEX.search(f)), None
    )
    _metal, geo, vcolor = _parse_vertex_colors(normalize_oin_for_comparison(body))
    if not vcolor:
        return {}, oin_string

    nverts = max(geometry_vertex_count(geo), max(vcolor) + 1)
    group = geometry_rotation_group(geo) or [tuple(range(nverts))]

    # v0.4.11: widen the candidate set by the within-fragment donor swaps, if enabled. With
    # the lever off this is exactly ``[{}]``, so every key below keeps its v0.4.5 value and
    # the emitted string and the returned map are both byte-identical to the old behaviour.
    from .levers import lever_enabled

    donor_perms = (
        _donor_swap_permutations(frags, vcolor)
        if lever_enabled("OIN_CANONICAL_DONOR_FOLD")
        else [{}]
    )

    best_key = None
    best_map: dict[int, int] = {}
    best_out = oin_string
    for perm in group:
        for dperm in donor_perms:
            # Donor swap first, then the rotation: ``mapping[s] = p(d(s))``. ``d`` only ever
            # exchanges same-coloured slots, so ``arr`` -- and hence the signature -- is
            # identical for every ``d``, which is why this widening is invisible to the
            # comparison key and changes only which labeling is emitted.
            mapping = {}
            for slot in vcolor:
                s = dperm.get(slot, slot)
                mapping[slot] = perm[s] if s < len(perm) else s
            arr = [VERTEX_SENTINEL] * nverts
            for slot, color in vcolor.items():
                arr[mapping[slot]] = color
            out, order_key = _render(frags, metal_pos, mapping)
            # Third and fourth terms: a final tie-break so the *map* is deterministic too.
            # They only ever separate candidates that already render identically, so they
            # can never make the emitted string depend on the incoming labeling. ``dperm``
            # is last so that with the lever off the key is unchanged from v0.4.5.
            key = (tuple(arr), order_key, perm, tuple(sorted(dperm.items())))
            if best_key is None or key < best_key:
                best_key, best_map, best_out = key, mapping, out
    return best_map, best_out + suffix


def canonical_slot_map(oin_string: str) -> dict[int, int]:
    """Return ``{slot_in_this_string: canonical_slot}`` -- the Lane 5 / Lane 6 entry point.

    Answering "what is donor *d*'s canonical slot index?" is a two-step join, and both
    steps already exist:

    1. take *d*'s slot as it appears in the emitted OIN string (the aligner's
       ``item["slot"]``, or the integer inside the ``{n}`` marker on *d*'s atom);
    2. ``canonical_slot_map(oin_string)[that_slot]``.

    Do **not** re-derive this from ``canonical_slot_permutation`` or from a fresh vertex
    coloring: a second derivation is a second thing that can drift, and the naive one is
    subtly wrong (see that function's warning).

    The map is the **identity** when ``OIN_CANONICAL_SLOTS`` is on, because the encoder has
    already applied it -- which is the point: a caller written against this helper is
    correct with the lever either way, and stays correct when the lever is promoted to
    default-on. Calling it on an already-canonical string is idempotent.

    **What is and is not unique (read this before building a stereo descriptor on it).**
    The canonical *string* is unique, always. The per-donor map is unique only up to the
    colored polyhedron's rotational automorphism group -- and exactly up to it, no more.
    fac-M(ppy)3 has a real C3 axis through the three equivalent ligands, so three donor
    labelings are equally canonical (measured: 3, matching the automorphism order, while
    the emitted string stays single-valued). That is not a defect to work around:

    * every member of that group is a **proper rotation**, so a descriptor that flips under
      reflection -- metal Delta/Lambda, an eta winding, an axial sign -- takes the *same*
      value on all of them. Folding here cannot destroy chirality the way the Y2 wave's
      sign-sorted axial token did, because that token was folded over a *reflection*.
    * a descriptor that is NOT invariant under the automorphism is not a property of the
      molecule in the first place; it is a property of which interchangeable ligand you
      decided to call "first".

    So: derive from the canonical arrangement, and let this map join a donor to it. Do not
    ask this map for an identity it cannot have.
    """
    mapping, _out = canonical_slot_relabeling(oin_string)
    return mapping


def canonicalize_oin_slots(oin_string: str) -> str:
    """Rewrite an inline OIN string with canonical slot labels and canonical fragment order.

    This is the encoder post-pass (``OIN_CANONICAL_SLOTS``), applied to the finished inline
    string in ``utils.perception_tmc.get_oin_string``. Idempotent.
    """
    _mapping, out = canonical_slot_relabeling(oin_string)
    return out
