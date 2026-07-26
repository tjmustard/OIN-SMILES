"""Order-independent tetrahedral stereo for the ligand-fragment rebuild (v0.4.5 Lane 8).

Why this exists
===============
``get_oin_string`` does not extract a ligand with ``Chem.PathToSubmol`` or
``RWMol.RemoveAtom``; it rebuilds every ligand from scratch:

1. heavy atoms are copied with ``mw.AddAtom(atom)``, which copies the parent
   atom's chiral tag **verbatim**;
2. hydrogens are dropped and folded into ``SetNumExplicitHs``;
3. bonds are re-added by iterating ``heavy_indices`` in **ascending parent index**
   and taking each neighbour with ``nbr_idx > old_idx``.

An RDKit chiral tag is not an absolute configuration. It is a parity *relative to
the order the neighbours appear on that atom*. Steps 2 and 3 both change that
order, and both are functions of the input atom numbering -- so the verbatim tag
copied in step 1 is re-interpreted against a different permutation. When that
permutation's parity differs, the same stereocentre in the same 3D structure
emits ``@`` from one input ordering and ``@@`` from another.

Measured on the v0.4.5 renumbering probe: the parent mol's perceived graph and
chiral tags were **identical** between the two orderings for 26 of 29 affected
molecules -- the divergence appeared for the first time inside this rebuild.

The fix
=======
Do not translate the tag; re-derive it. Copy the parent conformer's coordinates
onto the rebuilt fragment and let ``AssignAtomChiralTagsFromStructure`` stamp the
tag against the fragment's **own** neighbour order. Geometry is not a function of
atom numbering, so the result is order-independent by construction -- and, unlike
a descriptor made stable by being constant, it still inverts under reflection.

Scope guard
===========
Only atoms that **already** carry ``CHI_TETRAHEDRAL_CW``/``CCW`` are restamped,
and no other atom is touched. This can never add or remove stereochemistry: the
number of tagged atoms going in equals the number going out. Making the string
stable by dropping the tag is the failure mode this lane exists to prevent.
"""

import numpy as np
from rdkit import Chem

_TETRAHEDRAL = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)

# Below this |triple product| of the three *unit* neighbour vectors the centre is
# too close to planar to read a handedness from, so the existing tag is kept rather
# than resolved by numerical noise. RDKit declines at a similar threshold: the only
# three atoms where the rule below disagreed with RDKit across the whole fixture set
# were Cp ring carbons at |t| = 0.011-0.016.
_PLANARITY_TOL = 0.05


def _tag_from_geometry(mol, conf, idx):
    """Read a tetrahedral tag straight off the coordinates.

    Same convention as RDKit's ``assignChiralTypesFrom3D``: take the first three
    neighbours **in the atom's own bond order** and sign the triple product of the
    vectors to them; the implicit fourth substituent (a hydrogen, or a lone pair)
    sits opposite. Positive triple product is ``CHI_TETRAHEDRAL_CCW``.

    Validated against RDKit on every tagged atom of the project's 48 XYZ fixtures:
    255 agreements, 3 disagreements, all three at |unit triple product| < 0.016 --
    i.e. planar centres neither method should be reading, which ``_PLANARITY_TOL``
    now excludes.

    Needed because RDKit itself declines a whole class this encoder produces: a
    degree-3 atom with no hydrogen. After the metal bond is cut, every phosphine or
    amine donor looks exactly like that, and ``AssignAtomChiralTagsFromStructure``
    returns ``CHI_UNSPECIFIED`` for it -- which would leave the order-dependent
    copied tag in place on precisely the donor atoms this notation cares most about.
    """
    atom = mol.GetAtomWithIdx(idx)
    nbrs = [n.GetIdx() for n in atom.GetNeighbors()][:3]
    if len(nbrs) < 3:
        return None
    origin = np.asarray(conf.GetAtomPosition(idx))
    vecs = []
    for n in nbrs:
        v = np.asarray(conf.GetAtomPosition(n)) - origin
        norm = float(np.linalg.norm(v))
        if norm < 1e-6:
            return None
        vecs.append(v / norm)
    t = float(np.dot(np.cross(vecs[0], vecs[1]), vecs[2]))
    if abs(t) < _PLANARITY_TOL:
        return None
    return Chem.ChiralType.CHI_TETRAHEDRAL_CCW if t > 0 else Chem.ChiralType.CHI_TETRAHEDRAL_CW


def restamp_fragment_chirality(frag_rw, frag_to_parent, parent_conf):
    """Re-derive tetrahedral tags on a rebuilt fragment from the parent geometry.

    Parameters
    ----------
    frag_rw:
        The ``RWMol`` under construction. Modified in place.
    frag_to_parent:
        ``{fragment atom index: parent mol atom index}``. Must cover every
        fragment atom; a missing entry aborts the restamp (fail-safe, no change).
    parent_conf:
        The parent mol's 3D ``Conformer``. Its positions are indexed by parent
        atom index.

    Returns:
    -------
    int
        Number of tags actually flipped. ``0`` means the rebuild happened to
        preserve every parity for this ordering -- which is exactly what makes
        the defect intermittent.
    """
    prior = {}
    for atom in frag_rw.GetAtoms():
        tag = atom.GetChiralTag()
        if tag in _TETRAHEDRAL:
            prior[atom.GetIdx()] = tag
    if not prior:
        return 0

    n = frag_rw.GetNumAtoms()
    if any(i not in frag_to_parent for i in range(n)):
        return 0

    try:
        # Work on a copy: a perception failure must never corrupt the fragment
        # that is about to be serialized.
        work = Chem.RWMol(frag_rw)
        work.RemoveAllConformers()
        conf = Chem.Conformer(n)
        conf.Set3D(True)
        for i in range(n):
            conf.SetAtomPosition(i, parent_conf.GetAtomPosition(frag_to_parent[i]))
        work.AddConformer(conf, assignId=True)

        probe = work.GetMol()
        probe.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(probe)
        # replaceExistingTags=True: the copied tags are precisely what we distrust.
        Chem.AssignAtomChiralTagsFromStructure(probe)
        probe_conf = probe.GetConformer()
    except Exception:  # noqa: BLE001 - fail-safe: leave the fragment untouched
        return 0

    flipped = 0
    for idx, old_tag in prior.items():
        new_tag = probe.GetAtomWithIdx(idx).GetChiralTag()
        if new_tag not in _TETRAHEDRAL:
            # RDKit declined. Overwhelmingly this is a degree-3 heteroatom donor
            # with a lone pair rather than a hydrogen, which it skips by design --
            # read the handedness off the coordinates directly instead.
            try:
                new_tag = _tag_from_geometry(probe, probe_conf, idx)
            except Exception:  # noqa: BLE001
                new_tag = None
        # Still nothing means the centre is too planar to read. Keep the tag we have
        # rather than dropping stereochemistry -- an absent tag is a worse answer
        # than an uncertain one.
        if new_tag not in _TETRAHEDRAL:
            continue
        if new_tag != old_tag:
            frag_rw.GetAtomWithIdx(idx).SetChiralTag(new_tag)
            flipped += 1
    return flipped
