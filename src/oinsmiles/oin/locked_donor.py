"""Metal-locked donor stereochemistry — v0.4.5 Lane 6, injectivity blind spot P3.

Why this exists
===============
A secondary amine that is stereogenic **only because the metal occupies its fourth
position** loses its configuration completely. ``POJJOP`` (square-planar Pd, whose sole
stereocentre is the Pd-bound 2° amine) and its mirror encoded byte-identically::

    [Pd_SPL].Cc1ccc([NH]{0}Cc2ccccn{1}2)cc1.[Cl]{2}.[Cl]{3}

The amine appears as a bare ``[NH]{0}``. That is total *encoder* blindness, not merely
key blindness — the notation had no way to say which enantiomer it was.

Two rules cause it, and both are correct in general:

* ``ChiralityRecoveryUtility.recover``'s Zone-A clear (``core/chirality.py``) drops the
  chiral tag of any P/N whose ``total_degree < 4`` once the metal has been stripped, and
* RDKit clears trivalent-nitrogen chirality unconditionally, because a free amine
  really does invert.

Neither is true for a **metal-locked** donor: the metal bond is what makes the centre
stereogenic and it is also what stops the inversion. So the exemption is carved out for
exactly that case and no wider — the general rule is what keeps roughly a fifth of the
corpus from acquiring spurious sp3 handedness (see ``core/translator.py``'s
``with_stereo=False`` note).

The same argument applies verbatim to a metal-locked **phosphorus** donor whose tag was
already cleared by the Zone-A rule before Lane 8's ``restamp_fragment_chirality`` could
correct it — ``stable_stereo`` only fixes tags that already exist, so making them exist
is this module's job. One mechanism covers N and P; only the *survivability* of the
resulting tag differs (see `Notation`_).

Recovering the configuration
============================
Straight from the input 3D. The signed volume of the tetrahedron on the donor's four
neighbours inverts between enantiomers (``POJJOP`` −9.4 vs mirror +9.4, per
``tools/injectivity/config_oracle.py::bound_amine_centers``).

*Expressing* it reproducibly is the hard part, and the obvious route is a trap.
``bound_amine_centers`` orders neighbours with ``CanonicalRankAtoms(breakTies=True)``,
measured **not** invariant under input renumbering (2–11 distinct rank vectors over 20
renumberings) — fine as an oracle on one fixed molecule, useless as an emission
mechanism. Nor is "use the fragment's own bond order and let the canonical SMILES writer
normalise it" (Lane 8's argument) universally right here. Measured on rdkit 2025.09.3,
by holding one molecule fixed and feeding the writer two different bond orders at the
stereocentre:

============================  ==========================================================
fragment-atom class           does the emitted ``@``/``@@`` depend on bond order?
============================  ==========================================================
3 or 4 real bonds             **yes** — the tag is a parity relative to bond order
2 real bonds + implicit H     **no** — the same tag emits the same character either way
============================  ==========================================================

So the two classes need different frames, and a metal-locked donor lands in *both*: a
PR₃ phosphine donor keeps three heavy neighbours after the metal is cut, a secondary
amine or phosphine keeps only two.

* **3+ real bonds** — parity over the first three neighbours *in fragment bond order*,
  which is RDKit's own ``assignChiralTypesFrom3D`` rule. The numbering-dependence
  cancels: reorder the bonds and both this parity and the writer's interpretation of it
  flip together. Same frame as ``oin/stable_stereo.py``.
* **2 real bonds + implicit H** — the tag is an absolute label, so it must come from a
  frame that does not depend on numbering at all. Drop the metal (the neighbour that is
  *absent* from the fragment, standing in for the lone pair) and order the remaining
  three by ``CanonicalRankAtoms(breakTies=False)`` symmetry-class rank — which the
  stereogenicity gate below has already guaranteed to be three distinct values, so no
  index-dependent tiebreak is ever consulted. Measured on POJJOP: identical signed
  volume (+9.405) across 6 random renumberings, negated (−9.405) for the z-mirror.

Both frames then apply the *same* formula and sign convention as RDKit
(``t = (u₀ × u₁) · u₂`` over unit vectors to the three reference neighbours, positive →
``CCW``, the fourth substituent opposite), so the only difference between them is where
the ordering comes from.

**Honest limit on absolute sense.** For the 2-bond class RDKit has no answer to match:
``AssignAtomChiralTagsFromStructure`` returns ``CHI_UNSPECIFIED`` for a 3-coordinate
nitrogen even with all three neighbours explicit (measured), which is exactly why this
module has to do the work. So "which of ``@``/``@@`` is R" is an OIN convention fixed by
the rule above, not an inherited RDKit or CIP one. That is enough for injectivity — the
descriptor is reproducible, orientation-invariant, and inverts under reflection — and the
tests pin all three properties so the convention cannot drift silently.

Stereogenicity gate
===================
A metal-bound ammine (M,H,H,H) or primary amine (M,H,H,R) has symmetry-equivalent
hydrogens and is **not** a stereocentre. Emitting a descriptor for those is the
over-sensitivity failure the axial lane already had to fix once, so all four neighbours
must sit in four distinct symmetry classes before anything is stamped.

.. _Notation:

Notation
========
The standard SMILES tag: ``[N@]``/``[N@@]`` and ``[P@]``/``[P@@]``. Measured on rdkit
2025.09.3:

* the *writer* emits ``[N@@H]`` for a 2-heavy-neighbour nitrogen without complaint;
* ``AssignStereochemistry(cleanIt=True)`` clears nitrogen unconditionally, so a
  sanitising ``MolFromSmiles`` drops it too (``sanitize=False`` keeps it);
* phosphorus with three distinct substituents survives both, and is correctly cleared
  when two substituents are identical.

So the tag is emittable as long as it is applied **after** the last sanitising step,
which is why this module stamps a *property* during the fragment rebuild and
``recover()`` converts it to a tag as its final action. The honest limitation: for
nitrogen the descriptor is encoder-authoritative but not RDKit-round-trippable — a
sanitising re-parse of the OIN drops it, so the generator cannot rebuild it and the
round trip reports a loud false negative instead of a silent collision. Same trade the
axial lever makes, and the reason this lever is default-OFF.

An out-of-band token (``|amine:0+|``, following ``|ax:±|``) would survive re-parse, but
it needs its own atom identity, sign convention and multi-token ordering — three fresh
opportunities to re-run the Y2 "sorted by sign, therefore reflection-invariant" mistake —
and it would give phosphorus a *second* representation on top of the ``[P@]`` the
existing Zone-A lone-pair path already emits. Rejected for that, not for effort.

Gating
======
Nothing here runs unless ``OIN_EMIT_LOCKED_DONOR`` is set. With it unset no property is
stamped, ``recover()``'s restore loop finds nothing, and output is byte-identical.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from rdkit import Chem

from ..core.constants import TRANSITION_METALS_NUM
from .levers import lever_enabled as _lever_enabled

#: Environment lever. Default OFF -> byte-identical output.
ENV_LEVER = "OIN_EMIT_LOCKED_DONOR"

#: Atom property carrying the recovered configuration across the fragment rebuild and
#: the sanitising steps that would otherwise clear a nitrogen tag. ``"CW"`` or ``"CCW"``,
#: matching the ``Chem.ChiralType`` names.
LOCKED_TAG_PROP = "_OIN_LOCKED_STEREO"

#: Donors this applies to: nitrogen and phosphorus.
_DONOR_ATOMIC_NUMS: frozenset[int] = frozenset({7, 15})

#: Below this |triple product| of the three *unit* reference vectors the centre is too
#: close to planar to read a handedness from, so nothing is stamped rather than resolving
#: it by numerical noise. Same value and rationale as
#: ``oin/stable_stereo.py::_PLANARITY_TOL``.
_PLANARITY_TOL = 0.05

_TAG_BY_NAME = {
    "CW": Chem.ChiralType.CHI_TETRAHEDRAL_CW,
    "CCW": Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
}


class LockedDonorPlan(NamedTuple):
    """Per-molecule eligibility, computed once so the ranking is not repeated per fragment.

    Attributes:
    ----------
    indices:
        Parent atom indices of the metal-locked, genuinely stereogenic N/P donors.
    sym_ranks:
        ``CanonicalRankAtoms(parent, breakTies=False)`` — symmetry-class labels, used both
        for the gate and for the canonical frame of the 2-real-bond case.
    """

    indices: tuple[int, ...]
    sym_ranks: tuple[int, ...]

    def __bool__(self) -> bool:
        """A plan is truthy when it found at least one eligible donor."""
        return bool(self.indices)


_EMPTY_PLAN = LockedDonorPlan((), ())


def lever_enabled() -> bool:
    """True when the ``OIN_EMIT_LOCKED_DONOR`` lever is enabled.

    Routed through the central registry so ``OIN_EMIT_LOCKED_DONOR=0`` DISABLES it. The
    bare ``bool(os.environ.get(...))`` this replaced did the opposite -- "0" is a non-empty
    string, so opting out the obvious way switched the descriptor on.
    """
    return _lever_enabled(ENV_LEVER)


def plan_locked_donors(parent: Chem.Mol) -> LockedDonorPlan:
    """Find *parent*'s metal-locked, genuinely stereogenic N/P donors.

    Four conditions, all narrow on purpose:

    1. the atom is N or P and is **not** aromatic (an aromatic nitrogen cannot be a
       tetrahedral stereocentre, and a pyridine donor must never be touched);
    2. it is bonded to **exactly one** transition metal — a bridging donor's
       configuration is not a property of one centre;
    3. it has **exactly four** neighbours in the metal-present mol, i.e. the metal plus
       three, which is what "the metal locks the fourth position" means; and
    4. those four neighbours occupy **four distinct symmetry classes**, so a metal-bound
       ammine (M,H,H,H) or primary amine (M,H,H,R) is excluded.

    Returns an empty plan on any perception failure, so a molecule whose canonical
    ranking cannot be computed degrades to today's behaviour rather than raising.
    """
    try:
        # Symmetry-aware ranks: equivalent atoms SHARE a value. Compared only for equality
        # (the gate) and for ordering a set already known to be pairwise distinct (the
        # frame), so no index-dependent tiebreak is ever involved -- unlike
        # breakTies=True, which is not renumbering-invariant.
        sym_ranks = tuple(Chem.CanonicalRankAtoms(parent, breakTies=False))
    except Exception:  # noqa: BLE001 - guarded: no eligibility rather than a crash
        return _EMPTY_PLAN

    found: list[int] = []
    for atom in parent.GetAtoms():
        if atom.GetAtomicNum() not in _DONOR_ATOMIC_NUMS:
            continue
        if atom.GetIsAromatic():
            continue
        nbrs = list(atom.GetNeighbors())
        if len(nbrs) != 4:
            continue
        if sum(1 for n in nbrs if n.GetAtomicNum() in TRANSITION_METALS_NUM) != 1:
            continue
        if len({sym_ranks[n.GetIdx()] for n in nbrs}) != 4:
            continue
        found.append(atom.GetIdx())
    return LockedDonorPlan(tuple(found), sym_ranks)


def _reference_neighbours(
    parent: Chem.Mol,
    frag: Chem.Mol,
    frag_idx: int,
    frag_to_parent: dict[int, int],
    parent_idx: int,
    sym_ranks: tuple[int, ...],
) -> list[int] | None:
    """Three parent atom indices whose triple product defines the configuration.

    Which frame is used depends on how many real bonds the *fragment* atom kept — see
    the module docstring's table, which was measured rather than assumed.
    """
    frag_atom = frag.GetAtomWithIdx(frag_idx)
    # Bond order, not GetNeighbors() order: this is the ordering an RDKit chiral tag on a
    # 3-or-more-bond atom is a parity against.
    bonded: list[int] = []
    for bond in frag_atom.GetBonds():
        mapped = frag_to_parent.get(bond.GetOtherAtomIdx(frag_idx))
        if mapped is None:
            return None
        bonded.append(mapped)

    if len(bonded) >= 3:
        return bonded[:3]

    # 2 real bonds + implicit H: the tag is an absolute label, so the frame must not
    # depend on numbering. Drop the metal -- the neighbour the fragment is missing, which
    # is standing in for the lone pair, and therefore the one RDKit's own rule would put
    # opposite the other three -- and order the rest by symmetry-class rank.
    rest = [
        n.GetIdx()
        for n in parent.GetAtomWithIdx(parent_idx).GetNeighbors()
        if n.GetAtomicNum() not in TRANSITION_METALS_NUM
    ]
    if len(rest) != 3:
        return None
    return sorted(rest, key=lambda i: sym_ranks[i])


def _tag_name_from_geometry(
    parent: Chem.Mol,
    parent_conf: Chem.Conformer,
    frag: Chem.Mol,
    frag_idx: int,
    frag_to_parent: dict[int, int],
    parent_idx: int,
    sym_ranks: tuple[int, ...],
) -> str | None:
    """Read the fragment atom's tag name off the parent's coordinates.

    RDKit's own convention (``assignChiralTypesFrom3D``): sign the triple product of the
    unit vectors to three reference neighbours; positive is ``CCW``, and the fourth
    substituent sits opposite. Only the *source of the ordering* differs between the two
    frames -- see :func:`_reference_neighbours`.
    """
    ref = _reference_neighbours(parent, frag, frag_idx, frag_to_parent, parent_idx, sym_ranks)
    if ref is None:
        return None

    origin = np.asarray(parent_conf.GetAtomPosition(parent_idx))
    vecs = []
    for p in ref:
        v = np.asarray(parent_conf.GetAtomPosition(p)) - origin
        norm = float(np.linalg.norm(v))
        if norm < 1e-6:
            return None
        vecs.append(v / norm)

    t = float(np.dot(np.cross(vecs[0], vecs[1]), vecs[2]))
    if abs(t) < _PLANARITY_TOL:
        return None
    return "CCW" if t > 0 else "CW"


def stamp_locked_donor_stereo(
    frag: Chem.Mol,
    parent: Chem.Mol,
    parent_conf: Chem.Conformer,
    old_to_new: dict[int, int],
    plan: LockedDonorPlan | None = None,
) -> int:
    """Stamp :data:`LOCKED_TAG_PROP` on *frag*'s metal-locked stereogenic N/P donors.

    Parameters
    ----------
    frag:
        The rebuilt ligand fragment, modified in place. Only atom *properties* are set —
        no tag, bond or atom is touched — so on its own this cannot change the emitted
        string; ``ChiralityRecoveryUtility.recover`` does that.
    parent:
        The metal-present mol the fragment was rebuilt from.
    parent_conf:
        *parent*'s conformer. Must be the **pristine** input conformer, not one that has
        been through ``_align_to_pai`` — that alignment can reflect the coordinates,
        which would invert the recovered sign.
    old_to_new:
        ``{parent atom index: fragment atom index}`` from the rebuild.
    plan:
        Pre-computed :func:`plan_locked_donors` result, so the parent is ranked once per
        molecule rather than once per fragment. Computed here when omitted.

    Returns:
    -------
    int
        Number of atoms stamped.
    """
    if plan is None:
        plan = plan_locked_donors(parent)
    if not plan:
        return 0

    frag_to_parent = {v: k for k, v in old_to_new.items()}
    stamped = 0
    for parent_idx in plan.indices:
        frag_idx = old_to_new.get(parent_idx)
        if frag_idx is None:
            continue
        try:
            name = _tag_name_from_geometry(
                parent, parent_conf, frag, frag_idx, frag_to_parent, parent_idx, plan.sym_ranks
            )
        except Exception:  # noqa: BLE001 - guarded: stamp nothing rather than crash
            name = None
        if name is None:
            continue
        frag.GetAtomWithIdx(frag_idx).SetProp(LOCKED_TAG_PROP, name)
        stamped += 1
    return stamped


def restore_locked_donor_tags(rw: Chem.RWMol) -> int:
    """Turn :data:`LOCKED_TAG_PROP` back into a chiral tag, in place.

    Called as the **last** stereo action of ``ChiralityRecoveryUtility.recover`` because
    everything before it — ``AssignStereochemistry(cleanIt=True)`` in particular — clears
    a trivalent nitrogen's tag by design.

    Only atoms whose tag is currently ``CHI_UNSPECIFIED`` are written. A donor that
    already carries a tag has been handled by an existing, measured path (the Zone-A
    lone-pair branch for P, or the >=4-neighbour verify-and-flip) and must not be
    second-guessed here: this lane's remit is the donors that end up with *no* tag.

    Returns:
    -------
    int
        Number of tags set.
    """
    written = 0
    for atom in rw.GetAtoms():
        if not atom.HasProp(LOCKED_TAG_PROP):
            continue
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        tag = _TAG_BY_NAME.get(atom.GetProp(LOCKED_TAG_PROP))
        if tag is None:
            continue
        atom.SetChiralTag(tag)
        written += 1
    return written
