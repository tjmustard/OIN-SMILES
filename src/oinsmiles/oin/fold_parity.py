"""Reflection-parity veto for the within-fragment donor fold (v0.4.12 Lane 1).

WHAT v0.4.11 PROVED, AND WHY THE OBVIOUS FIX IS WRONG
=====================================================
``OIN_CANONICAL_DONOR_FOLD`` exchanges donors that are (a) in one
``CanonicalRankAtoms(breakTies=False)`` symmetry class of their fragment and (b) the same
vertex colour. It works: 393 molecules move ``key_equal/slot_renumber -> byte_exact``, none
in any other direction, and the comparison key moves on none of them. It also **collapses
enantiomers in 221 of those same 393 gains**.

The reason is one sentence, and it is the transferable part:

    A fragment's automorphism says nothing about the PARITY of the vertex permutation it
    induces on the coordination sphere.

The tempting repair -- require each swap to be a proper rotation of the polyhedron, reusing
``canonical_slots.derive_rotation_group``'s ``det > 0`` test -- is **wrong before it is
built**. A donor swap is a transposition that fixes every other vertex, and for a rank-3
polyhedron such a permutation does not preserve the Gram matrix at all, so that test rejects
*every* swap and the fold degenerates to the identity. The fold is not justified by a symmetry
of the polyhedron; it is justified by a symmetry of the LIGAND realized as a proper rotation of
the whole complex, which is not a property of the coordination graph and cannot be read off it.

SO THE PARITY IS TESTED WHERE THE ANSWER LIVES: THE GEOMETRY
===========================================================
This module lifts ``tools/mirror_audit_donor_fold.py``'s own verdict into the encoder as a
per-molecule veto. That tool is the instrument that caught the defect, so agreeing with it by
construction -- rather than by a second, independently-reasoned predicate -- is the point.

    S_rot        = canonicalize(inline,        fold OFF)
    S_fold       = canonicalize(inline,        fold ON)
    S_rot_m      = canonicalize(inline_mirror, fold OFF)
    S_fold_m     = canonicalize(inline_mirror, fold ON)

    veto   <=>   (S_rot != S_rot_m)  and  (S_fold == S_fold_m)
    emit    =    S_rot if veto else S_fold

**The left conjunct is load-bearing.** Without it every achiral molecule is vetoed -- its
mirror encodes identically with the fold ON *and* OFF, and that is not this lever's doing. The
same conjunct is what correctly declines to blame this lever for a metal-centred Delta/Lambda
pair, whose descriptor the shipped encoder already folds because ``OIN_EMIT_METAL_CONFIG`` is
held off. That is a pre-existing gap and it belongs to v0.4.16. This is byte-for-byte the audit
tool's ``REGRESSION <=> OFF_distinct and not ON_distinct`` implication.

COST: ONE EXTRA ENCODE, ON THE ~7.5% THAT ACTUALLY FOLD
=======================================================
``canonicalize_oin_slots`` is a pure string operation, so three of the four strings above are
free once the structure has been encoded once. Only ``inline_mirror`` needs perception, and it
is computed **only when the fold actually changes the string** -- when ``S_rot == S_fold`` there
is nothing to veto and the mirror is never built.

PRESENTATION-INVARIANCE
=======================
The veto predicate is a function of the STRUCTURE, not of the incoming atom numbering: both
arms of the comparison are re-derived from coordinates, and ``canonicalize_oin_slots`` is
already presentation-invariant on each. So two presentations of one complex reach the same
verdict and emit the same string. This is asserted, not assumed --
``tests/unit/test_fold_parity.py`` renumbers the input and requires the verdict to hold.
"""

from __future__ import annotations

import logging
import os
import tempfile
import threading

logger = logging.getLogger(__name__)

FOLD_LEVER = "OIN_CANONICAL_DONOR_FOLD"
VETO_LEVER = "OIN_FOLD_PARITY_VETO"

#: Re-entrancy guard. The veto encodes a mirrored structure, and that encode runs through the
#: very post-pass that calls this module. Without the guard it would recurse until the stack
#: gave out. Thread-local rather than a module global because the harness runs molecules
#: concurrently and one molecule's mirror pass must not disarm another molecule's veto.
_state = threading.local()


def _in_mirror_pass() -> bool:
    return getattr(_state, "in_mirror", False)


def last_outcome() -> str | None:
    """Why the most recent :func:`resolve` call decided as it did, on this thread.

    Five outcomes, and keeping them DISTINCT is the load-bearing part. Three of them
    (``declined_*``) emit the rotation-only labeling because the veto had no usable evidence;
    ``vetoed_collapse`` emits it because the evidence said the fold destroys a mirror pair.
    Both look identical in the output string, and a mirror audit cannot tell them apart --
    which is exactly how the first implementation of this module scored 18/18 "successes"
    while its self-check was declining on every single molecule and the veto proper never ran.
    The three fixture tests passed too, because declining also separates a mirror pair.

    Consumers: ``tools/fold_transition_sim.py`` and the lane's write-up. Never a control flow
    input -- it is an observation, not a decision.
    """
    return getattr(_state, "outcome", None)


def _note(outcome: str) -> None:
    _state.outcome = outcome
    logger.debug("fold parity: %s", outcome)


class _mirror_pass:
    """Mark the dynamic extent of the veto's own encode, so it cannot re-enter the veto."""

    def __enter__(self):
        _state.in_mirror = True
        return self

    def __exit__(self, *exc):
        _state.in_mirror = False
        return False


class _forced_lever:
    """Force one lever on/off for a dynamic extent, restoring the prior value exactly.

    Restores by DELETING the key when it was previously unset, rather than writing ``"0"``.
    Writing ``"0"`` would be a different state from unset for any code that still reads the
    environment directly, and this project has already been bitten twice by ``"0"`` not
    meaning what it looks like.
    """

    def __init__(self, name: str, on: bool):
        self.name = name
        self.on = on

    def __enter__(self):
        self.prior = os.environ.get(self.name)
        os.environ[self.name] = "1" if self.on else "0"
        return self

    def __exit__(self, *exc):
        if self.prior is None:
            os.environ.pop(self.name, None)
        else:
            os.environ[self.name] = self.prior
        return False


def _mirror_coords(coords):
    """Reflect through the xy-plane (negate z) -- the same mirror the audit tool applies.

    Any reflection works: they differ from one another by a proper rotation, which the
    encoder's own orientation canonicalization removes. Matching the audit tool's choice keeps
    the in-encoder verdict and the offline verdict comparable line by line.
    """
    import numpy as np

    out = np.array(coords, dtype=float, copy=True)
    out[:, 2] *= -1.0
    return out


def _atom_coord_pairs(tmc_mol, xyz_coords):
    """``[(symbol, coordinate)]`` for ``tmc_mol``, or ``None`` if they cannot be aligned.

    ⚠ ``tmc_mol``'s atom order is NOT the coordinate order. Perception rebuilds the molecule
    and stamps each atom with ``__origIdx``, its index into ``xyz_coords``; ``get_oin_string``
    itself reads coordinates as ``xyz_coords[orig_i]`` for exactly this reason. Zipping
    ``GetAtoms()`` against ``xyz_coords`` positionally instead produces a chemically different
    molecule that still encodes cleanly -- on ``BIWDIV_comp_0`` it yielded ``[Co_TBP]`` with
    invented bond orders in place of ``[Co_OCT]``. That is the failure mode this whole module
    exists to prevent, so it is handled here rather than trusted.

    Falls back to positional order only when NO atom carries the property (a mol that did not
    come from ``_get_tmc_mol_impl`` at all, e.g. the generator's contract mol on the scored
    path). ``resolve``'s self-check below is what makes that fallback safe.
    """
    atoms = list(tmc_mol.GetAtoms())
    if len(atoms) != len(xyz_coords):
        logger.debug("fold parity: atom/coordinate count mismatch, declining to fold")
        return None

    if all(a.HasProp("__origIdx") for a in atoms):
        try:
            idx = [a.GetIntProp("__origIdx") for a in atoms]
        except Exception:  # noqa: BLE001 -- a malformed property is a decline, not a crash
            return None
        if sorted(idx) != list(range(len(atoms))):
            logger.debug("fold parity: __origIdx is not a permutation, declining to fold")
            return None
    else:
        idx = list(range(len(atoms)))

    return [(a.GetSymbol(), xyz_coords[i]) for a, i in zip(atoms, idx)]


def _encode_pairs(pairs, mirror, fold=False):
    """Encode a reconstructed structure through the PUBLIC converter. ``None`` on any failure.

    Going through ``XYZToSMILES().convert()`` rather than reaching into ``get_oin_string``
    with a hand-mirrored mol is deliberate. ``tmc_mol`` carries chiral tags already perceived
    from the ORIGINAL coordinates; negating z without redoing perception would leave those
    tags describing the un-mirrored structure, and the veto would compare a mirror against
    itself and silently under-fire. Re-perceiving from written coordinates is the only way to
    be sure every stereo-bearing feature is re-read. It is also what makes this agree with
    ``tools/mirror_audit_donor_fold.py``, which encodes a written mirrored file the same way.

    A failure returns ``None``. A structure whose mirror cannot be encoded yields no evidence
    that the fold is safe, and the caller treats "no evidence" as "do not fold" -- the
    conservative direction, matching how an unknown geometry degrades to identity folding.
    """
    fh = None
    try:
        fh = tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False)
        fh.write(f"{len(pairs)}\n\n")
        for sym, xyz in pairs:
            z = -float(xyz[2]) if mirror else float(xyz[2])
            fh.write(f"{sym:<3} {float(xyz[0]):>16.10f} {float(xyz[1]):>16.10f} {z:>16.10f}\n")
        fh.close()

        from ..core.translator import XYZToSMILES

        # ⚠ The fold MUST be forced for this encode rather than inherited. Inheriting the
        # ambient (ON) setting returns the mirror's FOLDED string, which makes `s_rot_m`
        # identical to `s_fold_m`; the left conjunct then degenerates to "did the fold fire?",
        # which is true by construction at this point, and the achiral guard is silently
        # disarmed. Measured: with that bug the self-check declined on 18 of 18 movers and the
        # three fixture tests still PASSED -- declining to fold also separates a mirror pair,
        # so a dead veto and a working one are indistinguishable from the outside.
        with _mirror_pass(), _forced_lever(FOLD_LEVER, fold):
            return XYZToSMILES().convert(fh.name)
    except Exception:  # noqa: BLE001 -- an unencodable structure is a "do not fold", not a crash
        logger.debug("fold parity: reconstruction encode failed", exc_info=True)
        return None
    finally:
        if fh is not None:
            try:
                os.unlink(fh.name)
            except OSError:
                pass


def resolve(inline_oin: str, tmc_mol, xyz_coords) -> str:
    """Apply the slot post-pass, vetoing the donor fold when it would collapse a mirror pair.

    Returns the string to emit. With ``OIN_FOLD_PARITY_VETO`` off this is exactly
    ``canonicalize_oin_slots(inline_oin)`` -- one call, no mirror, byte-identical to v0.4.11.
    """
    from .canonical_slots import canonicalize_oin_slots
    from .levers import lever_enabled

    # Off, or already inside the veto's own mirror encode -> the plain post-pass. The second
    # condition is what stops infinite recursion, and it must be checked before anything
    # expensive happens.
    if not lever_enabled(VETO_LEVER) or _in_mirror_pass():
        return canonicalize_oin_slots(inline_oin)

    # The veto only has an opinion about the donor fold. With the fold itself off there is no
    # widened candidate set to police, so this is a no-op by construction.
    if not lever_enabled(FOLD_LEVER):
        return canonicalize_oin_slots(inline_oin)

    _note("pending")

    with _forced_lever(FOLD_LEVER, False):
        s_rot = canonicalize_oin_slots(inline_oin)
    s_fold = canonicalize_oin_slots(inline_oin)

    if s_rot == s_fold:
        # The fold did not fire on this molecule. Nothing to veto, and -- the reason this
        # check is here rather than after the mirror -- nothing to pay for either.
        _note("fold_inactive")
        return s_fold

    if xyz_coords is None or tmc_mol is None:
        _note("declined_no_conformer")
        return s_rot

    pairs = _atom_coord_pairs(tmc_mol, xyz_coords)
    if pairs is None:
        _note("declined_no_pairs")
        return s_rot

    # SELF-CHECK, and it is not ceremony. The veto's verdict is only as good as its
    # reconstruction of the structure, and a mis-paired symbol/coordinate list does not raise
    # -- it encodes cleanly as a DIFFERENT molecule, whose mirror comparison is then pure
    # fiction. This project has shipped that exact shape before: v0.4.7's first attachment
    # check was a silent no-op whose complete A/B run "reported what a genuine null result
    # looks like". So: re-encode the reconstruction UNMIRRORED and require it to reproduce the
    # labeling we are actually deciding about. If it does not, we have no instrument here and
    # we decline to fold rather than guess.
    # Both reconstruction encodes run with the fold FORCED OFF, so they yield the
    # rotation-only labelings the four-way comparison is defined on.
    self_oin = _encode_pairs(pairs, mirror=False, fold=False)
    if self_oin is None:
        _note("declined_no_self_encode")
        return s_rot
    with _forced_lever(FOLD_LEVER, False):
        if canonicalize_oin_slots(self_oin) != s_rot:
            _note("declined_reconstruction_drift")
            return s_rot

    mirror_oin = _encode_pairs(pairs, mirror=True, fold=False)
    if mirror_oin is None:
        # No evidence the fold is safe here. Decline to fold.
        _note("declined_no_mirror")
        return s_rot

    with _forced_lever(FOLD_LEVER, False):
        s_rot_m = canonicalize_oin_slots(mirror_oin)
    s_fold_m = canonicalize_oin_slots(mirror_oin)

    if s_rot == s_rot_m:
        # The shipped encoder ALREADY folds this pair -- an achiral molecule, or a
        # metal-centred Delta/Lambda whose descriptor is held off (v0.4.16's gap, not this
        # lever's). Nothing here is this fold's doing, so it is allowed.
        _note("allowed_preexisting_fold")
        return s_fold
    if s_fold == s_fold_m:
        _note("vetoed_collapse")
        return s_rot
    _note("allowed_separation_survives")
    return s_fold
