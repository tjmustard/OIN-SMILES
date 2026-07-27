"""Canonical ligand-body emission — v0.4.5 Lane 1.

WHY
===
Every OIN ligand body is already serialized with ``MolToSmiles(..., canonical=True)``, so
RDKit's *serializer* is canonical. The defect is one level up: **the graph handed to the
serializer is not**. Perception from 3D distances makes two non-unique choices --
``get_UA_pairs`` uses ``networkx.max_weight_matching`` (maximum matchings are not unique,
``utils/perception_core.py``) and ``AC2BO`` returns, in its own words, "an arbitrary
resonance form". Two geometries of one molecule therefore become two different graph
objects and ``canonical=True`` faithfully turns each into a different string:
``O=c1[cH]{2>}[cH]{2}c(=O)o1`` from the crystal, ``O=C1[CH]{2>}=[CH]{2}C(=O)O1`` from the
regenerated structure.

``oin.compare`` already repairs this **at compare time** -- that is why the round-trip
*key* is canonical while the emitted *string* is not. This module moves that machinery
upstream so the encoder emits the canonical form directly. It reuses
``compare._parse_fragment`` and ``compare.canonical_fragment_body`` rather than copying
them, so the encoder and the key cannot drift apart.

THE HARD PART: SLOT IDENTITY MUST SURVIVE THE REPARSE
=====================================================
``{n}`` markers are placed from ``_smilesAtomOutputOrder``, and a reparse changes that
order. Re-deriving the position by substructure match is *not* safe -- on a near-symmetric
ligand ``GetSubstructMatch`` can return a wrong automorphism, putting the marker on a CH
instead of the deprotonated X-type carbon. So slot identity is **carried through** the
reparse instead of re-derived: each donor gets an atom map number, the map numbers survive
``MolToSmiles`` -> ``MolFromSmiles``, the donors are recovered by map number, the map
numbers are cleared, and only then is the final canonical SMILES written and its output
order read.

Every failure returns ``None``, and the caller keeps the un-reparsed body for that
**whole** fragment. There is no partial application: a misplaced ``{n}`` silently corrupts
the coordination sphere, which is far worse than the notation drift this closes. Three
guards are load-bearing:

1. **Composition.** A map number forces brackets, and a bracket *changes* implicit-H
   semantics -- ``n`` in a five-ring means one implicit H, ``[n:1]`` means none. So the
   reparse is only accepted when the element counts, total H, and total formal charge come
   out unchanged. This is the phantom-hydrogen class that ``generate_robust_smiles``
   stage 1b exists to prevent, and it must not be re-opened here.
2. **Donor identity.** The atom recovered by map number ``k+1`` must match the donor it
   claims to be in element and heavy degree, and no map number may be lost or duplicated.
3. **Idempotence.** The emitted body must be its own canonical form, which is what makes
   the result a genuine canonical representative rather than one step of a walk. Measured
   over 6062 distinct bodies from the capstone corpus, one pass already reaches the fixed
   point for 6056; the other 6 **oscillate with period two** because RDKit flips ``@``/
   ``@@`` on adamantane-cage carbons across a parse/write cycle (a degenerate cage
   "stereocentre" that is not one). Those bail and keep the un-reparsed body. That RDKit
   instability also affects ``compare.canonical_fragment_body`` itself, so it is a
   pre-existing hazard in the comparison key, not something this module introduces.
"""

from rdkit import Chem

from .compare import _parse_fragment, canonical_fragment_body
from .hydrogen import h_faithful_smiles

__all__ = ["canonical_body", "canonical_body_emit"]

#: Above this many donors the chelate-ring probe is not worth building; mirrors the
#: identical guard in ``compare._chelate_locked_fragment_key``.
_MAX_PROBE_DONORS = 10

#: Reparse passes allowed while looking for the fixed point. One pass suffices for
#: 6056/6062 corpus bodies; the headroom costs one extra ``MolToSmiles`` and catches the
#: rare fragment that needs two. Anything still moving after this is an oscillator.
_MAX_REPARSE_PASSES = 4


def canonical_body(body_smiles: str) -> str:
    """Canonical form of one **slot-stripped** ligand body.

    This is the function Lane 2 should call for its vertex colors, and the predicate the
    Lane 1 acceptance test asserts against: an emitted body ``b`` is canonical iff
    ``canonical_body(b) == b``. It is ``compare.canonical_fragment_body`` under a name the
    encoder side can import without implying the comparison layer -- one implementation,
    three consumers, so they cannot drift apart.
    """
    return canonical_fragment_body(body_smiles)


def _output_order(mol):
    """``_smilesAtomOutputOrder`` from the last ``MolToSmiles``, as a list, or ``None``."""
    if mol is None or not mol.HasProp("_smilesAtomOutputOrder"):
        return None
    try:
        raw = mol.GetProp("_smilesAtomOutputOrder")
        return [int(x) for x in raw.strip("[]").rstrip(",").split(",") if x != ""]
    except Exception:
        return None


def _composition(mol):
    """``(element counts, total H, total formal charge)`` -- what the reparse must preserve.

    Radical electrons are deliberately excluded: RDKit re-derives them from the valence on
    every parse (a carbene written ``[CH2]`` comes back with two radical electrons even
    when the source mol carried none), so comparing them would reject every carbene while
    telling us nothing about whether the molecule changed. A radical that silently became
    a hydrogen *does* show up, in the total-H term.
    """
    counts: dict[int, int] = {}
    n_h = 0
    charge = 0
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        counts[z] = counts.get(z, 0) + 1
        n_h += atom.GetTotalNumHs()
        charge += atom.GetFormalCharge()
    return tuple(sorted(counts.items())), n_h, charge


def _n_aromatic(mol) -> int:
    return sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())


def _clear_chelate_locked_stereo(mol, donor_indices) -> None:
    r"""Clear E/Z on double bonds held rigid by a ring that closes through the metal.

    Mutates ``mol`` in place. Same semantics as
    ``compare._chelate_locked_fragment_key``: build a probe in which a dummy metal is
    bonded to every donor, take the rings through that metal, and clear ``BondStereo``
    plus every incident ``BondDir`` on the DOUBLE bonds those rings lock. A double bond in
    *no* metal ring -- a pendant, genuinely flippable alkene or imine -- is untouched, so
    real diastereomers still separate.

    Without this, a salicylaldiminate / beta-diketiminate / eta-alkene bond acquires a
    ``/`` or ``\`` whose sign depends on SMILES traversal, and the two round-trip
    directions disagree about a bond that was never freely E/Z.
    """
    if not donor_indices or len(donor_indices) > _MAX_PROBE_DONORS:
        return
    probe = Chem.RWMol(mol)
    dummy = probe.AddAtom(Chem.Atom(26))  # Fe, matching the compare-layer probe
    for d in donor_indices:
        if probe.GetBondBetweenAtoms(d, dummy) is None:
            probe.AddBond(d, dummy, Chem.BondType.SINGLE)
    p = probe.GetMol()
    Chem.FastFindRings(p)

    locked: set[int] = set()
    for ring in Chem.GetSymmSSSR(p):
        ring_set = set(ring)
        if dummy in ring_set:
            locked |= ring_set
    if not locked:
        return

    for bond in p.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a == dummy or b == dummy or a not in locked or b not in locked:
            continue
        rbond = mol.GetBondBetweenAtoms(a, b)
        if rbond is None or rbond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        rbond.SetStereo(Chem.BondStereo.STEREONONE)
        for end in (rbond.GetBeginAtom(), rbond.GetEndAtom()):
            for nb in end.GetBonds():
                if nb.GetBondDir() != Chem.BondDir.NONE:
                    nb.SetBondDir(Chem.BondDir.NONE)


def _heavy_degree(atom) -> int:
    return sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() != 1)


def _reparse_once(mol, donors):
    """One map-number-carrying reparse pass.

    Returns ``(smiles, new_donor_indices, reparsed_mol)`` -- the reparsed mol still carries
    ``_smilesAtomOutputOrder`` for ``smiles`` and has no atom map numbers left, so it can
    be fed straight back in for the next pass -- or ``None`` if any guard trips.
    """
    work = Chem.RWMol(mol)
    for k, d in enumerate(donors):
        work.GetAtomWithIdx(d).SetAtomMapNum(k + 1)
    # h_faithful_smiles, NOT Chem.MolToSmiles. This write is the one that decides the emitted
    # body, so a plain MolToSmiles here silently discards the OIN_H_FAITHFUL repair applied
    # upstream at perception_tmc.py -- the caller overwrites its H-faithful string with
    # whatever this
    # function returns. That is why the two levers were mutually exclusive.
    #
    # Safe by construction: h_faithful_smiles returns plain MolToSmiles output verbatim when
    # OIN_H_FAITHFUL is off, so with the shipped default this line is byte-identical. When the
    # lever IS on it can only change a string whose atom count was already measurably wrong.
    #
    # It also makes the composition guard below reachable in the useful direction: a fragment
    # whose bare 0-H symbol used to re-read one hydrogen heavier failed `_composition`, so the
    # whole reparse was discarded and the fragment kept a non-canonical body. Repairing the H
    # first lets those fragments canonicalize instead of being rejected.
    mapped = h_faithful_smiles(work, isomericSmiles=True, canonical=True)
    if not mapped:
        return None

    reparsed = _parse_fragment(mapped)
    if reparsed is None:
        return None

    # Guard 1: the reparse must not have changed the molecule. A map number forces a
    # bracket, and a bracket redefines the implicit-H count of the atom it wraps.
    if _composition(mol) != _composition(reparsed):
        return None

    # Guard 1b: the reparse must not DE-AROMATIZE. A fragment's aromatic flags can carry
    # information the free fragment no longer supports on its own: a metallo-porphyrin's
    # macrocycle is perceived aromatic in the complex, but once the metal is stripped its
    # four N's have no hydrogen and RDKit's default model calls the free base non-aromatic,
    # so a naive reparse would emit `C1=C2C=CC(=N2)...` where the encoder had `c1c2nc(...`.
    # That is fine for the COMPARISON key (both sides lose it identically) but at emit time
    # it is a loss of chemical fidelity -- exactly what tests/unit/test_aromatic_reencode.py
    # pins. Measured over the 713 distinct ligand bodies in the capstone `rdkit_canonical`
    # population, this guard rejects NOTHING (0 of 713 lose aromaticity) and all 500 rows
    # still converge, so it costs no folding power at all.
    if _n_aromatic(reparsed) < _n_aromatic(mol):
        return None

    by_map: dict[int, int] = {}
    for atom in reparsed.GetAtoms():
        num = atom.GetAtomMapNum()
        if not num:
            continue
        if num in by_map:
            return None  # duplicated label: the donors are no longer distinguishable
        by_map[num] = atom.GetIdx()
    if len(by_map) != len(donors):
        return None

    # Guard 2: every marker must still land on the donor it was stamped on.
    new_donors: list[int] = []
    for k, d in enumerate(donors):
        j = by_map.get(k + 1)
        if j is None:
            return None
        old_atom, new_atom = mol.GetAtomWithIdx(d), reparsed.GetAtomWithIdx(j)
        if old_atom.GetAtomicNum() != new_atom.GetAtomicNum():
            return None
        if _heavy_degree(old_atom) != _heavy_degree(new_atom):
            return None
        new_donors.append(j)

    _clear_chelate_locked_stereo(reparsed, new_donors)

    # ⚠ DO NOT restore the metal-locked donor tag here. Tried in v0.4.6 and MEASURED WRONG.
    #
    # The idea was obvious and the correspondence is available: Guard 2 above has already proven
    # donors[k] <-> new_donors[k] (same element, same heavy degree), so copying the chiral tag
    # from `mol` onto `reparsed` looks safe and does make P3 emit under OIN_CANONICAL_BODY.
    # POJJOP passes. It is still wrong.
    #
    # Setting a chiral tag AFTER the sanitize adds a stereocentre the canonical ranker did not
    # know about, which changes the canonical WRITE ORDER -- and `@`/`@@` is a parity relative to
    # that order, not an absolute label. On RIFGUJ_comp_2 (three Cu-bound amines on one
    # cyclohexane) the three ring-CARBON tags then flip between a structure and its mirror.
    # The geometry says they must not: AssignStereochemistryFrom3D + rdCIPLabeler give those
    # carbons lowercase `s` -- pseudo-asymmetric, i.e. a RELATIVE (all-cis) descriptor -- and
    # they read `s` identically for the structure AND its reflection. So the restoration
    # silently rewrote carbon stereochemistry that should have been untouched.
    #
    # tests/unit/test_locked_donor.py::TestMultiCentreDescriptor::test_flips_under_reflection and
    # ::test_three_locked_amines_all_invert_together are what caught it; the single-centre POJJOP
    # fixture could not, which is the Y2 lesson intact.
    #
    # A correct fix has to make the reparse preserve the tag WITHOUT perturbing the ranking --
    # e.g. keep the donor bracketed through the sanitize, or re-derive parity from the parent
    # geometry after the write order is fixed. Until then OIN_EMIT_LOCKED_DONOR stays
    # incompatible with OIN_CANONICAL_BODY, as recorded in oin/levers.py::_HELD_OFF.

    # Clear BEFORE the final emit: the body must be a pure canonical SMILES of the
    # reparsed graph, with no `:k` residue for the inline handler to trip over.
    for atom in reparsed.GetAtoms():
        atom.SetAtomMapNum(0)

    # h_faithful_smiles, NOT Chem.MolToSmiles: THIS is the write whose output becomes the
    # emitted body, so a plain MolToSmiles here discards the OIN_H_FAITHFUL repair no matter
    # what the intermediate pass did. Byte-identical when that lever is off.
    smiles = h_faithful_smiles(reparsed, isomericSmiles=True, canonical=True)
    if not smiles:
        return None
    return smiles, new_donors, reparsed


def canonical_body_emit(mol, donor_indices):
    """Canonical body SMILES for one ligand fragment, with the donors' SMILES positions.

    Args:
        mol: the sanitized fragment mol the encoder is about to serialize (the output of
            ``OINSanitizer.generate_robust_smiles`` after chirality recovery and
            ``SetDoubleBondNeighborDirections``). Not mutated.
        donor_indices: atom indices **in ``mol``'s own index space** that bond the metal,
            already routed through ``canonical_donor_representative`` /
            ``canonical_eta_set_representative`` by the caller. That choice is preserved:
            the map number is stamped on exactly the atom the caller named, so the
            canonical-representative remap is carried through, never undone.

    Returns:
        ``(smiles, {donor_index_in_mol: position_in_smiles}, reparsed_mol)``, or ``None``
        if anything at all fails -- in which case the caller must keep its existing body
        for the whole fragment. ``smiles`` is a fixed point of :func:`canonical_body`.
    """
    try:
        n_atoms = mol.GetNumAtoms()
        donors: list[int] = []
        for d in donor_indices:
            if not (0 <= d < n_atoms):
                return None
            if d not in donors:
                donors.append(d)

        # A map number already on the input would be indistinguishable from ours after the
        # reparse, so refuse rather than risk mislabelling a donor.
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum():
                return None

        # Iterate to the fixed point of the reparse, so the emitted body satisfies
        # canonical_body(body) == body rather than merely being one step along a walk.
        cur_mol, cur_donors = mol, donors
        last_smiles = None
        result = None
        for _ in range(_MAX_REPARSE_PASSES):
            out = _reparse_once(cur_mol, cur_donors)
            if out is None:
                return None
            smiles, new_donors, new_mol = out
            result = out
            if smiles == last_smiles:
                break
            last_smiles = smiles
            cur_mol, cur_donors = new_mol, new_donors
        else:
            return None  # oscillating (RDKit adamantyl @/@@ flip) -- keep the input body

        smiles, new_donors, reparsed = result
        order = _output_order(reparsed)
        if order is None:
            return None
        position = {frag_idx: pos for pos, frag_idx in enumerate(order)}

        remap: dict[int, int] = {}
        for d, j in zip(donors, new_donors):
            if j not in position:
                return None
            remap[d] = position[j]

        return smiles, remap, reparsed
    except Exception:
        return None
