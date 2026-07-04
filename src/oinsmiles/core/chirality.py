"""Chirality encoding utilities for the XYZ->OIN pipeline.

Two classes handle the full lifecycle of P/N stereocenter encoding:
  - CIPAssigner: pre-fragmentation -- reads 3D conformer, stores
    ``_OIN_CIPCode`` (metal-present convention) on P/N atoms, and
    ``_OIN_CIPCode_LP`` (fragment-local / lone-pair convention) on eligible
    Zone-A P atoms (Stereo Phase 4).
  - ChiralityRecoveryUtility: post-fragmentation -- verifies/corrects chiral
    tags after OINSanitizer strips the metal context. Zone-A P atoms carrying
    ``_OIN_CIPCode_LP`` are verified-and-flipped (never cleared); everything
    else keeps the original degree-keyed behaviour.

See ``spec/compiled/SuperPRD_StereoPhase4_ZoneA_P.md`` Sec.5.1 for the
dummy-metal-copy recipe and the lone-pair CIP convention this module
implements.
"""

from __future__ import annotations

import warnings

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from .constants import TRANSITION_METALS_NUM

#: RDKit atomic numbers for nitrogen (N) and phosphorus (P).
_PN_ATOMIC_NUMS: frozenset[int] = frozenset({7, 15})

#: RDKit atomic number for phosphorus (P) alone -- Zone-A LP-CIP is P-only
#: (Stereo Phase 4 Sec.3.2: nitrogen is explicitly out of scope, trivalent
#: [N@] is cleared by RDKit as non-stereogenic amine inversion).
_P_ATOMIC_NUM: int = 15

#: Atom property name for the fragment-local (lone-pair convention) CIP
#: label, distinct from the existing metal-present '_OIN_CIPCode'.
_LP_CIP_PROP: str = "_OIN_CIPCode_LP"


class OINStereoWarning(UserWarning):
    """All Phase-4 Zone-A P stereo diagnostics.

    The message ALWAYS embeds the atom index (defeats Python's warning
    dedup-by-message-text, so repeated warnings for different atoms are not
    silently collapsed into one).
    """


def _eligible_zone_a_p(mol: Chem.Mol) -> list[int]:
    """Return indices of P atoms bonded to exactly one metal atom.

    "Metal" is decided by the single existing predicate source
    (``TRANSITION_METALS_NUM``) -- never duplicated (TD-005).

    A P atom bonded to >=2 metal atoms (bridging phosphide) is excluded and
    triggers an ``OINStereoWarning`` naming the atom index (B7); this
    degrades that atom to today's clearing behaviour in ``recover()``.
    """
    eligible: list[int] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _P_ATOMIC_NUM:
            continue

        metal_neighbor_count = sum(
            1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() in TRANSITION_METALS_NUM
        )

        if metal_neighbor_count == 1:
            eligible.append(atom.GetIdx())
        elif metal_neighbor_count >= 2:
            warnings.warn(
                OINStereoWarning(
                    f"atom {atom.GetIdx()}: P bonded to {metal_neighbor_count} "
                    "metal atoms (bridging phosphide) -- Zone-A lone-pair CIP "
                    "not computed; degrades to clearing behaviour."
                ),
                stacklevel=2,
            )
    return eligible


def _build_dummy_metal_copy(mol: Chem.Mol, p_idx: int) -> Chem.Mol | None:
    """Build the dummy-metal copy of *mol* for P atom *p_idx*'s lone-pair CIP label.

    Normative recipe (SuperPRD Sec.5.1, B4): the single metal bonded to this
    P is turned into a Z=0 dummy atom (formal charge zeroed, isotope
    cleared); every OTHER metal-ligand bond -- including bonds belonging to
    any OTHER metal atom entirely -- is removed, leaving only the
    dummy--P bond. A Z=0 dummy is lowest CIP priority, same as a lone pair,
    and it sits exactly where the metal (and therefore the lone pair) was,
    so this reproduces the fragment-local CIP sense while the 3D conformer
    is still intact.

    Implementation note (found empirically, not spelled out by the
    normative recipe text): the metal--P bond in ``get_tmc_mol()`` output is
    a ``DATIVE`` bond. A dative bond counts toward the ACCEPTOR's valence
    but not the DONOR's (P is the donor here per ``xyz2mol.py``'s
    ``AddBond(ligand_atom, metal_atom, DATIVE)`` convention) -- so RDKit's
    from-3D chirality perception sees only 3 valence-counting substituents
    on P, tops up the "missing" 4th with a phantom implicit H, and silently
    refuses to assign a chiral tag (verified: ``AssignAtomChiralTagsFromStructure``
    returns ``CHI_UNSPECIFIED`` on the untouched dummy copy even though
    ``Chem.FindMolChiralCenters(..., includeUnassigned=True)`` correctly
    flags the atom as a stereocentre). The fix: turn the dummy--P bond into
    an ordinary valence-counting ``SINGLE`` bond and set ``NoImplicit(True)``
    on both the dummy atom and P, so RDKit sees exactly 4 real substituents
    and no phantom H. This is purely a graph/valence-bookkeeping device for
    this scratch copy -- it does not change the real mol or the OIN grammar.

    The entire body runs inside try/except: ANY failure (missing metal
    neighbour, bond-surgery error, ``Chem.SanitizeMol`` exception -- e.g. an
    eta-ligand aromatic ring stranded by bond removal) returns ``None`` and
    emits an ``OINStereoWarning`` naming *p_idx*. Callers degrade to
    store-nothing; ``convert()`` must never crash because of this helper
    (B4, RISK-7).
    """
    try:
        p_atom = mol.GetAtomWithIdx(p_idx)
        metal_idx = None
        for nbr in p_atom.GetNeighbors():
            if nbr.GetAtomicNum() in TRANSITION_METALS_NUM:
                metal_idx = nbr.GetIdx()
                break
        if metal_idx is None:
            raise ValueError(f"no metal neighbour found on P atom {p_idx}")

        rw = Chem.RWMol(mol)

        # Drop every metal-ligand bond except metal_idx--p_idx, for ALL
        # metal atoms present (not just metal_idx) -- normative recipe.
        bonds_to_remove: list[tuple[int, int]] = []
        for atom in rw.GetAtoms():
            if atom.GetAtomicNum() not in TRANSITION_METALS_NUM:
                continue
            for nbr in atom.GetNeighbors():
                if atom.GetIdx() == metal_idx and nbr.GetIdx() == p_idx:
                    continue  # keep the one bond we need
                pair = (atom.GetIdx(), nbr.GetIdx())
                if pair not in bonds_to_remove and pair[::-1] not in bonds_to_remove:
                    bonds_to_remove.append(pair)

        for a_idx, b_idx in bonds_to_remove:
            if rw.GetBondBetweenAtoms(a_idx, b_idx) is not None:
                rw.RemoveBond(a_idx, b_idx)

        dummy_atom = rw.GetAtomWithIdx(metal_idx)
        dummy_atom.SetAtomicNum(0)
        dummy_atom.SetFormalCharge(0)
        dummy_atom.SetIsotope(0)
        dummy_atom.SetNoImplicit(True)

        # Valence-bookkeeping fix (see docstring): make the dummy--P bond an
        # ordinary bond so it counts toward P's valence, and suppress the
        # phantom implicit H on P that would otherwise mask the 4th
        # substituent from chirality perception.
        dummy_bond = rw.GetBondBetweenAtoms(metal_idx, p_idx)
        dummy_bond.SetBondType(Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(p_idx).SetNoImplicit(True)

        dummy_mol = rw.GetMol()
        Chem.SanitizeMol(dummy_mol)
        return dummy_mol
    except Exception as exc:  # noqa: BLE001 - guarded degradation, B4
        warnings.warn(
            OINStereoWarning(
                f"atom {p_idx}: dummy-metal copy construction failed "
                f"({exc!r}) -- Zone-A lone-pair CIP not computed; degrades "
                "to today's clearing behaviour."
            ),
            stacklevel=2,
        )
        return None


def _lp_cip_label(dummy_mol: Chem.Mol, p_idx: int) -> str | None:
    """Compute the lone-pair-convention CIP label for P *p_idx* in *dummy_mol*.

    ``rdCIPLabeler`` is the ONLY label source for ``_OIN_CIPCode_LP`` (B5) --
    the legacy ``Chem.AssignStereochemistry`` ``_CIPCode`` is never used for
    this property (legacy may still run elsewhere for tag *perception*).
    Chiral tags must be (re)derived from the 3D conformer first because the
    dummy atom's identity differs from the original metal's, which can
    change the CIP priority ordering at this centre.

    Guarded: any exception, or no label assignable (e.g. P is not a CIP
    stereocentre -- two identical substituents), returns ``None``. This is
    the expected, silent outcome for symmetric P and must NOT warn (BDPP/
    BDNN negative controls rely on this).
    """
    try:
        Chem.AssignAtomChiralTagsFromStructure(dummy_mol)
        Chem.AssignStereochemistry(dummy_mol, cleanIt=True, force=True)
        rdCIPLabeler.AssignCIPLabels(dummy_mol)
        atom = dummy_mol.GetAtomWithIdx(p_idx)
        return atom.GetPropsAsDict().get("_CIPCode")
    except Exception:  # noqa: BLE001 - guarded, store-nothing degradation
        return None


def _metal_present_cip_label(mol: Chem.Mol, p_idx: int) -> str | None:
    """Diagnostic-only helper: metal-present-convention CIP label for the HITL print/table.

    NOT part of the normative dummy-metal recipe, NOT stored on the atom,
    NEVER compared to ``_OIN_CIPCode_LP`` (B1).

    Discovered empirically while building ``_build_dummy_metal_copy``: the
    legacy ``_OIN_CIPCode`` set by ``assign_all()``'s existing (unchanged)
    block is silently ``None`` for a Zone-A P/N stereocentre, for the exact
    same reason the dummy-metal copy needed a valence-bookkeeping fix -- the
    metal--P bond is DATIVE, so RDKit's from-3D perception never sees a real
    4th substituent on P and never assigns a tag. This is a PRE-EXISTING gap
    in the legacy metal-present computation (masked until now because
    ``recover()`` always cleared Zone-A tags before anyone could notice the
    stored value was already empty), not something Stereo Phase 4 changes
    the meaning of. Fixing the general legacy computation is out of scope
    for this MiniPRD; this helper exists ONLY so the HITL candidate-output
    table (Task 13 / SuperPRD Sec.9) has a real metal-present label to
    cross-check the literature (R,R) configuration against, instead of
    printing an uninformative ``None``.

    Applies the same bond-type/NoImplicit fix as ``_build_dummy_metal_copy``
    but keeps the metal's real identity (no Z=0 swap). Guarded: returns
    ``None`` on any failure.
    """
    try:
        p_atom = mol.GetAtomWithIdx(p_idx)
        metal_idx = None
        for nbr in p_atom.GetNeighbors():
            if nbr.GetAtomicNum() in TRANSITION_METALS_NUM:
                metal_idx = nbr.GetIdx()
                break
        if metal_idx is None:
            return None

        rw = Chem.RWMol(mol)
        bonds_to_remove: list[tuple[int, int]] = []
        for atom in rw.GetAtoms():
            if atom.GetAtomicNum() not in TRANSITION_METALS_NUM:
                continue
            for nbr in atom.GetNeighbors():
                if atom.GetIdx() == metal_idx and nbr.GetIdx() == p_idx:
                    continue
                pair = (atom.GetIdx(), nbr.GetIdx())
                if pair not in bonds_to_remove and pair[::-1] not in bonds_to_remove:
                    bonds_to_remove.append(pair)
        for a_idx, b_idx in bonds_to_remove:
            if rw.GetBondBetweenAtoms(a_idx, b_idx) is not None:
                rw.RemoveBond(a_idx, b_idx)

        bond = rw.GetBondBetweenAtoms(metal_idx, p_idx)
        bond.SetBondType(Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(metal_idx).SetNoImplicit(True)
        rw.GetAtomWithIdx(p_idx).SetNoImplicit(True)

        copy_mol = rw.GetMol()
        Chem.SanitizeMol(copy_mol)
        Chem.AssignAtomChiralTagsFromStructure(copy_mol)
        Chem.AssignStereochemistry(copy_mol, cleanIt=True, force=True)
        rdCIPLabeler.AssignCIPLabels(copy_mol)
        return copy_mol.GetAtomWithIdx(p_idx).GetPropsAsDict().get("_CIPCode")
    except Exception:  # noqa: BLE001 - diagnostic-only, never raises
        return None


def _attach_dummy_metal(mol: Chem.RWMol, p_idx: int) -> int:
    """Append a Z=0 wildcard as stereogenic P *p_idx*'s 4th substituent, in place.

    Generation-side counterpart to ``_build_dummy_metal_copy`` (Stereo Phase
    4 MiniPRD-C): that helper CONVERTS an EXISTING metal neighbour on an
    already-coordinated P (and raises when none exists) -- this helper
    ATTACHES a brand-new dummy to a metal-free fragment atom, before ETKDG
    ever embeds it, so the P's ``[P@]``/``[P@@]`` tag genuinely controls a
    4-coordinate tetrahedron instead of describing an under-determined
    3-coordinate one. New code, not a revival of the deleted
    ``PseudoAtomStrategy``, and never calls ``_build_dummy_metal_copy``.

    Mirrors ``_build_dummy_metal_copy``'s valence-bookkeeping fix (SINGLE
    bond + ``NoImplicit(True)`` on both atoms, see :func:`_build_dummy_metal_copy`
    docstring) so RDKit perceives 4 real substituents on P and honours the
    chiral tag through embedding instead of topping up a phantom implicit H.

    Mutates *mol* in place; the dummy is appended as the highest atom index.

    Returns:
    -------
    int
        The new dummy atom's index.
    """
    dummy_idx = mol.AddAtom(Chem.Atom(0))
    mol.AddBond(p_idx, dummy_idx, Chem.BondType.SINGLE)
    mol.GetAtomWithIdx(dummy_idx).SetNoImplicit(True)
    mol.GetAtomWithIdx(p_idx).SetNoImplicit(True)
    return dummy_idx


class CIPAssigner:
    """Assigns 3D-derived CIP codes and chiral tags to P/N stereocenters.

    Must be called on the **full** TMC mol (pre-fragmentation) so that the
    metal coordination context is available for correct CIP assignment.

    Preconditions
    -------------
    * *mol* must have a valid embedded 3D conformer (as returned by
      ``get_tmc_mol()``).
    * ``Chem.SanitizeMol(mol)`` must succeed -- any sanitisation exception
      propagates to the caller.

    Postconditions
    --------------
    * Each P/N atom with an assignable CIP code gains the atom property
      ``_OIN_CIPCode`` ('R' or 'S') -- the metal-present convention.
    * Chiral tags (``CHI_TETRAHEDRAL_CW`` / ``CHI_TETRAHEDRAL_CCW``) are set
      from the 3D conformer geometry.
    * Each eligible Zone-A P atom (bonded to exactly one metal) gains the
      atom property ``_OIN_CIPCode_LP`` ('R' or 'S') -- the fragment-local
      (lone-pair) convention, computed via the dummy-metal copy. Any
      pre-existing ``_OIN_CIPCode_LP`` is cleared first (idempotence: a
      second call never leaves a stale tag from a previous run).
    """

    def assign_all(self, mol: Chem.Mol, diagnostics: bool = True) -> Chem.Mol:
        """Assign stereo tags and store ``_OIN_CIPCode`` / ``_OIN_CIPCode_LP``.

        Parameters
        ----------
        mol:
            RDKit ``Mol`` with an embedded 3D conformer. Must not be ``None``.
        diagnostics:
            When ``True`` (default), re-run the lone-pair CIP computation on
            the same dummy-metal copy and warn (``OINStereoWarning``) on a
            same-convention conflict, and print the metal-present label for
            HITL visibility. Set ``False`` to skip this for batch users
            concerned about ``rdCIPLabeler`` runtime (RISK-10).

        Returns:
        -------
        Chem.Mol
            The same mol object with updated atom properties and chiral tags.

        Raises:
        ------
        ValueError
            If *mol* is ``None``.
        Chem.SanitizeMol exceptions
            Any sanitisation failure propagates unchanged so the caller can
            decide how to handle it.
        """
        if mol is None:
            raise ValueError("mol must not be None")

        # Hard precondition — exception propagates to caller.
        Chem.SanitizeMol(mol)

        # MUST precede AssignStereochemistry: sets CHI_TETRAHEDRAL_CW/CCW from
        # the 3D conformer geometry.  Without this call, AssignStereochemistry
        # computes _CIPCode but does NOT write @/@@ into subsequent SMILES.
        Chem.AssignAtomChiralTagsFromStructure(mol)

        # Compute _CIPCode from the chiral tags set above.
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in _PN_ATOMIC_NUMS:
                cip = atom.GetPropsAsDict().get("_CIPCode")
                if cip:
                    atom.SetProp("_OIN_CIPCode", cip)

        # --- Zone-A P lone-pair CIP (Stereo Phase 4) ---
        # Idempotence (B9): clear any stale _OIN_CIPCode_LP before
        # recomputing, so repeated calls never leave a tag stale.
        for atom in mol.GetAtoms():
            if atom.HasProp(_LP_CIP_PROP):
                atom.ClearProp(_LP_CIP_PROP)

        for p_idx in _eligible_zone_a_p(mol):
            dummy_mol = _build_dummy_metal_copy(mol, p_idx)
            if dummy_mol is None:
                continue  # _build_dummy_metal_copy already warned

            label = _lp_cip_label(dummy_mol, p_idx)
            if label is not None:
                mol.GetAtomWithIdx(p_idx).SetProp(_LP_CIP_PROP, label)

                # Seed a SPECIFIED chiral tag on the REAL atom so it survives
                # fragmentation (AddAtom copies ChiralTag) and gives
                # recover()'s verify-and-flip something to act on. The exact
                # handedness copied here does not need to already match the
                # fragment's own bond order -- recover() recomputes CIP on
                # the fragment mol and flips to match the stored label
                # regardless (same pattern as the pre-existing >=4-neighbour
                # branch). It only needs to be non-UNSPECIFIED to begin with.
                dummy_tag = dummy_mol.GetAtomWithIdx(p_idx).GetChiralTag()
                if dummy_tag == Chem.ChiralType.CHI_UNSPECIFIED:
                    dummy_tag = Chem.ChiralType.CHI_TETRAHEDRAL_CW
                mol.GetAtomWithIdx(p_idx).SetChiralTag(dummy_tag)

            if diagnostics:
                # Same-convention cross-check ONLY (B1): re-run on the SAME
                # dummy-metal copy. Never compare against the metal-present
                # convention -- that comparison is legitimately divergent by
                # construction (Sec.5.1), not a bug.
                recheck = _lp_cip_label(dummy_mol, p_idx)
                if recheck != label:
                    warnings.warn(
                        OINStereoWarning(
                            f"atom {p_idx}: Zone-A P lone-pair CIP cross-check "
                            f"mismatch on the dummy-metal copy (stored="
                            f"{label!r}, recheck={recheck!r})."
                        ),
                        stacklevel=2,
                    )
                # Metal-present label: print-only for HITL visibility, NEVER
                # compared against the lone-pair label (B1). The stored
                # legacy _OIN_CIPCode is frequently None here (see
                # _metal_present_cip_label docstring: dative-bond valence
                # perception gap, pre-existing, out of scope to fix
                # generally) -- fall back to the diagnostic recompute so the
                # HITL table has a real label to check (R,R) against.
                metal_present_cip = mol.GetAtomWithIdx(p_idx).GetPropsAsDict().get("_OIN_CIPCode")
                if metal_present_cip is None:
                    metal_present_cip = _metal_present_cip_label(mol, p_idx)
                print(
                    f"[Zone-A P diagnostic] atom {p_idx}: "
                    f"lone-pair convention={label!r}  "
                    f"metal-present convention={metal_present_cip!r} "
                    "(informational only, never compared)"
                )

        return mol


class ChiralityRecoveryUtility:
    """Verifies and corrects @/@@ on P/N atoms in post-fragmentation fragment mols.

    After ``OINSanitizer.generate_robust_smiles()`` the metal neighbour is
    absent. Branch order (Stereo Phase 4, B9):

    1. Zone-A P atoms carrying ``_OIN_CIPCode_LP`` (set by ``CIPAssigner``):
       KEEP the chiral tag and verify-and-flip it against a fresh
       ``rdCIPLabeler`` recompute on this fragment's own geometry-derived
       tags -- regardless of total degree. Multi-P fragments recompute CIP
       after each flip (bounded at 2 full passes, B6) because one P's flip
       can change another's CIP priority ranking.
    2. Everything else falls through to the existing degree-keyed behaviour
       unchanged: Zone A P/N atoms without the property (including every
       symmetric P and all N -- N is out of scope, Sec.3.2) have their
       chiral tag cleared (@/@@ is undefined without the full coordination
       sphere or a lone-pair tag); atoms with >=4 neighbours use the stored
       (metal-present) ``_OIN_CIPCode`` to flip/keep the tag; 4-neighbour
       atoms with no ``_OIN_CIPCode`` (non-standard valence) have the stray
       tag cleared.
    """

    def recover(self, mol: Chem.Mol) -> Chem.Mol:
        """Re-apply correct chiral tags to P/N atoms in a fragment mol.

        Parameters
        ----------
        mol:
            Fragment ``Mol`` returned by ``OINSanitizer.generate_robust_smiles()``.

        Returns:
        -------
        Chem.Mol
            Updated mol with corrected P/N chiral tags.
        """
        if mol is None:
            return mol

        rw = Chem.RWMol(mol)

        # --- Zone-A P: verify-and-flip keyed on _OIN_CIPCode_LP ---
        # Runs BEFORE the degree-keyed branches below: this is the whole
        # point of Stereo Phase 4 -- a Zone-A P (total_degree < 4 in the
        # fragment, metal excluded) with a lone-pair tag must be KEPT, not
        # cleared by the old unconditional Zone-A clear.
        tagged_p_indices = [
            atom.GetIdx()
            for atom in rw.GetAtoms()
            if atom.GetAtomicNum() == _P_ATOMIC_NUM and atom.HasProp(_LP_CIP_PROP)
        ]

        if tagged_p_indices:
            # Bounded fixed-point (B6): flipping one P's tag can change
            # another P's CIP priority ranking in the SAME fragment, so we
            # recompute rdCIPLabeler before evaluating each subsequent P,
            # across at most 2 full passes over the tagged-P set.
            for _pass in range(2):
                any_changed = False
                for p_idx in tagged_p_indices:
                    Chem.AssignStereochemistry(rw, cleanIt=True, force=True)
                    try:
                        rdCIPLabeler.AssignCIPLabels(rw)
                    except Exception:  # noqa: BLE001 - guarded recompute
                        continue

                    atom = rw.GetAtomWithIdx(p_idx)
                    ctag = atom.GetChiralTag()
                    if ctag == Chem.ChiralType.CHI_UNSPECIFIED:
                        continue  # nothing to flip

                    stored_lp = atom.GetPropsAsDict().get(_LP_CIP_PROP)
                    current_cip = atom.GetPropsAsDict().get("_CIPCode")

                    if stored_lp and current_cip != stored_lp:
                        if ctag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                        else:
                            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                        any_changed = True

                if not any_changed:
                    break

        # --- Existing degree-keyed branches (unchanged behaviour) ---
        # Compute CIP codes from existing chiral tags (no 3D conformer needed).
        Chem.AssignStereochemistry(rw, cleanIt=True, force=True)

        for atom in rw.GetAtoms():
            if atom.GetAtomicNum() not in _PN_ATOMIC_NUMS:
                continue
            if atom.GetIdx() in tagged_p_indices:
                continue  # already handled by the Zone-A lone-pair branch above

            stored_cip: str | None = atom.GetPropsAsDict().get("_OIN_CIPCode")
            current_cip: str | None = atom.GetPropsAsDict().get("_CIPCode")
            total_deg: int = atom.GetTotalDegree()

            if total_deg < 4:
                # Zone A atom without a lone-pair tag (symmetric P, or any
                # N -- N is out of scope, Sec.3.2). Metal removed, 3
                # neighbours remain; @/@@ is undefined without the full
                # coordination sphere or a lone-pair tag.
                atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

            elif stored_cip and total_deg >= 4:
                if current_cip != stored_cip:
                    # Flip the tag so SMILES output matches the stored CIP code.
                    ctag = atom.GetChiralTag()
                    if ctag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                        atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                    elif ctag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                        atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                    # CHI_UNSPECIFIED: cannot flip — leave as-is

            else:
                # 4 neighbours but no _OIN_CIPCode — non-standard valence.
                # Neutral fallback: clear the stray chiral tag.
                atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

        return rw.GetMol()
