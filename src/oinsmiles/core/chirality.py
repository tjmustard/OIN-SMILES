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

import logging
import warnings

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ..utils.aromaticity import sanitize_allowing_boron_cage
from .constants import TRANSITION_METALS_NUM

logger = logging.getLogger(__name__)

#: RDKit atomic numbers for nitrogen (N) and phosphorus (P).
_PN_ATOMIC_NUMS: frozenset[int] = frozenset({7, 15})

#: RDKit atomic number for phosphorus (P) alone -- Zone-A LP-CIP is P-only
#: (Stereo Phase 4 Sec.3.2: nitrogen is explicitly out of scope, trivalent
#: [N@] is cleared by RDKit as non-stereogenic amine inversion).
_P_ATOMIC_NUM: int = 15

#: Atom property name for the fragment-local (lone-pair convention) CIP
#: label, distinct from the existing metal-present '_OIN_CIPCode'.
_LP_CIP_PROP: str = "_OIN_CIPCode_LP"

#: Atom property name for the metal-free (fragment-convention) rdCIPLabeler CIP
#: label of a specified sp3 C/Si/S stereocentre, stamped by ``build_contract_mol``
#: so ``recover()`` can re-orient it on the metal-free fragment. A metal-/eta-adjacent
#: centre's CIP label flips between the metal-present contract mol and the metal-free
#: fragment, so the generator's metal-present flip loop mis-orients it.
_SP3_CIP_PROP: str = "_OIN_CIPCode_SP3"


#: Periodic table handle for valence-deficit capping (see _fill_open_valence_with_h).
_PERIODIC_TABLE = Chem.GetPeriodicTable()


def _fill_open_valence_with_h(mol: Chem.Mol) -> None:
    """Cap valence-deficient (metal-stripped donor) atoms with explicit H, in place.

    A donor atom in the fragment ``recover()`` sees has lost its dative metal bond, so
    it reads as an open-valence atom (donor ``O{0}`` -> bare ``[O]``). ``rdCIPLabeler``
    declines to assign CIP to a stereocentre bearing such a substituent, so a genuine
    sulfonimidoyl S(VI) (``JEKQAS``/``REPZUJ``/``ZORCOA``) gets no label and its
    re-orientation is skipped. Filling the deficit with H makes the neighbourhood
    labelable in the SAME convention ``_template_sp3_label`` reads. This runs only on a
    throwaway CIP-probe copy -- it is never emitted, and never alters the final SMILES
    (which caps donors with a formal charge, per ``OINDiscreteAligner``). Aromatic atoms
    are left untouched (an aromatic pyridine-N donor must not become ``[nH]``).
    """
    for a in mol.GetAtoms():
        if a.GetIsAromatic():
            continue
        default_val = _PERIODIC_TABLE.GetDefaultValence(a.GetAtomicNum())
        if default_val <= 0:
            continue
        deficit = default_val - a.GetTotalValence()
        if deficit > 0 and a.GetFormalCharge() == 0 and a.GetNumExplicitHs() == 0:
            a.SetNoImplicit(False)
            a.SetNumExplicitHs(deficit)
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:  # noqa: BLE001 - guarded
        pass


def _reparse_cip_label_once(smiles: str, probe: int, fill_deficit: bool) -> "str | None":
    """One rdCIPLabeler pass over *smiles*, optionally H-filling open valences first."""
    try:
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is None:
            return None
        try:
            Chem.SanitizeMol(m, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        except Exception:  # noqa: BLE001 - lenient perception
            m.UpdatePropertyCache(strict=False)
        if fill_deficit:
            _fill_open_valence_with_h(m)
        Chem.AssignStereochemistry(m, cleanIt=True, force=True)
        rdCIPLabeler.AssignCIPLabels(m)
        target = next((a for a in m.GetAtoms() if a.GetAtomMapNum() == probe), None)
        if target is not None and target.HasProp("_CIPCode"):
            return target.GetProp("_CIPCode")
        return None
    except Exception:  # noqa: BLE001 - guarded
        return None


def _reparse_aromatic_cip_label(mol: Chem.Mol, idx: int) -> "str | None":
    """Aromatic-preserving rdCIPLabeler label for atom *idx* on a fresh re-parse.

    rdCIPLabeler gives OPPOSITE R/S for a stereocentre bonded to an aromatic haptic
    ring (Cp, indenyl, fluorenyl) depending on the ring's aromatic/kekulized state,
    and a processed mol object can carry a corrupted aromatic state. Re-parsing from
    SMILES with kekulization skipped normalises that, so both the ``build_contract_mol``
    stamp (``_template_sp3_label``) and this ``recover()`` comparison read the label in
    ONE convention. An atom-map probe survives the SMILES atom re-ordering.

    Two passes: the first fills metal-stripped donor atoms' open valences with H so
    ``rdCIPLabeler`` can rank a centre bearing a radical donor substituent
    (sulfonimidoyl S ``JEKQAS``; carbene-/alkene-donor-adjacent carbon ``ORIHUU``/
    ``XILZID``) in the SAME convention ``_template_sp3_label`` reads. H-fill skips
    aromatic atoms, so an eta-Cp/arene-adjacent centre (``AHEBEV``/``BABWAD``/
    ``KAGXUM``) is unaffected. The second pass (no fill) is a fallback for the rare
    case where filling breaks perception. Returns None on any failure (caller degrades
    to leaving the tag as-is).
    """
    probe = 99
    try:
        tagged = Chem.Mol(mol)
        tagged.GetAtomWithIdx(idx).SetAtomMapNum(probe)
        smiles = Chem.MolToSmiles(tagged)
    except Exception:  # noqa: BLE001 - guarded
        return None
    for fill_deficit in (True, False):
        label = _reparse_cip_label_once(smiles, probe, fill_deficit)
        if label is not None:
            return label
    return None


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


def _clear_spurious_high_coordination_stereo(mol):
    """Drop octahedral / trigonal-bipyramidal tags on non-stereogenic main-group centres.

    ``AssignAtomChiralTagsFromStructure`` stamps a permutation tag on any atom whose
    3D neighbourhood is asymmetric, but a pentafluorosulfanyl (-SF5) sulfur is
    octahedral with **five identical terminal fluorines** and one substituent, so every
    arrangement is superimposable -- it is achiral. RDKit's stereo perception (legacy
    *and* modern ``FindPotentialStereo``) does not reduce the five equivalent F, so the
    tag survives into the OIN as a spurious ``[S@OH..]`` the 3D generator cannot
    reproduce (the embed places the F set symmetrically and the re-encode drops it),
    breaking the round trip. Clear the tag wherever a non-metal centre's stereochemistry
    rests on a set of identical terminal ligands. Transition metals are excluded -- their
    octahedral / TBP handedness is the coordination geometry and is carried elsewhere.
    """
    from collections import Counter

    from .constants import TRANSITION_METALS_NUM

    high_coord = (
        Chem.ChiralType.CHI_OCTAHEDRAL,
        Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL,
    )
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() not in high_coord:
            continue
        if atom.GetAtomicNum() in TRANSITION_METALS_NUM:
            continue
        terminal = Counter(n.GetAtomicNum() for n in atom.GetNeighbors() if n.GetDegree() == 1)
        # Non-stereogenic when only one ligand differs from a set of >=(degree-1)
        # identical terminal atoms (SF5: five terminal F on a degree-6 sulfur).
        if terminal and max(terminal.values()) >= atom.GetDegree() - 1:
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)


def _genuine_stereocentre_indices(mol) -> "frozenset[int] | None":
    """Atom indices RDKit's *modern* stereo perception treats as (potential) stereocentres.

    Computed on a sanitised copy (lenient property cache + ring perception) so a
    non-standard-valence fragment (C#O, charge-less Cp) does not break perception.
    ``includeUnassigned=True`` keeps a genuine centre whose handedness is not yet
    assigned. Returns ``None`` on any failure so the caller can degrade to leaving
    tags untouched rather than mistake a perception error for "no stereocentres".
    """
    try:
        probe = Chem.Mol(mol)
        probe.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(probe)
        centres = Chem.FindMolChiralCenters(
            probe, force=True, includeUnassigned=True, useLegacyImplementation=False
        )
        return frozenset(idx for idx, _ in centres)
    except Exception:  # noqa: BLE001 - guarded; None => leave tags as-is
        return None


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
        # Routed through the boron-cage-aware wrapper: with OIN_BORON_CAGE unset,
        # or on any mol without a deltahedral cage, this is exactly
        # Chem.SanitizeMol(mol) and propagates the same exceptions.
        sanitize_allowing_boron_cage(mol)

        # MUST precede AssignStereochemistry: sets CHI_TETRAHEDRAL_CW/CCW from
        # the 3D conformer geometry.  Without this call, AssignStereochemistry
        # computes _CIPCode but does NOT write @/@@ into subsequent SMILES.
        Chem.AssignAtomChiralTagsFromStructure(mol)

        # Compute _CIPCode from the chiral tags set above.
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        # Drop spurious -SF5-style high-coordination stereo (achiral, but tagged from
        # geometry) so the OIN does not carry an @ the generator cannot reproduce.
        _clear_spurious_high_coordination_stereo(mol)

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
                logger.debug(
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
                    atom = rw.GetAtomWithIdx(p_idx)
                    ctag = atom.GetChiralTag()
                    if ctag == Chem.ChiralType.CHI_UNSPECIFIED:
                        continue  # nothing to flip

                    stored_lp = atom.GetPropsAsDict().get(_LP_CIP_PROP)
                    # Read on the same aromatic-preserving re-parse the LP stamp uses,
                    # so a P donor bonded to an aromatic ring is compared in one
                    # convention (GUXPIA).
                    current_cip = _reparse_aromatic_cip_label(rw, p_idx)

                    if stored_lp and current_cip and current_cip != stored_lp:
                        if ctag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                        else:
                            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                        any_changed = True

                if not any_changed:
                    break

        # --- sp3 C/Si/S: verify-and-flip keyed on _OIN_CIPCode_SP3 ---
        # The generator's metal-present flip loop (build_contract_mol) orients a
        # backbone carbon correctly, but a metal-/eta-adjacent sp3 centre's CIP
        # label FLIPS between the metal-present contract mol and the metal-free
        # fragment (AHEBEV's benzylic C: R with the eta-arene bound to Cr, S once
        # the metal is stripped), so that loop mis-orients it. Re-orient here on
        # the metal-free fragment -- the same graph get_oin_string emits from --
        # against the rdCIPLabeler label taken from the (also metal-free) OIN
        # template. rdCIPLabeler, not legacy, for the same reason as Zone-A P.
        # No-op on the forward-encode path: only build_contract_mol stamps the prop.
        # N(7) is included for a genuine quaternary ammonium N+, which
        # build_contract_mol stamps with _SP3_CIP_PROP only when total_degree == 4
        # (POYJIX). A trivalent amine N never carries the prop, so it is untouched
        # here and still deferred by the degree-keyed branch / RDKit inversion clear.
        tagged_sp3_indices = [
            atom.GetIdx()
            for atom in rw.GetAtoms()
            if atom.GetAtomicNum() in (6, 7, 14, 16) and atom.HasProp(_SP3_CIP_PROP)
        ]
        if tagged_sp3_indices:
            for _pass in range(2):
                any_changed = False
                for idx in tagged_sp3_indices:
                    atom = rw.GetAtomWithIdx(idx)
                    ctag = atom.GetChiralTag()
                    if ctag == Chem.ChiralType.CHI_UNSPECIFIED:
                        continue
                    stored = atom.GetPropsAsDict().get(_SP3_CIP_PROP)
                    # Read the label on the SAME aromatic-preserving re-parse the
                    # stamp used, so a haptic-ring-adjacent centre is compared in one
                    # convention (BEPXEA broke when the two sides diverged).
                    current = _reparse_aromatic_cip_label(rw, idx)
                    if stored and current and current != stored:
                        atom.SetChiralTag(
                            Chem.ChiralType.CHI_TETRAHEDRAL_CCW
                            if ctag == Chem.ChiralType.CHI_TETRAHEDRAL_CW
                            else Chem.ChiralType.CHI_TETRAHEDRAL_CW
                        )
                        any_changed = True
                if not any_changed:
                    break

        # --- Spurious donor-S: clear a geometry-derived tag on a sulfur that is
        # NOT a genuine stereocentre once the metal is gone. build_contract_mol's
        # AssignStereochemistryFrom3D stamps a permutation tag ([S@SP3]/[S@SP1]/
        # [S@TB9H]) on a metal-donor thioether S from the metal-present geometry;
        # the tag survives fragmentation because legacy AssignStereochemistry(
        # cleanIt=True) does not scrub a pre-set permutation tag, so get_oin_string
        # emits an [S@..] the input crystal geometry never produced (BAZMOH, HUGSEI,
        # LUSKIV, YUMPIH, CIDDAU). recover() is on BOTH encode paths, so clearing
        # here makes the two sides symmetric. Gate on modern stereo perception so a
        # genuine chiral sulfonimidoyl S(VI) (degree 4) -- re-oriented by the
        # _SP3_CIP_PROP branch above, which also excludes it below -- is never masked.
        s_tagged = [
            atom.GetIdx()
            for atom in rw.GetAtoms()
            if atom.GetAtomicNum() == 16
            and atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            and not atom.HasProp(_SP3_CIP_PROP)
        ]
        if s_tagged:
            genuine = _genuine_stereocentre_indices(rw)
            if genuine is not None:
                for idx in s_tagged:
                    if idx not in genuine:
                        rw.GetAtomWithIdx(idx).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

        # --- Existing degree-keyed branches (unchanged behaviour) ---
        # Compute CIP codes from existing chiral tags (no 3D conformer needed).
        Chem.AssignStereochemistry(rw, cleanIt=True, force=True)

        for atom in rw.GetAtoms():
            if atom.GetAtomicNum() not in _PN_ATOMIC_NUMS:
                continue
            if atom.GetIdx() in tagged_p_indices:
                continue  # already handled by the Zone-A lone-pair branch above
            if atom.GetIdx() in tagged_sp3_indices:
                continue  # quaternary N+ handled by the sp3 branch above; do not clear

            stored_cip: str | None = atom.GetPropsAsDict().get("_OIN_CIPCode")
            current_cip = atom.GetPropsAsDict().get("_CIPCode")
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
