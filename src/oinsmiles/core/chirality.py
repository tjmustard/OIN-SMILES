"""Chirality encoding utilities for the XYZ→OIN pipeline.

Three classes handle the full lifecycle of P/N stereocenter encoding:
  - CIPAssigner: pre-fragmentation — reads 3D conformer, stores _OIN_CIPCode on P/N atoms.
  - ChiralityRecoveryUtility: post-fragmentation — verifies/corrects chiral tags after
    OINSanitizer strips the metal context.
  - PseudoAtomStrategy: fallback for P/N atoms with 4 neighbors but no computable CIP code.
"""

from __future__ import annotations

from rdkit import Chem

#: Sentinel atomic number used by PseudoAtomStrategy.
#: Z=0 is the RDKit wildcard atom ('*').
PSEUDO_ATOMIC_NUM: int = 0

#: RDKit atomic numbers for nitrogen (N) and phosphorus (P).
_PN_ATOMIC_NUMS: frozenset[int] = frozenset({7, 15})


class PseudoAtomStrategy:
    """Fallback for P/N atoms with 4 neighbours but no computable CIP code.

    Wildcard atoms (Z=0, '*') must be stripped before OIN serialization via
    ``strip_pseudo_atoms()``.
    """

    PSEUDO_ATOMIC_NUM: int = PSEUDO_ATOMIC_NUM

    @staticmethod
    def strip_pseudo_atoms(mol: Chem.Mol) -> Chem.Mol:
        """Remove all wildcard (*) atoms (Z=0) from *mol* before OIN output.

        Atoms are removed in reverse index order to keep remaining indices
        stable during iteration.
        """
        rw = Chem.RWMol(mol)
        indices = sorted(
            [a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 0],
            reverse=True,
        )
        for idx in indices:
            rw.RemoveAtom(idx)
        return rw.GetMol()


class CIPAssigner:
    """Assigns 3D-derived CIP codes and chiral tags to P/N stereocenters.

    Must be called on the **full** TMC mol (pre-fragmentation) so that the
    metal coordination context is available for correct CIP assignment.

    Preconditions
    -------------
    * *mol* must have a valid embedded 3D conformer (as returned by
      ``get_tmc_mol()``).
    * ``Chem.SanitizeMol(mol)`` must succeed — any sanitisation exception
      propagates to the caller.

    Postconditions
    --------------
    * Each P/N atom with an assignable CIP code gains the atom property
      ``_OIN_CIPCode`` ('R' or 'S').
    * Chiral tags (``CHI_TETRAHEDRAL_CW`` / ``CHI_TETRAHEDRAL_CCW``) are set
      from the 3D conformer geometry.
    """

    def assign_all(self, mol: Chem.Mol) -> Chem.Mol:
        """Assign stereo tags and store ``_OIN_CIPCode`` on all P/N atoms.

        Parameters
        ----------
        mol:
            RDKit ``Mol`` with an embedded 3D conformer.  Must not be ``None``.

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

        return mol


class ChiralityRecoveryUtility:
    """Verifies and corrects @/@@ on P/N atoms in post-fragmentation fragment mols.

    After ``OINSanitizer.generate_robust_smiles()`` the metal neighbour is
    absent.  This utility compares the CIP code computed from the surviving
    chiral tags against the pre-fragmentation ``_OIN_CIPCode`` stored by
    ``CIPAssigner``, and flips the tag when they disagree.

    Zone A P/N atoms (direct metal binders, fewer than 4 total neighbours in
    the fragment) have their chiral tags cleared — per the MiniPRD constraint
    that @/@@ on Zone A atoms requires the full coordination sphere.
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

        # Compute CIP codes from existing chiral tags (no 3D conformer needed).
        Chem.AssignStereochemistry(rw, cleanIt=True, force=True)

        for atom in rw.GetAtoms():
            if atom.GetAtomicNum() not in _PN_ATOMIC_NUMS:
                continue

            stored_cip: str | None = atom.GetPropsAsDict().get("_OIN_CIPCode")
            current_cip: str | None = atom.GetPropsAsDict().get("_CIPCode")
            total_deg: int = atom.GetTotalDegree()

            if total_deg < 4:
                # Zone A atom — metal removed, 3 neighbours remain.
                # @/@@ is undefined without the full coordination sphere.
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
                # PseudoAtomStrategy fallback: clear the stray chiral tag.
                atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

        return rw.GetMol()
