def _clear_chelate_locked_bond_stereo(tmc_mol):
    """Drop E/Z from double bonds that a metal-containing ring holds rigid.

    A chelating ligand closes its ring *through the metal* (acac, salicylaldiminate,
    beta-diketiminate, nacnac, bis-imine). ``get_oin_string`` then fragments the
    complex and strips the metal, which opens that ring -- so a C=C or C=N that is
    physically ring-locked looks acyclic to RDKit and gets a directional marker
    whose sign falls out of the SMILES traversal order rather than the chemistry.
    Re-encoding the generated 3D structure traverses differently and flips the
    marker, so the round trip fails on a bond that never had a free E/Z to begin
    with. An eta-bound alkene forms a metal-C-C ring and is suppressed for the
    same reason.

    This is the encoder-side counterpart of the coordination-sphere filter the
    generator already applies in ``generator3d/ligand.py`` when deciding which
    double bonds are worth constraining at embed time. A pendant, freely rotatable
    alkene lies in no metal ring and keeps its marker.

    Ring perception runs on a copy whose DATIVE bonds are upgraded to SINGLE:
    RDKit's SSSR **ignores dative bonds**, so every chelate ring is invisible to it
    (VOacac2: 30 atoms, 31 bonds, one fragment -- two cycles by Euler -- and
    ``GetSymmSSSR`` returns 0). That blindness is what lets DetectBondStereochemistry
    treat a ring-locked alkene as acyclic in the first place.
    """
    from rdkit import Chem

    from .constants import TRANSITION_METALS_NUM

    metal_indices = {
        a.GetIdx() for a in tmc_mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM
    }
    if not metal_indices:
        return

    try:
        probe = Chem.RWMol(tmc_mol)
        for bond in probe.GetBonds():
            if bond.GetBondType() == Chem.BondType.DATIVE:
                bond.SetBondType(Chem.BondType.SINGLE)
        probe = probe.GetMol()
        Chem.FastFindRings(probe)
        rings = Chem.GetSymmSSSR(probe)
    except Exception:
        return

    locked_atoms = set()
    for ring in rings:
        ring = set(ring)
        if ring & metal_indices:
            locked_atoms |= ring

    for bond in tmc_mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        if bond.GetBeginAtomIdx() in locked_atoms and bond.GetEndAtomIdx() in locked_atoms:
            bond.SetStereo(Chem.BondStereo.STEREONONE)
            bond.SetBondDir(Chem.BondDir.NONE)
            for nbr in (bond.GetBeginAtom(), bond.GetEndAtom()):
                for b in nbr.GetBonds():
                    if b.GetBondType() == Chem.BondType.SINGLE:
                        b.SetBondDir(Chem.BondDir.NONE)


class XYZToSMILES:
    """Convert an XYZ structure file into an OIN-SMILES string."""

    def convert(self, xyz_file_path: str) -> str:
        """Converts an XYZ file to an OIN-SMILES string."""
        from pathlib import Path

        from rdkit import Chem

        from ..core.chirality import CIPAssigner
        from ..utils.xyz2mol import get_oin_string, get_tmc_mol

        charge = 0  # Default
        path = Path(xyz_file_path)

        try:
            tmc_mol, xyz_coords = get_tmc_mol(path, charge, with_stereo=False)
        except Exception as e:
            raise ValueError(f"xyz2mol failed: {e}")

        # Assign 3D-derived CIP codes to P/N stereocenters before fragmentation.
        # tmc_mol from get_tmc_mol() is already sanitized and has a valid 3D conformer.
        Chem.SanitizeMol(tmc_mol)

        # Perceive C=C (cis/trans) stereo from the 3D geometry so get_oin_string can
        # carry it into the OIN string (its E/Z carry reads bond.GetStereo() on this
        # mol). Deliberately scoped to double bonds only -- this does NOT run
        # AssignAtomChiralTagsFromStructure, so it does not enforce sp3 handedness the
        # generator does not reproduce (which would regress ~20% of complexes). This
        # is what makes the OIN round trip lossless for, and able to verify, cis/trans.
        Chem.DetectBondStereochemistry(tmc_mol, -1)
        Chem.AssignStereochemistry(tmc_mol, force=True)

        _clear_chelate_locked_bond_stereo(tmc_mol)

        CIPAssigner().assign_all(tmc_mol)

        # 2. Generate OIN (ChiralityRecoveryUtility is applied inside get_oin_string)
        oin_string = get_oin_string(tmc_mol, xyz_coords)

        return oin_string
