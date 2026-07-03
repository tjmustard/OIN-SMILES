"""Module for the xyz2mol functionality for TMCs."""

import argparse
import logging
import signal
import subprocess
from itertools import combinations
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import GetPeriodicTable, rdchem, rdEHTTools, rdmolops
from rdkit.Chem.MolStandardize import rdMolStandardize

from ..core.chirality import ChiralityRecoveryUtility
from .oin_aligner import OINDiscreteAligner, OINSanitizer
from .xyz2mol_local import (
    AC2mol,
    chiral_stereo_check,
    read_xyz_file,
    xyz2AC_obabel,
)

# fmt: off
TRANSITION_METALS = ["Sc","Ti","V","Cr","Mn","Fe","Co","La","Ni","Cu","Zn",
                     "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","Lu",
                     "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
]

TRANSITION_METALS_NUM = [21,22,23,24,25,26,27,57,28,29,30,39,40,41,
                         42,43,44,45,46,47,48,71,72,73,74,75,76,77,78,79,80,
]


ALLOWED_OXIDATION_STATES = {
    "Sc": [3],
    "Ti": [3, 4],
    "V": [2, 3, 4, 5],
    "Cr": [2, 3, 4, 6],
    "Mn": [2, 3, 4, 6, 7],
    "Fe": [2, 3],
    "Co": [2, 3],
    "Ni": [2],
    "Cu": [1, 2],
    "Zn": [2],
    "Y": [3],
    "Zr": [4],
    "Nb": [3, 4, 5],
    "Mo": [2, 3, 4, 5, 6],
    "Tc": [2, 3, 4, 5, 6, 7],
    "Ru": [2, 3, 4, 5, 6, 7, 8],
    "Rh": [1, 3],
    "Pd": [2, 4],
    "Ag": [1],
    "Cd": [2],
    "La": [3],
    "Hf": [4],
    "Ta": [3, 4, 5],
    "W": [2, 3, 4, 5, 6],
    "Re": [2, 3, 4, 5, 6, 7],
    "Os": [3, 4, 5, 6, 7, 8],
    "Ir": [1, 3],
    "Pt": [2, 4],
    "Au": [1, 3],
    "Hg": [1, 2],
}
# fmt: on

logger = logging.getLogger(__name__)

params = Chem.MolStandardize.rdMolStandardize.MetalDisconnectorOptions()
params.splitAromaticC = True
params.splitGrignards = True
params.adjustCharges = False

MetalNon_Hg = (
    "[#3,#11,#12,#19,#13,#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,#39,#40,#41,"
    "#42,#43,#44,#45,#46,#47,#48,#57,#72,#73,#74,#75,#76,#77,#78,#79,#80]"
    "~[B,#6,#14,#15,#33,#51,#16,#34,#52,Cl,Br,I,#85,#1;!$([#1]~[#6])]"
)

pt = GetPeriodicTable

global atomic_valence_electrons

atomic_valence_electrons = {}
atomic_valence_electrons[1] = 1
atomic_valence_electrons[5] = 3
atomic_valence_electrons[6] = 4
atomic_valence_electrons[7] = 5
atomic_valence_electrons[8] = 6
atomic_valence_electrons[9] = 7
atomic_valence_electrons[13] = 3
atomic_valence_electrons[14] = 4
atomic_valence_electrons[15] = 5
atomic_valence_electrons[16] = 6
atomic_valence_electrons[17] = 7
atomic_valence_electrons[18] = 8
atomic_valence_electrons[32] = 4
atomic_valence_electrons[33] = 5  # As
atomic_valence_electrons[35] = 7
atomic_valence_electrons[34] = 6
atomic_valence_electrons[53] = 7

# TMs
atomic_valence_electrons[21] = 3  # Sc
atomic_valence_electrons[22] = 4  # Ti
atomic_valence_electrons[23] = 5  # V
atomic_valence_electrons[24] = 6  # Cr
atomic_valence_electrons[25] = 7  # Mn
atomic_valence_electrons[26] = 8  # Fe
atomic_valence_electrons[27] = 9  # Co
atomic_valence_electrons[28] = 10  # Ni
atomic_valence_electrons[29] = 11  # Cu
atomic_valence_electrons[30] = 12  # Zn

atomic_valence_electrons[39] = 3  # Y
atomic_valence_electrons[40] = 4  # Zr
atomic_valence_electrons[41] = 5  # Nb
atomic_valence_electrons[42] = 6  # Mo
atomic_valence_electrons[43] = 7  # Tc
atomic_valence_electrons[44] = 8  # Ru
atomic_valence_electrons[45] = 9  # Rh
atomic_valence_electrons[46] = 10  # Pd
atomic_valence_electrons[47] = 11  # Ag
atomic_valence_electrons[48] = 12  # Cd

atomic_valence_electrons[57] = 3  # La
atomic_valence_electrons[72] = 4  # Hf
atomic_valence_electrons[73] = 5  # Ta
atomic_valence_electrons[74] = 6  # W
atomic_valence_electrons[75] = 7  # Re
atomic_valence_electrons[76] = 8  # Os
atomic_valence_electrons[77] = 9  # Ir
atomic_valence_electrons[78] = 10  # Pt
atomic_valence_electrons[79] = 11  # Au
atomic_valence_electrons[80] = 12  # Hg


def shell(cmd, shell=False):
    """Run a shell command and return its captured stdout."""
    if shell:
        p = subprocess.Popen(
            cmd,
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    else:
        cmd = cmd.split()
        p = subprocess.Popen(
            cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

    output, err = p.communicate()
    return output


def fix_NO2(mol):
    """Localize nitro groups mis-assigned a charge of -2.

    Such groups are a neutral nitrogen bound to two negatively charged oxygen
    atoms. They are changed to reflect the correct neutral configuration of a
    nitro group. The oxidation state on the transition metal is changed
    accordingly.
    """
    # Create RWMol if not already
    if isinstance(mol, Chem.RWMol):
        emol = mol
    else:
        emol = Chem.RWMol(mol)

    patt = Chem.MolFromSmarts(
        "[#8-]-[#7+0]-[#8-].[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,#39,#40,#41,#42,#43,#44,#45,#46,#47,#48,#57,#72,#73,#74,#75,#76,#77,#78,#79,#80]"
    )
    matches = emol.GetSubstructMatches(patt)
    for a1, a2, a3, a4 in matches:
        if not emol.GetBondBetweenAtoms(a1, a4) and not emol.GetBondBetweenAtoms(a3, a4):
            tm = emol.GetAtomWithIdx(a4)
            o1 = emol.GetAtomWithIdx(a1)
            n = emol.GetAtomWithIdx(a2)
            tm_charge = tm.GetFormalCharge()
            new_charge = tm_charge - 2
            tm.SetFormalCharge(new_charge)
            n.SetFormalCharge(+1)
            o1.SetFormalCharge(0)
            emol.RemoveBond(a1, a2)
            emol.AddBond(a1, a2, rdchem.BondType.DOUBLE)

    Chem.SanitizeMol(emol)
    return emol


def fix_equivalent_Os(mol):
    """Fix a neutral coordinating atom linked to a charged atom via resonance.

    The charge is moved to the coordinating atom and charges fixed accordingly.
    """
    if isinstance(mol, Chem.RWMol):
        emol = mol
    else:
        emol = Chem.RWMol(mol)

    patt = Chem.MolFromSmarts("[#6-,#7-,#8-,#15-,#16-]-[*]=[#6,#7,#8,#15,#16]")

    matches = emol.GetSubstructMatches(patt)
    used_atom_ids_1 = []
    used_atom_ids_3 = []
    for atom in emol.GetAtoms():
        if atom.GetAtomicNum() in TRANSITION_METALS_NUM:
            neighbor_idxs = [a.GetIdx() for a in atom.GetNeighbors()]
            for a1, a2, a3 in matches:
                if (
                    a3 in neighbor_idxs
                    and a1 not in neighbor_idxs
                    and a1 not in used_atom_ids_1
                    and a3 not in used_atom_ids_3
                ):
                    used_atom_ids_1.append(a1)
                    used_atom_ids_3.append(a3)

                    emol.RemoveBond(a1, a2)
                    emol.AddBond(a1, a2, Chem.rdchem.BondType.DOUBLE)
                    emol.RemoveBond(a2, a3)
                    emol.AddBond(a2, a3, Chem.rdchem.BondType.SINGLE)
                    emol.GetAtomWithIdx(a1).SetFormalCharge(0)
                    emol.GetAtomWithIdx(a3).SetFormalCharge(-1)

    Chem.SanitizeMol(emol)
    return emol


def get_proposed_ligand_charge(ligand_mol, cutoff=-10):
    """Run an extended Hückel calculation for the ligand in ligand_mol.

    A suggested charge is found by filling electrons in orbitals <-10eV and
    comparing with total number of valence electrons. If charge is >= 1 (<-1)
    and the LUMO (HOMO) is low (high) in energy, two additional electrons are
    added (removed). The suggested charge is returned.
    """
    valence_electrons = 0
    passed, result = rdEHTTools.RunMol(ligand_mol)
    for a in ligand_mol.GetAtoms():
        valence_electrons += atomic_valence_electrons[a.GetAtomicNum()]

    passed, result = rdEHTTools.RunMol(ligand_mol)
    N_occ_orbs = sum(1 for i in result.GetOrbitalEnergies() if i < cutoff)
    charge = valence_electrons - 2 * N_occ_orbs
    percieved_homo = result.GetOrbitalEnergies()[N_occ_orbs - 1]
    if N_occ_orbs == len(result.GetOrbitalEnergies()):
        percieved_lumo = np.nan
    else:
        percieved_lumo = result.GetOrbitalEnergies()[N_occ_orbs]
    while charge >= 1 and percieved_lumo < -9:
        N_occ_orbs += 1
        charge += -2
        logger.debug("added two more electrons:", charge, percieved_lumo)
        percieved_lumo = result.GetOrbitalEnergies()[N_occ_orbs]
    while charge < -1 and percieved_homo > -10.2:
        N_occ_orbs -= 1
        charge += 2
        logger.debug("removed two electrons:", charge, percieved_homo)
        percieved_homo = result.GetOrbitalEnergies()[N_occ_orbs - 1]

    return charge


def get_basic_mol(xyz_file, overall_charge):
    """Build a basic mol object for an extended Hückel calculation.

    The object is constructed from the adjacency matrix evaluated from the
    xyz-coordinates. All bonds are single bonds, and charges are only assigned
    if necessary to work with it, i.e. a nitrogen with four neighbors gets a
    +1 charge, boron with 4 neighbors gets a -1 charge and oxygen with three
    neighbors gets a +1 charge.
    """
    atoms, _, xyz_coords = read_xyz_file(xyz_file)

    # AC, mol = xyz2AC_huckel(atoms, xyz_coords, overall_charge)
    AC, mol = xyz2AC_obabel(
        atoms, xyz_coords, tolerance=0.5
    )  # Modified tolerance to capture haptic bonds
    tm_indxs = [atoms.index(tm) for tm in TRANSITION_METALS_NUM if tm in atoms]

    rwMol = Chem.RWMol(mol)
    length_ac = len(AC)

    bondTypeDict = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    for i in range(length_ac):
        for j in range(i + 1, length_ac):
            bo = int(round(AC[i, j]))
            if bo == 0:
                continue
            bt = bondTypeDict.get(bo, Chem.BondType.SINGLE)
            rwMol.AddBond(i, j, bt)

    mol = rwMol.GetMol()

    for i, a in enumerate(mol.GetAtoms()):
        if a.GetAtomicNum() == 7:
            # explicit_valence = np.sum(AC[i])
            explicit_valence = sum([ele for idx, ele in enumerate(AC[i]) if idx not in tm_indxs])
            if explicit_valence == 4:
                a.SetFormalCharge(1)
        if a.GetAtomicNum() == 5:
            # Boron with 4 explicit bonds should be negative
            explicit_valence = sum([ele for idx, ele in enumerate(AC[i]) if idx not in tm_indxs])
            if explicit_valence == 4:
                a.SetFormalCharge(-1)
        if a.GetAtomicNum() == 8:
            explicit_valence = sum([ele for idx, ele in enumerate(AC[i]) if idx not in tm_indxs])
            if explicit_valence == 3:
                a.SetFormalCharge(1)

    return mol, xyz_coords


def lig_checks(lig_mol, coordinating_atoms):
    """Sending proposed ligand mol object through series of checks.

    - neighbouring coordinating atoms must be connected by pi-bond, aromatic
      bond (haptic), conjugated system
    - If I have two neighbouring identical charges -> fail, I would rather
      change the charge and make a bond
     -> suggest new charge adding/subtracting electrons based on these neighbouring charges
    - count partial charges: partial charges that are not negative on ligand
      coordinating atoms count against this ligand
      -> loop through resonance forms to see if any live up to this, then choose that one.
      -> partial positive charge on coordinating atom is big red flag
      -> If "bad" partial charges still exists suggest a new charge:
         add/subtract electrons based on the values of the partial charges
    """
    res_mols = rdchem.ResonanceMolSupplier(lig_mol)
    if len(res_mols) == 0:
        res_mols = rdchem.ResonanceMolSupplier(lig_mol, flags=Chem.ALLOW_INCOMPLETE_OCTETS)
    # Check for neighbouring coordinating atoms:
    possible_lig_mols = []

    for res_mol in res_mols:
        positive_atoms = []
        negative_atoms = []
        N_aromatic = 0
        for a in res_mol.GetAtoms():
            if a.GetIsAromatic():
                N_aromatic += 1
            if a.GetFormalCharge() > 0:
                positive_atoms.append(a.GetIdx())
            if a.GetFormalCharge() < 0 and a.GetIdx() not in coordinating_atoms:
                negative_atoms.append(a.GetIdx())

        possible_lig_mols.append((res_mol, len(positive_atoms), len(negative_atoms), N_aromatic))
    return possible_lig_mols


def get_lig_mol(mol, charge, coordinating_atoms):
    """Create a sanitizable mol object for the ligand.

    The checks defined in lig_checks are taken into account.
    We try different charge settings and settings where carbenes are
    allowed/not allowed in case no perfect solution (no partial charges
    on other than the coordinating atoms) can be found. Finally best
    found solution based on criteria in lig_checks is returned.
    """
    atoms = [a.GetAtomicNum() for a in mol.GetAtoms()]
    AC = Chem.rdmolops.GetAdjacencyMatrix(mol)
    lig_mol = AC2mol(mol, AC, atoms, charge, allow_charged_fragments=True, use_atom_maps=False)
    if not lig_mol and charge >= 0:
        charge += -2
        lig_mol = AC2mol(mol, AC, atoms, charge, allow_charged_fragments=True, use_atom_maps=False)
        if not lig_mol:
            return None, charge
    if not lig_mol and charge < 0:
        charge += 2
        lig_mol = AC2mol(mol, AC, atoms, charge, allow_charged_fragments=True, use_atom_maps=False)
        if not lig_mol:
            charge += -4
            lig_mol = AC2mol(
                mol,
                AC,
                atoms,
                charge,
                allow_charged_fragments=True,
                use_atom_maps=False,
            )
            if not lig_mol:
                return None, charge

    possible_res_mols = lig_checks(lig_mol, coordinating_atoms)
    best_res_mol, lowest_pos, lowest_neg, highest_aromatic = possible_res_mols[0]
    for res_mol, N_pos_atoms, N_neg_atoms, N_aromatic in possible_res_mols:
        if N_aromatic > highest_aromatic:
            best_res_mol, lowest_pos, lowest_neg, highest_aromatic = (
                res_mol,
                N_pos_atoms,
                N_neg_atoms,
                N_aromatic,
            )
        if N_aromatic == highest_aromatic and N_pos_atoms + N_neg_atoms < lowest_pos + lowest_neg:
            best_res_mol, lowest_pos, lowest_neg = res_mol, N_pos_atoms, N_neg_atoms
    if lowest_pos + lowest_neg == 0:
        return best_res_mol, charge

    lig_mol_no_carbene = AC2mol(
        mol,
        AC,
        atoms,
        charge,
        allow_charged_fragments=True,
        use_atom_maps=False,
        allow_carbenes=False,
    )
    allow_carbenes = True

    if lig_mol_no_carbene:
        res_mols_no_carbenes = lig_checks(lig_mol_no_carbene, coordinating_atoms)
        for res_mol, N_pos_atoms, N_neg_atoms, N_aromatic in res_mols_no_carbenes:
            if (
                N_aromatic > highest_aromatic
                and N_pos_atoms + N_neg_atoms <= lowest_pos + lowest_neg
            ):
                best_res_mol, lowest_pos, lowest_neg, highest_aromatic = (
                    res_mol,
                    N_pos_atoms,
                    N_neg_atoms,
                    N_aromatic,
                )
            if (
                N_aromatic == highest_aromatic
                and N_pos_atoms + N_neg_atoms < lowest_pos + lowest_neg
            ):
                best_res_mol, lowest_pos, lowest_neg = res_mol, N_pos_atoms, N_neg_atoms
                allow_carbenes = False

    if lowest_pos + lowest_neg == 0:
        logger.debug("found opt solution without carbenes")
        return best_res_mol, charge

    if lowest_pos - lowest_neg + charge < 0:
        new_charge = charge + 2
    else:
        new_charge = charge - 2  # if 0 maybe I should try both

    new_lig_mol = AC2mol(
        mol,
        AC,
        atoms,
        new_charge,
        allow_charged_fragments=True,
        use_atom_maps=False,
        allow_carbenes=allow_carbenes,
    )
    if not new_lig_mol:
        return best_res_mol, charge
    new_possible_res_mols = lig_checks(new_lig_mol, coordinating_atoms)
    for res_mol, N_pos_atoms, N_neg_atoms, N_aromatic in new_possible_res_mols:
        if N_aromatic > highest_aromatic:
            best_res_mol, lowest_pos, lowest_neg, highest_aromatic = (
                res_mol,
                N_pos_atoms,
                N_neg_atoms,
                N_aromatic,
            )
            charge = new_charge
        if N_aromatic == highest_aromatic and N_pos_atoms + N_neg_atoms < lowest_pos + lowest_neg:
            best_res_mol, lowest_pos, lowest_neg = res_mol, N_pos_atoms, N_neg_atoms
            charge = new_charge

    return best_res_mol, charge


def get_tmc_mol(xyz_file, overall_charge, with_stereo=False):
    """Get TMC mol object from given xyz file.

    Args:
        xyz_file (str) : Path to TMC xyz file
        overall_charge (int): Overall charge of TMC
        with_stereo (bool): Whether to percieve stereochemistry from the 3D data

    Returns:
        tmc_mol (rdkit.Chem.rdchem.Mol): TMC mol object
    """
    mol, xyz_coords = get_basic_mol(xyz_file, overall_charge)

    tmc_idx = None
    for a in mol.GetAtoms():
        a.SetIntProp("__origIdx", a.GetIdx())
        if a.GetAtomicNum() in TRANSITION_METALS_NUM:
            # tm_atom = a.GetSymbol()
            tmc_idx = a.GetIdx()

    if tmc_idx is None:
        raise Exception("Found no TM in the input file. Please supply an xyz file with a TM")

    coordinating_atoms = np.nonzero(Chem.rdmolops.GetAdjacencyMatrix(mol)[tmc_idx, :])[0]

    # Exclude aromatic ring carbons that are backbone atoms in bidentate chelates
    # (e.g., the phenylene bridge in bidentate phosphine ligands).
    # If a carbon is in an aromatic ring AND has a neighbor that is a heteroatom donor
    # (P, N, S, O) ALSO in coordinating_atoms, exclude the carbon — the heteroatom is the
    # real donor. This is safe for eta-ligands (Cp, arene) where no ring carbons neighbor
    # heteroatom donors that are coordinating.
    HETEROATOM_DONORS = {7, 8, 15, 16}  # N, O, P, S
    coordinating_set = set(int(idx) for idx in coordinating_atoms)
    filtered_atoms = []
    for atom_idx in coordinating_atoms:
        atom = mol.GetAtomWithIdx(int(atom_idx))
        if atom.GetAtomicNum() == 6 and atom.IsInRing():  # Carbon in a ring
            neighbors = atom.GetNeighbors()
            has_heteroatom_donor_neighbor = any(
                n.GetAtomicNum() in HETEROATOM_DONORS and n.GetIdx() in coordinating_set
                for n in neighbors
            )
            if has_heteroatom_donor_neighbor:
                continue  # Skip this ring carbon; the neighboring heteroatom is the donor
        filtered_atoms.append(atom_idx)
    coordinating_atoms = np.array(filtered_atoms, dtype=int)

    # frags = rdMolStandardize.DisconnectOrganometallics(mol, params)
    mdis = rdMolStandardize.MetalDisconnector(params)
    mdis.SetMetalNon(Chem.MolFromSmarts(MetalNon_Hg))
    frags = mdis.Disconnect(mol)
    frag_mols = rdmolops.GetMolFrags(frags, asMols=True)

    total_lig_charge = 0
    tm_idx = None
    lig_list = []
    for i, f in enumerate(frag_mols):
        m = Chem.Mol(f)
        atoms = m.GetAtoms()
        for atom in atoms:
            if atom.GetAtomicNum() in TRANSITION_METALS_NUM:
                tm_idx = i
                break
        else:
            lig_charge = get_proposed_ligand_charge(f)

            lig_coordinating_atoms = [
                a.GetIdx() for a in m.GetAtoms() if a.GetIntProp("__origIdx") in coordinating_atoms
            ]
            lig_mol, lig_charge = get_lig_mol(m, lig_charge, lig_coordinating_atoms)
            if not lig_mol:
                return None

            # Restore __origIdx from m to lig_mol
            if lig_mol.GetNumAtoms() == m.GetNumAtoms():
                for a_lig, a_orig in zip(lig_mol.GetAtoms(), m.GetAtoms()):
                    if a_orig.HasProp("__origIdx"):
                        a_lig.SetIntProp("__origIdx", a_orig.GetIntProp("__origIdx"))

            total_lig_charge += lig_charge
            lig_list.append(lig_mol)

    if tm_idx is None:
        raise Exception("Found no TM in the input file. Please supply an xyz file with a TM")

    tm = Chem.RWMol(frag_mols[tm_idx])
    tm_ox = overall_charge - total_lig_charge

    len(tm.GetAtoms())

    for a in tm.GetAtoms():
        if a.GetAtomicNum() in TRANSITION_METALS_NUM:
            a.SetFormalCharge(tm_ox)

    for lmol in lig_list:
        tm = Chem.CombineMols(tm, lmol)

    emol = Chem.RWMol(tm)
    coordinating_atoms_idx = [
        a.GetIdx() for a in emol.GetAtoms() if a.GetIntProp("__origIdx") in coordinating_atoms
    ]
    tm_idx = [a.GetIdx() for a in emol.GetAtoms() if a.GetIntProp("__origIdx") == tmc_idx][0]
    dMat = Chem.Get3DDistanceMatrix(emol)
    cut_atoms = []
    for i, j in combinations(coordinating_atoms_idx, 2):
        bond = emol.GetBondBetweenAtoms(int(i), int(j))
        if bond and abs(dMat[i, tm_idx] - dMat[j, tm_idx]) >= 0.4:
            logger.debug(
                "Haptic bond pattern with too great distance:",
                dMat[i, tm_idx],
                dMat[j, tm_idx],
            )
            if dMat[i, tm_idx] > dMat[j, tm_idx] and i in coordinating_atoms_idx:
                coordinating_atoms_idx.remove(i)
                cut_atoms.append(i)
            if dMat[j, tm_idx] > dMat[i, tm_idx] and j in coordinating_atoms_idx:
                coordinating_atoms_idx.remove(j)
                cut_atoms.append(j)
    for j in cut_atoms:
        for i in coordinating_atoms_idx:
            bond = emol.GetBondBetweenAtoms(int(i), int(j))
            if bond and dMat[i, tm_idx] - dMat[j, tm_idx] >= -0.1 and i in coordinating_atoms_idx:
                coordinating_atoms_idx.remove(i)

    for i in coordinating_atoms_idx:
        if emol.GetBondBetweenAtoms(i, tm_idx):
            continue
        emol.AddBond(i, tm_idx, Chem.BondType.DATIVE)

    # Fix specific cases
    # Operate on emol directly to preserve properties
    emol = fix_equivalent_Os(emol)
    emol = fix_NO2(emol)

    tmc_mol = emol.GetMol()
    Chem.SanitizeMol(tmc_mol)
    if with_stereo:
        chiral_stereo_check(tmc_mol)
    return tmc_mol, xyz_coords


def _align_to_pai(tmc_mol, xyz_coords, metal_idx):
    """Canonicalize the orientation of the molecule.

    1. Translates so the metal is at (0,0,0).
    2. Rotates so the Principal Axes of Inertia (PAI) align with the Cartesian axes.
       - Highest Moment of Inertia -> Z
       - Lowest Moment of Inertia -> X
       - Enforce Right-Handed System
    """
    import numpy as np

    coords = np.array(xyz_coords)
    masses = np.array([a.GetMass() for a in tmc_mol.GetAtoms()])

    # 1. Translate Metal to Origin
    metal_pos = coords[metal_idx]
    coords -= metal_pos

    # 2. Calculate Inertia Tensor relative to Origin (Metal)
    # I_jk = sum( m_i * (r_i^2 * delta_jk - r_i_j * r_i_k) )
    # Simplified: We can compute the covariance matrix weighted by mass?
    # Actually PAI axes are eigenvectors of the Inertia Tensor.
    # I = sum_i m_i * ( (r_i.r_i)I - r_i outer r_i )

    inertia = np.zeros((3, 3))
    for i in range(len(coords)):
        m = masses[i]
        pos = coords[i]
        sq_norm = np.dot(pos, pos)

        # Diagonal elements
        inertia[0, 0] += m * (sq_norm - pos[0] * pos[0])
        inertia[1, 1] += m * (sq_norm - pos[1] * pos[1])
        inertia[2, 2] += m * (sq_norm - pos[2] * pos[2])

        # Off-diagonal elements (symmetric, negative product)
        inertia[0, 1] -= m * (pos[0] * pos[1])
        inertia[0, 2] -= m * (pos[0] * pos[2])
        inertia[1, 2] -= m * (pos[1] * pos[2])

    inertia[1, 0] = inertia[0, 1]
    inertia[2, 0] = inertia[0, 2]
    inertia[2, 1] = inertia[1, 2]

    # 3. Diagonalize
    evals, evecs = np.linalg.eigh(inertia)
    # evals are sorted ascending: w1 <= w2 <= w3
    # v1 corresponds to w1 (Lowest Inertia)
    # v3 corresponds to w3 (Highest Inertia)

    # User Spec:
    # Lowest (w1) -> X
    # Highest (w3) -> Z

    x_axis = evecs[:, 0]  # v1
    # y_axis = evecs[:, 1] # v2
    z_axis = evecs[:, 2]  # v3

    # Enforce Right-Handed System: Y = Z x X
    # (Note: X x Z would be -Y)
    y_axis = np.cross(z_axis, x_axis)

    # Construct Rotation Matrix (Rows are new basis vectors)
    # We want to project existing coordinates onto these new axes.
    # New_X = Old_Vec . X_Axis
    R = np.vstack([x_axis, y_axis, z_axis])

    # Apply Rotation
    new_coords = coords @ R.T

    # 4. Handle Degeneracy (Axial Rotation)
    # PAI leaves X/Y rotation arbitrary if Ix ~ Iy (e.g. Ferrocene).
    # We fix this by aligning a "Pivot Atom" to the +X axis.
    # Pivot = Atom with max distance from Z-axis.
    # Tie-breaker: Lowest Index (stable across random rotations if input order is preserved).

    # Calculate Distances from Origin (Metal)
    dists_sq = np.sum(new_coords**2, axis=1)

    # Identify Candidates (within tolerance of max distance)
    # This prevents numerical noise from flipping the pivot among symmetric atoms.
    max_dist_sq = np.max(dists_sq)
    tolerance = 1e-5
    candidates = np.where(dists_sq >= max_dist_sq - tolerance)[0]

    # Tie-breaker: Choose lowest index among candidates
    # Relies on input atom order being preserved (which it is)
    pivot_idx = np.min(candidates)

    # Verify pivot is not on Z-axis (unlikely for max-dist atoms in 3D, unless linear)
    # If it is, we need to pick the next shell?
    # For now, assume molecule isn't linear along Z if we are doing XY alignment.
    # If linear, rotation doesn't matter anyway.

    pivot_pos = new_coords[pivot_idx]
    # Angle in XY plane
    angle = np.arctan2(pivot_pos[1], pivot_pos[0])

    # We want to rotate by -angle around Z to bring pivot to (r, 0, z)
    # R_z = [[c, -s, 0], [s, c, 0], [0, 0, 1]]
    # We want new_pos = R_z(-angle) @ pos
    # cos(-a) = cos(a), sin(-a) = -sin(a)
    c = np.cos(-angle)
    s = np.sin(-angle)

    R_pivot = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    canonical_coords = new_coords @ R_pivot.T

    # 5. Handle Z-Axis Sign Ambiguity
    # PAI Eigenvectors have arbitrary sign. Z vs -Z is random.
    # This affects which ligand is "Top" vs "Bottom".
    # Metric: sum(z_i * (i+1)**3) - Super-linear weighting to break symmetry
    # If negative, flip Z (and Y to maintain right-hand).

    z_moment_idx = 0.0
    for i in range(len(canonical_coords)):
        z_moment_idx += canonical_coords[i][2] * (i + 1) ** 3

    if z_moment_idx < 0:
        # Flip Z -> -Z
        # To maintain RHS (X x Y = Z), if Z flips, we must flip either X or Y.
        # But we determined X is fixed by Pivot.
        # So we must flip Y.
        # Transformation: x->x, y->-y, z->-z
        canonical_coords[:, 1] *= -1
        canonical_coords[:, 2] *= -1

    return canonical_coords.tolist()


def get_oin_string(tmc_mol, xyz_coords):
    """Generates the Open Isomer Notation (OIN) string for the molecule (V2.4).

    Pipeline:
    1. Identify Metal and Connections.
    2. CANONICALIZE ORIENTATION (Translation + PAI Alignment).
    3. Fragment Molecule into Ligands.
    4. Sort Ligands (Mass-First: Component Molecular Weight -> Binding Atom Mass).
    5. Sanitization-First: Force Explicit Hydrogens on Zone A atoms to prevent SMILES drift.
    6. Generate Canonical SMILES for each fragment.
    7. Align Geometry (OINDiscreteAligner V2.4).
    8. Serialize output.
    """
    # 1. Identify Metal and Connections
    metal_idx = -1
    for atom in tmc_mol.GetAtoms():
        if atom.GetAtomicNum() in TRANSITION_METALS_NUM:
            metal_idx = atom.GetIdx()
            break

    # 2. CANONICALIZE ORIENTATION (Translation + PAI Alignment)
    if metal_idx != -1:
        xyz_coords = _align_to_pai(tmc_mol, xyz_coords, metal_idx)

    if metal_idx == -1:
        raise ValueError("No transition metal found in molecule!")

    metal_atom = tmc_mol.GetAtomWithIdx(metal_idx)
    _zone_a_indices = [nbr.GetIdx() for nbr in metal_atom.GetNeighbors()]

    # 2. Fragment Molecule
    mol = Chem.RWMol(tmc_mol)

    # 2. Identify Coordinating Atoms & Bonds
    metal_bonds = mol.GetAtomWithIdx(metal_idx).GetBonds()
    coordinating_atoms = []
    bonds_to_remove = []

    for bond in metal_bonds:
        other_atom = bond.GetOtherAtomIdx(metal_idx)
        coordinating_atoms.append(other_atom)
        bonds_to_remove.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    # 3. Disconnect Metal-Ligand Bonds to define fragments
    for u, v in bonds_to_remove:
        mol.RemoveBond(u, v)

    # 3b. Neutralize Charges for OIN representation
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)  # Start neutral, Sanitize will adjust Zone A

    # 4. Fragment Identification
    frags_indices = Chem.GetMolFrags(mol, asMols=False)

    fragments_data = []  # Will hold dicts with info

    # Extract Coordinates for Alignment from original xyz_coords
    # Map from Mol Idx -> Original XYZ Index
    # We assume 'mol' atoms are in same order as 'tmc_mol' which matches 'xyz_coords'
    # (except maybe if hydrogens were added/removed, but usually 1-to-1).
    # tmc_mol comes from get_tmc_mol which builds from AC2mol but maintains ordering usually?
    # Let's verify: get_tmc_mol returns tmc_mol, xyz_coords.
    # The atom index in tmc_mol should map to xyz_coords index if no hydrogens added implicitly?
    # AC2mol might add H? "use_atom_maps=False".
    # Usually xyz2mol structures match input XYZ atoms unless H added.
    # We can rely on __origIdx if available.

    atom_map_to_xyz = {}
    for atom in mol.GetAtoms():
        if atom.HasProp("__origIdx"):
            atom_map_to_xyz[atom.GetIdx()] = atom.GetIntProp("__origIdx")
        else:
            # Fallback if property missing (shouldn't happen with get_tmc_mol)
            atom_map_to_xyz[atom.GetIdx()] = atom.GetIdx()

    for i, indices in enumerate(frags_indices):
        is_metal = metal_idx in indices

        # Calculate Masses for Sorting
        mass = 0.0
        binding_mass = 0.0

        frag_binding_atoms = []  # List of (global_idx, mass, coords, local_idx_placeholder)

        for idx in indices:
            atom = mol.GetAtomWithIdx(idx)
            m = atom.GetMass()
            mass += m
            if idx in coordinating_atoms:
                binding_mass = max(binding_mass, m)
                # Store info for alignment
                # We need coords
                orig_i = atom_map_to_xyz[idx]
                coords = np.array(xyz_coords[orig_i])
                frag_binding_atoms.append((idx, m, coords))

        if frag_binding_atoms:
            pass

        # Capture metal coordinates if this is the metal fragment
        metal_coords = None
        if is_metal and indices:
            # Use the first atom of the metal fragment as the center
            m_idx = indices[0]
            orig_m_idx = atom_map_to_xyz.get(m_idx, m_idx)
            metal_coords = np.array(xyz_coords[orig_m_idx])

        fragments_data.append(
            {
                "indices": indices,
                "is_metal": is_metal,
                "mass": mass,
                "binding_mass": binding_mass,
                "metal_coords": metal_coords,
                "smiles": "",  # Will be filled later
                "binding_atoms": frag_binding_atoms,  # Will be updated later
            }
        )

    # Now process each fragment to generate SMILES and refine binding_atoms
    for i, item in enumerate(fragments_data):
        indices = item["indices"]  # Old indices
        is_metal = item["is_metal"]
        frag_binding_atoms = item["binding_atoms"]  # Retrieve correct binding atoms

        # 1. Identify Heavy Atoms and Hydrogen Counts
        heavy_indices = []
        h_indices = set()
        atom_h_count = {}  # Old Idx -> Num H neighbors

        for old_idx in indices:
            atom = mol.GetAtomWithIdx(old_idx)
            if atom.GetAtomicNum() == 1:
                h_indices.add(old_idx)
            else:
                heavy_indices.append(old_idx)
                # Count H neighbors within this fragment
                h_count = 0
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 1 and nbr.GetIdx() in indices:
                        h_count += 1
                atom_h_count[old_idx] = h_count

        # 2. Extract Fragment (Heavy Atoms Only, or H atoms for pure-H fragments)
        mw = Chem.RWMol()
        old_to_new = {}

        # For pure-H fragments (hydride ligands like Fe-H after bond cutting),
        # heavy_indices is empty. Fall back to using H atoms directly so that
        # isolated hydrides are serialized as [H] in the OIN rather than "".
        atoms_to_include = heavy_indices if heavy_indices else sorted(h_indices)

        for old_idx in atoms_to_include:
            atom = mol.GetAtomWithIdx(old_idx)
            new_idx = mw.AddAtom(atom)

            # Apply calculated H count from stripped hydrogens
            # This is critical for OINSanitizer to know the original H count
            # since we set NoImplicit=True globally earlier.
            if old_idx in atom_h_count:
                mw.GetAtomWithIdx(new_idx).SetNumExplicitHs(atom_h_count[old_idx])

            old_to_new[old_idx] = new_idx
        # Loop over atoms in fragment, check neighbors
        # Add Bonds between heavy atoms
        for old_idx in heavy_indices:
            atom = mol.GetAtomWithIdx(old_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Check if neighbor is in heavy_indices and we haven't added bond yet
                if nbr_idx in heavy_indices and nbr_idx > old_idx:
                    bond = mol.GetBondBetweenAtoms(old_idx, nbr_idx)
                    if bond:
                        mw.AddBond(old_to_new[old_idx], old_to_new[nbr_idx], bond.GetBondType())

        frag_mol = mw.GetMol()

        # Identify binding atoms in new local indices
        frag_binding_indices_local = []

        # We also need to construct the full 'binding_atoms' list expected by Aligner
        # which needs 'local_idx' (atom index in the generated SMILES fragment)
        # But we don't know the SMILES order yet!
        # OINSanitizer generates SMILES.
        # OINSanitizer takes ligand_mol and binding_indices_in_ligand.

        for binding_item in frag_binding_atoms:
            g_idx = binding_item[0]
            if g_idx in old_to_new:
                frag_binding_indices_local.append(old_to_new[g_idx])

        sanitized_smiles = ""
        sanitized_mol = frag_mol  # Default fallback
        if not is_metal:
            sanitized_smiles, sanitized_mol = OINSanitizer.generate_robust_smiles(
                frag_mol, frag_binding_indices_local
            )
            # Re-apply correct @/@@ on P/N atoms using pre-fragmentation _OIN_CIPCode.
            sanitized_mol = ChiralityRecoveryUtility().recover(sanitized_mol)
            sanitized_smiles = Chem.MolToSmiles(sanitized_mol, isomericSmiles=True, canonical=True)
        else:
            sanitized_smiles = f"[{mol.GetAtomWithIdx(metal_idx).GetSymbol()}]"

        # Now we need the mapping from SMILES order to Fragment Atom order to
        # get 'local_idx' correctly.
        # RDKit's MolToSmiles canonicalization reorders atoms.

        smiles_mol = Chem.MolFromSmiles(sanitized_smiles, sanitize=False)

        if smiles_mol is None:
            # Fallback/Debug: Sanitization produced invalid SMILES?
            logger.error(f"Failed to parse generated SMILES: {sanitized_smiles}")
            frag_to_smiles_idx = {}
        else:
            # We need to map `frag_mol` (original indices) to `smiles_mol`.
            # We use `sanitized_mol` (which has same indices as frag_mol but compatible properties)
            # to match against `smiles_mol`.
            try:
                # Use sanitized_mol as Query? Or Target?
                # SubstructMatch(query). sanitized_mol is the "Source",
                # smiles_mol is "Target/Query structure".
                # match[i] maps query_atom_i to target_atom_j?
                # Mol.GetSubstructMatch(query) -> tuple of atom indices in Mol
                # that match query atoms.
                # So if we use `sanitized_mol.GetSubstructMatch(smiles_mol)`:
                # match[0] is the index in `sanitized_mol` that corresponds to
                # Atom 0 in `smiles_mol`.
                # match[1] is ... Atom 1 in `smiles_mol`.
                # This is exactly what we need.
                match = sanitized_mol.GetSubstructMatch(smiles_mol)
            except Exception as e:
                logger.warning(f"SubstructMatch failed for {sanitized_smiles}: {e}")
                match = None
        # match[i] = index of atom in frag_mol that corresponds to atom i in smiles_mol?
        # No, match is "indices of atoms in frag_mol that match atoms 0,1,2...
        # in query(smiles_mol)".
        # So match[0] is the atom index in frag_mol that corresponds to Atom 0 in SMILES.
        # match[1] is Atom 1 in SMILES, etc.

        # We want: Given atom in frag_mol (which maps to global), what is its SMILES index?
        frag_to_smiles_idx = {}
        if match:
            for s_idx, f_idx in enumerate(match):
                frag_to_smiles_idx[f_idx] = s_idx
        else:
            # Fallback (should not happen if SMILES generated from mol)
            # Maybe due to stereochem differences or sanitization?
            # Safe fallback: linear mapping? No.
            # If match fails, we might have issues.
            pass

        # Update frag_binding_atoms with local_idx
        # stored as: (global_idx, mass, coords)
        # We want a list for Aligner: [(global_idx, mass, coords, local_idx)]

        final_binding_atoms = []
        for g_idx, m, coords in frag_binding_atoms:
            if g_idx in old_to_new:
                l_idx_in_frag = old_to_new[g_idx]
                s_idx = frag_to_smiles_idx.get(l_idx_in_frag, 0)  # Default 0 if fail
                final_binding_atoms.append((g_idx, m, coords, s_idx))

        # Update the item in-place
        item["smiles"] = sanitized_smiles
        item["binding_atoms"] = final_binding_atoms

        # We don't need to append, we are modifying the dict referenced by 'item'
        # fragments_data.append({...}) -> OMITTED to avoid infinite loop

    # 5. Sort Fragments (Deterministic Baseline: Input Order)
    # We sort by the minimum original atom index in the fragment.
    # This ensures that 'fragments_data' order is deterministic based on input file.
    # Metal is Rank 0 (First).

    def get_input_order_key(item):
        if item["is_metal"]:
            return -1  # Metal first
        # Find min original index to ensure deterministic input order
        valid_indices = [atom_map_to_xyz.get(idx, idx) for idx in item["indices"]]
        return min(valid_indices) if valid_indices else float("inf")

    def get_canonical_sort_key(item):
        if item["is_metal"]:
            return (-float("inf"),)  # Metal always first

        # 1. Fragment Molecular Weight (Descending) -> Negate
        mw = item["mass"]

        # 2. Binding Atom Mass (Descending) -> Negate
        # Use binding_mass (max of binding atoms)
        b_mass = item["binding_mass"]

        # 3. SMILES (Ascending)
        smiles = item["smiles"]

        # 4. Tie-Breaker: Input Order (Ascending) for perfect duplicates
        # Matches PRD Step 2 requirement generally, with explicit tie-break
        input_order = get_input_order_key(item)

        return (-mw, -b_mass, smiles, input_order)

    fragments_data.sort(key=get_canonical_sort_key)

    # 6. Run Aligner (On Input-Ordered Fragments)
    # We pass ALL fragments. Metal is Rank 0.
    aligner = OINDiscreteAligner(0, fragments_data)
    geometry_string_raw = aligner.generate_canonical_vectors()

    # geometry_string_raw looks like: "g:SPL|w:1.0:0;2.0:1"
    # We need to parse this to get Slot Assignments for re-sorting.

    # Default if something failed
    if "w:NON" in geometry_string_raw or "error" in geometry_string_raw:
        # Fallback: Just use the current input order.
        full_smiles_parts = [f["smiles"] for f in fragments_data]
        full_smiles = ".".join(full_smiles_parts)
        sidecar_oin = f"{full_smiles} |{geometry_string_raw}|"

    else:
        # Parse Geometry and W-Tag
        # Format: g:GEO|w:Rank.Idx:Slot;...
        parts = geometry_string_raw.split("|")
        geo_tag = parts[0]  # g:GEO
        w_tag = parts[1]  # w:Rank.Idx:Slot

        w_content = w_tag[2:]  # Remove w: prefix

        # Parse assignments: Rank -> List of Slots
        rank_to_slots = {}
        # Also need to track detailed assignments to reconstruct the tag
        pair_data = []  # (Rank, LocalIdx, Slot)

        if w_content and w_content != "NON":
            for entry in w_content.split(";"):
                if not entry:
                    continue
                if ":" not in entry:
                    continue  # Fix: Skip malformed entries
                # entry: "1.0:0"
                # Rank refers to the index in 'fragments_data' (which is currently Input-Ordered)
                parts = entry.split(":")
                if len(parts) != 2:
                    logger.warning(f"Malformed w-tag entry: {entry}")
                    continue
                left, slot_str = parts

                # Need to strip heading chars ^, >, <
                heading_char = ""
                for char in ["^", ">", "<"]:
                    if char in slot_str:
                        heading_char = char
                        break

                slot_str_clean = slot_str.replace("^", "").replace(">", "").replace("<", "")

                slot = int(slot_str_clean)

                if "." in left:
                    r_str, l_str = left.split(".")
                    rank = int(float(r_str))
                    l_idx = int(l_str)
                else:
                    rank = int(float(left))
                    l_idx = 0

                if rank not in rank_to_slots:
                    rank_to_slots[rank] = []
                rank_to_slots[rank].append(slot)
                rank_to_slots[rank].append(slot)
                pair_data.append((rank, l_idx, slot, heading_char))

        # Assign primary slot to each fragment for sorting
        for r, frag in enumerate(fragments_data):
            if frag["is_metal"]:
                frag["_sort_slot"] = -1
            elif r in rank_to_slots:
                # Use minimum slot index as sort key
                frag["_sort_slot"] = min(rank_to_slots[r])
            else:
                frag["_sort_slot"] = float("inf")

        # 7. Re-Sort Fragments by Slot
        # Store original rank to map back for w-tag
        for r, frag in enumerate(fragments_data):
            frag["_orig_rank"] = r

        # Stable sort: Primary=Slot, Secondary=InputOrder
        fragments_data.sort(key=lambda x: (x["_sort_slot"], get_input_order_key(x)))

        # Build Map: OldRank -> NewRank
        old_to_new_rank = {}
        for new_r, frag in enumerate(fragments_data):
            old_r = frag["_orig_rank"]
            old_to_new_rank[old_r] = new_r

        # 8. Re-Construct W-Tag with New Ranks
        # pair_data has (OldRank, LocalIdx, Slot)
        # We need (NewRank, LocalIdx, Slot)

        new_pair_data = []
        for old_r, l_idx, slot, heading_char in pair_data:
            if old_r in old_to_new_rank:
                new_r = old_to_new_rank[old_r]
                new_pair_data.append((new_r, l_idx, slot, heading_char))

        # Sort W-Tag entries by NewRank (then LocalIdx)
        new_pair_data.sort(key=lambda x: (x[0], x[1]))

        new_w_parts = []
        for nr, li, sl, hd in new_pair_data:
            tag = f"{nr}.{li}:{sl}"
            new_pair_data.sort(key=lambda x: (x[0], x[1]))

        new_w_parts = []
        for nr, li, sl, hd_char in new_pair_data:
            tag = f"{nr}.{li}:{sl}"
            tag += hd_char
            new_w_parts.append(tag)
        new_w_tag = "w:" + ";".join(new_w_parts)

        final_geometry_string = f"{geo_tag}|{new_w_tag}"

        # Assemble Final String
        full_smiles_parts = [f["smiles"] for f in fragments_data]
        full_smiles = ".".join(full_smiles_parts)
        sidecar_oin = f"{full_smiles} |{final_geometry_string}|"

    # V3.0 Experimental: Convert to Inline Topology
    from ..oin.inline import OINInlineHandler

    inline_oin = OINInlineHandler.generate_inline_string(sidecar_oin)

    return inline_oin


# Re-export necessary components

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script takes a TMC xyzfile as input and returns a TMC SMILES"
    )
    parser.add_argument("--xyz_file", type=Path, help="The path to a TMC xyz file", required=True)
    parser.add_argument(
        "--charge",
        type=int,
        help="The overall charge of the TMC",
        required=True,
    )
    parser.add_argument(
        "--log_level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "DISABLE"],
        default="INFO",
        help="Set the logging level",
    )

    # Parse arguments
    args = parser.parse_args()

    if args.log_level == "DISABLE":
        logging.disable(logging.CRITICAL)
    else:
        # logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s")
        logging.basicConfig(format="")
        logger.setLevel(getattr(logging, args.log_level))

    # Stop the function if it runs too long.
    def timeout_handler(num, stack):
        """Raise an exception when the SIGALRM watchdog fires."""
        print("Received SIGALRM, terminating")
        raise Exception("Timeout")

    signal.signal(signal.SIGALRM, timeout_handler)

    # Set timeout length
    signal.alarm(300)

    tmc_mol, xyz_coords = get_tmc_mol(args.xyz_file, args.charge, with_stereo=False)

    # Generate OIN
    oin_smiles = get_oin_string(tmc_mol, xyz_coords)

    with open(args.xyz_file.stem + ".txt", "w") as _f:
        _f.write(oin_smiles)

    logger.info(f"Output SMILES: {oin_smiles}")
