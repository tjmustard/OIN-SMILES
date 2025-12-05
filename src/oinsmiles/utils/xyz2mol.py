"Module for the xyz2mol functionality for TMCs"

import argparse
import logging
import signal
import subprocess
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation
from rdkit import Chem
from rdkit.Chem import GetPeriodicTable, rdchem, rdEHTTools, rdmolops
from rdkit.Chem.MolStandardize import rdMolStandardize

from .xyz2mol_local import (
    AC2mol,
    chiral_stereo_check,
    read_xyz_file,
    xyz2AC_obabel,
)
from .oin_aligner import OINCanonicalAligner

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

MetalNon_Hg = "[#3,#11,#12,#19,#13,#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,#39,#40,#41,#42,#43,#44,#45,#46,#47,#48,#57,#72,#73,#74,#75,#76,#77,#78,#79,#80]~[B,#6,#14,#15,#33,#51,#16,#34,#52,Cl,Br,I,#85]"

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
    """Localizes nitro groups that have been assigned a charge of -2 (neutral
    Nitrogen bound to two negatively charged Oxygen atoms).

    These groups are changed to reflect the correct neutral
    configuration of a nitro group. The oxidation state on the
    transition metal is changed accordingly.
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
        if not emol.GetBondBetweenAtoms(a1, a4) and not emol.GetBondBetweenAtoms(
            a3, a4
        ):
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
    """Localizes and fixes where a neutral atom is coordinating to the metal
    but connected ro a negatively charged atom through resonane.

    The charge is moved to the coordinating atom and charges fixed
    accordingly.
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
    """Runs an extended Hückel calculation for the ligand defined in
    ligand_mol.

    A suggested charge is found by filling electrons in orbitals <-10eV
    and omparing with total number of valence electrons. If charge is >=
    1 (<-1) and the LUMO (HOMO) is low (high) in energy, two additional
    electrons are added (removed). The suggested charge is returned.
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
    """A basic mol-object (that can be usedto do an extended Hückel calculation
    is constructed based on the adjacency matrix evaluated from the xyz-
    coordinates.

    All bonds are single bonds, and charges are only asigned if
    necessary to work with it, i.e. a Nitrogen with four neihbors gets a
    +1 charge, Boron with 4 neighbors gets a -1 charge and oxygen with
    three neighbors gets a +1 charge.
    """
    atoms, _, xyz_coords = read_xyz_file(xyz_file)

    # AC, mol = xyz2AC_huckel(atoms, xyz_coords, overall_charge)
    AC, mol = xyz2AC_obabel(atoms, xyz_coords)
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
            explicit_valence = sum(
                [ele for idx, ele in enumerate(AC[i]) if idx not in tm_indxs]
            )
            if explicit_valence == 4:
                a.SetFormalCharge(1)
        if a.GetAtomicNum() == 5:
            # Boron with 4 explicit bonds should be negative
            explicit_valence = sum(
                [ele for idx, ele in enumerate(AC[i]) if idx not in tm_indxs]
            )
            if explicit_valence == 4:
                a.SetFormalCharge(-1)
        if a.GetAtomicNum() == 8:
            explicit_valence = sum(
                [ele for idx, ele in enumerate(AC[i]) if idx not in tm_indxs]
            )
            if explicit_valence == 3:
                a.SetFormalCharge(1)

    return mol, xyz_coords


def lig_checks(lig_mol, coordinating_atoms):
    """Sending proposed ligand mol object through series of checks.

    - neighbouring coordinating atoms must be connected by pi-bond, aromatic bond (haptic), conjugated system
    - If I have two neighbouring identical charges -> fail, I would rather change the charge and make a bond
     -> suggest new charge adding/subtracting electrons based on these neighbouring charges
    - count partial charges: partial charges that are not negative on ligand coordinating atoms count against this ligand
      -> loop through resonance forms to see if any live up to this, then choose that one.
      -> partial positive charge on coordinating atom is big red flag
      -> If "bad" partial charges still exists suggest a new charge: add/subtract electrons based on the values of the partial charges
    """
    res_mols = rdchem.ResonanceMolSupplier(lig_mol)
    if len(res_mols) == 0:
        res_mols = rdchem.ResonanceMolSupplier(
            lig_mol, flags=Chem.ALLOW_INCOMPLETE_OCTETS
        )
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

        possible_lig_mols.append(
            (res_mol, len(positive_atoms), len(negative_atoms), N_aromatic)
        )
    return possible_lig_mols


def get_lig_mol(mol, charge, coordinating_atoms):
    """A sanitizable mol object is created for the ligand, taking into account
    the checks defined in lig_checks.

    We try different charge settings and settings where carbenes are
    allowed/not allowed in case no perfect solution (no partial charges
    on other than the coordinating atoms) can be found. Finally best
    found solution based on criteria in lig_checks is returned.
    """

    atoms = [a.GetAtomicNum() for a in mol.GetAtoms()]
    AC = Chem.rdmolops.GetAdjacencyMatrix(mol)
    lig_mol = AC2mol(
        mol, AC, atoms, charge, allow_charged_fragments=True, use_atom_maps=False
    )
    if not lig_mol and charge >= 0:
        charge += -2
        lig_mol = AC2mol(
            mol, AC, atoms, charge, allow_charged_fragments=True, use_atom_maps=False
        )
        if not lig_mol:
            return None, charge
    if not lig_mol and charge < 0:
        charge += 2
        lig_mol = AC2mol(
            mol, AC, atoms, charge, allow_charged_fragments=True, use_atom_maps=False
        )
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
        if (
            N_aromatic == highest_aromatic
            and N_pos_atoms + N_neg_atoms < lowest_pos + lowest_neg
        ):
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
        if (
            N_aromatic == highest_aromatic
            and N_pos_atoms + N_neg_atoms < lowest_pos + lowest_neg
        ):
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
        raise Exception(
            "Found no TM in the input file. Please supply an xyz file with a TM"
        )

    coordinating_atoms = np.nonzero(Chem.rdmolops.GetAdjacencyMatrix(mol)[tmc_idx, :])[
        0
    ]

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
                a.GetIdx()
                for a in m.GetAtoms()
                if a.GetIntProp("__origIdx") in coordinating_atoms
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
        raise Exception(
            "Found no TM in the input file. Please supply an xyz file with a TM"
        )

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
        a.GetIdx()
        for a in emol.GetAtoms()
        if a.GetIntProp("__origIdx") in coordinating_atoms
    ]
    tm_idx = [
        a.GetIdx() for a in emol.GetAtoms() if a.GetIntProp("__origIdx") == tmc_idx
    ][0]
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
            if (
                bond
                and dMat[i, tm_idx] - dMat[j, tm_idx] >= -0.1
                and i in coordinating_atoms_idx
            ):
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


def get_oin_string(tmc_mol, xyz_coords):
    """Generates the Open Isomer Notation (OIN) string for the molecule.
    
    This function implements OIN v1.3 specification:
    1. Identifies the metal center.
    2. Identifies coordinating atoms (Zone A - bonding atoms).
    3. Groups coordinating atoms into ligands to detect haptic bonds.
    4. Disconnects metal-ligand bonds.
    5. Sorts fragments by Mass (Mass-First Canonicalization).
    6. Generates SMILES string with Selective Explicit Hydrogens:
       - Zone A (bonding atoms): Explicit H required (e.g., [NH3], [Cl])
       - Zone B (backbone atoms): Implicit H for readability (e.g., CC)
    7. Performs Principal Axis Alignment (PAI) on 3D coordinates.
    8. Generates 'w' tag (Zone A unit vectors), 'd' tag, and 'm' tag.
    """
    
    # Work on a copy to avoid modifying the input
    mol = Chem.RWMol(tmc_mol)
    
    # 1. Identify Metal
    metal_idx = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in TRANSITION_METALS_NUM:
            metal_idx = atom.GetIdx()
            metal_orig_idx = atom.GetIntProp("__origIdx")
            break
            
    if metal_idx is None:
        return Chem.MolToSmiles(tmc_mol) + " |w:|"

    # 2. Identify Coordinating Atoms & Bonds
    metal_bonds = mol.GetAtomWithIdx(metal_idx).GetBonds()
    coordinating_atoms = []
    bonds_to_remove = []
    
    for bond in metal_bonds:
        other_atom = bond.GetOtherAtomIdx(metal_idx)
        coordinating_atoms.append(other_atom)
        bonds_to_remove.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    # 3. Detect Haptic vs Dative
    haptic_groups = []
    dative_atoms = []
    
    visited = set()
    for idx in coordinating_atoms:
        if idx in visited:
            continue
            
        cluster = {idx}
        queue = [idx]
        visited.add(idx)
        
        while queue:
            curr = queue.pop(0)
            curr_atom = mol.GetAtomWithIdx(curr)
            for nbr in curr_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in coordinating_atoms and nbr_idx not in visited:
                    visited.add(nbr_idx)
                    cluster.add(nbr_idx)
                    queue.append(nbr_idx)
        
        if len(cluster) >= 3:
            haptic_groups.append(sorted(list(cluster)))
        else:
            for atom_idx in cluster:
                dative_atoms.append(atom_idx)

    # 4. Disconnect Metal-Ligand Bonds
    for u, v in bonds_to_remove:
        mol.RemoveBond(u, v)
        
    # 4b. Neutralize Charges (OIN Style)
    # OIN typically represents components as neutral fragments (e.g. [Cl] not [Cl-]).
    # We set formal charge to 0 for all atoms.
    # We also disable implicit hydrogens to prevent RDKit from adding them to satisfy valence (e.g. Cl -> ClH).
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
    
    # NOTE: In v1.2, we do NOT neutralize charges to preserve explicit H info.
    # CORRECTION: We DO neutralize to match PRD examples (e.g. [Cl], [CH3]).
    # Explicit H count is preserved by the graph structure.
    
    # 5. Mass-First Sorting & Canonicalization
    # We need to sort the disconnected fragments.
    # GetMolFrags returns tuple of atom indices for each fragment
    frags_indices = Chem.GetMolFrags(mol, asMols=False)
    
    # We need to identify which fragment is the metal
    metal_frag_idx = -1
    frag_info = []
    
    for i, indices in enumerate(frags_indices):
        if metal_idx in indices:
            metal_frag_idx = i
            # Metal is always first (Rank 0)
            frag_info.append({'indices': indices, 'mass': float('inf'), 'smi': '', 'is_metal': True})
        else:
            # Calculate mass
            mass = 0.0
            binding_atom_mass = 0.0
            for idx in indices:
                atom = mol.GetAtomWithIdx(idx)
                mass += atom.GetMass()
                if idx in coordinating_atoms:
                    binding_atom_mass = max(binding_atom_mass, atom.GetMass())
            
            # Generate temporary SMILES for tie-breaking
            # We need a mol for this fragment to generate SMILES
            # But we can just use the indices to extract it? 
            # Actually, let's just use the mass and binding atom mass for now.
            # For strict tie-breaking we'd need canonical SMILES of the fragment.
            # Let's extract the fragment as a Mol
            # Note: This is a bit expensive but necessary for strict sorting
            # However, we can't easily extract just by indices without renumbering issues
            # Let's rely on mass and binding atom mass first.
            
            frag_info.append({
                'indices': indices, 
                'mass': mass, 
                'binding_mass': binding_atom_mass,
                'is_metal': False
            })
    
    # Sort fragments: 
    # 1. Metal (already handled by putting it first or giving inf mass)
    # 2. Mass (Descending)
    # 3. Binding Atom Mass (Descending)
    # 4. Indices (as stable tie breaker if needed)
    
    def sort_key(item):
        if item['is_metal']:
            return (float('inf'), float('inf'))
        return (item['mass'], item['binding_mass'])

    frag_info.sort(key=sort_key, reverse=True)
    
    # Construct the new atom order
    new_order = []
    old_to_new = {}
    current_new_idx = 0
    
    for item in frag_info:
        # Within each fragment, we should also canonicalize?
        # Standard SMILES canonicalization within fragment is fine.
        # But we need to map old indices to new indices.
        # Let's just append them in their current internal order for now, 
        # relying on RDKit's MolToSmiles to handle intra-fragment canonicalization later?
        # NO, if we reorder atoms in the Mol, we must ensure MolToSmiles respects that order 
        # OR we let MolToSmiles do its thing and we map back.
        # But MolToSmiles on a disconnected mol might reorder fragments?
        # Yes, standard MolToSmiles canonicalizes the whole string.
        # We want to ENFORCE our fragment order.
        # So we should generate SMILES for each fragment individually and join them.
        
        # For the OIN string, we need the SMILES.
        # For the tags, we need the indices in that SMILES string.
        # So we must know exactly how the final SMILES is ordered.
        
        # Strategy:
        # 1. Extract each fragment as a separate Mol.
        # 2. Generate SMILES for each fragment (Explicit H).
        # 3. Join SMILES with '.'.
        # 4. Track atom mapping from Original -> Fragment -> Final Combined Index.
        pass

    # Let's implement the Strategy
    final_smiles_parts = []
    final_atom_mapping = {} # Old Idx -> New Idx in the virtual combined molecule
    
    # We need to process the Metal fragment first, then the sorted ligands
    
    current_offset = 0
    
    # Helper to extract fragment mol
    def get_frag_mol(mol, indices):
        # Create a mapping from old idx to new frag idx
        mapping = {old: new for new, old in enumerate(indices)}
        
        # Create new empty mol
        frag = Chem.RWMol()
        for old_idx in indices:
            atom = mol.GetAtomWithIdx(old_idx)
            new_idx = frag.AddAtom(atom)
            # Copy properties if needed
        
        # Add bonds
        for old_idx in indices:
            atom = mol.GetAtomWithIdx(old_idx)
            for bond in atom.GetBonds():
                nbr_idx = bond.GetOtherAtomIdx(old_idx)
                if nbr_idx in indices and nbr_idx > old_idx: # Add each bond once
                    bond_type = bond.GetBondType()
                    frag.AddBond(mapping[old_idx], mapping[nbr_idx], bond_type)
        
        return frag, mapping

    for item in frag_info:
        indices = item['indices'] # Old indices
        # We want to generate canonical SMILES for this fragment
        # AND know the mapping from Old Index -> Index in that SMILES string
        
        # Extract fragment
        frag_mol, frag_local_map = get_frag_mol(mol, indices)
        
        # OIN V1.3: Selective Explicit Hydrogen Handling
        # Zone A (bonding atoms): Use bracket notation with H subscript (e.g., [NH3], [Cl])
        # Zone B (backbone atoms): Use standard implicit H (e.g., CC, c)
        
        # Key insight: ALL H atoms from XYZ are explicit in the graph
        # Solution: Remove ALL H atoms, then set explicit H count ONLY for Zone A
        # This gives [NH2] for Zone A N, but CC for Zone B carbons
        
        # Identify all H atoms in this fragment
        all_h_atoms = []
        h_parent_map = {}  # H local idx -> parent local idx
        zone_a_heavy_atoms = []
        
        for old_idx in indices:
            frag_local_idx = frag_local_map[old_idx]
            atom = mol.GetAtomWithIdx(old_idx)
            
            if atom.GetAtomicNum() == 1:  # Hydrogen
                # Find what it's bonded to
                for bond in atom.GetBonds():
                    parent_idx = bond.GetOtherAtomIdx(old_idx)
                    if parent_idx in indices:  # Parent is in this fragment
                        all_h_atoms.append(old_idx)
                        h_parent_map[frag_local_idx] = frag_local_map[parent_idx]
                        break
            elif old_idx in coordinating_atoms:  # Zone A heavy atom
                zone_a_heavy_atoms.append(old_idx)
        
        # Remove ALL H atoms from fragment
        atoms_to_remove = [frag_local_map[h_idx] for h_idx in all_h_atoms]
        
        new_frag_mol = Chem.RWMol(frag_mol)
        for h_local_idx in sorted(atoms_to_remove, reverse=True):
            new_frag_mol.RemoveAtom(h_local_idx)
        
        frag_mol = new_frag_mol.GetMol()
        
        # Build mapping after H removal
        old_to_new_local = {}
        new_idx = 0
        for old_local_idx in range(len(indices)):
            if old_local_idx not in atoms_to_remove:
                old_to_new_local[old_local_idx] = new_idx
                new_idx += 1
        
        # Set explicit H count ONLY for Zone A atoms
        for old_idx in zone_a_heavy_atoms:
            old_local_idx = frag_local_map[old_idx]
            if old_local_idx not in old_to_new_local:
                continue
            new_local_idx = old_to_new_local[old_local_idx]
            atom = frag_mol.GetAtomWithIdx(new_local_idx)
            
            # Count H atoms that were bonded to this atom
            h_count = sum(1 for h_local in atoms_to_remove if h_parent_map.get(h_local) == old_local_idx)
            atom.SetNumExplicitHs(h_count)
            # Force bracket notation
            atom.SetAtomMapNum(new_local_idx + 1)
            # Set NoImplicit only for Zone A
            atom.SetNoImplicit(True)
        
        # Reset NoImplicit for Zone B atoms (Backbone)
        for old_idx in indices:
            if old_idx in zone_a_heavy_atoms:
                continue
            if old_idx in all_h_atoms:
                continue
                
            old_local_idx = frag_local_map[old_idx]
            if old_local_idx in old_to_new_local:
                new_local_idx = old_to_new_local[old_local_idx]
                atom = frag_mol.GetAtomWithIdx(new_local_idx)
                atom.SetNoImplicit(False)
        
        # Zone B atoms now have no H atoms and no explicit H set
        # RDKit will add implicit H automatically when generating SMILES
        
        try:
            Chem.SanitizeMol(frag_mol)
        except Exception as e:
            logger.warning(f"Sanitization failed for fragment: {e}")
        
        # Generate SMILES with map numbers to ensure unique mapping
        frag_smiles_with_maps = Chem.MolToSmiles(frag_mol, canonical=True, allHsExplicit=False)
        
        # Remove atom map numbers for the final string
        import re
        frag_smiles = re.sub(r':(\d+)', '', frag_smiles_with_maps)
        
        final_smiles_parts.append(frag_smiles)
        
        # Update global mapping (accounting for removed H atoms)
        remaining_local_to_old = {}
        for old_idx in indices:
            old_local_idx = frag_local_map[old_idx]
            if old_local_idx in old_to_new_local:
                new_local_idx = old_to_new_local[old_local_idx]
                remaining_local_to_old[new_local_idx] = old_idx
        
        # Create a temporary mol from the SMILES to get the SMILES order
        # Note: We use the version WITH maps to ensure robust matching
        temp_mol = Chem.MolFromSmiles(frag_smiles_with_maps)
        
        if temp_mol is None:
             # Fallback if something goes wrong (shouldn't happen)
             logger.warning(f"Failed to parse generated SMILES: {frag_smiles_with_maps}")
             # Use rank based mapping as fallback
             ranks = list(Chem.CanonicalRankAtoms(frag_mol, breakTies=True))
             rank_to_local = {r: i for i, r in enumerate(ranks)}
             for r in range(len(frag_mol.GetAtoms())):
                local_idx = rank_to_local[r]
                old_idx = remaining_local_to_old[local_idx]
                new_global_idx = current_offset + r
                old_to_new[old_idx] = new_global_idx
        else:
            # Map frag_mol atoms to temp_mol atoms (SMILES order)
            # match[i] is the index in frag_mol that corresponds to atom i in temp_mol
            match = frag_mol.GetSubstructMatch(temp_mol)
            
            if not match:
                logger.warning(f"Failed to match frag_mol to SMILES mol: {frag_smiles_with_maps}")
                 # Fallback
                ranks = list(Chem.CanonicalRankAtoms(frag_mol, breakTies=True))
                rank_to_local = {r: i for i, r in enumerate(ranks)}
                for r in range(len(frag_mol.GetAtoms())):
                    local_idx = rank_to_local[r]
                    old_idx = remaining_local_to_old[local_idx]
                    new_global_idx = current_offset + r
                    old_to_new[old_idx] = new_global_idx
            else:
                # match[smiles_idx] = frag_local_idx
                for smiles_idx, frag_local_idx in enumerate(match):
                    old_idx = remaining_local_to_old[frag_local_idx]
                    new_global_idx = current_offset + smiles_idx
                    old_to_new[old_idx] = new_global_idx
            
        current_offset += len(frag_mol.GetAtoms())

    final_smiles = ".".join(final_smiles_parts)
    
    # 7. Canonical Alignment (PBCA - OIN V1.4)
    # Center on Metal
    metal_coords = np.array(xyz_coords[metal_orig_idx])
    centered_coords = np.array(xyz_coords) - metal_coords
    
    # Prepare data for Aligner
    ligands_for_aligner = []
    
    # frag_info is sorted, final_smiles_parts corresponds to it.
    for i, item in enumerate(frag_info):
        if item['is_metal']:
            continue
            
        binding_atoms_data = []
        frag_indices = set(item['indices'])
        
        # We need to find binding atoms for this specific ligand
        # Coordinating atoms are global indices.
        # We also need to respect the order of sorting binding atoms if multiple?
        # The Aligner expects: 'binding_atoms': list of (original_index, atomic_mass, coords_xyz)
        # xyz2mol uses all binding atoms for PAI previously.
        # Here we just gather them. The Aligner will pick the best one (highest mass).
        
        # Filter coordinating atoms that belong to this fragment
        lig_binding_indices = [idx for idx in coordinating_atoms if idx in frag_indices]
        
        # Sort them by mass descending (as per PRD implication for "first binding atom")
        # PRD Step 1: "Sort Ligands ... Binding Atom Mass".
        # PRD Step 3 Phase A: "Select Zone A atom... highest atomic mass".
        # So we should provide them sorted or let Aligner sort.
        # The Aligner code I wrote: `p1_atom_coords = p1_ligand['binding_atoms'][0][2]`
        # So it takes the FIRST one. Thus I MUST sort them here by mass descending.
        
        lig_binding_indices.sort(key=lambda idx: mol.GetAtomWithIdx(idx).GetMass(), reverse=True)
        
        for atom_idx in lig_binding_indices:
            atom = mol.GetAtomWithIdx(atom_idx)
            mass = atom.GetMass()
            if atom.HasProp("__origIdx"):
                orig_idx = atom.GetIntProp("__origIdx")
            else:
                orig_idx = atom_idx # Fallback
            
            coords = centered_coords[orig_idx]
            binding_atoms_data.append((orig_idx, mass, coords))
            
        ligands_for_aligner.append({
            'smiles': final_smiles_parts[i],
            'mass': item['mass'],
            'binding_atoms': binding_atoms_data
        })
        
    if ligands_for_aligner:
        aligner = OINCanonicalAligner(ligands_for_aligner)
        rotation_matrix = aligner.get_best_alignment()
        
        # Apply rotation: v_new = v_old . R^T
        aligned_coords = np.dot(centered_coords, rotation_matrix.T)
    else:
        # No ligands (e.g. naked metal ion?), identity
        aligned_coords = centered_coords
    
    # 8. Generate Tags (V1.4)
    
    # v tag: Unified Vector Tags (Connectivity + Geometry)
    # Format: v:MetalIdx.LigandIdx:x,y,z;...
    v_entries = []
    
    new_metal_idx = old_to_new[metal_idx]
    
    # Sort by new index for canonical output
    sorted_coordinating = sorted(coordinating_atoms, key=lambda x: old_to_new[x])
    
    for old_idx in sorted_coordinating:
        new_ligand_idx = old_to_new[old_idx]
        # Use __origIdx to access aligned_coords (which matches xyz_coords order)
        atom_orig_idx = mol.GetAtomWithIdx(old_idx).GetIntProp("__origIdx")
        vec = aligned_coords[atom_orig_idx]
        norm = np.linalg.norm(vec)
        if norm > 1e-6:
            u_vec = vec / norm
        else:
            u_vec = vec # Should not happen for bonded atoms
            
        u_vec = np.round(u_vec, 3)
        u_vec[u_vec == -0.0] = 0.0
            
        v_entries.append(f"{new_metal_idx}.{new_ligand_idx}:{u_vec[0]:.3f},{u_vec[1]:.3f},{u_vec[2]:.3f}")
        
    v_tag = "v:" + ";".join(v_entries)
    
    # d tag is REMOVED in V1.4
        
    # m tag: Haptic Groups
    m_tags = []
    # Sort groups by smallest new index in group
    sorted_haptic = []
    for group in haptic_groups:
        new_indices = sorted([old_to_new[idx] for idx in group])
        sorted_haptic.append(new_indices)
        
    sorted_haptic.sort(key=lambda x: x[0])
    
    for group in sorted_haptic:
        atom_list = ".".join(map(str, group))
        m_tags.append(f"m:{new_metal_idx}:{atom_list}")
        
    # Assemble OIN
    tags = [v_tag]
    # d tag removed
    if m_tags:
        tags.extend(m_tags)
        
    # Add g tag placeholder if needed, but PRD says it's descriptive.
    # We'll omit it for now as we don't have geometry classification logic here.
        
    oin = f"{final_smiles} |{'|'.join(tags)}|"
    return oin


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script takes a TMC xyzfile as input and returns a TMC SMILES"
    )
    parser.add_argument(
        "--xyz_file", type=Path, help="The path to a TMC xyz file", required=True
    )
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
