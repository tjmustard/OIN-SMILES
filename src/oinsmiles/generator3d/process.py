"""---generate_molecule.py---.

Generate 2D graph molecule in ACE-Reaction format from xyz files or SMILES.
"""

import copy

import numpy as np
from rdkit import Chem
from scipy import spatial

from . import chem
from .utils import compute_chg_and_bo_pulp as compute_pulp


def get_block_diagonal_adj_from_fragments(fragment_adj_list, total_num_atom=0):
    """It returns a single block diagonal matrix by concatenating fragment adjacency matrices."""
    if total_num_atom == 0:
        for fragment_adj in fragment_adj_list:
            total_num_atom += len(fragment_adj)
    total_diagonal_matrix = np.zeros((total_num_atom, total_num_atom))
    cnt = 0
    for i in range(len(fragment_adj_list)):
        fragment_adj = fragment_adj_list[i]
        size = len(fragment_adj)
        for j in range(size):
            for k in range(size):
                total_diagonal_matrix[cnt + j][cnt + k] = fragment_adj[j][k]
        cnt += size
    return total_diagonal_matrix


def locate_molecule(ace_molecule, coordinate_list, update=False):
    """Locates atoms according to coordinate_list, be cautious on ordering of atoms."""
    atom_list = ace_molecule.atom_list
    for i in range(len(atom_list)):
        atom = atom_list[i]
        atom.set_coordinate(coordinate_list[i])
    if update:
        ace_molecule.set_adj_matrix(None)


def translate_molecule(ace_molecule, vector):
    """Translate the molecule."""
    atom_list = ace_molecule.atom_list
    for atom in atom_list:
        translate_atom(atom, vector)


def locate_atom(atom, coordinate):
    """Locates a single atom to input 'coordinate'."""
    atom.x = coordinate[0]
    atom.y = coordinate[1]
    atom.z = coordinate[2]


def translate_atom(atom, vector):
    """Translate the atom."""
    atom.x += vector[0]
    atom.y += vector[1]
    atom.z += vector[2]


def read_molecule(f, extension="xyz"):
    """Read the molecule."""
    molecule = chem.Molecule()
    atom_list = []
    info = []
    if extension == "com":
        try:
            chg, multiplicity = f.readline().strip().split()
            chg = int(chg)
            multiplicity = int(multiplicity)
            molecule.chg = chg
            molecule.multiplicity = multiplicity
        except Exception:
            print("Wrong format ! Should start with molecular charge and multiplicity !!!")
            return molecule, info
        while True:
            try:
                atom_line = f.readline().strip().split()
                print(atom_line)
                element = atom_line[0]
                x = float(atom_line[1])
                y = float(atom_line[2])
                z = float(atom_line[3])
                atom = chem.Atom(element)
                atom.x = x
                atom.y = y
                atom.z = z
                atom_list.append(atom)
            except Exception:
                break
        molecule.atom_list = atom_list
        return molecule, []

    elif extension == "xyz":
        try:
            atom_num = int(f.readline().strip())
            info = f.readline().strip().split()
        except Exception:
            atom_num = 0
        for i in range(atom_num):
            atom_line = f.readline().strip().split()
            element = atom_line[0]
            x = float(atom_line[1])
            y = float(atom_line[2])
            z = float(atom_line[3])
            atom = chem.Atom(element)
            atom.x = x
            atom.y = y
            atom.z = z
            atom_list.append(atom)
        molecule.atom_list = atom_list
        return molecule, info


def get_rd_mol_from_np_array(z_list, chg_list, bo_matrix):
    """Return the rd mol from np array."""
    n = len(chg_list)
    bond_types = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    # First generate atom_list
    rd_mol = Chem.Mol()
    rde_mol = Chem.EditableMol(rd_mol)
    for i in range(n):
        rd_atom = Chem.Atom(int(z_list[i]))
        rd_atom.SetFormalCharge(int(chg_list[i]))
        rde_mol.AddAtom(rd_atom)
    # Get bond list from np.array bo_matrix (Never use), you can expand this part
    bond_order_list = [1, 2, 3]
    for bond_order in bond_order_list:
        np_array_list = np.where(bo_matrix == bond_order)
        np_array_list = np.stack(np_array_list, axis=1)
        for array in np_array_list:
            begin = array[0]
            end = array[1]
            if begin < end:
                rde_mol.AddBond(int(begin), int(end), bond_types[bond_order])
    rd_mol = rde_mol.GetMol()
    return rd_mol


def get_ace_mol_with_coordinate(smiles):
    """Return the ace mol with coordinate."""
    molecule = chem.Molecule(smiles)
    coordinate_list = molecule.make_3d_coordinate()
    locate_molecule(molecule, coordinate_list)
    return molecule


def _dearomatize_stuck_rings(rd_molecule, add_hydrogen):
    """Kekulize a molecule whose aromatic system has no global Kekule structure.

    A quinoid ring -- an aromatic ring atom carrying an exocyclic double bond,
    e.g. a 2-iminopyridine / amidinate donor (the Ti/Hf ``no_conformers`` cases)
    -- makes ``Chem.Kekulize`` raise. The blanket fallback of treating EVERY
    remaining aromatic bond as single also degrades the ligand's *other*,
    well-behaved rings (a 2,6-dimethylphenyl bridge), so MetalloGen embeds those
    arenes non-planar and they re-encode as quinoid -> a spurious string
    mismatch. Instead, clear aromaticity on ONLY the stuck ring(s) and kekulize
    the rest normally, preserving correct alternating bond orders on the good
    rings. Returns a best-effort mol (the original on any failure); the
    bond-order lookup below still guards against a stray AROMATIC bond.
    """
    try:
        mol = Chem.AddHs(rd_molecule) if add_hydrogen else Chem.RWMol(rd_molecule).GetMol()
        mol = Chem.RWMol(mol)
        mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        except Exception:
            pass
        stuck = set()
        for ring in mol.GetRingInfo().AtomRings():
            ring_set = set(ring)
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if not atom.GetIsAromatic():
                    continue
                for b in atom.GetBonds():
                    if (
                        b.GetOtherAtomIdx(idx) not in ring_set
                        and b.GetBondType() == Chem.BondType.DOUBLE
                    ):
                        stuck.update(ring)
        if not stuck:
            return rd_molecule
        for idx in stuck:
            mol.GetAtomWithIdx(idx).SetIsAromatic(False)
        for b in mol.GetBonds():
            if b.GetBeginAtomIdx() in stuck and b.GetEndAtomIdx() in stuck:
                b.SetIsAromatic(False)
                if b.GetBondType() == Chem.BondType.AROMATIC:
                    b.SetBondType(Chem.BondType.SINGLE)
        out = mol.GetMol()
        out.UpdatePropertyCache(strict=False)
        Chem.Kekulize(out, clearAromaticFlags=True)
        return out
    except Exception:
        return rd_molecule


def get_ace_mol_from_rd_mol(rd_molecule, add_hydrogen=True, include_stereo=False):
    """It converts rd_molecule type info ace_molecule type."""
    # Kekulize molecule
    rd_molecule_copy = copy.deepcopy(rd_molecule)
    try:
        if add_hydrogen:
            rd_molecule = Chem.AddHs(rd_molecule)
        Chem.rdmolops.Kekulize(rd_molecule)
    except Exception:
        # Some aromatic ring has no valid Kekule structure (a quinoid ring: an
        # aromatic ring atom with an exocyclic double bond). Dearomatize ONLY that
        # ring so the rest of the molecule keeps correct bond orders, rather than
        # reverting the whole ligand to aromatic (which then hits the AROMATIC
        # KeyError below and aborts generation with "failed to generate any
        # conformers").
        rd_molecule = _dearomatize_stuck_rings(rd_molecule_copy, add_hydrogen)
    bond_types = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}
    n = rd_molecule.GetNumAtoms()
    atom_list = []
    chg_list = []
    atom_feature = dict()
    # Make atom_list
    for i in range(n):
        rd_atom = rd_molecule.GetAtomWithIdx(i)
        ace_atom = chem.Atom()
        chg_list.append(rd_atom.GetFormalCharge())
        """
        position = rd_molecule.GetAtomPosition(i)
        if position!=None:
            ace_atom.x = position[0]
            ace_atom.y = position[1]
            ace_atom.z = position[2]
        """
        ace_atom.atomic_number = rd_atom.GetAtomicNum()
        atom_list.append(ace_atom)
    atom_feature["chg"] = np.array(chg_list)
    # Make bond order matrix
    bonds = rd_molecule.GetBonds()
    bo_matrix = np.zeros((n, n))
    for bond in bonds:
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        # Safety net: _dearomatize_stuck_rings above normally leaves no AROMATIC
        # bonds, but if one survives (a ring that still would not kekulize),
        # approximate it as single order rather than KeyError-ing and aborting
        # generation with "failed to generate any conformers".
        bond_order = bond_types.get(bond.GetBondType(), 1)
        bo_matrix[begin][end] = bo_matrix[end][begin] = bond_order
    ace_molecule = chem.Molecule()
    ace_molecule.atom_list = atom_list
    ace_molecule.bo_matrix = bo_matrix
    ace_molecule.atom_feature = atom_feature
    return ace_molecule


def get_total_atom_list_from_molecule_list(molecule_list):
    """Returns total_atom_list by merging atom_lists in molecule_list, indexing is also."""
    atom_list = []
    cnt = 0
    for molecule in molecule_list:
        adding_atom_list = molecule.atom_list
        for atom in adding_atom_list:
            neighbors_info = atom.neighbors_info
            if neighbors_info is None:
                print("neighbors are not prepared!!! IM gneneration...")
            else:
                new_neighbors_info = dict()
                for neighbor_atom_idx in neighbors_info:
                    new_neighbors_info[neighbor_atom_idx + cnt] = neighbors_info[neighbor_atom_idx]
                neighbors_info = new_neighbors_info
        atom_list += adding_atom_list
        cnt += len(adding_atom_list)
    return atom_list


def get_block_diagonal_adj_by_distance(molecule_list, coeff=1.10):
    """From molecule_list, generate adjacency as block diagonal matrix by appending molecules."""
    total_atom = 0
    block_len = []
    for molecule in molecule_list:
        total_atom += len(molecule.atom_list)
        block_len.append(len(molecule.atom_list))
    block_diagonal_adj = np.zeros((total_atom, total_atom))
    index = 0
    for molecule in molecule_list:
        adj = get_adj_matrix_from_distance(molecule, coeff)
        block_diagonal_adj[index : index + len(adj)][index : index + len(adj)] = adj[:][:]
        index += len(adj)
    return block_diagonal_adj, block_len


def get_condensed_adj(molecule_list):
    """From molecule list, it returns adjacency of block diagonal matrix without considering."""


######### Not need to be ready... Just directly frag.py
'''
def append_neighbor(atom,bonded_atom_idx,bond_order):
    """ Modifies add new neighbor of class atom
    Args:
        |  atom (<class 'Atom'>): class instance of atom
        |  bonded_atom_idx(int): The index of a new atom that 'atom' is bonded to
        |  bond_order(int): bond_order
    Returns:
        | No return! It modifies atom!
    """
    if atom.neighbors_info == None:
        atom.neighbors_info = dict()
        atom.neighbors_info[bonded_atom_idx] = bond_order
    else:
        atom.neighbors_info[bonded_atom_idx] = bond_order

'''


def remove_neighbor(atom, bonded_atom_idx, bond_order):
    """Modifies add new neighbor of class atom."""
    if atom.neighbors_info is None:
        print("Wrong input!")
    else:
        del atom.neighbors_info[bonded_atom_idx]


# ########## From here, we have some operations that deal with multiple objects: intersection,
# different set, union,


def get_desired_set_from_two_molecule_list(molecule_list1, molecule_list2, desired_set):
    """Returns obtainable sets on molecule_lists."""
    set_list = dict()
    n1 = len(molecule_list1)
    n2 = len(molecule_list2)
    molecule_idx_list1 = set(range(n1))
    molecule_idx_list2 = set(range(n2))
    intersection_list = [
        set(),
        set(),
    ]  # Index for molecule_list1, molecule_list2 respectively
    # For every loop, obtain those functions
    for i in molecule_idx_list1:
        molecule1 = molecule_list1[i]
        for j in molecule_idx_list2:
            molecule2 = molecule_list2[j]
            if molecule1.is_same_molecule(molecule2):
                intersection_list[0].add(i)
                intersection_list[1].add(j)
    difference_set1 = molecule_idx_list1 - intersection_list[0]  # A - (A&B)
    difference_set2 = molecule_idx_list2 - intersection_list[1]  # B - (A&B)
    molecule_intersection_list = list(map(lambda x: molecule_list1[x], intersection_list[0]))
    molecule_union_list = list(map(lambda x: molecule_list1[x], molecule_idx_list1)) + list(
        map(lambda x: molecule_list2[x], difference_set2)
    )  # A|B = A | (B-A)
    molecule_difference_list1 = list(map(lambda x: molecule_list1[x], difference_set1))
    molecule_difference_list2 = list(map(lambda x: molecule_list2[x], difference_set2))
    for desired_type in desired_set:
        if desired_type == "union":
            set_list[desired_type] = molecule_union_list
        elif desired_type == "difference":
            set_list["1-2"] = molecule_difference_list1
            set_list["2-1"] = molecule_difference_list2

        elif desired_type == "intersection":
            set_list[desired_type] = molecule_intersection_list
    return set_list


def get_atom_list_from_z_list(z_list):
    """Return the atom list from z list."""
    atom_list = []
    for atomic_number in z_list:
        atom = chem.Atom(atomic_number)
        atom_list.append(atom)
    return atom_list


def get_atom_list_from_element_list(element_list):
    """Return the atom list from element list."""
    atom_list = []
    for element in element_list:
        atom = chem.Atom(element)
        atom_list.append(atom)
    return atom_list


def get_z_list_from_atom_list(atom_list):
    """Return the z list from atom list."""
    z_list = []
    for atom in atom_list:
        z_list.append(atom.get_atomic_number())
    return z_list


def get_element_list_from_atom_list(atom_list):
    """Return the element list from atom list."""
    element_list = []
    for atom in atom_list:
        element_list.append(atom.get_element())
    return element_list


def copy_atom_list(atom_list, include_geometry=False):
    """Copy the atom list."""
    new_atom_list = []
    for atom in atom_list:
        new_atom = copy.deepcopy(atom)
        new_atom_list.append(new_atom)
    return new_atom_list


def copy_molecule(molecule, shallow=True, include_geometry=False):
    """Copy the molecule."""
    new_molecule = chem.Molecule()
    if molecule.bo_matrix is not None:
        new_molecule.bo_matrix = np.copy(molecule.bo_matrix)
        new_molecule.adj_matrix = np.copy(molecule.adj_matrix)
    elif molecule.adj_matrix is not None:
        molecule.adj_matrix = np.copy(molecule.adj_matrix)
    new_molecule.chg = molecule.chg
    if shallow:
        if "chg" in molecule.atom_feature:
            new_molecule.atom_feature["chg"] = copy.deepcopy(molecule.atom_feature["chg"])
    else:
        new_molecule.atom_feature = copy.deepcopy(molecule.atom_feature)
    new_molecule.atom_list = copy_atom_list(molecule.atom_list, include_geometry)
    return new_molecule


def copy_molecule_list(molecule_list, include_geometry=False):
    """Copy the molecule list."""
    new_molecule_list = []
    for molecule in molecule_list:
        new_molecule_list.append(copy_molecule(molecule))
    return new_molecule_list


def copy_intermediate(intermediate, include_geometry=False):
    """Copy the intermediate."""
    atom_list = intermediate.atom_list
    adj_matrix = intermediate.get_matrix("adj")
    bo_matrix = intermediate.get_matrix("bo")
    chg_list = intermediate.get_chg_list()
    new_intermediate = chem.Intermediate((atom_list, adj_matrix, bo_matrix, chg_list))
    return new_intermediate


def compare_atom_list(atom_list1, atom_list2):
    """Compare atom list."""
    n = len(atom_list1)
    m = len(atom_list2)
    element_list1 = []
    element_list2 = []
    if n != m:
        return False
    else:
        for i in range(n):
            element_list1.append(atom_list1[i].get_element())
            element_list2.append(atom_list2[i].get_element())
            if not atom_list1[i].is_same_atom(atom_list2[i]):
                return False
        return True


def compare_molecule_list(molecule_list1, molecule_list2):
    """Compare molecule list."""
    n = len(molecule_list1)
    m = len(molecule_list2)
    if n != m:
        return False
    for i in range(n):
        molecule1 = molecule_list1[i]
        appeared = False
        for j in range(m):
            molecule2 = molecule_list2[j]
            if molecule1.is_same_molecule(molecule2):
                appeared = True
                break
        if not appeared:
            return False
    return True


def get_molecule_list_without_repitition(molecule_list):
    """Return the molecule list without repitition."""
    new_molecule_list = []
    for molecule in molecule_list:
        new_molecule = True
        for molecule_prime in new_molecule_list:
            if molecule.is_same_molecule(molecule_prime):
                new_molecule = False
        if new_molecule:
            new_molecule_list.append(molecule)
    return new_molecule_list


def get_permuted_molecule(molecule, permutation):
    """Return the permuted molecule."""
    atom_list = molecule.atom_list
    bo_matrix = molecule.bo_matrix
    adj_matrix = molecule.adj_matrix
    atom_feature = molecule.atom_feature
    new_molecule = chem.Molecule()
    if bo_matrix is not None:
        bo_matrix = get_permuted_matrix(bo_matrix, permutation)
        adj_matrix = np.where(bo_matrix > 0, 1, 0)
    elif adj_matrix is not None:
        adj_matrix = get_permuted_matrix(adj_matrix, permutation)
    new_molecule.adj_matrix = adj_matrix
    new_molecule.bo_matrix = bo_matrix
    new_molecule.atom_list = get_permuted_atom_list(atom_list, permutation)
    new_molecule.atom_feature = get_permuted_atom_feature(atom_feature, permutation)
    return new_molecule


def get_permuted_atom_list(atom_list, permutation):
    """Return the permuted atom list."""
    n = len(permutation)
    new_atom_list = [None] * n
    for i in range(n):
        new_atom_list[permutation[i]] = atom_list[i]
    return new_atom_list


def get_permuted_atom_feature(atom_feature, permutation):
    """Return the permuted atom feature."""
    new_atom_feature = dict()
    n = len(permutation)
    if atom_feature is None:
        return atom_feature
    for feature in atom_feature:
        feature_value = atom_feature[feature]
        new_feature_value = copy.deepcopy(feature_value)
        for i in range(n):
            value = permutation[i]
            new_feature_value[value] = feature_value[i]
        new_atom_feature[feature] = new_feature_value
    return new_atom_feature


def get_permuted_matrix(matrix, permutation):
    """Return the permuted matrix."""
    n = len(matrix)
    permuted_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            permuted_matrix[permutation[i]][permutation[j]] = matrix[i][j]
    return permuted_matrix


def get_adj_matrix_from_distance(molecule, coeff=1.10, criteria=1.0):
    """Returns adj_matrix from 3d coordinate of given molecule.

    It recognizes bond between two atoms, if the sum of radius * coeff is less than distance between
    two atoms.

    :param coeff(float):
        criteria for recognizing bond. If criteria gets higher, more and more bonds are generated
        between atoms,
        since criteria distance for bond distance gets higher.
        Appropriate criteria value is between 0.8 ~ 1.3, here we set default value as 1.10

    :return adj(pyclass 'numpy.ndarray'):
        connectivity matrix between atoms

    """
    MetalElements = [
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Y",
        "Zr",
        "Nb",
        "Mo",
        "Tc",
        "Ru",
        "Rh",
        "Pd",
        "Ag",
        "Cd",
        "Lu",
        "Hf",
        "Ta",
        "W",
        "Re",
        "Os",
        "Ir",
        "Pt",
        "Au",
        "Hg",
    ]

    atom_list = molecule.atom_list
    n = len(atom_list)
    radius_list = molecule.get_radius_list()
    radius_matrix_flatten = np.repeat(radius_list, n)
    radius_matrix = radius_matrix_flatten.reshape((n, n))
    radius_sum_matrix = radius_matrix + radius_matrix.T
    coordinate_list = molecule.get_coordinate_list()
    distance_matrix = spatial.distance_matrix(coordinate_list, coordinate_list)
    ratio_matrix = distance_matrix / radius_sum_matrix
    adj1 = np.where(distance_matrix < criteria, 1, 0)
    adj2 = np.where(ratio_matrix < coeff, 1, 0)
    adj3 = np.zeros((n, n))
    for i in range(n):
        atom = atom_list[i]
        a = atom.element
        if a in MetalElements:
            dist = radius_list[i]
            max_dist = radius_list[i] + 1.0
            for j in range(i + 1, n):
                b = atom_list[j].element
                if b != "H" and distance_matrix[i][j] < max_dist:
                    adj3[i][j] = adj3[j][i] = 1
                elif b == "H" and distance_matrix[i][j] < dist + 0.4:
                    adj3[i][j] = adj3[j][i] = 1

    adj = np.where(adj1 + adj2 + adj3 > 0, 1, 0)
    # adj = np.where(adj1+adj2>0,1,0)
    # adj = adj2
    np.fill_diagonal(adj, 0)
    return adj


def get_bo_matrix_from_adj_matrix(
    molecule, chg=None, method="SumofFragments", obtain_all_resonance=False
):
    """Returns bo_matrix from adj_matrix stored in pyclass 'Molecule'.

    :param chg(int):
        total charge of the molecule

    :param method(str):
        Mainly 'SumofFragments' and 'Ionic' is possible
        'SumofFragments' use user defined fragment charge to evaluate bo_matrix
        (also uses some chemical heuristics)
        'Ionic' method uses chemical heuristics to evaluate bo_matrix

    :param obtain_all_resonance(boolean):
        If True, it returns multiple possible resonance structures, therefore return as list of
        numpy.ndarray
        (Normally this function is not used)

    :return bo_matrix(pyclass 'numpy.ndarray' or list of pyclass 'numpy.ndarray'):
        possible bo_matrices for obtain_all_resonance=True, otherwise just single bo_matrix
    """
    if chg is None:
        chg = molecule.get_chg()
        if chg is None:
            print("Total charge is not specified! Provide charge information!")
            return None
    # The legacy fragment-based perceiver (``frag.AdjtoBO``) is not vendored in this
    # subset; route through the PuLP charge/bond-order solver, which is the engine's
    # active backend (see ``get_chg_list_and_bo_matrix_pulp``).
    _, bo_matrix = compute_pulp.compute_chg_and_bo(
        molecule, chg, resolve=True, cleanUp=True, HalogenConstraint=True
    )
    return [bo_matrix] if obtain_all_resonance else bo_matrix


def get_chg_list_from_bo_matrix(molecule, chg, bo_matrix, method="SumofFragments"):
    """Returns chg_list from a given bo_matrix stored in pyclass 'Molecule'.

    :param chg(int):
        total charge of the molecule

    :param bo_matrix(pyclass 'numpy.ndarray'):
        Possible bond order matrix of given molecule

    :param method(str):
        Mainly 'SumofFragments' and 'Ionic' is possible
        'SumofFragments' use user defined fragment charge to evaluate bo_matrix
        (also uses some chemical heuristics)
        'Ionic' method uses chemical heuristics to evaluate bo_matrix

    :return chg_list(pyclass 'numpy.ndarray' (1D-array)):
        formal charge of each atom
    """
    # Route through the PuLP solver (the ``frag.getFC`` perceiver is not vendored
    # in this subset); the solver recomputes charges consistently with the bonds.
    chg_list, _ = compute_pulp.compute_chg_and_bo(
        molecule, chg, resolve=True, cleanUp=True, HalogenConstraint=True
    )
    return np.array(chg_list)


def get_chg_and_bo(molecule, chg=None, method="SumofFragments"):
    """Return the chg and bo."""
    if chg is None:
        chg = molecule.get_chg()
        if chg is None:
            print("Total charge is not specified! Provide charge information!")
            return None, None
    # Unified charge/bond-order perception via the PuLP solver (the legacy
    # ``frag``/``compute_scipy`` backends are not vendored in this subset).
    chg_list, bo_matrix = compute_pulp.compute_chg_and_bo(
        molecule, chg, resolve=True, cleanUp=True, HalogenConstraint=True
    )
    return chg_list, bo_matrix


def get_chg_list_and_bo_matrix_pulp(molecule, chg=None, **kwargs):
    """Return the chg list and bo matrix pulp."""
    return compute_pulp.compute_chg_and_bo(
        molecule, chg, resolve=True, cleanUp=True, HalogenConstraint=True, **kwargs
    )


def get_ace_mol_from_minimal_data(minimal_data, object_type="molecule"):
    """Return the ace mol from minimal data."""
    z_list = minimal_data["z"]
    adj_matrix = minimal_data["adj"]
    bo_matrix = minimal_data["bo"]
    chg = minimal_data["chg"]
    chg_list = minimal_data["atom chg"]
    coordinate_list = minimal_data["coords"]
    if object_type == "molecule":
        molecule = chem.Molecule((z_list, adj_matrix, bo_matrix, chg_list))
    else:
        molecule = chem.Intermediate((z_list, adj_matrix, bo_matrix, chg_list))
    if chg is not None:
        molecule.chg = chg
    if coordinate_list is not None and coordinate_list[0][0] is not None:  # x value of zeroth atom
        locate_molecule(molecule, coordinate_list)
    return molecule


def add_atoms(ace_mol, new_atom_list):
    """Add the atoms."""
    atom_list = ace_mol.atom_list
    atom_feature = ace_mol.atom_feature
    new_chg_list = np.zeros((len(new_atom_list)))
    n = len(atom_list)
    m = len(new_atom_list)
    bo_matrix = ace_mol.get_matrix("bo")
    if bo_matrix is None:
        adj_matrix = ace_mol.get_matrix("adj")
        if adj_matrix is None:
            print("Cannot add atoms, since adjacency is not given")
            exit()
        new_adj_matrix = np.zeros((n + m, n + m))
        new_adj_matrix[:n, :n] = adj_matrix
        ace_mol.adj_matrix = new_adj_matrix
    else:
        bo_matrix = ace_mol.get_matrix("bo")
        new_bo_matrix = np.zeros((n + m, n + m))
        new_bo_matrix[:n, :n] = bo_matrix
        ace_mol.bo_matrix = new_bo_matrix
    for key in atom_feature:
        if key != "chg":
            atom_feature[key] = None
        else:
            atom_feature[key] = np.concatenate((atom_feature[key], new_chg_list), axis=0)
    for atom_type in new_atom_list:
        new_atom = chem.Atom(atom_type)
        atom_list.append(new_atom)


def add_bonds(ace_mol, bond_list):
    """Add the bonds."""
    bo_matrix = ace_mol.get_matrix("bo")
    if bo_matrix is None:
        print("impossible to add bond, since bond order is not given!")
        exit()
    else:
        for bond in bond_list:
            start = bond[0]
            end = bond[1]
            bond_order = bond[2]
            bo_matrix[start][end] = bo_matrix[end][start] = bond_order


def add_atoms_with_bonds(ace_mol, add_info_list):
    """Add the atoms with bonds."""
    n = len(ace_mol.atom_list)
    new_bond_list = []
    new_atom_list = []
    cnt = 0
    for add_info in add_info_list:
        atom_type = add_info[0]
        bond_info = add_info[1]
        new_atom_list.append(atom_type)
        for bond in bond_info:
            start = n + cnt
            end = bond[0]
            bond_order = bond[1]
            new_bond_list.append((start, end, bond_order))
        cnt += 1
    add_atoms(ace_mol, new_atom_list)
    add_bonds(ace_mol, new_bond_list)


def molecule_to_intermediate(ace_mol):
    """Molecule to intermediate."""
    intermediate = chem.Intermediate()
    intermediate.atom_list = ace_mol.atom_list
    intermediate.bo_matrix = ace_mol.bo_matrix
    intermediate.adj_matrix = ace_mol.adj_matrix
    intermediate.atom_feature = ace_mol.atom_feature
    return intermediate


def molecule_to_ase_atoms(molecule):
    """Molecule to ase atoms."""
    from ase import Atoms

    element_list = molecule.get_element_list()
    symbols = "".join(element_list)
    positions = molecule.get_coordinate_list()
    # bo_matrix = molecule.get_matrix('bo')
    # also give original charge information
    # charge_list = molecule.get_chg_list().tolist()
    ase_atoms = Atoms(
        symbols,
        positions=positions,
        charges=None,
    )
    # ase_atoms.set_initial_charges(charge_list)
    return ase_atoms


def read_geometries(directory):
    """Read the geometries."""
    conformers = []
    with open(directory) as f:
        while True:
            try:
                atom_num = int(f.readline().strip())
                energy = float(f.readline().strip())
                atom_list = []
                for i in range(atom_num):
                    line = f.readline().strip().split()
                    element = line[0]
                    atom = chem.Atom(element)
                    atom.x = float(line[1])
                    atom.y = float(line[2])
                    atom.z = float(line[3])
                    atom_list.append(atom)
                molecule = chem.Molecule()
                molecule.atom_list = atom_list
                molecule.energy = energy
                conformers.append(molecule)
            except Exception:
                break
    return conformers


def check_atom_validity(group, bo, chg, octet=4):
    # lone pair inequality
    """Check the atom validity."""
    if group - bo - chg < 0:
        return False
    # Octet rule inequality
    if group + bo - chg > 2 * octet:
        return False
    return True


def get_molecule_group(adj_matrix, index=0):
    """Return the molecule group."""
    current_list = set([index])
    total_list = set([index])
    while len(current_list) > 0:
        new_current_list = set([])
        for i in current_list:
            neighbor_list = np.where(adj_matrix[i] > 0)[0].tolist()
            new_current_list = new_current_list | set(neighbor_list)
        current_list = new_current_list - total_list
        total_list = total_list | new_current_list
    return total_list


def group_molecules(adj_matrix):
    """Group molecules."""
    n = len(adj_matrix)
    all_indices = set(range(n))
    groups = []
    index = 0
    while len(all_indices) > 0:
        indices = get_molecule_group(adj_matrix, index)
        all_indices = all_indices - indices
        groups.append(list(indices))
        if len(all_indices) > 0:
            index = min(all_indices)
        else:
            break
    return groups


def check_geometry(coordinate_list, criteria=0.5):
    """Check the geometry."""
    distance_matrix = spatial.distance_matrix(coordinate_list, coordinate_list)
    np.fill_diagonal(distance_matrix, 100)
    check_distance_matrix = np.where(distance_matrix < criteria, 1, 0)
    return np.sum(check_distance_matrix) < 1


def get_rmsd(molecule1, molecule2):
    """Return the rmsd."""
    n = len(molecule1.atom_list)
    if len(molecule2.atom_list) != n:
        print("Cannot calculate RMSD!!!")
        return None
    coordinate_list1 = np.array(molecule1.get_coordinate_list())
    coordinate_list2 = np.array(molecule2.get_coordinate_list())
    rmsd = np.sqrt(np.sum((coordinate_list1 - coordinate_list2) ** 2) / n)
    return rmsd


def is_same_connectivity(
    original_molecule, new_molecule, max_coeff=1.3, min_coeff=0.95, space=0.01
):
    """Return whether same connectivity."""
    coeff = min_coeff
    is_same = False
    if len(original_molecule.atom_list) == 1:
        return True, 0.95
    while coeff < max_coeff:
        adj_matrix = get_adj_matrix_from_distance(new_molecule, coeff)
        new_molecule.set_adj_matrix(adj_matrix)
        is_same = new_molecule.is_same_molecule(original_molecule, False)
        if is_same:
            break
        coeff += space
    return is_same, coeff


def minimize_rmsd(reference_molecule, changing_molecule):
    """Minimize rmsd."""
    from ase.build.rotate import minimize_rotation_and_translation

    reference_ase_atoms = molecule_to_ase_atoms(reference_molecule)
    target_ase_atoms = molecule_to_ase_atoms(changing_molecule)
    minimize_rotation_and_translation(reference_ase_atoms, target_ase_atoms)
    coordinate_list = target_ase_atoms.get_positions()
    locate_molecule(changing_molecule, coordinate_list)


def get_molecule_info_from_sdf(sdf_directory):
    """Return the molecule info from sdf."""
    from . import globalvars as gv

    with open(sdf_directory, "r") as f:
        lines = f.readlines()

    # First, get the number of atoms
    n_atoms = int(lines[3].split()[0])
    n_bonds = int(lines[3].split()[1])

    z_list = []
    coords = []
    adj_matrix = np.zeros((n_atoms, n_atoms))
    chg_list = np.zeros(n_atoms)
    metal_index = None

    for i in range(n_atoms):
        line = lines[4 + i].split()
        x = float(line[0])
        y = float(line[1])
        z = float(line[2])
        element = line[3]
        coords.append([x, y, z])
        z_list.append(chem.Atom(element).get_atomic_number())
        if element in gv.metal:
            metal_index = i
    for i in range(n_bonds):
        line = lines[4 + n_atoms + i].split()
        s, e = line[0], line[1]
        s_ = int(s.strip())
        if s_ > n_atoms:
            tmp = s
            s = tmp[:3].strip()
            e = tmp[3:].strip()
        if not s.isdigit() or not e.isdigit():
            print("WRONG SDF; Error occurs during parsing bond block ...")
        s = int(s) - 1
        e = int(e) - 1
        adj_matrix[s][e] = 1
        adj_matrix[e][s] = 1
    chg_line = lines[4 + n_atoms + n_bonds]
    if "CHG" in chg_line:
        chg_line = chg_line.strip().split()
        n_chg = int(chg_line[2])
        for i in range(n_chg):
            idx = int(chg_line[3 + 2 * i]) - 1
            chg = int(chg_line[4 + 2 * i])
            chg_list[idx] = chg

    coords = np.array(coords)

    return z_list, coords, adj_matrix, chg_list, metal_index
