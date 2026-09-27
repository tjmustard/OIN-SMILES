import copy
import logging

import numpy as np

from . import chem, embed, process
from . import globalvars as gv
from . import ligand as ligand_module

logger = logging.getLogger(__name__)


class Geometry:
    """Geometry."""

    def __init__(self, geometry_name):
        """Initialize the Geometry."""
        self.geometry_name = geometry_name  # str
        self.direction_vector = (
            gv.known_geometries_vector_dict[geometry_name]
            if geometry_name in gv.known_geometries_vector_dict.keys()
            else []
        )
        self.permutations = (
            gv.known_geometries_permutation_dict[geometry_name]
            if geometry_name in gv.known_geometries_permutation_dict.keys()
            else []
        )

    def get_steric_number(self):
        """Return the steric number."""
        return len(self.direction_vector)


class MetalComplex:
    """Metal complex."""

    def __init__(self, geometry_name, center_atom, ligands, chg, multiplicity):
        """Initialize the Metal complex."""
        self.geometry_type = Geometry(geometry_name)  # Geometry Object
        self.center_atom = center_atom  # chem.Atom
        self.ligands = ligands  # [Ligand Object]

        self.atom_indices_for_each_ligand = []  # [[int]]
        self.metal_index = None  # int
        self.adj_matrix = None  # np.array
        self.bo_matrix = None  # np.array
        self.is_actinide = False

        num_atom = 0
        if center_atom is not None:
            num_atom += 1
        for ligand in ligands:
            num_atom += len(ligand.molecule.atom_list)

        self.num_atom = num_atom
        self.chg = chg
        self.multiplicity = multiplicity
        self.name = None

    def get_atom_indices_for_each_ligand(self):
        """Return the atom indices for each ligand."""
        if len(self.atom_indices_for_each_ligand) > 0:
            return self.atom_indices_for_each_ligand
        else:
            ligands = self.ligands
            atom_indices_for_each_ligand = []
            n = 1
            for ligand in ligands:
                m = len(ligand.molecule.atom_list)
                atom_indices = [i for i in range(n, n + m)]
                atom_indices_for_each_ligand.append(atom_indices)
                n += m
            self.atom_indices_for_each_ligand = atom_indices_for_each_ligand
            return atom_indices_for_each_ligand

    def get_binding_groups(self):
        """Return the binding groups."""
        ligands = self.ligands
        atom_indices_for_each_ligand = self.get_atom_indices_for_each_ligand()
        steric_number = self.geometry_type.get_steric_number()
        n = len(ligands)
        binding_groups = [None] * steric_number
        for i in range(n):
            ligand = ligands[i]
            atom_indices = atom_indices_for_each_ligand[i]
            binding_infos = ligand.binding_infos
            for binding_info in binding_infos:
                binding_indices, binding_site = binding_info
                if binding_site is None:
                    logger.debug("Coordination Information is not determined !!!")
                    return []
                binding_groups[binding_site - 1] = [
                    atom_indices[index] for index in binding_indices
                ]

        return binding_groups

    def get_adj_matrix(self):
        """Return the adj matrix."""
        metal_index = self.metal_index
        ligands = self.ligands
        atom_indices_for_each_ligand = self.get_atom_indices_for_each_ligand()
        m = len(ligands)

        n = 1
        for ligand in ligands:
            n += len(ligand.molecule.atom_list)

        adj_matrix = np.zeros((n, n))

        for i in range(m):
            ligand = ligands[i]
            atom_indices = atom_indices_for_each_ligand[i]
            mol_adj = ligand.molecule.get_adj_matrix()
            for j in range(len(mol_adj)):
                for k in range(len(mol_adj)):
                    start = atom_indices[j]
                    end = atom_indices[k]
                    adj_matrix[start][end] = mol_adj[j][k]
                    adj_matrix[end][start] = mol_adj[k][j]
        binding_groups = self.get_binding_groups()
        for group in binding_groups:
            if group is None:
                # Unfilled coordination slot (fewer ligands than the geometry has
                # vertices) -- no metal-donor bond to add. Skipping keeps every
                # fully-coordinated complex byte-identical while letting an
                # under-coordinated one (e.g. a lone eta-arene in a 4-slot shell)
                # build instead of crashing on ``for i in None``.
                continue
            for i in group:
                adj_matrix[metal_index][i] = 1
                adj_matrix[i][metal_index] = 1
        return adj_matrix

    def get_atom_list(self):
        """Return the atom list."""
        metal_index = self.metal_index
        ligands = self.ligands
        atom_indices_for_each_ligand = self.get_atom_indices_for_each_ligand()

        n = 1
        for ligand in ligands:
            n += len(ligand.molecule.atom_list)

        atom_list = [None] * n
        atom_list[metal_index] = self.center_atom.copy()

        for i in range(len(ligands)):
            ligand = ligands[i]
            atom_indices = atom_indices_for_each_ligand[i]
            ligand_atom_list = ligand.molecule.atom_list
            for j in range(len(ligand_atom_list)):
                atom_list[atom_indices[j]] = ligand_atom_list[j].copy()

        return atom_list

    def get_molecule(self):
        """Return the molecule."""
        molecule = chem.Molecule()
        chg = self.chg
        mult = self.multiplicity

        atom_list = self.get_atom_list()
        adj_matrix = self.get_adj_matrix()

        molecule.atom_list = atom_list
        molecule.atom_feature = dict()
        molecule.adj_matrix = adj_matrix
        molecule.chg = chg
        molecule.multiplicity = mult
        if hasattr(self, "energy"):
            molecule.energy = self.energy

        return molecule

    def get_position(self):
        """Return the position."""
        atom_list = self.get_atom_list()
        positions = [[atom.x, atom.y, atom.z] for atom in atom_list]
        return np.array(positions)

    def set_position(self, positions):
        # CAUTION: positions should be in the order of atom_list
        """Set the position."""
        center_atom = self.center_atom
        ligands = self.ligands
        atom_indices_for_each_ligand = self.get_atom_indices_for_each_ligand()
        atom_list = self.get_atom_list()
        n = len(atom_list)
        if len(positions) != n:
            logger.debug("Number of atoms does not match ...")
            return

        process.locate_atom(center_atom, positions[self.metal_index])
        for i, atom_indices in enumerate(atom_indices_for_each_ligand):
            ligand = ligands[i]
            for j in range(len(atom_indices)):
                process.locate_atom(ligand.molecule.atom_list[j], positions[atom_indices[j]])

    def copy(self):
        """Copy."""
        geometry_name = self.geometry_type.geometry_name
        center_atom = self.center_atom.copy()
        ligands = [ligand.copy() for ligand in self.ligands]
        chg = self.chg
        multiplicity = self.multiplicity
        atom_indices_for_each_ligand = copy.deepcopy(self.atom_indices_for_each_ligand)

        new_complex = MetalComplex(geometry_name, center_atom, ligands, chg, multiplicity)
        new_complex.metal_index = self.metal_index
        new_complex.atom_indices_for_each_ligand = atom_indices_for_each_ligand
        new_complex.name = self.name

        return new_complex

    def get_stereoisomers(self):
        """Return the stereoisomers."""
        geometry_type = self.geometry_type
        permutations = geometry_type.permutations
        isomers = []
        if len(permutations) == 0:
            logger.debug("Not supported geometry type ...")
            return isomers
        else:
            for permutation in permutations:
                isomer = self.copy()
                ligands = isomer.ligands
                for i, ligand in enumerate(ligands):
                    for j, binding_info in enumerate(ligand.binding_infos):
                        _, _ = binding_info
                        ligand.binding_infos[j] = (binding_info[0], permutation[i])
                isomers.append(isomer)
        return isomers

    def get_embedding(self, num_conformer=10, d_criteria=0.5, align=True):
        """Return the embedding."""
        options = [0, 1]
        num_conf_per_option = int((num_conformer + 1) / 2)
        if num_conf_per_option == 0:
            num_conf_per_option = 1
        if num_conf_per_option > 1:
            scale_size = min(0.1, 0.4 / (num_conf_per_option - 1))
            start = max(0.8, 1 - scale_size * (num_conf_per_option - 1) / 2)
            end = min(1.2, 1 + scale_size * (num_conf_per_option - 1) / 2)
            scales = np.arange(start, end, scale_size)
        else:
            scales = np.array([1.0])
        candidate_positions = []
        for option in options:
            for scale in scales:
                if True:
                    positions = embed.get_embedding(self, scale, option, align=align)
                    if positions is not None:
                        candidate_positions.append(positions)
                else:
                    continue

            if len(candidate_positions) == num_conformer:
                break
        if len(candidate_positions) == 0:
            logger.debug("No valid embedding found ...")
            raise RuntimeError("No valid embedding found: all conformer embed attempts failed.")
        return candidate_positions

    def print_coordinate_list(self):
        """Print the coordinate list."""
        atom_list = self.get_atom_list()
        len(atom_list)
        for atom in atom_list:
            element = atom.get_element()
            coordinate = atom.get_coordinate()
            print_x = f"{coordinate[0]:>12.8f}"
            print_y = f"{coordinate[1]:>12.8f}"
            print_z = f"{coordinate[2]:>12.8f}"
            logger.debug(f"{element:<3} {print_x} {print_y} {print_z}")
        logger.debug("")

    def get_distances_from_center(self):
        """Get the distances of all atoms from the center atom."""
        metal_index = self.metal_index
        if metal_index is None:
            logger.debug("Metal index is not determined ...")
            return []
        adj_matrix = self.get_adj_matrix()
        neighbor_list = [-1] * len(adj_matrix)
        atom_set = set([metal_index])
        neighbor_list[metal_index] = 0
        distance = 1
        while len(atom_set) > 0:
            next_set = set()
            for atom in atom_set:
                for i in range(len(adj_matrix)):
                    if adj_matrix[atom][i] == 1 and neighbor_list[i] == -1:
                        neighbor_list[i] = distance
                        next_set.add(i)
            atom_set = next_set
            distance += 1
        return neighbor_list


def replace_actinide(metal_complex):
    """Replace actinide."""
    metal_atom = metal_complex.center_atom
    metal_atom = metal_atom.get_element().lower().capitalize()
    if metal_atom in gv.actinide_metal:
        metal_complex.is_actinide = True
        index = gv.actinide_metal.index(metal_atom)
        corresponding_lanthanide = gv.lanthanide_metal[index]
        metal_complex.center_atom.set_element(corresponding_lanthanide)


def group_binding_sites(binding_indices, adj_matrix):
    """Group binding sites."""
    all_set = set(binding_indices)
    group_list = []
    index = binding_indices[0]
    while len(all_set) > 0:
        molecule_indices = set([index])
        current_indices = set([index])
        while len(current_indices) > 0:
            next_indices = set()
            for i in current_indices:
                for j in binding_indices:
                    if adj_matrix[i][j] == 1:
                        next_indices.add(j)
            current_indices = next_indices - molecule_indices
            molecule_indices = molecule_indices.union(current_indices)
        all_set = all_set - molecule_indices
        group_list.append(list(molecule_indices))
        if len(all_set) > 0:
            index = list(all_set)[0]
        else:
            break
    return group_list


def construct_metal_complex(z_list, adj_matrix, geometry_name=None):
    """Construct metal complex."""
    metal_index = None
    for i, atomic_num in enumerate(z_list):
        if atomic_num in gv.metal_z_list:
            metal_index = i
            break
    if metal_index is None:
        logger.debug("Metal was not found !!!")
        raise RuntimeError("Metal atom was not found in the molecule.")

    broken_adj_matrix = np.copy(adj_matrix)
    ligands = []
    binding_sites = []

    n = len(z_list)

    # Cut all the bonds between metals and obtain binding sites ...
    for i in range(n):
        if broken_adj_matrix[metal_index][i] > 0:
            binding_sites.append(i)
            broken_adj_matrix[metal_index][i] = broken_adj_matrix[i][metal_index] = 0.0

    groups = process.group_molecules(broken_adj_matrix)
    order = 0
    atom_indices_for_each_ligand = []

    for group in groups:
        if metal_index in group:
            continue
        group.sort()
        atom_indices_for_each_ligand.append(group)
        ligand_atom_list = [chem.Atom(z_list[i]) for i in group]
        chg = 0
        if len(group) > 1:
            reduce_function = {group[i]: i for i in range(len(group))}
            index_function = np.ix_(group, group)
            ligand_adj_matrix = adj_matrix[index_function]
            # Find binding indices for given ligand
            binding_indices = [i for i in group if i in binding_sites]
            index_function = np.ix_(binding_indices, binding_indices)
            sub_adj_matrix = adj_matrix[index_function]
            ligand_binding_groups = process.group_molecules(sub_adj_matrix)
            chg = -len(ligand_binding_groups)
            order += len(ligand_binding_groups)
            binding_indices = [reduce_function[i] for i in binding_indices]
            binding_indices = group_binding_sites(binding_indices, ligand_adj_matrix)
        else:
            # Single atom binding ligand ...
            ligand_adj_matrix = np.zeros((1, 1))
            chg = -1
            binding_indices = [[0]]
            order += 1

        # Identify order for binding sites
        # Make ligand molecule
        ligand_molecule = chem.Molecule()
        ligand_molecule.atom_list = ligand_atom_list
        ligand_molecule.adj_matrix = ligand_adj_matrix
        ligand_molecule.chg = chg
        ligand_molecule.multiplicity = 1

        chg_list, bo_matrix = process.get_chg_and_bo(ligand_molecule, ligand_molecule.chg)

        ligand_molecule.set_atom_feature(chg_list, "chg")
        ligand_molecule.chg_list = chg_list
        ligand_molecule.bo_matrix = bo_matrix

        binding_infos = [[i, None] for i in binding_indices]

        ligands.append(ligand_module.Ligand(ligand_molecule, binding_infos))

    center_atom = chem.Atom(z_list[metal_index])

    chg = None
    multiplicity = None
    metal_complex = None

    if geometry_name is None:
        logger.debug("Geometry not specified")
        raise RuntimeError("Cannot build MetalComplex: geometry_name is None.")

    metal_complex = MetalComplex(geometry_name, center_atom, ligands, chg, multiplicity)
    metal_complex.metal_index = metal_index
    metal_complex.atom_indices_for_each_ligand = atom_indices_for_each_ligand

    return metal_complex


# High-spin ground-state assignment for weak-field first-row transition-metal
# complexes (see _spin_multiplicity). ``_HS_GROUP`` is the neutral-atom valence-
# electron count (group number); with an oxidation state ``ox`` the d-electron count
# is ``group - ox``, and ``_HS_UNPAIRED`` indexes the Hund's-rule high-spin unpaired
# count by d (d0..d10). ``_HS_COMMON_OX`` lists the metal's common oxidation states;
# the drawn electron-count parity selects a unique one (see _spin_multiplicity).
_HS_GROUP = {"V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11}
_HS_COMMON_OX = {
    "V": (3, 2),
    "Cr": (3, 2),
    "Mn": (2, 3),
    "Fe": (3, 2),
    "Co": (2, 3),
    "Ni": (2,),
    "Cu": (2,),
}
_HS_UNPAIRED = (0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)  # unpaired e- for high-spin d0..d10
_WEAK_FIELD_DONORS = frozenset({"F", "Cl", "Br", "I", "O"})


def _spin_multiplicity(metal_element, z_sum, chg, ligands):
    """Return the ground-state spin multiplicity for the (neutral-drawn) complex.

    Defaults to the minimum multiplicity consistent with the drawn electron-count
    parity -- the long-standing low-spin assumption, correct for closed-shell and
    strong-field complexes. For a first-row transition metal whose donor set is
    unambiguously weak-field (only halide/oxygen donors), assign the Hund's-rule
    high-spin ground state instead: the common oxidation state whose d-electron
    count is parity-consistent with the drawn electron count fixes the number of
    unpaired electrons.

    The generation m-SMILES is drawn all-neutral, so oxidation state / d-count are
    not directly recoverable. The parity constraint -- which the assigned unpaired
    count must satisfy for the ``.UHF`` an optimizer consumes -- selects a unique
    oxidation state among the consecutive common ones. Anything ambiguous (mixed or
    strong-field donors, a non-cohort metal, or no parity-consistent d-count) keeps
    the low-spin default, so this never lowers spin and leaves every strong-field /
    photocatalyst / organometallic complex unchanged.
    """
    parity = (z_sum - chg) % 2
    low_spin = parity + 1

    group = _HS_GROUP.get(metal_element)
    if group is None:
        return low_spin

    donors = [
        lig.molecule.atom_list[idx].get_element()
        for lig in ligands
        for binding_indices, _site in lig.binding_infos
        for idx in binding_indices
    ]
    if not donors or any(d not in _WEAK_FIELD_DONORS for d in donors):
        return low_spin

    for ox in _HS_COMMON_OX[metal_element]:
        d = group - ox
        if 0 <= d <= 10 and _HS_UNPAIRED[d] % 2 == parity:
            return _HS_UNPAIRED[d] + 1

    return low_spin


def get_om_from_modified_smiles(smiles):
    """Return the om from modified smiles."""
    from rdkit import Chem

    from . import ligand

    smiles_list = smiles.split("|")

    n = len(smiles_list)
    metal_atom = Chem.MolFromSmiles(smiles_list[0])
    metal_chg = Chem.GetFormalCharge(metal_atom)
    metal_atom = chem.Atom(metal_atom.GetAtomWithIdx(0).GetSymbol())
    z_sum = metal_atom.get_atomic_number()
    ligands = [ligand.get_ligand_from_smiles(smiles_list[i]) for i in range(1, n - 1)]

    ligand_chg = sum([lig.molecule.get_chg() for lig in ligands])
    z_sum += sum([np.sum(lig.molecule.get_z_list()) for lig in ligands])

    chg = ligand_chg + metal_chg
    geometry_name = smiles_list[-1]
    multiplicity = _spin_multiplicity(metal_atom.get_element(), z_sum, chg, ligands)

    metal_complex = MetalComplex(geometry_name, metal_atom, ligands, chg, multiplicity)
    metal_complex.metal_index = 0
    metal_complex.multiplicity = multiplicity
    replace_actinide(metal_complex)

    return metal_complex


def get_om_from_parsed(metal_smiles, ligand_specs, geometry_name):
    """Build a ``MetalComplex`` directly from ``ParsedOIN``-derived fragments (SL2).

    The default generation path round-trips through the ``metal|lig|...|geo``
    m-SMILES string, which cannot express eta-ring winding (the sign is recovered
    stochastically by a wide embed pool). This direct constructor takes the same
    per-fragment ``mapped_smiles`` the m-SMILES path builds -- so the load-bearing
    ligand chemistry (``get_ligand_from_smiles``: stereo recovery, AddHs, bare-donor
    H-stripping, Cp kekulization) is byte-identical -- but additionally attaches the
    OIN winding target to each haptic ligand for deterministic winding construction
    in the rigid placer.

    ``ligand_specs`` is a list (ordered exactly as the m-SMILES join) of
    ``(mapped_smiles, winding)`` where ``winding`` is ``None`` or
    ``(output_order, star_frag_idx, char)``: ``output_order[canonical_pos]`` is the
    atom's original fragment index (RDKit ``_smilesAtomOutputOrder`` from the
    canonicalizing ``MolToSmiles``), ``star_frag_idx`` is the heading atom's original
    fragment index, and ``char`` is the requested ``'>'``/``'<'``. Metal ``@SPn``
    chirality stays deferred (as on the m-SMILES path); the encoder can't reproduce
    it yet.
    """
    from rdkit import Chem

    from . import ligand as ligand_mod

    metal_atom = Chem.MolFromSmiles(metal_smiles)
    metal_chg = Chem.GetFormalCharge(metal_atom)
    metal_atom = chem.Atom(metal_atom.GetAtomWithIdx(0).GetSymbol())
    z_sum = metal_atom.get_atomic_number()

    ligands = []
    for mapped_smiles, winding in ligand_specs:
        lig = ligand_mod.get_ligand_from_smiles(mapped_smiles)
        # Attach the deterministic winding target only for a PURE haptic ligand --
        # a single binding_info whose face has >1 atom -- which is exactly the case
        # the rigid placer routes through _place_haptic. (A haptic arm inside a
        # chelate takes the chelate branch and keeps the stochastic behaviour.)
        if winding is not None and len(lig.binding_infos) == 1 and len(lig.binding_infos[0][0]) > 1:
            output_order, star_frag_idx, char = winding
            face_locals = lig.binding_infos[0][0]
            # Re-order the ring atoms into the encoder's SMILES/fragment order
            # (ascending ORIGINAL fragment index) so signed_circulation's "next
            # after star" is the same physical atom the encoder used.
            ordered = sorted(face_locals, key=lambda p: output_order[p])
            star_canon = next((p for p in face_locals if output_order[p] == star_frag_idx), None)
            if star_canon is not None:
                lig.winding = (ordered, ordered.index(star_canon), char)
        ligands.append(lig)

    ligand_chg = sum([lig.molecule.get_chg() for lig in ligands])
    z_sum += sum([np.sum(lig.molecule.get_z_list()) for lig in ligands])

    chg = ligand_chg + metal_chg
    multiplicity = _spin_multiplicity(metal_atom.get_element(), z_sum, chg, ligands)

    metal_complex = MetalComplex(geometry_name, metal_atom, ligands, chg, multiplicity)
    metal_complex.metal_index = 0
    metal_complex.multiplicity = multiplicity
    replace_actinide(metal_complex)

    return metal_complex
