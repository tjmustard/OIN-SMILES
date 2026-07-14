import logging

from rdkit import Chem

logger = logging.getLogger(__name__)


def get_ligand_from_smiles(mapped_smiles):
    """Return the ligand from smiles."""
    from . import process

    rd_mol = Chem.MolFromSmiles(mapped_smiles, sanitize=False)
    rd_mol.UpdatePropertyCache(strict=False)
    # Recover C=C (cis/trans) stereo from the parsed bond directions so it can be
    # enforced deterministically at embed time. MolFromSmiles(sanitize=False)
    # leaves the '/' and '\' markers as BondDir only; SetBondStereoFromDirections
    # materializes STEREOCIS/STEREOTRANS plus the two stereo reference atoms. The
    # ace_mol is otherwise stereo-blind (get_ace_mol_from_rd_mol copies only
    # atomic number / charge / bond order), so without this the embed picks a
    # random dihedral. Indices are ligand-local and survive AddHs (heavy-atom
    # order is preserved), matching the ace_mol atom order built below.
    stereo_bonds = []
    try:
        Chem.SetBondStereoFromDirections(rd_mol)
        for b in rd_mol.GetBonds():
            if (
                b.GetBondType() == Chem.BondType.DOUBLE
                and b.GetStereo() != Chem.BondStereo.STEREONONE
            ):
                sa = list(b.GetStereoAtoms())
                if len(sa) == 2:
                    stereo_bonds.append(
                        (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetStereo(), sa[0], sa[1])
                    )
    except Exception:
        stereo_bonds = []
    # Restrict enforcement to geometrically-free double bonds. A double bond whose
    # atoms coordinate the metal (atom-map number > 0) or neighbour a coordinating
    # atom is rigidly fixed by chelation/coordination -- its E/Z is already
    # reproduced by the ligand topology, so a distance-geometry stereo constraint
    # there is unnecessary AND over-constrains the metal embed (measured to cut the
    # conformer yield of metal-bound imine chelates, e.g. AGULIX 9/9 -> 3/9). Keep
    # only bonds clear of the coordination sphere; the target case -- a pendant,
    # freely-rotatable alkene -- is unaffected.
    #
    # This proxy is deliberately broader than the encoder's chelate-ring test in
    # core/translator._clear_chelate_locked_bond_stereo, so the encoder can emit an
    # E/Z the generator does not enforce (a free C=N hanging off a monodentate
    # donor: AFECIZ, XIZXAG -- both still round-trip-fail on that bond). Narrowing
    # it to the ring test is correct in principle but detonates a latent bug in
    # _apply_double_bond_stereo, which force-sets the bond back to DOUBLE without
    # touching formal charges: on AFECIZ the PuLP charge assignment wants that C=N
    # single with a charged N, so SanitizeMol then rejects a 4-valent N and every
    # ff_clean fails (565s -> full attempt budget). Fix that first, then narrow.
    if stereo_bonds:
        near_donor = {a.GetIdx() for a in rd_mol.GetAtoms() if a.GetAtomMapNum() > 0}
        for d in list(near_donor):
            for nb in rd_mol.GetAtomWithIdx(d).GetNeighbors():
                near_donor.add(nb.GetIdx())
        stereo_bonds = [
            sb for sb in stereo_bonds if sb[0] not in near_donor and sb[1] not in near_donor
        ]
    rd_mol = Chem.AddHs(rd_mol, explicitOnly=False, addCoords=False)
    # Recover sp3 atom chirality for the same reason as the C=C stereo above: the
    # '@'/'@@' survive MolFromSmiles(sanitize=False) as chiral tags, but the
    # ace_mol drops them, and the embed then returns a random enantiomer for every
    # stereocentre. Captured AFTER AddHs so the hydrogen is an explicit neighbour
    # and the indices line up with the ace_mol atom order built below.
    #
    # Skip metal-binding centres (atom map > 0): the complex adds a metal bond, so
    # their neighbour set -- and hence the meaning of the tag -- changes. Zone-A
    # donor stereo is re-asserted from the template by build_contract_mol instead.
    chiral_centers = []
    for atom in rd_mol.GetAtoms():
        tag = atom.GetChiralTag()
        if tag not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            continue
        if atom.GetAtomMapNum() > 0:
            continue
        nbrs = tuple(b.GetOtherAtomIdx(atom.GetIdx()) for b in atom.GetBonds())
        if len(nbrs) != 4:
            continue
        chiral_centers.append((atom.GetIdx(), nbrs, tag))
    ace_mol = process.get_ace_mol_from_rd_mol(rd_mol)
    ace_mol.stereo_bonds = stereo_bonds
    ace_mol.chiral_centers = chiral_centers
    n = len(ace_mol.atom_list)
    binding_infos = dict()
    for i in range(n):
        rd_atom = rd_mol.GetAtomWithIdx(i)
        mapping = rd_atom.GetAtomMapNum()
        if mapping > 0:
            if mapping in binding_infos:
                binding_infos[mapping].append(i)
            else:
                binding_infos[mapping] = [i]
    final_binding_infos = []
    for position in binding_infos:
        final_binding_infos.append([binding_infos[position], position])

    ligand = Ligand(ace_mol, final_binding_infos)
    return ligand


class Ligand:
    """Ligand."""

    def __init__(self, molecule, binding_infos):
        """Initialize the Ligand."""
        self.molecule = molecule  # chem.Molecule
        self.binding_infos = binding_infos  # [[[int],int]]

    def get_smiles(self):
        """Return the smiles."""
        return self.molecule.get_smiles("ace")

    def get_denticity(self):
        """Return the denticity."""
        return len(self.binding_infos)

    def update(self):
        """Update."""
        self.coord_list = self.molecule.get_coordinate_list()

    def get_adj_matrix(self):
        """Return the adj matrix."""
        adj_matrix = self.molecule.adj_matrix
        if adj_matrix is None:
            adj_matrix = self.molecule.get_adj_matrix()
        return adj_matrix

    def print_coordinate_list(self):
        """Print the coordinate list."""
        coordinate_list = self.coordinate_list()
        if coordinate_list is None:
            coordinate_list = self.get_coordinate_list()
        n = len(coordinate_list)

        logger.debug(n)
        logger.debug("")
        for i in range(n):
            symbol, x, y, z = self.coordinate_list[i]
            print_x = f"{x:.6f}"
            print_y = f"{y:.6f}"
            print_z = f"{z:.6f}"
            logger.debug(f"{symbol} {print_x} {print_y} {print_z}")
        logger.debug("")

    def copy(self):
        """Copy."""
        import copy

        molecule = self.molecule.copy()
        molecule.adj_matrix = self.molecule.get_adj_matrix()
        binding_infos = copy.deepcopy(self.binding_infos)
        new_ligand = Ligand(molecule, binding_infos)
        return new_ligand

    def __str__(self):
        """Return ``str(self)``."""
        pass
