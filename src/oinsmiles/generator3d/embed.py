import itertools
import logging
import random

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from rdkit.Geometry import Point3D
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation as R

from ..generation import _telemetry
from ..oin.winding import signed_circulation
from . import chem, clash, process
from .utils import ic

logger = logging.getLogger(__name__)


# get_alternative_molecule is a pure function of the (structurally fixed) complex
# and ``option`` -- it ignores ``scale`` -- and its result is consumed read-only
# (``get_rd_mol`` builds a fresh Mol per attempt). ``get_embedding`` rebuilds
# ``new_complex = metal_complex.copy()`` and recomputes it every attempt, so a
# per-generation memo removes those redundant rebuilds. The cache dict is OWNED by
# ``generate_3d_structures`` (one fresh dict per generation) and threaded in via the
# ``alt_cache`` parameter: direct callers pass None and always recompute, so there
# is NO cross-molecule staleness (a module-global keyed by option is unsafe -- a
# direct call would read another molecule's entry). Same spirit as the
# ``compute_chg_and_bo`` topology memo.
def _alt_mol_cached(new_complex, option, cache):
    """Memoized get_alternative_molecule via a caller-supplied per-generation dict."""
    if cache is not None and option in cache:
        return cache[option]
    val = get_alternative_molecule(new_complex, option)
    if cache is not None:
        cache[option] = val
    return val


def get_transition_metal_center(geometry_name):
    """Return the transition metal center."""
    tm_center = {
        "6_octahedral": 26,
        "4_tetrahedral": 30,
        "4_square_planar": 78,
    }
    return tm_center.get(geometry_name, 26)


def get_dummy_center_for_valid(geometry_name):
    # Get the dummy center to get valid ace molecule
    """Return the dummy center for valid."""
    dummy_centers = {
        "2": chem.Atom("O"),
        "3": chem.Atom("P"),
        "4": chem.Atom("Si"),
        "5": chem.Atom("P"),
        "6": chem.Atom("S"),
    }
    steric_number = geometry_name.split("_")[0]
    return dummy_centers[steric_number] if steric_number in dummy_centers else chem.Atom("Fe")


def get_dummy_center(geometry_name):
    # Get the dummy center to embed the molecule
    """Return the dummy center."""
    dummy_centers = {
        "2_linear": chem.Atom("Fe"),
        "2_bent_135": chem.Atom("S"),
        "2_bent_90": chem.Atom("S"),
        "3_t_shaped": chem.Atom("P"),
        "3_trigonal_planar": chem.Atom("Al"),
        "3_trigonal_pyramidal": chem.Atom("P"),
        "4_seesaw": chem.Atom("P"),
        "4_square_planar": chem.Atom("Pt"),
        "4_tetrahedral": chem.Atom("Si"),
        "5_pentagonal_planar": chem.Atom("Pd"),
        "5_square_pyramidal": chem.Atom("Fe"),
        "5_trigonal_bipyramidal": chem.Atom("P"),
        "6_hexagonal_planar": chem.Atom("Fe"),
        "6_octahedral": chem.Atom("Fe"),
        "6_pentagonal_pyramidal": chem.Atom("Fe"),
    }
    return (
        dummy_centers[geometry_name] if geometry_name in dummy_centers.keys() else chem.Atom("Fe")
    )


def get_dummy_atom_list(bond_num):
    # Get the dummy atom for the conformer
    """Return the dummy atom list."""
    dummy_atom = {
        1: [chem.Atom("O"), chem.Atom("Ag")],
        2: [chem.Atom("Al"), chem.Atom("P")],
        3: [chem.Atom("Si")],
        4: [chem.Atom("P")],
        5: [chem.Atom("Fe"), chem.Atom("S")],
        6: [chem.Atom("Fe")],
    }
    if bond_num in dummy_atom:
        return dummy_atom[bond_num]
    else:
        return [chem.Atom("Fe")]


def check_dummy_atom(binding_indices, binding_infos, option):
    """Check the dummy atom."""
    if option == 0:
        return len(binding_indices) > 1  # Only for haptic

    elif option == 1:
        return len(binding_indices) > 1 or len(binding_infos) == 1  # Haptic + single atom bonded

    else:
        return True


def initialize_molecule_properties(metal_complex, option):
    """Initialize the molecule properties."""
    total_atom_num = 1
    total_wo_dummy = 1
    total_chg = 0
    total_mult = metal_complex.multiplicity

    ligands = metal_complex.ligands
    atom_indices_for_each_ligand = metal_complex.get_atom_indices_for_each_ligand()

    for ligand in ligands:
        total_atom_num += len(ligand.molecule.atom_list)
        total_wo_dummy += len(ligand.molecule.atom_list)

    total_atom_list = [None] * total_atom_num

    for i in range(len(ligands)):
        ligand = ligands[i]
        atom_indices = atom_indices_for_each_ligand[i]
        atom_list = ligand.molecule.atom_list
        binding_infos = ligand.binding_infos

        total_chg += ligand.molecule.chg

        for j in range(len(atom_list)):
            total_atom_list[atom_indices[j]] = atom_list[j].copy()

        for j in range(len(binding_infos)):
            binding_indices = binding_infos[j][0]
            if check_dummy_atom(binding_indices, binding_infos, option):
                total_atom_num += 1

    return total_atom_num, total_wo_dummy, total_chg, total_mult, total_atom_list


def get_alternative_molecule(metal_complex, option):
    """Returns the alternative molecule for the conformer, with dummy atoms and dummy center,.

    params: metal_complex: complex_geometry.MetalComplex object
            distance: float
            option: True: set dummy atoms selectively, False: set dummy atoms uniformly

    returns: ace_mol: chem.Molecule object
    """
    metal_index = metal_complex.metal_index
    ligands = metal_complex.ligands
    atom_indices_for_each_ligand = metal_complex.get_atom_indices_for_each_ligand()
    metal_complex.get_adj_matrix()
    metal_complex.get_binding_groups()

    geometry = metal_complex.geometry_type
    geometry_name = geometry.geometry_name

    # Initialize molecule properties
    dummy_center = get_dummy_center(geometry_name)
    dummy_indices = []  # Indices of dummy atom
    # The final atom connecting to the metal atom (because non-haptic ligands)
    metal_binding_infos = []

    total_atom_num, total_wo_dummy, total_chg, total_mult, total_atom_list = (
        initialize_molecule_properties(metal_complex, option)
    )
    total_adj = np.zeros((total_atom_num, total_atom_num))

    # Accumulate each ligand's carried C=C stereo and sp3 chirality, remapped from
    # ligand-local to complex-global atom indices, so both can be enforced on the
    # embed rd_mol.
    total_stereo_bonds = []
    total_chiral_centers = []

    m = 0

    dummy_atom_cn_list = []

    for i in range(len(ligands)):
        ligand = ligands[i]
        atom_indices = atom_indices_for_each_ligand[i]

        binding_infos = ligand.binding_infos
        atom_list = ligand.molecule.atom_list
        mol_adj = ligand.molecule.adj_matrix

        for si, sj, stereo, sra, srb in getattr(ligand.molecule, "stereo_bonds", []):
            total_stereo_bonds.append(
                (
                    atom_indices[si],
                    atom_indices[sj],
                    stereo,
                    atom_indices[sra],
                    atom_indices[srb],
                )
            )

        for center, nbrs, tag in getattr(ligand.molecule, "chiral_centers", []):
            total_chiral_centers.append(
                (atom_indices[center], tuple(atom_indices[k] for k in nbrs), tag)
            )

        for j in range(len(atom_list)):
            for k in range(len(atom_list)):
                total_adj[atom_indices[j]][atom_indices[k]] = mol_adj[j][k]
                total_adj[atom_indices[k]][atom_indices[j]] = mol_adj[k][j]

        for j in range(len(binding_infos)):
            binding_info = binding_infos[j]
            binding_indices, binding_site = binding_info
            is_multidentate = len(binding_infos) > 1
            is_haptic = len(binding_indices) > 1
            if check_dummy_atom(binding_indices, binding_infos, option):
                for k in range(len(binding_indices)):
                    idx = atom_indices[binding_indices[k]]
                    total_adj[idx][total_wo_dummy + m] = 1
                    total_adj[total_wo_dummy + m][idx] = 1
                dummy_atom_cn_list.append(len(binding_indices))
                dummy_indices.append(total_wo_dummy + m)
                metal_binding_infos.append(
                    [total_wo_dummy + m, binding_site, is_multidentate, is_haptic]
                )
                m += 1
            else:
                idx = atom_indices[binding_indices[0]]
                metal_binding_infos.append([idx, binding_site, is_multidentate, is_haptic])

    total_atom_list[metal_index] = get_dummy_center_for_valid(geometry_name)

    # Connect center_binding atoms to the center
    for metal_binding_info in metal_binding_infos:
        binding_index = metal_binding_info[0]
        total_adj[binding_index][metal_index] = 1
        total_adj[metal_index][binding_index] = 1

    len(total_adj)

    # Make alternative molecules ...

    ace_mol_list = []
    total_dummy_atom_list = []
    if len(dummy_atom_cn_list) > 0:
        cn_dict = dict()
        dummy_atom_candidate_list = []
        for cn in dummy_atom_cn_list:
            if cn not in cn_dict:
                cn_dict[cn] = len(cn_dict)
                dummy_atom_candidate_list.append(get_dummy_atom_list(cn))
        dummy_atom_candidate_list = list(itertools.product(*dummy_atom_candidate_list))
        for candidate_list in dummy_atom_candidate_list:
            total_dummy_atom_list.append([candidate_list[cn_dict[cn]] for cn in dummy_atom_cn_list])
    else:
        total_dummy_atom_list.append([])

    for dummy_atom_list in total_dummy_atom_list:
        ace_mol = chem.Molecule()
        ace_mol.chg = total_chg
        ace_mol.multiplicity = total_mult
        ace_mol.atom_list = total_atom_list + [atom.copy() for atom in dummy_atom_list]
        ace_mol.adj_matrix = total_adj

        process.group_molecules(total_adj)

        metal_indices = [metal_index]

        try:
            valid_ace_mol = ace_mol.get_valid_molecule(
                False, method="pulp", MetalCenters=metal_indices
            )
        except Exception:
            valid_ace_mol = None

        if valid_ace_mol is None:
            try:
                valid_ace_mol = ace_mol.get_valid_molecule(
                    False, method="xyz2mol", MetalCenters=metal_indices
                )
            except Exception:
                pass

        if valid_ace_mol is None:
            _telemetry.record("embed.pulp_and_xyz2mol_both_failed", n_atoms=len(ace_mol.atom_list))
            valid_ace_mol = ace_mol
            valid_ace_mol.bo_matrix = ace_mol.adj_matrix.copy()
            valid_ace_mol.chg_list = np.zeros(len(ace_mol.atom_list))

        valid_ace_mol.chg = total_chg
        valid_ace_mol.multiplicity = total_mult

        valid_ace_mol.atom_list[metal_index] = dummy_center

        valid_ace_mol.stereo_bonds = list(total_stereo_bonds)
        valid_ace_mol.chiral_centers = list(total_chiral_centers)

        ace_mol_list.append(valid_ace_mol)

    return ace_mol_list, dummy_indices, metal_binding_infos


def get_repulsive_potential(coordinate_list, d_criteria=0.5, p=6):
    """Return the repulsive potential."""
    distance_matrix = cdist(coordinate_list, coordinate_list)
    potential_matrix = 1 / (distance_matrix - d_criteria) ** p
    return np.sum(potential_matrix)


def kabsch(P, Q):
    """Least-squares rigid transform mapping points ``P`` onto reference ``Q``.

    Returns ``(rot, trans)`` -- a proper rotation matrix (``det == +1``) and a
    translation vector minimizing ``sum_i || rot @ P[i] + trans - Q[i] ||**2`` over
    the correspondence ``P[i] <-> Q[i]``. Scale is fixed to 1: bond lengths placed by
    this transform must stay physical, so the fit may only rotate/translate, never
    stretch.

    The ``d`` sign term below is the Umeyama reflection guard. ``svd`` alone would
    return the improper transform (a reflection, ``det == -1``) whenever that fits a
    mirrored correspondence better -- which would silently flip a chiral ligand's
    handedness and, with it, metal-centered stereochemistry (Delta/Lambda, cis/trans).
    Forcing ``det == +1`` means a mirror image is fitted with a large residual instead
    of being accepted as the same structure.

    ``P`` and ``Q`` are ``(n, 3)`` array-likes with ``n >= 1`` rows in correspondence.
    Both orthogonal factors have ``det == +-1``, so ``d`` is always exactly ``+-1``
    (never 0); a single point (``n == 1``) yields identity rotation + pure translation.
    """
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    p_bar = P.mean(axis=0)
    q_bar = Q.mean(axis=0)
    H = (P - p_bar).T @ (Q - q_bar)
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    rot = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    trans = q_bar - rot @ p_bar
    return rot, trans


def align_double_single_ligand(metal_complex, positions, d_criteria=1.7):
    """Align double single ligand."""
    ligands = metal_complex.ligands
    atom_indices_for_each_ligand = metal_complex.get_atom_indices_for_each_ligand()

    # ratio_criteria = 0.6
    # d_criteria = 0.5

    for i in range(len(ligands)):
        tmp_positions = positions.copy()

        ligand = ligands[i]
        atom_indices = atom_indices_for_each_ligand[i]
        binding_infos = ligand.binding_infos

        if len(binding_infos) > 1 or len(binding_infos[0][0]) > 1:
            continue

        binding_indices, binding_site = binding_infos[0]
        adj_matrix = ligand.get_adj_matrix()
        binding_neighbors = np.where(adj_matrix[binding_indices[0], :] == 1)[0]

        l_positions = tmp_positions[atom_indices, :].copy()

        # First align
        v1 = l_positions[binding_indices[0], :]
        v2 = metal_complex.geometry_type.direction_vector[binding_site - 1]

        v1 = v1 / np.linalg.norm(v1)
        v2 = v2 / np.linalg.norm(v2)

        angle = np.arccos(np.dot(v1, v2))
        axis = np.cross(v1, v2)

        if np.linalg.norm(axis) > 0:
            axis = axis / np.linalg.norm(axis)
            r = R.from_rotvec(angle * axis)
            l_positions = r.apply(l_positions)

        # Shift to zero
        d = l_positions[binding_indices[0], :]
        l_positions = l_positions - l_positions[binding_indices[0], :]

        # Second align
        v1 = np.mean(l_positions[binding_neighbors, :], axis=0)
        v2 = metal_complex.geometry_type.direction_vector[binding_site - 1]

        if np.linalg.norm(v1) < 1e-6:
            continue
        v1 = v1 / np.linalg.norm(v1)
        v2 = v2 / np.linalg.norm(v2)

        angle = np.arccos(np.dot(v1, v2)) - 0.03  # Nearly 180o
        axis = np.cross(v1, v2)

        if np.linalg.norm(axis) > 0:
            axis = axis / np.linalg.norm(axis)
            r = R.from_rotvec(angle * axis)
            l_positions = r.apply(l_positions)

        # Shift back
        l_positions = l_positions + d
        tmp_positions[atom_indices, :] = l_positions

        # Check the validity of the new positions
        # Check the distance between the ligands
        distance_matrix = cdist(tmp_positions, tmp_positions)
        other_indices = np.setdiff1d(np.arange(len(tmp_positions)), atom_indices)
        for j in atom_indices:
            for k in other_indices:
                if distance_matrix[j, k] < d_criteria:
                    continue

        # Check potential
        current_potential = get_repulsive_potential(positions)
        new_potential = get_repulsive_potential(tmp_positions)
        if new_potential > 10 * current_potential:
            continue

        positions = tmp_positions.copy()

    return positions


def _promotion_keeps_valence(rd_mol):
    """True if ``rd_mol`` still passes RDKit's valence check.

    ``SANITIZE_PROPERTIES`` runs the same ``UpdatePropertyCache(strict=True)``
    valence check that a downstream ``SanitizeMol`` / ``MolToSmiles`` performs, on
    a throwaway copy so the caller's mol is left untouched. Used to decide whether
    promoting a bond to DOUBLE would over-valence an atom.
    """
    probe = Chem.Mol(rd_mol)
    try:
        Chem.SanitizeMol(probe, sanitizeOps=Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        return True
    except Exception:
        return False


def _charge_fix_promotion(rd_mol, i, j):
    r"""Restore valence validity after bond ``i``-``j`` was force-set to DOUBLE.

    Promoting a carried bond back to DOUBLE adds one to the bond order of exactly
    its two endpoints, so only ``i`` or ``j`` can become over-valent. PuLP demoted
    the bond to a lower order and put the electrons elsewhere; the encoder's
    charged Lewis form (e.g. a nitrone/N-oxide ``C=[N+]([O-])``, an aldiminate
    ``/N=C(\\[O-])``) instead carries a positive formal charge on the atom that
    gains the bond. Mirror that: bump an endpoint's formal charge by +1 (N 3->4 ->
    N+, O 2->3 -> O+, ...) and accept the first combination that sanitizes. A wrong
    bump (e.g. +1 on a 4-bond carbon -> C+ valence 5) fails the valence probe and
    is rejected, so trying each endpoint then both is self-correcting.

    The generator's internal formal charges never reach the output OIN -- the round
    trip re-encodes the final XYZ from geometry -- so the bump only has to keep the
    mol valence-valid for the embed and the FF cleanup; net charge need not be
    conserved. Returns True with the charge left in place on success; otherwise
    restores the original charges and returns False so the caller can decline.
    """
    a = rd_mol.GetAtomWithIdx(int(i))
    b = rd_mol.GetAtomWithIdx(int(j))
    saved = [(a.GetIdx(), a.GetFormalCharge()), (b.GetIdx(), b.GetFormalCharge())]
    for combo in ([a], [b], [a, b]):
        for atom in combo:
            atom.SetFormalCharge(atom.GetFormalCharge() + 1)
        if _promotion_keeps_valence(rd_mol):
            return True
        for idx, chg in saved:
            rd_mol.GetAtomWithIdx(idx).SetFormalCharge(chg)
    return False


def _apply_double_bond_stereo(rd_mol, stereo_bonds):
    """Set carried C=C (cis/trans) stereo on an embed-ready rd_mol.

    ``stereo_bonds`` is a list of ``(i, j, BondStereo, refA, refB)`` in the mol's
    own atom-index space (``get_rd_mol`` preserves atom order 1:1). Setting the
    stereo plus its two reference atoms makes RDKit distance geometry reproduce
    the requested E/Z; the ace_mol pipeline is otherwise stereo-blind, so the
    double-bond dihedral would embed at random (cis-biased). Skips silently for a
    bond that is absent or no longer double (e.g. a PuLP bond-order re-perception
    relocated it), degrading to the previous behavior rather than raising.
    """
    if not stereo_bonds:
        return
    n = rd_mol.GetNumAtoms()
    changed = False
    for i, j, stereo, ref_a, ref_b in stereo_bonds:
        if max(i, j, ref_a, ref_b) >= n:
            continue
        bond = rd_mol.GetBondBetweenAtoms(int(i), int(j))
        if bond is None:
            continue
        # A carried stereo bond is always a genuine C=C/C=N from the input SMILES.
        # The dummy-metal PuLP re-perception can drop it to SINGLE in the embed mol,
        # which would leave nothing to constrain and let the dihedral (hence the
        # scan target that seeds the FF cleanup) embed at random. Restore the
        # double bond so distance geometry enforces the requested E/Z. Promoting it
        # can over-valence an endpoint, because PuLP paid for the lower bond order
        # with a formal charge the encoder instead put on that atom:
        #   * a nitrone/N-oxide or aldiminate C=N (AHAZOZ, AFECIZ, XIZXAG) -- the N
        #     is +1 in the correct Lewis form, so bumping its charge (via
        #     _charge_fix_promotion) restores validity and the E/Z is enforced;
        #   * FIXYER's C#6 -> valence 5 -- no single +1 bump fixes a doubly-overfull
        #     bond, so the promotion is declined and the bond degrades to random
        #     (the documented fallback) rather than emit an invalid mol.
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            original_type = bond.GetBondType()
            bond.SetBondType(Chem.BondType.DOUBLE)
            if _promotion_keeps_valence(rd_mol):
                changed = True
            elif _charge_fix_promotion(rd_mol, i, j):
                changed = True
            else:
                bond.SetBondType(original_type)
                continue
        try:
            bond.SetStereoAtoms(int(ref_a), int(ref_b))
            bond.SetStereo(stereo)
        except Exception:
            continue
    if changed:
        try:
            rd_mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass


def _permutation_is_odd(source, target):
    """Parity of the permutation taking ``source`` order to ``target`` order."""
    pos = {v: k for k, v in enumerate(source)}
    perm = [pos[v] for v in target]
    inversions = sum(
        1 for a in range(len(perm)) for b in range(a + 1, len(perm)) if perm[a] > perm[b]
    )
    return inversions % 2 == 1


def _apply_atom_chirality(rd_mol, chiral_centers):
    """Set carried sp3 chirality on an embed-ready rd_mol.

    ``chiral_centers`` is a list of ``(center, (n0, n1, n2, n3), ChiralType)``
    where the neighbour tuple records the bond order the tag was read against.
    A chiral tag is only meaningful relative to that ordering, and the mol is
    torn down and rebuilt from the ace_mol (and again by PuLP re-perception), so
    the bond order around the centre can differ. Re-derive the tag by comparing
    the stored order with this mol's actual order: an odd permutation inverts
    the apparent handedness, so flip the tag to keep the same 3D configuration.

    ``EmbedParameters.enforceChirality`` already defaults to True, so the tags set
    here are what the embed honors. Skips silently when the centre's neighbour set
    changed (a re-perceived bond),
    degrading to the previous unconstrained behavior rather than raising.
    """
    if not chiral_centers:
        return
    n = rd_mol.GetNumAtoms()
    flip = {
        Chem.ChiralType.CHI_TETRAHEDRAL_CW: Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
        Chem.ChiralType.CHI_TETRAHEDRAL_CCW: Chem.ChiralType.CHI_TETRAHEDRAL_CW,
    }
    for center, nbrs, tag in chiral_centers:
        if center >= n or max(nbrs) >= n:
            continue
        atom = rd_mol.GetAtomWithIdx(int(center))
        current = [b.GetOtherAtomIdx(int(center)) for b in atom.GetBonds()]
        if sorted(current) != sorted(int(x) for x in nbrs):
            continue
        atom.SetChiralTag(flip[tag] if _permutation_is_odd(nbrs, current) else tag)


def _finalize_positions(
    positions,
    metal_complex,
    new_complex,
    radius_list,
    R,
    metal_index,
    align,
    atom_d_criteria,
    ratio_criteria,
    adj_ratio_criteria,
):
    """Bond-length-correct one embedded conformer and run the three validity criteria.

    Returns ``(good, candidate)`` with exactly one non-None:
      - ``good`` = ``positions[:num_atom]`` when all three criteria pass;
      - ``candidate`` = ``(positions, score)`` when a criterion fails, scored exactly
        as the original inline code did, so the caller can keep it as a best-effort
        fallback.
    Extracted verbatim from ``get_embedding`` so the serial and batched
    (``get_embeddings_batch``) paths share identical geometry acceptance. The serial
    path's byte-identity to pristine is pinned by the golden A/B.
    """
    # Update bond distances that are too long ... (Because of atom replacement)
    q_updates = dict()
    bond_scale = 1.0
    bond_list = np.stack(
        np.where(new_complex.get_adj_matrix() > 0), axis=1
    )  # Not for dummy atom ...
    for bond in bond_list:
        s, e = bond
        if s > e:
            continue
        if s == 0 or e == 0:  # Not for metal ...
            continue
        d = ic.get_distance(positions[: metal_complex.num_atom], s, e)
        # Check ratio ...
        r_sum = radius_list[s] + radius_list[e]
        ratio = d / r_sum
        delta_ratio = ratio - bond_scale
        if delta_ratio > 0.2:  # If too long ...
            q_updates[(s, e)] = -delta_ratio * r_sum
        else:
            q_updates[(s, e)] = 0.0
    ic.update_xyz(positions[: metal_complex.num_atom], q_updates)  # Others are the same ...

    if align is True:
        positions = align_double_single_ligand(new_complex, positions)

    # Check the validity of initial embedding ...
    # Check the distance matrix between different atoms ... (first criterion)
    distance_matrix = cdist(
        positions[: metal_complex.num_atom], positions[: metal_complex.num_atom]
    )
    np.fill_diagonal(distance_matrix, 1e6)
    if np.any(distance_matrix < atom_d_criteria):
        logger.debug("Atoms are too close ...")
        return None, (positions, -100000)

    # Check the collapse between ligands with ratio method (second criterion)
    ratio_matrix = distance_matrix / R
    min_ratio = np.min(ratio_matrix)
    if min_ratio < ratio_criteria:
        logger.debug("Atoms are too close ...")
        return None, (positions, -50000)

    # Finally, check the ratios that are ambiguous, between the ligands ...
    adj_matrix = np.where(distance_matrix / R < adj_ratio_criteria, 1, 0)
    original_adj_matrix = metal_complex.get_adj_matrix()
    # Remove bonds between the metal ...
    adj_matrix[metal_index, :] = 0.0
    adj_matrix[:, metal_index] = 0.0
    original_adj_matrix[metal_index, :] = 0.0
    original_adj_matrix[:, metal_index] = 0.0

    diff = np.sum(np.abs(adj_matrix - original_adj_matrix))

    if diff > 0:
        logger.debug("Undesired bond is detected ...")
        return None, (positions, -diff)

    # Fourth criterion (vdW steric clash) -- gated OFF by default (clash.VDW_ACCEPTANCE_ENABLED).
    # The three checks above forbid only atomic *fusion* (covalent overlap); a conformer can
    # clear them yet carry real inter-fragment vdW clashes -- the release's 53%-vs-5% finding,
    # from a gate blind to steric overlap. When enabled, reject a clashing conformer so the loop
    # keeps searching for a vdW-clean one; score it into (-1, 0) -- above any topology-broken
    # (-diff <= -2) candidate, least-clashing nearest 0 -- so the caller's unchanged best-rejected
    # fallback surfaces the least-clashing candidate when none is clean. Off by default because on
    # the current placement it loosens coordination (see clash.py); A4/A5 evaluate it on the
    # tight Kabsch pool. Disabled -> this returns success exactly as pre-A3 (byte-identical).
    if clash.VDW_ACCEPTANCE_ENABLED:
        atomic_numbers = [atom.get_atomic_number() for atom in new_complex.get_atom_list()]
        clash_vdw, _clash_severe, worst_overlap = clash.vdw_clash_count(
            positions[: metal_complex.num_atom], atomic_numbers
        )
        if clash_vdw > 0:
            logger.debug(f"vdW steric clash detected ... ({clash_vdw} pair(s))")
            penalty = clash_vdw + (1.0 - worst_overlap)
            score = -1.0 + 1.0 / (1.0 + penalty)
            return None, (positions, score)

    logger.debug("Embedding success!\n")
    return positions[: metal_complex.num_atom], None


# ---------------------------------------------------------------------------
# Rigid-placement embed (opt-in ``option=3``) -- v0.4.3 A4.
#
# Instead of the single dummy-metal distance-geometry embed (which disguises the
# metal and stretches ligand interiors to satisfy one global CoordMap), build each
# ligand independently with a clean ETKDG conformer, then RIGIDLY place it onto the
# ideal coordination vectors with the reflection-guarded ``kabsch`` fit above. The
# metal is fixed at the origin; the assembled guess is then relaxed by the SAME
# ``_finalize_positions`` criteria + downstream FF the DG path uses. This is gated
# off by default (``ff_params["use_kabsch"]``) so the default pool stays byte-
# identical until A5 decides whether to promote it.
# ---------------------------------------------------------------------------


def _unit(v):
    """Return ``v`` normalized; a zero vector is returned unchanged."""
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


def _rotation_align(a, b):
    """Proper rotation mapping unit vector ``a`` onto unit vector ``b``.

    Handles the parallel (identity) and antiparallel (180 deg about any
    perpendicular axis) degeneracies that a bare cross-product axis-angle misses.
    """
    a = _unit(a)
    b = _unit(b)
    c = float(np.dot(a, b))
    if c > 1 - 1e-9:
        return np.eye(3)
    if c < -1 + 1e-9:
        # 180 deg: rotate about any axis perpendicular to ``a``.
        perp = np.cross(a, [1.0, 0.0, 0.0])
        if np.linalg.norm(perp) < 1e-6:
            perp = np.cross(a, [0.0, 1.0, 0.0])
        return R.from_rotvec(np.pi * _unit(perp)).as_matrix()
    axis = _unit(np.cross(a, b))
    angle = np.arccos(np.clip(c, -1.0, 1.0))
    return R.from_rotvec(angle * axis).as_matrix()


def _embed_free_ligand(ligand, seed):
    """Return clean ligand-local coordinates ``(n, 3)`` from a lone ETKDG embed.

    Builds the ligand fragment ALONE (no dummy metal, no CoordMap), so its interior
    is undistorted -- the direct fix for the dummy-metal bond-length patch the DG
    path needs. Atom order matches ``ligand.molecule.atom_list`` (``get_rd_mol``
    preserves it), so local atom ``j`` maps to global slot ``atom_indices[j]``.
    Single-atom donors (a lone halide) skip the embed. Returns ``None`` if the
    fragment will not embed, so the caller falls through to the next attempt.
    """
    atom_list = ligand.molecule.atom_list
    n = len(atom_list)
    if n == 1:
        return np.zeros((1, 3))
    try:
        rd_mol = ligand.molecule.get_rd_mol()
    except Exception:
        return None
    rd_mol = Chem.Mol(rd_mol)
    params = rdDistGeom.srETKDGv3()
    params.randomSeed = int(seed)
    cid = AllChem.EmbedMolecule(rd_mol, params)
    if cid < 0:
        params.useRandomCoords = True
        cid = AllChem.EmbedMolecule(rd_mol, params)
    if cid < 0:
        return None
    try:
        return rd_mol.GetConformer(cid).GetPositions()
    except Exception:
        return None


def _place_monodentate(P, donor_local, neighbor_locals, v_slot, bond_len):
    """Place a single-donor ligand: donor exactly on ``v_slot*bond_len``, body out.

    Orients the donor->neighbor-centroid axis along ``+v_slot`` (substituents point
    AWAY from the metal at the origin) via the shared ``kabsch`` fit, then pins the
    donor exactly on the ideal vector. The remaining spin about the M-donor axis is
    left free here and resolved later by the packing pass.
    """
    target_donor = v_slot * bond_len
    if len(neighbor_locals) == 0:
        return P - P[donor_local] + target_donor
    nbr_centroid = P[neighbor_locals].mean(axis=0)
    src = np.vstack([P[donor_local], nbr_centroid])
    axis_len = np.linalg.norm(nbr_centroid - P[donor_local])
    tgt = np.vstack([target_donor, target_donor + v_slot * axis_len])
    rot, _ = kabsch(src, tgt)
    placed = (rot @ P.T).T
    return placed - placed[donor_local] + target_donor


def _place_haptic(P, face_locals, v_slot, metal_r, atom_r, scale, winding=None):
    """Place a haptic (eta) face: ring plane normal along ``v_slot``, centroid on it.

    Rotation only (the ring's own geometry is preserved, never distorted). The
    centroid height is chosen so the metal-to-face-atom distance is ~physical for
    the requested ``scale``. The SVD ring normal has an undetermined sign, so the
    initial placement presents an arbitrary face to the metal.

    ``winding`` (SL2 oin-direct-winding, ``None`` on the default kabsch path ->
    byte-identical) is ``(face_order, star_rank, char)`` from ``Ligand.winding``:
    the requested OIN winding for a *load-bearing* ring (indenyl, ansa, substituted
    Cp whose rac/meso must be reproduced). When given, measure the placed ring's
    winding with the encoder's shared ``signed_circulation`` -- SAME star atom and
    SAME ascending-fragment-order the encoder uses, so the sign convention cannot
    diverge -- and if it disagrees, turn the ring over with ONE 180 deg PROPER
    rotation about an in-plane axis through the ring centroid. That flip reverses
    the winding sign (``cross(v_star, v_next) || ring normal``, and the rotation
    sends normal -> -normal), preserves every metal-ring distance (the centroid
    lies on the axis; in-plane radii keep their length), and -- being a proper
    rotation, never a reflection -- leaves substituent / rac-vs-meso chirality
    intact. An orientation-free ring (Cp/arene) also passes through here harmlessly:
    the flip realizes a graph automorphism, giving an equivalent conformer.
    """
    face = P[face_locals]
    centroid = face.mean(axis=0)
    centered = face - centroid
    # Ring normal = smallest-variance singular direction of the face points.
    _u, _s, vt = np.linalg.svd(centered)
    normal = vt[2]
    r_ring = float(np.mean(np.linalg.norm(centered, axis=1)))
    target_atom_dist = (metal_r + atom_r) * scale
    h2 = target_atom_dist**2 - r_ring**2
    h = np.sqrt(h2) if h2 > 1e-6 else 0.1
    rot = _rotation_align(normal, v_slot)
    placed = (rot @ (P - centroid).T).T + v_slot * h

    if winding is not None:
        face_order, star_rank, target_char = winding
        ring_coords = placed[face_order]
        # Metal is at the origin, so the ring centroid IS the outward
        # metal->centroid axis -- exactly the encoder's convention
        # (oin_aligner._determine_winding uses the actual ring centroid).
        face_c = ring_coords.mean(axis=0)
        if signed_circulation(ring_coords, star_rank, face_c) != target_char:
            # In-plane flip axis: any unit vector perpendicular to the ring
            # normal (== v_slot after alignment). Deterministic choice.
            u = np.cross(v_slot, [1.0, 0.0, 0.0])
            if np.linalg.norm(u) < 1e-6:
                u = np.cross(v_slot, [0.0, 1.0, 0.0])
            flip = R.from_rotvec(np.pi * _unit(u)).as_matrix()
            placed = (flip @ (placed - face_c).T).T + face_c

    return placed


def _place_chelate(P, donor_locals, slot_vectors, bond_lens):
    """Place a polydentate chelate: least-squares fit its donors onto their slots.

    ``kabsch`` (reflection-guarded, so Delta/Lambda handedness is preserved) fits the
    ``m`` donor atoms onto their ``m`` ideal slot targets. The donor-donor spacing of
    a free ligand rarely equals the polyhedron's inter-vertex spacing, so a nonzero
    **bite residual** remains by construction -- the downstream constrained FF closes
    it. For a bidentate (2 donors, which underdetermine the spin about the donor-donor
    axis) the ligand centroid is added as a third correspondence pointing AWAY from
    the metal, so the body is never folded back over the metal.
    """
    donor_src = P[donor_locals]
    donor_tgt = np.array(
        [slot_vectors[k] * bond_lens[k] for k in range(len(donor_locals))], dtype=float
    )
    src, tgt = donor_src, donor_tgt
    if len(donor_locals) == 2:
        lig_centroid = P.mean(axis=0)
        donor_mid = donor_tgt.mean(axis=0)
        outward = _unit(donor_mid)
        body_len = np.linalg.norm(lig_centroid - donor_src.mean(axis=0))
        src = np.vstack([donor_src, lig_centroid])
        tgt = np.vstack([donor_tgt, donor_mid + outward * body_len])
    rot, trans = kabsch(src, tgt)
    return (rot @ P.T).T + trans


def _bite_residual(placed, donor_locals, slot_vectors, bond_lens):
    """RMS distance of placed donors from their ideal slot targets (chelate strain)."""
    tgt = np.array([slot_vectors[k] * bond_lens[k] for k in range(len(donor_locals))], dtype=float)
    got = placed[list(donor_locals)]
    return float(np.sqrt(np.mean(np.sum((got - tgt) ** 2, axis=1))))


def _rotate_about_axis(coords, pivot, axis, angle):
    """Rotate ``coords`` about the line through ``pivot`` along unit ``axis`` by ``angle`` (rad).

    ``axis`` must already be a unit vector (callers pass ``_unit(axis)``) so the
    rotation matrix stays bit-identical to the pre-factoring inline math.
    """
    rot = R.from_rotvec(angle * axis).as_matrix()
    return (rot @ (coords - pivot).T).T + pivot


def _pack_ligands(positions, ligand_axes, metal_index, n_steps=12, n_passes=3, score_indices=None):
    """Phase-3 packing: torsionally de-clash monodentate bodies about their M-donor axes.

    Independent placement can interpenetrate bulky ligands on adjacent vertices.
    For each ligand with a well-defined M-donor axis (monodentate), sweep a grid of
    rotations of its atoms about that axis and keep the angle that lowers the total
    repulsive potential. This is STRICTLY guarded -- a rotation is accepted only when
    it reduces ``get_repulsive_potential`` -- so packing can never ship a worse
    (more clashing) structure than plain placement.

    ``ligand_axes`` maps each ligand's global atom indices to its ``(pivot, axis)``
    (or ``None`` to skip that ligand). A monodentate rotates about its M-donor axis;
    a bidentate rotates about its donor-donor axis -- both keep every donor on its
    slot. Returns the (possibly improved) positions.

    ``n_passes`` sequential sweeps let a ligand react to neighbours moved earlier in
    the same pass (mutually-clashing tris-chelates need more than one look).

    ``score_indices`` (default None) restricts the repulsive-potential score to a
    subset of atoms. None scores over the whole complex -- byte-identical to before;
    a list scores over only those atoms, used by the greedy path to settle the
    already-placed *prefix* while the not-yet-placed atoms still sit at the origin.
    """
    positions = positions.copy()
    angles = [2 * np.pi * k / n_steps for k in range(1, n_steps)]

    def _score(pos):
        return get_repulsive_potential(pos if score_indices is None else pos[score_indices])

    for _pass in range(n_passes):
        improved_any = False
        for atom_indices, axis_info in ligand_axes:
            if axis_info is None:
                continue
            pivot, axis = axis_info
            axis = _unit(axis)
            if np.linalg.norm(axis) < 1e-9:
                continue
            best_potential = _score(positions)
            best_positions = None
            for angle in angles:
                trial = positions.copy()
                block = positions[atom_indices]
                trial[atom_indices] = _rotate_about_axis(block, pivot, axis, angle)
                pot = _score(trial)
                if pot < best_potential:
                    best_potential = pot
                    best_positions = trial
            if best_positions is not None:
                positions = best_positions
                improved_any = True
        if not improved_any:
            break
    return positions


def _gyration_radius(coords):
    """Radius of gyration of a point set -- a cheap steric-bulk proxy."""
    coords = np.asarray(coords, dtype=float)
    if len(coords) < 2:
        return 0.0
    centroid = coords.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((coords - centroid) ** 2, axis=1))))


def _ligand_difficulty_key(ligand, free_coords):
    """Placement-order key: most-constrained ligand FIRST (ascending sort).

    denticity desc -> haptic-before-monodentate -> steric bulk (radius of gyration)
    desc -> atom-count desc. A chelate/pincer anchors the sphere (it constrains the
    most slots); a haptic ring and bulky monodentates fit next; a lone halide
    (Rg 0, one atom) drops to the end, into whatever space remains. (SL3 ordering,
    user-confirmed.)
    """
    binding_infos = ligand.binding_infos
    denticity = len(binding_infos)
    is_haptic = 1 if (denticity == 1 and len(binding_infos[0][0]) > 1) else 0
    n_atoms = len(ligand.molecule.atom_list)
    rg = _gyration_radius(free_coords)
    return (-denticity, -is_haptic, -rg, -n_atoms)


def _place_one_ligand(ligand, P, direction_vector, metal_r, scale, greedy=False):
    """Rigidly place one free-embedded ligand ``P`` onto its ideal slot(s).

    Shared by the independent (``greedy=False``, byte-identical) and greedy paths.
    Returns ``(placed, axis_info)`` where ``axis_info`` is the ``(pivot, axis)`` for
    the guarded torsional de-clash, or ``None`` when no donor-preserving axis exists.
    A haptic face gets a de-clash axis (its slot-normal -- an in-plane spin that
    preserves winding sign AND every metal-ring distance) ONLY in greedy mode; on the
    default path it keeps ``None`` so the existing option=3 output is byte-identical.
    """
    atom_list = ligand.molecule.atom_list
    binding_infos = ligand.binding_infos
    if len(binding_infos) == 1 and len(binding_infos[0][0]) == 1:
        # Monodentate single donor.
        donor_local = binding_infos[0][0][0]
        slot = binding_infos[0][1]
        v_slot = _unit(direction_vector[slot - 1])
        bond_len = (metal_r + atom_list[donor_local].get_radius()) * scale
        adj_matrix = ligand.get_adj_matrix()
        neighbor_locals = np.where(adj_matrix[donor_local, :] == 1)[0]
        placed = _place_monodentate(P, donor_local, neighbor_locals, v_slot, bond_len)
        axis_info = (v_slot * bond_len, v_slot)
    elif len(binding_infos) == 1 and len(binding_infos[0][0]) > 1:
        # Haptic (eta) face.
        face_locals = binding_infos[0][0]
        slot = binding_infos[0][1]
        v_slot = _unit(direction_vector[slot - 1])
        atom_r = atom_list[face_locals[0]].get_radius()
        placed = _place_haptic(
            P, face_locals, v_slot, metal_r, atom_r, scale, winding=ligand.winding
        )
        axis_info = None
        if greedy:
            centroid = placed[list(face_locals)].mean(axis=0)
            axis_info = (centroid, v_slot)
    else:
        # Polydentate chelate: one donor per binding_info.
        donor_locals = [bi[0][0] for bi in binding_infos]
        slots = [bi[1] for bi in binding_infos]
        slot_vectors = [_unit(direction_vector[s - 1]) for s in slots]
        bond_lens = [(metal_r + atom_list[dl].get_radius()) * scale for dl in donor_locals]
        placed = _place_chelate(P, donor_locals, slot_vectors, bond_lens)
        logger.debug(
            "kabsch chelate bite residual: %.3f",
            _bite_residual(placed, donor_locals, slot_vectors, bond_lens),
        )
        axis_info = None
        if len(donor_locals) == 2:
            d1 = placed[donor_locals[0]]
            d2 = placed[donor_locals[1]]
            axis_info = (d1, d2 - d1)
    return placed, axis_info


def _interligand_clash_cost(ligand_coords, other_coords):
    """Cross-set repulsive potential between a ligand and the already-placed atoms.

    Same soft ``1/(d-0.5)**6`` form as the intra-set ``get_repulsive_potential``.
    Intra-ligand pairs are rigid under a rotation about the donor axis, so scoring
    only the cross term is equivalent to -- and cheaper than -- scoring the merged set.
    """
    d = cdist(np.asarray(ligand_coords, dtype=float), np.asarray(other_coords, dtype=float))
    return float(np.sum(1.0 / (d - 0.5) ** 6))


def _place_ligand_collision_aware(placed, axis_info, positions, committed_indices, n_steps=12):
    """Spin ``placed`` about its donor axis to the pose that least clashes with others.

    Donors lie ON the axis, so every trial keeps them exactly on-slot -- only the free
    rotational DOF is searched (the on-slot fidelity is what fixes fac/mer). Returns
    the chosen pose (the base pose if there is no axis or nothing to avoid).
    """
    if axis_info is None or len(committed_indices) == 0:
        return placed
    pivot, axis = axis_info
    axis = _unit(axis)
    if np.linalg.norm(axis) < 1e-9:
        return placed
    other = positions[committed_indices]
    best_pose = placed
    best_cost = _interligand_clash_cost(placed, other)
    for k in range(1, n_steps):
        angle = 2 * np.pi * k / n_steps
        trial = _rotate_about_axis(placed, pivot, axis, angle)
        cost = _interligand_clash_cost(trial, other)
        if cost < best_cost:
            best_cost = cost
            best_pose = trial
    return best_pose


def _kabsch_embedding(
    metal_complex,
    new_complex,
    scale,
    align,
    seed,
    atom_d_criteria,
    ratio_criteria,
    adj_ratio_criteria,
    greedy=False,
):
    """Independent-build + rigid-placement embed for ``option=3`` (see module note).

    Returns accepted ``positions[:num_atom]`` (metal-first ordering, the downstream
    contract) or a best-effort fallback via the shared ``_finalize_positions`` -- the
    exact same return contract as ``get_embedding`` -- or ``None`` if a ligand will
    not build.
    """
    total_atom_list = new_complex.get_atom_list()
    radius_list = [atom.get_radius() for atom in total_atom_list]
    n = len(radius_list)
    R_sum = np.repeat(np.array(radius_list), n).reshape((n, n))
    R_sum = R_sum + R_sum.T

    direction_vector = np.asarray(new_complex.geometry_type.direction_vector, dtype=float)
    metal_index = new_complex.metal_index
    metal_r = new_complex.center_atom.get_radius()
    num_atom = new_complex.num_atom
    atom_indices_for_each_ligand = new_complex.get_atom_indices_for_each_ligand()

    positions = np.zeros((num_atom, 3))
    ligand_axes = []

    if seed is None:
        seed = random.randint(0, 1000000)

    ligands = list(new_complex.ligands)

    if not greedy:
        # Independent placement in parse order (byte-identical option=3 / A4 path).
        for i, ligand in enumerate(ligands):
            P = _embed_free_ligand(ligand, seed + i)
            if P is None:
                return None
            global_indices = atom_indices_for_each_ligand[i]
            placed, axis_info = _place_one_ligand(
                ligand, P, direction_vector, metal_r, scale, greedy=False
            )
            for j, g in enumerate(global_indices):
                positions[g] = placed[j]
            ligand_axes.append((list(global_indices), axis_info))
        positions[metal_index] = np.array([0.0, 0.0, 0.0])
        # Phase 3 -- guarded torsional packing (only kept if it lowers clash).
        positions = _pack_ligands(positions, ligand_axes, metal_index)
    else:
        # Greedy: difficulty-ordered, collision-aware sequential placement (SL3).
        # Pre-embed every ligand FIRST, seeded by its ORIGINAL index so the free
        # conformer -- and thus byte-level reproducibility -- is independent of the
        # placement order we are about to choose.
        free_coords = []
        for i, ligand in enumerate(ligands):
            P = _embed_free_ligand(ligand, seed + i)
            if P is None:
                return None
            free_coords.append(P)
        order = sorted(
            range(len(ligands)),
            key=lambda i: _ligand_difficulty_key(ligands[i], free_coords[i]),
        )
        # Metal anchors the origin first; each ligand is then placed into the space
        # the higher-priority ligands already occupy, and the committed prefix is
        # re-settled so earlier ligands can yield to the newcomer. Atoms are always
        # written to their canonical global indices, so metal-first / slot ordering
        # is preserved -- only the *iteration* order changes.
        positions[metal_index] = np.array([0.0, 0.0, 0.0])
        committed_indices = [metal_index]
        committed_axes = []
        for i in order:
            ligand = ligands[i]
            global_indices = atom_indices_for_each_ligand[i]
            placed, axis_info = _place_one_ligand(
                ligand, free_coords[i], direction_vector, metal_r, scale, greedy=True
            )
            placed = _place_ligand_collision_aware(
                placed, axis_info, positions, committed_indices
            )
            for j, g in enumerate(global_indices):
                positions[g] = placed[j]
            committed_indices.extend(global_indices)
            committed_axes.append((list(global_indices), axis_info))
            ligand_axes.append((list(global_indices), axis_info))
            # Guarded between-placement settle of the committed prefix only (the
            # not-yet-placed atoms still sit at the origin, so score on the prefix).
            positions = _pack_ligands(
                positions,
                committed_axes,
                metal_index,
                n_passes=1,
                score_indices=committed_indices,
            )
        # Final guarded torsional packing over the whole assembled complex.
        positions = _pack_ligands(positions, ligand_axes, metal_index)

    good, candidate = _finalize_positions(
        positions,
        metal_complex,
        new_complex,
        radius_list,
        R_sum,
        metal_index,
        False,  # already rigidly placed; do NOT re-run the monodentate axis-orient
        atom_d_criteria,
        ratio_criteria,
        adj_ratio_criteria,
    )
    if good is not None:
        return good
    return candidate[0][:num_atom]


def get_embedding(
    metal_complex,
    scale=1.0,
    option=0,
    align=False,
    use_random=True,
    seed=None,
    alt_cache=None,
    greedy=False,
):
    """Return the embedding.

    ``seed`` fixes the distance-geometry random seed. Pass an int for a
    reproducible embed; pass None (with ``use_random``) to draw one, which
    makes every sp3 stereocentre's handedness a fresh coin flip per call.

    ``alt_cache`` is an optional per-generation dict (see ``_alt_mol_cached``);
    None means recompute the alternative molecule every call.
    """
    atom_d_criteria = 0.5
    ratio_criteria = 0.65
    adj_ratio_criteria = 1.4

    # Make a conformer based on the representation
    new_complex = metal_complex.copy()
    total_atom_list = new_complex.get_atom_list()

    # Opt-in rigid-placement embed (A4). option=3 bypasses the whole dummy-metal
    # distance-geometry mechanism below -- it builds each ligand independently and
    # places it with kabsch -- so it is dispatched here, before the alternative-mol
    # machinery. Never reached unless the caller opts in via the options list.
    if option == 3:
        return _kabsch_embedding(
            metal_complex,
            new_complex,
            scale,
            align,
            seed,
            atom_d_criteria,
            ratio_criteria,
            adj_ratio_criteria,
            greedy=greedy,
        )

    # new_complex gives the geometry ... (metal complex remains no change)

    alternative_ace_mol_list, dummy_indices, metal_binding_infos = _alt_mol_cached(
        new_complex, option, alt_cache
    )

    metal_index = new_complex.metal_index
    metal_r = new_complex.center_atom.get_radius()

    radius_list = [atom.get_radius() for atom in total_atom_list]  # Before alteration of atoms ...
    n = len(radius_list)
    R = np.repeat(np.array(radius_list), n).reshape((n, n))
    R = R + R.T

    # Make cmap/params ...
    direction_vector = new_complex.geometry_type.direction_vector
    params = rdDistGeom.EmbedParameters()
    # useRandomCoords=True stays the PRIMARY start: it is byte-identical to pristine on
    # every embed that already succeeds, and -- crucially -- flexible ligands (e.g. a
    # Pt-SPL amino-alcohol arm) need the random exploration to fold the donor onto the
    # metal; seeding those from the metric matrix instead lands the donor off-metal and
    # silently drops its binding. The metric-matrix start is applied ONLY as a retry on
    # a -1 return (see the `useRandomCoords=False` block after the primary embed below),
    # so it rescues the OCT loose-scale failures without disturbing anything that works.
    params.useRandomCoords = True
    params.maxIterations = 100
    params.useBasicKnowledge = False
    params.ignoreSmoothingFailures = True

    cmap = dict()

    if seed is not None:
        params.randomSeed = int(seed)
    elif use_random is True:
        params.randomSeed = random.randint(0, 1000000)

    candidate_list = []
    scales_for_haptic = [0.4, 0.5, 0.6, 0.7]

    haptic_exist = False

    for alternative_ace_mol in alternative_ace_mol_list:
        alternative_ace_mol_list.index(alternative_ace_mol)
        rd_mol = alternative_ace_mol.get_rd_mol()
        # Carried stereo persists across the repeated EmbedMolecule calls below
        # (both are atom/bond properties, not cleared by embedding), so set once here.
        stereo_bonds = getattr(alternative_ace_mol, "stereo_bonds", [])
        chiral_centers = getattr(alternative_ace_mol, "chiral_centers", [])
        # No enforceChirality here: rdDistGeom.EmbedParameters() already defaults it
        # to True. It had nothing to act on before, because the ace_mol carried no
        # chiral tags -- _apply_atom_chirality is what gives it something to enforce.
        _apply_double_bond_stereo(rd_mol, stereo_bonds)
        _apply_atom_chirality(rd_mol, chiral_centers)
        logger.debug("%s %s", "Trying ", Chem.MolToSmiles(rd_mol))

        positions = None

        # Try many scale values for haptic (haptic embedding seems to work differently)
        for haptic_scale in scales_for_haptic:
            failed = False
            embed_raised = False  # True only if the primary EmbedMolecule *raised*
            # Make cmap for embedding
            for metal_binding_info in metal_binding_infos:
                binding_index, binding_site, is_multidentate, is_haptic = (
                    metal_binding_info  # list, int, bool, bool
                )
                atom_r = alternative_ace_mol_list[0].atom_list[binding_index].get_radius()
                if is_haptic:
                    distance = (metal_r + atom_r) * haptic_scale
                    haptic_exist = True
                else:
                    distance = (metal_r + atom_r) * scale
                x, y, z = direction_vector[binding_site - 1] * distance
                cmap[binding_index] = Point3D(x, y, z)
            cmap[metal_index] = Point3D(0.0, 0.0, 0.0)
            params.SetCoordMap(cmap)

            # EmbedMolecule returns -1 on failure; it does NOT raise. Capture the
            # return code and treat -1 as failed directly, instead of waiting for
            # GetConformer() to throw. This also lets us skip the rebuild retry
            # below on a plain -1: that retry re-runs the identical seeded embed on
            # the identical PuLP mol (same params.randomSeed, unchanged ace_mol), so
            # it reproduces the same -1 deterministically -- measured 8/8 dead on
            # fac-Ir(ppy)3, and via the `continue` it also suppressed the
            # `if not haptic_exist: break` below, forcing 3 more redundant primary
            # embeds per failing non-haptic combo. The rebuild is only meaningful
            # when the embed RAISED, because the except-branch then swaps the metal
            # to a real element and re-solves, giving the retry a genuinely
            # different molecule.
            try:
                rc = AllChem.EmbedMolecule(rd_mol, params)
            except Exception:
                embed_raised = True
                rc = -1
                try:
                    geometry_name = new_complex.geometry_type.geometry_name
                    temp_metal_num = get_transition_metal_center(geometry_name)
                    alternative_ace_mol.atom_list[metal_index].set_atomic_number(temp_metal_num)
                    rc = AllChem.EmbedMolecule(rd_mol, params)
                except Exception:
                    logger.debug("Embedding failed ...")
                    failed = True
            if rc == -1 and not embed_raised:
                # P9: a plain -1 (no exception) means the random-coords start left the
                # distance-geometry first minimization infeasible. On the OCT
                # tris-chelates the loose-scale (1.1/1.2) cmap over-stretches the
                # fused-aromatic chelate bite, so every one of the 100 iterations fails
                # FIRST_MINIMIZATION (params.GetFailureCounts()) and the embed returns
                # -1 having built nothing. Retry ONCE from the metric-matrix eigen-
                # decomposition (useRandomCoords=False, RDKit's default start), which
                # begins near-feasible and flips these to a valid conformer. This is a
                # MATERIALLY different embed (different initial-coord method) -- not the
                # dead same-params re-run P3 removed, which reproduced the identical -1 --
                # and it stays deterministic (same params.randomSeed). It fires only on a
                # -1, so every embed that already succeeds keeps its exact pristine
                # geometry, including flexible donors whose arm the random start folds
                # onto the metal (a metric-matrix start would drop that binding).
                params.useRandomCoords = False
                try:
                    rc = AllChem.EmbedMolecule(rd_mol, params)
                except Exception:
                    rc = -1
                params.useRandomCoords = True  # restore for the next combo's primary
            if rc == -1:
                failed = True
            try:
                conformer = rd_mol.GetConformer()
            except Exception:
                if not failed:
                    logger.debug("Conformer not obtained ...")
                    failed = True
            try:
                positions = conformer.GetPositions()
            except Exception:
                if not failed:
                    logger.debug("Position not obtained ...")
                    failed = True

            if embed_raised and (failed or positions is None):
                # Try different molecule ...  (raised-embed path only -- see above)
                alternative_ace_mol = alternative_ace_mol.get_valid_molecule(
                    False, method="pulp", MetalCenters=[metal_index]
                )
                rd_mol = alternative_ace_mol.get_rd_mol()
                # get_valid_molecule preserves atom order, so the captured indices
                # still apply; re-set the stereo on the rebuilt rd_mol.
                _apply_double_bond_stereo(rd_mol, stereo_bonds)
                _apply_atom_chirality(rd_mol, chiral_centers)
                logger.debug("%s %s", "Trying ", Chem.MolToSmiles(rd_mol))

                try:
                    rc = AllChem.EmbedMolecule(rd_mol, params)
                except Exception:
                    logger.debug("Embedding failed ...")
                    continue
                if rc == -1:
                    continue
                try:
                    conformer = rd_mol.GetConformer()
                except Exception:
                    logger.debug("Conformer not obtained ...")
                    continue
                try:
                    positions = conformer.GetPositions()
                except Exception:
                    logger.debug("Position not obtained ...")
                    continue

            if not failed:
                break

            # No need to perform multiple embedding ...
            if not haptic_exist:
                break

        # If position does not exist, move on to the next alternative ace mol
        if positions is None:
            continue

        good, candidate = _finalize_positions(
            positions,
            metal_complex,
            new_complex,
            radius_list,
            R,
            metal_index,
            align,
            atom_d_criteria,
            ratio_criteria,
            adj_ratio_criteria,
        )
        if good is not None:
            return good
        candidate_list.append(candidate)
        continue

    logger.debug("Embedding failed ...")

    if len(candidate_list) == 0:
        return None
    else:
        # Get the best position among the position list ... (Check by value ...)
        maximum_value = -100000000
        final_positions = candidate_list[0][0]
        for candidate in candidate_list:
            value = candidate[1]
            if value > maximum_value:
                maximum_value = value
                final_positions = candidate[0]
        logger.debug(f"Returning the best position ... maximum value {maximum_value}")
        _telemetry.record(
            "embed.best_rejected_returned",
            score=maximum_value,
            n_candidates=len(candidate_list),
        )
        return final_positions[: metal_complex.num_atom]


def get_embeddings_batch(
    metal_complex,
    scale=1.0,
    option=0,
    num_confs=4,
    num_threads=0,
    align=False,
    seed=None,
    alt_cache=None,
):
    """Batched, C++-parallel counterpart to :func:`get_embedding` for NON-haptic complexes.

    Runs ``num_confs`` distance-geometry embeds in a single ``EmbedMultipleConfs``
    call, which releases the GIL and parallelizes across ``num_threads`` cores
    (``0`` = all). Every conformer is bond-corrected and criteria-checked by the
    same :func:`_finalize_positions` the serial path uses, so geometry acceptance
    is identical; only the *sampling* differs. Returns a list of accepted positions
    (each ``positions[:num_atom]``), or a single best-effort candidate if no
    conformer passed all criteria (mirrors ``get_embedding``'s fallback so the batch
    is non-regressive).

    Returns ``None`` to signal "not batchable -- fall back to serial
    ``get_embedding``": a haptic complex, whose cmap must sweep ``scales_for_haptic``
    the way the serial path does. The caller keeps that path serial.

    NOT byte-identical to the serial path: ``EmbedMultipleConfs`` derives its own
    per-conformer seed sequence from ``seed``, so the conformers differ from
    ``num_confs`` serial ``EmbedMolecule`` calls. This is why the parallel path is
    gated behind ``num_threads != 1`` (default 1 = serial, byte-identical) and
    validated by the accuracy gate rather than byte-identity.
    """
    atom_d_criteria = 0.5
    ratio_criteria = 0.65
    adj_ratio_criteria = 1.4

    new_complex = metal_complex.copy()
    total_atom_list = new_complex.get_atom_list()

    alternative_ace_mol_list, dummy_indices, metal_binding_infos = _alt_mol_cached(
        new_complex, option, alt_cache
    )

    # Haptic binding needs the serial scales_for_haptic sweep; signal fallback.
    if any(metal_binding_info[3] for metal_binding_info in metal_binding_infos):
        return None

    metal_index = new_complex.metal_index
    metal_r = new_complex.center_atom.get_radius()

    radius_list = [atom.get_radius() for atom in total_atom_list]
    n = len(radius_list)
    R = np.repeat(np.array(radius_list), n).reshape((n, n))
    R = R + R.T

    direction_vector = new_complex.geometry_type.direction_vector
    params = rdDistGeom.EmbedParameters()
    params.useRandomCoords = True
    params.maxIterations = 100
    params.useBasicKnowledge = False
    params.ignoreSmoothingFailures = True
    params.numThreads = int(num_threads)
    if seed is not None:
        params.randomSeed = int(seed)
    else:
        params.randomSeed = random.randint(0, 1000000)

    accepted = []
    candidates = []
    for alternative_ace_mol in alternative_ace_mol_list:
        rd_mol = alternative_ace_mol.get_rd_mol()
        stereo_bonds = getattr(alternative_ace_mol, "stereo_bonds", [])
        chiral_centers = getattr(alternative_ace_mol, "chiral_centers", [])
        _apply_double_bond_stereo(rd_mol, stereo_bonds)
        _apply_atom_chirality(rd_mol, chiral_centers)

        cmap = dict()
        for metal_binding_info in metal_binding_infos:
            binding_index, binding_site, is_multidentate, is_haptic = metal_binding_info
            atom_r = alternative_ace_mol_list[0].atom_list[binding_index].get_radius()
            distance = (metal_r + atom_r) * scale
            x, y, z = direction_vector[binding_site - 1] * distance
            cmap[binding_index] = Point3D(x, y, z)
        cmap[metal_index] = Point3D(0.0, 0.0, 0.0)
        params.SetCoordMap(cmap)

        try:
            conformer_ids = list(AllChem.EmbedMultipleConfs(rd_mol, int(num_confs), params))
        except Exception:
            conformer_ids = []

        for conformer_id in conformer_ids:
            try:
                positions = rd_mol.GetConformer(conformer_id).GetPositions()
            except Exception:
                continue
            good, candidate = _finalize_positions(
                positions,
                metal_complex,
                new_complex,
                radius_list,
                R,
                metal_index,
                align,
                atom_d_criteria,
                ratio_criteria,
                adj_ratio_criteria,
            )
            if good is not None:
                accepted.append(good)
            elif candidate is not None:
                candidates.append(candidate)

        # Prefer the first alternative ace mol that yields any valid conformer,
        # matching get_embedding's "return on first success" bias.
        if accepted:
            return accepted

    if candidates:
        # No fully-valid conformer from any alt mol: hand back the single best-effort
        # position (least-bad score), exactly as get_embedding's candidate fallback.
        best = max(candidates, key=lambda c: c[1])
        return [best[0][: metal_complex.num_atom]]

    return accepted
