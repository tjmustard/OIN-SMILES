"""
Metal-complex RMSD computation for coordination-sphere verification.

This module provides RMSD calculation specifically for transition metal complexes,
where only the coordination sphere (metal + directly-bonded donors) is meaningful
for geometric validation. Ligand backbone conformations may differ due to
rotational flexibility and DG sampling, so full-molecule RMSD is not useful.

Algorithm (coordination-sphere RMSD):
1. Locate metal atom in both molecules by atomic number
2. Extract atoms within dynamic cutoff (largest gap in metal-atom distances)
3. Group donor atoms by element (C, N, Cl, etc.)
4. Enumerate all permutations within each element group
5. For each permutation, apply Kabsch rotation to mol2 and compute RMSD
6. Return minimum RMSD across all permutations

This handles ligand reordering (e.g., Ir(ppy)3 with 3 equivalent ppy ligands)
and finds the rotation that optimally aligns the coordination frame.
"""

import itertools

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation


def calculate_tmc_rmsd(mol1, mol2, mol2_bonded=None):
    """
    Calculate RMSD for transition metal complexes using coordination-sphere atoms only.

    Uses bond information when available (mol2_bonded with connectivity), otherwise
    falls back to dynamic cutoff detection on distances.

    Args:
        mol1: RDKit mol object (input structure, usually from XYZ)
        mol2: RDKit mol object (generated structure, usually from XYZ)
        mol2_bonded: Optional RDKit mol with bonds (from Molassembler). If provided,
                     uses bonds to define coordination sphere for mol2.

    Returns:
        float: RMSD of coordination sphere atoms (Angstrom), or 999.0 on failure
    """
    try:
        # Extract coordinates
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()

        coords1 = np.array(
            [conf1.GetAtomPosition(i) for i in range(mol1.GetNumAtoms())], dtype=float
        )
        coords2 = np.array(
            [conf2.GetAtomPosition(i) for i in range(mol2.GetNumAtoms())], dtype=float
        )

        # Find metal atom in both molecules (lookup by atomic number)
        metal_idx1, metal_idx2 = _find_metal_atoms(mol1, mol2)
        if metal_idx1 is None or metal_idx2 is None:
            return 999.0

        # Extract coordination spheres
        sphere1 = _extract_coordination_sphere(mol1, coords1, metal_idx1, use_bonds=False)

        # For mol2, try to use bonded information if available
        mol2_for_sphere = mol2_bonded if mol2_bonded is not None else mol2
        conf2_for_coords = mol2_bonded.GetConformer() if mol2_bonded is not None else conf2
        coords2_for_sphere = (
            np.array(
                [conf2_for_coords.GetAtomPosition(i) for i in range(mol2_for_sphere.GetNumAtoms())],
                dtype=float,
            )
            if mol2_bonded is not None
            else coords2
        )

        sphere2 = _extract_coordination_sphere(
            mol2_for_sphere, coords2_for_sphere, metal_idx2, use_bonds=True
        )

        if not sphere1 or not sphere2:
            return 998.0

        # Check element composition match
        if set(sphere1.keys()) != set(sphere2.keys()):
            return 997.0

        for element in sphere1:
            if len(sphere1[element]) != len(sphere2[element]):
                print(
                    f"[RMSD DEBUG] Count mismatch for {element}: {len(sphere1[element])} vs {len(sphere2[element])}"
                )
                return 996.0

        # Center both at the metal
        metal_pos1 = coords1[metal_idx1]
        metal_pos2 = coords2[metal_idx2]

        sphere1_centered = {k: v - metal_pos1 for k, v in sphere1.items()}
        sphere2_centered = {k: v - metal_pos2 for k, v in sphere2.items()}

        # Find minimum RMSD over all permutations within element groups
        return _compute_permutation_rmsd(sphere1_centered, sphere2_centered)

    except Exception as e:
        print(f"RMSD calculation failed: {e}")
        import traceback

        traceback.print_exc()
        return 995.0


def _find_metal_atoms(mol1, mol2):
    """
    Find the metal atom in each molecule by atomic number.
    Common metals: Fe=26, Pt=78, Pd=46, Ir=77, Ru=44, Rh=45, Co=27, Ni=28

    Returns:
        tuple: (metal_idx_mol1, metal_idx_mol2) or (None, None) if not found
    """
    METAL_ATOMIC_NUMBERS = {
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,  # Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn (3d)
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        48,  # Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd (4d)
        72,
        73,
        74,
        75,
        76,
        77,
        78,
        79,
        80,  # Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg (5d)
        57,
        58,
        59,
        60,  # La, Ce, Pr, Nd (lanthanides)
    }

    metal_idx1 = None
    metal_idx2 = None

    for i in range(mol1.GetNumAtoms()):
        if mol1.GetAtomWithIdx(i).GetAtomicNum() in METAL_ATOMIC_NUMBERS:
            metal_idx1 = i
            break

    for i in range(mol2.GetNumAtoms()):
        if mol2.GetAtomWithIdx(i).GetAtomicNum() in METAL_ATOMIC_NUMBERS:
            metal_idx2 = i
            break

    return metal_idx1, metal_idx2


def _extract_coordination_sphere(mol, coords, metal_idx, use_bonds=True, cutoff=None):
    """
    Extract metal + directly-bonded donor atoms.

    If use_bonds=True and molecule has bonding info, uses direct metal neighbors.
    Otherwise uses dynamic distance cutoff.

    Returns:
        dict: {atom_symbol: [pos1, pos2, ...]} grouped by element
    """
    metal_pos = coords[metal_idx]

    # Try bond-based extraction first (more reliable when bonds are available)
    if use_bonds and mol.GetNumBonds() > 0:
        sphere = {}
        metal_atom = mol.GetAtomWithIdx(metal_idx)

        # Add metal itself
        atom_sym = metal_atom.GetSymbol()
        if atom_sym not in sphere:
            sphere[atom_sym] = []
        sphere[atom_sym].append(coords[metal_idx])

        # Add all atoms bonded to the metal (excluding H)
        for neighbor in metal_atom.GetNeighbors():
            if neighbor.GetSymbol() == "H":  # Skip H atoms
                continue
            neighbor_idx = neighbor.GetIdx()
            atom_sym = neighbor.GetSymbol()
            if atom_sym not in sphere:
                sphere[atom_sym] = []
            sphere[atom_sym].append(coords[neighbor_idx])

        return sphere

    # Fall back to connectivity-derived cutoff (covalent radius sum per atom pair).
    # This avoids a single fixed cutoff that is simultaneously too tight for long
    # M-I bonds (~2.75 Å) and too loose for ring-C backbone atoms (~2.9 Å).
    # Cutoff per ligand atom = r_cov(metal) + r_cov(ligand) + 0.45 Å tolerance.
    COVALENT_RADII = {
        1: 0.31,
        6: 0.76,
        7: 0.71,
        8: 0.66,
        9: 0.57,  # H C N O F
        15: 1.07,
        16: 1.05,
        17: 1.02,
        35: 1.20,
        53: 1.39,  # P S Cl Br I
        22: 1.36,
        23: 1.22,
        24: 1.18,
        25: 1.17,
        26: 1.25,  # Ti V Cr Mn Fe
        27: 1.16,
        28: 1.21,
        29: 1.38,
        30: 1.31,  # Co Ni Cu Zn
        40: 1.45,
        41: 1.34,
        42: 1.30,
        43: 1.27,
        44: 1.25,  # Zr Nb Mo Tc Ru
        45: 1.25,
        46: 1.28,
        47: 1.44,
        48: 1.48,  # Rh Pd Ag Cd
        72: 1.44,
        73: 1.34,
        74: 1.30,
        75: 1.28,
        76: 1.26,  # Hf Ta W Re Os
        77: 1.27,
        78: 1.28,
        79: 1.44,
        80: 1.32,  # Ir Pt Au Hg
    }
    DEFAULT_RADIUS = 1.50  # fallback for unlisted elements
    TOLERANCE = 0.55  # generous tolerance covers bond elongation in generated structures

    metal_atom = mol.GetAtomWithIdx(metal_idx)
    r_metal = COVALENT_RADII.get(metal_atom.GetAtomicNum(), DEFAULT_RADIUS)

    # Extract atoms within per-element bonding cutoff (exclude H)
    sphere = {}

    # Add metal itself
    atom_sym = metal_atom.GetSymbol()
    sphere[atom_sym] = [coords[metal_idx]]

    for i in range(mol.GetNumAtoms()):
        if i == metal_idx:
            continue
        atom = mol.GetAtomWithIdx(i)
        if atom.GetSymbol() == "H":
            continue
        r_ligand = COVALENT_RADII.get(atom.GetAtomicNum(), DEFAULT_RADIUS)
        bond_cutoff = cutoff if cutoff is not None else (r_metal + r_ligand + TOLERANCE)
        dist = np.linalg.norm(coords[i] - metal_pos)
        if dist <= bond_cutoff:
            sym = atom.GetSymbol()
            if sym not in sphere:
                sphere[sym] = []
            sphere[sym].append(coords[i])

    return sphere


def _compute_greedy_rmsd(sphere1_centered, sphere2_centered):
    """
    Compute RMSD using greedy optimal assignment when exhaustive permutation is infeasible.
    Uses scipy's linear_sum_assignment to find the minimum-cost matching.

    Args:
        sphere1_centered: {element: [pos1, pos2, ...]} centered at metal
        sphere2_centered: {element: [pos1, pos2, ...]} centered at metal

    Returns:
        float: RMSD with optimal atom matching
    """
    coords1_matched = []
    coords2_matched = []

    for element in sorted(sphere1_centered.keys()):
        coords1_list = np.array(sphere1_centered[element])
        coords2_list = np.array(sphere2_centered[element])

        # Compute pairwise distances between atoms in the two sets
        # Use distance matrix for optimal assignment
        n1, n2 = len(coords1_list), len(coords2_list)
        dist_matrix = np.zeros((n1, n2))
        for i in range(n1):
            for j in range(n2):
                dist_matrix[i, j] = np.linalg.norm(coords1_list[i] - coords2_list[j])

        # Find optimal assignment
        row_ind, col_ind = linear_sum_assignment(dist_matrix)

        coords1_matched.append(coords1_list[row_ind])
        coords2_matched.append(coords2_list[col_ind])

    coords1_all = np.vstack(coords1_matched)
    coords2_all = np.vstack(coords2_matched)

    # Compute Kabsch rotation and RMSD
    try:
        rot, rmsd_scipy = Rotation.align_vectors(coords1_all, coords2_all)
        coords2_rotated = rot.apply(coords2_all)
        rmsd = np.sqrt(np.mean(np.sum((coords1_all - coords2_rotated) ** 2, axis=1)))
        return rmsd
    except Exception:
        return 994.0


def _compute_permutation_rmsd(sphere1_centered, sphere2_centered):
    """
    Compute minimum RMSD over all element-group permutations using Kabsch alignment.

    Args:
        sphere1_centered: {element: [pos1, pos2, ...]} centered at metal
        sphere2_centered: {element: [pos1, pos2, ...]} centered at metal

    Returns:
        float: Minimum RMSD across all valid permutations
    """
    # Ensure both have the same element composition
    if set(sphere1_centered.keys()) != set(sphere2_centered.keys()):
        return 993.0

    for element in sphere1_centered:
        if len(sphere1_centered[element]) != len(sphere2_centered[element]):
            return 992.0

    # Generate all permutations within each element group
    # Cap at 6+ atoms per element to avoid combinatorial explosion (6! = 720; 7! = 5040)
    element_perms = {}
    for element in sorted(sphere1_centered.keys()):
        coords1_list = sphere1_centered[element]
        n_atoms = len(coords1_list)
        if n_atoms > 5:
            # Too many atoms for exhaustive permutation; use best-match assignment instead
            return _compute_greedy_rmsd(sphere1_centered, sphere2_centered)
        element_perms[element] = list(itertools.permutations(range(n_atoms)))

    # Cartesian product over all element groups
    element_order = sorted(sphere1_centered.keys())
    perm_iterators = [element_perms[el] for el in element_order]

    min_rmsd = float("inf")

    for perm_combo in itertools.product(*perm_iterators):
        # Build matched atom-pair arrays
        coords1_matched = []
        coords2_matched = []

        for elem_idx, element in enumerate(element_order):
            perm = perm_combo[elem_idx]
            coords1_list = np.array(sphere1_centered[element])
            coords2_list = np.array(sphere2_centered[element])

            coords1_matched.append(coords1_list)
            coords2_matched.append(coords2_list[list(perm)])

        coords1_all = np.vstack(coords1_matched)
        coords2_all = np.vstack(coords2_matched)

        # Compute Kabsch rotation and RMSD
        try:
            rot, rmsd_scipy = Rotation.align_vectors(coords1_all, coords2_all)
            coords2_rotated = rot.apply(coords2_all)
            rmsd = np.sqrt(np.mean(np.sum((coords1_all - coords2_rotated) ** 2, axis=1)))

            if rmsd < min_rmsd:
                min_rmsd = rmsd
        except Exception:
            continue

    return min_rmsd if min_rmsd != float("inf") else 991.0
