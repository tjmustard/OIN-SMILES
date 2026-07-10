"""
Metal-complex RMSD computation for coordination-sphere verification.

This module provides RMSD calculation specifically for transition metal complexes,
where only the coordination sphere (metal + directly-bonded donors) is meaningful
for geometric validation. Ligand backbone conformations may differ due to
rotational flexibility and DG sampling, so full-molecule RMSD is not useful.

Algorithm (coordination-sphere RMSD):
1. Locate the metal atom in both molecules (``TRANSITION_METALS_NUM``).
2. Build the *generated* sphere from real bonds: the metal plus its heavy neighbours.
   That connectivity is the chemical ground truth for what coordinates.
3. Select the *input* sphere to match that composition: for each element with count
   k, take the k heavy atoms of that element nearest the metal, rejecting the
   selection when the k-th is beyond a covalent-radius ceiling.
4. Enumerate permutations within each element group (Hungarian + ICP above 5 atoms).
5. For each permutation, apply Kabsch rotation to mol2 and compute the MEAN RMSD.
6. Return the minimum RMSD across all permutations.

Why the input sphere is composition-matched rather than distance-thresholded: the
input mol usually comes from ``Chem.MolFromXYZFile`` and carries no bonds, so a
covalent-radius cutoff has to stand in for connectivity -- and it errs in both
directions. A real long apical Pd-N bond at 2.57 A falls outside a 2.54 A cutoff
(DAPZIF), while a non-donor ring carbon at 2.19 A falls inside one (ROJXIY). No
single threshold fixes both. Matching to the bonded sphere's composition does.

This cannot silently pass a real coordination bug: a genuinely absent donor lies far
outside the ``CEILING_TOL`` ceiling and is reported as a mapping failure. (Callers
such as ``verify_ir_complexes.py`` invoke this without a preceding OIN string gate,
so the ceiling, not the gate, is what makes the selection safe.)

Two entry points:
- ``calculate_tmc_rmsd_detailed`` -> ``(rmsd, None)`` or ``(None, reason)``. Prefer it.
- ``calculate_tmc_rmsd`` -> a bare float, using >=991 sentinel codes for failure.
  Retained for callers that format the result directly.
"""

import itertools
import os
import sys

import numpy as np
from rdkit import Chem
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.core.constants import TRANSITION_METALS_NUM

_PERIODIC_TABLE = Chem.GetPeriodicTable()

# A donor may sit this far beyond r_cov(metal) + r_cov(donor) and still be accepted as
# coordinating. Sized from real long bonds: DAPZIF's apical Pd-N overshoots by 0.47 A,
# ABETIK's Zr-C(allyl) by 0.29 A. A genuinely absent donor overshoots by several A.
CEILING_TOL = 1.0

# Legacy sentinel codes returned by calculate_tmc_rmsd(). Any value >= 900 means
# "the metric could not run", not "the geometry is bad" -- see the module docstring.
_SENTINEL_NO_METAL = 999.0
_SENTINEL_METAL_MISMATCH = 998.0
_SENTINEL_ELEMENT_ABSENT = 997.0
_SENTINEL_MAPPING_FAILED = 996.0
_SENTINEL_EXCEPTION = 995.0


def calculate_tmc_rmsd(mol1, mol2, mol2_bonded=None):
    """
    Calculate RMSD for transition metal complexes using coordination-sphere atoms only.

    Thin wrapper over :func:`calculate_tmc_rmsd_detailed` for callers that consume the
    result as a bare float. Prefer the detailed variant, which distinguishes "the
    geometry is bad" from "the coordination spheres could not be mapped onto each
    other" instead of collapsing both into a number.

    Args:
        mol1: RDKit mol object (input structure, usually from XYZ)
        mol2: RDKit mol object (generated structure, usually from XYZ)
        mol2_bonded: Optional RDKit mol with bonds (from the generator). If provided,
                     uses bonds to define the coordination sphere for mol2.

    Returns:
        float: mean RMSD of coordination-sphere atoms (Angstrom), or a sentinel
        >= 991.0 when the spheres could not be mapped.
    """
    rmsd, _reason, code = _calculate(mol1, mol2, mol2_bonded)
    return rmsd if rmsd is not None else code


def calculate_tmc_rmsd_detailed(mol1, mol2, mol2_bonded=None):
    """
    Coordination-sphere RMSD, distinguishing bad geometry from an unmappable sphere.

    Args:
        mol1: RDKit mol object (input structure, usually from XYZ)
        mol2: RDKit mol object (generated structure, usually from XYZ)
        mol2_bonded: Optional RDKit mol with bonds (from the generator).

    Returns:
        tuple: ``(rmsd, None)`` where rmsd is the mean coordination-sphere RMSD in
        Angstrom, or ``(None, reason)`` when the two spheres cannot be mapped onto
        each other. A ``reason`` is never a statement about geometric quality.
    """
    rmsd, reason, _code = _calculate(mol1, mol2, mol2_bonded)
    return rmsd, reason


def _calculate(mol1, mol2, mol2_bonded=None):
    """Shared implementation. Returns ``(rmsd_or_None, reason_or_None, sentinel_code)``."""
    try:
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()

        coords1 = np.array(
            [conf1.GetAtomPosition(i) for i in range(mol1.GetNumAtoms())], dtype=float
        )

        # The bonded mol, when supplied, is the authority for both connectivity and
        # coordinates -- reading the metal index out of one mol and the positions out
        # of another only works while their atom orderings agree.
        mol2_for_sphere = mol2_bonded if mol2_bonded is not None else mol2
        conf2_for_coords = mol2_for_sphere.GetConformer() if mol2_bonded is not None else conf2
        coords2 = np.array(
            [conf2_for_coords.GetAtomPosition(i) for i in range(mol2_for_sphere.GetNumAtoms())],
            dtype=float,
        )

        metal_idx1 = _find_metal(mol1)
        metal_idx2 = _find_metal(mol2_for_sphere)
        if metal_idx1 is None:
            return None, "no transition metal found in the input structure", _SENTINEL_NO_METAL
        if metal_idx2 is None:
            return None, "no transition metal found in the generated structure", _SENTINEL_NO_METAL

        sym1 = mol1.GetAtomWithIdx(metal_idx1).GetSymbol()
        sym2 = mol2_for_sphere.GetAtomWithIdx(metal_idx2).GetSymbol()
        if sym1 != sym2:
            return (
                None,
                f"metal element differs: input {sym1}, generated {sym2}",
                _SENTINEL_METAL_MISMATCH,
            )

        # Generated side: real bonds when we have them, distance cutoff otherwise.
        sphere2 = _extract_coordination_sphere(mol2_for_sphere, coords2, metal_idx2, use_bonds=True)
        if not sphere2:
            return None, "generated coordination sphere is empty", _SENTINEL_METAL_MISMATCH

        # Input side: selected to match the generated sphere's composition.
        composition = {el: len(pos) for el, pos in sphere2.items()}
        sphere1, reason, code = _select_matched_sphere(mol1, coords1, metal_idx1, composition)
        if sphere1 is None:
            return None, reason, code

        metal_pos1 = coords1[metal_idx1]
        metal_pos2 = coords2[metal_idx2]

        sphere1_centered = {k: np.asarray(v) - metal_pos1 for k, v in sphere1.items()}
        sphere2_centered = {k: np.asarray(v) - metal_pos2 for k, v in sphere2.items()}

        rmsd = _compute_permutation_rmsd(sphere1_centered, sphere2_centered)
        if rmsd >= 900:
            return None, f"RMSD kernel found no valid alignment (code {rmsd:.0f})", rmsd
        return rmsd, None, 0.0

    except Exception as e:
        print(f"RMSD calculation failed: {e}")
        import traceback

        traceback.print_exc()
        return None, f"exception in RMSD metric: {type(e).__name__}: {e}", _SENTINEL_EXCEPTION


def _rcov(atomic_num):
    """Covalent radius in Angstrom, from RDKit's periodic table."""
    return _PERIODIC_TABLE.GetRcovalent(atomic_num)


def _find_metal(mol):
    """Index of the first transition-metal atom in ``mol``, or None.

    The metal list is imported, never copied: a second copy is exactly how Sc and Y
    went missing here and turned 32 correct round-trips into "RMSD 999" (TD-005).
    """
    for i in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(i).GetAtomicNum() in TRANSITION_METALS_NUM:
            return i
    return None


def _find_metal_atoms(mol1, mol2):
    """
    Find the metal atom in each molecule by atomic number.

    Returns:
        tuple: (metal_idx_mol1, metal_idx_mol2); either element is None if not found.
    """
    return _find_metal(mol1), _find_metal(mol2)


def _select_matched_sphere(mol, coords, metal_idx, composition):
    """
    Select an input-side coordination sphere matching ``composition`` element-for-element.

    For each element with count k, take the k heavy atoms of that element nearest the
    metal. Reject when the k-th nearest lies beyond r_cov(metal) + r_cov(el) +
    CEILING_TOL, which is what keeps a genuinely missing donor from being papered over
    with some distant atom of the right element.

    Returns:
        tuple: ``(sphere, None, 0.0)`` or ``(None, reason, sentinel_code)``.
    """
    metal_atom = mol.GetAtomWithIdx(metal_idx)
    metal_sym = metal_atom.GetSymbol()
    metal_pos = coords[metal_idx]
    r_metal = _rcov(metal_atom.GetAtomicNum())

    # Bucket every heavy non-metal atom by element, nearest the metal first.
    candidates = {}
    for i in range(mol.GetNumAtoms()):
        if i == metal_idx:
            continue
        atom = mol.GetAtomWithIdx(i)
        if atom.GetAtomicNum() == 1:
            continue
        dist = float(np.linalg.norm(coords[i] - metal_pos))
        candidates.setdefault(atom.GetSymbol(), []).append((dist, i))
    for bucket in candidates.values():
        bucket.sort()

    sphere = {}
    for element, count in composition.items():
        picked = []
        needed = count

        # The metal is in its own sphere; only extra atoms of the same element compete.
        if element == metal_sym:
            picked.append(metal_pos)
            needed -= 1
        if needed == 0:
            sphere[element] = np.array(picked)
            continue

        available = candidates.get(element, [])
        if len(available) < needed:
            return (
                None,
                f"input has only {len(available)} {element} atom(s) but the generated "
                f"coordination sphere needs {needed}",
                _SENTINEL_ELEMENT_ABSENT,
            )

        ceiling = r_metal + _rcov(_PERIODIC_TABLE.GetAtomicNumber(element)) + CEILING_TOL
        furthest = available[needed - 1][0]
        if furthest > ceiling:
            return (
                None,
                f"{element} donor #{needed} is {furthest:.2f} A from {metal_sym}, "
                f"beyond the {ceiling:.2f} A ceiling",
                _SENTINEL_MAPPING_FAILED,
            )

        picked.extend(coords[i] for _dist, i in available[:needed])
        sphere[element] = np.array(picked)

    return sphere, None, 0.0


def _extract_coordination_sphere(mol, coords, metal_idx, use_bonds=True, cutoff=None):
    """
    Extract metal + directly-bonded donor atoms.

    If use_bonds=True and the molecule has bonding info, uses direct metal neighbours --
    this is the normal path for the generated structure. The distance fallback below is
    reached only when the generator returned no bond topology (``gen_result.mol is
    None``, e.g. some eta cases), where a cutoff is the only option available.

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
    TOLERANCE = 0.55  # generous tolerance covers bond elongation in generated structures

    # Unlisted elements (Sc, Y, La, Lu, B, Si, ...) previously fell back to a blanket
    # 1.50 A, which is wildly wrong in both directions -- too small for Y (1.90), too
    # large for B (0.84). Defer to RDKit rather than extend the table by hand. Listed
    # values are left exactly as they were, so no complex that passes today can shift.
    def radius(atomic_num):
        return COVALENT_RADII.get(atomic_num) or _rcov(atomic_num)

    metal_atom = mol.GetAtomWithIdx(metal_idx)
    r_metal = radius(metal_atom.GetAtomicNum())

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
        r_ligand = radius(atom.GetAtomicNum())
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


def _compute_robust_rmsd(sphere1_centered, sphere2_centered):
    """
    Rotation-robust coordination-sphere RMSD for element groups too large for
    exhaustive permutation (>5 atoms).

    ``_compute_greedy_rmsd`` matches atoms in the *un-rotated* metal-centred frame,
    so it mis-pairs symmetric multi-ring ligands whenever the generated complex is
    rotated relative to the input (e.g. the tilted Cp rings of a bent ansa-
    metallocene) and badly over-estimates the RMSD -- a ferrocene whose generated
    orientation happens to match the input passes, but a tilted TiCat1 does not,
    even when its coordination sphere is geometrically correct.

    This finds the true minimum over both correspondence and rotation:
      1. enumerate candidate rotations from element-compatible anchor-atom pairs,
      2. Hungarian-assign every atom (per element) in each candidate frame,
      3. Kabsch-refine and score, keep the global best,
      4. ICP-polish the best frame to convergence.
    A correct minimum is always <= the greedy estimate, so this can only *lower* a
    complex's RMSD -- it never raises one, hence no regression risk for the
    currently-passing complexes. Falls back to the greedy estimate on any failure.
    """
    try:
        elements = sorted(sphere1_centered.keys())
        P = np.vstack([np.array(sphere1_centered[el]) for el in elements])
        Q = np.vstack([np.array(sphere2_centered[el]) for el in elements])
        labP = np.array([el for el in elements for _ in sphere1_centered[el]])
        labQ = np.array([el for el in elements for _ in sphere2_centered[el]])

        def assign(rot):
            qc = rot.apply(Q)
            colmap = np.empty(len(Q), dtype=int)
            for el in elements:
                pi = np.where(labP == el)[0]
                qi = np.where(labQ == el)[0]
                dist = np.linalg.norm(P[pi][:, None, :] - qc[qi][None, :, :], axis=2)
                r, c = linear_sum_assignment(dist)
                colmap[pi[r]] = qi[c]
            return colmap

        def refine(colmap):
            rot, _ = Rotation.align_vectors(P, Q[colmap])
            qr = rot.apply(Q)
            rmsd = np.sqrt(np.mean(np.sum((P - qr[colmap]) ** 2, axis=1)))
            return rmsd, rot

        def icp(seed, iters=50):
            rot = seed
            prev = None
            rmsd, _ = refine(assign(rot))
            for _ in range(iters):
                colmap = assign(rot)
                rmsd, rot = refine(colmap)
                if prev is not None and np.array_equal(colmap, prev):
                    break
                prev = colmap
            return rmsd

        norms = np.linalg.norm(P, axis=1)
        anchors = [i for i in np.argsort(-norms) if norms[i] > 1e-6]

        best = float("inf")
        best_seed = None
        for ai in anchors[:3]:
            for bi in anchors[:5]:
                if bi == ai:
                    continue
                cosang = abs(np.dot(P[ai], P[bi]) / (norms[ai] * norms[bi] + 1e-12))
                if cosang > 0.97:  # near-collinear anchors -> rotation underdetermined
                    continue
                qa = np.where(labQ == labP[ai])[0]
                qb = np.where(labQ == labP[bi])[0]
                for cj in qa:
                    for dj in qb:
                        if dj == cj:
                            continue
                        rot, _ = Rotation.align_vectors(
                            np.array([P[ai], P[bi]]), np.array([Q[cj], Q[dj]])
                        )
                        rmsd, rref = refine(assign(rot))
                        if rmsd < best:
                            best, best_seed = rmsd, rref
                break  # one non-collinear partner per anchor is sufficient

        greedy = _compute_greedy_rmsd(sphere1_centered, sphere2_centered)
        if best_seed is None:
            return greedy
        # Floor with the legacy greedy estimate: the anchor/ICP search covers a
        # strict superset of correspondences, so min(...) can never exceed greedy
        # and therefore cannot regress a currently-passing complex. (A greedy
        # failure sentinel ~994 is large, so it never wins the min over a real value.)
        return min(best, icp(best_seed), greedy)
    except Exception:
        return _compute_greedy_rmsd(sphere1_centered, sphere2_centered)


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
            # Too many atoms for exhaustive permutation; use a rotation-robust
            # anchor-pair search (greedy single-shot assignment mis-pairs tilted
            # symmetric rings and over-estimates the RMSD).
            return _compute_robust_rmsd(sphere1_centered, sphere2_centered)
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
