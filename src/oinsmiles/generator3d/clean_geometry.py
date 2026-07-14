import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from scipy.spatial.distance import cdist

from .bond_lengths import bond_length, sigma_table_applies
from .embed import _apply_atom_chirality, _apply_double_bond_stereo


def _binding_distance(metal_sym, ligand_sym, metal_r, atom_r, scale, is_sigma):
    """Per-atom metal–donor scan target (Å) for ``ff_clean``.

    σ donors (single-atom binding groups) for which
    :func:`~oinsmiles.generator3d.bond_lengths.sigma_table_applies` is true consult
    the curated :func:`~oinsmiles.generator3d.bond_lengths.bond_length` table,
    falling back to the covalent-radius sum when the pair is unlisted (a flat
    default would be worse than the covalent sum). Haptic groups (η-rings, >1 atom
    per slot), non-enabled metals, and short-pin-exempt pairs (e.g. Pd–P, whose
    short curated target over-strains ``ff_clean`` — see
    ``bond_lengths.SHORT_PIN_EXEMPT_PAIRS``) stay entirely on the covalent path — a
    table entry such as ``Fe: {C: 1.80}`` is a σ carbene/carbonyl distance, not an
    η⁵-Cp centroid distance, and Fe is a metal where the table regresses fidelity
    (see ``ENABLED_METALS``). ``scale`` (the swept pool-diversity factor) always
    multiplies, so a per-pair correction cannot collapse the conformer pool the way
    a global factor would.
    """
    if is_sigma and sigma_table_applies(metal_sym, ligand_sym):
        tabulated = bond_length(metal_sym, ligand_sym)
        base = tabulated if tabulated is not None else metal_r + atom_r
    else:
        base = metal_r + atom_r
    return base * scale


def print_rd_geometry(rd_mol, positions):
    """Print the rd geometry."""
    for i, rd_atom in enumerate(rd_mol.GetAtoms()):
        element = rd_atom.GetSymbol()
        x, y, z = positions[i]
        if abs(x) < 0.0001:
            x = 0.00
        if abs(y) < 0.0001:
            y = 0.00
        if abs(z) < 0.0001:
            z = 0.00
        print_x = f"{x:12.8f}"
        print_y = f"{y:12.8f}"
        print_z = f"{z:12.8f}"
        print(f"{element:<3} {print_x} {print_y} {print_z}")
    print()


class TMCOptimizer:
    """Tmc optimizer."""

    def __init__(
        self,
        step_size=0.1,
        ff_max_iters=200,
        ff_force_tol=1e-4,
        ff_energy_tol=1e-6,
        d_converge=0.05,
        num_relaxation=5,
        default_ff="uff",
    ):
        """Initialize the Tmc optimizer."""
        self.step_size = step_size
        self.num_relaxation = num_relaxation
        self.maximal_displacement = 0.5
        self.ratio_criteria = 0.6
        self.atom_d_criteria = 0.5
        self.bond_criteria = 1.2
        self.d_converge = d_converge  # Scan-level geometry convergence (Angstrom)
        self.default_ff = default_ff
        # RDKit ForceField.Minimize() convergence knobs (defaults match RDKit).
        self.ff_max_iters = ff_max_iters
        self.ff_force_tol = ff_force_tol
        self.ff_energy_tol = ff_energy_tol
        self.fix_value = 20.0
        self.binding_fix_value = 2000.0
        self.chk_file = "scan.chk"
        self.scale_factor = {
            1: 1,
            2: 1.1,
            3: 1.1,
            4: 1.1,
            5: 1.2,
            6: 1.2,
            7: 1.2,
            8: 1.2,
            9: 1.2,
            "else": 1.6,
        }

    def clean_geometry(self, metal_complex, scale=1.0):
        """Clean geometry."""
        print("Embedded geometry ...")
        metal_complex.print_coordinate_list()

        print("FF cleaning ...")
        final_energy = 0.0
        try:
            ff_success, final_energy = self.ff_clean(metal_complex, scale)
        except Exception as e:
            print(f"Internal failure for ff clean ... {e}")
            ff_success = False

        if ff_success:
            print("FF clean success!")
            print("FF cleaned geometry ...")
            metal_complex.print_coordinate_list()
            metal_complex.energy = final_energy
            return True
        else:
            print("FF clean has failed ...")
            print("FF cleaned geometry ...")
            metal_complex.print_coordinate_list()
            return False

    def ff_clean(self, metal_complex, scale=1.0):
        """Ff clean."""
        ligands = metal_complex.ligands
        atom_indices_for_each_ligand = metal_complex.get_atom_indices_for_each_ligand()
        center_atom = metal_complex.center_atom
        metal_r = center_atom.get_radius()
        metal_xyz = center_atom.get_coordinate()

        # Prepare FF setting ...
        rd_mol_list = []
        tmp_positions = []
        scanning_indices = []
        target_values = dict()
        binding_groups_infos = dict()
        ligand_binding_group_infos = dict()
        cnt = 0
        ligand_to_metal = dict()

        step_size = self.step_size
        ratio_criteria = self.ratio_criteria
        atom_d_criteria = self.atom_d_criteria
        bond_criteria = self.bond_criteria
        d_converge = self.d_converge
        fix_value = self.fix_value
        binding_fix_value = self.binding_fix_value

        final_positions = metal_complex.get_position()
        ligand_atom_indices = []
        radius_list = [center_atom.get_radius()]
        ligand_adj_matrices = []
        # Carried C=C stereo, offset into combined_rd_mol's index space.
        combined_stereo_bonds = []
        combined_chiral_centers = []

        # Gather ligand information for the scan ...
        for i in range(len(ligands)):
            ligand = ligands[i]
            ace_mol = ligand.molecule.copy()
            for atom in ligand.molecule.atom_list:
                radius_list.append(atom.get_radius())

            try:
                valid_ace_mol = ace_mol.get_valid_molecule(False)
            except Exception:
                valid_ace_mol = None

            if valid_ace_mol is None:
                valid_ace_mol = ace_mol.get_valid_molecule(False, "xyz2mol")
            if len(valid_ace_mol.atom_list) == 1:
                valid_ace_mol.atom_list[0].set_element("Cl")  # In case H- fail ...
                valid_ace_mol.atom_feature["chg"] = np.array([-1])
                valid_ace_mol.chg = -1
            bo_matrix = valid_ace_mol.get_bo_matrix()
            atom_list = ligand.molecule.atom_list
            if bo_matrix is None:
                valid_ace_mol.initialize()

            period_list, group_list = valid_ace_mol.get_period_group_list()

            # TODO: Change charge in each atom for octet deficient atoms (Remove if bad ...)
            n = len(valid_ace_mol.atom_list)
            for j in range(n):
                if period_list[j] == 1:  # Pass for hydrogen
                    continue
                else:
                    bo = np.sum(bo_matrix[j])
                    chg = valid_ace_mol.atom_feature["chg"][j]
                    g = group_list[j]
                    valence = g + bo - chg  # If valence less than 4 ...
                    if valence < 8:
                        valid_ace_mol.atom_feature["chg"][j] -= 8 - valence

            rd_mol = valid_ace_mol.get_rd_mol()
            # Enforce carried C=C (cis/trans) stereo on the per-ligand FF mol.
            # This is where the double bond genuinely exists: the dummy-metal embed
            # mol drops it to a single bond (PuLP re-perception), so this per-ligand
            # FF mol -- freshly embedded below with the double bond intact -- is the
            # reliable place to make distance geometry reproduce the requested E/Z.
            # Indices are ligand-local here; also record them offset into
            # combined_rd_mol for a post-sanitize re-assert (SanitizeMol can drop
            # double-bond stereo that lacks bond directions).
            _apply_double_bond_stereo(rd_mol, ligand.molecule.stereo_bonds)
            _apply_atom_chirality(rd_mol, getattr(ligand.molecule, "chiral_centers", []))
            for si, sj, stereo, sra, srb in getattr(ligand.molecule, "stereo_bonds", []):
                combined_stereo_bonds.append((cnt + si, cnt + sj, stereo, cnt + sra, cnt + srb))
            for center, nbrs, tag in getattr(ligand.molecule, "chiral_centers", []):
                combined_chiral_centers.append((cnt + center, tuple(cnt + k for k in nbrs), tag))
            rd_mol_list.append(rd_mol)
            atom_indices = atom_indices_for_each_ligand[i]
            binding_infos = ligand.binding_infos  # [[indices, geometric_idx]]

            tmp_indices = []
            for j in range(len(atom_indices)):
                tmp_positions.append(final_positions[atom_indices[j]])
                ligand_to_metal[len(tmp_positions) - 1] = atom_indices[j]
                tmp_indices.append(cnt + j + 1)  # Because metal goes to index 0
            ligand_atom_indices.append(tmp_indices)
            ligand_adj_matrices.append(valid_ace_mol.get_adj_matrix())
            total_binding_groups = []
            for info in binding_infos:
                sum_d = 0
                binding_groups = []
                # A single-atom binding group is a σ donor; >1 atom is an η/haptic
                # group. Only σ donors take the curated per-pair bond-length table.
                is_sigma = len(info[0]) == 1
                for idx in info[0]:
                    binding_groups.append(cnt + idx)
                    final_positions[atom_indices[idx]].tolist()
                    atom_r = atom_list[idx].get_radius()
                    sum_d += _binding_distance(
                        center_atom.get_element(),
                        atom_list[idx].get_element(),
                        metal_r,
                        atom_r,
                        scale,
                        is_sigma,
                    )
                ref_d = sum_d / len(info[0])
                if len(info[0]) < 10:
                    ref_d *= self.scale_factor[
                        len(info[0])
                    ]  # Consider elongation of haptic interaction ...
                else:
                    ref_d *= 1.6
                target_values[tuple(binding_groups)] = ref_d
                binding_groups_infos[tuple(binding_groups)] = len(atom_indices)
                total_binding_groups.append(tuple(binding_groups))
                scanning_indices += binding_groups
            ligand_binding_group_infos[tuple(total_binding_groups)] = list(
                range(cnt, cnt + len(atom_indices))
            )
            cnt += len(atom_indices)

        # Construct original_ligand_adj_matrix ...
        original_ligand_adj_matrix = np.zeros((cnt + 1, cnt + 1))
        for k, atom_indices in enumerate(ligand_atom_indices):
            reduce_function = np.ix_(atom_indices, atom_indices)
            original_ligand_adj_matrix[reduce_function] = ligand_adj_matrices[k]

        combined_rd_mol = rd_mol_list[0]
        for rd_mol in rd_mol_list[1:]:
            combined_rd_mol = Chem.CombineMols(combined_rd_mol, rd_mol)

        Chem.SanitizeMol(combined_rd_mol)
        # Re-assert stereo after sanitize (which can strip double-bond stereo
        # lacking bond directions) so the embed below reproduces the E/Z and the
        # sp3 handedness.
        _apply_double_bond_stereo(combined_rd_mol, combined_stereo_bonds)
        _apply_atom_chirality(combined_rd_mol, combined_chiral_centers)
        # The scan below overwrites every atom position from ``tmp_positions``
        # before the first Minimize (see the GetConformer() loop), so this
        # conformer's coordinates are never used — it exists only to attach the
        # conformer object GetConformer() requires at :GetConformer sites. An
        # explicit empty conformer is ~400× cheaper than ETKDG and, unlike
        # EmbedMolecule (which returns -1 → no conformer → GetConformer() raises →
        # the whole scale/option aborts), never fails on a ligand set that
        # sanitizes but will not embed. Stereo asserted just above is untouched
        # (AddConformer adds only coordinates).
        combined_rd_mol.RemoveAllConformers()
        combined_rd_mol.AddConformer(Chem.Conformer(combined_rd_mol.GetNumAtoms()), assignId=True)

        # Make force fields ...
        # MMFF atom typing (``MMFFGetMoleculeProperties``) depends only on the fixed
        # ligand graph, never on 3D coordinates, so compute it once here and reuse it
        # for every per-iteration MMFF rebuild in the scan below (byte-identical: the
        # same properties object yields the same force field).
        tmp_positions = np.array(tmp_positions)
        mmff = None
        uff = None
        mmff_props = None
        try:
            mmff_props = AllChem.MMFFGetMoleculeProperties(combined_rd_mol)
            mmff = AllChem.MMFFGetMoleculeForceField(combined_rd_mol, mmff_props)
        except Exception:
            pass
        try:
            uff = AllChem.UFFGetMoleculeForceField(combined_rd_mol)
        except Exception:
            pass

        if mmff is None and uff is None:
            print("Force field is not supported ...")
            print(Chem.MolToSmiles(combined_rd_mol))
            return False

        n = len(radius_list)
        R = np.repeat(np.array(radius_list), n).reshape((n, n))
        R = R + R.T

        # Sort binding_groups by number of atoms ... (Small to large)
        list(ligand_binding_group_infos.keys())
        sorted_ligand_binding_groups_list = sorted(
            ligand_binding_group_infos, key=lambda x: len(ligand_binding_group_infos[x])
        )

        final_success = True

        # Constant across the whole scan (fixed graph): O(1) donor lookup + atom count.
        scanning_set = set(scanning_indices)
        n_ff_atoms = combined_rd_mol.GetNumAtoms()

        def _make_constrained_ff(ff_kind):
            """Build one FF for the conformer's current coords, pinning atoms in place.

            Every atom is constrained to its present position (scanning donors held
            softer). The molecular graph is fixed across the scan, so building only the
            FF actually Minimized -- and reusing the coordinate-independent
            ``mmff_props`` -- is byte-identical to the old code that rebuilt BOTH force
            fields every iteration.
            """
            if ff_kind == "mmff":
                ff = AllChem.MMFFGetMoleculeForceField(combined_rd_mol, mmff_props)
            else:
                ff = AllChem.UFFGetMoleculeForceField(combined_rd_mol)
            if ff is None:
                return None
            ff.Initialize()
            for atom_idx in range(n_ff_atoms):
                force_constant = binding_fix_value if atom_idx in scanning_set else fix_value
                if ff_kind == "mmff":
                    ff.MMFFAddPositionConstraint(
                        atom_idx, maxDispl=0.00, forceConstant=force_constant
                    )
                else:
                    ff.UFFAddPositionConstraint(
                        atom_idx, maxDispl=0.00, forceConstant=force_constant
                    )
            return ff

        # Scan by each ligand ...
        for ligand_idx, ligand_binding_groups in enumerate(sorted_ligand_binding_groups_list):
            for k in range(100):
                # Also, get old adj matrix ...
                old_positions = np.copy(tmp_positions)
                positions_with_metal = np.vstack((metal_xyz, old_positions))
                distance_matrix = cdist(positions_with_metal, positions_with_metal)
                np.fill_diagonal(distance_matrix, 1e6)
                ratio_matrix = distance_matrix / R

                old_ligand_adj_matrix = np.where(ratio_matrix < bond_criteria, 1, 0)
                old_ligand_adj_matrix[0, :] = 0
                old_ligand_adj_matrix[:, 0] = 0

                # Translation ...
                abs_delta = 0
                current_binding_indices = []
                for binding_groups in list(ligand_binding_groups):
                    binding_vectors = []
                    for idx in binding_groups:
                        binding_vectors.append(tmp_positions[idx].tolist())
                    binding_vectors = np.array(binding_vectors)
                    ref_d = target_values[tuple(binding_groups)]
                    v = np.mean(binding_vectors, axis=0)
                    d = np.linalg.norm(v)
                    delta_d = d - ref_d
                    current_binding_indices += list(binding_groups)

                    # Adjust to ref_d
                    if delta_d > step_size:
                        delta_d = step_size
                    elif delta_d < -step_size:
                        delta_d = -step_size
                    if abs_delta < abs(delta_d):
                        abs_delta = delta_d
                    for idx in binding_groups:
                        tmp_positions[idx] -= delta_d * v / d

                if k > 0 and abs_delta < d_converge:  # Must perform at least one FF opt
                    break

                ff_success = False
                if len(scanning_indices) < len(tmp_positions):
                    conformer = combined_rd_mol.GetConformer()
                    for i, position in enumerate(tmp_positions):
                        x, y, z = position
                        conformer.SetAtomPosition(i, Point3D(x, y, z))

                    # Build+constrain each FF lazily -- only the one actually Minimized
                    # is constructed (see _make_constrained_ff). In the common case
                    # (default_ff="uff", UFF succeeds first) this skips the entire MMFF
                    # rebuild -- atom typing included -- every iteration.
                    ff_order = (
                        ["uff", "mmff"] if self.default_ff.lower() == "uff" else ["mmff", "uff"]
                    )
                    for ff_kind in ff_order:
                        ff = _make_constrained_ff(ff_kind)
                        if ff is None:
                            continue
                        try:
                            ff.Minimize(
                                maxIts=self.ff_max_iters,
                                forceTol=self.ff_force_tol,
                                energyTol=self.ff_energy_tol,
                            )
                        except Exception:
                            continue

                        conformer = combined_rd_mol.GetConformer()
                        tmp_positions = conformer.GetPositions()

                        # Check validity of the geometry (tmp_positions)
                        # Insert metal at zero ...
                        positions_with_metal = np.vstack((metal_xyz, tmp_positions))
                        distance_matrix = cdist(positions_with_metal, positions_with_metal)
                        np.fill_diagonal(distance_matrix, 1e6)
                        ratio_matrix = distance_matrix / R
                        min_ratio = np.min(ratio_matrix)

                        # Check the Collapse of geometry ...
                        if (
                            min_ratio < ratio_criteria
                            or not np.all(distance_matrix) > atom_d_criteria
                        ):
                            print(
                                "[FF Scan] Atoms are too close ... "
                                "Restoring to the original positions !"
                            )
                            print_rd_geometry(combined_rd_mol, tmp_positions)
                            tmp_positions = old_positions  # Restore to original ...
                            continue

                        # Check the ratio between metal and the binding indices ...
                        tmp_indices = [i + 1 for i in current_binding_indices]
                        min_ratio = np.min(ratio_matrix[0, tmp_indices])
                        min_distance = np.min(distance_matrix[0, tmp_indices])
                        ligand_indices = ligand_binding_group_infos[ligand_binding_groups]
                        tmp_indices = [i + 1 for i in ligand_indices]
                        total_min_ratio = np.min(ratio_matrix[0, tmp_indices])
                        total_min_distance = np.min(distance_matrix[0, tmp_indices])
                        if min_ratio < bond_criteria:
                            if (
                                min_ratio > total_min_ratio + 0.1
                                or min_distance > total_min_distance + 0.2
                            ):  # May change ...
                                print(
                                    "[FF scan] Other atoms are likely to bind to the metal "
                                    "... Using the previous positions !"
                                )
                                print_rd_geometry(combined_rd_mol, tmp_positions)
                                tmp_positions = old_positions  # Restore to original ...
                                continue

                        # Check the distance between ligands ...
                        ligand_adj_matrix = np.where(ratio_matrix < bond_criteria, 1, 0)
                        ligand_adj_matrix[0, :] = 0
                        ligand_adj_matrix[:, 0] = 0
                        delta_matrix = ligand_adj_matrix - old_ligand_adj_matrix
                        formed_bonds = np.stack(np.where(delta_matrix > 0), axis=1).tolist()
                        removed_bonds = np.stack(np.where(delta_matrix < 0), axis=1).tolist()
                        # Compare with the original ligand adj matrix
                        adj_change = False
                        for bond in formed_bonds:
                            s, e = bond
                            if original_ligand_adj_matrix[s][e] == 0:
                                adj_change = True
                                break

                        if not adj_change:
                            for bond in removed_bonds:
                                s, e = bond
                                if original_ligand_adj_matrix[s][e] > 0:
                                    adj_change = True
                                    break

                        if adj_change:
                            print(
                                "[FF Scan] Adjacent matrix has changed ... "
                                "Restoring to the original positions !"
                            )
                            print_rd_geometry(combined_rd_mol, tmp_positions)
                            tmp_positions = old_positions
                            for bond in formed_bonds + removed_bonds:
                                s, e = bond
                                if s < e:
                                    continue
                            continue
                        else:
                            ff_success = True
                            break
                else:
                    ff_success = True

                if not ff_success:
                    final_success = False
                    break
        # check move_dict to determine the success of FF clean
        # Less than int(success_criteria/step_size)+1 should be left for successful clean ...

        # If fine, update final positions
        for i in ligand_to_metal:
            x, y, z = tmp_positions[i]
            final_positions[ligand_to_metal[i]] = [x, y, z]
        # Update ligand ...
        metal_complex.set_position(final_positions)
        final_energy = 0.0
        if final_success:
            try:
                conformer = combined_rd_mol.GetConformer()
                for i, position in enumerate(tmp_positions):
                    x, y, z = position
                    conformer.SetAtomPosition(i, Point3D(x, y, z))
                if self.default_ff.lower() == "uff":
                    final_ff = AllChem.UFFGetMoleculeForceField(combined_rd_mol)
                else:
                    final_ff = AllChem.MMFFGetMoleculeForceField(
                        combined_rd_mol, AllChem.MMFFGetMoleculeProperties(combined_rd_mol)
                    )
                if final_ff is not None:
                    final_ff.Initialize()
                    final_energy = final_ff.CalcEnergy()
            except Exception as e:
                print(f"Failed to calculate final FF energy: {e}")

        return final_success, final_energy
