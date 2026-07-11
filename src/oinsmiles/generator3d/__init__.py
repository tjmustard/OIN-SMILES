from . import clean_geometry, embed, om


def calculate_heavy_atom_rmsd(mol1, mol2):
    """Calculate the heavy atom rmsd."""
    import numpy as np
    from scipy.spatial.transform import Rotation

    c1 = np.array([a.get_coordinate() for a in mol1.atom_list if a.get_atomic_number() > 1])
    c2 = np.array([a.get_coordinate() for a in mol2.atom_list if a.get_atomic_number() > 1])
    if len(c1) == 0 or len(c1) != len(c2):
        return float("inf")
    c1 -= c1.mean(axis=0)
    c2 -= c2.mean(axis=0)
    try:
        rot, rmsd = Rotation.align_vectors(c1, c2)
        return rmsd
    except Exception:
        return float("inf")


def _dihedral_deg(p0, p1, p2, p3):
    """Signed dihedral (degrees) of the p0-p1-p2-p3 sequence."""
    import math

    import numpy as np

    b1 = np.asarray(p2, float) - np.asarray(p1, float)
    n = np.linalg.norm(b1)
    if n == 0:
        return 0.0
    b1 /= n
    b0 = np.asarray(p0, float) - np.asarray(p1, float)
    b2 = np.asarray(p3, float) - np.asarray(p2, float)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def _complex_stereo_targets(metal_complex):
    """Map each ligand's carried C=C stereo to complex-global atom indices.

    The generator's embed uses a random seed, so a pendant alkene's dihedral is
    not guaranteed to land on the requested side every time (newer RDKit made
    distance geometry honor the constraint less reliably). Returning the target
    tuples lets ``generate_3d_structures`` reject conformers that embedded the
    wrong E/Z. Each tuple is ``(i, j, stereo, ref_a, ref_b)`` in the same index
    space as ``metal_complex.get_position()``.
    """
    targets = []
    try:
        idx_map = metal_complex.get_atom_indices_for_each_ligand()
        ligands = metal_complex.ligands
    except Exception:
        return targets
    for li, ligand in enumerate(ligands):
        if li >= len(idx_map):
            continue
        amap = idx_map[li]
        stereo_bonds = getattr(getattr(ligand, "molecule", None), "stereo_bonds", []) or []
        for si, sj, stereo, sra, srb in stereo_bonds:
            if max(si, sj, sra, srb) >= len(amap):
                continue
            targets.append((amap[si], amap[sj], stereo, amap[sra], amap[srb]))
    return targets


def _complex_chiral_targets(metal_complex):
    """Map each ligand's carried sp3 chirality to complex-global atom indices.

    ``enforceChirality`` makes the embed honor a chiral volume constraint, but a
    distorted or fallback embed can still land on the wrong enantiomer. These
    targets let ``generate_3d_structures`` reject such conformers. Each tuple is
    ``(center, (n0, n1, n2, n3), ChiralType)`` in ``get_position()`` index space.
    """
    targets = []
    try:
        idx_map = metal_complex.get_atom_indices_for_each_ligand()
        ligands = metal_complex.ligands
    except Exception:
        return targets
    for li, ligand in enumerate(ligands):
        if li >= len(idx_map):
            continue
        amap = idx_map[li]
        centers = getattr(getattr(ligand, "molecule", None), "chiral_centers", []) or []
        for center, nbrs, tag in centers:
            if center >= len(amap) or max(nbrs) >= len(amap):
                continue
            targets.append((amap[center], tuple(amap[k] for k in nbrs), tag))
    return targets


def _chiral_targets_satisfied(positions, targets):
    """True if every carried sp3 centre embedded with the requested handedness.

    For neighbours in the order the tag was recorded against, the signed volume
    ``(n1-n0) . ((n2-n0) x (n3-n0))`` is positive for CHI_TETRAHEDRAL_CW and
    negative for CHI_TETRAHEDRAL_CCW (measured against RDKit; pinned by
    ``test_generator_atom_chirality.test_signed_volume_sign_convention``). Like
    the E/Z check this is self-consistent with the stored neighbour order, so it
    needs no CIP perception. A near-planar centre (|volume| tiny) is treated as
    unresolved and rejected rather than guessed.
    """
    import numpy as np
    from rdkit import Chem

    for center, nbrs, tag in targets:
        try:
            p0, p1, p2, p3 = (np.asarray(positions[i], dtype=float) for i in nbrs)
        except Exception:
            continue
        volume = float(np.dot(p1 - p0, np.cross(p2 - p0, p3 - p0)))
        if abs(volume) < 0.1:  # degenerate/planar embed -- not a resolved centre
            return False
        wants_positive = tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW
        if wants_positive != (volume > 0):
            return False
    return True


def _stereo_targets_satisfied(positions, targets):
    """True if every carried C=C stereo is cleanly reproduced by ``positions``.

    CIS/STEREOZ means the two reference atoms sit on the same side of the double
    bond (dihedral ~0); TRANS/STEREOE means opposite sides (dihedral ~180). The
    check is self-consistent with the stored reference atoms, so it needs no CIP
    perception.

    A genuine alkene is planar, so a correct conformer's dihedral is close to 0
    or 180. We accept only clearly-resolved geometry (|dihedral| <= 60 for cis,
    >= 120 for trans) and reject the ambiguous middle band: a distorted embed
    that lands near 90 reads as cis to a naive sign test but is perceived as the
    opposite isomer by downstream CIP (from the H positions), which is exactly
    how a wrong E/Z used to leak through.
    """
    from rdkit import Chem

    for i, j, stereo, ref_a, ref_b in targets:
        try:
            ang = abs(_dihedral_deg(positions[ref_a], positions[i], positions[j], positions[ref_b]))
        except Exception:
            continue
        wants_same_side = stereo in (Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOZ)
        if wants_same_side:
            if ang > 30.0:  # not a clean, near-planar cis
                return False
        else:
            if ang < 150.0:  # not a clean, near-planar trans
                return False
    return True


def generate_3d_structures(
    m_smiles,
    num_conformers=1,
    optimizer=None,
    pool_size=5,
    ff_params=None,
    uff_pool_size=50,
    rmsd_threshold=0.5,
    energy_threshold=2.0,
    timeout=None,
    seed=42,
    embed_time_budget=None,
):
    """Generate 3D structures from an m-SMILES string.

    Returns a list of successful geometries as molecule objects, sorted by energy if an optimizer is
    used.

    ff_params: optional dict of TMCOptimizer convergence knobs
    (ff_max_iters, ff_force_tol, ff_energy_tol, d_converge, num_relaxation).

    seed: base distance-geometry seed. The embed was previously seeded from an
    unseeded ``random.randint``, so repeated runs of the same m-SMILES returned
    different structures -- and every sp3 stereocentre landed on a random
    enantiomer. Defaulting to 42 matches the rest of the project.

    embed_time_budget: optional wall-clock cap (seconds) on the attempt loop. The
    FF-only path had no bound at all -- ``timeout`` was consumed only by the ASE
    optimizer -- so a molecule whose embed never satisfies the validity checks ran
    the full ``max_attempts`` budget (ZIHGEE_comp_0: 250 attempts, ~1696 s) before
    returning nothing. The deadline is checked between attempts (never mid-embed;
    the in-flight attempt always finishes), so a molecule that DOES embed is
    unaffected -- it fills the pool and breaks first. On exhaustion the pool built
    so far is returned; an empty pool becomes the same ``[]`` as before, only fast.
    None (the default) preserves the prior unbounded behavior for direct callers.
    """
    import time

    try:
        metal_complex = om.get_om_from_modified_smiles(m_smiles)
    except Exception as e:
        print(f"Failed to parse m-SMILES: {e}")
        return []

    clean_ff_params = {k: v for k, v in (ff_params or {}).items() if k != "max_attempts"}
    cleaner = clean_geometry.TMCOptimizer(**clean_ff_params)
    options = [0, 1, 2]  # Added one more option to increase pool variety
    scales = [0.8, 0.9, 1.0, 1.1, 1.2]

    # Target number of initial structures to generate
    target_pool = uff_pool_size
    successful_mols = []
    # Carried C=C (E/Z) and sp3 chirality targets, and a pool of otherwise-valid
    # conformers that embedded the wrong stereochemistry (kept only as a
    # last-resort fallback so a stubborn embed never hard-fails generation).
    stereo_targets = _complex_stereo_targets(metal_complex)
    chiral_targets = _complex_chiral_targets(metal_complex)
    stereo_rejects = []

    import itertools

    combinations = list(itertools.product(scales, options))

    default_max = max(target_pool * 5, 250)
    max_attempts = ff_params.get("max_attempts", default_max) if ff_params else default_max

    # Bound the work the stereo filters can add, keyed on REJECTIONS rather than on
    # attempts. A rejected conformer does not fill the pool, so without a cap a
    # molecule whose embed rarely satisfies its constraints never reaches
    # target_pool and runs the FULL attempt budget: AFECIZ (a chelate imine whose
    # C=N the embed seldom lands planar, because useBasicKnowledge=False gives
    # distance geometry no double-bond planarity term) went from 565s unconstrained
    # to >27 min, all of it spent rejecting.
    #
    # Counting rejections, not attempts, leaves healthy molecules untouched -- they
    # reject nothing, so the cap never fires -- and does not penalise a molecule
    # that burns attempts for unrelated reasons (a valence exception, a haptic scale
    # that will not embed). Once the cap is hit we keep whatever satisfied the
    # constraints; if nothing did, the stereo_rejects fallback below returns the same
    # conformer the unfiltered code would have returned, so this is never worse.
    reject_budget = max(target_pool * 2, 25)

    # Wall-clock deadline for the attempt loop (see docstring). Checked between
    # attempts so a molecule that embeds cleanly is never interrupted; it exists
    # only to keep a molecule whose embed never validates from running the full
    # max_attempts budget (ZIHGEE ~1696 s) before giving up.
    deadline = (time.monotonic() + embed_time_budget) if embed_time_budget else None

    for i in range(max_attempts):
        if len(successful_mols) >= target_pool:
            break
        if deadline is not None and time.monotonic() > deadline:
            print(
                f"Embed wall-clock budget ({embed_time_budget}s) reached after "
                f"{i} attempt(s) with {len(successful_mols)} conformer(s) "
                f"({len(stereo_rejects)} stereo-reject fallback(s)); stopping rather "
                f"than exhausting {max_attempts} attempts."
            )
            break
        if len(stereo_rejects) >= reject_budget:
            print(
                f"Stereo reject budget ({reject_budget}) reached with "
                f"{len(successful_mols)} conformer(s) satisfying the requested stereo; "
                f"stopping rather than exhausting {max_attempts} attempts."
            )
            break
        scale, option = combinations[i % len(combinations)]

        # A single scale/option combo can raise inside the embed (e.g. an
        # RDKit valence exception on a dative donor); skip it rather than
        # letting one bad combo abort the whole pool.
        try:
            # Distinct-but-reproducible seed per attempt: the pool keeps its
            # variety, but the same m-SMILES always yields the same conformers.
            # (Stride idiom borrowed from molassembler_adapter's retry loop.)
            positions = embed.get_embedding(
                metal_complex, scale, option, align=True, seed=seed + i * 1009
            )
        except Exception as e:
            print(f"Embedding failed (scale={scale}, option={option}): {e}")
            positions = None
        if positions is not None:
            tmp_complex = metal_complex.copy()
            tmp_complex.set_position(positions)

            # cleaner.clean_geometry will print logs, could be silenced later
            success = cleaner.clean_geometry(tmp_complex, scale)

            if success:
                position = tmp_complex.get_position()
                wrong_ez = stereo_targets and not _stereo_targets_satisfied(
                    position, stereo_targets
                )
                wrong_chirality = chiral_targets and not _chiral_targets_satisfied(
                    position, chiral_targets
                )
                if wrong_ez or wrong_chirality:
                    # Right topology, wrong stereochemistry -- keep as fallback only.
                    stereo_rejects.append(tmp_complex.get_molecule())
                    continue
                successful_mols.append(tmp_complex.get_molecule())

    if not successful_mols and stereo_rejects:
        # No embed reproduced the requested stereochemistry within the attempt
        # budget; return the best available rather than nothing (non-regressive
        # vs. the prior, unfiltered behavior).
        print("WARNING: no conformer reproduced the requested stereo; using best available.")
        successful_mols = stereo_rejects

    if not successful_mols:
        return []

    # Sort by UFF energy if available (handle None values safely)
    successful_mols.sort(
        key=lambda m: (
            getattr(m, "energy", None) if getattr(m, "energy", None) is not None else float("inf")
        )
    )

    # Deduplicate
    dedup_mols = []
    for mol in successful_mols:
        is_unique = True
        for acc_mol in dedup_mols:
            rmsd = calculate_heavy_atom_rmsd(mol, acc_mol)
            e1 = getattr(mol, "energy", None)
            e1 = e1 if e1 is not None else 0.0
            e2 = getattr(acc_mol, "energy", None)
            e2 = e2 if e2 is not None else 0.0
            if rmsd < rmsd_threshold and abs(e1 - e2) <= energy_threshold:
                is_unique = False
                break
        if is_unique:
            dedup_mols.append(mol)

    successful_mols = dedup_mols

    if optimizer:
        from .ml_optimizer import ASEOptimizer

        # We explicitly do NOT catch initialization exceptions here.
        # If the user asks for an optimizer and it fails to load,
        # we want to fail loudly rather than silently falling back to FF.
        opt = ASEOptimizer(method=optimizer, timeout=timeout)

        if opt:
            optimized_mols = []
            mols_to_optimize = successful_mols[:num_conformers]
            for mol in mols_to_optimize:
                success, energy, new_mol = opt.optimize(mol)
                if success:
                    optimized_mols.append((energy, new_mol))
                else:
                    # Keep the original if optimization fails but penalize its rank
                    optimized_mols.append((float("inf"), mol))

            # Sort by energy
            optimized_mols.sort(key=lambda x: x[0])
            successful_mols = [m[1] for m in optimized_mols]

    return successful_mols[:num_conformers]


def get_xyz_string(molecule):
    """Returns the XYZ string format for a given generated molecule."""
    atom_list = molecule.atom_list
    lines = [
        str(len(atom_list)),
        f"Generated by MetalloGen-3D, Charge: {molecule.chg}, "
        f"Multiplicity: {molecule.multiplicity}",
    ]
    for atom in atom_list:
        element = atom.get_element()
        coord = atom.get_coordinate()
        lines.append(f"{element:<3} {coord[0]:>12.8f} {coord[1]:>12.8f} {coord[2]:>12.8f}")
    return "\n".join(lines)
