"""Molassembler 3D generation adapter for OIN3DGenerator.

Replaces the ArchitectorAdapter + ArchitectorWrapper pair with a single
adapter that uses SCINE Molassembler for distance-geometry (DG) conformer
generation.  The C++-backed DG call executes in a subprocess via
ProcessPoolExecutor so that GIL-holding C++ cannot block the Python process.

API surface confirmed in: .agents/memory/molassembler_spike_results.md
  - import scine_molassembler as masm
  - masm.io.experimental.from_smiles(smiles: str) -> masm.Molecule
  - masm.dg.generate_conformation(mol: Molecule, seed: int) -> ndarray | dg.Error
  - masm.io.write(filename: str, mol: Molecule, positions: ndarray)
  - positions are in Angstrom (despite docs saying bohr — verified empirically)

Template-based generator (added for improved geometric fidelity):
  For complexes whose geometry is known (geo_code in TEMPLATES), the primary
  generation path places each ligand by aligning its binding atom(s) to the
  OIN template slot vectors using Kabsch rotation.  This gives near-perfect
  metal–ligand geometry without relying on Molassembler DG quality.  Eta-type
  ligands (multiple binding atoms at the same slot) are routed to Molassembler
  DG, which handles them correctly.
"""
from __future__ import annotations

import os
import tempfile
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import TimeoutError as FuturesTimeout
from dataclasses import dataclass
from typing import Optional

import numpy as np
from rdkit import Chem

from .oin_parser import ParsedOIN, TEMPLATES


@dataclass
class GeneratedStructure:
    """Result of 3D structure generation: XYZ block + optional bonded RDKit mol."""
    xyz: str
    mol: Optional[Chem.Mol] = None


# ===========================================================================
# Typical metal–ligand bond lengths (Angstrom) used by template generator
# ===========================================================================

_BOND_LENGTHS: dict[str, dict[str, float]] = {
    "Ir": {"C": 2.00, "N": 2.12, "O": 2.05, "P": 2.30, "Cl": 2.35, "Br": 2.50, "I": 2.65, "S": 2.30, "F": 1.95},
    "Pt": {"Cl": 2.31, "N": 2.10, "P": 2.25, "C": 2.00, "O": 2.05, "S": 2.30, "Br": 2.40, "I": 2.60, "F": 2.05},
    "Pd": {"Cl": 2.30, "N": 2.10, "P": 2.30, "C": 2.00, "O": 2.05, "S": 2.30, "Br": 2.40, "As": 2.40, "F": 2.10},
    "Fe": {"C": 1.80, "N": 2.00, "O": 2.00, "Cl": 2.20, "Br": 2.35, "I": 2.50, "H": 1.55, "P": 2.20},
    "Cr": {"C": 1.85, "N": 2.05, "O": 1.95, "Cl": 2.30},
    "Re": {"C": 1.90, "N": 2.15, "O": 2.00, "Cl": 2.40},
    "Au": {"C": 2.00, "N": 2.10, "Cl": 2.30, "P": 2.25, "S": 2.30},
    "Rh": {"C": 2.00, "N": 2.10, "O": 2.05, "Cl": 2.35, "P": 2.25},
    "Cu": {"C": 1.90, "N": 2.00, "O": 1.95, "Cl": 2.25},
    "Hg": {"I": 2.65, "Cl": 2.40, "Br": 2.50, "N": 2.30},
    "Ti": {"Cl": 2.30, "N": 2.05, "O": 1.80, "C": 2.05},
    "V":  {"O": 1.60, "N": 2.00, "Cl": 2.30},
    "Zn": {"N": 2.05, "O": 2.00, "Cl": 2.25, "Br": 2.35, "S": 2.25, "Se": 2.40},
    "Ni": {"N": 2.00, "O": 2.00, "Cl": 2.20, "P": 2.20, "C": 1.90, "S": 2.20},
    "Ag": {"N": 2.20, "O": 2.30, "S": 2.40, "Cl": 2.40},
    "Cd": {"N": 2.30, "O": 2.30, "S": 2.55, "Cl": 2.50},
    "Ru": {"C": 1.90, "N": 2.10, "O": 2.00, "Cl": 2.35, "P": 2.25},
    "La": {"N": 2.60, "O": 2.50, "S": 2.80, "Se": 2.90},
}
_DEFAULT_BOND_LENGTH = 2.10  # Fallback for unknown (metal, ligand) pairs


def _bond_length(metal_sym: str, ligand_sym: str) -> float:
    """Return typical bond length (Å) for the given metal–ligand pair."""
    return _BOND_LENGTHS.get(metal_sym, {}).get(ligand_sym, _DEFAULT_BOND_LENGTH)


# ===========================================================================
# Template-based 3D generator (primary + eta-ligand paths)
# ===========================================================================


def _analytic_ring_geometry(
    n: int,
    symbols: list[str],
    cc_bond: float = 1.40,
    ch_bond: float = 1.08,
) -> tuple[np.ndarray, list[str]]:
    """Return positions and symbols for a regular n-membered ring in the xy-plane.

    Ring atoms are placed at angles 0, 2π/n, 4π/n, … at circumradius
    ``cc_bond / (2*sin(π/n))``.  One H is added radially outward from each
    ring atom.  The ring centroid is at the origin; the plane normal is +z.
    """
    circumradius = cc_bond / (2.0 * np.sin(np.pi / n))
    angles = [2.0 * np.pi * k / n for k in range(n)]
    positions: list[np.ndarray] = []
    out_symbols: list[str] = []
    for k, ang in enumerate(angles):
        pos = np.array([circumradius * np.cos(ang), circumradius * np.sin(ang), 0.0])
        positions.append(pos)
        out_symbols.append(symbols[k] if k < len(symbols) else "C")
        # Outward-pointing H
        h_pos = pos * (1.0 + ch_bond / circumradius)
        positions.append(h_pos)
        out_symbols.append("H")
    return np.array(positions, dtype=float), out_symbols


def _stitch_multi_eta_fragment(
    frag_smiles: str,
    vectors: list,
    metal_sym: str,
) -> tuple[np.ndarray, list[str], "Chem.Mol | None"] | None:
    """Place an ansa-metallocene fragment with eta groups at multiple slot directions.

    For fragments where multiple distinct slot directions bind eta atoms (e.g.
    SiMe2-bridged Cp2), groups binding atoms by slot, ETKDG-embeds the full
    organic fragment, then uses Rotation.align_vectors to simultaneously align
    all eta-group centroids to their respective slot units.

    Returns
    -------
    (positions, symbols, mol) or None on failure.
    mol is the RDKit Mol with bond connectivity.
    """
    from rdkit.Chem import AllChem
    from scipy.spatial.transform import Rotation

    slot_groups: dict[tuple, list[int]] = {}
    for v in vectors:
        key = tuple(round(x, 4) for x in v.vector)
        slot_groups.setdefault(key, []).append(v.atom_in_fragment_idx)

    if len(slot_groups) < 2:
        return None

    mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            pass

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    try:
        r = AllChem.EmbedMolecule(mol, params)
    except Exception:
        r = -1
    if r != 0:
        return None

    n_atoms = mol.GetNumAtoms()
    positions = np.array(
        [list(mol.GetConformer().GetAtomPosition(i)) for i in range(n_atoms)],
        dtype=float,
    )
    symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_atoms)]

    source_vecs = []
    target_vecs = []
    for slot_key, bidxs in slot_groups.items():
        valid = [i for i in bidxs if i < n_atoms]
        if not valid:
            return None
        centroid = np.mean([positions[i] for i in valid], axis=0)
        c_norm = np.linalg.norm(centroid)
        if c_norm < 1e-6:
            return None
        source_vecs.append(centroid / c_norm)

        slot_unit = np.array(slot_key, dtype=float)
        slot_norm = np.linalg.norm(slot_unit)
        if slot_norm < 1e-9:
            return None
        slot_unit = slot_unit / slot_norm
        target_vecs.append(slot_unit)

    try:
        rot, _ = Rotation.align_vectors(target_vecs, source_vecs)
    except Exception:
        return None

    rotated = rot.apply(positions)

    T = np.zeros(3)
    for slot_key, bidxs in slot_groups.items():
        valid = [i for i in bidxs if i < n_atoms]
        centroid = np.mean([rotated[i] for i in valid], axis=0)
        slot_unit = np.array(slot_key, dtype=float)
        slot_norm = np.linalg.norm(slot_unit)
        slot_unit = slot_unit / slot_norm

        ring_radius = float(np.mean([np.linalg.norm(rotated[i] - centroid) for i in valid]))
        binding_sym = mol.GetAtomWithIdx(valid[0]).GetSymbol()
        d_mc_bond = _bond_length(metal_sym, binding_sym)
        if d_mc_bond > ring_radius:
            d_mc = float(np.sqrt(d_mc_bond**2 - ring_radius**2))
        else:
            d_mc = d_mc_bond * 0.80
        T += slot_unit * d_mc - centroid

    T /= len(slot_groups)
    final_positions = rotated + T
    return final_positions, symbols, mol


def _stitch_eta_fragment(
    frag_smiles: str,
    binding_idxs: list[int],
    slot_unit: np.ndarray,
    metal_sym: str,
) -> tuple[np.ndarray, list[str], "Chem.Mol | None"] | None:
    """Place an eta-type ligand (Cp, arene) by centroid-plane alignment.

    For eta-n ligands where *binding_idxs* all share the same slot direction:
    1. Tries ETKDGv3 embedding of the organic fragment first.
    2. Falls back to analytic regular-ring geometry if ETKDG fails (e.g. Cp
       anion `[cH]1[cH][cH][cH][cH]1` which RDKit cannot kekulize).
    3. Computes centroid of the binding atoms and estimates M–centroid distance
       from the M–C bond length and the ring circumradius.
    4. Rotates the fragment so its binding-atom plane normal aligns with
       *slot_unit* (pointing away from the metal).
    5. Translates the centroid to ``slot_unit * centroid_dist``.

    Returns
    -------
    (positions, symbols, mol) or None on failure.
    mol is the RDKit Mol with bond connectivity, or None for analytic-geometry fallback.
    """
    from rdkit.Chem import AllChem  # noqa: PLC0415
    from scipy.spatial.transform import Rotation  # noqa: PLC0415

    n_binding = len(binding_idxs)
    if n_binding < 2:
        return None

    positions: np.ndarray | None = None
    symbols: list[str] = []
    valid_idxs: list[int] = []
    etkdg_mol: Chem.Mol | None = None
    smiles_mol: Chem.Mol | None = None

    # ── Attempt 1: ETKDGv3 ───────────────────────────────────────────────────
    mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
    if mol is not None:
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            try:
                Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
            except Exception:
                for atom in mol.GetAtoms():
                    pass  # valence already calculated during partial sanitization
        mol = Chem.AddHs(mol)
        smiles_mol = mol
        _params = AllChem.ETKDGv3()
        _params.randomSeed = 42
        try:
            r = AllChem.EmbedMolecule(mol, _params)
        except Exception:
            r = -1
        if r == 0:
            n_atoms = mol.GetNumAtoms()
            positions = np.array(
                [list(mol.GetConformer().GetAtomPosition(i)) for i in range(n_atoms)],
                dtype=float,
            )
            symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_atoms)]
            valid_idxs = [i for i in binding_idxs if i < n_atoms]
            etkdg_mol = mol

    # ── Attempt 2: analytic regular-ring geometry ────────────────────────────
    if positions is None or len(valid_idxs) < 2:
        # Determine ring atom symbols from the fragment SMILES (best effort).
        mol_ns = Chem.MolFromSmiles(frag_smiles, sanitize=False)
        if mol_ns is not None:
            ring_syms = [
                mol_ns.GetAtomWithIdx(i).GetSymbol()
                for i in range(mol_ns.GetNumAtoms())
                if i in binding_idxs
            ]
        else:
            ring_syms = ["C"] * n_binding

        positions, symbols = _analytic_ring_geometry(n_binding, ring_syms)
        # heavy-atom indices in analytic geometry are 0, 2, 4, … (alternating with H)
        valid_idxs = list(range(0, 2 * n_binding, 2))

    binding_pos = np.array([positions[i] for i in valid_idxs], dtype=float)
    centroid = binding_pos.mean(axis=0)

    # Ring circumradius from binding atoms
    ring_radius = float(np.mean([np.linalg.norm(p - centroid) for p in binding_pos]))
    if ring_radius < 1e-6:
        return None

    # M–C bond length; binding atoms are usually C in Cp/arene ligands.
    # Determine binding atom symbol from the fragment SMILES (best effort).
    _mol_ns_sym = Chem.MolFromSmiles(frag_smiles, sanitize=False)
    if _mol_ns_sym is not None and valid_idxs and valid_idxs[0] < _mol_ns_sym.GetNumAtoms():
        try:
            binding_sym_str = _mol_ns_sym.GetAtomWithIdx(valid_idxs[0]).GetSymbol()
        except Exception:
            binding_sym_str = "C"
    else:
        binding_sym_str = "C"
    m_c_bl = _bond_length(metal_sym, binding_sym_str)

    # M–centroid distance via Pythagorean theorem.
    if m_c_bl > ring_radius:
        centroid_dist = float(np.sqrt(m_c_bl**2 - ring_radius**2))
    else:
        centroid_dist = m_c_bl * 0.80

    target_centroid = slot_unit * centroid_dist

    # Plane normal: smallest SVD singular vector of centred binding positions.
    centered_bp = binding_pos - centroid
    _, _, vh = np.linalg.svd(centered_bp)
    plane_normal = vh[-1]

    # Ensure normal points in the same half-space as slot_unit.
    if float(np.dot(plane_normal, slot_unit)) < 0:
        plane_normal = -plane_normal

    # Rotate fragment so plane_normal → slot_unit (rotation around centroid).
    try:
        rot, _ = Rotation.align_vectors([slot_unit], [plane_normal])
    except Exception:
        return None

    positions = rot.apply(positions - centroid) + centroid
    # Translate centroid to target.
    positions = positions + (target_centroid - centroid)

    # ── Recover bond topology from SMILES mol when ETKDG path failed ─────────
    # When analytic geometry was used, positions are in interleaved [C,H,C,H,…]
    # order (even indices = heavy atoms, odd indices = H).  The smiles_mol from
    # Chem.AddHs() has heavy atoms at indices 0…n_binding-1, then H at
    # n_binding…2*n_binding-1.  Reorder positions to match, then use smiles_mol
    # as the returned mol so that _template_generate can assemble bonds.
    if etkdg_mol is None and smiles_mol is not None and smiles_mol.GetNumAtoms() == len(positions):
        heavy_indices = list(range(0, 2 * n_binding, 2))
        h_indices = list(range(1, 2 * n_binding, 2))
        heavy_pos = positions[heavy_indices]
        h_pos = positions[h_indices]
        positions = np.vstack([heavy_pos, h_pos])
        symbols = [smiles_mol.GetAtomWithIdx(i).GetSymbol() for i in range(smiles_mol.GetNumAtoms())]
        etkdg_mol = smiles_mol

    return positions, symbols, etkdg_mol


def _stitch_fragment(
    frag_smiles: str,
    binding_idxs: list[int],
    target_positions: list[np.ndarray],
    slot_units: list[np.ndarray] | None = None,
    forbidden_positions: list[np.ndarray] | None = None,
) -> tuple[np.ndarray, list[str], "Chem.Mol"] | None:
    """Generate organic fragment and Kabsch-align binding atoms to targets.

    Parameters
    ----------
    frag_smiles:
        SMILES of the ligand fragment (slot markers already stripped).
    binding_idxs:
        Atom indices (in *frag_smiles* parse order) of the binding atoms.
    target_positions:
        Desired 3D positions for each binding atom (in global frame where
        metal is at the origin).
    slot_units:
        Unit vectors from metal toward each binding atom.  Used only for
        monodentate orientation: non-binding atoms are rotated to face away
        from the metal (in the +slot direction) to prevent spurious M-H bonds.

    Returns
    -------
    (positions, symbols) or None on failure.
        positions: (N_atoms, 3) float array in Angstrom.
        symbols: list[str] of element symbols, length N_atoms.

    Compatibility check for bidentate
    ----------------------------------
    For bidentate (and higher) ligands, the Kabsch alignment is only applied
    when the distance between binding atoms in the organic ETKDG conformer
    is within 2 Å of the distance between the target positions.  If the
    distances differ by more, the organic structure cannot be rigidly
    superposed onto the chelate geometry without severe distortion (e.g.
    ppy: organic C…N ≈ 7 Å vs chelated ≈ 2.7 Å).  In that case, return
    ``None`` so the caller falls back to Molassembler DG.
    """
    from rdkit.Chem import AllChem  # noqa: PLC0415
    from scipy.spatial.transform import Rotation  # noqa: PLC0415

    mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
    if mol is None:
        return None

    # Sanitize ligand, but allow unusual valences (e.g. C#O for carbonyl).
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            for atom in mol.GetAtoms():
                pass  # valence already calculated during partial sanitization

    # Set NoImplicit on binding atoms to avoid spurious H addition.
    # When the binding atom forms an M-L bond in the complex its valence is
    # already fully used (e.g. C in C#N: triple bond + 1 metal = 4 = max).
    # Use +2 (not +1) to account for double bonds to the metal (e.g. V=O:
    # O has _bsum=0, _dv=2 → 0+2=2 >= 2 → NoImplicit, preventing H2O gen).
    rw = Chem.RWMol(mol)
    try:
        _pt = Chem.GetPeriodicTable()
        for bidx in binding_idxs:
            if bidx >= rw.GetNumAtoms():
                continue
            _a = rw.GetAtomWithIdx(bidx)
            _dv = _pt.GetDefaultValence(_a.GetSymbol())
            if _dv > 0:
                _bsum = int(sum(b.GetBondTypeAsDouble() for b in _a.GetBonds()))
                if _bsum + 2 >= _dv:
                    _a.SetNoImplicit(True)
        mol = Chem.AddHs(rw.GetMol())
    except Exception:
        mol = Chem.AddHs(mol)

    _etkdg_params = AllChem.ETKDGv3()
    _etkdg_params.randomSeed = 42
    try:
        result = AllChem.EmbedMolecule(mol, _etkdg_params)
    except Exception:
        result = -1
    if result == -1:
        if mol.GetNumAtoms() == 1:
            # Single-atom fragment (e.g. hydride [H]) — ETKDG needs ≥ 2 atoms.
            # Create a trivial conformer at origin; the monodentate translation
            # below will move it to the correct target position.
            conf = Chem.Conformer(1)
            conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
            mol.AddConformer(conf, assignId=True)
        else:
            return None

    conf = mol.GetConformer()
    n_atoms = mol.GetNumAtoms()
    positions = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(n_atoms)], dtype=float
    )
    symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_atoms)]

    # Validate binding indices (after AddHs, n_atoms may be larger than the
    # original heavy-atom count but binding_idxs were computed from heavy atoms
    # so they are still valid as long as < original heavy-atom count).
    for bidx in binding_idxs:
        if bidx >= n_atoms:
            return None

    if len(binding_idxs) == 1:
        # Monodentate: rigid translation so binding atom lands on target.
        t = target_positions[0] - positions[binding_idxs[0]]
        positions = positions + t

        # Orientation fix: rotate fragment around the binding atom so that
        # non-binding atoms face away from the metal (in the +slot direction).
        # Without this, H atoms from NH3 / [CH] groups may end up between the
        # metal and the binding atom, causing xyz2mol to form spurious M-H bonds.
        if slot_units and len(slot_units) >= 1 and n_atoms > 1:
            slot_u = np.array(slot_units[0], dtype=float)
            bidx = binding_idxs[0]
            binding_pos = positions[bidx]
            other_idxs = [i for i in range(n_atoms) if i != bidx]
            if other_idxs:
                other_pos = positions[other_idxs]
                outward = other_pos.mean(axis=0) - binding_pos
                outward_norm = float(np.linalg.norm(outward))
                if outward_norm > 1e-6:
                    outward_unit = outward / outward_norm
                    if float(np.dot(outward_unit, slot_u)) < 0.9:
                        try:
                            rot_o, _ = Rotation.align_vectors(
                                [slot_u], [outward_unit]
                            )
                            centered = positions - binding_pos
                            positions = rot_o.apply(centered) + binding_pos
                        except Exception:
                            pass  # Keep current orientation on failure

        return positions, symbols, mol

    # ── Bidentate or higher ──────────────────────────────────────────────────
    # Compatibility check: if organic binding-atom distances differ from target
    # distances by more than 2 Å, Kabsch would severely distort the ligand.
    current_ba = np.array([positions[i] for i in binding_idxs], dtype=float)
    target_ba = np.array(target_positions, dtype=float)

    for i in range(len(binding_idxs)):
        for j in range(i + 1, len(binding_idxs)):
            org_dist = float(np.linalg.norm(current_ba[i] - current_ba[j]))
            tgt_dist = float(np.linalg.norm(target_ba[i] - target_ba[j]))
            if abs(org_dist - tgt_dist) > 2.0:
                return None  # Incompatible geometry → caller uses DG

    c_center = current_ba.mean(axis=0)
    t_center = target_ba.mean(axis=0)
    c_centered = current_ba - c_center
    t_centered = target_ba - t_center

    if np.linalg.norm(c_centered[0]) < 1e-6 or np.linalg.norm(t_centered[0]) < 1e-6:
        return positions + (t_center - c_center), symbols, mol

    rot, _ = Rotation.align_vectors(t_centered, c_centered)
    positions_aligned = rot.apply(positions - c_center) + t_center

    # ── Bite-axis rotation optimisation (bidentate/polydentate) ──────────────
    # scipy's align_vectors with ≤2 vectors in 3D leaves one rotation DOF
    # (around the bite axis) undetermined.  It picks an arbitrary rotation,
    # which can place chelate backbone atoms or H atoms dangerously close to
    # the metal centre (origin) or to already-placed ligand atoms.
    #
    # Fix: search 360° in 5° increments around the bite axis and keep the
    # orientation that (1) avoids clashes with forbidden_positions (already-
    # placed atoms) and (2) maximises min distance from the metal.
    # Clash threshold: 1.5 Å between any non-binding fragment atom and any
    # already-placed atom.  If all angles clash, fall back to best by dist.
    bite_axis = t_centered[0] - t_centered[-1]
    bite_norm = float(np.linalg.norm(bite_axis))
    if len(binding_idxs) >= 2 and bite_norm > 1e-6:
        bite_axis_unit = bite_axis / bite_norm

        # Precompute non-binding atom indices and forbidden positions array.
        _binding_set = set(binding_idxs)
        _nb_idxs = [i for i in range(n_atoms) if i not in _binding_set]
        _forb_np = (
            np.array(forbidden_positions, dtype=float)
            if forbidden_positions
            else None
        )

        def _has_clash(pos_arr: np.ndarray) -> bool:
            """Return True if any non-binding atom is within 1.5 Å of any
            already-placed atom (forbidden_positions)."""
            if _forb_np is None or len(_nb_idxs) == 0:
                return False
            nb = pos_arr[_nb_idxs]
            dists = np.sqrt(
                ((nb[:, None, :] - _forb_np[None, :, :]) ** 2).sum(axis=-1)
            )
            return bool(dists.min() < 1.5)

        centered = positions_aligned - t_center

        initial_min_dist = float(np.linalg.norm(positions_aligned, axis=1).min())
        initial_no_clash = not _has_clash(positions_aligned)

        best_angle = 0.0
        best_no_clash = initial_no_clash
        best_no_clash_dist = initial_min_dist if initial_no_clash else -1.0
        best_any_dist = initial_min_dist  # Fallback when all angles clash

        for deg in range(5, 360, 5):
            rot_try = Rotation.from_rotvec(bite_axis_unit * np.radians(deg))
            pos_try = rot_try.apply(centered) + t_center
            min_dist = float(np.linalg.norm(pos_try, axis=1).min())
            no_clash = not _has_clash(pos_try)

            if no_clash:
                if not best_no_clash or min_dist > best_no_clash_dist:
                    best_no_clash = True
                    best_no_clash_dist = min_dist
                    best_angle = float(np.radians(deg))
            elif not best_no_clash and min_dist > best_any_dist:
                best_any_dist = min_dist
                best_angle = float(np.radians(deg))

        if best_angle != 0.0:
            rot_best = Rotation.from_rotvec(bite_axis_unit * best_angle)
            positions_aligned = rot_best.apply(centered) + t_center

    # Reject bidentate placement if any non-binding heavy atom ends up within
    # 1.7 Å of the metal centre (origin).  This catches cases where Kabsch
    # alignment with a large bite-distance mismatch (e.g. ppy: free C–N ~4.5 Å
    # vs chelated ~2.83 Å) folds inner ring atoms into the metal sphere, causing
    # XYZToSMILES to misidentify the topology on the round-trip.
    if len(binding_idxs) >= 2:
        _binding_set = set(binding_idxs)
        _nb_heavy_idxs = [
            i for i in range(n_atoms)
            if i not in _binding_set and symbols[i] != "H"
        ]
        if _nb_heavy_idxs:
            _min_d_metal = float(
                np.linalg.norm(positions_aligned[_nb_heavy_idxs], axis=1).min()
            )
            if _min_d_metal < 1.7:
                return None  # Bidentate distortion too severe → fall back to DG

    return positions_aligned, symbols, mol


def _template_generate(parsed_oin: "ParsedOIN") -> "tuple[str, Chem.Mol | None] | None":
    """Generate a 3D XYZ block using OIN template slot vectors.

    Returns the XYZ string on success, or ``None`` if template generation
    cannot be applied (unknown geometry, eta ligands, bidentate ligand with
    incompatible organic geometry, or embedding failure).

    Eta-ligand detection
    --------------------
    Eta-ligand handling
    -------------------
    A fragment is an eta-type ligand when two or more OINVectors share the
    same slot direction (same rounded vector tuple).  These are handled by
    ``_stitch_eta_fragment`` (centroid-plane alignment) rather than the
    atom-by-atom Kabsch alignment used for mono/bidentate ligands.  If
    eta-fragment placement fails, falls back to Molassembler DG.

    Bidentate compatibility
    -----------------------
    For bidentate chelates where the binding-atom distance in the isolated
    organic molecule differs by > 2 Å from the chelated distance (e.g. ppy:
    ~7 Å organic vs ~2.7 Å chelated), ``_stitch_fragment`` returns ``None``
    and this function propagates ``None`` to trigger DG fallback.
    """
    geo_code = parsed_oin.geo_code
    if not geo_code:
        return None

    template = TEMPLATES.get(geo_code)
    if template is None:
        return None  # NON or unknown geometry

    # Extract metal symbol from the metal fragment SMILES.
    metal_frag = parsed_oin.fragments[parsed_oin.metal_fragment_idx]
    metal_mol = Chem.MolFromSmiles(metal_frag, sanitize=False)
    if metal_mol is None:
        return None
    metal_sym = metal_mol.GetAtomWithIdx(0).GetSymbol()

    # Group OINVectors by fragment index.
    frag_vecs: dict[int, list] = {}
    for vec in parsed_oin.vectors:
        frag_vecs.setdefault(vec.fragment_idx, []).append(vec)

    # Build the XYZ atom list.
    all_syms: list[str] = [metal_sym]
    all_pos: list[np.ndarray] = [np.zeros(3)]
    all_frag_idxs: list[int] = [-1]  # -1 = metal centre
    # Track eta fragments for post-placement ring-rotation optimisation:
    #   [(atom_start_idx, atom_end_idx_exclusive, slot_unit_vector), …]
    eta_frag_ranges: list[tuple[int, int, np.ndarray]] = []
    # Track fragment mol objects and their binding indices for combined mol assembly.
    # Each entry: (frag_mol, binding_idxs_in_frag)
    fragment_mol_parts: list[tuple[Chem.Mol, list[int]]] = []
    has_all_mols: bool = True
    has_multi_eta: bool = False  # True if any fragment has eta at multiple slots

    for frag_idx, frag_smiles in enumerate(parsed_oin.fragments):
        if frag_idx == parsed_oin.metal_fragment_idx:
            continue

        vecs = frag_vecs.get(frag_idx)
        if not vecs:
            continue  # Uncoordinated ligand — skip.

        # ── Eta-ligand detection (per fragment) ──────────────────────────────
        vec_tuples = [tuple(round(x, 4) for x in v.vector) for v in vecs]
        is_eta = len(vec_tuples) != len(set(vec_tuples))

        if is_eta:
            # All binding atoms share a single slot direction (eta-n Cp/arene).
            # Use centroid-plane alignment instead of per-atom Kabsch.
            unique_dirs = list({vt: None for vt in vec_tuples}.keys())
            if len(unique_dirs) != 1:
                # Multi-eta-slot fragment (ansa-metallocene: one fragment,
                # multiple eta groups at different slots).
                has_multi_eta = True
                result = _stitch_multi_eta_fragment(
                    frag_smiles,
                    vecs,
                    metal_sym,
                )
                if result is None:
                    return None
                frag_positions, frag_symbols, frag_mol = result
                all_binding_idxs = [v.atom_in_fragment_idx for v in vecs]
                eta_start = len(all_pos)
                all_pos.extend(frag_positions)
                all_syms.extend(frag_symbols)
                all_frag_idxs.extend([frag_idx] * len(frag_positions))
                eta_frag_ranges.append((eta_start, len(all_pos), np.array(unique_dirs[0], dtype=float)))
                if frag_mol is not None:
                    fragment_mol_parts.append((frag_mol, all_binding_idxs))
                else:
                    has_all_mols = False
                continue

            eta_slot_vec = np.array(unique_dirs[0], dtype=float)
            eta_slot_norm = float(np.linalg.norm(eta_slot_vec))
            if eta_slot_norm < 1e-9:
                return None
            eta_slot_unit = eta_slot_vec / eta_slot_norm

            eta_binding_idxs = [v.atom_in_fragment_idx for v in vecs]

            result = _stitch_eta_fragment(
                frag_smiles,
                eta_binding_idxs,
                eta_slot_unit,
                metal_sym,
            )
            if result is None:
                return None  # Eta fragment failed → use DG
            frag_positions, frag_symbols, frag_mol = result
            eta_start = len(all_pos)
            all_pos.extend(frag_positions)
            all_syms.extend(frag_symbols)
            all_frag_idxs.extend([frag_idx] * len(frag_positions))
            eta_frag_ranges.append((eta_start, len(all_pos), eta_slot_unit))
            if frag_mol is not None:
                fragment_mol_parts.append((frag_mol, eta_binding_idxs))
            else:
                has_all_mols = False
            continue  # Next fragment

        # ── Normal (monodentate / bidentate) path ────────────────────────────
        # Collect binding atom indices and their target 3D positions.
        binding_idxs: list[int] = []
        target_positions: list[np.ndarray] = []
        slot_units_list: list[np.ndarray] = []
        for v in vecs:
            slot_vec = np.array(v.vector, dtype=float)
            slot_norm = float(np.linalg.norm(slot_vec))
            if slot_norm < 1e-9:
                continue
            slot_unit = slot_vec / slot_norm
            binding_mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
            if binding_mol is None or v.atom_in_fragment_idx >= binding_mol.GetNumAtoms():
                continue
            binding_sym = binding_mol.GetAtomWithIdx(v.atom_in_fragment_idx).GetSymbol()
            bl = _bond_length(metal_sym, binding_sym)
            binding_idxs.append(v.atom_in_fragment_idx)
            target_positions.append(slot_unit * bl)
            slot_units_list.append(slot_unit)

        if not binding_idxs:
            continue

        result = _stitch_fragment(
            frag_smiles,
            binding_idxs,
            target_positions,
            slot_units=slot_units_list,
            forbidden_positions=list(all_pos),
        )
        if result is None:
            # Fragment incompatible with template (e.g. bidentate with large
            # organic vs chelated binding-atom distance).  Fall back to DG.
            return None

        frag_positions, frag_symbols, frag_mol = result
        all_pos.extend(frag_positions)
        all_syms.extend(frag_symbols)
        all_frag_idxs.extend([frag_idx] * len(frag_positions))
        fragment_mol_parts.append((frag_mol, list(binding_idxs)))

    n = len(all_syms)
    if n < 2:
        return None  # Only metal — not useful.

    # ── Inter-fragment collision check (with eta ring-rotation optimisation) ──
    # For eta-ligands placed by centroid-plane alignment the SVD-based rotation
    # picks an arbitrary ring-rotation angle around the metal→centroid axis.
    # When two eta rings collide (inter-ring C-C < 1.60 Å, e.g. TET bent
    # metallocene Cp2TiMe2) we attempt to recover by rotating each eta ring
    # around its centroid–metal axis to maximise minimum inter-fragment distance.
    # Intra-fragment pairs (same frag_idx → bonded atoms, C-C ~1.40–1.54 Å)
    # are excluded from the check.
    def _inter_frag_min(pos_list: list, frag_list: list) -> float:
        heavy = [(i, s) for i, s in enumerate(all_syms) if i > 0 and s != "H"]
        if len(heavy) < 2:
            return 999.0
        hi = [i for i, _ in heavy]
        hp = np.array([pos_list[i] for i in hi], dtype=float)
        hf = np.array([frag_list[i] for i in hi])
        diffs = hp[:, None, :] - hp[None, :, :]
        dists = np.sqrt((diffs**2).sum(axis=-1))
        np.fill_diagonal(dists, 999.0)
        same = hf[:, None] == hf[None, :]
        return float(np.where(same, 999.0, dists).min())

    if eta_frag_ranges:
        # Attempt ring-rotation optimisation for each eta fragment in turn.
        # Always run — not just when clashes exist — because the analytic ring
        # geometry uses an arbitrary rotational phase around the metal→centroid
        # axis.  Optimising to maximise inter-fragment separation gives the most
        # physically reasonable placement even when no hard clash is present
        # (e.g. TiCp2Me2 TET geometry where Cp rings and methyl groups start far
        # apart but the ring phase still strongly affects RMSD).
        from scipy.spatial.transform import Rotation as _Rot  # noqa: PLC0415

        for eta_start, eta_end, slot_u in eta_frag_ranges:
            best_angle = 0.0
            best_min = -1.0
            eta_atoms = list(range(eta_start, eta_end))
            centroid = np.mean([all_pos[i] for i in eta_atoms], axis=0)
            orig_rel = [all_pos[i] - centroid for i in eta_atoms]

            for deg in range(0, 360, 5):
                angle_rad = np.radians(deg)
                rot = _Rot.from_rotvec(slot_u * angle_rad)
                rotated = [rot.apply(r) + centroid for r in orig_rel]
                # Temporarily update positions for this eta ring
                old_pos = [all_pos[i] for i in eta_atoms]
                for k, i in enumerate(eta_atoms):
                    all_pos[i] = rotated[k]
                md = _inter_frag_min(all_pos, all_frag_idxs)
                if md > best_min:
                    best_min = md
                    best_angle = angle_rad
                # Restore
                for k, i in enumerate(eta_atoms):
                    all_pos[i] = old_pos[k]

            # Apply the best rotation found
            rot_best = _Rot.from_rotvec(slot_u * best_angle)
            for k, i in enumerate(eta_atoms):
                all_pos[i] = rot_best.apply(orig_rel[k]) + centroid

        # Final collision check after ring-rotation optimisation
        final_min = _inter_frag_min(all_pos, all_frag_idxs)
        if final_min < 1.60:
            return None

    lines = [str(n), f"Template-generated from OIN ({geo_code})"]
    for sym, pos in zip(all_syms, all_pos):
        lines.append(f"{sym:<2}  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}")
    xyz_str = "\n".join(lines) + "\n"

    # Build combined RDKit mol with bond connectivity + 3D conformer from all_pos.
    # CombineMols preserves existing ETKDG conformers from fragments; we strip
    # all of them and set a single conformer from the final all_pos array so
    # the written MOL/SDF has the correct template-placed positions.
    combined_mol: Chem.Mol | None = None
    if has_all_mols and fragment_mol_parts:
        try:
            combined_rw = Chem.RWMol(metal_mol)
            for frag_mol, frag_binding_idxs in fragment_mol_parts:
                frag_start = combined_rw.GetNumAtoms()
                # Strip conformers from fragment so CombineMols doesn't carry
                # old ETKDG geometry into the combined mol's conformer list.
                frag_no_conf = Chem.RWMol(frag_mol)
                frag_no_conf.RemoveAllConformers()
                combined_rw = Chem.RWMol(Chem.CombineMols(combined_rw.GetMol(), frag_no_conf.GetMol()))
                for bidx in frag_binding_idxs:
                    global_bidx = frag_start + bidx
                    if global_bidx < combined_rw.GetNumAtoms():
                        combined_rw.AddBond(0, global_bidx, Chem.BondType.DATIVE)
            # Set the single conformer from the final (collision-checked) all_pos.
            conf = Chem.Conformer(combined_rw.GetNumAtoms())
            for i, pos in enumerate(all_pos):
                conf.SetAtomPosition(i, pos.tolist())
            combined_rw.AddConformer(conf, assignId=True)
            combined_mol = combined_rw.GetMol()
        except Exception:
            combined_mol = None

    return xyz_str, combined_mol


# ===========================================================================
# OIN geo_code → Molassembler shape name mapping
# ===========================================================================

_SHAPE_MAP: dict[str, str] = {
    "LIN": "Line",
    "TPL": "EquilateralTriangle",
    "TET": "Tetrahedron",
    "TPY": "TrigonalPyramid",
    "SPL": "Square",
    "SPY": "SquarePyramid",
    "TBP": "TrigonalBipyramid",
    "OCT": "Octahedron",
    "PBP": "PentagonalBipyramid",
}

# Symbols that are never the metal centre (organic + common non-metals)
_ORGANIC_SYMS: frozenset[str] = frozenset(
    {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I", "H", "Si", "Se", "As"}
)


# ===========================================================================
# Custom exceptions
# ===========================================================================


class MolassemblerTimeoutError(RuntimeError):
    """Raised when Molassembler DG conformer generation exceeds the timeout."""


# ===========================================================================
# Helper: reconstruct bonded RDKit mol from connected SMILES + XYZ block
# (used for the DG path where the subprocess only returns an XYZ string)
# ===========================================================================


def _reconstruct_mol_from_smiles_and_xyz(smiles: str, xyz_block: str) -> "Chem.Mol | None":
    """Build an RDKit mol with bonds (from SMILES) and 3D coords (from XYZ).

    Returns None if the atom counts or element order do not match, or on any
    parsing error.  The caller should treat None as "no bonded mol available"
    and fall back to coordinate-only output.
    """
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            try:
                Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
            except Exception:
                for atom in mol.GetAtoms():
                    pass  # valence already calculated during partial sanitization
        mol = Chem.AddHs(mol)

        lines = xyz_block.strip().split("\n")
        atom_count = int(lines[0].strip())
        positions = []
        xyz_syms = []
        for line in lines[2 : 2 + atom_count]:
            parts = line.split()
            if len(parts) >= 4:
                xyz_syms.append(parts[0])
                positions.append((float(parts[1]), float(parts[2]), float(parts[3])))

        if len(positions) != mol.GetNumAtoms():
            return None

        rdkit_syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]
        if xyz_syms != rdkit_syms:
            return None  # Atom ordering mismatch between masm and RDKit

        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, pos in enumerate(positions):
            conf.SetAtomPosition(i, pos)
        mol.AddConformer(conf, assignId=True)
        return mol
    except Exception:
        return None


# ===========================================================================
# RDKit ETKDG fallback — used when Molassembler DG returns GraphImpossible
# ===========================================================================


def _rdkit_etkdg_fallback(smiles: str) -> dict:
    """Generate an XYZ block via RDKit ETKDGv3 for topologies that Molassembler
    DG cannot embed (e.g. octahedral bidentate complexes).

    Returns the same dict schema as ``_molassembler_worker``.
    """
    from rdkit import Chem  # noqa: PLC0415
    from rdkit.Chem import AllChem  # noqa: PLC0415

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return {"error": "RDKit could not parse SMILES for ETKDG fallback", "ok": False}
    # Try full sanitization first; if that fails (e.g. kekulization issues with
    # aromatic eta ligands), fall back to softer flags that skip problematic steps.
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            # Skip both KEKULIZE (can't kekulize unkekulizable aromatic systems)
            # and SETAROMATICITY (aromatic atoms from SMILES are already set).
            Chem.SanitizeMol(
                mol,
                Chem.SanitizeFlags.SANITIZE_ALL
                ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
            )
        except Exception:
            pass  # best effort; keep whatever state we have

    # Try to add Hs; if it fails (e.g. aromatic C atoms that can't be kekulized),
    # proceed with ETKDG on the structure without explicit Hs.
    try:
        mol = Chem.AddHs(mol)
    except Exception:
        pass  # keep mol as-is; ETKDG can still embed unkekulized aromatic systems

    _fb_params = AllChem.ETKDGv3()
    _fb_params.randomSeed = 42
    try:
        result = AllChem.EmbedMolecule(mol, _fb_params)
    except Exception as e:
        return {"error": f"RDKit ETKDG embedding failed: {e}", "ok": False}
    if result == -1:
        return {"error": "RDKit ETKDG embedding failed", "ok": False}


    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    lines = [str(n), ""]
    for i in range(n):
        sym = mol.GetAtomWithIdx(i).GetSymbol()
        pos = conf.GetAtomPosition(i)
        lines.append(f"{sym:<2}  {pos.x:12.6f}  {pos.y:12.6f}  {pos.z:12.6f}")
    xyz_block = "\n".join(lines) + "\n"
    return {"xyz_block": xyz_block, "ok": True}


# ===========================================================================
# DG strategy helpers — MUST be at module level for ProcessPoolExecutor pickle
# ===========================================================================


def _min_inter_atomic_dist(positions: np.ndarray) -> float:
    """Return minimum pairwise inter-atomic distance (Å). Higher = fewer clashes."""
    if len(positions) < 2:
        return 999.0
    from scipy.spatial.distance import pdist  # noqa: PLC0415
    return float(pdist(positions).min())


def _best_from_ensemble(mol, seed: int, n: int = 10) -> np.ndarray | None:
    """Generate multiple conformers and return the one with best inter-atomic distances.

    Parameters
    ----------
    mol:
        SCINE Molassembler Molecule object with stereopermutators already assigned.
    seed:
        Seed for deterministic ensemble generation.
    n:
        Number of conformers to generate (default 10).

    Returns
    -------
    np.ndarray of shape (N_atoms, 3) or None if all conformers failed.
    """
    import scine_molassembler as masm  # noqa: PLC0415

    results = masm.dg.generate_ensemble(mol, n, seed)
    scored = [
        (pos, _min_inter_atomic_dist(pos))
        for pos in results
        if not isinstance(pos, masm.dg.Error)
    ]
    if not scored:
        return None
    return max(scored, key=lambda x: x[1])[0]


def _best_from_directed(mol, seed: int, max_size: int = 50) -> np.ndarray | None:
    """Generate conformers via exhaustive dihedral enumeration, pick best by distance.

    If ideal_ensemble_size exceeds max_size, falls back to ensemble generation to
    avoid combinatorial explosion.

    Parameters
    ----------
    mol:
        SCINE Molassembler Molecule object with stereopermutators already assigned.
    seed:
        Seed for deterministic enumeration.
    max_size:
        Ensemble cap — if ideal_ensemble_size > max_size, falls back to ensemble.

    Returns
    -------
    np.ndarray of shape (N_atoms, 3) or None if enumeration failed/skipped.
    """
    import scine_molassembler as masm  # noqa: PLC0415

    try:
        gen = masm.DirectedConformerGenerator(mol)
    except Exception:
        return None

    if gen.ideal_ensemble_size > max_size:
        return _best_from_ensemble(mol, seed, n=max_size)

    best: list = []
    best_score: list = [-1.0]

    def _cb(dl, positions):  # noqa: ARG001
        score = _min_inter_atomic_dist(positions)
        if score > best_score[0]:
            best_score[0] = score
            best[:] = [positions]

    try:
        gen.enumerate(_cb, seed)
    except Exception:
        return None
    return best[0] if best else None


# ===========================================================================
# Module-level worker — MUST be at module level for ProcessPoolExecutor pickle
# ===========================================================================


def _molassembler_worker(args: dict) -> dict:
    """Module-level DG worker for ProcessPoolExecutor.

    Parameters (via args dict)
    -------------------------
    smiles : str
        Connected SMILES of the full TMC (metal bonded to all ligands).
    seed : int
        Seed for deterministic DG generation.

    Returns
    -------
    dict
        On success: ``{"xyz_block": str, "ok": True}``
        On failure: ``{"error": str, "ok": False}``
    """
    import scine_molassembler as masm  # noqa: PLC0415

    smiles: str = args["smiles"]
    seed: int = args.get("seed", 42)
    geo_code: str = args.get("geo_code", "")
    perm_idx: int = args.get("perm_idx", 0)
    expected_trans_sym_pairs: list | None = args.get("expected_trans_sym_pairs")
    expected_bindings: list = args.get("expected_bindings", [])
    dg_strategy: str = args.get("dg_strategy", "single")
    ensemble_size: int = args.get("ensemble_size", 10)
    max_directed_size: int = args.get("max_directed_size", 50)

    try:
        mol = masm.io.experimental.from_smiles(smiles)
    except Exception:
        # Aromatic η-ligands (Cp, arene) and other SMILES that Molassembler
        # cannot interpret natively.  Delegate to RDKit ETKDG which handles
        # arbitrary SMILES topologies.
        return _rdkit_etkdg_fallback(smiles)

    # Enforce correct geometry at the metal centre when geo_code is available.
    # Without this, Molassembler defaults all 4-coordinate centres to Tetrahedron.
    result = None  # will be set below
    if geo_code:
        shape_name = _SHAPE_MAP.get(geo_code)
        if shape_name:
            try:
                from rdkit import Chem as _Chem  # noqa: PLC0415

                rdmol = _Chem.MolFromSmiles(smiles, sanitize=False)
                if rdmol is not None:
                    metal_idx = None
                    for i, atom in enumerate(rdmol.GetAtoms()):
                        if atom.GetSymbol() not in _ORGANIC_SYMS:
                            metal_idx = i
                            break
                    if metal_idx is not None:
                        shape = getattr(masm.shapes.Shape, shape_name)
                        # set_shape_at_atom raises if the Molassembler-inferred
                        # coordination number does not match the target shape's
                        # vertex count (e.g. OCT bidentate chelates where
                        # Molassembler assigns 'triangle' instead of
                        # 'octahedron').  Catch this early and fall through to
                        # ETKDG rather than generating with the wrong shape.
                        shape_ok = False
                        try:
                            mol.set_shape_at_atom(metal_idx, shape)
                            shape_ok = True
                        except Exception:
                            pass

                        if not shape_ok:
                            # Shape size mismatch → Molassembler cannot handle
                            # this topology (e.g. bidentate chelate).  Delegate
                            # to ETKDG which handles arbitrary connectivities.
                            result = _rdkit_etkdg_fallback(smiles)
                        else:
                            sp = mol.stereopermutators.option(metal_idx)
                            n_perms = sp.num_stereopermutations if sp is not None else 1

                            if expected_trans_sym_pairs:
                                # Feedback loop: find the stereopermutation(s)
                                # whose 3D geometry matches the OIN topology.
                                #
                                # Phase 1 — trans-pair topology check:
                                #   _check_trans_sym_pairs tests whether the set
                                #   of anti-parallel binding-atom symbol pairs
                                #   matches what the OIN slot vectors predict.
                                #   This distinguishes CIS/TRANS (SPL) and
                                #   fac/mer (OCT) at the topology level.
                                #
                                # Phase 2 — exact slot tiebreaker (binding_info):
                                #   When multiple permutations pass Phase 1 (e.g.
                                #   fac perm 2 and perm 3 for Ir(ppy)3 both have
                                #   the same all-C-N trans pairs), we run
                                #   _check_exact_slot_match to verify each
                                #   binding atom's actual direction is closest to
                                #   its OIN-assigned template slot vector.
                                #   Only triggers when len(candidates) > 1.
                                #
                                # Inner seed retry (seeds seed..seed+4) handles
                                # seed-dependent DG failures such as
                                # RefinedStructureInacceptable for OCT bidentate
                                # complexes (e.g. Ir(ppy)3 perm 0 at seed 42).
                                #
                                # The try/except around assign_stereopermutator
                                # guards against num_stereopermutations > actual
                                # valid range (Molassembler bug for some topologies).
                                candidates: list = []  # (positions,) passing trans check
                                for try_perm in range(n_perms):
                                    try:
                                        mol.assign_stereopermutator(
                                            metal_idx, try_perm
                                        )
                                    except Exception:
                                        break  # no more valid perms
                                    test_res = None
                                    for try_seed in range(seed, seed + 5):
                                        r = masm.dg.generate_conformation(
                                            mol, try_seed
                                        )
                                        if not isinstance(r, masm.dg.Error):
                                            test_res = r
                                            break
                                    if test_res is None:
                                        continue  # all seeds failed this perm
                                    if _check_trans_sym_pairs(
                                        test_res,
                                        smiles,
                                        metal_idx,
                                        expected_trans_sym_pairs,
                                    ):
                                        candidates.append(test_res)

                                if candidates:
                                    if len(candidates) == 1 or not expected_bindings:
                                        # Unambiguous topology match — use it.
                                        result = candidates[0]
                                    else:
                                        # Multiple topology-matching perms —
                                        # use exact slot alignment as tiebreaker.
                                        # _check_exact_slot_match uses element-
                                        # symbol + direction matching so it does
                                        # not require canonical atom-index
                                        # agreement between RDKit and Molassembler.
                                        for cand in candidates:
                                            if _check_exact_slot_match(
                                                cand,
                                                smiles,
                                                metal_idx,
                                                expected_bindings,
                                                geo_code,
                                            ):
                                                result = cand
                                                break
                                        if result is None:
                                            # No exact match — fall back to first
                                            # topology-matching candidate.
                                            result = candidates[0]
                                # If no candidate matched, fall through to
                                # standard generation below (result stays None).
                            elif (
                                expected_bindings
                                and geo_code not in ("SPL",)
                                and n_perms > 1
                            ):
                                # Exact-slot-only feedback for geometries without
                                # anti-parallel trans pairs (e.g. TPY).  SPL is
                                # excluded because _pick_masm_permutation already
                                # handles CIS/TRANS for square-planar complexes.
                                for try_perm in range(n_perms):
                                    try:
                                        mol.assign_stereopermutator(
                                            metal_idx, try_perm
                                        )
                                    except Exception:
                                        break
                                    test_res = None
                                    for try_seed in range(seed, seed + 5):
                                        r = masm.dg.generate_conformation(
                                            mol, try_seed
                                        )
                                        if not isinstance(r, masm.dg.Error):
                                            test_res = r
                                            break
                                    if test_res is None:
                                        continue
                                    if _check_exact_slot_match(
                                        test_res,
                                        smiles,
                                        metal_idx,
                                        expected_bindings,
                                        geo_code,
                                    ):
                                        result = test_res
                                        break
                                    if result is None:
                                        # Keep first valid conformation as fallback
                                        result = test_res
                                # If no perm matched, fall through to standard
                                # generation below (result may already be set
                                # to the first valid conformation as fallback).
                            else:
                                safe_perm = min(perm_idx, n_perms - 1)
                                try:
                                    mol.assign_stereopermutator(
                                        metal_idx, safe_perm
                                    )
                                except Exception:
                                    pass  # invalid perm; proceed with default
            except Exception:
                pass  # shape enforcement is best-effort; proceed with defaults

    if isinstance(result, dict):
        # Already an ETKDG fallback result dict — return it directly.
        return result

    # Apply the chosen DG strategy. For ensemble/directed, replace the perm-finding
    # result with the best conformer from a wider search. For single, fall through to
    # the seed-retry loop below (which handles result is None).
    if dg_strategy == "ensemble":
        better = _best_from_ensemble(mol, seed, ensemble_size)
        if better is not None:
            result = better
    elif dg_strategy == "directed":
        better = _best_from_directed(mol, seed, max_directed_size)
        if better is not None:
            result = better

    if result is None:
        # Try a few seeds to handle seed-dependent DG failures
        # (RefinedStructureInacceptable is not always deterministic).
        for _try_seed in range(seed, seed + 5):
            result = masm.dg.generate_conformation(mol, _try_seed)
            if not isinstance(result, masm.dg.Error):
                break

    if isinstance(result, masm.dg.Error):
        # Fall back to RDKit ETKDG for any Molassembler DG failure.
        # GraphImpossible: topology cannot be embedded (e.g. bidentate chelates).
        # RefinedStructureInacceptable: DG ran but refinement failed (e.g. tight
        # chelate ring systems like OCT Ir(ppy)3).
        return _rdkit_etkdg_fallback(smiles)

    positions = result  # ndarray (N_atoms, 3) in Angstrom

    # Write XYZ to a temp file then read back as a string.
    tmp_path: str = ""
    try:
        with tempfile.NamedTemporaryFile(suffix=".xyz", delete=False) as f:
            tmp_path = f.name
        masm.io.write(tmp_path, mol, positions)
        with open(tmp_path) as f:
            xyz_block = f.read()
    finally:
        if tmp_path:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass

    return {"xyz_block": xyz_block, "ok": True}


# ===========================================================================
# Helper: build connected SMILES from ParsedOIN
# ===========================================================================


def _build_connected_smiles(parsed_oin: ParsedOIN) -> str:
    """Construct a single connected SMILES for Molassembler from ParsedOIN.

    Uses RDKit CombineMols + AddBond to join the metal atom to each ligand's
    binding atom.  Falls back to the dot-disconnected ParsedOIN.smiles if
    the RDKit mol operations fail.

    Notes
    -----
    Binding atom indices come from ``OINVector.atom_in_fragment_idx``.  With
    the current regex-based ``parse_inline_string()`` these are always 0
    (first atom in the fragment SMILES).  This is correct for monodentate
    ligands; for bidentate ligands the second binding atom is also approximated
    as 0 — a known limitation awaiting a richer vector parser.
    """
    metal_smiles = parsed_oin.fragments[parsed_oin.metal_fragment_idx]

    # Use sanitize=False so transition metals with unusual valence don't raise.
    metal_mol = Chem.MolFromSmiles(metal_smiles, sanitize=False)
    if metal_mol is None:
        return parsed_oin.smiles  # fallback

    # Group vectors by fragment index.
    frag_vectors: dict[int, list] = {}
    for vec in parsed_oin.vectors:
        frag_vectors.setdefault(vec.fragment_idx, []).append(vec)

    rw = Chem.RWMol(metal_mol)
    metal_atom_idx = 0  # metal is the sole atom in its fragment mol

    for frag_idx, frag_smiles in enumerate(parsed_oin.fragments):
        if frag_idx == parsed_oin.metal_fragment_idx:
            continue

        vectors = frag_vectors.get(frag_idx)
        if not vectors:
            continue  # skip uncoordinated ligands

        # Load ligand fragment from original SMILES to preserve atom ordering.
        # atom_in_fragment_idx values come from _count_smiles_atoms_before which
        # counts atoms in the OIN inline string in appearance order — identical
        # to the order Chem.MolFromSmiles assigns atom indices.  We must NOT
        # round-trip through Chem.MolToSmiles (which produces canonical SMILES
        # with a different atom ordering), or binding atom indices will point at
        # wrong atoms (e.g. for ppy, idx 11 is pyridine-N in the original SMILES
        # but may be a C in canonical SMILES → Ir bonds to C instead of N).
        lig_mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
        if lig_mol is None:
            continue

        # Full sanitization to compute valences (needed for NoImplicit check
        # and for correct GetBondTypeAsDouble below).
        try:
            Chem.SanitizeMol(lig_mol)
        except Exception:
            try:
                Chem.SanitizeMol(lig_mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
            except Exception:
                for atom in lig_mol.GetAtoms():
                    pass  # valence already calculated during partial sanitization

        # Explicitly add all implicit H before combining so that ligand
        # H counts are not altered when the metal-ligand bond is added.
        # e.g. N (NH3, 3 implicit H) must stay NH3 after Pt-N bond forms.
        # Exception: binding atoms whose valence will be fully consumed by
        # (existing bonds + 1 metal bond) must NOT receive extra H from
        # AddHs, or the combined mol will have a valence overflow
        # (e.g. C in C#N: bond_order 3 + 1 metal = 4 = C default valence
        # → no room for H; without this guard, AddHs adds [CH]#N-Cu).
        try:
            _pt = Chem.GetPeriodicTable()
            lig_rw_pre = Chem.RWMol(lig_mol)
            for _vec in vectors:
                _bidx = _vec.atom_in_fragment_idx
                if _bidx < lig_rw_pre.GetNumAtoms():
                    _a = lig_rw_pre.GetAtomWithIdx(_bidx)
                    _dv = _pt.GetDefaultValence(_a.GetSymbol())
                    if _dv > 0:  # skip variable-valence metals (dv == -1)
                        _bsum = int(sum(
                            b.GetBondTypeAsDouble() for b in _a.GetBonds()
                        ))
                        if _bsum + 1 >= _dv:
                            _a.SetNoImplicit(True)
            lig_mol = Chem.AddHs(lig_rw_pre.GetMol())
        except Exception:
            pass  # keep lig_mol as-is if AddHs fails

        # Kekulize for Molassembler compatibility (aromatic lowercase atoms are
        # not understood by Molassembler).  Done in-place after AddHs — no
        # SMILES round-trip — so atom_in_fragment_idx indices remain valid.
        try:
            Chem.Kekulize(lig_mol, clearAromaticFlags=True)
        except Exception:
            pass  # keep aromatic form if kekulization fails

        lig_atom_offset = rw.GetNumAtoms()
        # Merge ligand atoms into combined mol
        rw = Chem.RWMol(Chem.CombineMols(rw.GetMol(), lig_mol))

        bonded_pairs: set[tuple[int, int]] = set()
        for vec in vectors:
            binding_idx = lig_atom_offset + vec.atom_in_fragment_idx
            pair = (metal_atom_idx, binding_idx)
            if binding_idx < rw.GetNumAtoms() and pair not in bonded_pairs:
                rw.AddBond(metal_atom_idx, binding_idx, Chem.BondType.SINGLE)
                bonded_pairs.add(pair)

    # Minimal sanitization: skip SETAROMATICITY so atoms already kekulized
    # at the fragment level are not re-aromaticized in the combined mol
    # (which would prevent kekulization of the metal-chelate ring system).
    _SKIP_AROM = (
        Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
    )
    Chem.SanitizeMol(rw, _SKIP_AROM, catchErrors=True)
    try:
        smiles = Chem.MolToSmiles(rw.GetMol())
        # Molassembler does not understand dative bond arrows (-> / <-) that
        # RDKit emits when a ligand atom exceeds its standard valence after
        # forming a metal-ligand bond (e.g. NH3 coordinating to Pt).
        # Replace them with plain single bonds before returning.
        smiles = smiles.replace("->", "-").replace("<-", "-")
        return smiles
    except Exception:
        return parsed_oin.smiles  # ultimate fallback


# ===========================================================================
# Helper: get binding atom symbol from a fragment SMILES
# (module-level so it is available in the subprocess)
# ===========================================================================


def _get_binding_sym(frag_smiles: str, atom_idx: int) -> str | None:
    """Return the element symbol of the binding atom in *frag_smiles*."""
    from rdkit import Chem as _Chem  # noqa: PLC0415

    mol = _Chem.MolFromSmiles(frag_smiles, sanitize=False)
    if mol is None or atom_idx >= mol.GetNumAtoms():
        return None
    return mol.GetAtomWithIdx(atom_idx).GetSymbol()


def _check_trans_sym_pairs(
    positions,
    smiles: str,
    metal_idx: int,
    expected_pairs: list,
) -> bool:
    """Return True if the 3D *positions* have the expected trans symbol pairs.

    Computes unit vectors from metal to each heavy neighbour and checks
    antiparallel pairs (dot < -0.5) against *expected_pairs*, a list of
    (symbol_A, symbol_B) tuples.

    Both same-element pairs (C-C, N-N) and cross-element pairs (C-N) are
    detected — this is essential for distinguishing fac from mer in OCT
    bidentate complexes like Ir(ppy)3.
    """
    import numpy as _np  # noqa: PLC0415
    from rdkit import Chem as _Chem  # noqa: PLC0415

    rdmol = _Chem.MolFromSmiles(smiles, sanitize=False)
    if rdmol is None:
        return False

    metal_pos = positions[metal_idx]
    # Build a flat list of (symbol, unit_vector) for every neighbour.
    # Includes H atoms bonded directly to the metal (hydride ligands),
    # which must be detected to distinguish trans-H pairs (e.g. in FeH2).
    nbr_list: list[tuple[str, "_np.ndarray"]] = []
    for nbr in rdmol.GetAtomWithIdx(metal_idx).GetNeighbors():
        sym = nbr.GetSymbol()
        v = positions[nbr.GetIdx()] - metal_pos
        norm = float(_np.linalg.norm(v))
        if norm > 1e-6:
            nbr_list.append((sym, v / norm))

    # Detect all anti-parallel pairs — including same-element (C-C, N-N).
    actual_pairs: set = set()
    for i in range(len(nbr_list)):
        for j in range(i + 1, len(nbr_list)):
            sym_i, vec_i = nbr_list[i]
            sym_j, vec_j = nbr_list[j]
            if float(_np.dot(vec_i, vec_j)) < -0.5:
                actual_pairs.add(tuple(sorted([sym_i, sym_j])))

    expected_set = {tuple(sorted(p)) for p in expected_pairs}
    return actual_pairs == expected_set


def _check_exact_slot_match(
    positions,
    smiles: str,
    metal_idx: int,
    expected_bindings: list[tuple[str, tuple[float, float, float]]],
    geo_code: str,
) -> bool:
    """Return True if actual binding-atom directions best match OIN slot vectors.

    Uses element-symbol matching: for each expected ``(element, slot_vec)``
    in *expected_bindings*, finds the best-matching unassigned heavy neighbour
    of the metal that has the same element symbol and whose actual 3D direction
    is closest to *slot_vec*.  Rejects if any other template slot is a better
    match for that neighbour (margin 0.1, ≈6°).

    This approach is robust to atom-ordering differences between RDKit and
    Molassembler (no canonical atom-index tracking required) and works for any
    geometry with distinct binding-atom element symbols (e.g. C and N in ppy).

    It is strictly more discriminating than ``_check_trans_sym_pairs``: it can
    distinguish orientationally distinct permutations with the same trans-pair
    signature — e.g. fac perm 2 vs fac perm 3 for Ir(ppy)3.

    Parameters
    ----------
    positions:
        Numpy array (N_atoms × 3) of Molassembler DG positions.
    smiles:
        Connected SMILES of the full TMC (same string passed to Molassembler).
    metal_idx:
        Atom index of the metal in *smiles* (and in *positions*).
    expected_bindings:
        List of ``(element_symbol, slot_unit_vector_tuple)`` for every
        coordinated binding atom, derived from the OIN slot assignments.
    geo_code:
        Three-letter OIN geometry code (e.g. ``'OCT'``).

    Returns
    -------
    bool
        ``True`` only if every expected binding atom can be uniquely matched
        to a metal neighbour that is closest to the expected slot.
    """
    import numpy as _np  # noqa: PLC0415
    from rdkit import Chem as _Chem  # noqa: PLC0415

    if not expected_bindings:
        return False

    all_slots = TEMPLATES.get(geo_code)
    if all_slots is None:
        return False

    # Pre-normalise all template slot vectors.
    slot_norms = _np.linalg.norm(all_slots, axis=1, keepdims=True)
    all_slots_norm = all_slots / _np.where(slot_norms > 1e-9, slot_norms, 1.0)

    rdmol = _Chem.MolFromSmiles(smiles, sanitize=False)
    if rdmol is None:
        return False

    metal_pos = positions[metal_idx]

    # Build flat list of (list_index, element, actual_unit_vector) for every
    # neighbour of the metal. Includes H atoms bonded directly to the metal
    # (hydride ligands), which must be detected in slot matching.
    nbr_list: list[tuple[int, str, "_np.ndarray"]] = []
    for nbr in rdmol.GetAtomWithIdx(metal_idx).GetNeighbors():
        sym = nbr.GetSymbol()
        idx = nbr.GetIdx()
        if idx >= len(positions):
            continue
        v = positions[idx] - metal_pos
        v_norm = float(_np.linalg.norm(v))
        if v_norm < 1e-6:
            continue
        nbr_list.append((len(nbr_list), sym, v / v_norm))

    # Optimal assignment (Hungarian algorithm) grouped by element.
    #
    # Greedy matching fails for chelate ligands (e.g. ppy in Ir(ppy)3) where
    # a binding atom sits between two adjacent OCT template slots: the greedy
    # pick assigns it to whichever expected slot happens to come first, leaving
    # the other expected slot without a good match.  The optimal assignment
    # simultaneously maximises the total dot-product sum, ensuring the globally
    # best correspondence is found before the per-pair slot check is applied.
    from collections import defaultdict as _defaultdict  # noqa: PLC0415
    from scipy.optimize import linear_sum_assignment as _lsa  # noqa: PLC0415

    # Normalise expected slot vectors.
    exp_units: list["_np.ndarray"] = []
    for _, ev in expected_bindings:
        eu = _np.array(ev, dtype=float)
        en = float(_np.linalg.norm(eu))
        exp_units.append(eu / en if en > 1e-9 else eu)

    # Group indices into expected_bindings and nbr_list by element symbol.
    exp_by_elem: dict[str, list[int]] = _defaultdict(list)
    for k, (es, _) in enumerate(expected_bindings):
        exp_by_elem[es].append(k)

    act_by_elem: dict[str, list[tuple[int, "_np.ndarray"]]] = _defaultdict(list)
    for i, sym, actual_unit in nbr_list:
        act_by_elem[sym].append((i, actual_unit))

    # --- Phase 1: element-grouped Hungarian assignment — collect all pairs ---
    all_actual_paired: list["_np.ndarray"] = []
    all_expected_paired: list["_np.ndarray"] = []

    for sym, exp_indices in exp_by_elem.items():
        act_entries = act_by_elem.get(sym, [])
        if len(exp_indices) != len(act_entries):
            return False  # count mismatch

        n = len(exp_indices)
        if n == 0:
            continue

        # Build cost matrix (rows = actual neighbours, cols = expected slots).
        cost = _np.array(
            [
                [float(_np.dot(act_entries[r][1], exp_units[c])) for c in exp_indices]
                for r in range(n)
            ]
        )

        # linear_sum_assignment minimises; negate to maximise total dot product.
        row_ind, col_ind = _lsa(-cost)

        for row, col in zip(row_ind, col_ind):
            all_actual_paired.append(act_entries[row][1])
            all_expected_paired.append(exp_units[exp_indices[col]])

    if not all_actual_paired:
        return False

    # --- Phase 2: Kabsch rotation for frame-invariant slot check ---
    #
    # Molassembler may generate structures in an arbitrary orientation that
    # differs from the OIN template frame.  For example, TrigonalPyramid
    # always generates the apex at −Z while the OIN TPY template places
    # slot 0 at +Z.  A direct dot-product check fails in that case.
    #
    # Fix: find the Kabsch rotation R such that R.apply(expected) ≈ actual,
    # then rotate both the expected vectors and all template slots before
    # the slot-closest check.  The resulting check is frame-invariant.
    all_slots_check = all_slots_norm  # default: original template frame
    expected_check = _np.array(all_expected_paired)

    if len(all_actual_paired) >= 2:
        try:
            from scipy.spatial.transform import Rotation as _Rotation  # noqa: PLC0415

            rot, _ = _Rotation.align_vectors(
                _np.array(all_actual_paired),
                _np.array(all_expected_paired),
            )
            all_slots_check = rot.apply(all_slots_norm)
            expected_check = rot.apply(expected_check)
        except Exception:
            pass  # keep original frame if alignment fails

    # --- Phase 3: slot-closest check in the (possibly rotated) frame ---
    for act_unit, exp_unit in zip(all_actual_paired, expected_check):
        exp_dot = float(_np.dot(act_unit, exp_unit))
        # Reject this permutation if any template slot is a better match
        # than the expected slot (margin 0.1 ≈ 6°; tolerates chelate strain).
        for sv in all_slots_check:
            if float(_np.dot(act_unit, sv)) > exp_dot + 0.1:
                return False

    return True


# ===========================================================================
# Helper: pick stereopermutation index from OIN slot assignments
# ===========================================================================


def _pick_masm_permutation(parsed_oin: ParsedOIN) -> int:
    """Return the Molassembler stereopermutation index to assign at the metal.

    For SPL (square-planar) geometry only: if two ligands with identical
    fragment SMILES are assigned to slots whose template vectors have a
    dot-product ≈ -1 (i.e. *trans* to each other), return 1 (the TRANS
    permutation).  Otherwise return 0 (the CIS / default permutation).

    For all other geometries returns 0 unconditionally.
    """
    if parsed_oin.geo_code != "SPL":
        return 0

    import numpy as np
    from collections import defaultdict

    # Map fragment SMILES → list of template vectors assigned to it.
    smiles_to_vectors: dict[str, list] = defaultdict(list)
    for vec in parsed_oin.vectors:
        frag_smiles = parsed_oin.fragments[vec.fragment_idx]
        smiles_to_vectors[frag_smiles].append(np.array(vec.vector))

    # If any identical-ligand pair has vectors pointing in opposite directions
    # (dot product < -0.5), the complex is TRANS → use permutation 1.
    for vecs in smiles_to_vectors.values():
        if len(vecs) < 2:
            continue
        for i in range(len(vecs)):
            for j in range(i + 1, len(vecs)):
                dot = float(np.dot(vecs[i], vecs[j]))
                if dot < -0.5:
                    return 1  # TRANS arrangement detected

    return 0  # Default: CIS


def _compute_expected_trans_sym_pairs(parsed_oin: ParsedOIN) -> list | None:
    """Compute expected trans atom-symbol pairs from OIN slot vectors.

    Returns a list of (symbol_A, symbol_B) tuples where A and B are the element
    symbols of ligand binding atoms that should be geometrically *trans* to each
    other (their OIN slot vectors have dot product ≈ -1).

    Used to distinguish isomers in the Molassembler stereopermutator feedback
    loop.  For SPL complexes this separates CIS/TRANS for all-different-ligand
    cases.  For OCT bidentate complexes (e.g. Ir(ppy)3) this distinguishes
    fac (all C-N trans pairs) from mer (contains C-C and N-N trans pairs).

    Returns ``None`` if:
    - SPL with identical ligands — handled by ``_pick_masm_permutation`` instead.
    - No anti-parallel slot vectors exist (e.g. TET, TPL, LIN with same element).
    """
    import numpy as _np  # noqa: PLC0415

    geo = parsed_oin.geo_code

    # SPL with identical ligands is handled by _pick_masm_permutation
    # (perm 0 = CIS, perm 1 = TRANS for same-element pairs).
    if geo == "SPL":
        ligand_smiles = [parsed_oin.fragments[v.fragment_idx] for v in parsed_oin.vectors]
        if len(set(ligand_smiles)) < len(ligand_smiles):
            return None

    trans_pairs: list = []
    vecs = parsed_oin.vectors
    for i in range(len(vecs)):
        for j in range(i + 1, len(vecs)):
            dot = float(_np.dot(vecs[i].vector, vecs[j].vector))
            if dot < -0.5:
                sym_i = _get_binding_sym(
                    parsed_oin.fragments[vecs[i].fragment_idx],
                    vecs[i].atom_in_fragment_idx,
                )
                sym_j = _get_binding_sym(
                    parsed_oin.fragments[vecs[j].fragment_idx],
                    vecs[j].atom_in_fragment_idx,
                )
                if sym_i and sym_j:
                    trans_pairs.append(tuple(sorted([sym_i, sym_j])))

    return trans_pairs if trans_pairs else None


# ===========================================================================
# Adapter class
# ===========================================================================


class MolassemblerAdapter:
    """3D structure generator using SCINE Molassembler via ProcessPoolExecutor.

    Replaces the ArchitectorAdapter + ArchitectorWrapper pair.  Calls the
    module-level ``_molassembler_worker`` in an isolated subprocess so that
    C++ GIL-holding code cannot block the main Python process.
    """

    def __init__(
        self,
        timeout: int = 60,
        dg_strategy: str = "single",
        ensemble_size: int = 10,
        max_directed_size: int = 50,
    ) -> None:
        self.timeout = timeout
        self.dg_strategy = dg_strategy
        self.ensemble_size = ensemble_size
        self.max_directed_size = max_directed_size

    def generate(self, parsed_oin: ParsedOIN, seed: int = 42) -> str:
        """Generate a 3D XYZ block for the TMC described by *parsed_oin*.

        Generation strategy (in priority order)
        ----------------------------------------
        1. **Template-based** (primary): places each ligand's binding atoms at
           the OIN slot-vector positions using Kabsch alignment.  Gives near-
           perfect metal–ligand geometry for rigid ligands (CO, ppy, bipy,
           en, …) without DG roughness.  Skipped for eta-type ligands (where
           ≥2 binding atoms share the same slot) and for NON/unknown geometry.
        2. **Molassembler DG** (fallback): used when template generation is
           unavailable (eta ligands, DG needed for stereopermutation selection).

        Parameters
        ----------
        parsed_oin:
            Parsed OIN with fragment SMILES and coordination vectors.
        seed:
            Seed for deterministic DG conformer generation (default 42).

        Returns
        -------
        str
            XYZ block string of the generated conformer.

        Raises
        ------
        MolassemblerTimeoutError
            If DG generation exceeds ``self.timeout`` seconds.
        RuntimeError
            If all generation strategies fail.
        """
        # ── Primary path: template-based Kabsch alignment ──────────────────
        # Bypasses Molassembler DG for all non-eta geometries, giving better
        # geometric fidelity for rigid ligands and eliminating DG roughness
        # for chelate ring systems (e.g. Ir(ppy)3).
        if parsed_oin.geo_code and parsed_oin.geo_code != "NON":
            template_result = _template_generate(parsed_oin)
            if template_result is not None:
                xyz_str, template_mol = template_result
                return GeneratedStructure(xyz=xyz_str, mol=template_mol)

        # ── Fallback: Molassembler DG ───────────────────────────────────────
        connected_smiles = _build_connected_smiles(parsed_oin)
        perm_idx = _pick_masm_permutation(parsed_oin)
        expected_trans_sym_pairs = _compute_expected_trans_sym_pairs(parsed_oin)
        # Compute expected (element, slot_vector) pairs for each coordinated
        # binding atom.  Used by _check_exact_slot_match as a tiebreaker when
        # multiple permutations share the same trans-pair topology.
        expected_bindings: list[tuple[str, tuple[float, float, float]]] = []
        for vec in parsed_oin.vectors:
            if vec.fragment_idx == parsed_oin.metal_fragment_idx:
                continue
            sym = _get_binding_sym(
                parsed_oin.fragments[vec.fragment_idx],
                vec.atom_in_fragment_idx,
            )
            if sym is not None:
                expected_bindings.append((sym, tuple(vec.vector)))  # type: ignore[arg-type]
        args = {
            "smiles": connected_smiles,
            "seed": seed,
            "geo_code": parsed_oin.geo_code,
            "perm_idx": perm_idx,
            "expected_trans_sym_pairs": expected_trans_sym_pairs,
            "expected_bindings": expected_bindings,
            "dg_strategy": self.dg_strategy,
            "ensemble_size": self.ensemble_size,
            "max_directed_size": self.max_directed_size,
        }

        with ProcessPoolExecutor(max_workers=1) as ex:
            fut = ex.submit(_molassembler_worker, args)
            try:
                result = fut.result(timeout=self.timeout)
            except FuturesTimeout:
                raise MolassemblerTimeoutError(
                    f"Molassembler timed out after {self.timeout}s"
                )

        if not result.get("ok"):
            raise RuntimeError(
                f"Molassembler error: {result.get('error', 'unknown')}"
            )

        xyz_block = result["xyz_block"]
        bonded_mol = _reconstruct_mol_from_smiles_and_xyz(connected_smiles, xyz_block)
        return GeneratedStructure(xyz=xyz_block, mol=bonded_mol)
