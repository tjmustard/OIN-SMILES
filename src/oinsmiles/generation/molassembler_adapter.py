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
import sys
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import TimeoutError as FuturesTimeout
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ..core.chirality import (
    OINStereoWarning,
    _attach_dummy_metal,
    _build_dummy_metal_copy,
    _lp_cip_label,
)
from ..core.constants import TRANSITION_METALS_NUM
from ..oin.winding import signed_circulation
from .oin_parser import TEMPLATES, ParsedOIN


@dataclass
class GeneratedStructure:
    """Result of 3D structure generation: XYZ block + optional bonded RDKit mol."""

    xyz: str
    mol: Optional[Chem.Mol] = None
    # Stereo Phase 3: per-ring haptic-face correction decisions (`fired` /
    # `skipped` / `conflict` / `no-op`), one dict per eta ring encountered
    # during template-based generation. Empty for non-eta / DG-fallback
    # generation. Inspectable metadata -- see `_stitch_eta_fragment` and
    # `_stitch_multi_eta_fragment` for the decision vocabulary.
    haptic_face_decisions: list = field(default_factory=list)


# ===========================================================================
# Typical metal–ligand bond lengths (Angstrom) used by template generator
# ===========================================================================

_BOND_LENGTHS: dict[str, dict[str, float]] = {
    "Ir": {
        "C": 2.00,
        "N": 2.12,
        "O": 2.05,
        "P": 2.30,
        "Cl": 2.35,
        "Br": 2.50,
        "I": 2.65,
        "S": 2.30,
        "F": 1.95,
    },
    "Pt": {
        "Cl": 2.31,
        "N": 2.10,
        "P": 2.25,
        "C": 2.00,
        "O": 2.05,
        "S": 2.30,
        "Br": 2.40,
        "I": 2.60,
        "F": 2.05,
    },
    "Pd": {
        "Cl": 2.30,
        "N": 2.10,
        "P": 2.30,
        "C": 2.00,
        "O": 2.05,
        "S": 2.30,
        "Br": 2.40,
        "As": 2.40,
        "F": 2.10,
    },
    "Fe": {
        "C": 1.80,
        "N": 2.00,
        "O": 2.00,
        "Cl": 2.20,
        "Br": 2.35,
        "I": 2.50,
        "H": 1.55,
        "P": 2.20,
    },
    "Cr": {"C": 1.85, "N": 2.05, "O": 1.95, "Cl": 2.30},
    "Re": {"C": 1.90, "N": 2.15, "O": 2.00, "Cl": 2.40},
    "Au": {"C": 2.00, "N": 2.10, "Cl": 2.30, "P": 2.25, "S": 2.30},
    "Rh": {"C": 2.00, "N": 2.10, "O": 2.05, "Cl": 2.35, "P": 2.25},
    "Cu": {"C": 1.90, "N": 2.00, "O": 1.95, "Cl": 2.25},
    "Hg": {"I": 2.65, "Cl": 2.40, "Br": 2.50, "N": 2.30},
    "Ti": {"Cl": 2.30, "N": 2.05, "O": 1.80, "C": 2.05},
    "V": {"O": 1.60, "N": 2.00, "Cl": 2.30},
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


# ===========================================================================
# Stereo Phase 3: haptic-face correction helpers
# ===========================================================================
#
# Shared by `_stitch_eta_fragment` (single-slot eta ring) and
# `_stitch_multi_eta_fragment` (bridged ansa-metallocene, multi-slot). These
# helpers build the *geometry* of the proper 180-degree in-plane correction;
# the winding-sign math itself always goes through `signed_circulation`
# (`oin/winding.py`) -- never duplicated here.


def _eta_binding_signature(mol: "Chem.Mol", idx: int, binding_set: set) -> tuple:
    """Local substituent signature of one ring (binding) atom.

    (element symbol, sorted tuple of non-ring neighbour symbols). Used to
    detect exocyclic substituents and ring symmetry -- never geometry.
    """
    atom = mol.GetAtomWithIdx(idx)
    subs = tuple(
        sorted(nbr.GetSymbol() for nbr in atom.GetNeighbors() if nbr.GetIdx() not in binding_set)
    )
    return (atom.GetSymbol(), subs)


def _eta_ring_has_exocyclic_substituents(mol: "Chem.Mol | None", binding_idxs: list[int]) -> bool:
    """True if any ring (binding) atom carries a non-H, non-ring neighbour.

    Plain Cp/arene rings (ferrocene) have none; substituted rings (e.g. the
    Ferrocene-halide-face fixture) do. Used to gate the de-aromatized ETKDG
    embedding attempt (Attempt 1b) so unsubstituted rings keep their exact
    pre-Phase-3 analytic-fallback placement (byte-identical, US-005).
    """
    if mol is None:
        return False
    binding_set = set(binding_idxs)
    for idx in binding_idxs:
        if idx >= mol.GetNumAtoms():
            continue
        _, subs = _eta_binding_signature(mol, idx, binding_set)
        if any(s != "H" for s in subs):
            return True
    return False


def _eta_ring_is_symmetric(mol: "Chem.Mol | None", binding_idxs: list[int]) -> bool:
    """True if every binding atom has an identical local substituent signature.

    When true, no OIN winding marker can encode a geometrically observable
    difference (SuperPRD Stereo Phase 3, US-003: "winding is not a geometric
    observable for symmetric rings") -- the correction must be an identity
    no-op regardless of the measured/target windings, e.g. ferrocene's plain
    Cp rings.
    """
    if mol is None:
        return True  # Fail safe toward no-op, never toward an unverified flip.
    binding_set = set(binding_idxs)
    sigs = set()
    for idx in binding_idxs:
        if idx >= mol.GetNumAtoms():
            return True  # Can't verify asymmetry -- be conservative.
        sigs.add(_eta_binding_signature(mol, idx, binding_set))
    return len(sigs) <= 1


def _extract_ring_winding_marker(
    vecs: list, binding_order: list[int]
) -> tuple[Optional[str], Optional[int]]:
    """Return (target_winding, star_atom_idx) for one eta ring's OINVectors.

    Exactly one ``OINVector`` in *vecs* may carry a non-``None`` ``.winding``
    -- the heading/star atom. More than one is a canonical-form violation
    (SuperPRD Stereo Phase 3, "multi-marker same slot") and raises
    ``ValueError`` rather than silently picking a winner. Zero markers is the
    legitimate "zero-marker eta ring" case (legacy/hand-authored OIN) and
    returns ``(None, None)``.
    """
    markers = [v for v in vecs if v.winding is not None]
    if len(markers) > 1:
        raise ValueError(
            "Multi-marker haptic slot: more than one winding marker "
            f"({[m.winding for m in markers]!r}) on ring atoms "
            f"{[m.atom_in_fragment_idx for m in markers]!r} of {binding_order!r}; "
            "canonical OIN form allows exactly one heading atom per haptic ring."
        )
    if not markers:
        return None, None
    star = markers[0]
    return star.winding, star.atom_in_fragment_idx


def _find_slot_index_for_direction(
    template: list, direction: tuple, tol: float = 1e-3
) -> Optional[int]:
    """Return the template slot index whose vector matches *direction*, if any."""
    for i, vec in enumerate(template):
        if all(abs(float(a) - float(b)) < tol for a, b in zip(vec, direction)):
            return i
    return None


def _in_plane_correction_axis(
    binding_pos: np.ndarray,
    centroid: np.ndarray,
    axis_unit: np.ndarray,
) -> "np.ndarray | None":
    """Build the in-plane 180-degree correction rotation axis.

    Axis = centroid -> binding-atom[0], projected into the plane
    perpendicular to *axis_unit* (the metal->centroid outward axis).
    Epsilon-guarded fallback to binding-atom[1] + Gram-Schmidt against
    *axis_unit* when the projection is ill-conditioned (SuperPRD Stereo
    Phase 3, R7). Returns a unit vector, or ``None`` if degenerate even
    after the fallback.
    """

    def _project(v: np.ndarray) -> np.ndarray:
        return v - np.dot(v, axis_unit) * axis_unit

    proj = _project(binding_pos[0] - centroid)
    if float(np.linalg.norm(proj)) < 1e-6 and len(binding_pos) > 1:
        proj = _project(binding_pos[1] - centroid)
    norm = float(np.linalg.norm(proj))
    if norm < 1e-9:
        return None
    return proj / norm


def _proper_180_rotation(rot_axis_unit: np.ndarray):
    """Build a proper (det +1) 180-degree rotation about *rot_axis_unit*.

    Asserts ``det(R) ~= +1`` before returning (SuperPRD Stereo Phase 3, R8 /
    US-006): a reflection would invert pendant substituent chirality and
    must never reach placement.
    """
    from scipy.spatial.transform import Rotation as _Rot  # noqa: PLC0415

    rot = _Rot.from_rotvec(rot_axis_unit * np.pi)
    det = float(np.linalg.det(rot.as_matrix()))
    assert abs(det - 1.0) < 1e-6, (
        f"Haptic-face correction rotation is not proper (det={det}); refusing "
        "to apply what would be a reflection."
    )
    return rot


# ===========================================================================
# Stereo Phase 4 (MiniPRD-B): Zone-A P verify-and-re-embed helpers
# ===========================================================================
#
# Consumes the `_OIN_CIPCode_LP` / `[P@]`/`[P@@]` contract established by
# MiniPRD_ZoneA_P_Encode.md (core/chirality.py). The dummy-metal-copy +
# rdCIPLabeler recipe is REUSED from that module (`_build_dummy_metal_copy`,
# `_lp_cip_label`), never reimplemented (SuperPRD Stereo Phase 4, negative
# constraint). Enforcement lives here -- the adapter's assembled-complex
# stage (`_template_generate`) -- and NOT in `OIN3DGenerator`/engine.py
# (Resolved Q2/RISK-2: `GeneratedStructure.mol` is `Optional` and post-hoc).


def _fresh_fragment_mol(frag_smiles: str) -> "Chem.Mol | None":
    """Parse+sanitize *frag_smiles* independent of any embedding attempt.

    Used only to inspect graph-level properties (chiral tags, CIP labels)
    that do not depend on 3D coordinates -- mirrors the sanitize-with-
    fallback-flags pattern used throughout this module.
    """
    mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            Chem.SanitizeMol(
                mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
            )
        except Exception:
            return None
    return mol


def _zone_a_p_expected_labels(frag_smiles: str) -> list[tuple[int, str]]:
    """Return ``[(local_p_idx, expected_lp_label), ...]`` for *frag_smiles*.

    Graph-based ONLY (no 3D, no embedding) -- computes `rdCIPLabeler`
    directly off the parsed+sanitized fragment mol's EXISTING chiral tag,
    exactly mirroring the last step of ``ChiralityRecoveryUtility.recover()``
    (``core/chirality.py``) that baked this same tag into the OIN string in
    the first place (verify-and-flip, keyed on ``_OIN_CIPCode_LP``). This is
    the "input tag's expected label" the MiniPRD-B verify step checks the
    assembled complex against: it reflects the OIN string's own encoded
    intent, independent of whatever a (possibly mis-embedded) 3D conformer
    later produces. Never derives a tag from 3D perception of a trivalent P
    (perception fails -- SuperPRD spike 3); this is the permitted graph-based
    recompute from an EXISTING tag, same as ``recover()``'s verify step.
    """
    mol = _fresh_fragment_mol(frag_smiles)
    if mol is None:
        return []
    out: list[tuple[int, str]] = []
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        rdCIPLabeler.AssignCIPLabels(mol)
    except Exception:  # noqa: BLE001 - guarded, no expected labels found
        return []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:  # phosphorus only (Sec.3.2: N out of scope)
            continue
        if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        label = atom.GetPropsAsDict().get("_CIPCode")
        if label:
            out.append((atom.GetIdx(), label))
    return out


def _graph_cip_label(mol: "Chem.Mol", p_idx: int) -> "str | None":
    """Graph-based CIP label for P *p_idx* from *mol*'s OWN existing chiral tag.

    No 3D conformer needed -- same technique as ``_zone_a_p_expected_labels``
    (``AssignStereochemistry`` + ``rdCIPLabeler`` off an existing tag, never
    ``AssignAtomChiralTagsFromStructure``), but operates on an in-memory
    ``Mol`` rather than re-parsing a SMILES string. Used as the Task-3
    pre-placement parity guard's reference in ``_stitch_fragment``: it must
    reflect whatever chiral tag the mol CURRENTLY carries -- including a
    TEST-ONLY ``_test_flip_chiral_idx`` seam mutation -- never the original
    OIN-encoded tag, so the guard only ever fires on a genuine ETKDG
    embedding anomaly and never on the test seam itself (that mismatch is
    caught downstream by ``_verify_zone_a_p`` + the re-embed loop, Task 7).
    """
    try:
        copy_mol = Chem.Mol(mol)
        Chem.AssignStereochemistry(copy_mol, cleanIt=True, force=True)
        rdCIPLabeler.AssignCIPLabels(copy_mol)
        return copy_mol.GetAtomWithIdx(p_idx).GetPropsAsDict().get("_CIPCode")
    except Exception:  # noqa: BLE001 - guarded, no reference label available
        return None


def _zone_a_p_tags_in_parsed_oin(parsed_oin: "ParsedOIN") -> list[tuple[int, int]]:
    """Return every tagged Zone-A P in *parsed_oin*, as ``(frag_idx, local_p_idx)``.

    Independent of which generation strategy is used. Lets
    ``MolassemblerAdapter.generate()`` warn on fallback paths (eta
    fallback with ``mol=None``, Molassembler DG fallback) that never run the
    verify-and-re-embed enforcement pinned to `_template_generate`'s
    assembled-complex stage (Task 5, RISK-9 -- "no assembled mol" must never
    be a silent skip).
    """
    out: list[tuple[int, int]] = []
    for frag_idx, frag_smiles in enumerate(parsed_oin.fragments):
        if frag_idx == parsed_oin.metal_fragment_idx:
            continue
        for local_p_idx, _label in _zone_a_p_expected_labels(frag_smiles):
            out.append((frag_idx, local_p_idx))
    return out


def _warn_zone_a_p_fallback(parsed_oin: "ParsedOIN", context: str) -> None:
    """Emit one ``OINStereoWarning`` per Zone-A P tag found in *parsed_oin*.

    Used when no assembled RDKit mol exists to verify/enforce it against
    (Task 5). *context* names the fallback path for the message.
    """
    for frag_idx, local_idx in _zone_a_p_tags_in_parsed_oin(parsed_oin):
        warnings.warn(
            OINStereoWarning(
                f"atom {local_idx} (fragment {frag_idx}): stereo unenforced "
                f"on fallback path ({context}) -- no assembled RDKit mol "
                "available to verify the Zone-A P lone-pair CIP tag; "
                "structure emitted without enforcement for this atom."
            ),
            stacklevel=2,
        )


def _verify_zone_a_p(
    assembled_mol: "Chem.Mol", fragment_inputs: list[tuple[int, str]]
) -> list[int]:
    """Verify each Zone-A P atom in *assembled_mol* against its expected OIN-encoded label.

    *fragment_inputs* is a list of ``(global_atom_idx, expected_lp_label)``
    pairs -- ``global_atom_idx`` indexes into *assembled_mol* (the metal is
    present, bonded DATIVE to each binding atom -- same convention as
    ``xyz2mol.get_tmc_mol()``, which is what ``_build_dummy_metal_copy``
    requires), ``expected_lp_label`` is the lone-pair-convention CIP label
    ('R'/'S') the OIN fragment SMILES encodes for that P atom (see
    ``_zone_a_p_expected_labels``).

    Reuses the MiniPRD-A dummy-metal-copy + `rdCIPLabeler` recipe
    (``core.chirality._build_dummy_metal_copy`` / ``_lp_cip_label``) --
    never reimplemented. Both sides of the comparison are the SAME
    lone-pair convention (never cross-convention -- SuperPRD B1).

    Returns the list of global atom indices whose measured label disagrees
    with the expected one. A dummy-copy construction/label failure (guarded
    internally by the reused helpers, which already warn) is NOT counted as
    a mismatch -- there is nothing actionable to re-embed toward.
    """
    mismatched: list[int] = []
    for global_idx, expected_label in fragment_inputs:
        dummy_mol = _build_dummy_metal_copy(assembled_mol, global_idx)
        if dummy_mol is None:
            continue  # already warned by _build_dummy_metal_copy
        measured = _lp_cip_label(dummy_mol, global_idx)
        if measured is None:
            continue  # non-stereogenic in this copy; nothing to enforce
        if measured != expected_label:
            mismatched.append(global_idx)
    return mismatched


def _zone_a_p_clash_offending_frags(
    normal_frag_meta: list[dict],
    all_pos: list[np.ndarray],
    all_syms: list[str],
    heavy_threshold: float = 1.7,
    h_threshold: float = 1.8,
) -> set[int]:
    """Return ``frag_idx`` of Zone-A-P fragments with an inter-fragment clash.

    Stereo Phase 4 (MiniPRD-C spike finding): a Zone-A-P fragment's
    residual-DOF orientation sweep (``_sweep_rotation_for_clash_avoidance``
    in ``_stitch_fragment``) only sees ``forbidden_positions`` as of ITS OWN
    placement time -- fragments processed LATER in ``parsed_oin.fragments``
    order don't exist yet, so an EARLIER Zone-A-P fragment cannot avoid
    clashing with one placed after it (confirmed empirically: a bulky
    co-ligand placed second folded a methyl H of a Zone-A-P ligand placed
    first to 1.7 A of the metal -- not a literal atom overlap, but close
    enough for xyz2mol to misperceive an M-H bond on round-trip). Checked
    here, AFTER every fragment has been placed, so the SAME re-embed loop
    that already retries on a stereo mismatch can also retry with the
    NOW-COMPLETE ``forbidden_positions`` -- purely geometric, never a
    stereo/mirror operation. Thresholds mirror the existing bidentate
    bite-distortion guard (1.7 A heavy / 1.8 A H) for consistency.
    """
    offending: set[int] = set()
    for meta in normal_frag_meta:
        if not _zone_a_p_expected_labels(meta["frag_smiles"]):
            continue
        start, end = meta["atom_start"], meta["atom_end"]
        binding_global = {start + b for b in meta["binding_idxs"]}
        for gi in range(start, end):
            if gi in binding_global:
                continue  # the binding atom is SUPPOSED to be near the metal
            threshold = h_threshold if all_syms[gi] == "H" else heavy_threshold
            for gj in range(len(all_pos)):
                if start <= gj < end:
                    continue  # intra-fragment (bonded-length) distances, skip
                if float(np.linalg.norm(all_pos[gi] - all_pos[gj])) < threshold:
                    offending.add(meta["frag_idx"])
                    break
            if meta["frag_idx"] in offending:
                break
    return offending


def _zone_a_p_measured_labels_dg(bonded_mol: "Chem.Mol") -> list[str]:
    """Lone-pair-convention CIP labels for every metal-bound P in *bonded_mol*.

    Stereo Phase 4 (MiniPRD-C, Task 5, C2): used to verify Zone-A-P
    enforcement on the DG-produced mol. Unlike ``_verify_zone_a_p`` (which
    matches specific global atom indices against specific expected labels,
    using index bookkeeping only the TEMPLATE path can provide), the DG
    path's ``bonded_mol`` is reconstructed from a CANONICALISED connected
    SMILES (``_build_connected_smiles`` calls ``Chem.MolToSmiles``), so a
    ligand-fragment-local atom index cannot be mapped to a global index here
    without deeper surgery. This returns the SORTED list of measured labels
    instead, for a set-based comparison against the OIN-encoded expected
    labels: correct for every fixture this MiniPRD's own test suite exercises
    (same-label bidentates, e.g. DIPAMP (R,R)) but -- documented honestly --
    it cannot in principle distinguish a coincidental same-set "swap" between
    two DIFFERENTLY-labelled Zone-A-P atoms in the same fragment. Exact
    per-atom index tracking through the DG pipeline is a follow-up item.
    """
    labels: list[str] = []
    for atom in bonded_mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        if not any(nbr.GetAtomicNum() in TRANSITION_METALS_NUM for nbr in atom.GetNeighbors()):
            continue
        dummy_mol = _build_dummy_metal_copy(bonded_mol, atom.GetIdx())
        if dummy_mol is None:
            continue
        label = _lp_cip_label(dummy_mol, atom.GetIdx())
        if label is not None:
            labels.append(label)
    return sorted(labels)


def _assemble_combined_mol(
    metal_mol: "Chem.Mol",
    fragment_mol_parts: list[tuple["Chem.Mol", list[int]]],
    all_pos: list[np.ndarray],
) -> "Chem.Mol | None":
    """Build the single combined RDKit mol (metal + all fragments + one conformer).

    Metal bonds are DATIVE, the conformer comes from *all_pos*. Factored out
    of `_template_generate` (unchanged logic, moved verbatim)
    so the Stereo-Phase-4 enforcement loop (Task 3) can rebuild it after
    re-embedding an offending fragment without duplicating this logic.
    """
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
        return combined_rw.GetMol()
    except Exception:
        return None


def _stitch_multi_eta_fragment(
    frag_smiles: str,
    vectors: list,
    metal_sym: str,
) -> tuple[np.ndarray, list[str], "Chem.Mol | None", list] | None:
    """Place an ansa-metallocene fragment by decomposing into independent rings.

    For bridged fragments (e.g. SiMe2-bridged bis-Cp), the algorithm:
    1. Generates an ETKDG conformer of the complete fragment (both rings + Si + H)
    2. Extracts ring positions directly from the conformer (no SMILES re-parsing)
    3. Applies centroid-plane alignment to transform each ring to its target slot
    4. Reconstructs the bridge atom (Si) and its methyl substituents geometrically

    This avoids SMILES re-parsing failures (e.g. standalone Cp anion kekulization)
    and correctly includes all H atoms in the coordinate transformation.

    Stereo Phase 3 (haptic-face correction): after both rings are placed, a
    single coherent whole-fragment 180-degree correction is applied only when
    both rings' measured circulation disagrees with their target winding in
    the SAME sense; if they disagree (one wants a flip, the other does not),
    placement is left unchanged and the disagreement is reported as a
    `conflict` (never an independent per-ring correction inside a bridged
    fragment -- see the module-level haptic-face-correction helpers).

    Returns (positions, symbols, mol, decisions) or None on failure, where
    `decisions` is a list of one dict per ring (`kind`, `slot`, `status` in
    {"fired", "skipped", "conflict", "no-op"}, `target`, `measured_before`,
    `symmetric`).
    """
    from rdkit.Chem import AllChem as _AllChem  # noqa: PLC0415
    from scipy.spatial.transform import Rotation as _Rotation  # noqa: PLC0415

    print(
        f"[DEBUG] _stitch_multi_eta_fragment called with frag_smiles={frag_smiles}, "
        f"{len(vectors)} vectors",
        file=sys.stderr,
    )

    # ── Phase 1: Setup ─────────────────────────────────────────────────────
    slot_groups: dict[tuple, list[int]] = {}
    # Stereo Phase 3: parallel per-slot OINVector list (preserves SMILES
    # order, like slot_groups) so each ring's own winding marker can be
    # recovered later without routing through winding_by_slot.
    slot_group_vecs: dict[tuple, list] = {}
    for v in vectors:
        key = tuple(round(x, 4) for x in v.vector)
        slot_groups.setdefault(key, []).append(v.atom_in_fragment_idx)
        slot_group_vecs.setdefault(key, []).append(v)

    print(f"[DEBUG] slot_groups: {len(slot_groups)} groups", file=sys.stderr)

    if len(slot_groups) < 2:
        print(
            f"[DEBUG] Not multi-eta (only {len(slot_groups)} slot group(s)), returning None",
            file=sys.stderr,
        )
        return None  # Not multi-eta

    # Parse fragment mol
    mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            Chem.SanitizeMol(
                mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
            )
        except Exception:
            pass

    # ── Phase 2: Find bridge atom and ipso carbons ─────────────────────────
    all_binding: set[int] = set(i for idxs in slot_groups.values() for i in idxs)
    bridge_idx: int | None = None
    ipso_per_slot: dict[tuple, int] = {}

    for atom in mol.GetAtoms():
        if atom.GetIdx() in all_binding:
            continue
        connected_slots: dict[tuple, int] = {}
        for nbr in atom.GetNeighbors():
            for slot_key, bidxs in slot_groups.items():
                if nbr.GetIdx() in bidxs:
                    connected_slots[slot_key] = nbr.GetIdx()
                    break
        if len(connected_slots) >= 2:
            bridge_idx = atom.GetIdx()
            ipso_per_slot = connected_slots
            break

    if bridge_idx is None:
        return None

    # ── Phase 3: BFS to find ring atom sets ────────────────────────────────
    def _bfs_atoms(mol, seeds: set, exclude_idx: int) -> set:
        visited = set(seeds)
        queue = list(seeds)
        while queue:
            cur = queue.pop()
            for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
                ni = nbr.GetIdx()
                if ni != exclude_idx and ni not in visited:
                    visited.add(ni)
                    queue.append(ni)
        return visited

    # ── Phase 4: ETKDG on full fragment mol ────────────────────────────────
    def _embed_fragment(smiles: str) -> "Chem.Mol | None":
        """Generate 3D coordinates for fragment SMILES using ETKDG.

        Handles aromatic rings (Cp, indenyl) that cannot be kekulized by
        converting aromatic bonds to SINGLE and clearing aromatic flags
        before distance-geometry embedding.
        """
        _params = _AllChem.ETKDGv3()
        _params.randomSeed = 42

        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is None:
            return None

        # Compute implicit valences (needed for AddHs)
        try:
            Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            return None

        # Convert aromatic bonds to SINGLE to avoid kekulization failures
        # (Cp rings have 5π electrons, can't be kekulized with standard valences)
        rw = Chem.RWMol(m)
        for bond in rw.GetBonds():
            if bond.GetIsAromatic() or bond.GetBondTypeAsDouble() == 1.5:
                bond.SetBondType(Chem.BondType.SINGLE)
                bond.SetIsAromatic(False)
        for atom in rw.GetAtoms():
            atom.SetIsAromatic(False)
        m = rw.GetMol()

        # Add explicit hydrogens
        try:
            m = Chem.AddHs(m)
        except Exception:
            return None

        # Embed with ETKDG
        if _AllChem.EmbedMolecule(m, _params) == 0:
            return m

        return None

    mol_h = _embed_fragment(frag_smiles)
    etkdg_ok = mol_h is not None

    if etkdg_ok:
        n_atoms_h = mol_h.GetNumAtoms()
        etkdg_pos = np.array(
            [list(mol_h.GetConformer().GetAtomPosition(i)) for i in range(n_atoms_h)],
            dtype=float,
        )
        print(f"[DEBUG] ETKDG embedding succeeded: {n_atoms_h} atoms", file=sys.stderr)
    else:
        print("[DEBUG] ETKDG embedding failed, will use analytic fallback", file=sys.stderr)
        etkdg_pos = None

    def _ring_atoms_with_H(mol_h: Chem.Mol, heavy_idxs: set) -> list[int]:
        """Return heavy_idxs ∪ all H neighbors of those heavy atoms."""
        result = list(heavy_idxs)
        for hi in heavy_idxs:
            for nbr in mol_h.GetAtomWithIdx(hi).GetNeighbors():
                ni = nbr.GetIdx()
                if mol_h.GetAtomWithIdx(ni).GetAtomicNum() == 1:
                    result.append(ni)
        return result

    # ── Phase 5: Transform ring positions to target slot ──────────────────
    placed_results: dict[tuple, tuple] = {}
    # Stereo Phase 3: per-ring haptic-face bookkeeping (target/measured
    # winding, in SMILES order) -- populated below, consumed after Phase 8.
    ring_haptic_info: dict[tuple, dict] = {}

    for slot_idx, (slot_key, bidxs) in enumerate(slot_groups.items()):
        ring_heavy_idxs = _bfs_atoms(mol, set(bidxs), bridge_idx)

        slot_unit = np.array(slot_key, dtype=float)
        slot_norm = float(np.linalg.norm(slot_unit))
        if slot_norm < 1e-9:
            return None
        slot_unit = slot_unit / slot_norm

        if etkdg_ok:
            # Direct ETKDG path
            ring_all_idxs = _ring_atoms_with_H(mol_h, ring_heavy_idxs)
            ring_etkdg_pos = etkdg_pos[ring_all_idxs]
            ring_syms = [mol_h.GetAtomWithIdx(i).GetSymbol() for i in ring_all_idxs]

            binding_etkdg_pos = np.array([etkdg_pos[i] for i in bidxs], dtype=float)
            centroid = binding_etkdg_pos.mean(axis=0)
            ring_radius = float(np.mean([np.linalg.norm(p - centroid) for p in binding_etkdg_pos]))
            print(
                f"[DEBUG] Ring {slot_idx} (ETKDG): {len(ring_all_idxs)} atoms, "
                f"ring_radius={ring_radius:.4f}",
                file=sys.stderr,
            )
        else:
            # Analytic fallback: place binding atoms on a circle with correct H counts
            n_bind = len(bidxs)
            n_H_per_heavy = {i: mol.GetAtomWithIdx(i).GetTotalNumHs() for i in ring_heavy_idxs}
            sorted_heavy = sorted(ring_heavy_idxs)
            cc_bond = 1.40
            circumradius = cc_bond / (2.0 * np.sin(np.pi / n_bind))
            ring_etkdg_pos_list = []
            ring_syms = []
            for k, hidx in enumerate(sorted_heavy):
                ang = 2.0 * np.pi * k / len(sorted_heavy)
                pos = np.array([circumradius * np.cos(ang), circumradius * np.sin(ang), 0.0])
                ring_etkdg_pos_list.append(pos)
                ring_syms.append(mol.GetAtomWithIdx(hidx).GetSymbol())
                n_H = n_H_per_heavy[hidx]
                for _ in range(n_H):
                    h_pos = pos * (1.0 + 1.08 / circumradius)
                    ring_etkdg_pos_list.append(h_pos)
                    ring_syms.append("H")
            ring_etkdg_pos = np.array(ring_etkdg_pos_list, dtype=float)
            ring_all_idxs = list(range(len(ring_etkdg_pos_list)))

            bidx_to_local = {hidx: k for k, hidx in enumerate(sorted_heavy)}
            binding_local_idxs = [bidx_to_local[i] for i in bidxs if i in bidx_to_local]
            binding_etkdg_pos = ring_etkdg_pos[binding_local_idxs]
            centroid = binding_etkdg_pos.mean(axis=0)
            ring_radius = float(np.mean([np.linalg.norm(p - centroid) for p in binding_etkdg_pos]))
            print(
                f"[DEBUG] Ring {slot_idx} (analytic): {len(ring_etkdg_pos)} atoms, "
                f"ring_radius={ring_radius:.4f}",
                file=sys.stderr,
            )

        if ring_radius < 1e-6:
            return None

        # M-C bond length → centroid distance
        binding_sym_str = mol.GetAtomWithIdx(list(bidxs)[0]).GetSymbol()
        m_c_bl = _bond_length(metal_sym, binding_sym_str)
        if m_c_bl > ring_radius:
            centroid_dist = float(np.sqrt(m_c_bl**2 - ring_radius**2))
        else:
            centroid_dist = m_c_bl * 0.80

        target_centroid = slot_unit * centroid_dist

        # Plane normal from SVD on centered binding positions
        centered_bp = binding_etkdg_pos - centroid
        _, _, vh = np.linalg.svd(centered_bp)
        plane_normal = vh[-1]
        if float(np.dot(plane_normal, slot_unit)) < 0:
            plane_normal = -plane_normal

        # Rotate and translate
        try:
            rot, _ = _Rotation.align_vectors([slot_unit], [plane_normal])
        except Exception:
            return None
        ring_pos = rot.apply(ring_etkdg_pos - centroid) + target_centroid

        # Ipso position
        ipso_old = ipso_per_slot[slot_key]
        if etkdg_ok:
            ipso_local = ring_all_idxs.index(ipso_old)
        else:
            sorted_heavy = sorted(ring_heavy_idxs)
            ipso_local = sorted_heavy.index(ipso_old)
        ipso_placed = ring_pos[ipso_local]

        placed_results[slot_key] = (ring_pos, ring_syms, ipso_placed)

        # ── Stereo Phase 3: measure this ring's circulation ─────────────────
        # bidxs is already SMILES/fragment order (a filtered subsequence of
        # `vectors`, which arrives in ascending SMILES-atom-index order --
        # see `_extract_ring_winding_marker` for the same assumption used by
        # `_stitch_eta_fragment`). Map each bidxs atom to its placed position
        # within ring_pos (mirrors the existing `ipso_local` lookup above).
        if etkdg_ok:
            binding_local_positions = [ring_all_idxs.index(b) for b in bidxs]
        else:
            _sorted_heavy_for_bind = sorted(ring_heavy_idxs)
            binding_local_positions = [_sorted_heavy_for_bind.index(b) for b in bidxs]

        ring_binding_pos = ring_pos[binding_local_positions]
        target_winding, star_atom_idx = _extract_ring_winding_marker(
            slot_group_vecs[slot_key], bidxs
        )
        star_local_idx = bidxs.index(star_atom_idx) if star_atom_idx is not None else None
        measured_winding = (
            signed_circulation(ring_binding_pos, star_local_idx, slot_unit)
            if star_local_idx is not None
            else None
        )
        ring_haptic_info[slot_key] = {
            "binding_local_positions": binding_local_positions,
            "axis": slot_unit,
            "target": target_winding,
            "measured": measured_winding,
            "symmetric": _eta_ring_is_symmetric(mol, bidxs),
        }

    # ── Phase 6: Place Si bridge atom ──────────────────────────────────────
    def _place_bridge_atom(p1: np.ndarray, p2: np.ndarray, bond_length: float = 1.87) -> np.ndarray:
        """Find Si position given two ipso C positions and bond length constraints."""
        mid = (p1 + p2) / 2.0
        half_d = float(np.linalg.norm(p2 - p1) / 2.0)
        h_sq = bond_length**2 - half_d**2
        if h_sq < 0:
            h = 0.0
        else:
            h = float(np.sqrt(h_sq))

        seg = (p2 - p1) / float(np.linalg.norm(p2 - p1))
        mid_perp = mid - float(np.dot(mid, seg)) * seg
        norm_mid_perp = float(np.linalg.norm(mid_perp))

        if norm_mid_perp < 1e-9:
            arb = (
                np.array([1.0, 0.0, 0.0], dtype=float)
                if abs(seg[0]) < 0.9
                else np.array([0.0, 1.0, 0.0], dtype=float)
            )
            direction = arb - float(np.dot(arb, seg)) * seg
            direction = direction / float(np.linalg.norm(direction))
        else:
            direction = mid_perp / norm_mid_perp

        return mid + h * direction

    slot_keys = list(slot_groups.keys())
    ipso1 = placed_results[slot_keys[0]][2]  # ipso_placed from phase 5
    ipso2 = placed_results[slot_keys[1]][2]  # ipso_placed from phase 5

    print(f"[DEBUG] ipso1={ipso1}, ipso2={ipso2}", file=sys.stderr)

    si_pos = _place_bridge_atom(ipso1, ipso2)
    print(f"[DEBUG] si_pos={si_pos}", file=sys.stderr)
    print(
        f"[DEBUG] Si-ipso1 dist={np.linalg.norm(si_pos - ipso1):.4f}, "
        f"Si-ipso2 dist={np.linalg.norm(si_pos - ipso2):.4f}",
        file=sys.stderr,
    )

    # ── Phase 7: Place Si substituents (methyl groups) ──────────────────────
    def _place_tetrahedral_methyls(
        si_pos: np.ndarray,
        ipso1: np.ndarray,
        ipso2: np.ndarray,
        si_c_bond: float = 1.87,
        c_h_bond: float = 1.09,
    ) -> list[tuple[str, np.ndarray]]:
        """Place methyl C and H atoms tetrahedral around Si."""
        v1 = (ipso1 - si_pos) / float(np.linalg.norm(ipso1 - si_pos))
        v2 = (ipso2 - si_pos) / float(np.linalg.norm(ipso2 - si_pos))

        sum12 = v1 + v2
        cross12 = np.cross(v1, v2)
        cross_norm = float(np.linalg.norm(cross12))

        if cross_norm < 1e-9:
            arb = (
                np.array([1.0, 0.0, 0.0], dtype=float)
                if abs(v1[0]) < 0.9
                else np.array([0.0, 1.0, 0.0], dtype=float)
            )
            perp = arb - float(np.dot(arb, v1)) * v1
            perp = perp / float(np.linalg.norm(perp))
        else:
            perp = cross12 / cross_norm

        mid_v34 = -sum12 / 2.0
        mid_norm = float(np.linalg.norm(mid_v34))

        if mid_norm < 1e-9:
            t = 1.0 / np.sqrt(2.0)
        else:
            t_sq = 1.0 - mid_norm**2
            t = float(np.sqrt(max(0.0, t_sq)))

        # BUG FIX: Remove the * 2.0 scaling that breaks tetrahedral angles
        v3 = mid_v34 + t * perp
        v4 = mid_v34 - t * perp

        v3_norm = float(np.linalg.norm(v3))
        v4_norm = float(np.linalg.norm(v4))

        if v3_norm > 1e-9:
            v3 = v3 / v3_norm
        else:
            v3 = np.array([0.0, 0.0, 1.0], dtype=float)

        if v4_norm > 1e-9:
            v4 = v4 / v4_norm
        else:
            v4 = np.array([0.0, 0.0, -1.0], dtype=float)

        me1_pos = si_pos + si_c_bond * v3
        me2_pos = si_pos + si_c_bond * v4

        results: list[tuple[str, np.ndarray]] = []
        tet_angle = float(np.arccos(-1.0 / 3.0))

        for me_pos, v_me in [(me1_pos, v3), (me2_pos, v4)]:
            results.append(("C", me_pos))
            v_to_si = -v_me
            arb2 = (
                np.array([1.0, 0.0, 0.0], dtype=float)
                if abs(v_to_si[0]) < 0.9
                else np.array([0.0, 1.0, 0.0], dtype=float)
            )
            u1 = arb2 - float(np.dot(arb2, v_to_si)) * v_to_si
            u1 = u1 / float(np.linalg.norm(u1))
            u2 = np.cross(v_to_si, u1)

            for k in range(3):
                phi = 2.0 * np.pi * k / 3.0
                # BUG FIX: H atoms should be at tet_angle from v_to_si (away from Si),
                # not from v_me. Change: cos(tet_angle)*v_me → -cos(tet_angle)*v_me
                h_dir = -float(np.cos(tet_angle)) * v_me + float(np.sin(tet_angle)) * (
                    float(np.cos(phi)) * u1 + float(np.sin(phi)) * u2
                )
                h_dir = h_dir / float(np.linalg.norm(h_dir))
                results.append(("H", me_pos + c_h_bond * h_dir))

        return results

    me_atoms = _place_tetrahedral_methyls(si_pos, ipso1, ipso2)

    # Debug: print methyl positions
    for i, (sym, pos) in enumerate(me_atoms):
        if sym == "C":
            dist_to_si = np.linalg.norm(pos - si_pos)
            print(f"[DEBUG] Methyl C{i}: pos={pos}, dist_to_si={dist_to_si:.4f}", file=sys.stderr)

    # ── Phase 8: Assemble and return ───────────────────────────────────────
    all_positions: list[np.ndarray] = []
    all_symbols: list[str] = []

    # Assemble positions from all rings
    for i, slot_key in enumerate(slot_keys):
        ring_pos, ring_syms, ipso_placed = placed_results[slot_key]
        all_positions.extend(ring_pos)
        all_symbols.extend(ring_syms)

    # Add Si
    si_atom_idx = len(all_positions)
    all_positions.append(si_pos)
    all_symbols.append("Si")

    # Add methyls (2 C atoms + 6 H atoms)
    for sym, pos in me_atoms:
        all_positions.append(pos)
        all_symbols.append(sym)

    # Heavy-atom tracking: map input fragment heavy atom indices to output global indices
    heavy_atom_map: dict[int, int] = {}
    output_idx = 0
    for slot_idx, slot_key in enumerate(slot_keys):
        ring_pos, ring_syms, _ = placed_results[slot_key]
        ring_heavy_idxs = _bfs_atoms(mol, set(slot_groups[slot_key]), bridge_idx)
        if etkdg_ok:
            ring_all_idxs = _ring_atoms_with_H(mol_h, ring_heavy_idxs)
            for local_i, global_heavy_idx in enumerate(ring_all_idxs):
                if mol_h.GetAtomWithIdx(global_heavy_idx).GetAtomicNum() != 1:
                    heavy_atom_map[global_heavy_idx] = output_idx + local_i
        else:
            sorted_heavy = sorted(ring_heavy_idxs)
            for k, hidx in enumerate(sorted_heavy):
                heavy_atom_map[hidx] = output_idx + k
        output_idx += len(ring_syms)

    heavy_atom_map[bridge_idx] = si_atom_idx
    print(f"[DEBUG] Heavy atom map: {heavy_atom_map}", file=sys.stderr)

    # ── Stereo Phase 3: coherent whole-fragment correction (or conflict) ────
    # Never an independent per-ring correction inside a bridged fragment
    # (SuperPRD Stereo Phase 3, negative constraint): apply ONE proper
    # rotation to the whole assembled fragment only when both rings agree
    # they need a flip; if they disagree, leave placement unchanged.
    decisions: list[dict] = []
    key0, key1 = slot_keys[0], slot_keys[1]
    info0, info1 = ring_haptic_info[key0], ring_haptic_info[key1]

    def _has_preference(info: dict) -> bool:
        return info["target"] is not None and not info["symmetric"]

    def _wants_flip(info: dict) -> bool:
        return _has_preference(info) and info["measured"] != info["target"]

    if not _has_preference(info0) and not _has_preference(info1):
        status = "no-op"
    else:
        want0, want1 = _wants_flip(info0), _wants_flip(info1)
        if want0 != want1:
            status = "conflict"
        elif want0 and want1:
            status = "fired"
        else:
            status = "skipped"

    if status == "fired":
        axis0, axis1 = info0["axis"], info1["axis"]
        cross_ax = np.cross(axis0, axis1)
        cross_norm = float(np.linalg.norm(cross_ax))
        if cross_norm > 1e-6:
            # cross(axis0, axis1) is perpendicular to BOTH rings' own
            # metal->centroid axes, so a single 180 deg rotation about it
            # inverts each ring's axis (axis_i -> -axis_i) simultaneously --
            # exactly the per-ring face-flip operation, applied coherently.
            rot_axis = cross_ax / cross_norm
        else:
            # axis0 / axis1 are (anti)parallel (the common near-linear
            # sandwich case): any axis perpendicular to axis0 is valid;
            # reuse ring0's own in-plane construction (R7 fallback).
            ring0_pos = placed_results[key0][0]
            b0 = ring0_pos[info0["binding_local_positions"]]
            cen0 = b0.mean(axis=0)
            rot_axis = _in_plane_correction_axis(b0, cen0, axis0)
            if rot_axis is None:
                rot_axis = np.array([1.0, 0.0, 0.0])

        rot = _proper_180_rotation(rot_axis)
        pivot = np.mean(all_positions, axis=0)
        all_positions = [rot.apply(p - pivot) + pivot for p in all_positions]

    for key, info in ((key0, info0), (key1, info1)):
        decisions.append(
            {
                "kind": "multi-eta",
                "slot": key,
                "status": status,
                "target": info["target"],
                "measured_before": info["measured"],
                "symmetric": info["symmetric"],
            }
        )

    # Return XYZ positions and symbols; mol=None for now (XYZ is the deliverable)
    return np.array(all_positions, dtype=float), all_symbols, None, decisions


def _stitch_eta_fragment(
    frag_smiles: str,
    binding_idxs: list[int],
    slot_unit: np.ndarray,
    metal_sym: str,
    winding: Optional[str] = None,
    star_atom_idx: Optional[int] = None,
) -> tuple[np.ndarray, list[str], "Chem.Mol | None", dict] | None:
    """Place an eta-type ligand (Cp, arene) by centroid-plane alignment.

    For eta-n ligands where *binding_idxs* all share the same slot direction:
    1. Tries ETKDGv3 embedding of the organic fragment first.
    2. Falls back to a de-aromatized ETKDGv3 embedding when the ring carries
       exocyclic substituents that block kekulization (Attempt 1b).
    3. Falls back to analytic regular-ring geometry if both fail (e.g. plain
       Cp anion `[cH]1[cH][cH][cH][cH]1`, which RDKit cannot kekulize and
       which carries nothing worth preserving beyond the ring itself).
    4. Computes centroid of the binding atoms and estimates M–centroid distance
       from the M–C bond length and the ring circumradius.
    5. Rotates the fragment so its binding-atom plane normal aligns with
       *slot_unit* (pointing away from the metal).
    6. Translates the centroid to ``slot_unit * centroid_dist``.
    7. (Stereo Phase 3) Measures the placed ring's circulation and, if it
       disagrees with *winding*, applies a proper 180-degree in-plane
       correction (see module-level haptic-face-correction helpers).

    Parameters
    ----------
    winding:
        Target OIN winding character (``'>'``/``'<'``) for this ring, i.e.
        the single non-``None`` ``OINVector.winding`` in the ring's slot
        group. ``None`` means no marker was present (legacy/zero-marker
        ring) -- correction is skipped and recorded as a ``no-op``.
    star_atom_idx:
        The *original* fragment atom index (as found in *binding_idxs*) of
        the heading/star atom that carries *winding*. Required whenever
        *winding* is not ``None``.

    Returns:
    -------
    (positions, symbols, mol, decision) or None on failure.
    mol is the RDKit Mol with bond connectivity, or None for analytic-geometry
    fallback. decision is a dict describing the haptic-face-correction
    outcome for this ring (`status` in {"fired", "skipped", "conflict" (never
    for a single ring), "no-op"}).
    """
    from rdkit.Chem import AllChem  # noqa: PLC0415

    n_binding = len(binding_idxs)
    if n_binding < 2:
        return None

    positions: np.ndarray | None = None
    symbols: list[str] = []
    valid_idxs: list[int] = []
    etkdg_mol: Chem.Mol | None = None
    smiles_mol: Chem.Mol | None = None
    reordered_for_smiles_mol = False

    # Fresh, untouched parse used only for substituent/symmetry inspection --
    # independent of whatever mutations the embedding attempts below apply.
    sig_mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
    has_substituents = _eta_ring_has_exocyclic_substituents(sig_mol, binding_idxs)

    # ── Attempt 1: ETKDGv3 (plain aromatic-preserving sanitize) ───────────────
    mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
    if mol is not None:
        try:
            Chem.SanitizeMol(mol)
            Chem.SetAromaticity(mol)
        except Exception:
            try:
                Chem.SanitizeMol(
                    mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                )
                Chem.SetAromaticity(mol)
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

    # ── Attempt 1b: anionic-ring ETKDGv3 (substituted rings only) ────────────
    # Attempt 1 fails to kekulize many substituted Cp/arene rings: a 5- (or
    # other odd-) membered aromatic ring where every ring atom already has an
    # explicit non-H substituent has no valid neutral alternating-bond Kekule
    # structure -- it needs exactly one ring atom to carry the anionic
    # (cyclopentadienide-like) formal charge that a plain `[cH]` ring gets
    # for free from RDKit's implicit-H/aromaticity model. Which position
    # works is parity-dependent, so try each ring atom in turn and keep the
    # first that sanitizes. This keeps bonds AROMATIC (unlike a de-aromatized
    # single-bond embed), so ETKDG's force field gives genuinely planar,
    # aromatic-length (~1.40 A) geometry natively -- required for the
    # independent XYZ->OIN bond-order perception (xyz2mol) to recognize the
    # ring at all. For a PLAIN ring (ferrocene, TiCp2Me2, …) this attempt is
    # never reached (`has_substituents` is False) -- the analytic fallback
    # (Attempt 2) reproduces the exact pre-Phase-3 placement unchanged
    # (byte-identical, US-005).
    if (positions is None or len(valid_idxs) < 2) and has_substituents:
        mol_b = None
        for charge_idx in binding_idxs:
            candidate = Chem.MolFromSmiles(frag_smiles, sanitize=False)
            if candidate is None or charge_idx >= candidate.GetNumAtoms():
                continue
            candidate.GetAtomWithIdx(charge_idx).SetFormalCharge(-1)
            try:
                Chem.SanitizeMol(candidate)
            except Exception:
                continue
            mol_b = candidate
            break

        if mol_b is not None:
            try:
                mol_b = Chem.AddHs(mol_b)
                _params_b = AllChem.ETKDGv3()
                _params_b.randomSeed = 42
                r_b = AllChem.EmbedMolecule(mol_b, _params_b)
            except Exception:
                r_b = -1
            if r_b == 0:
                n_atoms = mol_b.GetNumAtoms()
                positions = np.array(
                    [list(mol_b.GetConformer().GetAtomPosition(i)) for i in range(n_atoms)],
                    dtype=float,
                )
                symbols = [mol_b.GetAtomWithIdx(i).GetSymbol() for i in range(n_atoms)]
                valid_idxs = [i for i in binding_idxs if i < n_atoms]
                etkdg_mol = mol_b
                smiles_mol = mol_b

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
    from scipy.spatial.transform import Rotation  # noqa: PLC0415

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
        symbols = [
            smiles_mol.GetAtomWithIdx(i).GetSymbol() for i in range(smiles_mol.GetNumAtoms())
        ]
        etkdg_mol = smiles_mol
        reordered_for_smiles_mol = True

    # ── Stereo Phase 3: haptic-face correction ───────────────────────────────
    # Map each SMILES-order binding atom (binding_idxs order) to its index in
    # the FINAL `positions` array, tracking the heavy-first reorder above.
    if reordered_for_smiles_mol:
        final_binding_local_idxs = list(range(n_binding))
    else:
        final_binding_local_idxs = valid_idxs

    decision: dict = {
        "kind": "eta",
        "status": "no-op",
        "target": winding,
        "measured_before": None,
        "measured_after": None,
        "symmetric": _eta_ring_is_symmetric(sig_mol, binding_idxs),
    }

    if winding is not None and len(final_binding_local_idxs) == n_binding:
        if star_atom_idx is None:
            raise ValueError("_stitch_eta_fragment: winding target given without star_atom_idx")
        try:
            star_local_idx = binding_idxs.index(star_atom_idx)
        except ValueError as exc:
            raise ValueError(
                f"star_atom_idx={star_atom_idx} is not one of this ring's binding_idxs "
                f"{binding_idxs!r}"
            ) from exc

        binding_final = positions[final_binding_local_idxs]
        centroid_final = binding_final.mean(axis=0)
        # Metal is always at the origin in this local frame (see caller);
        # by construction the ring centroid lies exactly along slot_unit.
        measured = signed_circulation(binding_final, star_local_idx, slot_unit)
        decision["measured_before"] = measured

        if decision["symmetric"]:
            decision["status"] = "no-op"
        elif measured == winding:
            decision["status"] = "skipped"
            decision["measured_after"] = measured
        else:
            rot_axis = _in_plane_correction_axis(binding_final, centroid_final, slot_unit)
            if rot_axis is None:
                decision["status"] = "no-op"
            else:
                correction = _proper_180_rotation(rot_axis)
                positions = correction.apply(positions - centroid_final) + centroid_final
                binding_final2 = positions[final_binding_local_idxs]
                centroid_final2 = binding_final2.mean(axis=0)
                axis2 = centroid_final2
                axis2_norm = float(np.linalg.norm(axis2))
                axis2_unit = axis2 / axis2_norm if axis2_norm > 1e-9 else slot_unit
                decision["measured_after"] = signed_circulation(
                    binding_final2, star_local_idx, axis2_unit
                )
                decision["status"] = "fired"

    return positions, symbols, etkdg_mol, decision


def _strip_dummy_atoms(
    positions: np.ndarray,
    symbols: list[str],
    mol: "Chem.Mol",
    dummy_indices: list[int],
) -> tuple[np.ndarray, list[str], "Chem.Mol"]:
    """Delete *dummy_indices* from *positions*/*symbols*/*mol*, re-deriving indices.

    Stereo Phase 4 (MiniPRD-C, Task 4/C4.1): every Z=0 dummy attached by
    ``_attach_dummy_metal`` must be gone before this fragment leaves
    ``_stitch_fragment`` -- it must never reach ``combined_mol`` or the
    written XYZ block. Dummies were appended as the highest original-heavy-
    atom indices (before ``AddHs``), so ``binding_idxs`` (all < those
    indices) are unaffected by their removal; only trailing H-atom indices
    shift down, which downstream code never assumes are stable.
    """
    if not dummy_indices:
        return positions, symbols, mol
    keep_mask = np.ones(len(symbols), dtype=bool)
    for di in dummy_indices:
        keep_mask[di] = False
    new_positions = positions[keep_mask]
    new_symbols = [s for i, s in enumerate(symbols) if keep_mask[i]]
    rw = Chem.RWMol(mol)
    for di in sorted(dummy_indices, reverse=True):
        rw.RemoveAtom(di)
    return new_positions, new_symbols, rw.GetMol()


def _sweep_rotation_for_clash_avoidance(
    positions: np.ndarray,
    binding_pos: np.ndarray,
    axis: np.ndarray,
    other_idxs: list[int],
    forbidden_positions: "list | None",
    require_existing_clash: bool = False,
    clash_threshold: float = 1.5,
) -> np.ndarray:
    """Sweep rotation about *axis* (through *binding_pos*) to maximise min distance to forbidden.

    Resolves the one residual rotational DOF a single-vector alignment
    leaves undetermined: aligning one vector (e.g. P->dummy onto a slot
    direction, or the pre-existing "outward" vector onto a slot direction)
    fixes only 2 of 3 rotational DOF -- the remaining rotation about that
    axis is arbitrary and can swing the rest of a bulky monodentate fragment
    into an already-placed neighbour (confirmed empirically during Stereo
    Phase 4 MiniPRD-C: an unconstrained residual rotation folded a methyl
    group to a sub-Angstrom clash with a neighbouring ligand).

    When *require_existing_clash* is True, only searches -- and only ever
    changes *positions* -- if the CURRENT orientation already clashes
    (below *clash_threshold*). This is a monotonic guard: it can only turn an
    already-bad placement into a better one, never move an already-clash-free
    placement, so every existing tag-free golden fixture's byte-identical
    output is unaffected.
    """
    from scipy.spatial.transform import Rotation as _Rotation  # noqa: PLC0415

    if not forbidden_positions or not other_idxs:
        return positions
    forb_np = np.array(forbidden_positions, dtype=float)

    def _min_dist(pos_arr: np.ndarray) -> float:
        other = pos_arr[other_idxs]
        dists = np.sqrt(((other[:, None, :] - forb_np[None, :, :]) ** 2).sum(axis=-1))
        return float(dists.min())

    current_dist = _min_dist(positions)
    if require_existing_clash and current_dist >= clash_threshold:
        return positions

    axis_norm = float(np.linalg.norm(axis))
    if axis_norm < 1e-9:
        return positions
    axis_unit = axis / axis_norm
    centered = positions - binding_pos
    best_angle = 0.0
    best_dist = current_dist
    for deg in range(5, 360, 5):
        rot_try = _Rotation.from_rotvec(axis_unit * np.radians(deg))
        pos_try = rot_try.apply(centered) + binding_pos
        dist = _min_dist(pos_try)
        if dist > best_dist:
            best_dist = dist
            best_angle = float(np.radians(deg))
    if best_angle == 0.0:
        return positions
    rot_best = _Rotation.from_rotvec(axis_unit * best_angle)
    return rot_best.apply(centered) + binding_pos


def _stitch_fragment(
    frag_smiles: str,
    binding_idxs: list[int],
    target_positions: list[np.ndarray],
    slot_units: list[np.ndarray] | None = None,
    forbidden_positions: list[np.ndarray] | None = None,
    seed: int = 42,
    _test_flip_chiral_idx: int | None = None,
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
    seed:
        ETKDG ``randomSeed``. Stereo Phase 4 (MiniPRD-B) re-embeds the SAME
        fragment with a NEW seed here when the Zone-A P verify step
        (`_verify_zone_a_p`) detects a mismatch -- never a mirror/improper
        transform (SuperPRD B2/B3).
    _test_flip_chiral_idx:
        TEST-ONLY injection seam (Stereo Phase 4 MiniPRD-B, Task 7). When
        set, the chiral tag of this ONE fragment-local atom is flipped
        in-place BEFORE ETKDG embeds, so the resulting conformer is
        genuinely mis-embedded at exactly that atom (a realistic simulation
        of an ETKDG error) while every other stereocenter in the fragment
        (e.g. a co-resident Zone-A P in a bidentate fragment) embeds
        normally. Never used by production code paths -- exists solely so
        the bounded re-embed enforcement loop in `_template_generate` has a
        deterministic way to be exercised by tests.

    Returns:
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
            Chem.SanitizeMol(
                mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
            )
        except Exception:
            for atom in mol.GetAtoms():
                pass  # valence already calculated during partial sanitization

    # Stereo Phase 4 (MiniPRD-B, Task 7): TEST-ONLY forced mis-embed. Flips
    # ONE atom's chiral tag before ETKDG embeds, so this fragment's own
    # `randomSeed` picks the correct handedness for every OTHER stereocenter
    # but the wrong one for this atom -- unlike a whole-fragment mirror, this
    # never touches a co-resident stereocenter's configuration.
    if _test_flip_chiral_idx is not None and _test_flip_chiral_idx < mol.GetNumAtoms():
        _flip_atom = mol.GetAtomWithIdx(_test_flip_chiral_idx)
        _flip_tag = _flip_atom.GetChiralTag()
        if _flip_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            _flip_atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif _flip_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            _flip_atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

    # Set NoImplicit on binding atoms to avoid spurious H addition.
    # When the binding atom forms an M-L bond in the complex its valence is
    # already fully used (e.g. C in C#N: triple bond + 1 metal = 4 = max).
    # Use +2 (not +1) to account for double bonds to the metal (e.g. V=O:
    # O has _bsum=0, _dv=2 → 0+2=2 >= 2 → NoImplicit, preventing H2O gen).
    rw = Chem.RWMol(mol)

    # Stereo Phase 4 (MiniPRD-C): Zone-A P dummy-metal embed. Gate on ANY
    # denticity (Task 1, C2) -- never co-conditioned on len(binding_idxs)
    # == 1. `_zone_a_p_expected_labels` is the same graph-based recompute
    # `_verify_zone_a_p` trusts, so local_p_idx here matches this fragment's
    # own atom ordering with no re-derivation needed.
    zone_a_p_targets = _zone_a_p_expected_labels(frag_smiles)
    dummy_by_p_idx: dict[int, int] = {}
    pre_embed_expected: dict[int, "str | None"] = {}
    if zone_a_p_targets:
        for p_idx, _expected_label in zone_a_p_targets:
            if p_idx >= rw.GetNumAtoms():
                continue
            # Parity-guard reference (Task 3): the CIP implied by THIS mol's
            # own CURRENT chiral tag -- test-seam-aware (see
            # _graph_cip_label), NOT the original OIN-encoded tag. Computed
            # before the dummy attach so it reflects only a real chiral-tag
            # mutation (e.g. _test_flip_chiral_idx), never an artifact of
            # adding the dummy substituent itself.
            pre_embed_expected[p_idx] = _graph_cip_label(rw.GetMol(), p_idx)
            # Pinned order (C4.2): SetNoImplicit(P) -> attach dummy, both
            # before AddHs/embed.
            rw.GetAtomWithIdx(p_idx).SetNoImplicit(True)
            dummy_by_p_idx[p_idx] = _attach_dummy_metal(rw, p_idx)

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
        # Fall back to rw's current state (preserves the dummy attach above
        # and any partial NoImplicit edits) -- never the pre-edit `mol`,
        # which would silently drop the just-attached dummy.
        mol = Chem.AddHs(rw.GetMol())

    if dummy_by_p_idx:
        for _p_idx, _dummy_idx in dummy_by_p_idx.items():
            assert mol.GetAtomWithIdx(_dummy_idx).GetAtomicNum() == 0, (
                f"dummy atom at idx {_dummy_idx} lost its Z=0 identity after AddHs"
            )

    _etkdg_params = AllChem.ETKDGv3()
    _etkdg_params.randomSeed = seed
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
    positions = np.array([list(conf.GetAtomPosition(i)) for i in range(n_atoms)], dtype=float)
    symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_atoms)]

    # Validate binding indices (after AddHs, n_atoms may be larger than the
    # original heavy-atom count but binding_idxs were computed from heavy atoms
    # so they are still valid as long as < original heavy-atom count).
    for bidx in binding_idxs:
        if bidx >= n_atoms:
            return None

    # Stereo Phase 4 (MiniPRD-C, Task 3/C4.2): pre-placement parity guard.
    # A rigid rotation (the orientation step below) can never change a CIP
    # label, so this only checks that ETKDG itself embedded the tetrahedron
    # with the handedness the OIN string encodes. A mismatch means ETKDG
    # embedded the wrong-handed conformer -- a hard failure, never silently
    # masked by the (unrelated) re-embed loop upstream.
    if zone_a_p_targets:
        for p_idx, _expected_label in zone_a_p_targets:
            if p_idx not in dummy_by_p_idx:
                continue
            reference = pre_embed_expected.get(p_idx)
            if reference is None:
                continue  # nothing meaningful to check against
            measured = _lp_cip_label(mol, p_idx)
            assert measured == reference, (
                f"Zone-A P atom {p_idx}: pre-placement lone-pair CIP parity "
                f"mismatch (this mol's own chiral tag implies {reference!r}, "
                f"embedded conformer gives {measured!r}) -- ETKDG embedded "
                "the wrong-handed tetrahedron; hard failure, never a silent "
                "loop-mask."
            )

    if len(binding_idxs) == 1:
        # Monodentate: rigid translation so binding atom lands on target.
        t = target_positions[0] - positions[binding_idxs[0]]
        positions = positions + t

        bidx = binding_idxs[0]
        _zone_a_dummy_idx = dummy_by_p_idx.get(bidx)
        if _zone_a_dummy_idx is not None and slot_units and len(slot_units) >= 1:
            # US-C2 (monodentate): replace the "outward" centroid heuristic
            # below with the exact embedded P->dummy (metal-facing) vector
            # for Zone-A-P fragments only -- the dummy IS the faithful
            # metal-facing reference the tag was embedded against.
            slot_u = np.array(slot_units[0], dtype=float)
            binding_pos = positions[bidx]
            dummy_vec = positions[_zone_a_dummy_idx] - binding_pos
            dummy_norm = float(np.linalg.norm(dummy_vec))
            if dummy_norm > 1e-9:
                dummy_unit = dummy_vec / dummy_norm
                metal_dir = -slot_u  # P->metal direction (metal sits at the origin)
                if float(np.dot(dummy_unit, metal_dir)) < 0.999999:
                    try:
                        rot_o, _ = Rotation.align_vectors([metal_dir], [dummy_unit])
                        centered = positions - binding_pos
                        positions = rot_o.apply(centered) + binding_pos
                    except Exception:
                        pass  # Keep current orientation on failure

                # Aligning a single vector (P->dummy onto metal_dir) fixes
                # only 2 of 3 rotational DOF -- the remaining rotation about
                # that (now-fixed) axis is arbitrary and can swing the rest
                # of the ligand into an already-placed fragment (confirmed
                # empirically: an unconstrained residual rotation folded a
                # methyl group to a 0.27 A H...H clash with a neighbouring
                # ligand). Sweep that residual angle to maximise distance
                # from forbidden_positions -- never touches the fixed
                # metal-facing vector itself, so the CIP-determining
                # handedness is unaffected by this search. Unconditional
                # (not gated on an existing clash): this is new Zone-A-P-only
                # code, no pre-existing golden output to preserve.
                other_idxs_for_sweep = [i for i in range(n_atoms) if i != bidx]
                positions = _sweep_rotation_for_clash_avoidance(
                    positions,
                    binding_pos,
                    metal_dir,
                    other_idxs_for_sweep,
                    forbidden_positions,
                    require_existing_clash=False,
                )
        # Orientation fix: rotate fragment around the binding atom so that
        # non-binding atoms face away from the metal (in the +slot direction).
        # Without this, H atoms from NH3 / [CH] groups may end up between the
        # metal and the binding atom, causing xyz2mol to form spurious M-H bonds.
        elif slot_units and len(slot_units) >= 1 and n_atoms > 1:
            slot_u = np.array(slot_units[0], dtype=float)
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
                            rot_o, _ = Rotation.align_vectors([slot_u], [outward_unit])
                            centered = positions - binding_pos
                            positions = rot_o.apply(centered) + binding_pos
                        except Exception:
                            pass  # Keep current orientation on failure

        positions, symbols, mol = _strip_dummy_atoms(
            positions, symbols, mol, list(dummy_by_p_idx.values())
        )
        assert not any(
            mol.GetAtomWithIdx(i).GetAtomicNum() == 0 for i in range(mol.GetNumAtoms())
        ), "dummy atom leaked past _stitch_fragment (monodentate path)"
        assert len(positions) == mol.GetNumAtoms() == len(symbols), (
            "positions/mol/symbols length mismatch after dummy strip (monodentate path)"
        )
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
        _degenerate_positions, _degenerate_symbols, _degenerate_mol = _strip_dummy_atoms(
            positions + (t_center - c_center), symbols, mol, list(dummy_by_p_idx.values())
        )
        assert not any(
            _degenerate_mol.GetAtomWithIdx(i).GetAtomicNum() == 0
            for i in range(_degenerate_mol.GetNumAtoms())
        ), "dummy atom leaked past _stitch_fragment (bidentate degenerate path)"
        return _degenerate_positions, _degenerate_symbols, _degenerate_mol

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
        _forb_np = np.array(forbidden_positions, dtype=float) if forbidden_positions else None

        def _has_clash(pos_arr: np.ndarray) -> bool:
            """Return True if a non-binding atom clashes with a placed atom.

            A clash is any non-binding atom within 1.5 Å of an
            already-placed atom (forbidden_positions).
            """
            if _forb_np is None or len(_nb_idxs) == 0:
                return False
            nb = pos_arr[_nb_idxs]
            dists = np.sqrt(((nb[:, None, :] - _forb_np[None, :, :]) ** 2).sum(axis=-1))
            return bool(dists.min() < 1.5)

        centered = positions_aligned - t_center

        # Stereo Phase 4 (MiniPRD-C, Q-C1): for Zone-A-P bidentate fragments
        # this bite-axis rotation is the ONE residual DOF scipy's 2-vector
        # align_vectors leaves undetermined -- and it is EXACTLY the DOF
        # that decides which face each P's dummy (metal-facing reference)
        # ends up on. Resolve it by maximising a least-squares-style
        # alignment score of the embedded P->dummy vector(s) onto the
        # metal-facing direction(s) (-slot_units[binding_rank]), overriding
        # the clash-avoidance objective used below for non-Zone-A-P
        # fragments -- stereochemical correctness takes priority here. The
        # existing post-search clash/distance rejection further below still
        # applies unchanged and routes an incompatible-bite result to DG
        # (Task 5) if this angle distorts the geometry too much.
        _p_dummy_ranks = (
            [
                (rank, dummy_by_p_idx[bidx])
                for rank, bidx in enumerate(binding_idxs)
                if bidx in dummy_by_p_idx
            ]
            if slot_units
            else []
        )

        if _p_dummy_ranks:

            def _dummy_alignment_score(pos_arr: np.ndarray) -> float:
                score = 0.0
                for rank, dummy_idx in _p_dummy_ranks:
                    p_bidx = binding_idxs[rank]
                    vec = pos_arr[dummy_idx] - pos_arr[p_bidx]
                    vnorm = float(np.linalg.norm(vec))
                    if vnorm < 1e-9:
                        continue
                    target_dir = -np.array(slot_units[rank], dtype=float)
                    score += float(np.dot(vec / vnorm, target_dir))
                return score

            best_angle = 0.0
            best_score = _dummy_alignment_score(positions_aligned)
            for deg in range(5, 360, 5):
                rot_try = Rotation.from_rotvec(bite_axis_unit * np.radians(deg))
                pos_try = rot_try.apply(centered) + t_center
                score = _dummy_alignment_score(pos_try)
                if score > best_score:
                    best_score = score
                    best_angle = float(np.radians(deg))

            if best_angle != 0.0:
                rot_best = Rotation.from_rotvec(bite_axis_unit * best_angle)
                positions_aligned = rot_best.apply(centered) + t_center
        else:
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
    #
    # Also reject on non-binding H atoms within 1.8 Å of the metal centre.
    # A genuine metal hydride's H is a BINDING atom (excluded from this check
    # by definition), so hydride complexes are unaffected. A non-binding H
    # this close to the metal is a spurious clash (e.g. DIPAMP: isolated
    # P···P 4.408 Å vs target 3.182 Å bite-delta of 1.226 Å folds backbone H
    # atoms to 1.39–1.65 Å) that gets misperceived as a Rh–H hydride on the
    # XYZ→OIN round-trip; the DG fallback places these H atoms cleanly
    # (~3.19 Å) instead.
    # Zone-A-P dummy atoms are intentionally placed metal-facing (that is
    # their whole purpose) and are stripped before this fragment is ever
    # combined with the metal -- exclude them from the clash checks below,
    # which police REAL non-binding atoms only.
    _dummy_idx_set = set(dummy_by_p_idx.values())

    if len(binding_idxs) >= 2:
        _binding_set = set(binding_idxs)
        _nb_heavy_idxs = [
            i
            for i in range(n_atoms)
            if i not in _binding_set and i not in _dummy_idx_set and symbols[i] != "H"
        ]
        if _nb_heavy_idxs:
            _min_d_metal = float(np.linalg.norm(positions_aligned[_nb_heavy_idxs], axis=1).min())
            if _min_d_metal < 1.7:
                return None  # Bidentate distortion too severe → fall back to DG

        _nb_h_idxs = [
            i
            for i in range(n_atoms)
            if i not in _binding_set and i not in _dummy_idx_set and symbols[i] == "H"
        ]
        if _nb_h_idxs:
            _min_d_metal_h = float(np.linalg.norm(positions_aligned[_nb_h_idxs], axis=1).min())
            if _min_d_metal_h < 1.8:
                return None  # Non-binding H too close to metal → fall back to DG

    positions_aligned, symbols, mol = _strip_dummy_atoms(
        positions_aligned, symbols, mol, list(dummy_by_p_idx.values())
    )
    assert not any(mol.GetAtomWithIdx(i).GetAtomicNum() == 0 for i in range(mol.GetNumAtoms())), (
        "dummy atom leaked past _stitch_fragment (bidentate path)"
    )
    assert len(positions_aligned) == mol.GetNumAtoms() == len(symbols), (
        "positions/mol/symbols length mismatch after dummy strip (bidentate path)"
    )
    return positions_aligned, symbols, mol


def _template_generate(
    parsed_oin: "ParsedOIN",
) -> "tuple[str, Chem.Mol | None, list] | None":
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
    # Stereo Phase 3: per-ring haptic-face correction decisions, collected
    # across all fragments and surfaced via GeneratedStructure (see
    # OIN3DGenerator.generate()'s sole call site of this function).
    haptic_decisions: list[dict] = []
    # Stereo Phase 4 (MiniPRD-B): per-"normal" (non-eta) fragment placement
    # metadata, enough to re-run `_stitch_fragment` with a new seed and
    # re-splice its atom range if Zone-A P verification (`_verify_zone_a_p`)
    # finds a mismatch. Zone-A P atoms are always metal-bound directly (never
    # eta/shared-slot), so only this path needs tracking.
    normal_frag_meta: list[dict] = []

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
                has_multi_eta = True  # noqa: F841  (state flag, kept for readability)
                result = _stitch_multi_eta_fragment(
                    frag_smiles,
                    vecs,
                    metal_sym,
                )
                if result is None:
                    return None
                frag_positions, frag_symbols, frag_mol, ring_decisions = result
                for rd in ring_decisions:
                    haptic_decisions.append({"fragment_idx": frag_idx, **rd})
                all_binding_idxs = [v.atom_in_fragment_idx for v in vecs]
                all_pos.extend(frag_positions)
                all_syms.extend(frag_symbols)
                all_frag_idxs.extend([frag_idx] * len(frag_positions))
                # Do NOT add to eta_frag_ranges: the Kabsch alignment in
                # _stitch_multi_eta_fragment already optimises placement
                # for all eta groups simultaneously. Adding a single-axis
                # ring-rotation entry would degrade the other ring's placement.
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

            # Stereo Phase 3: the single non-None OINVector.winding in this
            # ring's slot group is the target marker (never winding_by_slot
            # -- see negative constraints). Multi-marker same-slot is a
            # canonical-form violation -> ValueError, never a silent pick.
            target_winding, star_atom_idx = _extract_ring_winding_marker(vecs, eta_binding_idxs)

            # Guard the oin_parser.py template-gating hole (a SlotAssignment
            # can be dropped from `parsed_oin.vectors` -- and hence from
            # `vecs` here -- when its slot index falls outside the resolved
            # template, while still being recorded in `winding_by_slot`).
            # Fail loudly rather than silently treating a lost marker as a
            # legitimate zero-marker ring.
            matched_slot_idx = _find_slot_index_for_direction(template, unique_dirs[0])
            if matched_slot_idx is not None:
                slot_marker = parsed_oin.winding_by_slot.get(matched_slot_idx)
                if slot_marker is not None and target_winding is None:
                    raise AssertionError(
                        f"eta ring at fragment {frag_idx} (slot {matched_slot_idx}) "
                        f"lost its winding marker {slot_marker!r}: winding_by_slot "
                        "records a marker but no surviving OINVector for this ring "
                        "carries it (oin_parser.py template-gating hole)."
                    )

            result = _stitch_eta_fragment(
                frag_smiles,
                eta_binding_idxs,
                eta_slot_unit,
                metal_sym,
                winding=target_winding,
                star_atom_idx=star_atom_idx,
            )
            if result is None:
                return None  # Eta fragment failed → use DG
            frag_positions, frag_symbols, frag_mol, ring_decision = result
            haptic_decisions.append({"fragment_idx": frag_idx, **ring_decision})
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
        atom_start = len(all_pos)
        all_pos.extend(frag_positions)
        all_syms.extend(frag_symbols)
        all_frag_idxs.extend([frag_idx] * len(frag_positions))
        atom_end = len(all_pos)
        fragment_mol_parts.append((frag_mol, list(binding_idxs)))
        normal_frag_meta.append(
            {
                "frag_idx": frag_idx,
                "atom_start": atom_start,
                "atom_end": atom_end,
                "frag_smiles": frag_smiles,
                "binding_idxs": list(binding_idxs),
                "target_positions": list(target_positions),
                "slot_units_list": list(slot_units_list),
                "fragment_mol_parts_idx": len(fragment_mol_parts) - 1,
            }
        )

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

        # Final collision check after ring-rotation optimisation.
        #
        # Threshold note (Stereo Phase 3): lowered from the pre-Phase-3 1.60
        # to 1.45 Å. This is a *monotonic* loosening -- it can only convert a
        # previously-rejected (None -> DG fallback) placement into an
        # accepted one; it can never change the *geometry* chosen for any
        # placement that already cleared 1.60 (byte-identical for all
        # existing goldens, incl. ferrocene -- see test_winding_inertness.py).
        # Needed because retaining eta-ring substituents (Attempt 1b in
        # `_stitch_eta_fragment`) exposed a real, resolution-independent
        # steric floor for the heavily-substituted Ferrocene-halide-face
        # golden fixture (best achievable inter-ring separation ~1.494 Å,
        # confirmed identical under a much finer/multi-pass search -- this
        # is the true geometric optimum for this fixture's bond-length model,
        # not a search-quality artifact).
        final_min = _inter_frag_min(all_pos, all_frag_idxs)
        if final_min < 1.45:
            return None

    # Build combined RDKit mol with bond connectivity + 3D conformer from all_pos.
    # CombineMols preserves existing ETKDG conformers from fragments; we strip
    # all of them and set a single conformer from the final all_pos array so
    # the written MOL/SDF has the correct template-placed positions.
    combined_mol: Chem.Mol | None = None
    if has_all_mols and fragment_mol_parts:
        combined_mol = _assemble_combined_mol(metal_mol, fragment_mol_parts, all_pos)

    # ── Stereo Phase 4 (MiniPRD-B): Zone-A P verify-and-re-embed ───────────
    # Runs on the assembled-complex mol ONLY (Resolved Q2/RISK-2 -- never in
    # OIN3DGenerator/engine.py). `combined_mol is None` (eta fallback, or a
    # combine-step exception) is handled at the `generate()` call site via
    # `_warn_zone_a_p_fallback` -- Task 5 (RISK-9), never a silent skip.
    if combined_mol is not None and normal_frag_meta:
        zone_a_targets: list[tuple[int, str]] = []
        for meta in normal_frag_meta:
            for local_p_idx, expected_label in _zone_a_p_expected_labels(meta["frag_smiles"]):
                zone_a_targets.append((meta["atom_start"] + local_p_idx, expected_label))

        if zone_a_targets:
            mismatched = _verify_zone_a_p(combined_mol, zone_a_targets)
            # Stereo Phase 4 (MiniPRD-C spike finding): a Zone-A-P fragment's
            # residual-DOF clash sweep only sees fragments placed BEFORE it
            # (forbidden_positions is a snapshot at ITS placement time) -- so
            # it can clash with a fragment placed LATER. Reuse this SAME
            # retry loop (now that every fragment has been placed at least
            # once, forbidden_positions is complete) to re-run the sweep
            # with full information. Purely geometric; never touches stereo.
            clashing_frags = _zone_a_p_clash_offending_frags(normal_frag_meta, all_pos, all_syms)
            attempts = 0
            # Hard cap: 3 re-embed attempts total (SuperPRD B8/RISK-8) -- a
            # stereo mismatch must never become a timeout death. Each attempt
            # is a single new-seed ETKDG re-embed of only the offending
            # fragment(s); never a mirror/improper transform (B2/B3), so any
            # co-resident stereocenter in the SAME fragment (e.g. DIPAMP's
            # other Zone-A P) is never touched by this loop.
            while (mismatched or clashing_frags) and attempts < 3:
                attempts += 1
                offending_frag_idxs = {
                    meta["frag_idx"]
                    for meta in normal_frag_meta
                    if any(meta["atom_start"] <= gidx < meta["atom_end"] for gidx in mismatched)
                } | clashing_frags
                for meta in normal_frag_meta:
                    if meta["frag_idx"] not in offending_frag_idxs:
                        continue
                    redo = _stitch_fragment(
                        meta["frag_smiles"],
                        meta["binding_idxs"],
                        meta["target_positions"],
                        slot_units=meta["slot_units_list"],
                        forbidden_positions=list(all_pos),
                        seed=42 + attempts * 1009,
                    )
                    if redo is None:
                        continue  # keep previous placement for this fragment
                    new_positions, _new_symbols, new_frag_mol = redo
                    start, end = meta["atom_start"], meta["atom_end"]
                    if len(new_positions) == (end - start):
                        for k, gi in enumerate(range(start, end)):
                            all_pos[gi] = new_positions[k]
                        fragment_mol_parts[meta["fragment_mol_parts_idx"]] = (
                            new_frag_mol,
                            meta["binding_idxs"],
                        )
                combined_mol = _assemble_combined_mol(metal_mol, fragment_mol_parts, all_pos)
                if combined_mol is None:
                    break
                mismatched = _verify_zone_a_p(combined_mol, zone_a_targets)
                clashing_frags = _zone_a_p_clash_offending_frags(
                    normal_frag_meta, all_pos, all_syms
                )

            if mismatched:
                # Persistent mismatch after 3 attempts: emit the structure
                # anyway + warn (never a timeout death, never silent -- B8).
                for gidx in mismatched:
                    warnings.warn(
                        OINStereoWarning(
                            f"atom {gidx}: Zone-A P lone-pair CIP could not be "
                            f"enforced to match the OIN-encoded tag after "
                            f"{attempts} re-embed attempt(s) -- emitting the "
                            "structure as generated."
                        ),
                        stacklevel=2,
                    )

    lines = [str(n), f"Template-generated from OIN ({geo_code})"]
    for sym, pos in zip(all_syms, all_pos):
        lines.append(f"{sym:<2}  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}")
    xyz_str = "\n".join(lines) + "\n"

    return xyz_str, combined_mol, haptic_decisions


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
                Chem.SanitizeMol(
                    mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                )
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
    """Generate an XYZ block via RDKit ETKDGv3 as a Molassembler fallback.

    Used for topologies that Molassembler DG cannot embed (e.g. octahedral
    bidentate complexes). Returns the same dict schema as ``_molassembler_worker``.
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
                ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
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

    Returns:
    -------
    np.ndarray of shape (N_atoms, 3) or None if all conformers failed.
    """
    import scine_molassembler as masm  # noqa: PLC0415

    results = masm.dg.generate_ensemble(mol, n, seed)
    scored = [
        (pos, _min_inter_atomic_dist(pos)) for pos in results if not isinstance(pos, masm.dg.Error)
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

    Returns:
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

    Returns:
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
                                        mol.assign_stereopermutator(metal_idx, try_perm)
                                    except Exception:
                                        break  # no more valid perms
                                    test_res = None
                                    for try_seed in range(seed, seed + 5):
                                        r = masm.dg.generate_conformation(mol, try_seed)
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
                            elif expected_bindings and geo_code not in ("SPL",) and n_perms > 1:
                                # Exact-slot-only feedback for geometries without
                                # anti-parallel trans pairs (e.g. TPY).  SPL is
                                # excluded because _pick_masm_permutation already
                                # handles CIS/TRANS for square-planar complexes.
                                for try_perm in range(n_perms):
                                    try:
                                        mol.assign_stereopermutator(metal_idx, try_perm)
                                    except Exception:
                                        break
                                    test_res = None
                                    for try_seed in range(seed, seed + 5):
                                        r = masm.dg.generate_conformation(mol, try_seed)
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
                                    mol.assign_stereopermutator(metal_idx, safe_perm)
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

    Notes:
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
                Chem.SanitizeMol(
                    lig_mol,
                    Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
                )
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
                        _bsum = int(sum(b.GetBondTypeAsDouble() for b in _a.GetBonds()))
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
    _SKIP_AROM = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
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

    Returns:
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

    from collections import defaultdict

    import numpy as np

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
        """Initialize the adapter with subprocess timeout and DG generation options."""
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

        Returns:
        -------
        str
            XYZ block string of the generated conformer.

        Raises:
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
                xyz_str, template_mol, haptic_decisions = template_result
                if template_mol is None:
                    # Stereo Phase 4 (Task 5, RISK-9): no assembled mol means
                    # `_template_generate`'s Zone-A P verify-and-re-embed loop
                    # never ran (it is gated on a non-None combined_mol) --
                    # e.g. an eta fragment's mol couldn't be built. Warn if
                    # the OIN carries any Zone-A P tag so the gap is visible.
                    _warn_zone_a_p_fallback(parsed_oin, "template path, no assembled mol")
                return GeneratedStructure(
                    xyz=xyz_str, mol=template_mol, haptic_face_decisions=haptic_decisions
                )

        # ── Fallback: Molassembler DG ───────────────────────────────────────
        # Stereo Phase 4 (MiniPRD-C, Task 5, C2): this path keeps its own
        # placement (incompatible-bite chelates are routed here by ee0b3f0's
        # bite-distortion guard), but Zone-A-P enforcement now ALSO runs
        # here via a bounded reseed-and-verify loop -- see
        # _zone_a_p_measured_labels_dg for the set-based comparison this
        # requires (the canonicalised connected SMILES makes exact per-atom
        # index tracking unavailable here, unlike the template path).
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

        zone_a_expected_labels = sorted(
            label
            for frag_idx, frag_smiles in enumerate(parsed_oin.fragments)
            if frag_idx != parsed_oin.metal_fragment_idx
            for _local_idx, label in _zone_a_p_expected_labels(frag_smiles)
        )

        bonded_mol: "Chem.Mol | None" = None
        xyz_block = ""
        dg_attempts = 0
        # Hard cap 3 attempts (mirrors the template path's bounded re-embed
        # loop, B8/RISK-8): a stereo mismatch must never become a timeout
        # death. Only reseeds when there is something to enforce; a plain
        # (non-Zone-A-P) DG generation always runs exactly once, unchanged.
        max_attempts = 3 if zone_a_expected_labels else 1
        while dg_attempts < max_attempts:
            dg_attempts += 1
            args = {
                "smiles": connected_smiles,
                "seed": seed if dg_attempts == 1 else seed + dg_attempts * 1009,
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
                    raise MolassemblerTimeoutError(f"Molassembler timed out after {self.timeout}s")

            if not result.get("ok"):
                raise RuntimeError(f"Molassembler error: {result.get('error', 'unknown')}")

            xyz_block = result["xyz_block"]
            bonded_mol = _reconstruct_mol_from_smiles_and_xyz(connected_smiles, xyz_block)

            if not zone_a_expected_labels:
                break  # nothing to enforce, first (only) attempt stands
            if bonded_mol is None:
                continue  # try a fresh seed; nothing to verify without a mol
            measured = _zone_a_p_measured_labels_dg(bonded_mol)
            if measured == zone_a_expected_labels:
                break  # enforced -- no warning

        if zone_a_expected_labels:
            if bonded_mol is None:
                _warn_zone_a_p_fallback(
                    parsed_oin, "Molassembler DG fallback path, no reconstructed mol"
                )
            elif _zone_a_p_measured_labels_dg(bonded_mol) != zone_a_expected_labels:
                warnings.warn(
                    OINStereoWarning(
                        "Zone-A P lone-pair CIP could not be enforced to match the "
                        f"OIN-encoded tag(s) on the Molassembler DG fallback path "
                        f"after {dg_attempts} attempt(s) -- emitting the structure "
                        "as generated."
                    ),
                    stacklevel=2,
                )

        return GeneratedStructure(xyz=xyz_block, mol=bonded_mol)
