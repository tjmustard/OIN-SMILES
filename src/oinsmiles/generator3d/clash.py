"""Whole-complex van-der-Waals steric-clash perception.

One shared definition of a "vdW clash", used by every acceptance/ranking site in the
generator (``embed._finalize_positions``, ``clean_geometry.ff_clean``, the ``__init__``
pool dedup, and ``metallogen_adapter._select_by_geometry``). It is a faithful lift of the
release quality metric, ``tools/structure_distortion_report.py::geometry_metrics``, so the
gate that *rejects* a clash and the metric that *measures* it agree by construction: same
RDKit radii, same ``1.3*sum(R_cov)`` bond perception, same ``0.75`` overlap cutoff, same
1-2 / 1-3 / self exemptions. ``tests/unit/test_vdw_acceptance.py`` cross-checks the two on a
shared fixture to catch any drift (that tool is owned by another session, hence the copy).

A clash is a *non-bonded, non-geminal* atom pair inside vdW contact. Metal-donor and every
intra-ligand bond are exempt implicitly: they fall inside ``1.3*sum(R_cov)`` and read as
bonds. Fusion pre-checks upstream guarantee no two atoms sit within covalent-bonding
distance except real bonds, so the geometric perception never mislabels a genuine severe
overlap as a bond at these call sites.
"""

import os

import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist

_PT = Chem.GetPeriodicTable()

# vdW acceptance/selection term -- ON by default (promoted in A5, v0.4.3). When A3 first
# added it, the term reduced clash only by loosening coordination on the then-current pool
# (clash 62.5%->2.5% but round-trip 70%->42.5% -- real ligand ejection, e.g.
# [Ru_TBP] -> [Ru_TPL]), so it shipped gated OFF. A5 re-ran the four-arm A/B on the improved
# v0.4.3 pool (B1's electronic geometry prior supplies the tight-AND-clean conformers A3's
# pool lacked) over a 40-molecule worst-cohort sample: the term now drops the vdW-clash
# fraction 92.5%->5.1% with ZERO coordination-number regressions (0 ejections across all 39
# generatable molecules) and round-trip 92.5%->90.0% (the single delta is a TBP<->SPY label
# flip with all donors bound, not a loss). Goldens still round-trip under the gate. So the
# ejection mechanism is gone and the term is promoted to default.
# Disable per run via OIN_VDW_ACCEPTANCE=0, or by setting clash.VDW_ACCEPTANCE_ENABLED = False.
VDW_ACCEPTANCE_ENABLED = os.environ.get("OIN_VDW_ACCEPTANCE", "1") != "0"

# Stretched-bond acceptance term -- the *complement* of the vdW gate above (that flags
# non-bonded pairs too CLOSE; this flags bonded pairs too FAR). OFF by default (SL1, v0.4.4):
# it is an opt-in acceptance criterion for the generate-until-key-exact early-exit path, not a
# default-path change. A bond is "stretched" when its length exceeds STRETCH_RATIO * the sum of
# the two RDKit covalent radii -- 1.5 leaves normal covalent AND metal-donor/dative bonds well
# under the bar (Pt-Cl, Fe-C(eta) sit ~1.0) while catching the pathology it targets: a methyl
# threaded through a phenyl or interlocked rings tangling a bond out to ~4 Angstrom (ratio ~2.6).
# Enable per run via OIN_STRETCHED_BOND=1, or by setting clash.STRETCHED_BOND_ENABLED = True.
STRETCH_RATIO = 1.5
STRETCHED_BOND_ENABLED = os.environ.get("OIN_STRETCHED_BOND", "0") != "0"


def _rcov(z):
    """RDKit covalent radius for atomic number ``z`` (fallback 0.7 Angstrom)."""
    try:
        return _PT.GetRcovalent(int(z))
    except Exception:
        return 0.7


def _rvdw(z):
    """RDKit van-der-Waals radius for atomic number ``z`` (fallback 1.7 Angstrom)."""
    try:
        return _PT.GetRvdw(int(z))
    except Exception:
        return 1.7


def vdw_clash_count(
    positions,
    atomic_numbers,
    clash_cutoff=0.75,
    severe_cutoff=0.60,
    adj_factor=1.3,
):
    """Count non-bonded, non-geminal atom pairs inside van-der-Waals contact.

    ``positions`` is an ``(N, 3)`` array; ``atomic_numbers`` an ``(N,)`` sequence of Z,
    aligned index-for-index. Returns ``(clash_vdw, clash_severe, worst_overlap)`` where
    ``overlap = dist / (rvdw_i + rvdw_j)``: ``clash_vdw`` counts pairs below
    ``clash_cutoff``, ``clash_severe`` below ``severe_cutoff``, and ``worst_overlap`` is
    the smallest non-bonded overlap ratio (1.0 when there are no non-bonded pairs).
    """
    pos = np.asarray(positions, dtype=float)
    z = np.asarray(atomic_numbers)
    n = len(z)
    if n < 2:
        return 0, 0, 1.0

    dist = cdist(pos, pos)
    np.fill_diagonal(dist, 9e9)
    rcov = np.array([_rcov(zi) for zi in z])
    rvdw = np.array([_rvdw(zi) for zi in z])

    adj = dist < adj_factor * (rcov[:, None] + rcov[None, :])
    overlap = dist / (rvdw[:, None] + rvdw[None, :])

    a = adj.astype(int)
    geminal = (a @ a) > 0
    nonbond = ~(adj | geminal | np.eye(n, dtype=bool))

    iu = np.triu_indices(n, 1)
    nb = nonbond[iu]
    ov = overlap[iu]
    clash_vdw = int(((ov < clash_cutoff) & nb).sum())
    clash_severe = int(((ov < severe_cutoff) & nb).sum())
    worst_overlap = float(ov[nb].min()) if nb.any() else 1.0
    return clash_vdw, clash_severe, worst_overlap


def mol_clash_count(mol):
    """Whole-complex ``clash_vdw`` count for a generator ``Molecule``.

    Duck-typed: ``mol.atom_list`` of atoms exposing ``get_coordinate()`` and
    ``get_atomic_number()`` (metal + every ligand atom). Convenience for the pool-ranking
    and geometry-selection sites, which rank conformers by whole-complex steric clash.

    Returns 0 when the molecule's atoms/coordinates cannot be read (e.g. a bare RDKit
    ``Mol`` with no ``atom_list``), so a caller that ranks by clash degrades to its prior
    (energy / geometry-fit) order for that candidate rather than dropping it.
    """
    try:
        positions = [a.get_coordinate() for a in mol.atom_list]
        atomic_numbers = [a.get_atomic_number() for a in mol.atom_list]
    except AttributeError:
        return 0
    clash_vdw, _severe, _worst = vdw_clash_count(positions, atomic_numbers)
    return clash_vdw


def stretched_bond_count(positions, atomic_numbers, adjacency, stretch_ratio=STRETCH_RATIO):
    """Count *bonded* atom pairs stretched beyond ``stretch_ratio`` x sum of covalent radii.

    The complement of :func:`vdw_clash_count`: that flags *non-bonded* pairs inside vdW contact
    (too close); this flags *bonded* pairs pulled too far apart. ``positions`` is an ``(N, 3)``
    array, ``atomic_numbers`` an ``(N,)`` sequence of Z aligned index-for-index, and
    ``adjacency`` the molecule's TRUE connectivity as an ``(N, N)`` matrix (nonzero == bonded),
    e.g. a MetalloGen ``mol.get_adj_matrix()``. The geometric ``1.3*sum(R_cov)`` bond perception
    ``vdw_clash_count`` uses CANNOT be reused here: an over-stretched bond exceeds that cutoff and
    silently drops out of the perceived adjacency -- exactly the pair this metric must see. Both
    matrix orientations are honored (``adj | adj.T``) so a directed adjacency still works.

    Returns ``(n_stretched, worst_ratio)`` where ``ratio = dist / (rcov_i + rcov_j)``:
    ``n_stretched`` counts bonded pairs above ``stretch_ratio`` and ``worst_ratio`` is the
    largest bonded ratio (``0.0`` when there are no bonds).
    """
    pos = np.asarray(positions, dtype=float)
    z = np.asarray(atomic_numbers)
    n = len(z)
    if n < 2 or adjacency is None:
        return 0, 0.0
    adj = np.asarray(adjacency)
    if adj.shape != (n, n):
        return 0, 0.0

    dist = cdist(pos, pos)
    rcov = np.array([_rcov(zi) for zi in z])
    ratio = dist / (rcov[:, None] + rcov[None, :])

    bonded = (adj != 0) | (adj.T != 0)
    iu = np.triu_indices(n, 1)
    bp = bonded[iu]
    if not bp.any():
        return 0, 0.0
    rr = ratio[iu]
    n_stretched = int(((rr > stretch_ratio) & bp).sum())
    worst_ratio = float(rr[bp].max())
    return n_stretched, worst_ratio


def mol_stretched_bond_count(mol):
    """Whole-complex stretched-bond count for a generator ``Molecule``.

    Complement of :func:`mol_clash_count`. Duck-typed: reads ``mol.atom_list`` (atoms exposing
    ``get_coordinate()`` / ``get_atomic_number()``) and the true adjacency from
    ``mol.get_adj_matrix()``. Returns 0 when the atoms, coordinates, or adjacency cannot be read
    (a bare RDKit ``Mol``, or a ``Molecule`` with no perceived connectivity), so a caller that
    gates on stretched bonds degrades to accepting the candidate rather than dropping it.
    """
    try:
        positions = [a.get_coordinate() for a in mol.atom_list]
        atomic_numbers = [a.get_atomic_number() for a in mol.atom_list]
        adjacency = mol.get_adj_matrix()
    except AttributeError:
        return 0
    n_stretched, _worst = stretched_bond_count(positions, atomic_numbers, adjacency)
    return n_stretched
