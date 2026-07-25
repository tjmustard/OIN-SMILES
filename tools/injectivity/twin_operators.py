"""The twin operators the Y1 plan called for and never built: ``swap_donor``, ``invert_stereocenter``.

:mod:`tools.injectivity.twin_collision` ships exactly one operator, ``mirror_z``. A mirror is a
coordinate transform, so it is trivially valid: bond lengths, angles and clashes are all
preserved by construction. It also answers exactly one question -- *enantiomerism* -- and
whole-molecule enantiomerism at that. Two questions it cannot reach:

* **Diastereomerism.** A square-planar complex is planar and therefore achiral, so its mirror
  is never a new isomer however different its four donors are. What distinguishes ``@SP1`` from
  ``@SP2`` from ``@SP3`` is *which donors sit trans*, and only a donor exchange produces that.
  This is why the ``@SP`` half of the metal-stereo question had no instrument at all.
* **One axis at a time.** Mirroring a molecule with several stereo elements flips them all
  together, so it cannot say which one the encoder lost. Lane 4's multi-axis cases need a
  single-axis flip.

Both operators here are **structural edits**, so unlike a mirror they can produce nonsense --
a swapped bulky ligand driven into another, an inverted centre with a broken ring. Every
operator therefore returns a :class:`EditedTwin` carrying the vdW clash counts before and
after, and :func:`probe_operator` refuses to score a twin the clash gate rejects. An oracle
verdict on a bad geometry is meaningless.

Distinctness is certified the same way throughout: the twin is a distinct isomer when no
proper rotation, over the graph automorphisms **and over the structure's whole torsion orbit**,
superimposes it on the original (:mod:`tools.injectivity.torsion_oracle`). Using the torsion
orbit rather than a rigid fit matters here -- a donor swap moves whole ligands, and a rigid
comparison would call every swap "distinct" simply because the substituents landed differently.

Run:
  PYTHONPATH=$PWD/src python -m tools.injectivity.twin_operators <file.xyz>
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem

from oinsmiles.core.constants import TRANSITION_METALS_NUM
from oinsmiles.generator3d.clash import vdw_clash_count
from oinsmiles.oin.axial import detect_axial_axes

from .oracle import load_mol
from .torsion_oracle import (
    _automorphism_perms,
    _batch_proper_rmsd,
    _heavy_indices,
    _optimise,
    _Scorer,
    _side_without,
    dihedral_negating_theta,
    rotatable_torsions,
)

#: a twin is rejected when the edit introduces vdW clashes the original did not have.
CLASH_TOLERANCE = 0


@dataclass
class EditedTwin:
    """A structural edit, with everything needed to decide whether to trust it."""

    operator: str
    detail: str
    coords: np.ndarray = field(repr=False)
    clash_before: tuple = (0, 0)
    clash_after: tuple = (0, 0)
    error: str = ""
    #: bonds the distinctness search must NOT rotate, because the edit changed them
    freeze_bonds: tuple = ()

    @property
    def geometry_ok(self) -> bool:
        """The edit did not introduce clashes the original was free of."""
        if self.error:
            return False
        return self.clash_after[0] - self.clash_before[0] <= CLASH_TOLERANCE


def _clash(mol: Chem.Mol, coords: np.ndarray) -> tuple:
    c, s, _w = vdw_clash_count(coords, [a.GetAtomicNum() for a in mol.GetAtoms()])
    return (c, s)


def relax_torsions(
    mol: Chem.Mol, coords: np.ndarray, *, freeze_bonds=(), sweeps: int = 2, grid: int = 12
) -> np.ndarray:
    """Reduce vdW clashes by rotating bonds only -- never by moving a configuration.

    A structural edit lands the moved fragment at whatever rotamer it happened to have, so a
    twin can be rejected for a clash that a plain bond rotation would remove. Relaxing in
    *torsion space* is the one repair that is safe here: it cannot change any configuration, so
    it cannot quietly turn the twin back into the original. Bonds the edit itself changed are
    frozen. Objective is ``(clash_vdw, clash_severe, -worst_overlap)``, minimised greedily.
    """
    from .torsion_oracle import apply_torsions

    root = int(np.argmax([a.GetAtomicNum() for a in mol.GetAtoms()]))
    frozen = {frozenset(b) for b in freeze_bonds}
    tors = [t for t in rotatable_torsions(mol, root) if frozenset((t[0], t[1])) not in frozen]
    if not tors:
        return coords
    z = [a.GetAtomicNum() for a in mol.GetAtoms()]

    def score(theta):
        c, s, w = vdw_clash_count(apply_torsions(coords, tors, theta), z)
        return (c, s, -w)

    theta = np.zeros(len(tors))
    best = score(theta)
    for _ in range(sweeps):
        improved = False
        for t in range(len(tors)):
            keep = theta[t]
            for step in np.linspace(0.0, 360.0, grid, endpoint=False):
                theta[t] = step
                v = score(theta)
                if v < best:
                    best, keep, improved = v, step, True
            theta[t] = keep
        if not improved:
            break
    return apply_torsions(coords, tors, theta)


def _metal_idx(mol: Chem.Mol) -> int:
    metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM]
    if len(metals) != 1:
        raise ValueError(f"expected exactly one transition metal, found {len(metals)}")
    return metals[0]


def _rotation_taking(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """The minimal rotation matrix taking unit vector ``u`` onto unit vector ``v``."""
    c = float(np.clip(np.dot(u, v), -1.0, 1.0))
    if c > 1.0 - 1e-12:
        return np.eye(3)
    axis = np.cross(u, v)
    n = float(np.linalg.norm(axis))
    if n < 1e-12:  # antiparallel: any perpendicular axis, 180 deg
        axis = np.cross(u, np.array([1.0, 0.0, 0.0]))
        if np.linalg.norm(axis) < 1e-9:
            axis = np.cross(u, np.array([0.0, 1.0, 0.0]))
        n = float(np.linalg.norm(axis))
    k = axis / n
    s = float(np.sqrt(max(0.0, 1.0 - c * c))) if n > 1e-12 else 0.0
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + s * K + (1.0 - c) * (K @ K)


def _exchange_branches(coords, pivot, head_a, head_b, branch_a, branch_b):
    """Swap two branches around ``pivot`` by a RIGID motion of each.

    Each branch is rotated so its head points where the other's did, then *translated* along
    that new direction so the head sits at the other's distance from the pivot. Rotation plus
    translation is an isometry, so every distance inside a branch survives exactly.

    The obvious shortcut -- rotate and then **scale** by the distance ratio -- is a similarity,
    not an isometry, and it silently shrinks the branch: swapping an N-H against an N-Pt branch
    that way compressed every bond in the moved branch by half.
    """
    ua = coords[head_a] - pivot
    ub = coords[head_b] - pivot
    ra, rb = float(np.linalg.norm(ua)), float(np.linalg.norm(ub))
    ua, ub = ua / ra, ub / rb
    new = coords.copy()
    new[branch_a] = (coords[branch_a] - pivot) @ _rotation_taking(ua, ub).T + pivot + (rb - ra) * ub
    new[branch_b] = (coords[branch_b] - pivot) @ _rotation_taking(ub, ua).T + pivot + (ra - rb) * ua
    return new


# --- operator 1: swap_donor ----------------------------------------------------------


def donor_groups(mol: Chem.Mol) -> dict[int, list[int]]:
    """``{donor_atom_idx: [atom indices of the branch that hangs off it]}``.

    A branch is what detaches when the metal-donor **bond** is cut. For a chelate that is the
    whole chelate (its other donor bond still holds it), which is exactly why a chelate donor
    cannot be swapped independently -- :func:`swap_donor` refuses those.
    """
    m = _metal_idx(mol)
    n = mol.GetNumAtoms()
    out = {}
    for nb in mol.GetAtomWithIdx(m).GetNeighbors():
        d = nb.GetIdx()
        side = _side_without(mol, n, cut=(m, d), start=d)
        if side is not None and m not in side:
            out[d] = sorted(side)
    return out


def swap_donor(mol: Chem.Mol, coords: np.ndarray, donor_a: int, donor_b: int) -> EditedTwin:
    """Exchange the coordination sites of two donors, moving each ligand rigidly.

    Each ligand is rotated about the metal by the minimal rotation carrying its own donor
    direction onto the other's, and its metal-donor distance is preserved. Everything internal
    to each ligand -- bond lengths, angles, torsions, and any stereocentre it carries -- is
    untouched, so the only thing that changes is *which site each donor occupies*. That is a
    positional isomerisation, not a deformation.

    Refuses a donor whose branch is not detachable by cutting its own metal bond, i.e. a
    chelate donor: swapping one arm of a chelate is not a rigid motion of anything.
    """
    m = _metal_idx(mol)
    groups = donor_groups(mol)
    for d in (donor_a, donor_b):
        if d not in groups:
            return EditedTwin(
                "swap_donor",
                f"{donor_a}<->{donor_b}",
                coords,
                error=f"donor {d} belongs to a chelate: not independently movable",
            )
    if set(groups[donor_a]) & set(groups[donor_b]):
        return EditedTwin(
            "swap_donor",
            f"{donor_a}<->{donor_b}",
            coords,
            error="both donors belong to the same ligand",
        )

    new = _exchange_branches(coords, coords[m], donor_a, donor_b, groups[donor_a], groups[donor_b])
    return EditedTwin(
        "swap_donor",
        f"{mol.GetAtomWithIdx(donor_a).GetSymbol()}{donor_a}"
        f"<->{mol.GetAtomWithIdx(donor_b).GetSymbol()}{donor_b}",
        new,
        clash_before=_clash(mol, coords),
        clash_after=_clash(mol, new),
    )


def enumerate_donor_swaps(mol: Chem.Mol) -> list[tuple[int, int]]:
    """Every swappable, non-equivalent donor pair, cheapest-to-interpret first.

    Pairs of *symmetry-equivalent* donors are skipped: exchanging them is the identity on the
    isomer, so a "collision" there would be correct behaviour, not a blind spot. This is the
    trap the hand-written donor-swap probe fell into.
    """
    groups = donor_groups(mol)
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    ds = sorted(groups)
    out = []
    for i, a in enumerate(ds):
        for b in ds[i + 1 :]:
            if ranks[a] == ranks[b] and len(groups[a]) == len(groups[b]):
                continue
            out.append((a, b))
    return out


# --- operator 2: invert_stereocenter --------------------------------------------------


def invert_axial(mol: Chem.Mol, coords: np.ndarray, axis_index: int = 0) -> EditedTwin:
    """Negate the signed dihedral of ONE hindered biaryl axis, leaving every other alone.

    This is the single-axis flip a whole-molecule mirror cannot do. On a molecule carrying two
    axes it converts the like (R,R) diastereomer into the unlike (R,S) one -- so it, not
    ``mirror_z``, is the operator that isolates which axis an encoder lost.
    """
    axes = [ax for ax in detect_axial_axes(mol) if ax.hindered]
    if axis_index >= len(axes):
        return EditedTwin(
            "invert_axial",
            f"axis {axis_index}",
            coords,
            error=f"only {len(axes)} hindered axes present",
        )
    ax = axes[axis_index]
    root = int(np.argmax([a.GetAtomicNum() for a in mol.GetAtoms()]))
    tors = rotatable_torsions(mol, root)
    match = [t for t, (a, b, _m, _d) in enumerate(tors) if {a, b} == {ax.a1, ax.a2}]
    if not match:
        return EditedTwin(
            "invert_axial",
            f"axis {ax.a1}-{ax.a2}",
            coords,
            error="the axis bond is not an independently rotatable torsion (ring-locked)",
        )
    theta = np.zeros(len(tors))
    theta[match[0]] = dihedral_negating_theta(mol, coords, tors)[match[0]]
    from .torsion_oracle import apply_torsions

    new = apply_torsions(coords, tors, theta)
    return EditedTwin(
        "invert_axial",
        f"axis {ax.a1}-{ax.a2} (sign {ax.sign} -> {-ax.sign})",
        new,
        clash_before=_clash(mol, coords),
        clash_after=_clash(mol, new),
        freeze_bonds=((ax.a1, ax.a2),),
    )


def invert_tetrahedral(mol: Chem.Mol, coords: np.ndarray, center: int) -> EditedTwin:
    """Invert one tetrahedral centre by exchanging two of its substituent branches.

    Swapping two substituents is the textbook parity flip, and as a *rigid* exchange of two
    branches it leaves both branches' internal geometry intact -- so the edit changes the
    configuration at ``center`` and nothing else.

    Two refusals keep the operator from producing a twin that is not a twin.

    * The branches must be **detachable**: a ring bond cannot be moved rigidly.
    * The centre must be a **genuine stereocentre** -- every substituent carrying a distinct
      symmetry rank, the same test :func:`tools.injectivity.config_oracle.bound_amine_centers`
      applies. Without it the operator happily "inverts" a metal-bound ammine by exchanging one
      of its three equivalent hydrogens for the metal, producing the identical molecule and
      reporting success. That is the hand-written donor-swap mistake in a different costume.

    Every eligible exchange is generated and ranked by the clash it introduces, and the best is
    then relaxed in torsion space (:func:`relax_torsions`) before the gate rules on it.

    **Measured scope limit, stated because it decides when to reach for this operator.** A rigid
    substituent exchange only works when the centre has room. Every stereocentre in the current
    fixture set is *locked* -- inside a chelate ring (BDPP, DPDME) or bound to the metal
    (POJJOP) -- and for all of them the exchange drives a substituent into the coordination
    sphere: the inversion itself is correct (POJJOP's signed tetrahedral volume flips sign) but
    the clash gate rejects the geometry, torsion relaxation notwithstanding. Reaching those
    would need bond-angle relaxation, which a generator-free instrument must not do. For a
    locked centre use ``mirror_z`` (whole-molecule) or :func:`invert_axial` (single axis); this
    operator is for stereocentres on a freely-rotating pendant.
    """
    n = mol.GetNumAtoms()
    atom = mol.GetAtomWithIdx(center)
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    nbrs = list(atom.GetNeighbors())
    if len({ranks[nb.GetIdx()] for nb in nbrs}) != len(nbrs):
        return EditedTwin(
            "invert_tetrahedral",
            f"atom {center}",
            coords,
            error="substituents are not all symmetry-distinct: not a stereocentre",
        )
    branches = {}
    for nb in nbrs:
        v = nb.GetIdx()
        side = _side_without(mol, n, cut=(center, v), start=v)
        if side is not None and center not in side:
            branches[v] = sorted(side)
    if len(branches) < 2:
        return EditedTwin(
            "invert_tetrahedral",
            f"atom {center}",
            coords,
            error=f"only {len(branches)} detachable substituents (ring-locked centre)",
        )
    # generate every candidate exchange and let the clash gate choose, rather than committing
    # to one pair and hoping: the smallest two branches are the least disruptive on average but
    # not always, and an edit the gate rejects is a probe that cannot be scored at all.
    order = sorted(branches, key=lambda v: (len(branches[v]), v))
    before = _clash(mol, coords)
    candidates = []
    for i, x in enumerate(order):
        for y in order[i + 1 :]:
            new = _exchange_branches(coords, coords[center], x, y, branches[x], branches[y])
            candidates.append((_clash(mol, new)[0] - before[0], x, y, new))
    candidates.sort(key=lambda c: (c[0], len(branches[c[1]]) + len(branches[c[2]])))
    _delta, a, b, new = candidates[0]
    if _clash(mol, new)[0] > before[0]:
        # a clash may just be a bad rotamer of the moved branch; try to rotate it away before
        # giving up. Torsion-only, so the inversion itself cannot be undone.
        relaxed = relax_torsions(mol, new)
        if _clash(mol, relaxed)[0] < _clash(mol, new)[0]:
            new = relaxed
    return EditedTwin(
        "invert_tetrahedral",
        f"atom {center} ({atom.GetSymbol()}): swapped branches at {a} and {b}",
        new,
        clash_before=_clash(mol, coords),
        clash_after=_clash(mol, new),
    )


# --- distinctness + encoder probe ------------------------------------------------------


def torsion_orbit_distance(
    mol: Chem.Mol,
    coords: np.ndarray,
    target: np.ndarray,
    *,
    restarts: int = 8,
    sweeps: int = 4,
    grid: int = 12,
    seed: int = 42,
    freeze_bonds: tuple = (),
) -> tuple[float, float]:
    """``(d_target, d_control)``: how close ``coords``' torsion orbit gets to ``target``.

    Same machinery and the same paired positive control as
    :func:`tools.injectivity.torsion_oracle.configurational_verdict`, but against an arbitrary
    target instead of the mirror -- which is what a *diastereomer* comparison needs.

    ``freeze_bonds`` removes named bonds from the search, and an operator that edits a torsion
    **must** use it. Otherwise the orbit simply rotates the edit back: an atropisomer flip is a
    torsion change, so a free search reaches it in one step and every axial twin would be
    scored "not distinct". Atropisomerism is conformational *plus a rotational barrier*, and a
    barrier is exactly what this geometric search does not model -- see the caveat in
    :mod:`tools.injectivity.torsion_oracle`.
    """
    from .torsion_oracle import apply_torsions

    heavy = _heavy_indices(mol)
    perms = _automorphism_perms(mol, heavy, 4000)
    tors = rotatable_torsions(mol, int(np.argmax([a.GetAtomicNum() for a in mol.GetAtoms()])))
    frozen = {frozenset(b) for b in freeze_bonds}
    tors = [t for t in tors if frozenset((t[0], t[1])) not in frozen]
    if not tors:
        return float(_batch_proper_rmsd(target[heavy], coords[heavy], perms).min()), 0.0
    d, _theta = _optimise(
        _Scorer(coords, heavy, perms, target[heavy], tors),
        np.random.default_rng(seed),
        restarts=restarts,
        sweeps=sweeps,
        grid=grid,
        seeds=(dihedral_negating_theta(mol, coords, tors),),
    )
    ctrl = apply_torsions(coords, tors, np.random.default_rng(seed + 1).uniform(0, 360, len(tors)))
    dc, _ = _optimise(
        _Scorer(coords, heavy, perms, ctrl[heavy], tors),
        np.random.default_rng(seed + 2),
        restarts=restarts,
        sweeps=sweeps,
        grid=grid,
    )
    return d, dc


def _write_xyz(mol: Chem.Mol, coords: np.ndarray, dst: Path) -> Path:
    lines = [str(mol.GetNumAtoms()), "twin operator output"]
    for i, a in enumerate(mol.GetAtoms()):
        x, y, z = coords[i]
        lines.append(f"{a.GetSymbol():<3} {x:>14.8f} {y:>14.8f} {z:>14.8f}")
    dst.write_text("\n".join(lines) + "\n")
    return dst


@dataclass
class OperatorOutcome:
    name: str
    operator: str
    detail: str
    geometry_ok: bool
    clash_before: tuple
    clash_after: tuple
    oracle_distinct: bool | None
    d_twin: float
    d_control: float
    oin_base: str
    oin_twin: str
    raw_equal: bool | None
    key_equal: bool | None
    verdict: str

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["clash_before"] = list(self.clash_before)
        d["clash_after"] = list(self.clash_after)
        return d


def probe_operator(
    xyz_path: str | Path,
    twin: EditedTwin,
    *,
    charge: int = 0,
    threshold: float = 0.5,
    restarts: int = 8,
    name: str | None = None,
) -> OperatorOutcome:
    """Encode a structure and one structural-edit twin; classify the collision, if any."""
    import tempfile

    from oinsmiles import XYZToSMILES

    from .twin_collision import _key, _silence_fds, classify

    path = Path(xyz_path)
    name = name or path.stem
    mol, coords = load_mol(path, charge)

    if not twin.geometry_ok:
        return OperatorOutcome(
            name=name,
            operator=twin.operator,
            detail=twin.detail,
            geometry_ok=False,
            clash_before=twin.clash_before,
            clash_after=twin.clash_after,
            oracle_distinct=None,
            d_twin=float("nan"),
            d_control=float("nan"),
            oin_base="",
            oin_twin="",
            raw_equal=None,
            key_equal=None,
            verdict=twin.error or "rejected: edit introduced vdW clashes",
        )

    d_twin, d_control = torsion_orbit_distance(
        mol, coords, twin.coords, restarts=restarts, freeze_bonds=twin.freeze_bonds
    )
    distinct = d_twin > threshold and d_control <= threshold

    with tempfile.TemporaryDirectory() as tmp:
        dst = _write_xyz(mol, twin.coords, Path(tmp) / "twin.xyz")
        with _silence_fds():
            oin_base = XYZToSMILES().convert(str(path))
            oin_twin = XYZToSMILES().convert(str(dst))

    raw_equal = oin_base == oin_twin
    key_equal = _key(oin_base) == _key(oin_twin)
    return OperatorOutcome(
        name=name,
        operator=twin.operator,
        detail=twin.detail,
        geometry_ok=True,
        clash_before=twin.clash_before,
        clash_after=twin.clash_after,
        oracle_distinct=distinct,
        d_twin=round(d_twin, 4),
        d_control=round(d_control, 4),
        oin_base=oin_base,
        oin_twin=oin_twin,
        raw_equal=raw_equal,
        key_equal=key_equal,
        verdict=classify(raw_equal, key_equal, distinct)
        if d_control <= threshold
        else "inconclusive (torsion control failed)",
    )


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    for p in argv:
        path = Path(p)
        mol, coords = load_mol(path)
        print(f"\n=== {path.stem} ===")
        for a, b in enumerate_donor_swaps(mol)[:6]:
            o = probe_operator(path, swap_donor(mol, coords, a, b))
            print(
                f"  swap_donor {o.detail:20s} geom_ok={o.geometry_ok} "
                f"distinct={o.oracle_distinct} (d={o.d_twin}/ctl {o.d_control}) "
                f"raw_eq={o.raw_equal} key_eq={o.key_equal} -> {o.verdict}"
            )
        for k in range(2):
            t = invert_axial(mol, coords, k)
            if t.error:
                break
            o = probe_operator(path, t)
            print(
                f"  invert_axial {o.detail:24s} distinct={o.oracle_distinct} "
                f"raw_eq={o.raw_equal} key_eq={o.key_equal} -> {o.verdict}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
