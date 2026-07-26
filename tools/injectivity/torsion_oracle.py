"""Torsion-aware configurational oracle -- separate CONFORMATIONAL from CONFIGURATIONAL chirality.

The gap this closes. :func:`tools.injectivity.oracle.geometric_chirality` is a **rigid**
superposition test: it takes the minimum proper-rotation RMSD between a structure and its
mirror over the graph automorphisms. That is exactly right for a rigid species and exactly
wrong for a flexible one, because a floppy *achiral* molecule's mirror is a
non-superimposable **conformer** and reads as "chiral". Wave 1 hit this (a 30-structure
sample read ~97% chiral) and correctly refused to publish a rate; the Wave-3 unknown-unknown
hunt left 70 of 299 structures in an ambiguous residue for the same reason.

The distinction that matters:

* **conformational chirality** -- the mirror is reachable from the original by *rotating
  bonds*. The two are conformers of one isomer, so an encoder that gives them one string is
  **correct** (that is conformer invariance).
* **configurational chirality** -- no amount of bond rotation reaches the mirror. The two are
  genuinely different isomers, and collapsing them is a losslessness failure.

Why not just embed a conformer pool. The obvious implementation -- ``EmbedMultipleConfs`` on
the molecule's graph, then ask whether any conformer fits the mirror -- is **unsound here**,
and the reason is the whole point of this audit. The configurations under test (metal
Delta/Lambda, atropisomerism) are precisely the ones RDKit does *not* carry in the graph, so a
free embed generates BOTH handednesses and the pool fits the mirror every time. The test would
answer "conformational" unconditionally.

**What this module does instead: torsion-space search from the actual structure.** Rotating a
dihedral is a continuous deformation of the given structure; it cannot change any
configuration -- not a tetrahedral centre, not the metal's coordination handedness, not a
biaryl's *configuration* (only its torsion angle, which is the point). So the set of
structures reachable by torsion rotations from the input **is** the conformational orbit of
the input's configuration. We minimise the mirror RMSD over that orbit:

    d_mirror = min over torsion vectors THETA, and over graph automorphisms,
               of  proper-rotation RMSD( mirror(base), coords(THETA) )

If ``d_mirror`` falls below the threshold the mirror lies in the input's own conformational
orbit -> conformational chirality -> collapsing is correct. If it does not, the difference is
configurational.

**The positive control, and why an absolute threshold alone would not do.** A search that
fails proves nothing on its own: it may just be a search that did not converge. So every
verdict is paired with a control on the *same* molecule, the *same* optimiser and the *same*
budget: take the base structure, apply a random torsion vector, and ask the optimiser to
recover it. That target is reachable **by construction**. Three outcomes:

* control converges, mirror converges  -> ``conformational``
* control converges, mirror does not   -> ``configurational``   (the search demonstrably works)
* control does not converge            -> ``inconclusive``      (no evidence either way)

**"No match" is evidence, not proof.** The reported ``budget`` (restarts x sweeps x grid) and
``threshold`` are part of every verdict, and :class:`TorsionVerdict` carries them so a caller
can never quote the conclusion without the search size that produced it.

Two documented modelling choices:

* **Heavy atoms only.** Terminal rotors (-CH3, -OH, -NH2) move hydrogens without changing any
  configuration, and including them would add search dimensions that only add noise. Mirror
  chirality of these species is determined by the heavy-atom skeleton.
* **Automorphism subset during descent.** The full automorphism set is used to score the final
  torsion vector, but only the best ``AUTO_DESCENT_K`` are carried through the inner descent.
  More automorphisms can only *lower* an RMSD, so a subset biases the tool toward
  ``configurational`` -- the direction that costs hand inspection rather than the direction
  that hides a blind spot.

Run:  PYTHONPATH=$PWD/src python -m tools.injectivity.torsion_oracle <file.xyz> ...
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem

from .oracle import MAX_AUTOMORPHISMS, load_mol

#: RMSD (A) below which a torsion vector is called a match. Torsion rotations preserve bond
#: lengths and angles exactly, so a genuinely reachable target optimises to near zero; this is
#: the same 0.5 A the rigid oracle uses, kept identical so the two are directly comparable.
MATCH_THRESHOLD = 0.5

#: number of automorphisms carried through the inner torsion descent (see module docstring).
AUTO_DESCENT_K = 8

#: default search budget.
N_RESTARTS = 12
N_SWEEPS = 4
COARSE_GRID = 12  # 30 deg steps
REFINE_STEPS = (10.0, 3.0)


@dataclass
class TorsionVerdict:
    """A configurational verdict, inseparable from the search that produced it."""

    name: str
    n_torsions: int
    n_heavy: int
    rigid_mirror_rmsd: float
    d_mirror: float
    d_control: float
    threshold: float
    budget: int
    n_autos_total: int
    n_autos_descent: int
    verdict: str  # conformational | configurational | inconclusive | rigid_achiral
    note: str = ""
    extra: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "verdict": self.verdict,
            "d_mirror": round(self.d_mirror, 4),
            "d_control": round(self.d_control, 4),
            "rigid_mirror_rmsd": round(self.rigid_mirror_rmsd, 4),
            "threshold": self.threshold,
            "n_torsions": self.n_torsions,
            "n_heavy": self.n_heavy,
            "budget_evaluations": self.budget,
            "n_autos_total": self.n_autos_total,
            "n_autos_descent": self.n_autos_descent,
            "note": self.note,
            **self.extra,
        }


# --- torsion model -------------------------------------------------------------------


def _heavy_indices(mol: Chem.Mol) -> np.ndarray:
    return np.array([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1], dtype=int)


def rotatable_torsions(mol: Chem.Mol, root: int) -> list[tuple[int, int, np.ndarray, int]]:
    """Rotatable bonds as ``(a, b, moving_atom_indices, depth_from_root)``.

    A bond qualifies when it is acyclic, both ends carry another heavy neighbour (so the
    rotation actually moves heavy atoms), and removing it disconnects the graph. ``moving`` is
    the side that does **not** contain ``root``; ``depth`` orders the rotations so that a
    rotation nearer the root is applied before the ones it carries, which makes the composition
    well defined. Metal-ligand bonds are included -- rotating a monodentate ligand about its
    donor bond is genuine conformational freedom.
    """
    n = mol.GetNumAtoms()
    heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1}
    dist = Chem.GetDistanceMatrix(mol)
    out = []
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i not in heavy or j not in heavy:
            continue
        hi = [nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors() if nb.GetIdx() in heavy]
        hj = [nb.GetIdx() for nb in mol.GetAtomWithIdx(j).GetNeighbors() if nb.GetIdx() in heavy]
        if len(hi) < 2 or len(hj) < 2:
            continue  # terminal on one side: no heavy atom would move
        # side of j reachable without traversing the i-j BOND
        side = _side_without(mol, n, cut=(i, j), start=j)
        if side is None or root in side:
            side = _side_without(mol, n, cut=(i, j), start=i)
            if side is None or root in side:
                continue
            i, j = j, i
        depth = int(dist[root][i]) if np.isfinite(dist[root][i]) else n
        out.append((i, j, np.fromiter(sorted(side), dtype=int), depth))
    out.sort(key=lambda t: (t[3], t[0], t[1]))
    return out


def _side_without(mol: Chem.Mol, n: int, *, cut: tuple[int, int], start: int) -> set[int] | None:
    """Atoms reachable from ``start`` without traversing the ``cut`` bond; ``None`` if all.

    The cut must remove the **bond**, not either atom. Blocking the atom instead is wrong in
    exactly the case that matters here: a metal is a cut vertex, so blocking it detaches every
    ligand from every other, and each metal-donor bond of a *chelate* then looks rotatable even
    though the chelate ring holds it fixed. That bug made fac-Ir(ppy)3's Delta/Lambda mirror
    reachable by "rotating" whole ligands off their own chelate rings.
    """
    a, b = cut
    seen = {start}
    stack = [start]
    while stack:
        u = stack.pop()
        for nb in mol.GetAtomWithIdx(u).GetNeighbors():
            v = nb.GetIdx()
            if v in seen or {u, v} == {a, b}:
                continue
            seen.add(v)
            stack.append(v)
    return None if len(seen) >= n else seen


def _rotate(coords: np.ndarray, a: int, b: int, moving: np.ndarray, deg: float) -> None:
    """Rotate ``moving`` about the current ``a->b`` axis by ``deg``, in place."""
    if deg == 0.0:
        return
    axis = coords[b] - coords[a]
    norm = float(np.linalg.norm(axis))
    if norm < 1e-8:
        return
    k = axis / norm
    t = np.radians(deg)
    c, s = np.cos(t), np.sin(t)
    v = coords[moving] - coords[b]
    coords[moving] = v * c + np.cross(k, v) * s + np.outer(v @ k, k) * (1.0 - c) + coords[b]


def apply_torsions(base: np.ndarray, tors, thetas: np.ndarray) -> np.ndarray:
    """``base`` deformed by the torsion vector ``thetas`` (degrees), root-outward."""
    coords = base.copy()
    for (a, b, moving, _d), deg in zip(tors, thetas, strict=True):
        _rotate(coords, a, b, moving, float(deg))
    return coords


def _signed_dihedral(coords: np.ndarray, i: int, j: int, k: int, ll: int) -> float:
    b0 = coords[i] - coords[j]
    b1 = coords[k] - coords[j]
    b2 = coords[ll] - coords[k]
    b1n = b1 / max(float(np.linalg.norm(b1)), 1e-12)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def dihedral_negating_theta(mol: Chem.Mol, coords: np.ndarray, tors) -> np.ndarray:
    """The torsion vector that negates every rotatable signed dihedral.

    This is the search's most important starting point, and it is not a heuristic. A
    reflection is exactly the operation that preserves all bond lengths and bond angles and
    **negates every signed dihedral**. So if a structure's chirality is purely conformational
    -- every rigid fragment individually achiral, every local configuration non-stereogenic --
    then negating the rotatable dihedrals reproduces the mirror up to a proper rotation and a
    graph automorphism, and the search finds it immediately. If it does not, what is left over
    is precisely the *configurational* part: ring puckers, sp3 centres, the metal's coordination
    handedness -- the dihedrals a torsion rotation cannot touch.

    Random restarts alone are not a substitute. Reaching the mirror generally requires several
    torsions to flip *together*, and one-at-a-time coordinate descent has to climb before it
    descends; on an unsubstituted biphenyl (``EDOQIZ``) 12 random restarts miss it entirely.
    """
    heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1}
    out = np.zeros(len(tors))
    for t, (a, b, _moving, _d) in enumerate(tors):
        ra = [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() != b]
        rb = [n.GetIdx() for n in mol.GetAtomWithIdx(b).GetNeighbors() if n.GetIdx() != a]
        ra = sorted(i for i in ra if i in heavy) or sorted(ra)
        rb = sorted(i for i in rb if i in heavy) or sorted(rb)
        if not ra or not rb:
            continue
        out[t] = -2.0 * _signed_dihedral(coords, ra[0], a, b, rb[0])
    return out


# --- scoring -------------------------------------------------------------------------


def heavy_skeleton(mol: Chem.Mol, heavy: np.ndarray) -> Chem.Mol:
    """``mol`` with every hydrogen deleted, heavy-atom order preserved."""
    rw = Chem.RWMol(mol)
    for idx in sorted((a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1), reverse=True):
        rw.RemoveAtom(idx)
    assert rw.GetNumAtoms() == len(heavy)
    return rw.GetMol()


def _automorphism_perms(mol: Chem.Mol, heavy: np.ndarray, max_autos: int) -> np.ndarray:
    """Graph automorphisms as heavy-atom index permutations, ``(n_auto, n_heavy)``.

    Enumerated on the **heavy skeleton**, not on the H-explicit graph. Enumerating on the
    H-explicit graph is what the rigid oracle does, and on anything bearing methyls it spends
    the whole ``maxMatches`` budget on H permutations that leave the heavy atoms fixed: on
    ``EDOQIZ`` (two tert-butyls) 4000 full-graph matches collapse to **6** distinct heavy
    images, while the heavy skeleton has **864**. Since more automorphisms can only lower an
    RMSD, that starvation makes the mirror look unreachable -- it is what made a freely
    rotating biphenyl read as configurational.

    Heavy-skeleton automorphisms are then filtered to those that also preserve each atom's
    hydrogen count and formal charge, so a -CH2- can never be matched onto a -CH3.
    """
    skel = heavy_skeleton(mol, heavy)
    matches = skel.GetSubstructMatches(
        skel, uniquify=False, useChirality=False, maxMatches=max_autos
    )
    nh = np.array(
        [
            sum(1 for x in mol.GetAtomWithIdx(int(i)).GetNeighbors() if x.GetAtomicNum() == 1)
            for i in heavy
        ]
    )
    chg = np.array([mol.GetAtomWithIdx(int(i)).GetFormalCharge() for i in heavy])
    n = len(heavy)
    perms, seen = [], set()
    for m in matches:
        if len(m) != n:
            continue
        img = np.fromiter(m, dtype=int, count=n)
        if not (np.array_equal(nh[img], nh) and np.array_equal(chg[img], chg)):
            continue
        key = img.tobytes()
        if key in seen:
            continue
        seen.add(key)
        perms.append(img)
    if not perms:
        perms = [np.arange(n)]
    return np.asarray(perms, dtype=int)


def _batch_proper_rmsd(target: np.ndarray, moved: np.ndarray, perms: np.ndarray) -> np.ndarray:
    """Proper-rotation RMSD of ``moved[perm]`` onto ``target`` for every ``perm``.

    Batched Kabsch with the reflection forbidden (``diag(1, 1, det)``), matching
    :func:`tools.injectivity.oracle._kabsch_proper_rmsd` term for term.
    """
    P = target - target.mean(0)
    Q = moved[perms]  # (k, n, 3)
    Q = Q - Q.mean(1, keepdims=True)
    H = np.einsum("kni,nj->kij", Q, P)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(
        np.linalg.det(np.einsum("kij,kjl->kil", Vt.transpose(0, 2, 1), U.transpose(0, 2, 1)))
    )
    D = np.zeros((len(perms), 3, 3))
    D[:, 0, 0] = D[:, 1, 1] = 1.0
    D[:, 2, 2] = d
    R = np.einsum("kij,kjl,klm->kim", Vt.transpose(0, 2, 1), D, U.transpose(0, 2, 1))
    Qr = np.einsum("kni,kji->knj", Q, R)
    return np.sqrt(((Qr - P) ** 2).sum(-1).mean(-1))


class _Scorer:
    """Min proper-rotation RMSD of a torsion-deformed structure onto a fixed target."""

    def __init__(self, base, heavy, perms, target_heavy, tors):
        self.base = base
        self.heavy = heavy
        self.perms_all = perms
        self.target = target_heavy
        self.tors = tors
        # the automorphisms that fit the undeformed structure best lead the descent
        d0 = _batch_proper_rmsd(self.target, base[heavy], perms)
        order = np.argsort(d0)[:AUTO_DESCENT_K]
        self.perms_descent = perms[order]
        self.calls = 0

    def _rmsds(self, thetas, perms):
        coords = apply_torsions(self.base, self.tors, thetas)
        self.calls += 1
        return _batch_proper_rmsd(self.target, coords[self.heavy], perms)

    def descent_score(self, thetas):
        return float(self._rmsds(thetas, self.perms_descent).min())

    def final_score(self, thetas):
        return float(self._rmsds(thetas, self.perms_all).min())

    def residual_profile(self, mol, thetas):
        """Per-atom leftover deviation at ``thetas``, worst first.

        For a ``configurational`` verdict this localises *what* could not be reached, which is
        what a hand inspection needs: a residual concentrated on one motif names the candidate
        axis, one spread over the whole skeleton usually means the search simply lost.
        """
        coords = apply_torsions(self.base, self.tors, thetas)
        moved = coords[self.heavy]
        d = _batch_proper_rmsd(self.target, moved, self.perms_all)
        perm = self.perms_all[int(np.argmin(d))]
        P = self.target - self.target.mean(0)
        Q = moved[perm] - moved[perm].mean(0)
        H = Q.T @ P
        U, _, Vt = np.linalg.svd(H)
        det = float(np.sign(np.linalg.det(Vt.T @ U.T)))
        R = Vt.T @ np.diag([1.0, 1.0, det]) @ U.T
        dev = np.linalg.norm(Q @ R.T - P, axis=1)
        order = np.argsort(-dev)[:6]
        syms = [mol.GetAtomWithIdx(int(self.heavy[i])).GetSymbol() for i in order]
        return {
            "worst_atom_deviation": round(float(dev.max()), 3),
            "worst_atoms": [
                {"element": s, "heavy_pos": int(i), "deviation": round(float(dev[i]), 3)}
                for s, i in zip(syms, order, strict=True)
            ],
        }


def _optimise(scorer: _Scorer, rng, *, restarts, sweeps, grid, seeds=()) -> float:
    """Coordinate descent over the torsion vector from ``seeds`` plus random restarts.

    ``seeds`` carries the informed starting points (see :func:`dihedral_negating_theta`); the
    random restarts stop a seed from being the only thing tried.
    """
    tors = scorer.tors
    k = len(tors)
    if k == 0:
        return scorer.final_score(np.zeros(0)), np.zeros(0)
    best_theta, best = None, np.inf
    coarse = np.linspace(0.0, 360.0, grid, endpoint=False)
    starts = [np.zeros(k), *[np.asarray(s, dtype=float) for s in seeds]]
    for r in range(len(starts) + restarts):
        theta = starts[r].copy() if r < len(starts) else rng.uniform(0.0, 360.0, k)
        cur = scorer.descent_score(theta)
        for _ in range(sweeps):
            improved = False
            for t in range(k):
                keep, base_v = theta[t], cur
                for step in (*coarse, *[keep + d for d in (-10.0, 10.0, -3.0, 3.0)]):
                    theta[t] = step
                    v = scorer.descent_score(theta)
                    if v < cur - 1e-6:
                        cur, keep = v, step
                theta[t] = keep
                improved |= cur < base_v - 1e-6
            if not improved:
                break
        if cur < best:
            best, best_theta = cur, theta.copy()
    return scorer.final_score(best_theta), best_theta


# --- top level -----------------------------------------------------------------------


def configurational_verdict(
    xyz_path: str | Path,
    *,
    charge: int = 0,
    threshold: float = MATCH_THRESHOLD,
    restarts: int = N_RESTARTS,
    sweeps: int = N_SWEEPS,
    grid: int = COARSE_GRID,
    seed: int = 42,
    max_autos: int = MAX_AUTOMORPHISMS,
    mol=None,
    coords=None,
    name: str | None = None,
) -> TorsionVerdict:
    """Is ``mirror(structure)`` reachable from it by bond rotation alone?

    See the module docstring for the method and for what each verdict does and does not
    license. ``mol``/``coords`` let a caller reuse an already-loaded structure.
    """
    path = Path(xyz_path)
    name = name or path.stem
    if mol is None:
        mol, coords = load_mol(path, charge)
    coords = np.asarray(coords, dtype=float)

    heavy = _heavy_indices(mol)
    perms = _automorphism_perms(mol, heavy, max_autos)
    mirror = coords.copy()
    mirror[:, 2] *= -1.0

    rigid = float(_batch_proper_rmsd(coords[heavy], mirror[heavy], perms).min())

    root = int(np.argmax([a.GetAtomicNum() for a in mol.GetAtoms()]))  # heaviest atom
    tors = rotatable_torsions(mol, root)
    # +3 starts that are not random: the undeformed structure and the two signs of the
    # dihedral-negating seed. Counted so the reported budget is not an under-report.
    budget = (restarts + 3) * max(1, sweeps) * max(1, len(tors)) * (grid + 4)

    if rigid <= threshold:
        return TorsionVerdict(
            name=name,
            n_torsions=len(tors),
            n_heavy=len(heavy),
            rigid_mirror_rmsd=rigid,
            d_mirror=rigid,
            d_control=0.0,
            threshold=threshold,
            budget=0,
            n_autos_total=len(perms),
            n_autos_descent=min(AUTO_DESCENT_K, len(perms)),
            verdict="rigid_achiral",
            note="mirror superimposes without any torsion change; no search needed",
        )

    rng = np.random.default_rng(seed)

    # the question: is the mirror in the input's own conformational orbit?
    # the dihedral-negating vector is the physically motivated seed; +/- covers the sign
    # convention of ``_rotate`` without depending on it.
    neg = dihedral_negating_theta(mol, coords, tors) if tors else np.zeros(0)
    mirror_scorer = _Scorer(coords, heavy, perms, mirror[heavy], tors)
    d_mirror, best_theta = _optimise(
        mirror_scorer,
        np.random.default_rng(seed),
        restarts=restarts,
        sweeps=sweeps,
        grid=grid,
        seeds=(neg, -neg),
    )

    # positive control: a target that IS reachable by construction, same optimiser, same budget
    if tors:
        ctrl_target = apply_torsions(coords, tors, rng.uniform(0.0, 360.0, len(tors)))
        d_control, _ = _optimise(
            _Scorer(coords, heavy, perms, ctrl_target[heavy], tors),
            np.random.default_rng(seed + 1),
            restarts=restarts,
            sweeps=sweeps,
            grid=grid,
        )
    else:
        d_control = 0.0

    if d_mirror <= threshold:
        verdict, note = "conformational", "the mirror lies in the structure's own torsion orbit"
    elif not tors:
        verdict = "configurational"
        note = (
            "no rotatable torsion exists (every acyclic single bond is terminal or held by a "
            "chelate ring), so the structure is rigid and the rigid oracle's verdict stands"
        )
    elif d_control <= threshold:
        verdict = "configurational"
        note = (
            "no torsion vector reaches the mirror, while the optimiser recovered a "
            "reachable control target on the same budget"
        )
    else:
        verdict = "inconclusive"
        note = (
            f"the optimiser failed its own positive control (d_control={d_control:.2f} A > "
            f"{threshold} A), so a failure to reach the mirror is not evidence"
        )

    return TorsionVerdict(
        name=name,
        n_torsions=len(tors),
        n_heavy=len(heavy),
        rigid_mirror_rmsd=rigid,
        d_mirror=d_mirror,
        d_control=d_control,
        threshold=threshold,
        budget=budget,
        n_autos_total=len(perms),
        n_autos_descent=min(AUTO_DESCENT_K, len(perms)),
        verdict=verdict,
        note=note,
        extra=(
            mirror_scorer.residual_profile(mol, best_theta)
            if verdict == "configurational" and tors
            else {}
        ),
    )


def _fmt(v: TorsionVerdict) -> str:
    return (
        f"\n## {v.name}\n"
        f"  heavy atoms      : {v.n_heavy}   rotatable torsions: {v.n_torsions}\n"
        f"  automorphisms    : {v.n_autos_total} total, {v.n_autos_descent} led the descent\n"
        f"  rigid mirror RMSD: {v.rigid_mirror_rmsd:.3f} A  (the OLD oracle's number)\n"
        f"  d_mirror         : {v.d_mirror:.3f} A   <- best over the torsion orbit\n"
        f"  d_control        : {v.d_control:.3f} A   <- reachable-by-construction target\n"
        f"  budget           : {v.budget} evaluations, threshold {v.threshold} A\n"
        f"  VERDICT          : {v.verdict.upper()}  -- {v.note}"
    )


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    for p in argv:
        print(_fmt(configurational_verdict(p)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
