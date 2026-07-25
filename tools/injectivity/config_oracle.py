"""Configurational oracle -- recover the *configuration* of each named blind-spot axis
straight from the 3D geometry, independent of the OIN encoder.

Wave 1's :mod:`oracle` answered only *is the mirror a distinct isomer?* (a scalar RMSD).
Wave 2 needs the next thing up: *what is the configuration?* -- a per-axis descriptor that
FLIPS between enantiomers. This is the measurement that reframes the three confirmed blind
spots from "documented-permanent" to "recoverable-but-unemitted": for every one of them the
information is sitting in the input coordinates; the encoder simply discards it.

Three recovery routes, one per axis, matching what the geometry actually supports:

* **P1 metal Δ/Λ** -- RDKit *does* perceive non-tetrahedral stereo from 3D
  (``AssignStereochemistryFrom3D`` -> ``@OH`` / ``@SP`` + a permutation index). The
  permutation flips between enantiomers (fac-Ir(ppy)3: 10 vs 8). We read it directly.
* **P2 axial / atropisomer** -- RDKit will NOT perceive this from pure 3D: per the RDKit
  Book, an atropisomer bond is only marked when a neighbor bond is *wedged*, and the
  configuration is then read from 3D. With no wedge, ``FindPotentialStereo`` returns
  nothing. So we recover it ourselves: the signed biaryl dihedral (BINAP: -70 vs +70).
* **P3 metal-bound 2° amine** -- RDKit clears it as an invertible amine
  (``CHI_UNSPECIFIED`` even after 3D assignment). We recover it as the signed tetrahedral
  volume of the four neighbours of the metal-locked N (POJJOP: -9.4 vs +9.4).

Every descriptor here is a **distinguisher** (it flips for an enantiomer) but is NOT yet a
**canonical token**: the metal permutation is relative to RDKit's neighbour order, the
axial/amine signs are relative to a canonical-rank reference neighbour. Making these
orientation-invariant (so two orientations of one isomer give one string) is the remaining
fix work -- exactly the canonicalization problem v0.4.5 tackles for the achiral string. This
module is measurement only; it changes no encoder output.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem

from oinsmiles.core.constants import TRANSITION_METALS_NUM
from oinsmiles.oin.axial import AxialAxis
from oinsmiles.oin.axial import detect_axial_axes as axial_axes
from oinsmiles.utils.xyz2mol import get_tmc_mol

#: RDKit non-tetrahedral chiral tags that carry a metal Δ/Λ permutation.
_NONTET_TAGS = {
    Chem.ChiralType.CHI_SQUAREPLANAR: "SP",
    Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL: "TB",
    Chem.ChiralType.CHI_OCTAHEDRAL: "OH",
}


@dataclass(frozen=True)
class MetalStereo:
    """Metal-centred non-tetrahedral configuration recovered by RDKit from 3D."""

    atom_idx: int
    shape: str  # "OH" | "SP" | "TB"
    permutation: int  # RDKit's permutation index (flips for the enantiomer)


@dataclass(frozen=True)
class AmineCenter:
    """A metal-locked secondary/tertiary amine N and its signed tetrahedral volume."""

    atom_idx: int
    sign: int  # +1 / -1: the recoverable configuration
    volume: float


@dataclass
class ConfigSignature:
    """The full configurational signature recovered from one structure."""

    metal: list[MetalStereo] = field(default_factory=list)
    axial: list[AxialAxis] = field(default_factory=list)
    bound_amine: list[AmineCenter] = field(default_factory=list)


# --- P1: metal Δ/Λ via RDKit non-tetrahedral perception -----------------------------


def metal_stereo_descriptors(mol: Chem.Mol) -> list[MetalStereo]:
    """Metal-centred non-tetrahedral configurations RDKit perceives from the 3D coords."""
    probe = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(probe)
    out: list[MetalStereo] = []
    for a in probe.GetAtoms():
        if a.GetAtomicNum() not in TRANSITION_METALS_NUM:
            continue
        shape = _NONTET_TAGS.get(a.GetChiralTag())
        if shape is None:
            continue
        perm = a.GetPropsAsDict().get("_chiralPermutation")
        if perm in (None, 0):
            continue
        out.append(MetalStereo(atom_idx=a.GetIdx(), shape=shape, permutation=int(perm)))
    return out


# --- P2: axial / atropisomer via signed biaryl dihedral ------------------------------
# Detection lives in ``oinsmiles.oin.axial`` (single source of truth, shared with the
# encoder's opt-in emit); ``axial_axes`` above is that module's ``detect_axial_axes``.


# --- P3: metal-bound amine via signed tetrahedral volume -----------------------------


def bound_amine_centers(mol: Chem.Mol) -> list[AmineCenter]:
    """Signed tetrahedral volume at each metal-locked 4-coordinate N.

    The metal bond raises the amine to degree 4, freezing the inversion the encoder's
    Zone-A rule assumes; the configuration is the sign of the neighbour tetrahedron's
    volume, neighbours ordered by canonical rank (mirror-consistent)."""
    conf = mol.GetConformer()
    # symmetry-aware ranks (equivalent atoms SHARE a rank) gate stereogenicity;
    # tie-broken ranks give a stable, mirror-consistent neighbour ordering.
    sym_ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    ord_ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=True))
    out: list[AmineCenter] = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 7:
            continue
        nbrs = list(a.GetNeighbors())
        if not any(n.GetAtomicNum() in TRANSITION_METALS_NUM for n in nbrs):
            continue
        if len(nbrs) != 4:
            continue
        # stereogenic only if the four substituents are all distinct: a metal-bound
        # ammine (M,H,H,H) or primary amine (M,H,H,R) has symmetry-equivalent H's and is
        # NOT a stereocentre -- exclude by requiring four distinct *symmetry* ranks.
        if len({sym_ranks[n.GetIdx()] for n in nbrs}) != 4:
            continue
        order = sorted(nbrs, key=lambda n: ord_ranks[n.GetIdx()])
        p = [np.array(conf.GetAtomPosition(n.GetIdx())) for n in order]
        vol = float(np.dot(np.cross(p[1] - p[0], p[2] - p[0]), p[3] - p[0]))
        out.append(
            AmineCenter(atom_idx=a.GetIdx(), sign=1 if vol > 0 else -1, volume=round(vol, 3))
        )
    return out


# --- top level -----------------------------------------------------------------------


def configurational_signature(mol: Chem.Mol) -> ConfigSignature:
    """Recover all three axis descriptors from one structure's 3D coordinates."""
    return ConfigSignature(
        metal=metal_stereo_descriptors(mol),
        axial=[ax for ax in axial_axes(mol) if ax.hindered],
        bound_amine=bound_amine_centers(mol),
    )


def load_mol(xyz_path: str | Path, charge: int = 0) -> Chem.Mol:
    mol, _ = get_tmc_mol(Path(xyz_path), charge, with_stereo=False)
    Chem.SanitizeMol(mol)
    return mol


def _mirror_mol(mol: Chem.Mol) -> Chem.Mol:
    """A z-reflected copy (enantiomer) with identical atom order and graph."""
    m = Chem.Mol(mol)
    conf = m.GetConformer()
    for i in range(m.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(p.x, p.y, -p.z))
    return m


@dataclass
class MirrorFlipReport:
    """Per-axis check that the recovered configuration FLIPS for the mirror enantiomer."""

    metal_flips: bool
    axial_flips: bool
    amine_flips: bool
    base: ConfigSignature
    mirror: ConfigSignature

    @property
    def any_recovered(self) -> bool:
        return self.metal_flips or self.axial_flips or self.amine_flips


def mirror_flip_report(xyz_path: str | Path, charge: int = 0) -> MirrorFlipReport:
    """Load a structure + its z-mirror; report which axis descriptors flip.

    A flip proves the axis is *recoverable* from 3D -- the encoder had the information and
    dropped it. This is the constructive core of the Wave 2 feasibility claim.
    """
    mol = load_mol(xyz_path, charge)
    mir = _mirror_mol(mol)
    b = configurational_signature(mol)
    m = configurational_signature(mir)

    metal_flips = any(
        bm.atom_idx == mm.atom_idx and bm.permutation != mm.permutation
        for bm in b.metal
        for mm in m.metal
    )
    axial_flips = any(
        ba.a1 == ma.a1 and ba.a2 == ma.a2 and ba.sign != ma.sign for ba in b.axial for ma in m.axial
    )
    amine_flips = any(
        bn.atom_idx == mn.atom_idx and bn.sign != mn.sign
        for bn in b.bound_amine
        for mn in m.bound_amine
    )
    return MirrorFlipReport(metal_flips, axial_flips, amine_flips, base=b, mirror=m)


if __name__ == "__main__":
    import sys

    for path in sys.argv[1:]:
        rep = mirror_flip_report(path)
        print(f"{path}")
        print(f"  metal Δ/Λ   flips: {rep.metal_flips}  base={rep.base.metal}")
        print(f"  axial       flips: {rep.axial_flips}  base={rep.base.axial}")
        print(f"  bound amine flips: {rep.amine_flips}  base={rep.base.bound_amine}")
