"""Guards for the two metal-stereo fixtures Lane 5 (metal-centred Δ/Λ) is validated on.

Why this module exists. The Y2 wave shipped an axial token that was accidentally
*reflection-invariant*; every guard passed because the only fixture exercised the easy
single-axis case. ``fac-Ir(ppy)3`` is the analogous easy case for metal Δ/Λ: its three
chelates are **unsymmetric** (C,N), so a descriptor that actually encodes *fac/mer* rather
than *helicity* would still appear to work on it.

The two fixtures here close that hole from both sides:

* ``ZUMNEC.xyz`` -- tris(catecholato)Mo. A homoleptic tris-bidentate whose two O donors are
  **symmetry-equivalent** (``CanonicalRankAtoms(breakTies=False)`` gives them one rank), so
  the complex has no fac/mer distinction at all and metal helicity (Δ/Λ) is its **sole**
  stereogenic element. A Lane-5 descriptor that leans on ligand asymmetry fails here.
* ``JEGKOW.xyz`` -- Rh(I) square planar with four **different** donor elements
  (N, P, C-carbonyl, I). Exercises RDKit's ``@SP`` path rather than ``@OH``.

Note the deliberate asymmetry in what is asserted. ZUMNEC's mirror **is** a distinct isomer
(Δ/Λ enantiomers). JEGKOW's mirror is **not**, and that is chemically correct rather than a
weak fixture: a square-planar complex is planar, so its coordination plane is a mirror plane
and four different *donors* produce diastereomers, not enantiomers. Reflection is therefore
the wrong distinctness operator for ``@SP`` -- the right one is a donor swap (which permutes
which donors sit trans). ``TestJegkowSquarePlanar`` asserts the achirality and the recovered
``@SP`` permutation; the swap-based distinctness probe rides on the ``swap_donor`` twin
operator in ``tools/injectivity/twin_operators.py``.
"""

import sys
import unittest
from pathlib import Path

from rdkit import Chem

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from oinsmiles.core.constants import TRANSITION_METALS_NUM  # noqa: E402
from oinsmiles.generator3d.clash import vdw_clash_count  # noqa: E402
from tools.injectivity.config_oracle import load_mol, mirror_flip_report  # noqa: E402
from tools.injectivity.oracle import is_distinct_enantiomer  # noqa: E402
from tools.injectivity.twin_collision import (  # noqa: E402
    VERDICT_INVARIANT_OK,
    VERDICT_KEY_BLIND,
    probe_mirror,
)

FIX = _ROOT / "tests" / "fixtures"
ZUMNEC = FIX / "ZUMNEC.xyz"  # Δ/Λ tris-bidentate, symmetry-equivalent donors (P1)
JEGKOW = FIX / "JEGKOW.xyz"  # 4-different-donor square planar (@SP)


def _donor_ranks(mol):
    """Symmetry ranks of the metal's donor atoms, grouped per ligand fragment."""
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM]
    assert len(metals) == 1
    m = metals[0]
    donors = [n.GetIdx() for n in mol.GetAtomWithIdx(m).GetNeighbors()]
    rw = Chem.RWMol(mol)
    rw.RemoveAtom(m)
    out = []
    for frag in Chem.GetMolFrags(rw.GetMol()):
        ds = [d for d in donors if (d if d < m else d - 1) in frag]
        if len(ds) >= 2:
            out.append({ranks[d] for d in ds})
    return out


class TestFixtureGeometryIsClean(unittest.TestCase):
    """The oracle's verdict is meaningless on a bad geometry, so gate the geometry first."""

    def test_no_vdw_clashes(self):
        for path in (ZUMNEC, JEGKOW):
            mol = load_mol(path)
            conf = mol.GetConformer()
            clash, severe, _worst = vdw_clash_count(
                conf.GetPositions(), [a.GetAtomicNum() for a in mol.GetAtoms()]
            )
            self.assertEqual((clash, severe), (0, 0), f"{path.name} has vdW clashes")

    def test_single_metal_center(self):
        for path in (ZUMNEC, JEGKOW):
            mol = load_mol(path)
            metals = [a for a in mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM]
            self.assertEqual(len(metals), 1, f"{path.name} must have exactly one metal")


class TestZumnecTrisBidentate(unittest.TestCase):
    """Δ/Λ tris-bidentate: helicity is the fixture's SOLE stereogenic element."""

    def test_three_chelates_with_symmetry_equivalent_donors(self):
        """Each chelate's two donors share one symmetry rank -> no fac/mer to lean on."""
        groups = _donor_ranks(load_mol(ZUMNEC))
        self.assertEqual(len(groups), 3, "expected three bidentate chelates")
        for g in groups:
            self.assertEqual(len(g), 1, "chelate donors must be symmetry-equivalent")

    def test_mirror_is_a_distinct_isomer(self):
        v = is_distinct_enantiomer(ZUMNEC)
        self.assertTrue(v.distinct, f"Δ/Λ mirror must be distinct (rmsd={v.rmsd})")
        # comfortably clear of the 0.5 A threshold and >10x the achiral controls (~0.05-0.10 A)
        self.assertGreater(v.rmsd, 1.0)

    def test_metal_descriptor_is_octahedral_and_flips(self):
        rep = mirror_flip_report(ZUMNEC)
        self.assertEqual([m.shape for m in rep.base.metal], ["OH"])
        self.assertTrue(rep.metal_flips, "the @OH permutation must flip for the mirror")
        self.assertNotEqual(rep.base.metal[0].permutation, rep.mirror.metal[0].permutation)

    def test_no_other_axis_is_implicated(self):
        """Neither axial nor bound-amine stereo is present, so a Lane-5 pass is unambiguous."""
        rep = mirror_flip_report(ZUMNEC)
        self.assertEqual(rep.base.axial, [])
        self.assertEqual(rep.base.bound_amine, [])

    def test_currently_key_blind(self):
        """Assert-current-behavior: the P1 collapse reproduces on a symmetric tris-chelate."""
        o = probe_mirror(ZUMNEC)
        self.assertTrue(o.oracle_distinct)
        self.assertTrue(o.key_equal, "P1 regression: ZUMNEC enantiomers now key-diverge")
        self.assertEqual(o.verdict, VERDICT_KEY_BLIND)


class TestZumnecAspirational(unittest.TestCase):
    """What Lane 5 must achieve. Fails today; flips to an unexpected success when it lands."""

    @unittest.expectedFailure
    def test_metal_chirality_should_diverge_at_key(self):
        o = probe_mirror(ZUMNEC)
        self.assertFalse(o.key_equal, "Δ/Λ must diverge at the key on a symmetric tris-chelate")


class TestJegkowSquarePlanar(unittest.TestCase):
    """4-different-donor square planar: the ``@SP`` path, and why its mirror is achiral."""

    def test_four_distinct_donor_elements(self):
        mol = load_mol(JEGKOW)
        metal = next(a for a in mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM)
        donors = sorted(n.GetSymbol() for n in metal.GetNeighbors())
        self.assertEqual(len(donors), 4, "square planar must have four donors")
        self.assertEqual(len(set(donors)), 4, f"donors must all differ, got {donors}")

    def test_square_planar_geometry_recovered(self):
        rep = mirror_flip_report(JEGKOW)
        self.assertEqual([m.shape for m in rep.base.metal], ["SP"])
        self.assertNotEqual(rep.base.metal[0].permutation, 0)

    def test_mirror_is_the_same_isomer(self):
        """Correct invariance, not a weak fixture: a square-planar complex is achiral.

        The coordination plane is a mirror plane, so reflection cannot generate a new
        isomer here however different the four donors are. The ``@SP`` permutation
        accordingly does NOT flip -- distinguishing ``@SP`` isomers needs a donor swap.
        """
        rep = mirror_flip_report(JEGKOW)
        self.assertFalse(rep.metal_flips, "@SP should not flip under reflection")
        o = probe_mirror(JEGKOW)
        self.assertFalse(o.oracle_distinct, f"mirror should superimpose (rmsd={o.oracle_rmsd})")
        self.assertTrue(o.raw_equal)
        self.assertEqual(o.verdict, VERDICT_INVARIANT_OK)


if __name__ == "__main__":
    unittest.main()
